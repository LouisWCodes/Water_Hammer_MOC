# Complete Water Hammer Analysis using Method of Characteristics (MOC)
# This method is a 1-Dimensional analysis technique for fluid transients using the MOC
# It is not intended to be a final answer rather a useful tool for engineers. Engineers should
# use it in conjunction with other methods and tools - such as Bentley's Hammer Software or AFT's Impulse
# for best results.
# Louis Walker, Process Engineer, 2025

from math import exp, sqrt, pi 
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------Pipe & Boundary Parameters--------------------------
length = 5500       # Pipe length (m)
diam = 0.4          # Pipe diameter (m)
wave_speed = 950  # Wave speed (m/s) - from paper
f = 0.009            # Friction factor
n_nodes = 215      # Number of nodes
H_upstream = 70 # Upstream reservoir head (m)
H_downstream = 60 # Downstream reservoir head (m)
end_time = 600  # Total simulation time (s)
xkvalve_0 = 0.1    # Initial valve resistance - fully open butterfly valve 

grav = 9.81         # Gravity (m/s²)
density = 1000      # Water density (kg/m³)

#-----------------------------Class Definition---------------------------------------
class WaterhammerMOC:
    def __init__(self, length, diam, wave_speed, f, n_nodes, store_interval = 1):
        # Pipe properties
        self.length = length
        self.diam = diam
        self.wave_speed = wave_speed
        self.f = f
        self.n_nodes = n_nodes
        self.TIME = 0.0  # Initialize time
        
        # Calculate cross-sectional area once and store
        self.area = np.pi * self.diam**2 / 4
        
        # Boundary conditions
        self.H_upstream = H_upstream
        self.H_downstream = H_downstream

        # Physical constants
        self.grav = grav
        self.density = density

        # Simulation time
        self.start_time = 0.0
        self.end_time = end_time
        self.dx = self.length / (self.n_nodes - 1) # Sets spatial (distance) step size
        self.dt = self.dx / self.wave_speed # Sets temporal (time) step size
        self.n_steps = int(self.end_time / self.dt) # Sets number of time steps
        
        # Courant number check
        self.courant = self.wave_speed * self.dt / self.dx
        if abs(self.courant - 1.0) > 0.01:
           print(f"Warning: Courant number = {self.courant:.4f} (should be ≈1.0)") #Courant number should be approximately 1.0
           # CFL is a necessary condition for convergence while solving hyperbolic PDEs numerically. CFL should be <= 1 for 1-D problems.
           print(f"Consider adjusting n_nodes for better accuracy")
            
        # Valve parameters
        self.xkvalve_0 = xkvalve_0
        self.ANGLE1 = 90.0
        self.ANGLE2 = 57
        self.ANGLE3 = 18
        self.ANGLE4 = 0.0
        self.T1 = 0.0
        self.T2 = 16.5
        self.T3 = 23
        self.T4 = 50

        # Numerical arrays
        self.H = np.zeros(self.n_nodes)          # Current head (m)
        self.H_old = np.zeros(self.n_nodes)      # Previous head (m)
        self.V = np.zeros(self.n_nodes)          # Current velocity (m/s)
        self.V_old = np.zeros(self.n_nodes)      # Previous velocity (m/s)
        self.Q_current = np.zeros(self.n_nodes)  # Current flow rate (m³/s)
        self.Q_old = np.zeros(self.n_nodes)      # Previous flow rate (m³/s)

        # Characteristic arrays
        self.CP = np.zeros(self.n_nodes)          # C+ characteristic (m)
        self.CM = np.zeros(self.n_nodes)          # C- characteristic (m)

        # Results storage
        self.time_history = [] 
        self.H_history = []
        self.V_history = []
        self.Q_history = []

        # Coefficients 
        self.B = self.wave_speed / self.grav # B 
        self.XF = self.f * self.dx / (2 * self.grav * self.diam * self.area**2) # Valve friction term

        # Store results at intervals for long simulations
        self.store_interval = store_interval
        self.store_counter = 0
        
        # Initialize steady state
        self.initialise_steady_state()

    #------------------Initial Steady-State------------------
    def initialise_steady_state(self):
        """Set up initial steady flow conditions"""
        K = self.xkvalve_0  # Valve resistance coefficient
        
        # Match Fortran calculation exactly: QX=(HIN-HOUT)/(XKVALVE/AREA/19.62+XF*XL/DELTX)
        valve_resistance_term = K / (self.area * 2 * self.grav)  # XKVALVE/AREA/19.62
        pipe_resistance_term = self.XF * self.length / self.dx   # XF*XL/DELTX
        
        total_resistance = valve_resistance_term + pipe_resistance_term
        
        if total_resistance <= 0:
            raise ValueError("Non-positive total resistance in initial flow calculation")
        
        # Calculate initial flow parameters
        delta_H = self.H_upstream - self.H_downstream
        
        # This calculates Q^2, then takes square root for Q 
        Q_squared = delta_H / total_resistance
        self.Q_0 = np.sqrt(Q_squared)  # Flow rate (m³/s)
        self.V_0 = self.Q_0 / self.area  # Velocity (m/s)
        
        # Set uniform velocity and flow rate
        self.V[:] = self.V_0
        self.V_old[:] = self.V_0
        self.Q_current[:] = self.Q_0
        self.Q_old[:] = self.Q_0

        # Linear head distribution accounting for pipe friction only
        # (valve drop occurs at boundary, not distributed along pipe)
        x = np.linspace(0, self.length, self.n_nodes)
        head_loss_gradient = self.f * self.V_0**2 / (2 * self.grav * self.diam)
        self.H = self.H_upstream - head_loss_gradient * x
        self.H_old[:] = self.H[:]

        # Store initial conditions (t=0)
        self.store_results()
        
        # Print initial conditions for validation
        print(f"Initial conditions set:")
        print(f"  Initial velocity: {self.V_0:.3f} m/s")
        print(f"  Initial flow rate: {self.Q_0:.4f} m³/s")
        print(f"  Valve resistance term: {valve_resistance_term:.6f} s²/m⁵")
        print(f"  Pipe resistance term: {pipe_resistance_term:.6f} s²/m⁵")
        print(f"  Head at upstream: {self.H[0]:.2f} m")
        print(f"  Head at downstream: {self.H[-1]:.2f} m")
        print(f"  Time step: {self.dt:.4f} s")
        print(f"  Number of time steps: {self.n_steps}")

    #------------------Calculate Characteristics------------------
    def calculate_characteristics(self):
        """Calculate C+ and C- characteristics using OLD timestep values"""
        for i in range(self.n_nodes):
            dynamic_head = self.B * self.Q_old[i] / self.area # This calculates dyamic head i.e. this is B*VP
            friction_head = self.XF * self.Q_old[i] * abs(self.Q_old[i]) 
            self.CP[i] = self.H_old[i] + dynamic_head - friction_head # Friction is negative in the characteristic + positive (C+) direction
            self.CM[i] = self.H_old[i] - dynamic_head + friction_head # Friction is positive in the characteristic - negative (C-) direction

    #----------------Solving the Interior Points------------------
    def solve_interior_points(self):
        """Solve interior nodes using the Method of Characteristics (MOC)
        
        Linearized 1D unsteady flow equations for water hammer:
            dH/dt + B * dQ/dx = 0
            dQ/dt + (1/B) * dH/dx + friction_term = 0

        Using C+ and C- characteristics:
            H_new = 0.5 * (H_C+ from upstream + H_C- from downstream)
            V_new = 0.5 * (H_C+ - H_C-) / B
            Q_new = V_new * AREA
        """
        for i in range(1, self.n_nodes-1):  # Loop over interior nodes
            # Compute new head at next time step
            self.H[i] = 0.5 * (self.CP[i-1] + self.CM[i+1])
            # Compute new flow rate ( Wylie & Streeter) update flow using the node itself, not i-1 for CP (standard linearized MOC)
            ''' Standardized linear MOC can produce unstable oscillations after many iterations where CFL is near 1'''
            self.Q_current[i] = 0.5 * (self.CP[i]- self.CM[i+1]) / self.B * self.area

    #-----------------Calculate valve coefficient-----------------
    def get_valve_coefficient(self):
        """Calculate valve resistance coefficient using vectorized piecewise-linear interpolation."""

        # Define valve times and angles as arrays
        times = np.array([0.0, self.T1, self.T2, self.T3, self.T4])
        angles = np.array([self.ANGLE1, self.ANGLE1, self.ANGLE2, self.ANGLE3, self.ANGLE4])

        # Interpolate current valve angle at self.TIME
        if self.TIME >= self.T4:
            return float('inf')  # Valve fully closed

        current_angle = np.interp(self.TIME, times, angles)

        # Calculate empirical valve coefficient
        if current_angle <= 0:
            return float('inf')
            
        valve_coefficient = exp((3.78 - 0.038 * current_angle) * 2.3)
        return valve_coefficient

    def apply_boundary_conditions(self):
        """Apply upstream and downstream boundary conditions"""
        
        # Upstream boundary - reservoir head 
        self.H[0] = self.H_upstream # Set initial head to head present in reservoir 
        self.Q_current[0] = (self.H_upstream - self.CM[1]) * self.area / self.B # Set initial flow as approximation from upstream head and downstream characteristic

        # Define Characteristics at last interior node
        CP_last = self.CP[-2]  # characteristic from last interior node
        CM_last = self.CM[-2]  # optional if needed

        # Downstream boundary - valve
        valve_coeff = self.get_valve_coefficient()
        
        if valve_coeff == float('inf'):
            # Valve fully closed
            self.Q_current[-1] = 0.0
            self.H[-1] = self.CP[-2]
        else:
            # Simplified valve boundary using CP characteristic
            # Based on: CP = H + BQ/A and valve equation H_loss = K*Q²/(2gA²)
            
            # Valve equation: H_downstream + K*Q²/(2gA²) = H_valve
            # CP equation: H_valve = CP - BQ/A
            # Combined: H_downstream + K*Q²/(2gA²) = CP - BQ/A
            
            # Rearranging: (K/(2gA²))Q² + (B/A)Q + (H_downstream - CP) = 0
            a = valve_coeff / (2 * self.grav * self.area**2)
            b = self.B / self.area
            c = self.H_downstream - self.CP[-2]
            
            # Standard quadratic formula discriminant
            discriminant = b**2 - 4*a*c
            # We should also handle when discriminant becomes very small else numerical instability will be introduced
            # tolerance depends on user preference but lets choose something suitable for now
            tolerance = 1e-10
            # If discriminant is less than tolerance sent the new flowrate to -b/2a (one distinct and repeated root)
            if abs (discriminant) < tolerance:
                self.Q_current[-1] = -b / (2*a)
                self.H[-1] = CP_last - self.B * self.Q_current[-1] / self.area
            # Continue simulation as normal when the discriminant is positive and reasonable
            elif discriminant >= tolerance:
                Q1 = (-b + sqrt(discriminant)) / (2*a)
                Q2 = (-b - sqrt(discriminant)) / (2*a)
                # choose root closest to previous flow (prevents overshoot)
                if abs(Q1 - self.Q_old[-1]) < abs(Q2 - self.Q_old[-1]):
                    self.Q_current[-1] = Q1
                else:
                    self.Q_current[-1] = Q2
                # update head
                self.H[-1] = CP_last - self.B * self.Q_current[-1] / self.area
            else:
                self.Q_current[-1] = 0.0
                self.H[-1] = CP_last

    #-------------Updates Arrays with New Values---------------
    def update_arrays(self):
        """Update arrays for next timestep"""
        # Calculate velocities
        self.V[:] = self.Q_current[:] / self.area
        
        # Store current values as old values for next iteration
        self.H_old[:] = self.H[:]
        self.V_old[:] = self.V[:]
        self.Q_old[:] = self.Q_current[:]
    #---------Function checks for large changes in head which may cause extreme oscillations---------------
    def check_for_jumps(self, step_name):
        """Check for sudden value changes that could cause oscillations"""
        if hasattr(self, '_last_H'):
            max_H_change = np.max(np.abs(self.H - self._last_H))
            max_Q_change = np.max(np.abs(self.Q_current - self._last_Q))
            
            if max_H_change > 5.0:  # Adjust threshold as needed
                print(f"LARGE HEAD JUMP at {step_name}: {max_H_change:.3f} m")
            if max_Q_change > 0.1:  # Adjust threshold as needed  
                print(f"LARGE FLOW JUMP at {step_name}: {max_Q_change:.4f} m³/s")
        
        self._last_H = self.H.copy()
        self._last_Q = self.Q_current.copy()
    #-------------Stores values for use in post-processing/plotting---------------
    def store_results(self, force_store=False):
        """Store results at current time step"""
        if force_store or self.store_counter % self.store_interval == 0:
            self.time_history.append(self.TIME)
            self.H_history.append(self.H.copy())
            self.V_history.append(self.V.copy())
            self.Q_history.append(self.Q_current.copy())
        self.store_counter += 1
          
    #-------------Main Simulation Loop---------------
    def run_simulation(self):
        """Main simulation loop"""
        print("Starting Water Hammer Simulation...")
        
        step_count = 0
        max_steps = int(self.end_time / self.dt)
        
        # Main time loop
        while step_count < max_steps and self.TIME < self.end_time:
            # Advance time
            self.TIME += self.dt
            step_count += 1
            
            # MOC solution steps
            self.calculate_characteristics()    # Use old values
            self.solve_interior_points()       # Solve interior nodes
            self.apply_boundary_conditions()   # Apply boundary conditions
            self.update_arrays()               # Update old arrays ONCE
            self.store_results()               # Store results
            
            # Progress indicator
            if step_count % 100 == 0:
                print(f"  Time: {self.TIME:.3f}s, Steps: {step_count}")
                print(f"    Outlet flow: {self.Q_current[-1]:.4f} m³/s")
                print(f"    Outlet head: {self.H[-1]:.2f} m")
        
        print(f"Simulation Complete!")
            # Force storage of final result if not already stored
        if (step_count - 1) % self.store_interval != 0:
            self.store_results(force_store=True)
        
        print(f"  Total time steps: {step_count}")
    # ... rest of completion message
        print(f"  Total time steps: {step_count}")
        print(f"  Final time: {self.TIME:.6f}s")
        print(f"  Results stored: {len(self.time_history)} points")
        
    #-------------Results Plotting----------------
    def plot_results(self):
        """Plot the simulation results"""
        
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
        
        # Plot head at downstream end
        downstream_heads = [h[-1] for h in self.H_history]
        ax1.plot(self.time_history, downstream_heads, 'b-', linewidth=2)
        ax1.set_ylabel('Head at Outlet (m)')
        ax1.set_title('Water Hammer Analysis Results')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, self.end_time)
        
        # Plot flow at downstream end
        downstream_flows = [q[-1] for q in self.Q_history]
        ax2.plot(self.time_history, downstream_flows, 'r-', linewidth=2)
        ax2.set_ylabel('Flow at Outlet (m³/s)')
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(0, self.end_time)
        
        # Plot valve angle
        valve_angles = []
        for t in self.time_history:
            times = np.array([0.0, self.T1, self.T2, self.T3, self.T4])
            angles = np.array([self.ANGLE1, self.ANGLE1, self.ANGLE2, self.ANGLE3, self.ANGLE4])
            if t >= self.T4:
                valve_angles.append(self.ANGLE4)
            else:
                valve_angles.append(np.interp(t, times, angles))
        
        ax3.plot(self.time_history, valve_angles, 'g-', linewidth=2)
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Valve Angle (degrees)')
        ax3.grid(True, alpha=0.3)
        ax3.set_xlim(0, self.end_time)
        
        plt.tight_layout()
        plt.show()

#------------------Test Setup------------------
# Create and run the simulation
if __name__ == "__main__":
    # Create MOC solver
    solver = WaterhammerMOC(length, diam, wave_speed, f, n_nodes)
    
    # Run the simulation
    solver.run_simulation()
    
    # Plot results
    solver.plot_results()