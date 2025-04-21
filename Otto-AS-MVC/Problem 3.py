import sys
import os
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from Air import Air  # Corrected import to match Air.py


# Debug the path (optional, retained for troubleshooting)
script_dir = os.path.dirname(os.path.abspath(__file__))
otto_dir = os.path.abspath(os.path.join(script_dir, '..', 'Otto-AS-MVC'))
print(f"Adding to sys.path: {otto_dir}")
sys.path.append(otto_dir)

# Check if Air.py exists
if os.path.exists(os.path.join(otto_dir, 'Air.py')):
    print("Air.py found!")
else:
    print("Air.py not found in the specified directory!")


class DieselModel:
    """Models a Diesel cycle for thermodynamic analysis using molar-based properties.

    Calculates the four states of the Diesel cycle (isentropic compression,
    constant-pressure heat addition, isentropic expansion, constant-volume heat
    rejection) using the Air class for thermodynamic properties.

    Attributes:
        r (float): Compression ratio (V1/V2).
        rc (float): Cutoff ratio (V3/V2).
        T1 (float): Initial temperature in Kelvin.
        P1 (float): Initial pressure in Pa.
        V1 (float): Initial specific volume in m³/mol (set by Air).
        air (Air): Air model for thermodynamic calculations.
        results (dict): Dictionary containing cycle results (temperatures, pressures, etc.).
    """

    def __init__(self):
        """Initializes the Diesel model with default parameters."""
        self.r = 18              # Default compression ratio
        self.rc = 2              # Default cutoff ratio
        self.T1 = 300            # Initial temperature in Kelvin
        self.P1 = 0.1e6          # Initial pressure in Pa
        self.V1 = 0.001          # Placeholder initial volume (overwritten)
        self.air = Air()         # Air model for thermodynamic calculations
        self.results = {}        # Store cycle results

    def calculate_cycle(self):
        """Calculates the Diesel cycle states, work, heat, and efficiency.

        Returns:
            tuple: (results, pv_points) where results is a dictionary of cycle
                   parameters and pv_points is a list of (volume, pressure) tuples
                   for the P-V diagram.
        """
        # State 1: Initial state
        self.air.set(T=self.T1, P=self.P1)  # Set initial state
        s1 = self.air.s                     # Entropy in J/mol·K
        u1 = self.air.u                     # Internal energy in J/mol
        v1 = self.air.v                     # Specific volume in m³/mol
        self.V1 = v1                        # Update V1 to specific volume

        # State 2: Isentropic compression
        V2 = self.V1 / self.r               # Compressed volume
        self.air.set(s=s1, v=V2)            # Constant entropy
        T2 = self.air.T                     # Temperature in K
        P2 = self.air.P                     # Pressure in Pa
        u2 = self.air.u                     # Internal energy in J/mol

        # State 3: Constant-pressure heat addition
        P3 = P2                             # Constant pressure
        V3 = self.rc * V2                   # Volume based on cutoff ratio
        self.air.set(P=P3, v=V3)            # Set state
        T3 = self.air.T                     # Temperature in K
        s3 = self.air.s                     # Entropy in J/mol·K
        u3 = self.air.u                     # Internal energy in J/mol

        # State 4: Isentropic expansion
        V4 = self.V1                        # Return to initial volume
        self.air.set(s=s3, v=V4)            # Constant entropy
        T4 = self.air.T                     # Temperature in K
        P4 = self.air.P                     # Pressure in Pa
        u4 = self.air.u                     # Internal energy in J/mol

        # Calculate work, heat, and efficiency (per mole)
        q_in = u3 - u2                      # Heat input (constant pressure)
        q_out = u1 - u4                     # Heat rejection (constant volume)
        w_net = q_in - q_out                # Net work
        efficiency = w_net / q_in * 100 if q_in != 0 else 0  # Efficiency in %

        # Store results
        self.results = {
            'T1': self.T1, 'P1': self.P1, 'V1': self.V1,
            'T2': T2, 'P2': P2, 'V2': V2,
            'T3': T3, 'P3': P3, 'V3': V3,
            'T4': T4, 'P4': P4, 'V4': V4,
            'q_in': q_in, 'q_out': q_out,
            'w_net': w_net, 'efficiency': efficiency
        }

        # Generate P-V diagram points
        pv_points = []
        # Isentropic compression (1-2)
        v = np.linspace(self.V1, V2, 50)
        for vol in v:
            self.air.set(s=s1, v=vol)  # Constant entropy
            pv_points.append((vol, self.air.P))
        # Constant-pressure heat addition (2-3)
        v = np.linspace(V2, V3, 50)
        for vol in v:
            self.air.set(P=P2, v=vol)  # Constant pressure
            pv_points.append((vol, self.air.P))
        # Isentropic expansion (3-4)
        v = np.linspace(V3, V4, 50)
        for vol in v:
            self.air.set(s=s3, v=vol)  # Constant entropy
            pv_points.append((vol, self.air.P))
        # Close the cycle (4-1)
        pv_points.append((V4, P4))
        pv_points.append((self.V1, self.P1))

        return self.results, pv_points


class DieselView:
    """Manages the Tkinter GUI for displaying Diesel cycle results and P-V diagram.

    Provides input fields for cycle parameters, displays results, and plots the P-V
    diagram using Matplotlib.

    Attributes:
        root (tk.Tk): Main Tkinter window.
        controller (DieselController): Controller for handling calculations.
        r_entry (ttk.Entry): Entry field for compression ratio.
        rc_entry (ttk.Entry): Entry field for cutoff ratio.
        t1_entry (ttk.Entry): Entry field for initial temperature.
        p1_entry (ttk.Entry): Entry field for initial pressure.
        result_labels (dict): Dictionary of readonly entry fields for results.
        figure (matplotlib.figure.Figure): Matplotlib figure for plotting.
        ax (matplotlib.axes.Axes): Matplotlib axes for P-V diagram.
        canvas (FigureCanvasTkAgg): Matplotlib canvas for Tkinter integration.
    """

    def __init__(self, root, controller):
        """Initializes the view with input fields, result displays, and plot canvas.

        Args:
            root (tk.Tk): Main Tkinter window.
            controller (DieselController): Controller for handling calculations.
        """
        self.root = root
        self.controller = controller
        self.root.title("Diesel Cycle Analysis")

        # Input frame
        input_frame = ttk.LabelFrame(self.root, text="Input Parameters")
        input_frame.grid(row=0, column=0, padx=10, pady=5, sticky="ew")

        # Compression ratio input
        ttk.Label(input_frame, text="Compression Ratio (r):").grid(row=0, column=0, padx=5, pady=5)
        self.r_entry = ttk.Entry(input_frame)
        self.r_entry.insert(0, "18")  # Default value
        self.r_entry.grid(row=0, column=1, padx=5, pady=5)

        # Cutoff ratio input
        ttk.Label(input_frame, text="Cutoff Ratio (rc):").grid(row=1, column=0, padx=5, pady=5)
        self.rc_entry = ttk.Entry(input_frame)
        self.rc_entry.insert(0, "2")  # Default value
        self.rc_entry.grid(row=1, column=1, padx=5, pady=5)

        # Initial temperature input
        ttk.Label(input_frame, text="T1 (K):").grid(row=2, column=0, padx=5, pady=5)
        self.t1_entry = ttk.Entry(input_frame)
        self.t1_entry.insert(0, "300")  # Default value
        self.t1_entry.grid(row=2, column=1, padx=5, pady=5)

        # Initial pressure input
        ttk.Label(input_frame, text="P1 (MPa):").grid(row=3, column=0, padx=5, pady=5)
        self.p1_entry = ttk.Entry(input_frame)
        self.p1_entry.insert(0, "0.1")  # Default value
        self.p1_entry.grid(row=3, column=1, padx=5, pady=5)

        # Calculate button
        ttk.Button(self.root, text="Calculate", command=self.controller.calculate).grid(row=1, column=0, pady=10)

        # Output frame
        output_frame = ttk.LabelFrame(self.root, text="Results")
        output_frame.grid(row=2, column=0, padx=10, pady=5, sticky="ew")

        # Initialize result display fields
        self.result_labels = {}
        row = 0
        for param in ['T1', 'T2', 'T3', 'T4']:
            ttk.Label(output_frame, text=f"{param} (K):").grid(row=row, column=0, padx=5, pady=2)
            self.result_labels[param] = ttk.Entry(output_frame, state='readonly')
            self.result_labels[param].grid(row=row, column=1, padx=5, pady=2)
            row += 1

        for param in ['P1', 'P2', 'P3', 'P4']:
            ttk.Label(output_frame, text=f"{param} (Pa):").grid(row=row, column=0, padx=5, pady=2)
            self.result_labels[param] = ttk.Entry(output_frame, state='readonly')
            self.result_labels[param].grid(row=row, column=1, padx=5, pady=2)
            row += 1

        for param in ['q_in', 'q_out', 'w_net']:
            ttk.Label(output_frame, text=f"{param} (J):").grid(row=row, column=0, padx=5, pady=2)
            self.result_labels[param] = ttk.Entry(output_frame, state='readonly')
            self.result_labels[param].grid(row=row, column=1, padx=5, pady=2)
            row += 1

        ttk.Label(output_frame, text="Efficiency (%):").grid(row=row, column=0, padx=5, pady=2)
        self.result_labels['efficiency'] = ttk.Entry(output_frame, state='readonly')
        self.result_labels['efficiency'].grid(row=row, column=1, padx=5, pady=2)

        # Initialize Matplotlib plot
        self.figure, self.ax = plt.subplots(figsize=(5, 4))
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=1, rowspan=3, padx=10, pady=5)

    def update_results(self, results):
        """Updates the result display fields with cycle calculations.

        Args:
            results (dict): Dictionary containing cycle parameters (T1, P1, q_in, etc.).
        """
        for param, value in results.items():
            self.result_labels[param].configure(state='normal')  # Enable editing
            self.result_labels[param].delete(0, tk.END)         # Clear field
            self.result_labels[param].insert(0, f"{value:.2f}") # Insert new value
            self.result_labels[param].configure(state='readonly') # Lock field

    def plot_pv_diagram(self, pv_points):
        """Plots the P-V diagram for the Diesel cycle.

        Args:
            pv_points (list): List of (volume, pressure) tuples for the cycle.
        """
        self.ax.clear()  # Clear previous plot
        v, p = zip(*pv_points)  # Unpack volume and pressure
        self.ax.plot(v, p, 'b-', label='Diesel Cycle')  # Plot cycle
        self.ax.set_xlabel('Volume (m³/mol)')           # Label x-axis
        self.ax.set_ylabel('Pressure (Pa)')             # Label y-axis
        self.ax.set_title('P-V Diagram')                # Set title
        self.ax.grid(True)                              # Add grid
        self.ax.legend()                                # Show legend
        self.canvas.draw()                              # Redraw canvas


class DieselController:
    """Controller for managing Diesel cycle calculations and GUI updates.

    Coordinates between the DieselModel for calculations and DieselView for
    GUI interactions, handling input validation and result display.

    Attributes:
        model (DieselModel): Diesel cycle model instance.
        view (DieselView): View instance for GUI updates.
    """

    def __init__(self, model, view):
        """Initializes the controller with a model and view.

        Args:
            model (DieselModel): Diesel cycle model.
            view (DieselView): Diesel cycle view.
        """
        self.model = model
        self.view = view

    def calculate(self):
        """Performs Diesel cycle calculations based on user inputs.

        Validates inputs, updates the model, and triggers view updates.

        Raises:
            ValueError: If inputs are invalid (non-numeric, non-positive).
        """
        try:
            # Validate and get inputs
            r = float(self.view.r_entry.get())
            rc = float(self.view.rc_entry.get())
            T1 = float(self.view.t1_entry.get())
            P1 = float(self.view.p1_entry.get()) * 1e6  # Convert MPa to Pa
            if any(x <= 0 for x in [r, rc, T1, P1]):
                raise ValueError("All inputs must be positive")

            # Update model parameters
            self.model.r = r
            self.model.rc = rc
            self.model.T1 = T1
            self.model.P1 = P1

            # Perform calculations and update view
            results, pv_points = self.model.calculate_cycle()
            self.view.update_results(results)
            self.view.plot_pv_diagram(pv_points)
        except ValueError as e:
            tk.messagebox.showerror("Input Error", f"Invalid input: {str(e)}")
        except Exception as e:
            tk.messagebox.showerror("Calculation Error", f"Calculation failed: {str(e)}")


if __name__ == "__main__":
    root = tk.Tk()  # Initialize Tkinter window
    model = DieselModel()  # Create model
    view = DieselView(root, None)  # Create view (controller set later)
    controller = DieselController(model, view)  # Create controller
    view.controller = controller  # Assign controller to view
    root.mainloop()  # Run Tkinter event loop