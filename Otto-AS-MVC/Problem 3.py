import sys
import os

# Debug the path
script_dir = os.path.dirname(os.path.abspath(__file__))
otto_dir = os.path.join(script_dir, '..', 'Otto-AS-MVC')
otto_dir = os.path.abspath(otto_dir)
print(f"Adding to sys.path: {otto_dir}")
sys.path.append(otto_dir)

# Check if Air.py exists
if os.path.exists(os.path.join(otto_dir, 'Air.py')):
    print("Air.py found!")
else:
    print("Air.py not found in the specified directory!")

import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from air import Air

# Model
class DieselModel:
    def __init__(self):
        self.r = 18
        self.rc = 2
        self.T1 = 300
        self.P1 = 0.1e6
        self.V1 = 0.001
        self.air = Air()
        self.results = {}

    def calculate_cycle(self):
        # State 1
        self.air.set(T=self.T1, P=self.P1)
        s1 = self.air.s
        u1 = self.air.u
        v1 = self.air.v
        self.V1 = v1

        # State 2: Isentropic compression
        V2 = self.V1 / self.r
        self.air.set(s=s1, v=V2)
        T2 = self.air.T
        P2 = self.air.P
        u2 = self.air.u

        # State 3: Constant pressure heat addition
        P3 = P2
        V3 = self.rc * V2
        self.air.set(P=P3, v=V3)
        T3 = self.air.T
        s3 = self.air.s
        u3 = self.air.u

        # State 4: Isentropic expansion
        V4 = self.V1
        self.air.set(s=s3, v=V4)
        T4 = self.air.T
        P4 = self.air.P
        u4 = self.air.u

        # Calculate work, heat, and efficiency
        q_in = u3 - u2
        q_out = u1 - u4
        w_net = q_in - q_out
        efficiency = w_net / q_in * 100

        self.results = {
            'T1': self.T1, 'P1': self.P1, 'V1': self.V1,
            'T2': T2, 'P2': P2, 'V2': V2,
            'T3': T3, 'P3': P3, 'V3': V3,
            'T4': T4, 'P4': P4, 'V4': V4,
            'q_in': q_in, 'q_out': q_out,
            'w_net': w_net, 'efficiency': efficiency
        }

        # P-V diagram data
        pv_points = []
        v = np.linspace(self.V1, V2, 50)
        for vol in v:
            self.air.set(s=s1, v=vol)
            pv_points.append((vol, self.air.P))
        v = np.linspace(V2, V3, 50)
        for vol in v:
            self.air.set(P=P2, v=vol)
            pv_points.append((vol, self.air.P))
        v = np.linspace(V3, V4, 50)
        for vol in v:
            self.air.set(s=s3, v=vol)
            pv_points.append((vol, self.air.P))
        pv_points.append((V4, P4))
        pv_points.append((self.V1, self.P1))

        return self.results, pv_points

# View
class DieselView:
    def __init__(self, root, controller):
        self.root = root
        self.controller = controller
        self.root.title("Diesel Cycle Analysis")

        # Input frame
        input_frame = ttk.LabelFrame(self.root, text="Input Parameters")
        input_frame.grid(row=0, column=0, padx=10, pady=5, sticky="ew")

        ttk.Label(input_frame, text="Compression Ratio (r):").grid(row=0, column=0, padx=5, pady=5)
        self.r_entry = ttk.Entry(input_frame)
        self.r_entry.insert(0, "18")
        self.r_entry.grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(input_frame, text="Cutoff Ratio (rc):").grid(row=1, column=0, padx=5, pady=5)
        self.rc_entry = ttk.Entry(input_frame)
        self.rc_entry.insert(0, "2")
        self.rc_entry.grid(row=1, column=1, padx=5, pady=5)

        ttk.Label(input_frame, text="T1 (K):").grid(row=2, column=0, padx=5, pady=5)
        self.t1_entry = ttk.Entry(input_frame)
        self.t1_entry.insert(0, "300")
        self.t1_entry.grid(row=2, column=1, padx=5, pady=5)

        ttk.Label(input_frame, text="P1 (MPa):").grid(row=3, column=0, padx=5, pady=5)
        self.p1_entry = ttk.Entry(input_frame)
        self.p1_entry.insert(0, "0.1")
        self.p1_entry.grid(row=3, column=1, padx=5, pady=5)

        # Calculate button
        ttk.Button(self.root, text="Calculate", command=self.controller.calculate).grid(row=1, column=0, pady=10)

        # Output frame
        output_frame = ttk.LabelFrame(self.root, text="Results")
        output_frame.grid(row=2, column=0, padx=10, pady=5, sticky="ew")

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

        # Matplotlib graph
        self.figure, self.ax = plt.subplots(figsize=(5, 4))
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=1, rowspan=3, padx=10, pady=5)

    def update_results(self, results):
        for param, value in results.items():
            self.result_labels[param].configure(state='normal')
            self.result_labels[param].delete(0, tk.END)
            self.result_labels[param].insert(0, f"{value:.2f}")
            self.result_labels[param].configure(state='readonly')

    def plot_pv_diagram(self, pv_points):
        self.ax.clear()
        v, p = zip(*pv_points)
        self.ax.plot(v, p, 'b-', label='Diesel Cycle')
        self.ax.set_xlabel('Volume (m^3)')
        self.ax.set_ylabel('Pressure (Pa)')
        self.ax.set_title('P-V Diagram')
        self.ax.grid(True)
        self.ax.legend()
        self.canvas.draw()

# Controller
class DieselController:
    def __init__(self, model, view):
        self.model = model
        self.view = view

    def calculate(self):
        self.model.r = float(self.view.r_entry.get())
        self.model.rc = float(self.view.rc_entry.get())
        self.model.T1 = float(self.view.t1_entry.get())
        self.model.P1 = float(self.view.p1_entry.get()) * 1e6

        results, pv_points = self.model.calculate_cycle()
        self.view.update_results(results)
        self.view.plot_pv_diagram(pv_points)

# Main application
if __name__ == "__main__":
    root = tk.Tk()
    model = DieselModel()
    view = DieselView(root, None)
    controller = DieselController(model, view)
    view.controller = controller
    root.mainloop()