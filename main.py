import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter

def generate_sierpinski(level):
    v0 = np.array([0.0, 0.0])
    v1 = np.array([1.0, 0.0])
    v2 = np.array([0.5, np.sqrt(3)/2])
    base_triangle = [v0, v1, v2]
    def F(i, pts):
        if i == 0:
            return [0.5 * p for p in pts]
        elif i == 1:
            return [0.5 * p + [0.5, 0.0] for p in pts]
        elif i == 2:
            return [0.5 * p + [0.25, np.sqrt(3)/4] for p in pts]
    triangles = [base_triangle]
    for _ in range(level):
        new_triangles = []
        for tri in triangles:
            for i in range(3):
                new_triangles.append(F(i, tri))
        triangles = new_triangles
    tol = 1e-8
    vertex_map = {}
    vertices = []
    def register_vertex(p):
        key = tuple(np.round(p / tol).astype(int))
        if key not in vertex_map:
            vertex_map[key] = len(vertices)
            vertices.append(p)
        return vertex_map[key]
    edges = set()
    for tri in triangles:
        idx = [register_vertex(p) for p in tri]
        for i in range(3):
            a, b = idx[i], idx[(i+1)%3]
            edges.add(tuple(sorted((a, b))))
    return np.array(vertices), list(edges)

class WaveApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Sierpinski Gasket Wave Simulation")
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        self.frame2d = ttk.Frame(self.notebook)
        self.frame3d = ttk.Frame(self.notebook)
        self.frame_energy = ttk.Frame(self.notebook)
        self.notebook.add(self.frame2d, text="2D Animation")
        self.notebook.add(self.frame3d, text="3D Animation")
        self.notebook.add(self.frame_energy, text="Energy vs Time")
        self.build_tab_2d()
        self.build_tab_3d()
        self.build_tab_energy()
        self.frame2d.rowconfigure(0, weight=1)
        self.frame2d.columnconfigure(1, weight=1)
        self.frame3d.rowconfigure(0, weight=1)
        self.frame3d.columnconfigure(0, weight=1)
        self.frame_energy.rowconfigure(0, weight=1)
        self.frame_energy.columnconfigure(0, weight=1)

    def build_tab_2d(self):
        controls = tk.Frame(self.frame2d)
        controls.grid(row=0, column=0, sticky='ns')

        # Defaults
        param_list = [
            ("Gasket Level", "level_var", 3),
            ("TFINAL", "tfinal_var", 2.0),
            ("Frames", "frames_var", 200),
            ("Displacement Vertex", "disp_vertex_var", 0),
            ("Displacement Amplitude", "disp_amp_var", 1.0),
            ("Displacement Sigma", "disp_sigma_var", 0.1),
            ("Velocity Vertex", "vel_vertex_var", 0),
            ("Velocity Amplitude", "vel_amp_var", 0.0),
            ("Velocity Sigma", "vel_sigma_var", 0.1)
        ]
        self.vars = {}
        for i, (label, varname, default) in enumerate(param_list):
            tk.Label(controls, text=label).grid(row=i, column=0, sticky='w')
            vtype = tk.DoubleVar if isinstance(default, float) else tk.IntVar
            self.vars[varname] = vtype(value=default)
            tk.Entry(controls, textvariable=self.vars[varname], width=6).grid(row=i, column=1)

        tk.Button(controls, text="Start Simulation", command=self.start_simulation).grid(row=10, column=0, columnspan=2, pady=10)
        tk.Button(controls, text="Reset", command=self.reset).grid(row=11, column=0, columnspan=2, pady=5)
        tk.Button(controls, text="Export 2D GIF", command=self.export_gif).grid(row=12, column=0, columnspan=2, pady=5)
        tk.Button(controls, text="Export 3D GIF", command=self.export_gif_3d).grid(row=13, column=0, columnspan=2, pady=5)

        self.fig2d, self.ax2d = plt.subplots(figsize=(5, 5))
        self.canvas2d = FigureCanvasTkAgg(self.fig2d, master=self.frame2d)
        self.canvas2d.get_tk_widget().grid(row=0, column=1, sticky='nsew')

    def build_tab_3d(self):
        self.fig3d = plt.figure(figsize=(5,5))
        self.ax3d = self.fig3d.add_subplot(111, projection='3d')
        self.canvas3d = FigureCanvasTkAgg(self.fig3d, master=self.frame3d)
        self.canvas3d.get_tk_widget().grid(row=0, column=0, sticky='nsew')

    def build_tab_energy(self):
        self.fig_energy, self.ax_energy = plt.subplots(figsize=(5, 3))
        self.canvas_energy = FigureCanvasTkAgg(self.fig_energy, master=self.frame_energy)
        self.canvas_energy.get_tk_widget().grid(row=0, column=0, sticky='nsew')

    def start_simulation(self):
        p = self.vars
        level = p["level_var"].get()
        tfinal = p["tfinal_var"].get()
        frames = p["frames_var"].get()
        disp_vertex = p["disp_vertex_var"].get()
        disp_amp = p["disp_amp_var"].get()
        disp_sigma = p["disp_sigma_var"].get()
        vel_vertex = p["vel_vertex_var"].get()
        vel_amp = p["vel_amp_var"].get()
        vel_sigma = p["vel_sigma_var"].get()

        positions, edges = generate_sierpinski(level)
        n = len(positions)
        if disp_vertex >= n or vel_vertex >= n:
            tk.messagebox.showerror("Error", f"Vertex index out of range (max index {n-1})")
            return

        A = np.zeros((n, n))
        for i, j in edges:
            A[i, j] = A[j, i] = 1
        D = np.diag(A.sum(axis=1))
        L = D - A
        r = 3/5
        L = (1/r**level) * L

        dpos = np.linalg.norm(positions - positions[disp_vertex], axis=1)
        f = disp_amp * np.exp(-dpos**2 / (2 * disp_sigma**2))

        vpos = np.linalg.norm(positions - positions[vel_vertex], axis=1)
        v = vel_amp * np.exp(-vpos**2 / (2 * vel_sigma**2))

        eigvals, eigvecs = np.linalg.eigh(L)
        f_hat = eigvecs.T @ f
        v_hat = eigvecs.T @ v
        times = np.linspace(0, tfinal, frames)
        sqrt_lambda = np.sqrt(eigvals)
        with np.errstate(divide='ignore', invalid='ignore'):
            v_term = np.where(sqrt_lambda != 0, v_hat / sqrt_lambda, 0)
        u_over_time = np.array([
            eigvecs @ (f_hat * np.cos(sqrt_lambda * t) + v_term * np.sin(sqrt_lambda * t))
            for t in times
        ])
        u_over_time[:, eigvals==0] += v_hat[eigvals==0] * times[:, None]

        u_t_over_time = np.array([
            eigvecs @ (-f_hat * sqrt_lambda * np.sin(sqrt_lambda * t) + v_hat * np.cos(sqrt_lambda * t))
            for t in times
        ])

        self.positions = positions
        self.edges = edges
        self.u_over_time = u_over_time
        self.u_t_over_time = u_t_over_time
        self.times = times
        self.L = L

        kinetic = np.sum(u_t_over_time ** 2, axis=1) / 2
        potential = np.array([u @ (L @ u) / 2 for u in u_over_time])
        total_energy = kinetic + potential

        self.ax_energy.clear()
        self.ax_energy.plot(times, total_energy, label="Total Energy")
        self.ax_energy.set_xlabel("Time")
        self.ax_energy.set_ylabel("Energy")
        self.ax_energy.legend()
        self.canvas_energy.draw()

        self.start_2d_animation()
        self.start_3d_animation()

    def start_2d_animation(self):
        x, y = self.positions[:, 0], self.positions[:, 1]
        self.ax2d.clear()
        for i, j in self.edges:
            self.ax2d.plot(x[[i, j]], y[[i, j]], color='gray', lw=0.8)
        sc = self.ax2d.scatter(x, y, c=self.u_over_time[0], cmap='coolwarm',
                               s=1000//self.vars["level_var"].get(), edgecolors='black', vmin=-1, vmax=1)
        self.ax2d.set_aspect('equal')
        def update(frame):
            sc.set_array(self.u_over_time[frame])
            self.ax2d.set_title(f"t = {self.times[frame]:.2f}")
            return sc,
        self.anim2d = animation.FuncAnimation(self.fig2d, update, frames=len(self.times), blit=False, interval=50)
        self.canvas2d.draw()

    def start_3d_animation(self):
        x, y = self.positions[:, 0], self.positions[:, 1]
        z = self.u_over_time[0]
        self.ax3d.clear()
        self.scatter3d = self.ax3d.scatter(x, y, z, c=z, cmap='coolwarm', s=50, edgecolor='k')
        self.lines3d = [self.ax3d.plot(x[[i,j]], y[[i,j]], z[[i,j]], color='gray', lw=0.5)[0] for i,j in self.edges]
        def update_3d(frame):
            z = self.u_over_time[frame]
            self.scatter3d.remove()
            self.scatter3d = self.ax3d.scatter(x, y, z, c=z, cmap='coolwarm', s=50, edgecolor='k')
            for k, (i,j) in enumerate(self.edges):
                self.lines3d[k].set_data(x[[i,j]], y[[i,j]])
                self.lines3d[k].set_3d_properties(z[[i,j]])
            self.ax3d.set_title(f"t = {self.times[frame]:.2f}")
            return [self.scatter3d] + self.lines3d
        self.anim3d = animation.FuncAnimation(self.fig3d, update_3d, frames=len(self.times), blit=False, interval=50)
        self.canvas3d.draw()

    def reset(self):
        if hasattr(self, 'anim2d') and self.anim2d:
            self.anim2d.event_source.stop()
            self.anim2d = None
        if hasattr(self, 'anim3d') and self.anim3d:
            self.anim3d.event_source.stop()
            self.anim3d = None
        self.ax2d.clear()
        self.canvas2d.draw()
        self.ax3d.clear()
        self.canvas3d.draw()
        self.ax_energy.clear()
        self.canvas_energy.draw()

    def export_gif(self):
        if self.anim2d is None:
            messagebox.showinfo("Info", "Run simulation first!")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".gif", filetypes=[("GIF files","*.gif")])
        if file_path:
            writer = PillowWriter(fps=20)
            self.anim2d.save(file_path, writer=writer)
            messagebox.showinfo("Done", f"Saved 2D GIF to:\n{file_path}")

    def export_gif_3d(self):
        if self.anim3d is None:
            messagebox.showinfo("Info", "Run simulation first!")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".gif", filetypes=[("GIF files","*.gif")])
        if file_path:
            writer = PillowWriter(fps=20)
            self.anim3d.save(file_path, writer=writer)
            messagebox.showinfo("Done", f"Saved 3D GIF to:\n{file_path}")

if __name__ == "__main__":
    root = tk.Tk()
    app = WaveApp(root)
    root.mainloop()
