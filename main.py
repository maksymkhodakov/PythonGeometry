import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import shapely.geometry as geometry
import shapely.ops as ops
import tkinter as tk
from tkinter import simpledialog
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def find_largest_inscribed_circle(points, hull):
    poly = geometry.Polygon(points[hull.vertices])
    center = poly.centroid.coords[0]
    min_dist_to_edge = poly.boundary.distance(geometry.Point(center))
    return center, min_dist_to_edge


class ConvexHullApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Convex Hull and Largest Inscribed Circle")
        self.geometry("1920x1080")
        self.create_widgets()

    def create_widgets(self):
        self.points_entry = tk.Entry(self)
        self.points_entry.pack(pady=10)
        self.plot_button = tk.Button(self, text="Generate and Plot", command=self.generate_and_plot)
        self.plot_button.pack(pady=10)
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def generate_and_plot(self):
        n = int(self.points_entry.get())
        points = np.random.rand(n, 2) * 100
        hull = ConvexHull(points)
        center, radius = find_largest_inscribed_circle(points, hull)

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.plot(points[:, 0], points[:, 1], 'o', markersize=3)
        for simplex in hull.simplices:
            ax.plot(points[simplex, 0], points[simplex, 1], 'k-')
        circle = plt.Circle(center, radius, color='r', fill=False)
        ax.add_artist(circle)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(f'Convex Hull and Inscribed Circle (Radius: {radius:.2f})')
        self.canvas.draw()


if __name__ == '__main__':
    app = ConvexHullApp()
    app.mainloop()
