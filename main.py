import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from tkinter import simpledialog, Toplevel
from matplotlib.figure import Figure
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def point_inside_polygon(point, polygon):
    x, y = point
    inside = False
    n = len(polygon)
    p1x, p1y = polygon[0]
    for i in range(n + 1):
        p2x, p2y = polygon[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside


def distance_point_to_segment(point, segment):
    """Calculate the minimum distance from a point to a line segment."""
    p1, p2 = np.array(segment[0]), np.array(segment[1])
    p = np.array(point)
    line_vec = p2 - p1
    point_vec = p - p1
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len
    p1_to_p = np.linalg.norm(point_vec)
    projection = point_vec.dot(line_unitvec)
    if projection < 0:
        return p1_to_p
    if projection > line_len:
        return np.linalg.norm(p - p2)
    projected_point = p1 + line_unitvec * projection
    return np.linalg.norm(p - projected_point)


def generate_star_polygon(N, inner_radius, outer_radius):
    angles = np.linspace(0, 2 * np.pi, 2 * N, endpoint=False)
    vertices = []
    for angle in angles:
        r = outer_radius if angle % (2 * np.pi / N) == 0 else inner_radius
        vertices.append((r * np.cos(angle), r * np.sin(angle)))
    return np.array(vertices)


def plot_polygon_and_circle(N, outer_radius=10, inner_radius=5):
    # Generate star polygon and calculate Voronoi diagram
    polygon_vertices = generate_star_polygon(N, inner_radius, outer_radius)
    vor = Voronoi(polygon_vertices)
    inside_vertices = [v for v in vor.vertices if point_inside_polygon(v, polygon_vertices)]

    # Find the maximum inscribed circle
    max_area = 0
    max_center = None
    for vertex in inside_vertices:
        distances = [distance_point_to_segment(vertex, (polygon_vertices[i], polygon_vertices[(i + 1) % N])) for i in
                     range(N)]
        radius = min(distances)
        area = np.pi * radius ** 2
        if area > max_area:
            max_area = area
            max_center = vertex

    # Create a new Tkinter window
    new_window = Toplevel(ROOT)
    new_window.title("Polygon Plot")

    # Create a Figure for Matplotlib
    fig = Figure(figsize=(6, 4))
    ax = fig.add_subplot(111)

    # Plotting
    ax.plot(polygon_vertices[:, 0], polygon_vertices[:, 1], 'bo-')
    voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors='orange', line_width=2)
    if max_center is not None:
        circle = plt.Circle(max_center, np.sqrt(max_area / np.pi), color='red', alpha=0.5)
        ax.add_patch(circle)
    ax.set_xlim(min(polygon_vertices[:, 0]) - 1, max(polygon_vertices[:, 0]) + 1)
    ax.set_ylim(min(polygon_vertices[:, 1]) - 1, max(polygon_vertices[:, 1]) + 1)
    ax.set_title("Star-shaped Polygon, Voronoi Diagram, and Maximum Inscribed Circle by Area")
    ax.grid(True)

    # Embedding the Matplotlib figure in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=new_window)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


# Creating the UI
ROOT = tk.Tk()
ROOT.withdraw()
# The input dialog
USER_INP = simpledialog.askinteger(title="Star Polygon Generator",
                                   prompt="Enter the number of points for the star polygon:")

# Call the function with the user input
if USER_INP:
    plot_polygon_and_circle(USER_INP)
    ROOT.mainloop()
