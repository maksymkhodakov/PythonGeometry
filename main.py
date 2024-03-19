import tkinter as tk
from tkinter import Toplevel, Label, Entry, Button
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib.figure import Figure
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


def create_plot_frame(N, plot_window=None):
    # Ensure N is an integer
    N = int(N)

    # Generate star polygon and calculate Voronoi diagram
    outer_radius, inner_radius = 10, 5
    polygon_vertices = generate_star_polygon(N, inner_radius, outer_radius)
    vor = Voronoi(polygon_vertices)

    # Create a Figure for Matplotlib
    fig = Figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # Plotting
    ax.plot(polygon_vertices[:, 0], polygon_vertices[:, 1], 'bo-', markersize=5)  # Polygon points
    # Close the polygon path
    ax.plot([polygon_vertices[-1, 0], polygon_vertices[0, 0]], [polygon_vertices[-1, 1], polygon_vertices[0, 1]], 'bo-', markersize=5)

    # Find the Voronoi vertices inside the polygon and calculate the maximum inscribed circle
    inside_vertices = np.array([v for v in vor.vertices if point_inside_polygon(v, polygon_vertices)])

    # Find the point inside the polygon that maximizes the distance to the closest edge
    # This will be an approximation to the problem and works well when the Voronoi vertex is close to the center of the polygon
    max_radius = 0
    max_vertex = None
    for vertex in inside_vertices:
        radius = np.min([distance_point_to_segment(vertex, (polygon_vertices[i], polygon_vertices[(i + 1) % N])) for i in range(N)])
        if radius > max_radius:
            max_radius = radius
            max_vertex = vertex

    # Draw the circle inside the polygon
    if max_vertex is not None:
        circle = plt.Circle(max_vertex, max_radius, color='red', fill=False, linestyle='--', linewidth=1.5)
        ax.add_patch(circle)
        # Annotate the circle's radius
        ax.annotate(f'Radius: {max_radius:.2f}', xy=(0, 1), xycoords='axes fraction', fontsize=9,
                    xytext=(5, -5), textcoords='offset points', ha='left', va='top')

    # Set plot limits
    ax.set_xlim([polygon_vertices[:, 0].min() - 1, polygon_vertices[:, 0].max() + 1])
    ax.set_ylim([polygon_vertices[:, 1].min() - 1, polygon_vertices[:, 1].max() + 1])

    ax.set_title("Star-shaped Polygon and Maximum Inscribed Circle")

    # Embedding the Matplotlib figure in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=plot_window)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


# Creating the UI
def setup_ui():
    root = tk.Tk()
    root.title("Star Polygon Generator")

    label = Label(root, text="Enter the number of points for the star polygon:")
    label.pack()

    n_entry = Entry(root)
    n_entry.pack()
    n_entry.insert(0, "5")  # Default starting value

    plot_window = None  # This will hold the Toplevel window for plotting

    def draw_plot():
        nonlocal plot_window
        N = n_entry.get()
        if plot_window and plot_window.winfo_exists():
            plot_window.destroy()  # Close existing plot window if open
        plot_window = Toplevel()
        plot_window.title("Polygon Plot")
        create_plot_frame(N, plot_window)

    draw_button = Button(root, text="Draw Plot", command=draw_plot)
    draw_button.pack()

    root.mainloop()


if __name__ == "__main__":
    setup_ui()
