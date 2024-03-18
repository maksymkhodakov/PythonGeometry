import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d


def point_inside_polygon(point, polygon):
    """
    Check if a point is inside a polygon.

    Parameters:
        point (tuple): Coordinates of the point (x, y).
        polygon (list of tuples): List of vertices of the polygon [(x1, y1), (x2, y2), ..., (xn, yn)].

    Returns:
        bool: True if the point is inside the polygon, False otherwise.
    """
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


def distance(p1, p2):
    """
    Compute the Euclidean distance between two points.

    Parameters:
        p1 (tuple): First point coordinates (x1, y1).
        p2 (tuple): Second point coordinates (x2, y2).

    Returns:
        float: Euclidean distance between the points.
    """
    return np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


def generate_star_polygon(N, inner_radius, outer_radius):
    angles = np.linspace(0, 2 * np.pi, 2 * N, endpoint=False)
    vertices = []
    for angle in angles:
        vertices.append((outer_radius * np.cos(angle), outer_radius * np.sin(angle)))
        vertices.append((inner_radius * np.cos(angle + np.pi / N), inner_radius * np.sin(angle + np.pi / N)))
    return np.array(vertices)


# User input for number of points
N = int(input("Enter the number of points for the star polygon: "))

# Generate random star-shaped polygon vertices
outer_radius = 10
inner_radius = 5
polygon_vertices = generate_star_polygon(N, inner_radius, outer_radius)

# Generate Voronoi diagram
vor = Voronoi(polygon_vertices)

# Find the Voronoi vertices within the polygon
inside_vertices = []
for vertex in vor.vertices:
    if point_inside_polygon(vertex, polygon_vertices):
        inside_vertices.append(vertex)

# Find the maximum inscribed circle
max_radius = 0
max_center = None
for vertex in inside_vertices:
    radius = min(distance(vertex, v) for v in polygon_vertices)
    if radius > max_radius:
        max_radius = radius
        max_center = vertex

# Plot the polygon, Voronoi diagram, and maximum inscribed circle
plt.figure(figsize=(8, 6))
plt.plot(polygon_vertices[:, 0], polygon_vertices[:, 1], 'bo-')  # Plot polygon
voronoi_plot_2d(vor, show_vertices=False, line_colors='orange', line_width=2)  # Plot Voronoi diagram
circle = plt.Circle(max_center, max_radius, color='red', alpha=0.5)  # Maximum inscribed circle
plt.gca().add_patch(circle)
plt.xlim(min(polygon_vertices[:, 0]) - 1, max(polygon_vertices[:, 0]) + 1)
plt.ylim(min(polygon_vertices[:, 1]) - 1, max(polygon_vertices[:, 1]) + 1)
plt.title("Star-shaped Polygon, Voronoi Diagram, and Maximum Inscribed Circle")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()
