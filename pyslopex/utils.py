"""Geometry and math utilities for slope stability analysis."""

import math


def line_intersection(p1, p2, p3, p4):
    """Find intersection of two lines defined by point pairs.

    Args:
        p1, p2: Points defining the first line.
        p3, p4: Points defining the second line.

    Returns:
        (x, y) tuple or None if lines are parallel.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1e-12:
        return None

    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom

    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)
    return (x, y)


def segment_intersection(p1, p2, p3, p4):
    """Find intersection of two line segments.

    Returns (x, y) if segments intersect, None otherwise.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1e-12:
        return None

    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom

    if 0 <= t <= 1 and 0 <= u <= 1:
        x = x1 + t * (x2 - x1)
        y = y1 + t * (y2 - y1)
        return (x, y)
    return None


def polygon_area(vertices):
    """Calculate area of a polygon using the shoelace formula.

    Args:
        vertices: List of (x, y) tuples.

    Returns:
        Absolute area of the polygon.
    """
    n = len(vertices)
    if n < 3:
        return 0.0
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
    return abs(area) / 2.0


def polygon_centroid(vertices):
    """Calculate the centroid of a polygon.

    Returns:
        (cx, cy) tuple.
    """
    n = len(vertices)
    if n == 0:
        return (0, 0)
    if n < 3:
        cx = sum(v[0] for v in vertices) / n
        cy = sum(v[1] for v in vertices) / n
        return (cx, cy)

    a = 0.0
    cx = 0.0
    cy = 0.0
    for i in range(n):
        j = (i + 1) % n
        cross = vertices[i][0] * vertices[j][1] - vertices[j][0] * vertices[i][1]
        a += cross
        cx += (vertices[i][0] + vertices[j][0]) * cross
        cy += (vertices[i][1] + vertices[j][1]) * cross
    a /= 2.0
    if abs(a) < 1e-12:
        cx = sum(v[0] for v in vertices) / n
        cy = sum(v[1] for v in vertices) / n
        return (cx, cy)
    cx /= (6.0 * a)
    cy /= (6.0 * a)
    return (cx, cy)


def clip_polygon_by_horizontal_line(vertices, y_level, keep_above=True):
    """Clip a polygon by a horizontal line using Sutherland-Hodgman algorithm.

    Args:
        vertices: Polygon vertices as list of (x, y).
        y_level: Y-coordinate of the clipping line.
        keep_above: If True, keep the portion above y_level.

    Returns:
        Clipped polygon vertices.
    """
    if not vertices:
        return []

    def inside(p):
        return p[1] >= y_level if keep_above else p[1] <= y_level

    def intersect_edge(p1, p2):
        if abs(p2[1] - p1[1]) < 1e-12:
            return p1
        t = (y_level - p1[1]) / (p2[1] - p1[1])
        x = p1[0] + t * (p2[0] - p1[0])
        return (x, y_level)

    output = []
    n = len(vertices)
    for i in range(n):
        current = vertices[i]
        next_pt = vertices[(i + 1) % n]

        curr_inside = inside(current)
        next_inside = inside(next_pt)

        if curr_inside:
            output.append(current)
            if not next_inside:
                output.append(intersect_edge(current, next_pt))
        elif next_inside:
            output.append(intersect_edge(current, next_pt))

    return output


def ground_elevation(x, slope_length, slope_height):
    """Get ground surface elevation at a given x-coordinate.

    The ground surface consists of:
    - Flat base (y=0) for x < 0
    - Slope face for 0 <= x <= slope_length
    - Flat top (y=H) for x > slope_length

    Args:
        x: Horizontal coordinate.
        slope_length: Horizontal projection of the slope face (m).
        slope_height: Slope height (m).

    Returns:
        Ground elevation y at position x.
    """
    if x <= 0:
        return 0.0
    elif x <= slope_length:
        return x * slope_height / slope_length
    else:
        return float(slope_height)
