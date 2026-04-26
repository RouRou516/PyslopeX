"""Plotting utilities for slope stability analysis.

Uses matplotlib for visualization of slope geometry, failure surfaces,
material layers, and water table.
"""

import math

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

from .utils import ground_elevation


# Default colors for materials
_MATERIAL_COLORS = [
    '#D2B48C',  # tan
    '#A0522D',  # sienna
    '#8B7355',  # burlywood
    '#BC8F8F',  # rosybrown
    '#C4A882',  # khaki-ish
    '#9C8A6E',  # olive-ish
    '#B8860B',  # darkgoldenrod
    '#DAA520',  # goldenrod
]


def _get_material_color(idx, material):
    """Get color for a material layer."""
    if material.color is not None:
        return material.color
    return _MATERIAL_COLORS[idx % len(_MATERIAL_COLORS)]


def plot_boundary(slope, show=True, figsize=(10, 6), title=None,
                  left_extend=None, right_extend=None):
    """Plot the slope boundary with material layers and water table.

    Args:
        slope: Slope object.
        show: Whether to call plt.show().
        figsize: Figure size.
        title: Plot title.
        left_extend: Distance to extend view left of toe.
        right_extend: Distance to extend view right of crest.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    H = slope.height
    L = slope.length

    if left_extend is None:
        left_extend = max(1.5 * H, 5.0)
    if right_extend is None:
        right_extend = max(1.5 * H, 5.0)

    # Draw material layers
    layers = slope.get_layer_boundaries()
    patches = []

    for idx, (y_top, y_bot, mat) in enumerate(layers):
        color = _get_material_color(idx, mat)
        x_max = L + right_extend

        # Layer extends across the entire model width
        # Calculate layer polygon considering slope face
        vertices = _layer_polygon(
            y_top, y_bot, -left_extend, x_max, L, H)
        if len(vertices) >= 3:
            poly = mpatches.Polygon(vertices, closed=True,
                                    facecolor=color, edgecolor='gray',
                                    alpha=0.6, linewidth=0.5)
            ax.add_patch(poly)
            patches.append((poly, mat))

    # Draw ground surface
    ground_pts = slope.get_ground_surface(left_extend, right_extend)
    gx = [p[0] for p in ground_pts]
    gy = [p[1] for p in ground_pts]
    ax.plot(gx, gy, 'k-', linewidth=2)

    # Draw water table
    if slope.water_table_depth is not None:
        y_wt = H - slope.water_table_depth
        ax.axhline(y=y_wt, color='dodgerblue', linestyle='--',
                    linewidth=1.5, label='Water Table')
        ax.fill_between([-left_extend, L + right_extend],
                        y_wt, [0, 0], alpha=0.15, color='dodgerblue')

    # Draw loads
    _draw_loads(ax, slope, L, H)

    # Draw material legend
    legend_handles = []
    for poly, mat in patches:
        legend_handles.append(
            mpatches.Patch(color=poly.get_facecolor(),
                           label=mat.name))
    if legend_handles:
        ax.legend(handles=legend_handles, loc='upper right', fontsize=8)

    # Labels and formatting
    ax.set_xlabel('Horizontal Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title(title or f'Slope: H={H}m, β={slope.angle:.1f}°')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-left_extend, L + right_extend)
    ax.set_ylim(-0.5, H + 1)

    plt.tight_layout()
    if show:
        plt.show()
    return fig


def _layer_polygon(y_top, y_bot, x_min, x_max, slope_length, slope_height):
    """Build polygon for a material layer clipped to the ground surface."""
    pts = []

    # Bottom-left corner
    pts.append((x_min, y_bot))
    # Bottom-right corner
    pts.append((x_max, y_bot))

    # Right side up to layer top
    pts.append((x_max, y_top))

    # Along the top of the layer, following ground surface where needed
    # Top-right going left
    # Check if the top of the layer intersects the slope face
    # Slope face: y = x * H / L for 0 <= x <= L
    if slope_length > 0:
        slope_rate = slope_height / slope_length
    else:
        slope_rate = float('inf')

    # Intersection of y_top with slope face
    if y_top <= slope_height and y_top >= 0:
        x_on_slope_at_ytop = y_top / slope_rate if slope_rate > 0 else 0
    else:
        x_on_slope_at_ytop = None

    if y_bot <= slope_height and y_bot >= 0:
        x_on_slope_at_ybot = y_bot / slope_rate if slope_rate > 0 else 0
    else:
        x_on_slope_at_ybot = None

    # Build top edge (left to right)
    top_pts = []

    # From x_min going right along y_top
    # If we cross the slope face, follow it
    if x_on_slope_at_ytop is not None and x_on_slope_at_ytop > x_min:
        # Ground is above y_top for x < x_on_slope_at_ytop on the slope face
        # and below or at y_top for x > x_on_slope_at_ytop
        top_pts.append((x_min, y_top))
    else:
        top_pts.append((x_min, y_top))

    # Left side down
    pts.append((x_min, y_top))

    # Reorder: build the full polygon properly
    # Simplified approach: use a rectangular clip with ground surface
    result = []
    n_pts = 100
    for i in range(n_pts + 1):
        x = x_min + (x_max - x_min) * i / n_pts
        y_ground = ground_elevation(x, slope_length, slope_height)
        y_layer_top = min(y_top, y_ground)
        y_layer_bot = max(y_bot, 0)
        if y_layer_top > y_layer_bot:
            result.append((x, y_layer_top))

    for i in range(n_pts, -1, -1):
        x = x_min + (x_max - x_min) * i / n_pts
        y_ground = ground_elevation(x, slope_length, slope_height)
        y_layer_top = min(y_top, y_ground)
        y_layer_bot = max(y_bot, 0)
        if y_layer_top > y_layer_bot:
            result.append((x, y_layer_bot))

    return result


def plot_planar_result(slope, result, show=True, figsize=(10, 6)):
    """Plot slope with the critical planar failure surface.

    Args:
        slope: Slope object.
        result: PlanarResult from analysis.
        show: Whether to call plt.show().

    Returns:
        matplotlib Figure.
    """
    fig = plot_boundary(slope, show=False, figsize=figsize,
                        title=f'Planar Sliding: FOS = {result.fos:.3f}, '
                              f'θ = {result.critical_angle:.1f}°')

    ax = fig.axes[0]

    # Draw critical failure plane
    if result.failure_plane:
        fx = [p[0] for p in result.failure_plane]
        fy = [p[1] for p in result.failure_plane]
        ax.plot(fx, fy, 'r-', linewidth=2.5, label='Failure Plane')

    # Draw sliding wedge
    if result.sliding_wedge:
        wedge = mpatches.Polygon(
            result.sliding_wedge, closed=True,
            facecolor='red', alpha=0.15, edgecolor='red',
            linewidth=1.5, label='Sliding Wedge')
        ax.add_patch(wedge)

    # FOS annotation
    if result.failure_plane:
        x_mid = (result.failure_plane[0][0]
                 + result.failure_plane[1][0]) / 2
        y_mid = (result.failure_plane[0][1]
                 + result.failure_plane[1][1]) / 2
        ax.annotate(f'FOS = {result.fos:.3f}',
                    xy=(x_mid, y_mid),
                    fontsize=12, fontweight='bold', color='red',
                    bbox=dict(boxstyle='round,pad=0.3',
                              facecolor='white', alpha=0.8))

    ax.legend(loc='upper right', fontsize=8)
    plt.tight_layout()
    if show:
        plt.show()
    return fig


def plot_polyline_result(slope, result, show=True, figsize=(10, 6)):
    """Plot slope with the polyline failure surface.

    Args:
        slope: Slope object.
        result: PolylineResult from analysis.
        show: Whether to call plt.show().

    Returns:
        matplotlib Figure.
    """
    fig = plot_boundary(slope, show=False, figsize=figsize,
                        title=f'Polyline Sliding ({result.method}): '
                              f'FOS = {result.fos:.3f}')

    ax = fig.axes[0]

    # Draw failure surface
    if result.failure_surface:
        fx = [p[0] for p in result.failure_surface]
        fy = [p[1] for p in result.failure_surface]
        ax.plot(fx, fy, 'r-', linewidth=2.5, label='Failure Surface')

        # Draw failure surface vertices
        ax.plot(fx, fy, 'ro', markersize=5)

    # Draw thrust arrows at slice boundaries
    if (result.thrust_distribution
            and result.failure_surface
            and len(result.thrust_distribution) > 1):
        max_thrust = max(abs(t) for t in result.thrust_distribution)
        if max_thrust > 1e-6:
            for i, thrust in enumerate(result.thrust_distribution):
                if i < len(result.failure_surface):
                    px, py = result.failure_surface[i]
                    scale = thrust / max_thrust * slope.height * 0.3
                    if abs(scale) > 0.1:
                        ax.annotate(
                            '', xy=(px, py + scale * 0.5),
                            xytext=(px, py),
                            arrowprops=dict(
                                arrowstyle='->', color='blue',
                                lw=1.5))

    # Draw slice boundaries
    if result.failure_surface:
        for pt in result.failure_surface:
            y_ground = ground_elevation(pt[0], slope.length, slope.height)
            ax.plot([pt[0], pt[0]], [pt[1], y_ground],
                    'b--', linewidth=0.8, alpha=0.5)

    # FOS annotation
    if result.failure_surface and len(result.failure_surface) >= 2:
        mid_idx = len(result.failure_surface) // 2
        px, py = result.failure_surface[mid_idx]
        y_ground = ground_elevation(px, slope.length, slope.height)
        ax.annotate(f'FOS = {result.fos:.3f}',
                    xy=(px, (py + y_ground) / 2),
                    fontsize=12, fontweight='bold', color='red',
                    bbox=dict(boxstyle='round,pad=0.3',
                              facecolor='white', alpha=0.8))

    ax.legend(loc='upper right', fontsize=8)
    plt.tight_layout()
    if show:
        plt.show()
    return fig


def _draw_loads(ax, slope, x_crest, H):
    """Draw load indicators on the plot."""
    for udl in slope.udls:
        x_start = x_crest + udl.offset
        x_end = x_start + (udl.length if udl.length else H)
        color = udl.color or 'green'
        ax.annotate('', xy=(x_start, H), xytext=(x_start, H + 0.8),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        ax.annotate('', xy=(x_end, H), xytext=(x_end, H + 0.8),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        ax.plot([x_start, x_end], [H + 0.8, H + 0.8],
                color=color, linewidth=2)
        ax.text((x_start + x_end) / 2, H + 1.1,
                f'{udl.magnitude} kPa', ha='center', fontsize=8,
                color=color)

    for ll in slope.line_loads:
        x_pos = x_crest + ll.offset
        color = ll.color or 'purple'
        ax.annotate('', xy=(x_pos, H), xytext=(x_pos, H + 1.2),
                    arrowprops=dict(arrowstyle='->', color=color, lw=2))
        ax.text(x_pos, H + 1.4, f'{ll.magnitude} kN/m',
                ha='center', fontsize=8, color=color)
