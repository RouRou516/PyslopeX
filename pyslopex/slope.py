"""Slope model definition for slope stability analysis.

Coordinate system:
    - x: horizontal, positive to the right
    - y: vertical, positive upward
    - Toe at origin (0, 0)
    - Crest at (slope_length, slope_height)
    - Slope face rises from toe (left) to crest (right)
"""

import math

from .material import Material
from .loads import Udl, LineLoad
from .utils import ground_elevation


class Slope:
    """Represents a 2D slope model for stability analysis.

    Args:
        height: Slope height (m).
        angle: Slope angle from horizontal (degrees). Use None if providing length.
        length: Horizontal projection of slope face (m). Use None if providing angle.

    Raises:
        ValueError: If neither or both of angle/length are provided.
    """

    def __init__(self, height, angle=None, length=None):
        if angle is None and length is None:
            raise ValueError("Either angle or length must be provided")
        if angle is not None and length is not None:
            raise ValueError("Only one of angle or length should be provided")

        self.height = height

        if angle is not None:
            self.angle = float(angle)
            self.length = height / math.tan(math.radians(angle))
        else:
            self.length = float(length)
            self.angle = math.degrees(math.atan2(height, length))

        self.materials = []
        self.water_table_depth = None
        self.udls = []
        self.line_loads = []
        self._analysis_left_limit = None
        self._analysis_right_limit = None

    @property
    def toe(self):
        """Toe (bottom) of the slope."""
        return (0.0, 0.0)

    @property
    def crest(self):
        """Crest (top) of the slope."""
        return (self.length, float(self.height))

    def get_top_coordinates(self):
        """Return the crest coordinates."""
        return self.crest

    def get_bottom_coordinates(self):
        """Return the toe coordinates."""
        return self.toe

    def set_materials(self, *materials):
        """Assign material layers to the slope.

        Materials are sorted by depth_to_bottom automatically.
        The last material extends infinitely downward.

        Args:
            *materials: Material objects.
        """
        mats = list(materials)
        for i, m in enumerate(mats):
            for j, m2 in enumerate(mats):
                if i != j and m.depth_to_bottom == m2.depth_to_bottom:
                    raise ValueError(
                        f"Duplicate depth_to_bottom: {m.depth_to_bottom}")
        self.materials = sorted(mats, key=lambda m: m.depth_to_bottom)

    def set_water_table(self, depth=None):
        """Set water table depth from the top of the slope.

        Args:
            depth: Depth from top surface to water table (m).
                   None removes the water table.
        """
        self.water_table_depth = depth

    def set_udls(self, *udls):
        """Assign uniform distributed loads.

        Args:
            *udls: Udl objects.
        """
        self.udls = list(udls)

    def remove_udls(self, *udls):
        """Remove specific UDLs from the slope."""
        for u in udls:
            if u in self.udls:
                self.udls.remove(u)

    def set_line_loads(self, *line_loads):
        """Assign line loads.

        Args:
            *line_loads: LineLoad objects.
        """
        self.line_loads = list(line_loads)

    def set_analysis_limits(self, left, right):
        """Set left and right limits for failure surface search.

        Args:
            left: Left (minimum x) limit.
            right: Right (maximum x) limit.
        """
        self._analysis_left_limit = left
        self._analysis_right_limit = right

    def get_ground_elevation(self, x):
        """Get ground surface elevation at x.

        Args:
            x: Horizontal coordinate.

        Returns:
            Elevation y.
        """
        return ground_elevation(x, self.length, self.height)

    def get_ground_surface(self, left_extend=None, right_extend=None):
        """Return ground surface as a list of (x, y) points.

        Args:
            left_extend: Distance to extend left of the toe.
            right_extend: Distance to extend right of the crest.

        Returns:
            List of (x, y) tuples.
        """
        if left_extend is None:
            left_extend = max(1.5 * self.height, 5.0)
        if right_extend is None:
            right_extend = max(1.5 * self.height, 5.0)

        return [
            (-left_extend, 0.0),
            (0.0, 0.0),
            (self.length, float(self.height)),
            (self.length + right_extend, float(self.height)),
        ]

    def get_layer_boundaries(self):
        """Return layer boundaries as list of (y_top, y_bottom, material).

        The last material extends below all defined layers.
        """
        H = self.height
        layers = []
        prev_depth = 0.0
        for mat in self.materials:
            y_top = H - prev_depth
            y_bot = H - mat.depth_to_bottom
            layers.append((y_top, y_bot, mat))
            prev_depth = mat.depth_to_bottom
        return layers

    def get_material_at_y(self, y):
        """Get the material at a given elevation y.

        Args:
            y: Elevation.

        Returns:
            Material object at that elevation.
        """
        for y_top, y_bot, mat in self.get_layer_boundaries():
            if y_bot - 1e-9 <= y <= y_top + 1e-9:
                return mat
        if self.materials:
            return self.materials[-1]
        return None

    def analyse_planar(self, **kwargs):
        """Run planar sliding analysis.

        Returns:
            PlanarResult with factor of safety and critical failure plane.
        """
        from .planar import PlanarAnalysis
        analysis = PlanarAnalysis(self)
        return analysis.analyse(**kwargs)

    def analyse_polyline(self, failure_surface_points=None, **kwargs):
        """Run polyline sliding analysis (transfer coefficient method).

        Args:
            failure_surface_points: List of (x, y) points defining the failure
                surface, ordered from entry (upper right) to exit (lower left/toe).

        Returns:
            PolylineResult with factor of safety and thrust distribution.
        """
        from .polyline import PolylineAnalysis
        analysis = PolylineAnalysis(self)
        return analysis.analyse(
            failure_surface_points=failure_surface_points, **kwargs)

    def plot_boundary(self, **kwargs):
        """Plot slope boundary with materials and water table."""
        from .plotting import plot_boundary
        return plot_boundary(self, **kwargs)

    def plot_critical_planar(self, result=None, **kwargs):
        """Plot slope with critical planar failure surface."""
        from .plotting import plot_planar_result
        if result is None:
            result = self.analyse_planar()
        return plot_planar_result(self, result, **kwargs)

    def plot_polyline(self, failure_surface_points=None, result=None, **kwargs):
        """Plot slope with polyline failure surface."""
        from .plotting import plot_polyline_result
        if result is None:
            result = self.analyse_polyline(
                failure_surface_points=failure_surface_points)
        return plot_polyline_result(self, result, **kwargs)

    def __repr__(self):
        return (f"Slope(height={self.height}, angle={self.angle:.1f}, "
                f"length={self.length:.2f})")
