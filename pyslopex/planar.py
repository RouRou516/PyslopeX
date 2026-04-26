"""Planar sliding method (平面滑动法) for slope stability analysis.

Analyzes slope stability assuming a planar (straight) failure surface
passing through the toe of the slope. Searches for the critical plane
angle that gives the minimum factor of safety.

The factor of safety is calculated as:
    FOS = (Σ resisting forces) / (Σ driving forces)

For a planar surface at angle θ through the toe:
    - Driving force: W * sin(θ)
    - Resisting force: c * L + (W * cos(θ) - U) * tan(φ)
"""

import math
import numpy as np

from .utils import ground_elevation


class PlanarResult:
    """Result from planar sliding analysis.

    Attributes:
        fos: Minimum factor of safety.
        critical_angle: Critical failure plane angle (degrees).
        failure_plane: List of [(x1,y1), (x2,y2)] defining the failure plane.
        all_results: List of (angle, fos) for all trial angles.
        sliding_wedge: Vertices of the critical sliding wedge.
    """

    def __init__(self):
        self.fos = float('inf')
        self.critical_angle = 0.0
        self.failure_plane = []
        self.all_results = []
        self.sliding_wedge = []

    def __repr__(self):
        return (f"PlanarResult(fos={self.fos:.4f}, "
                f"critical_angle={self.critical_angle:.2f}°)")


class PlanarAnalysis:
    """Planar sliding method analysis.

    The failure surface is a plane passing through the toe of the slope.
    The method searches over a range of plane angles to find the
    minimum factor of safety.

    Args:
        slope: Slope object.
    """

    def __init__(self, slope):
        self.slope = slope

    def analyse(self, min_angle=5.0, max_angle=None, num_angles=200,
                num_slices=50):
        """Run the planar sliding analysis.

        Args:
            min_angle: Minimum failure plane angle to search (degrees).
            max_angle: Maximum failure plane angle (degrees).
                Default is slope angle - 0.1.
            num_angles: Number of trial angles.
            num_slices: Number of vertical slices for FOS calculation.

        Returns:
            PlanarResult with the minimum FOS and critical failure plane.
        """
        slope = self.slope
        beta = slope.angle

        if max_angle is None:
            max_angle = beta - 0.1
        if min_angle >= max_angle:
            min_angle = max(1.0, max_angle - 1.0)

        result = PlanarResult()
        angles = np.linspace(min_angle, max_angle, num_angles)

        for theta in angles:
            fos = self._calculate_fos(float(theta), num_slices)
            result.all_results.append((float(theta), fos))

            if fos < result.fos:
                result.fos = fos
                result.critical_angle = float(theta)

        # Compute critical failure plane coordinates
        theta_cr = math.radians(result.critical_angle)
        H = slope.height
        beta_rad = math.radians(beta)
        x_crest = H / math.tan(beta_rad) if math.tan(beta_rad) > 1e-9 else 1e9
        x_fail = H / math.tan(theta_cr) if math.tan(theta_cr) > 1e-9 else 1e9

        result.failure_plane = [(0.0, 0.0), (x_fail, float(H))]
        result.sliding_wedge = [
            (0.0, 0.0),
            (x_crest, float(H)),
            (x_fail, float(H)),
        ]

        return result

    def _calculate_fos(self, theta_deg, num_slices=50):
        """Calculate FOS for a planar failure at angle theta.

        Uses vertical slices to handle layered soils properly.

        Args:
            theta_deg: Failure plane angle (degrees).
            num_slices: Number of vertical slices.

        Returns:
            Factor of safety.
        """
        slope = self.slope
        theta = math.radians(theta_deg)
        beta_rad = math.radians(slope.angle)
        H = slope.height

        tan_theta = math.tan(theta)
        if tan_theta < 1e-9:
            return float('inf')

        x_fail = H / tan_theta
        x_crest = H / math.tan(beta_rad) if math.tan(beta_rad) > 1e-9 else 1e9

        # Failure plane must exit beyond the crest
        if x_fail <= x_crest + 0.001:
            return float('inf')

        dx = x_fail / num_slices
        cos_theta = math.cos(theta)
        sin_theta = math.sin(theta)

        total_driving = 0.0
        total_resisting = 0.0

        for i in range(num_slices):
            x_left = i * dx
            x_right = (i + 1) * dx
            x_mid = (x_left + x_right) / 2.0

            y_top = ground_elevation(x_mid, slope.length, H)
            y_bot = x_mid * tan_theta

            h = y_top - y_bot
            if h <= 1e-9:
                continue

            weight = self._calculate_slice_weight(
                x_left, x_right, y_bot, y_top)

            # External loads on top surface (only if slice is on top)
            weight += self._calculate_slice_loads(
                x_left, x_right, x_crest)

            # Base length along failure plane
            base_length = dx / cos_theta

            # Shear strength at base
            mat = slope.get_material_at_y(y_bot)
            if mat is None:
                continue
            c_base = mat.cohesion
            phi_base = math.radians(mat.friction_angle)

            # Pore water pressure at base
            u = self._get_pore_pressure(y_bot)
            pore_force = u * base_length

            # Force equilibrium
            driving = weight * sin_theta
            normal = weight * cos_theta
            effective_normal = max(0.0, normal - pore_force)
            resisting = (c_base * base_length
                         + effective_normal * math.tan(phi_base))

            total_driving += driving
            total_resisting += resisting

        if total_driving <= 1e-9:
            return float('inf')

        return total_resisting / total_driving

    def _calculate_slice_weight(self, x_left, x_right, y_bot, y_top):
        """Calculate the weight of a vertical slice considering layers.

        Args:
            x_left, x_right: Horizontal bounds of the slice.
            y_bot: Elevation of the slice base.
            y_top: Elevation of the slice top.

        Returns:
            Weight per unit width (kN/m).
        """
        slope = self.slope
        dx = x_right - x_left
        layers = slope.get_layer_boundaries()

        weight = 0.0
        for y_ltop, y_lbot, mat in layers:
            # Clip layer to slice height
            y_clip_top = min(y_top, y_ltop)
            y_clip_bot = max(y_bot, y_lbot)
            if y_clip_top <= y_clip_bot:
                continue
            layer_height = y_clip_top - y_clip_bot
            weight += layer_height * dx * mat.unit_weight

        # Handle area below defined layers (use last material)
        if layers:
            lowest_y = layers[-1][1]
        else:
            lowest_y = float(slope.height)

        if y_bot < lowest_y:
            extra_height = min(lowest_y, y_top) - y_bot
            if extra_height > 0:
                mat = slope.materials[-1] if slope.materials else None
                if mat:
                    weight += extra_height * dx * mat.unit_weight

        # Handle area above defined layers (use first material)
        if layers:
            highest_y = layers[0][0]
        else:
            highest_y = 0.0

        if y_top > highest_y and layers:
            extra_height = y_top - max(highest_y, y_bot)
            if extra_height > 0:
                weight += extra_height * dx * slope.materials[0].unit_weight

        return weight

    def _calculate_slice_loads(self, x_left, x_right, x_crest):
        """Calculate external loads on a slice.

        Args:
            x_left, x_right: Slice horizontal bounds.
            x_crest: X-coordinate of the crest.

        Returns:
            Additional weight from loads (kN/m).
        """
        slope = self.slope
        extra_weight = 0.0

        for udl in slope.udls:
            load_start = x_crest + udl.offset
            load_end = (load_start + udl.length
                        if udl.length is not None else float('inf'))
            overlap_start = max(x_left, load_start)
            overlap_end = min(x_right, load_end)
            if overlap_end > overlap_start:
                extra_weight += udl.magnitude * (overlap_end - overlap_start)

        for ll in slope.line_loads:
            load_x = x_crest + ll.offset
            if x_left <= load_x <= x_right:
                extra_weight += ll.magnitude

        return extra_weight

    def _get_pore_pressure(self, y):
        """Calculate pore water pressure at elevation y.

        Args:
            y: Elevation.

        Returns:
            Pore pressure (kPa). Zero if no water table or above it.
        """
        slope = self.slope
        if slope.water_table_depth is None:
            return 0.0

        y_water = slope.height - slope.water_table_depth
        if y >= y_water:
            return 0.0

        gamma_w = 9.81  # Unit weight of water (kN/m³)
        return gamma_w * (y_water - y)
