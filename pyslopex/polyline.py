"""Polyline sliding method (折线滑动法 / 传递系数法) for slope stability.

Implements the transfer coefficient method (传递系数法), widely used in
Chinese geotechnical practice for analyzing slopes with non-circular
(piecewise linear) failure surfaces.

The method divides the sliding mass into blocks at each vertex of the
polyline failure surface. The thrust from upper blocks is transferred
to lower blocks via a transfer coefficient.

Two solution methods are provided:
    - 'implicit': FOS appears in the transfer coefficient; solved by
      bisection (standard Chinese code GB 50330 approach).
    - 'explicit': FOS calculated directly without iteration.
"""

import math
import numpy as np

from .utils import ground_elevation, polygon_area


class PolylineResult:
    """Result from polyline sliding analysis.

    Attributes:
        fos: Factor of safety.
        failure_surface: Failure surface points [(x,y), ...].
        slice_data: Detailed data for each slice/block.
        thrust_at_toe: Residual thrust at the exit point.
        thrust_distribution: Thrust at each slice boundary.
        method: Analysis method used ('implicit' or 'explicit').
    """

    def __init__(self):
        self.fos = float('inf')
        self.failure_surface = []
        self.slice_data = []
        self.thrust_at_toe = 0.0
        self.thrust_distribution = []
        self.method = 'implicit'

    def __repr__(self):
        return (f"PolylineResult(fos={self.fos:.4f}, "
                f"method='{self.method}')")


class _SliceData:
    """Data for a single slice/block in polyline analysis."""

    def __init__(self):
        self.alpha = 0.0        # Base angle (rad)
        self.alpha_deg = 0.0    # Base angle (deg)
        self.length = 0.0       # Base length (m)
        self.weight = 0.0       # Weight per unit width (kN/m)
        self.c = 0.0            # Cohesion at base (kPa)
        self.phi = 0.0          # Friction angle at base (deg)
        self.U = 0.0            # Pore pressure force on base (kN/m)
        self.driving = 0.0      # W * sin(alpha)
        self.resisting = 0.0    # c*L + (W*cos(alpha) - U)*tan(phi)
        self.p_top = None       # Upper vertex
        self.p_bot = None       # Lower vertex


class PolylineAnalysis:
    """Polyline sliding analysis using the transfer coefficient method.

    Args:
        slope: Slope object.
    """

    def __init__(self, slope):
        self.slope = slope

    def analyse(self, failure_surface_points=None, method='implicit',
                auto_search=False, num_internal_points=2,
                grid_resolution=20, entry_x_range=None):
        """Run polyline sliding analysis.

        Args:
            failure_surface_points: List of (x, y) tuples defining the
                failure surface, ordered from entry point (upper-right,
                on top surface) to exit point (lower-left, at toe).
            method: 'implicit' (iterative bisection) or 'explicit'
                (direct calculation).
            auto_search: If True, search for the critical polyline
                surface automatically.
            num_internal_points: Number of internal break points for
                auto search (1-3 recommended).
            grid_resolution: Grid resolution for auto search.
            entry_x_range: (x_min, x_max) range for entry point search.
                Default is from crest to crest + height.

        Returns:
            PolylineResult.
        """
        if failure_surface_points is not None:
            slices = self._prepare_slices(failure_surface_points)
            result = self._solve(slices, method)
            result.failure_surface = list(failure_surface_points)
            return result
        elif auto_search:
            return self._search_critical(
                method, num_internal_points, grid_resolution, entry_x_range)
        else:
            raise ValueError(
                "Either failure_surface_points must be provided "
                "or auto_search must be True")

    def _prepare_slices(self, points):
        """Prepare slice data from failure surface points.

        Points are ordered from upper-right (entry) to lower-left (exit).
        """
        slope = self.slope
        H = slope.height
        n_segments = len(points) - 1
        if n_segments < 1:
            raise ValueError(
                "At least 2 points needed for polyline failure surface")

        slices = []
        for i in range(n_segments):
            p_upper = points[i]      # Upper/right point
            p_lower = points[i + 1]  # Lower/left point

            sd = _SliceData()
            sd.p_top = p_upper
            sd.p_bot = p_lower

            # Base angle: positive when base slopes downward from right to left
            dx = p_upper[0] - p_lower[0]
            dy = p_upper[1] - p_lower[1]
            sd.alpha = math.atan2(dy, dx)
            sd.alpha_deg = math.degrees(sd.alpha)

            # Base length
            seg_dx = p_lower[0] - p_upper[0]
            seg_dy = p_lower[1] - p_upper[1]
            sd.length = math.sqrt(seg_dx ** 2 + seg_dy ** 2)

            # Block weight (area between ground surface and failure surface)
            sd.weight = self._calculate_block_weight(p_upper, p_lower)

            # External loads on this block
            sd.weight += self._calculate_block_loads(p_upper, p_lower)

            # Shear strength at base midpoint
            y_mid = (p_upper[1] + p_lower[1]) / 2.0
            mat = slope.get_material_at_y(y_mid)
            if mat is None:
                raise ValueError(f"No material at elevation {y_mid:.2f}")
            sd.c = mat.cohesion
            sd.phi = mat.friction_angle

            # Pore pressure on base
            sd.U = self._calculate_base_pore_pressure(p_upper, p_lower)

            # Driving and resisting components
            sd.driving = sd.weight * math.sin(sd.alpha)
            cos_alpha = math.cos(sd.alpha)
            effective_normal = sd.weight * cos_alpha - sd.U
            effective_normal = max(0.0, effective_normal)
            sd.resisting = (sd.c * sd.length
                            + effective_normal
                            * math.tan(math.radians(sd.phi)))

            slices.append(sd)

        return slices

    def _solve(self, slices, method='implicit'):
        """Solve for factor of safety.

        Args:
            slices: List of _SliceData objects.
            method: 'implicit' or 'explicit'.

        Returns:
            PolylineResult.
        """
        result = PolylineResult()
        result.slice_data = slices
        result.method = method

        if not slices:
            result.fos = float('inf')
            return result

        if method == 'implicit':
            self._solve_implicit(slices, result)
        else:
            self._solve_explicit(slices, result)

        return result

    def _solve_implicit(self, slices, result):
        """Solve using implicit method (bisection).

        FOS appears in the transfer coefficient:
            ψ_{i-1} = cos(α_{i-1} - α_i) - sin(α_{i-1} - α_i) * tan(φ_i) / F_s

        Binary search for F_s such that the thrust at the last slice P_n = 0.
        """
        fos_low = 0.01
        fos_high = 20.0
        tol = 1e-8
        max_iter = 200

        for _ in range(max_iter):
            fos_mid = (fos_low + fos_high) / 2.0
            p_last = self._calculate_thrust(slices, fos_mid)

            if abs(p_last) < tol:
                break
            if p_last > 0:
                fos_high = fos_mid
            else:
                fos_low = fos_mid

            if (fos_high - fos_low) < tol:
                break

        result.fos = (fos_low + fos_high) / 2.0
        result.thrust_at_toe = self._calculate_thrust(slices, result.fos)
        result.thrust_distribution = self._get_thrust_distribution(
            slices, result.fos)

    def _solve_explicit(self, slices, result):
        """Solve using explicit method.

        Transfer coefficient does not depend on FOS:
            ψ_{i-1} = cos(α_{i-1} - α_i) - sin(α_{i-1} - α_i) * tan(φ_i)

        When ψ < 0 (large angle difference), clamp to 0 per Chinese practice.
        If the denominator is negative, the explicit method is not applicable
        and we fall back to the implicit method.
        """
        n = len(slices)

        # Calculate transfer coefficients (clamp negative to 0)
        psi_list = []
        for i in range(n - 1):
            alpha_i = slices[i].alpha
            alpha_next = slices[i + 1].alpha
            phi_next = math.radians(slices[i + 1].phi)
            psi_i = (math.cos(alpha_i - alpha_next)
                     - math.sin(alpha_i - alpha_next) * math.tan(phi_next))
            # Clamp: negative transfer coefficient is physically unrealistic
            psi_i = max(0.0, psi_i)
            psi_list.append(psi_i)

        # Cumulative products from bottom to top:
        # Psi[i] = product of psi from i to n-2
        Psi = [1.0] * n
        for i in range(n - 2, -1, -1):
            Psi[i] = psi_list[i] * Psi[i + 1]

        # Calculate FOS
        numerator = 0.0
        denominator = 0.0
        for i in range(n):
            numerator += slices[i].resisting * Psi[i]
            denominator += slices[i].driving * Psi[i]

        if denominator <= 1e-12:
            # Explicit method not applicable, fall back to implicit
            result.method = 'explicit->implicit_fallback'
            self._solve_implicit(slices, result)
            return

        result.fos = numerator / denominator
        result.thrust_at_toe = 0.0
        result.thrust_distribution = self._get_thrust_distribution_explicit(
            slices, result.fos, psi_list)

    def _calculate_thrust(self, slices, fos):
        """Calculate the thrust at the last slice for a given FOS (implicit).

        P_i = F_s * T_i + P_{i-1} * ψ_{i-1} - R_i

        where:
            T_i = W_i * sin(α_i)
            R_i = c_i * L_i + (W_i * cos(α_i) - U_i) * tan(φ_i)
            ψ_{i-1} = cos(α_{i-1} - α_i) - sin(α_{i-1} - α_i) * tan(φ_i) / F_s
        """
        P = 0.0
        prev_alpha = 0.0

        for i, s in enumerate(slices):
            if i == 0:
                psi = 1.0
            else:
                psi = (math.cos(prev_alpha - s.alpha)
                       - math.sin(prev_alpha - s.alpha)
                       * math.tan(math.radians(s.phi)) / fos)

            P = fos * s.driving + P * psi - s.resisting
            # 国标做法: 推力不小于0(负推力表示该块体自身稳定，不向下方传递)
            if i < len(slices) - 1:
                P = max(0.0, P)
            prev_alpha = s.alpha

        return P

    def _get_thrust_distribution(self, slices, fos):
        """Get thrust at each slice boundary (implicit method)."""
        distribution = [0.0]
        P = 0.0
        prev_alpha = 0.0

        for i, s in enumerate(slices):
            if i == 0:
                psi = 1.0
            else:
                psi = (math.cos(prev_alpha - s.alpha)
                       - math.sin(prev_alpha - s.alpha)
                       * math.tan(math.radians(s.phi)) / fos)

            P = fos * s.driving + P * psi - s.resisting
            if i < len(slices) - 1:
                P = max(0.0, P)
            distribution.append(P)
            prev_alpha = s.alpha

        return distribution

    def _get_thrust_distribution_explicit(self, slices, fos, psi_list):
        """Get thrust at each slice boundary (explicit method)."""
        distribution = [0.0]
        P = 0.0

        for i, s in enumerate(slices):
            if i == 0:
                psi = 1.0
            else:
                psi = psi_list[i - 1]

            P = fos * s.driving + P * psi - s.resisting
            if i < len(slices) - 1:
                P = max(0.0, P)
            distribution.append(P)

        return distribution

    def _calculate_block_weight(self, p_upper, p_lower):
        """Calculate the weight of a block between ground surface and
        failure surface, bounded by vertical lines through p_upper and p_lower.

        Uses vertical strips to handle layered soils.
        """
        slope = self.slope
        H = slope.height

        x_left = min(p_upper[0], p_lower[0])
        x_right = max(p_upper[0], p_lower[0])
        if x_right - x_left < 1e-9:
            return 0.0

        num_strips = 20
        dx = (x_right - x_left) / num_strips
        total_weight = 0.0

        for j in range(num_strips):
            x_mid = x_left + (j + 0.5) * dx

            # Ground surface elevation
            y_top = ground_elevation(x_mid, slope.length, H)

            # Failure surface elevation (linear interpolation)
            if abs(p_upper[0] - p_lower[0]) < 1e-9:
                y_bot = min(p_upper[1], p_lower[1])
            else:
                t = ((x_mid - p_upper[0])
                     / (p_lower[0] - p_upper[0]))
                t = max(0.0, min(1.0, t))
                y_bot = p_upper[1] + t * (p_lower[1] - p_upper[1])

            if y_top <= y_bot:
                continue

            # Calculate strip weight by layers
            layers = slope.get_layer_boundaries()
            strip_weight = 0.0

            for y_ltop, y_lbot, mat in layers:
                y_clip_top = min(y_top, y_ltop)
                y_clip_bot = max(y_bot, y_lbot)
                if y_clip_top <= y_clip_bot:
                    continue
                strip_weight += (y_clip_top - y_clip_bot) * dx * mat.unit_weight

            # Below defined layers
            if layers:
                lowest_y = layers[-1][1]
                if y_bot < lowest_y:
                    h = min(lowest_y, y_top) - y_bot
                    if h > 0 and slope.materials:
                        strip_weight += h * dx * slope.materials[-1].unit_weight

            total_weight += strip_weight

        return total_weight

    def _calculate_block_loads(self, p_upper, p_lower):
        """Calculate external loads on a block."""
        slope = self.slope
        x_left = min(p_upper[0], p_lower[0])
        x_right = max(p_upper[0], p_lower[0])
        x_crest = slope.length
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

    def _calculate_base_pore_pressure(self, p_upper, p_lower):
        """Calculate total pore water pressure force on a slice base."""
        slope = self.slope
        if slope.water_table_depth is None:
            return 0.0

        H = slope.height
        y_water = H - slope.water_table_depth
        gamma_w = 9.81

        # Pore pressure at each endpoint
        u_top = max(0.0, gamma_w * (y_water - p_upper[1]))
        u_bot = max(0.0, gamma_w * (y_water - p_lower[1]))

        seg_dx = p_lower[0] - p_upper[0]
        seg_dy = p_lower[1] - p_upper[1]
        base_length = math.sqrt(seg_dx ** 2 + seg_dy ** 2)

        # Average pore pressure * base length
        return (u_top + u_bot) / 2.0 * base_length

    def _search_critical(self, method, num_internal_points,
                         grid_resolution, entry_x_range):
        """Search for the critical polyline failure surface.

        Strategy:
            1. Define entry points on the upper surface.
            2. Place internal break points at evenly spaced x-positions.
            3. Search over y-positions of internal points.
            4. Exit point is always at the toe.
        """
        slope = self.slope
        H = slope.height
        L = slope.length

        if entry_x_range is None:
            entry_x_range = (L + 0.5, L + max(H, 3.0))

        best_result = PolylineResult()
        best_result.fos = float('inf')
        best_points = None

        x_entries = np.linspace(
            entry_x_range[0], entry_x_range[1], grid_resolution)

        for x_entry in x_entries:
            # Internal x-positions evenly spaced between entry and toe
            x_positions = np.linspace(
                x_entry, 0, num_internal_points + 2)  # +2 for entry and exit

            # Create a grid of y-values for each internal point
            y_grid = np.linspace(0.1, H * 0.95, grid_resolution)

            # Build all combinations using iterative approach
            points_list = self._grid_search(
                x_entry, x_positions, y_grid, num_internal_points,
                H, method)

            for points in points_list:
                try:
                    slices = self._prepare_slices(points)
                    result = self._solve(slices, method)
                    result.failure_surface = points
                    if result.fos < best_result.fos:
                        best_result = result
                        best_points = points
                except (ValueError, ZeroDivisionError):
                    continue

        if best_points is not None:
            best_result.failure_surface = best_points

        return best_result

    def _grid_search(self, x_entry, x_positions, y_grid,
                      num_internal_points, H, method):
        """Generate trial failure surfaces for grid search.

        Returns list of point lists for valid failure surfaces.
        """
        slope = self.slope
        results = []

        # Internal point indices (exclude first=entry, last=exit)
        internal_indices = list(range(1, len(x_positions) - 1))

        if not internal_indices:
            # Direct line from entry to toe
            points = [(x_entry, H), (0.0, 0.0)]
            results.append(points)
            return results

        # Build combinations iteratively
        from itertools import product

        internal_y_combos = list(product(y_grid, repeat=len(internal_indices)))

        for y_combo in internal_y_combos:
            points = [(x_entry, H)]
            valid = True
            prev_y = H
            prev_x = x_entry

            for idx, yi in zip(internal_indices, y_combo):
                xi = x_positions[idx]
                # Ensure monotonically decreasing y and x
                if yi >= prev_y or xi >= prev_x or yi <= 0:
                    valid = False
                    break
                # Must be below ground surface
                y_ground = ground_elevation(xi, slope.length, H)
                if yi > y_ground:
                    valid = False
                    break
                points.append((float(xi), float(yi)))
                prev_y = yi
                prev_x = xi

            if not valid:
                continue

            # Exit at toe
            if prev_y <= 0.1:
                continue
            points.append((0.0, 0.0))

            if len(points) >= 2:
                results.append(points)

        return results
