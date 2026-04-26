"""Tests for polyline sliding analysis (transfer coefficient method)."""

import math
import pytest
from pyslopex import Slope, Material


class TestPolylineBasic:
    """Basic polyline analysis tests."""

    def test_simple_two_segment(self):
        """Test with a simple two-segment failure surface."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        # Failure surface: from upper surface through one bend to toe
        points = [(15, 10), (8, 5), (0, 0)]
        result = s.analyse_polyline(failure_surface_points=points)
        assert result.fos > 0
        assert result.fos < 20

    def test_implicit_vs_explicit(self):
        """Implicit and explicit methods should give similar results."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        points = [(15, 10), (8, 5), (0, 0)]

        r_implicit = s.analyse_polyline(
            failure_surface_points=points, method='implicit')
        r_explicit = s.analyse_polyline(
            failure_surface_points=points, method='explicit')

        # Should be within 10% of each other
        diff = abs(r_implicit.fos - r_explicit.fos)
        avg = (r_implicit.fos + r_explicit.fos) / 2
        assert diff / avg < 0.10

    def test_three_segments(self):
        """Test with three segments (4 points)."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        points = [(15, 10), (11, 7), (5, 3), (0, 0)]
        result = s.analyse_polyline(failure_surface_points=points)
        assert result.fos > 0

    def test_thrust_distribution(self):
        """Thrust distribution should have one more entry than slices."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        points = [(15, 10), (8, 5), (0, 0)]
        result = s.analyse_polyline(failure_surface_points=points)

        # 2 segments → 3 thrust values (including initial 0)
        assert len(result.thrust_distribution) == 3
        assert result.thrust_distribution[0] == 0.0

    def test_high_cohesion_gives_high_fos(self):
        """High cohesion should give high FOS."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=50, depth_to_bottom=10)
        s.set_materials(m)

        points = [(15, 10), (8, 5), (0, 0)]
        result = s.analyse_polyline(failure_surface_points=points)
        assert result.fos > 3.0

    def test_water_table_reduces_fos(self):
        """Water table should reduce FOS."""
        s1 = Slope(height=10, angle=45)
        m1 = Material(unit_weight=20, friction_angle=30,
                      cohesion=10, depth_to_bottom=10)
        s1.set_materials(m1)

        s2 = Slope(height=10, angle=45)
        m2 = Material(unit_weight=20, friction_angle=30,
                      cohesion=10, depth_to_bottom=10)
        s2.set_materials(m2)
        s2.set_water_table(depth=2)

        points = [(15, 10), (8, 5), (0, 0)]
        r1 = s1.analyse_polyline(failure_surface_points=points)
        r2 = s2.analyse_polyline(failure_surface_points=points)
        assert r1.fos > r2.fos

    def test_two_materials(self):
        """Test with two material layers."""
        s = Slope(height=10, angle=45)
        m1 = Material(unit_weight=20, friction_angle=35,
                      cohesion=15, depth_to_bottom=5)
        m2 = Material(unit_weight=18, friction_angle=25,
                      cohesion=10, depth_to_bottom=10)
        s.set_materials(m1, m2)

        points = [(15, 10), (8, 5), (0, 0)]
        result = s.analyse_polyline(failure_surface_points=points)
        assert result.fos > 0

    def test_udl_reduces_fos(self):
        """Adding a load should reduce FOS."""
        from pyslopex import Udl

        s1 = Slope(height=10, angle=45)
        m1 = Material(unit_weight=20, friction_angle=30,
                      cohesion=10, depth_to_bottom=10)
        s1.set_materials(m1)

        s2 = Slope(height=10, angle=45)
        m2 = Material(unit_weight=20, friction_angle=30,
                      cohesion=10, depth_to_bottom=10)
        s2.set_materials(m2)
        s2.set_udls(Udl(magnitude=50))

        points = [(15, 10), (8, 5), (0, 0)]
        r1 = s1.analyse_polyline(failure_surface_points=points)
        r2 = s2.analyse_polyline(failure_surface_points=points)
        assert r1.fos > r2.fos

    def test_no_points_raises_error(self):
        """Should raise error if no points provided and auto_search is False."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        with pytest.raises(ValueError):
            s.analyse_polyline()


class TestPolylineAutoSearch:
    """Test automatic critical surface search."""

    def test_auto_search_returns_result(self):
        """Auto search should return a valid result."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        result = s.analyse_polyline(
            auto_search=True, num_internal_points=1,
            grid_resolution=8)
        assert result.fos > 0
        assert len(result.failure_surface) >= 2

    def test_auto_search_fos_reasonable(self):
        """Auto search FOS should be reasonable."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        result = s.analyse_polyline(
            auto_search=True, num_internal_points=1,
            grid_resolution=8)
        # FOS should be between 0.5 and 5 for this configuration
        assert 0.5 < result.fos < 5.0


class TestPolylineSingleSegment:
    """Test polyline with a single straight segment (should match planar)."""

    def test_single_segment_matches_planar(self):
        """A single-segment polyline should approximate a planar failure."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        # Single segment from (14.14, 10) to (0, 0) → angle ≈ 35.26°
        theta = 35.0
        x_entry = 10 / math.tan(math.radians(theta))
        points = [(x_entry, 10), (0, 0)]

        r_poly = s.analyse_polyline(
            failure_surface_points=points, method='implicit')
        r_plan = s.analyse_planar(min_angle=34, max_angle=36, num_angles=5)

        # Should be within 15% (difference due to slice discretization)
        diff = abs(r_poly.fos - r_plan.fos)
        avg = (r_poly.fos + r_plan.fos) / 2
        assert diff / avg < 0.15
