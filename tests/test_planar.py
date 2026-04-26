"""Tests for planar sliding analysis."""

import math
import pytest
from pyslopex import Slope, Material


class TestPlanarBasic:
    """Basic planar analysis tests."""

    def test_single_material(self):
        """Test planar analysis with a single uniform soil layer."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        result = s.analyse_planar()
        assert result.fos > 0
        assert result.fos < 20
        assert 0 < result.critical_angle < 45
        assert len(result.failure_plane) == 2
        assert len(result.sliding_wedge) == 3

    def test_cohesionless_soil(self):
        """For pure friction soil (c=0), critical angle should be close
        to the slope angle since FOS = tan(φ) / tan(θ) is minimum at θ = β."""
        s = Slope(height=10, angle=30)
        m = Material(unit_weight=20, friction_angle=35,
                     cohesion=0, depth_to_bottom=10)
        s.set_materials(m)

        result = s.analyse_planar()
        # FOS at slope angle = tan(35)/tan(30) = 1.21
        assert result.fos < 1.5
        assert result.critical_angle > 25  # Close to slope angle

    def test_high_cohesion_gives_high_fos(self):
        """High cohesion should give high FOS."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=50, depth_to_bottom=10)
        s.set_materials(m)

        result = s.analyse_planar()
        assert result.fos > 3.0

    def test_steeper_slope_lower_fos(self):
        """Steeper slope should give lower FOS for the same soil."""
        s1 = Slope(height=10, angle=30)
        m1 = Material(unit_weight=20, friction_angle=30,
                      cohesion=10, depth_to_bottom=10)
        s1.set_materials(m1)

        s2 = Slope(height=10, angle=60)
        m2 = Material(unit_weight=20, friction_angle=30,
                      cohesion=10, depth_to_bottom=10)
        s2.set_materials(m2)

        r1 = s1.analyse_planar()
        r2 = s2.analyse_planar()
        assert r1.fos > r2.fos

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

        r1 = s1.analyse_planar()
        r2 = s2.analyse_planar()
        assert r1.fos > r2.fos

    def test_two_materials(self):
        """Test with two material layers."""
        s = Slope(height=10, angle=45)
        m1 = Material(unit_weight=20, friction_angle=35,
                      cohesion=15, depth_to_bottom=5)
        m2 = Material(unit_weight=18, friction_angle=25,
                      cohesion=10, depth_to_bottom=10)
        s.set_materials(m1, m2)

        result = s.analyse_planar()
        assert result.fos > 0
        assert result.fos < 20

    def test_udl_reduces_fos(self):
        """Adding a uniform load should reduce FOS."""
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

        r1 = s1.analyse_planar()
        r2 = s2.analyse_planar()
        assert r1.fos > r2.fos

    def test_all_results_populated(self):
        """Check that all_results is populated."""
        s = Slope(height=10, angle=45)
        m = Material(unit_weight=20, friction_angle=30,
                     cohesion=10, depth_to_bottom=10)
        s.set_materials(m)

        result = s.analyse_planar(num_angles=50)
        assert len(result.all_results) == 50
        for angle, fos in result.all_results:
            assert 0 < angle < 45
            assert fos > 0


class TestPlanarGeometry:
    """Test slope geometry for planar analysis."""

    def test_slope_by_length(self):
        """Test creating a slope by length."""
        s = Slope(height=10, angle=None, length=10)
        assert abs(s.angle - 45.0) < 0.01
        assert abs(s.length - 10.0) < 0.01

    def test_slope_by_angle(self):
        """Test creating a slope by angle."""
        s = Slope(height=10, angle=45, length=None)
        assert abs(s.length - 10.0) < 0.01

    def test_invalid_slope(self):
        """Test that neither or both angle/length raises error."""
        with pytest.raises(ValueError):
            Slope(height=10)
        with pytest.raises(ValueError):
            Slope(height=10, angle=45, length=10)

    def test_toe_and_crest(self):
        """Test toe and crest coordinates."""
        s = Slope(height=10, angle=45)
        assert s.toe == (0.0, 0.0)
        assert abs(s.crest[0] - 10.0) < 0.01
        assert abs(s.crest[1] - 10.0) < 0.01


class TestPlanarValidation:
    """Validation tests against known analytical solutions."""

    def test_dry_cohesionless_fos_formula(self):
        """For dry cohesionless soil, FOS = tan(φ) / tan(θ).
        The minimum FOS occurs at θ = slope angle."""
        phi = 35.0
        beta = 30.0
        s = Slope(height=10, angle=beta)
        m = Material(unit_weight=20, friction_angle=phi,
                     cohesion=0, depth_to_bottom=10)
        s.set_materials(m)

        result = s.analyse_planar()
        expected_fos = math.tan(math.radians(phi)) / math.tan(math.radians(beta))
        # Should be within 10% of analytical solution
        assert abs(result.fos - expected_fos) / expected_fos < 0.10
