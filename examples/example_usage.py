"""Example usage of the pyslopex library.

Demonstrates:
    1. Planar sliding analysis (平面滑动法)
    2. Polyline sliding analysis (折线滑动法 / 传递系数法)
    3. Effect of water table
    4. Effect of external loads
    5. Multi-layer soil
    6. Automatic critical surface search
"""

import sys
import os

# Fix Windows encoding for Chinese characters
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

# Add parent directory to path if running directly
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyslopex import Slope, Material, Udl, LineLoad


def example_1_planar_basic():
    """Example 1: Basic planar sliding analysis."""
    print("=" * 60)
    print("Example 1: Planar Sliding Analysis (平面滑动法)")
    print("=" * 60)

    s = Slope(height=10, angle=45)
    m = Material(unit_weight=20, friction_angle=30,
                 cohesion=10, depth_to_bottom=10,
                 name="Silty Clay")
    s.set_materials(m)

    result = s.analyse_planar()
    print(f"  Slope: H={s.height}m, angle={s.angle:.1f}°")
    print(f"  Soil: γ={m.unit_weight} kN/m³, φ={m.friction_angle}°, "
          f"c={m.cohesion} kPa")
    print(f"  Critical FOS: {result.fos:.3f}")
    print(f"  Critical angle: {result.critical_angle:.1f}°")
    print(f"  Failure plane: {result.failure_plane}")
    print()


def example_2_polyline_basic():
    """Example 2: Polyline sliding analysis with transfer coefficient method."""
    print("=" * 60)
    print("Example 2: Polyline Sliding (折线滑动法 / 传递系数法)")
    print("=" * 60)

    s = Slope(height=10, angle=45)
    m = Material(unit_weight=20, friction_angle=30,
                 cohesion=10, depth_to_bottom=10)
    s.set_materials(m)

    # Define failure surface points from entry to exit
    points = [(15, 10), (8, 5), (0, 0)]

    result = s.analyse_polyline(failure_surface_points=points)
    print(f"  Failure surface points: {points}")
    print(f"  FOS (implicit): {result.fos:.3f}")
    print(f"  Thrust at toe: {result.thrust_at_toe:.4f}")
    print(f"  Thrust distribution: "
          f"{[f'{t:.2f}' for t in result.thrust_distribution]}")
    print()

    # Also try explicit method
    result_exp = s.analyse_polyline(
        failure_surface_points=points, method='explicit')
    print(f"  FOS (explicit): {result_exp.fos:.3f}")
    print()


def example_3_water_table():
    """Example 3: Effect of water table on stability."""
    print("=" * 60)
    print("Example 3: Water Table Effect")
    print("=" * 60)

    s = Slope(height=10, angle=45)
    m = Material(unit_weight=20, friction_angle=30,
                 cohesion=10, depth_to_bottom=10)
    s.set_materials(m)

    # Without water table
    r1 = s.analyse_planar()
    print(f"  Without water table: FOS = {r1.fos:.3f}")

    # With water table at 2m depth
    s.set_water_table(depth=2)
    r2 = s.analyse_planar()
    print(f"  Water table at 2m depth: FOS = {r2.fos:.3f}")

    # With water table at surface (worst case)
    s.set_water_table(depth=0)
    r3 = s.analyse_planar()
    print(f"  Water table at surface: FOS = {r3.fos:.3f}")

    print(f"  FOS reduction: {(1 - r3.fos/r1.fos)*100:.1f}%")
    print()


def example_4_external_loads():
    """Example 4: Effect of external loads."""
    print("=" * 60)
    print("Example 4: External Load Effect")
    print("=" * 60)

    s = Slope(height=10, angle=45)
    m = Material(unit_weight=20, friction_angle=30,
                 cohesion=10, depth_to_bottom=10)
    s.set_materials(m)

    # No load
    r1 = s.analyse_planar()
    print(f"  No load: FOS = {r1.fos:.3f}")

    # With UDL
    s.set_udls(Udl(magnitude=50))
    r2 = s.analyse_planar()
    print(f"  UDL 50 kPa at crest: FOS = {r2.fos:.3f}")

    # With line load
    s.set_udls()
    s.set_line_loads(LineLoad(magnitude=30, offset=2))
    r3 = s.analyse_planar()
    print(f"  Line load 30 kN/m at 2m from crest: FOS = {r3.fos:.3f}")
    print()


def example_5_multilayer():
    """Example 5: Multi-layer soil."""
    print("=" * 60)
    print("Example 5: Multi-Layer Soil")
    print("=" * 60)

    s = Slope(height=10, angle=45)
    m1 = Material(unit_weight=20, friction_angle=35,
                  cohesion=15, depth_to_bottom=4,
                  name="Dense Sand")
    m2 = Material(unit_weight=19, friction_angle=28,
                  cohesion=12, depth_to_bottom=7,
                  name="Silty Clay")
    m3 = Material(unit_weight=18, friction_angle=22,
                  cohesion=8, depth_to_bottom=10,
                  name="Soft Clay")
    s.set_materials(m1, m2, m3)

    # Planar analysis
    r_plan = s.analyse_planar()
    print(f"  Planar FOS: {r_plan.fos:.3f}")
    print(f"  Critical angle: {r_plan.critical_angle:.1f}°")

    # Polyline analysis
    points = [(16, 10), (11, 7), (5, 3), (0, 0)]
    r_poly = s.analyse_polyline(failure_surface_points=points)
    print(f"  Polyline FOS: {r_poly.fos:.3f}")
    print()


def example_6_auto_search():
    """Example 6: Automatic critical polyline surface search."""
    print("=" * 60)
    print("Example 6: Auto Search for Critical Polyline Surface")
    print("=" * 60)

    s = Slope(height=10, angle=45)
    m = Material(unit_weight=20, friction_angle=30,
                 cohesion=10, depth_to_bottom=10)
    s.set_materials(m)

    result = s.analyse_polyline(
        auto_search=True, num_internal_points=2,
        grid_resolution=10)

    print(f"  Best FOS found: {result.fos:.3f}")
    print(f"  Critical surface: {result.failure_surface}")
    print(f"  Method: {result.method}")
    print()


def example_7_plotting():
    """Example 7: Generate plots (requires matplotlib)."""
    print("=" * 60)
    print("Example 7: Generate Plots")
    print("=" * 60)

    s = Slope(height=10, angle=45)
    m = Material(unit_weight=20, friction_angle=30,
                 cohesion=10, depth_to_bottom=10)
    s.set_materials(m)
    s.set_water_table(depth=3)

    # Planar result plot
    r_plan = s.analyse_planar()
    print(f"  Planar FOS: {r_plan.fos:.3f}")
    print(f"  Plot with: s.plot_critical_planar(r_plan)")

    # Polyline result plot
    points = [(15, 10), (8, 5), (0, 0)]
    r_poly = s.analyse_polyline(failure_surface_points=points)
    print(f"  Polyline FOS: {r_poly.fos:.3f}")
    print(f"  Plot with: s.plot_polyline(points)")

    # Uncomment to show plots:
    # s.plot_critical_planar(r_plan)
    # s.plot_polyline(points, r_poly)
    print()


if __name__ == "__main__":
    example_1_planar_basic()
    example_2_polyline_basic()
    example_3_water_table()
    example_4_external_loads()
    example_5_multilayer()
    example_6_auto_search()
    example_7_plotting()

    print("All examples completed.")
