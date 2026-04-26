"""pyslopex - Slope stability analysis using planar and polyline sliding methods.

Provides two methods for slope stability analysis:
    - Planar sliding method (平面滑动法): Assumes a planar failure surface.
    - Polyline sliding method (折线滑动法 / 传递系数法): Uses the transfer
      coefficient method for piecewise linear failure surfaces.

Usage:
    from pyslopex import Slope, Material, Udl, LineLoad

    s = Slope(height=10, angle=45)
    s.set_materials(Material(unit_weight=20, friction_angle=30,
                             cohesion=10, depth_to_bottom=10))

    result = s.analyse_planar()
    print(f"Planar FOS: {result.fos:.3f}")
"""

from .slope import Slope
from .material import Material
from .loads import Udl, LineLoad

__version__ = "0.1.0"
__all__ = ["Slope", "Material", "Udl", "LineLoad"]
