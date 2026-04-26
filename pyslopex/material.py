"""Material definition for slope stability analysis."""


class Material:
    """Represents a soil layer with geotechnical properties.

    Args:
        unit_weight: Unit weight of soil (kN/m³).
        friction_angle: Effective friction angle (degrees).
        cohesion: Effective cohesion (kPa).
        depth_to_bottom: Depth from top of slope to bottom of this layer (m).
        name: Optional name for the material.
        color: Optional color for plotting.
    """

    def __init__(self, unit_weight, friction_angle, cohesion, depth_to_bottom,
                 name=None, color=None):
        self.unit_weight = unit_weight
        self.friction_angle = friction_angle
        self.cohesion = cohesion
        self.depth_to_bottom = depth_to_bottom
        self.name = name or f"Soil(γ={unit_weight}, φ={friction_angle}, c={cohesion})"
        self.color = color

    def __repr__(self):
        return (f"Material(unit_weight={self.unit_weight}, "
                f"friction_angle={self.friction_angle}, "
                f"cohesion={self.cohesion}, "
                f"depth_to_bottom={self.depth_to_bottom})")
