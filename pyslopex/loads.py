"""Load definitions for slope stability analysis."""


class Udl:
    """Uniform distributed load applied on the top surface.

    Args:
        magnitude: Load magnitude (kPa).
        offset: Horizontal offset from crest of slope (m).
        length: Length of the load (m). None means infinite extent.
        dynamic_offset: If True, the load is moved in dynamic analysis.
        color: Optional color for plotting.
    """

    def __init__(self, magnitude, offset=0, length=None,
                 dynamic_offset=False, color=None):
        self.magnitude = magnitude
        self.offset = offset
        self.length = length
        self.dynamic_offset = dynamic_offset
        self.color = color

    def __repr__(self):
        return (f"Udl(magnitude={self.magnitude}, offset={self.offset}, "
                f"length={self.length})")


class LineLoad:
    """Concentrated line load applied on the top surface.

    Args:
        magnitude: Load magnitude (kN/m).
        offset: Horizontal offset from crest of slope (m).
        color: Optional color for plotting.
    """

    def __init__(self, magnitude, offset=0, color=None):
        self.magnitude = magnitude
        self.offset = offset
        self.color = color

    def __repr__(self):
        return f"LineLoad(magnitude={self.magnitude}, offset={self.offset})"
