/// Pixel coordinates (column, row), with the first pixel centred at `(0, 0)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PixelPos {
    pub x: f64,
    pub y: f64,
}

impl PixelPos {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PixelConvention {
    ZeroBased,
    Fits,
    /// Coordinates in (row, column) order, with the first pixel centred at `(0.5, 0.5)`.
    PixelEdgeRowMajor,
}

impl PixelConvention {
    pub fn to_internal(self, a: f64, b: f64) -> PixelPos {
        match self {
            Self::ZeroBased => PixelPos::new(a, b),
            Self::Fits => PixelPos::new(a - 1.0, b - 1.0),
            Self::PixelEdgeRowMajor => PixelPos::new(b - 0.5, a - 0.5),
        }
    }

    pub fn from_internal(self, p: PixelPos) -> (f64, f64) {
        match self {
            Self::ZeroBased => (p.x, p.y),
            Self::Fits => (p.x + 1.0, p.y + 1.0),
            Self::PixelEdgeRowMajor => (p.y + 0.5, p.x + 0.5),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conventions_round_trip_and_agree_on_first_pixel() {
        let first = PixelPos::new(0.0, 0.0);
        assert_eq!(PixelConvention::Fits.from_internal(first), (1.0, 1.0));
        assert_eq!(
            PixelConvention::PixelEdgeRowMajor.from_internal(first),
            (0.5, 0.5)
        );

        let p = PixelPos::new(12.25, 700.75);
        for c in [
            PixelConvention::ZeroBased,
            PixelConvention::Fits,
            PixelConvention::PixelEdgeRowMajor,
        ] {
            let (a, b) = c.from_internal(p);
            assert_eq!(c.to_internal(a, b), p, "{c:?}");
        }
        assert_eq!(
            PixelConvention::PixelEdgeRowMajor.from_internal(p),
            (701.25, 12.75)
        );
    }
}
