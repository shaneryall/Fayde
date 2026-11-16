//! Pixel positions and the conversions between pixel conventions used by other tools.

/// A pixel position in Fayde's internal convention: `(x, y)` = (column, row), with the centre of
/// the first pixel at `(0, 0)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PixelPos {
    /// Column coordinate.
    pub x: f64,
    /// Row coordinate.
    pub y: f64,
}

impl PixelPos {
    /// Construct from internal-convention coordinates.
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

/// Pixel conventions Fayde reads and writes at its I/O edges.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PixelConvention {
    /// Fayde / NumPy / Astropy origin 0: first pixel centre at `(0, 0)`, `(x, y)` order.
    ZeroBased,
    /// FITS (and Fayde source lists): first pixel centre at `(1, 1)`, `(x, y)` order.
    Fits,
    /// tetra3 / cedar-solve / cedar-detect: first pixel centre at `(0.5, 0.5)`, given as
    /// `(y, x)` — row first.
    PixelEdgeRowMajor,
}

impl PixelConvention {
    /// Convert a coordinate pair written in this convention into a [`PixelPos`].
    ///
    /// For [`PixelConvention::PixelEdgeRowMajor`] the pair is `(row, column)`.
    pub fn to_internal(self, a: f64, b: f64) -> PixelPos {
        match self {
            Self::ZeroBased => PixelPos::new(a, b),
            Self::Fits => PixelPos::new(a - 1.0, b - 1.0),
            Self::PixelEdgeRowMajor => PixelPos::new(b - 0.5, a - 0.5),
        }
    }

    /// Convert a [`PixelPos`] into a coordinate pair in this convention.
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
        // Row-major swap: (row, col) = (700.75 + 0.5, 12.25 + 0.5)
        assert_eq!(
            PixelConvention::PixelEdgeRowMajor.from_internal(p),
            (701.25, 12.75)
        );
    }
}
