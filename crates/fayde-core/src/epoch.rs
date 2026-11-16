//! Astrometric epochs.

/// An epoch as a Julian year (TT/TDB, 365.25-day years from J2000.0).
///
/// Sub-second time-scale differences (UTC vs TT) are irrelevant at catalogue proper-motion
/// precision: 69 s of TT−UTC moves a 1″/yr star by 2 µas.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Epoch(f64);

/// Modified Julian Date of J2000.0 (2000-01-01T12:00 TT).
const MJD_J2000: f64 = 51_544.5;
const DAYS_PER_JULIAN_YEAR: f64 = 365.25;

impl Epoch {
    /// J2000.0.
    pub const J2000: Epoch = Epoch(2000.0);

    /// Julian years an observation or index epoch may plausibly take. A value outside is a
    /// clock error (a Julian Date in an MJD keyword, an unset camera clock), and catalogue
    /// proper motions mean nothing that far from the catalogue's own epoch.
    pub const PLAUSIBLE_YEARS: std::ops::RangeInclusive<f64> = 1800.0..=2200.0;

    /// Whether the epoch lies in [`Epoch::PLAUSIBLE_YEARS`].
    pub fn is_plausible(self) -> bool {
        Self::PLAUSIBLE_YEARS.contains(&self.0)
    }

    /// From a Julian year such as `2026.5`. Returns `None` if not finite.
    pub fn julian_year(year: f64) -> Option<Self> {
        year.is_finite().then_some(Self(year))
    }

    /// From a Modified Julian Date.
    pub fn from_mjd(mjd: f64) -> Option<Self> {
        Self::julian_year(2000.0 + (mjd - MJD_J2000) / DAYS_PER_JULIAN_YEAR)
    }

    /// From a Julian Date.
    pub fn from_jd(jd: f64) -> Option<Self> {
        Self::from_mjd(jd - 2_400_000.5)
    }

    /// The Julian year.
    pub fn year(self) -> f64 {
        self.0
    }

    /// Years from `other` to `self`.
    pub fn years_since(self, other: Epoch) -> f64 {
        self.0 - other.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conversions() {
        assert_eq!(Epoch::from_mjd(51_544.5), Some(Epoch::J2000));
        assert_eq!(Epoch::from_jd(2_451_545.0), Some(Epoch::J2000));
        let e = Epoch::from_mjd(51_544.5 + 365.25 * 26.0).unwrap();
        assert!((e.year() - 2026.0).abs() < 1e-12);
        assert!(Epoch::julian_year(f64::NAN).is_none());
        assert!(e.is_plausible() && Epoch::J2000.is_plausible());
        // A Julian Date in an MJD keyword, and an unset clock.
        assert!(!Epoch::from_mjd(2_459_175.2).unwrap().is_plausible());
        assert!(!Epoch::julian_year(1.0).unwrap().is_plausible());
    }
}
