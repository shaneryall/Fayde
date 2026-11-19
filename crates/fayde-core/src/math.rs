#[cfg(all(target_arch = "x86", not(target_feature = "sse2")))]
compile_error!(
    "Fayde needs SSE2 on 32-bit x86: x87 arithmetic is not reproducible across machines"
);

#[doc(hidden)]
pub fn round_decimal(value: f64, places: usize) -> f64 {
    format!("{value:.places$}").parse().unwrap_or(value)
}

pub trait Portable: Copy {
    fn psin(self) -> Self;
    fn pcos(self) -> Self;
    fn psin_cos(self) -> (Self, Self);
    fn ptan(self) -> Self;
    fn pasin(self) -> Self;
    fn pacos(self) -> Self;
    fn patan(self) -> Self;
    fn patan2(self, x: Self) -> Self;
    fn pexp(self) -> Self;
    fn pln(self) -> Self;
    fn plog10(self) -> Self;
    fn ppowf(self, e: Self) -> Self;
    fn ppowi(self, n: i32) -> Self;
    fn phypot(self, other: Self) -> Self;
    fn perfc(self) -> Self;
}

macro_rules! portable {
    ($t:ty, $sin:ident, $cos:ident, $sincos:ident, $tan:ident, $asin:ident, $acos:ident,
     $atan:ident, $atan2:ident, $exp:ident, $ln:ident, $log10:ident, $pow:ident, $hypot:ident,
     $erfc:ident) => {
        impl Portable for $t {
            #[inline]
            fn psin(self) -> Self {
                libm::$sin(self)
            }
            #[inline]
            fn pcos(self) -> Self {
                libm::$cos(self)
            }
            #[inline]
            fn psin_cos(self) -> (Self, Self) {
                libm::$sincos(self)
            }
            #[inline]
            fn ptan(self) -> Self {
                libm::$tan(self)
            }
            #[inline]
            fn pasin(self) -> Self {
                libm::$asin(self)
            }
            #[inline]
            fn pacos(self) -> Self {
                libm::$acos(self)
            }
            #[inline]
            fn patan(self) -> Self {
                libm::$atan(self)
            }
            #[inline]
            fn patan2(self, x: Self) -> Self {
                libm::$atan2(self, x)
            }
            #[inline]
            fn pexp(self) -> Self {
                libm::$exp(self)
            }
            #[inline]
            fn pln(self) -> Self {
                libm::$ln(self)
            }
            #[inline]
            fn plog10(self) -> Self {
                libm::$log10(self)
            }
            #[inline]
            fn ppowf(self, e: Self) -> Self {
                libm::$pow(self, e)
            }
            #[inline]
            fn ppowi(self, n: i32) -> Self {
                let mut base = self;
                let mut e = n.unsigned_abs();
                let mut acc: Self = 1.0;
                while e > 0 {
                    if e & 1 == 1 {
                        acc *= base;
                    }
                    e >>= 1;
                    if e > 0 {
                        base *= base;
                    }
                }
                if n >= 0 {
                    acc
                } else if acc.is_finite() {
                    1.0 / acc
                } else {
                    (1.0 / self).ppowi(n.saturating_neg().max(1))
                        * if n == i32::MIN { 1.0 / self } else { 1.0 }
                }
            }
            #[inline]
            fn phypot(self, other: Self) -> Self {
                libm::$hypot(self, other)
            }
            #[inline]
            fn perfc(self) -> Self {
                libm::$erfc(self)
            }
        }
    };
}

portable!(f64, sin, cos, sincos, tan, asin, acos, atan, atan2, exp, log, log10, pow, hypot, erfc);
portable!(
    f32, sinf, cosf, sincosf, tanf, asinf, acosf, atanf, atan2f, expf, logf, log10f, powf, hypotf,
    erfcf
);

#[cfg(test)]
mod tests {
    use super::Portable;

    #[test]
    fn powi_matches_repeated_multiplication() {
        let x = 1.000_000_1_f64;
        assert_eq!(x.ppowi(0), 1.0);
        assert_eq!(x.ppowi(1), x);
        assert_eq!(x.ppowi(2), x * x);
        assert_eq!(x.ppowi(3), x * (x * x));
        assert_eq!(2f64.ppowi(-2), 0.25);
        assert_eq!(3f32.ppowi(4), 81.0);
        assert_eq!(2f64.ppowi(-1024), f64::from_bits(1 << 50));
        assert_eq!(2f32.ppowi(-128), f32::from_bits(1 << 21));
        assert_eq!(2f64.ppowi(i32::MIN), 0.0);
        assert_eq!(0.5f64.ppowi(i32::MIN), f64::INFINITY);
    }

    #[test]
    fn results_are_pinned_to_the_bit() {
        let x = 0.7_f64;
        let got = [
            x.psin().to_bits(),
            x.pcos().to_bits(),
            x.patan2(1.3).to_bits(),
            x.pln().to_bits(),
            x.pexp().to_bits(),
            x.phypot(1.3).to_bits(),
        ];
        assert_eq!(got, PINNED, "{got:#x?}");
    }

    const PINNED: [u64; 6] = [
        0x3fe4_9d6e_6946_19b8,
        0x3fe8_7996_529f_9d93,
        0x3fdf_9cbc_4269_ab2a,
        0xbfd6_d3c3_24e1_3f50,
        0x4000_1c2a_6126_8987,
        0x3ff7_9fab_e917_e2f5,
    ];
}
