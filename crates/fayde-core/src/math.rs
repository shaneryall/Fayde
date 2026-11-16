//! Portable floating-point maths: the same bits on every OS and CPU.
//!
//! `f64::sin`, `cos`, `atan2`, `ln`, `exp`, `hypot`, `powf` and friends call the platform's
//! maths library (Apple's libm, glibc, the MSVC CRT), which are not correctly rounded and differ
//! in the last bit. `powi`'s precision is documented as unspecified and differs between debug and
//! release builds. Index files are hashed and that hash is part of the solver's identity, so
//! Fayde must build the same bytes from the same catalogue on macOS, Linux and Windows.
//!
//! [`Portable`] routes every such call through the `libm` crate, a pure-Rust port of musl's libm
//! (itself from FreeBSD msun / Sun fdlibm), built without its `arch` feature so no
//! platform-specific instruction is used. `+ − × ÷` and `sqrt` are correctly rounded by IEEE 754
//! and need no wrapper; Rust never contracts `a * b + c` into a fused multiply-add on its own.
//!
//! The std methods are banned workspace-wide by `clippy.toml` (`disallowed-methods`).
//!
//! Supported targets: any with IEEE 754 binary64 arithmetic in SSE2-class registers, i.e. every
//! 64-bit target (x86_64, aarch64) and 32-bit x86 with SSE2 (the standard `i686` targets). On
//! x86 without SSE2 (`i586`) plain `+ − × ÷` use the x87 FPU's extended precision and `libm`
//! must use x87 instructions, so results cannot match other machines: the build refuses it.

#[cfg(all(target_arch = "x86", not(target_feature = "sse2")))]
compile_error!(
    "Fayde needs SSE2 on 32-bit x86: x87 arithmetic is not reproducible across machines"
);

/// Transcendental functions with platform-independent results. The `p` prefix marks the portable
/// version of the std method of the same name.
pub trait Portable: Copy {
    /// Sine.
    fn psin(self) -> Self;
    /// Cosine.
    fn pcos(self) -> Self;
    /// `(sin, cos)`.
    fn psin_cos(self) -> (Self, Self);
    /// Tangent.
    fn ptan(self) -> Self;
    /// Arcsine.
    fn pasin(self) -> Self;
    /// Arccosine.
    fn pacos(self) -> Self;
    /// Arctangent.
    fn patan(self) -> Self;
    /// Four-quadrant arctangent of `self / x`.
    fn patan2(self, x: Self) -> Self;
    /// `e^self`.
    fn pexp(self) -> Self;
    /// Natural logarithm.
    fn pln(self) -> Self;
    /// Base-10 logarithm.
    fn plog10(self) -> Self;
    /// `self^e`.
    fn ppowf(self, e: Self) -> Self;
    /// `self^n` by binary exponentiation in a fixed order.
    fn ppowi(self, n: i32) -> Self;
    /// `sqrt(self² + other²)` without undue overflow.
    fn phypot(self, other: Self) -> Self;
    /// Complementary error function, `1 − erf(self)`.
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
                    // |self|^|n| overflowed, but self^n may still be a (subnormal) number:
                    // exponentiate the reciprocal instead.
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
        // Overflowing positive power, representable negative power.
        assert_eq!(2f64.ppowi(-1024), f64::from_bits(1 << 50));
        assert_eq!(2f32.ppowi(-128), f32::from_bits(1 << 21));
        assert_eq!(2f64.ppowi(i32::MIN), 0.0);
        assert_eq!(0.5f64.ppowi(i32::MIN), f64::INFINITY);
    }

    /// Fixed bit patterns: any platform that disagrees is not running the portable code.
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
