//! Request values outside what the solver can honour are rejected before any search.

mod common;

use common::field;
use fayde_core::Epoch;
use fayde_solve::{solve, SolveError};
use fayde_wcs::Parity;

#[test]
fn implausible_epoch_and_sip_order_are_invalid_requests() {
    let (cat, f) = field((150.0, 2.0), 20.0, Parity::Negative, 3);
    assert!(solve(&f.request, &cat).is_ok());

    let mut req = f.request.clone();
    req.policy.epoch = Epoch::julian_year(8591.7);
    let e = solve(&req, &cat).unwrap_err();
    assert!(
        matches!(&e, SolveError::InvalidInput(m) if m.contains("epoch")),
        "{e}"
    );

    let mut req = f.request.clone();
    req.policy.max_sip_order = 10;
    let e = solve(&req, &cat).unwrap_err();
    assert!(
        matches!(&e, SolveError::InvalidInput(m) if m.contains("max_sip_order")),
        "{e}"
    );

    let mut req = f.request.clone();
    req.policy.max_sip_order = 9;
    req.policy.epoch = Epoch::julian_year(1890.0);
    assert!(solve(&req, &cat).is_ok());
}
