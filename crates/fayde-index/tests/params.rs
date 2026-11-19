use fayde_index::BuildParams;

fn check(edit: impl FnOnce(&mut BuildParams)) -> Result<(), String> {
    let mut p = BuildParams {
        fov_deg: 1.0,
        ..BuildParams::default()
    };
    edit(&mut p);
    p.validate().map_err(|e| e.to_string())
}

#[test]
fn defaults_and_the_recommended_set_pass() {
    check(|_| {}).unwrap();
    check(|p| {
        p.stars_per_fov = 23.7;
        p.pattern_stars_per_fov = 23.7;
        p.mag_limit = Some(15.5);
        p.min_edge_deg = 0.04;
    })
    .unwrap();
    check(|p| p.fov_deg = 0.2).unwrap();
}

#[test]
fn unbounded_lattices_are_rejected() {
    for v in [1e12, f64::INFINITY] {
        let e = check(|p| p.field_oversampling = v).unwrap_err();
        assert!(e.contains("field_oversampling"), "{e}");
    }
    let e = check(|p| {
        p.fov_deg = 0.01;
        p.min_edge_deg = 0.001;
    })
    .unwrap_err();
    assert!(e.contains("lattice fields"), "{e}");
}

#[test]
fn values_that_would_silently_disable_a_cut_are_rejected() {
    assert!(check(|p| p.mag_limit = Some(f32::NAN)).is_err());
    assert!(check(|p| p.max_motion_frac = f64::NAN).is_err());
    assert!(check(|p| p.motion_years = -1.0).is_err());
    assert!(check(|p| p.epoch_year = 99_999.0).is_err());
    assert!(check(|p| p.epoch_year = f64::NAN).is_err());
    assert!(check(|p| p.pattern_max_error = 1e-6).is_err());
}
