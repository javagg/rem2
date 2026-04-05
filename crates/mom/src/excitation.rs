//! Excitation vectors for MoM (plane wave, lumped port, etc.)

use crate::surface_mesh::SurfaceMesh;
use num_complex::Complex64;

/// Build the right-hand-side excitation vector for a +z-directed,
/// x-polarised incident plane wave:
///
/// E_inc(r) = x̂ exp(-jk·r·ẑ) = x̂ exp(-jkz)
///
/// For pulse basis: V[m] = -∫_Tm E_inc·dS (Galerkin testing with constant basis)
/// For RWG basis:   V[n] = -∫_{T⁺∪T⁻} f_n(r)·E_inc(r) dS
///
/// We use the face centroid as the quadrature point (1-point rule) which is
/// sufficient for the 1-point Gauss rule; for RWG we also use the centroid
/// of each half-support.
pub fn plane_wave_rhs(surf: &SurfaceMesh, k: f64, basis: &str) -> Vec<Complex64> {
    match basis.to_lowercase().as_str() {
        "pulse" => plane_wave_pulse(surf, k),
        _       => plane_wave_rwg(surf, k),
    }
}

/// Pulse-basis RHS: V[m] = -E_x(centroid_m) * area_m
fn plane_wave_pulse(surf: &SurfaceMesh, k: f64) -> Vec<Complex64> {
    surf.faces.iter().map(|face| {
        let z_c = face.centroid[2];
        let e_inc_x = Complex64::new(0.0, -k * z_c).exp(); // exp(-jkz)
        // V[m] = -∫_Tm x̂ · f_m dS; pulse: f_m = 1, testing = 1
        -e_inc_x * face.area
    }).collect()
}

/// RWG-basis RHS: V[n] = -∫_{T⁺∪T⁻} f_n(r)·E_inc(r) dS
/// Using centroid quadrature (1-point) on each triangle.
fn plane_wave_rwg(surf: &SurfaceMesh, k: f64) -> Vec<Complex64> {
    use crate::basis::rwg::generate_rwg_bases;
    let bases = generate_rwg_bases(surf);

    bases.iter().map(|b| {
        let mut val = Complex64::ZERO;
        for &(face_idx, in_plus) in &[(b.plus_face, true), (b.minus_face, false)] {
            let face = &surf.faces[face_idx];
            let r = &face.centroid;
            let e_inc_x = Complex64::new(0.0, -k * r[2]).exp();
            let fn_ = b.eval(r, surf, in_plus);
            // Only x-component of E_inc contributes
            val += e_inc_x * fn_[0] * face.area;
        }
        -val
    }).collect()
}
