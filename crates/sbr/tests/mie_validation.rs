//! Integration test: SBR+ PEC sphere RCS vs Mie series.
//!
//! Validates that the full SBR+ pipeline (mesh → BVH → ray tracing →
//! PO integral → RCS) produces results consistent with the Mie analytical
//! solution at high frequency (kα ≈ 31.4 >> 1).
//!
//! Tolerance: monostatic (θ=180°) RCS within 3 dB of Mie.

use rem_config::{load_config_from_str, ConfigFormat};
use rem_mesh::{gen::pec_sphere_msh, load_mesh_from_bytes};
use rem_parallel::NoComm;
use rem_sbr::run_with_mesh;

/// Compute monostatic PEC sphere RCS via Mie series.
/// Returns σ [m²] at θ=180° (backscatter).
fn mie_monostatic(radius: f64, freq: f64) -> f64 {
    use std::f64::consts::PI;
    let k = 2.0 * PI * freq / 2.997_924_58e8;
    let sigma_vec = rem_mom::mie::pec_sphere_rcs(radius, k, &[180.0_f64], None);
    sigma_vec[0]
}

#[test]
fn sbr_sphere_rcs_vs_mie() {
    // ── Geometry ────────────────────────────────────────────────────────────
    // a = 0.5 m, f = 1 GHz → ka ≈ 10.5  (optical regime, PO accurate)
    // λ = 0.3 m, λ/4 = 0.075 m. Sphere arc (pole to equator) = πa/2 = 0.785 m.
    // Need rings > 0.785/0.075 ≈ 11 → 24 rings is adequate (2× oversampled).
    let radius  = 0.5_f64;
    let freq    = 1.0e9_f64;

    // ── Generate surface mesh ───────────────────────────────────────────────
    // 24 lat × 48 lon → 2304 faces; well resolved at 1 GHz (ka≈10.5)
    let msh_str = pec_sphere_msh(radius, 24, 48, 1);

    // ── Minimal Palace config ───────────────────────────────────────────────
    let out_dir = std::env::temp_dir().join("sbr_sphere_test");
    std::fs::create_dir_all(out_dir.join("postpro")).unwrap();
    let out_str = out_dir.to_string_lossy().replace('\\', "/");

    let config_json = format!(r#"{{
        "Problem": {{ "Type": "SBR", "Output": "{out}" }},
        "Model": {{ "Mesh": "dummy.msh", "L0": 1.0 }},
        "Boundaries": {{ "PEC": {{ "Attributes": [1] }} }},
        "Solver": {{
            "SBR": {{
                "FreqMin": {f}, "FreqMax": {f}, "FreqStep": 0.0,
                "RayDensity": 3000.0,
                "MaxBounces": 2,
                "WeightThresh": 1e-4,
                "TargetType": "PEC",
                "ThetaInc": 0.0, "PhiInc": 0.0,
                "Polarization": "theta"
            }}
        }},
        "Postprocessing": {{
            "RCS": {{
                "ThetaDeg": [180.0],
                "PhiDeg":   [0.0]
            }}
        }}
    }}"#, out = out_str, f = freq);

    let config = load_config_from_str(&config_json, ConfigFormat::Json)
        .expect("config parse failed");

    let comm = NoComm;
    let mesh = load_mesh_from_bytes(&config, msh_str.as_bytes(), &comm)
        .expect("mesh load failed");

    let sbr_cfg = config.solver.sbr.as_ref().unwrap();
    run_with_mesh(&config, sbr_cfg, &mesh).expect("SBR+ solve failed");

    // ── Read RCS from CSV ───────────────────────────────────────────────────
    let csv_path = out_dir.join("postpro").join("rcs_sbr.csv");
    let csv = std::fs::read_to_string(&csv_path).expect("rcs_sbr.csv not found");

    // Parse last line (θ=180, φ=0)
    let rcs_dbsm: f64 = csv.lines()
        .filter(|l| !l.starts_with("Freq"))
        .last()
        .and_then(|l| l.split(',').nth(3))
        .and_then(|v| v.trim().parse().ok())
        .expect("could not parse RCS from CSV");

    // ── Compare to Mie ─────────────────────────────────────────────────────
    let mie_m2   = mie_monostatic(radius, freq);
    let mie_dbsm = 10.0 * mie_m2.log10();

    let error_db = (rcs_dbsm - mie_dbsm).abs();
    println!("SBR+ monostatic RCS = {:.2} dBsm", rcs_dbsm);
    println!("Mie monostatic RCS  = {:.2} dBsm", mie_dbsm);
    println!("Error               = {:.2} dB", error_db);

    // Clean up
    let _ = std::fs::remove_dir_all(&out_dir);

    assert!(
        error_db < 3.0,
        "SBR+ monostatic RCS {:.2} dBsm deviates from Mie {:.2} dBsm by {:.2} dB (limit: 3 dB)",
        rcs_dbsm, mie_dbsm, error_db
    );
}
