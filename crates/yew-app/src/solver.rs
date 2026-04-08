use crate::examples;
use js_sys::Promise;
use serde::Serialize;
use wasm_bindgen::prelude::*;
use wasm_bindgen_futures::JsFuture;

// ---------------------------------------------------------------------------
// JS bridge: remSim.runInWorker(configJson, meshBytes) → Promise
// ---------------------------------------------------------------------------

#[wasm_bindgen]
extern "C" {
    /// Call `globalThis.remSim.runInWorker(configJson, meshBytes)`.
    /// Returns a JS Promise that resolves to the simulation result JsValue.
    #[wasm_bindgen(js_namespace = remSim, js_name = runInWorker)]
    fn run_in_worker_js(config_json: &str, mesh_bytes: &[u8]) -> Promise;
}

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Serialize)]
pub struct SParamPoint {
    pub freq_hz: f64,
    pub s11_db:  f64,
    pub s11_re:  f64,
    pub s11_im:  f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct SimResult {
    pub energy: f64,
    pub node_count: usize,
    pub max_e: Option<f64>,
    pub max_b: Option<f64>,
    pub frequencies_hz: Option<Vec<f64>>,
    /// Driven: S11 vs frequency
    pub s_params: Option<Vec<SParamPoint>>,
    /// Transient: port voltage vs time
    pub time_points: Option<Vec<f64>>,
    pub port_voltages: Option<Vec<f64>>,
    /// Eigenmode: Q-factors per mode (dielectric loss perturbation)
    pub q_factors: Option<Vec<f64>>,
    /// MoM/SBR: RCS pattern, (freq_hz, theta_deg, phi_deg, rcs_dbsm)
    pub rcs_data: Option<Vec<(f64, f64, f64, f64)>>,
}

#[derive(Clone, Debug)]
pub struct OutputArtifact {
    pub file_name: String,
    pub content: String,
}

#[derive(Clone, Debug)]
pub struct SimRun {
    pub summary: SimResult,
    pub artifacts: Vec<OutputArtifact>,
}

fn phi_csv(phi: &[f64]) -> String {
    let mut out = String::from("index,phi\n");
    for (i, v) in phi.iter().enumerate() {
        out.push_str(&format!("{},{}\n", i, v));
    }
    out
}

fn vec3_csv(name: &str, field: &[[f64; 3]]) -> String {
    let mut out = format!("index,{}_x,{}_y,{}_z,{}_norm\n", name, name, name, name);
    for (i, v) in field.iter().enumerate() {
        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        out.push_str(&format!("{},{},{},{},{}\n", i, v[0], v[1], v[2], norm));
    }
    out
}

fn s_params_csv(pts: &[SParamPoint]) -> String {
    let mut out = String::from("freq_hz,freq_ghz,s11_db,s11_re,s11_im\n");
    for p in pts {
        out.push_str(&format!(
            "{},{:.6},{:.4},{:.6},{:.6}\n",
            p.freq_hz, p.freq_hz / 1e9, p.s11_db, p.s11_re, p.s11_im
        ));
    }
    out
}

fn time_series_csv(times: &[f64], voltages: &[f64]) -> String {
    let mut out = String::from("time_s,time_ns,port_voltage\n");
    for (t, v) in times.iter().zip(voltages.iter()) {
        out.push_str(&format!("{},{:.4},{:.6}\n", t, t * 1e9, v));
    }
    out
}

fn rcs_csv(data: &[(f64, f64, f64, f64)]) -> String {
    let mut out = String::from("freq_hz,freq_ghz,theta_deg,phi_deg,rcs_dbsm\n");
    for &(f, th, ph, db) in data {
        out.push_str(&format!("{},{:.6},{:.1},{:.1},{:.4}\n", f, f / 1e9, th, ph, db));
    }
    out
}

pub async fn run_example(key: &str) -> Result<SimRun, String> {
    let example = examples::find_example(key)
        .ok_or_else(|| format!("Unknown example: {}", key))?;

    let mesh_bytes = examples::get_mesh_bytes(key);

    // Run the simulation in a dedicated Web Worker so the main thread stays
    // responsive. `run_in_worker_js` returns a JS Promise; we await it here.
    let js_value = JsFuture::from(run_in_worker_js(example.config_json, &mesh_bytes))
        .await
        .map_err(|e| format!("{} solve failed: {:?}", example.problem_type, e))?;

    let result: rem_wasm::SimulationResult = serde_wasm_bindgen::from_value(js_value)
        .map_err(|e| format!("Failed to decode simulation result: {}", e))?;

    let max_e = result.e_field.as_ref()
        .map(|field| {
            field.iter()
                .map(|v| (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt())
                .fold(0.0f64, f64::max)
        });
    let max_b = result.b_field.as_ref()
        .map(|field| {
            field.iter()
                .map(|v| (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt())
                .fold(0.0f64, f64::max)
        });

    // Convert s_params from wasm type to UI type
    let s_params: Option<Vec<SParamPoint>> = result.s_params.as_ref().map(|pts| {
        pts.iter().map(|p| SParamPoint {
            freq_hz: p.freq_hz,
            s11_db:  p.s11_db,
            s11_re:  p.s11_re,
            s11_im:  p.s11_im,
        }).collect()
    });

    // Flatten rcs_data: Vec<(freq, Vec<RcsPoint>)> → Vec<(freq, theta, phi, dbsm)>
    let rcs_data: Option<Vec<(f64, f64, f64, f64)>> = result.rcs_data.as_ref().map(|freqs| {
        freqs.iter().flat_map(|(f, pts)| {
            pts.iter().map(move |p| (*f, p.theta_deg, p.phi_deg, p.rcs_dbsm))
        }).collect()
    });

    let node_count = if !result.phi.is_empty() {
        result.phi.len()
    } else if let Some(pts) = &s_params {
        pts.len()
    } else if let Some(rd) = &rcs_data {
        rd.len()
    } else {
        0
    };

    let mut artifacts = vec![];

    if !result.phi.is_empty() {
        artifacts.push(OutputArtifact {
            file_name: "phi.csv".to_string(),
            content: phi_csv(&result.phi),
        });
    }

    if let Some(freqs) = &result.frequencies_hz {
        let q_factors = result.q_factors.as_deref().unwrap_or(&[]);
        let mut csv = String::from("mode,frequency_hz,frequency_ghz,q_factor\n");
        for (i, &f) in freqs.iter().enumerate() {
            let q = q_factors.get(i).copied().unwrap_or(f64::INFINITY);
            let q_str = if q.is_infinite() { "inf".to_string() } else { format!("{:.1}", q) };
            csv.push_str(&format!("{},{},{:.6},{}\n", i + 1, f, f / 1e9, q_str));
        }
        artifacts.push(OutputArtifact {
            file_name: "eigenfrequencies.csv".to_string(),
            content: csv,
        });
    }

    if let Some(field) = &result.e_field {
        artifacts.push(OutputArtifact {
            file_name: "e_field.csv".to_string(),
            content: vec3_csv("e", field),
        });
    }

    if let Some(field) = &result.b_field {
        artifacts.push(OutputArtifact {
            file_name: "b_field.csv".to_string(),
            content: vec3_csv("b", field),
        });
    }

    if let Some(pts) = &s_params {
        artifacts.push(OutputArtifact {
            file_name: "s_params.csv".to_string(),
            content: s_params_csv(pts),
        });
    }

    if let (Some(times), Some(voltages)) = (&result.time_points, &result.port_voltages) {
        artifacts.push(OutputArtifact {
            file_name: "port_voltage.csv".to_string(),
            content: time_series_csv(times, voltages),
        });
    }

    if let Some(rd) = &rcs_data {
        artifacts.push(OutputArtifact {
            file_name: "rcs.csv".to_string(),
            content: rcs_csv(rd),
        });
    }

    Ok(SimRun {
        summary: SimResult {
            energy: result.energy,
            node_count,
            max_e,
            max_b,
            frequencies_hz: result.frequencies_hz,
            s_params,
            time_points: result.time_points,
            port_voltages: result.port_voltages,
            q_factors: result.q_factors,
            rcs_data,
        },
        artifacts,
    })
}
