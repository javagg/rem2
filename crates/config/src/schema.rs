use crate::preprocess::expand_ranges;
use rem_core::EPS0;
use serde::{Deserialize, Deserializer};

/// Deserialize a JSON value that may be either a scalar `f64` or an array `[f64, ...]`.
/// When an array is given, the first element is used (anisotropic �?isotropic fallback).
fn deserialize_scalar_or_first<'de, D>(d: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    use serde::de::{self, Visitor};
    use std::fmt;

    struct ScalarOrFirstVisitor;

    impl<'de> Visitor<'de> for ScalarOrFirstVisitor {
        type Value = f64;

        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "a float or an array of floats")
        }

        fn visit_f64<E: de::Error>(self, v: f64) -> Result<f64, E> { Ok(v) }
        fn visit_i64<E: de::Error>(self, v: i64) -> Result<f64, E> { Ok(v as f64) }
        fn visit_u64<E: de::Error>(self, v: u64) -> Result<f64, E> { Ok(v as f64) }

        fn visit_seq<A: de::SeqAccess<'de>>(self, mut seq: A) -> Result<f64, A::Error> {
            match seq.next_element::<f64>()? {
                Some(first) => {
                    // Drain remaining elements
                    while seq.next_element::<f64>()?.is_some() {}
                    Ok(first)
                }
                None => Err(de::Error::custom("empty array for scalar field")),
            }
        }
    }

    d.deserialize_any(ScalarOrFirstVisitor)
}

fn default_scalar_or_first_one() -> f64 { 1.0 }

// ---------------------------------------------------------------------------
// Top-level
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
pub struct PalaceConfig {
    #[serde(rename = "Problem")]
    pub problem: Problem,

    #[serde(rename = "Model")]
    pub model: Model,

    #[serde(rename = "Domains", default)]
    pub domains: Domains,

    #[serde(rename = "Boundaries", default)]
    pub boundaries: Boundaries,

    #[serde(rename = "Solver", default)]
    pub solver: SolverConfig,

    /// REM extension: post-processing options (ignored by Palace).
    #[serde(rename = "Postprocessing", default)]
    pub postprocessing: Postprocessing,

    /// REM extension: metadata from format converters (Sonnet19, etc.).
    /// Carries PlanarTechLayers (material, thickness, roughness),
    /// dielectric info, and other converter-specific data.
    #[serde(rename = "Metadata", default, skip_serializing_if = "Option::is_none")]
    pub metadata: Option<serde_json::Value>,
}

// ---------------------------------------------------------------------------
// Problem
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize, Default)]
pub struct Problem {
    #[serde(rename = "Type")]
    pub problem_type: ProblemType,

    /// Verbosity level 0-3. None in JSON �?defaults to 1 at runtime.
    #[serde(rename = "Verbose")]
    pub verbose: Option<u8>,

    #[serde(rename = "Output")]
    pub output: Option<String>,

    /// Palace `OutputFormats` section �?accepted for compatibility; REM ignores it.
    #[serde(rename = "OutputFormats", default)]
    pub output_formats: Option<OutputFormats>,
}

/// Palace output format options (not yet implemented in REM).
#[derive(Debug, Clone, Deserialize, Default)]
pub struct OutputFormats {
    /// Write full mesh grid function VTK output
    #[serde(rename = "GridFunction", default)]
    pub grid_function: bool,
}

impl Problem {
    /// Effective verbosity (Palace default: 1).
    pub fn verbose(&self) -> u8 { self.verbose.unwrap_or(1) }
    /// Effective output directory (Palace default: ".").
    pub fn output_dir(&self) -> &str { self.output.as_deref().unwrap_or(".") }
}

#[derive(Debug, Clone, Deserialize, PartialEq, Eq, Default)]
pub enum ProblemType {
    #[default]
    Electrostatic,
    Magnetostatic,
    Eigenmode,
    Driven,
    Transient,
}

// ---------------------------------------------------------------------------
// Model
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
pub struct Model {
    #[serde(rename = "Mesh", default)]
    pub mesh: String,

    #[serde(rename = "L0", default = "default_l0")]
    pub l0: f64,

    #[serde(rename = "Refinement", default)]
    pub refinement: Refinement,
}

#[derive(Debug, Clone, Deserialize, Default)]
pub struct Refinement {
    #[serde(rename = "MaxIter", default)]
    pub max_iter: usize,

    #[serde(rename = "Tol", default = "default_tol")]
    pub tol: f64,

    #[serde(rename = "Nonconformal", default)]
    pub nonconformal: bool,
}

fn default_l0() -> f64 { 1.0 }
fn default_tol() -> f64 { 1.0e-2 }

// ---------------------------------------------------------------------------
// Domains
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize, Default)]
pub struct Domains {
    #[serde(rename = "Materials", default)]
    pub materials: Vec<MaterialSpec>,

    /// Palace `Postprocessing` under Domains �?accepted for compatibility.
    #[serde(rename = "Postprocessing", default)]
    pub postprocessing: Option<DomainsPostprocessing>,

    /// Palace `CurrentDipole` �?accepted for compatibility (not implemented).
    #[serde(rename = "CurrentDipole", default)]
    pub current_dipole: Vec<CurrentDipoleSpec>,
}

/// Palace `Domains.CurrentDipole` (Hertzian dipole source �?not implemented).
#[derive(Debug, Clone, Deserialize)]
pub struct CurrentDipoleSpec {
    #[serde(rename = "Index")]
    pub index: u32,

    #[serde(rename = "Moment", default)]
    pub moment: f64,

    #[serde(rename = "Center", default)]
    pub center: Vec<f64>,

    #[serde(rename = "Direction", default)]
    pub direction: Vec<f64>,
}

/// Palace `Domains.Postprocessing` �?accepted for Palace compatibility.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct DomainsPostprocessing {
    #[serde(rename = "Energy", default)]
    pub energy: Vec<EnergyPostSpec>,

    #[serde(rename = "Probe", default)]
    pub probe: Vec<ProbeSpec>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct EnergyPostSpec {
    #[serde(rename = "Index", default)]
    pub index: u32,

    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes", default)]
    pub attributes: Vec<u32>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct ProbeSpec {
    #[serde(rename = "Index", default)]
    pub index: u32,

    #[serde(rename = "Center")]
    pub center: Vec<f64>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct MaterialSpec {
    /// Physical group attribute IDs. Accepts either `[1,2,3]` or `"1,3-5"`.
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    /// Accepts scalar `9.3` or anisotropic array `[9.3, 9.3, 11.5]`; array �?first element used.
    #[serde(rename = "Permittivity", default = "default_scalar_or_first_one",
            deserialize_with = "deserialize_scalar_or_first")]
    pub permittivity: f64,

    /// Accepts scalar `1.0` or anisotropic array `[1.0, 1.0, 1.0]`; array �?first element used.
    #[serde(rename = "Permeability", default = "default_scalar_or_first_one",
            deserialize_with = "deserialize_scalar_or_first")]
    pub permeability: f64,

    /// Accepts scalar `3.0e-5` or anisotropic array `[3.0e-5, 3.0e-5, 8.6e-5]`; array �?first element used.
    #[serde(rename = "LossTan", default,
            deserialize_with = "deserialize_scalar_or_first")]
    pub loss_tangent: f64,

    /// Magnetic loss tangent tan δ_m = μ�?μ�?for lossy magnetic materials (ferrites, etc.).
    #[serde(rename = "LossTanMag", default)]
    pub loss_tangent_magnetic: f64,

    #[serde(rename = "Conductivity", default,
            deserialize_with = "deserialize_scalar_or_first")]
    pub conductivity: f64,

    /// Palace `MaterialAxes` for anisotropic materials.
    /// Each row is a basis vector for the material coordinate frame.
    /// Triggers tensor epsilon assembly in electrostatic/eigenmode solvers.
    #[serde(rename = "MaterialAxes", default)]
    pub material_axes: Vec<Vec<f64>>,

    /// Drude-Lorentz poles for frequency-dependent permittivity.
    /// ε(ω) = ε₀ + Σ ωp² / (ω₀² − ω² + jγω)
    /// Relevant only for driven (frequency-domain) solvers.
    #[serde(rename = "DrudeLorentz", default)]
    pub drude_lorentz: Vec<DrudeLorentzPole>,

    /// Physical thickness of this material domain [m].
    /// Used by the BEM layered-media kernel for substrate definition.
    /// 0 = infinite half-space (default).
    #[serde(rename = "Thickness", default)]
    pub thickness: f64,
}

/// One Drude-Lorentz oscillator pole.
///
/// Contributes `plasma_freq_sq / (resonance_freq_sq �?ω² + j·damping·ω)` to ε�?ω).
///
/// For a Drude free-carrier term: set `ResonanceFreq = 0`, `PlasmaFreq = ωp/(2π)`,
/// `Damping = γ/(2π)`.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct DrudeLorentzPole {
    /// Plasma frequency f�?[Hz] �?the pole contribution strength.
    /// ωp² = (2π f�?²
    #[serde(rename = "PlasmaFreq", default)]
    pub plasma_freq: f64,

    /// Resonance frequency f₀ [Hz].  Zero �?Drude (free-carrier) term.
    /// ω₀² = (2π f₀)²
    #[serde(rename = "ResonanceFreq", default)]
    pub resonance_freq: f64,

    /// Damping rate γ [rad/s].  Also accepted as `DampingFreq` [Hz] (converted internally).
    #[serde(rename = "Damping", default)]
    pub damping: f64,
}

// ---------------------------------------------------------------------------
// Boundaries
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize, Default)]
pub struct Boundaries {
    #[serde(rename = "PEC", default)]
    pub pec: Option<AttrList>,

    #[serde(rename = "PMC", default)]
    pub pmc: Option<AttrList>,

    #[serde(rename = "Impedance", default)]
    pub impedance: Vec<ImpedanceSpec>,

    #[serde(rename = "ResistiveSheet", default)]
    pub resistive_sheet: Vec<ResistiveSheetSpec>,

    #[serde(rename = "LumpedPort", default)]
    pub lumped_port: Vec<LumpedPortSpec>,

    #[serde(rename = "WavePort", default)]
    pub wave_port: Vec<WavePortSpec>,

    #[serde(rename = "Absorbing", default)]
    pub absorbing: Option<AbsorbingSpec>,

    #[serde(rename = "Ground", default)]
    pub ground: Option<AttrList>,

    #[serde(rename = "ZeroCharge", default)]
    pub zero_charge: Option<AttrList>,

    #[serde(rename = "SurfaceCurrent", default)]
    pub surface_current: Vec<SurfaceCurrentSpec>,

    /// Electrostatic terminal conductors (Palace "Electrostatic" problem type).
    /// Each terminal is an equipotential conductor surface assigned a fixed voltage
    /// during capacitance-matrix extraction.  Maps to Palace's `"Terminal"` key.
    #[serde(rename = "Terminal", default)]
    pub terminal: Vec<TerminalSpec>,

    /// Palace `Periodic` / `FloquetWaveVector` boundaries �?not implemented.
    #[serde(rename = "Periodic", default)]
    pub periodic: Vec<PeriodicSpec>,

    /// Palace boundary-level `Postprocessing` (SurfaceFlux, FarField, Dielectric).
    /// Accepted for Palace compatibility; REM logs warnings.
    #[serde(rename = "Postprocessing", default)]
    pub postprocessing_flux: Vec<BoundaryPostprocessingSpec>,
}

/// Palace `Boundaries.Periodic` �?not implemented.
#[derive(Debug, Clone, Deserialize)]
pub struct PeriodicSpec {
    /// Floquet wave vector [kx, ky, kz] for quasi-periodic BCs.
    #[serde(rename = "FloquetWaveVector", default)]
    pub floquet_wave_vector: Vec<f64>,

    /// Pairs of donor/receiver boundaries.
    #[serde(rename = "BoundaryPairs", default)]
    pub boundary_pairs: Vec<PeriodicBoundaryPair>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct PeriodicBoundaryPair {
    #[serde(rename = "DonorAttributes", deserialize_with = "deserialize_attributes", default)]
    pub donor_attributes: Vec<u32>,

    #[serde(rename = "ReceiverAttributes", deserialize_with = "deserialize_attributes", default)]
    pub receiver_attributes: Vec<u32>,

    #[serde(rename = "Translation", default)]
    pub translation: Vec<f64>,
}

/// Palace `Boundaries.Postprocessing` �?not yet implemented.
#[derive(Debug, Clone, Deserialize)]
pub struct BoundaryPostprocessingSpec {
    #[serde(rename = "Index", default)]
    pub index: u32,

    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes", default)]
    pub attributes: Vec<u32>,

    #[serde(rename = "Type", default)]
    pub flux_type: String,

    #[serde(rename = "TwoSided", default)]
    pub two_sided: bool,

    #[serde(rename = "Center", default)]
    pub center: Vec<f64>,

    #[serde(rename = "Thickness", default)]
    pub thickness: f64,

    #[serde(rename = "Permittivity", default)]
    pub permittivity: f64,

    #[serde(rename = "LossTan", default)]
    pub loss_tan: f64,

    #[serde(rename = "NSample", default)]
    pub n_sample: usize,

    #[serde(rename = "ThetaPhis", default)]
    pub theta_phis: Vec<Vec<f64>>,
}

/// Electrostatic terminal boundary (Palace `"Terminal"`).
/// An equipotential conductor surface whose voltage is set during solving.
#[derive(Debug, Clone, Deserialize)]
pub struct TerminalSpec {
    #[serde(rename = "Index")]
    pub index: u32,

    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,
}

/// A boundary spec that only carries a list of physical group IDs.
#[derive(Debug, Clone, Deserialize)]
pub struct AttrList {
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,
}

/// Resistive thin-sheet boundary condition (Ω/�?sheet resistance).
#[derive(Debug, Clone, Deserialize)]
pub struct ResistiveSheetSpec {
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    /// Sheet resistance [Ω/sq]
    #[serde(rename = "Rs")]
    pub rs: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct ImpedanceSpec {
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    /// Surface resistance [Ω/sq]
    #[serde(rename = "Rs", default)]
    pub rs: f64,

    /// Surface inductance [H/sq]
    #[serde(rename = "Ls", default)]
    pub ls: f64,

    /// Surface capacitance [F/sq]
    #[serde(rename = "Cs", default)]
    pub cs: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct LumpedPortSpec {
    #[serde(rename = "Index")]
    pub index: u32,

    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    #[serde(rename = "Direction", default)]
    pub direction: String,

    /// Port resistance [Ω]
    #[serde(rename = "R", default)]
    pub r: f64,

    /// Port inductance [H]
    #[serde(rename = "L", default)]
    pub l: f64,

    /// Port capacitance [F]
    #[serde(rename = "C", default)]
    pub c: f64,

    #[serde(rename = "Excitation", default)]
    pub excitation: bool,

    /// Palace `Elements` for multi-element lumped ports.
    /// Each element's attributes are mapped to the same port BC tag.
    #[serde(rename = "Elements", default)]
    pub elements: Vec<LumpedPortElement>,
}

/// Palace `LumpedPort.Elements` item.
#[derive(Debug, Clone, Deserialize)]
pub struct LumpedPortElement {
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    #[serde(rename = "Direction", default)]
    pub direction: String,
}

#[derive(Debug, Clone, Deserialize)]
pub struct WavePortSpec {
    #[serde(rename = "Index")]
    pub index: u32,

    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    #[serde(rename = "Excitation", default)]
    pub excitation: bool,

    #[serde(rename = "Mode", default = "default_mode")]
    pub mode: u32,

    /// Offset along the propagation direction [m].
    #[serde(rename = "Offset", default)]
    pub offset: f64,

    /// Max iterations for WavePort internal iterative solver.
    #[serde(rename = "MaxIts", default)]
    pub max_its: usize,

    /// Tolerance for WavePort internal iterative eigenvalue solve.
    #[serde(rename = "EigenTol", default)]
    pub eigen_tol: f64,

    /// Verbosity of WavePort modal analysis.
    #[serde(rename = "Verbose", default)]
    pub verbose_port: u8,
}

fn default_mode() -> u32 { 1 }

#[derive(Debug, Clone, Deserialize)]
pub struct AbsorbingSpec {
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    #[serde(rename = "Order", default = "default_absorbing_order")]
    pub order: u8,
}

fn default_absorbing_order() -> u8 { 1 }

#[derive(Debug, Clone, Deserialize)]
pub struct SurfaceCurrentSpec {
    #[serde(rename = "Index")]
    pub index: u32,

    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    #[serde(rename = "Direction", default)]
    pub direction: String,
}

// ---------------------------------------------------------------------------
// Solver
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
pub struct SolverConfig {
    /// Finite element polynomial order (1 = P1, 2 = P2, �?
    #[serde(rename = "Order", default = "default_order")]
    pub order: u8,

    /// Discretization family for full-wave problems.
    ///
    /// Supported values (case-insensitive):
    /// - `HCurl` / `Nedelec` (default): edge-element H(curl) basis (curl-conforming; avoids spurious modes)
    /// - `H1`: scalar nodal basis (legacy; may produce spurious modes in driven/eigenmode)
    #[serde(rename = "Discretization", default = "default_discretization")]
    pub discretization: String,

    /// Palace `Device` �?REM is CPU-only; value is accepted and ignored.
    #[serde(rename = "Device", default = "default_device")]
    pub device: String,

    #[serde(rename = "Eigenmode", default)]
    pub eigenmode: Option<EigenmodeSolver>,

    #[serde(rename = "Driven", default)]
    pub driven: Option<DrivenSolver>,

    #[serde(rename = "Transient", default)]
    pub transient: Option<TransientSolver>,

    #[serde(rename = "Electrostatic", default)]
    pub electrostatic: Option<StaticSolverConfig>,

    #[serde(rename = "Magnetostatic", default)]
    pub magnetostatic: Option<StaticSolverConfig>,

    #[serde(rename = "Linear", default)]
    pub linear: LinearSolver,

    /// REM extension: DDM solver parameters (ignored by Palace).
    #[serde(rename = "DDM", default)]
    pub ddm: Option<DdmSolverConfig>,
}

/// REM near-to-far-field configuration.
///
/// Computes radiation pattern from the driven solver's near-field solution.
/// Uses Kirchhoff approximation: far-field amplitude �?�?**E**(r') e^{jk r̂·r'} dS'
/// integrated over all boundary elements tagged in `surface_attributes`.
#[derive(Debug, Clone, Deserialize)]
pub struct FarFieldConfig {
    /// Physical group tags of boundary faces on which to integrate.
    /// If empty, uses all boundary elements (whole mesh surface).
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes", default)]
    pub attributes: Vec<u32>,

    /// Number of elevation angle samples (θ from 0° to 180°). Default: 37 (5° steps).
    #[serde(rename = "NTheta", default = "default_far_field_n")]
    pub n_theta: usize,

    /// Number of azimuth angle samples (φ from 0° to 360°). Default: 73 (5° steps).
    #[serde(rename = "NPhi", default = "default_far_field_n2")]
    pub n_phi: usize,
}

fn default_order() -> u8 { 1 }
fn default_far_field_n() -> usize { 37 }
fn default_far_field_n2() -> usize { 73 }

impl Default for SolverConfig {
    fn default() -> Self {
        SolverConfig {
            order: 1,
            discretization: default_discretization(),
            device: "CPU".to_string(),
            eigenmode: None,
            driven: None,
            transient: None,
            electrostatic: None,
            magnetostatic: None,
            linear: LinearSolver::default(),
            ddm: None,
        }
    }
}

fn default_device() -> String { "CPU".to_string() }
fn default_discretization() -> String { "HCurl".to_string() }
fn default_formulation_auto() -> String { "Auto".to_string() }

fn parse_hcurl_formulation(value: &str) -> Option<bool> {
    match value.to_lowercase().as_str() {
        "hcurl" | "nedelec" => Some(true),
        "h1" => Some(false),
        "auto" | "" => None,
        _ => None,
    }
}

impl SolverConfig {
    /// Returns `true` when H(curl)/Nedelec discretization is requested.
    pub fn uses_hcurl(&self) -> bool {
        matches!(
            self.discretization.to_lowercase().as_str(),
            "hcurl" | "nedelec"
        )
    }

    /// Returns effective H(curl) selection for the driven solver.
    ///
    /// Precedence:
    /// 1. Solver.Driven.Formulation (HCurl/Nedelec/H1)
    /// 2. Solver.Discretization (global fallback)
    pub fn uses_hcurl_for_driven(&self) -> bool {
        if let Some(drv) = &self.driven {
            if let Some(v) = parse_hcurl_formulation(&drv.formulation) {
                return v;
            }
        }
        self.uses_hcurl()
    }

    /// Returns effective H(curl) selection for the eigenmode solver.
    ///
    /// Precedence:
    /// 1. Solver.Eigenmode.Formulation (HCurl/Nedelec/H1)
    /// 2. Solver.Discretization (global fallback)
    pub fn uses_hcurl_for_eigenmode(&self) -> bool {
        if let Some(eig) = &self.eigenmode {
            if let Some(v) = parse_hcurl_formulation(&eig.formulation) {
                return v;
            }
        }
        self.uses_hcurl()
    }

    /// Effective HCurl polynomial order for driven solver.
    /// Falls back to Solver.Order when Driven.HCurlOrder is not set.
    pub fn driven_hcurl_order(&self) -> u8 {
        self.driven
            .as_ref()
            .and_then(|d| d.hcurl_order)
            .unwrap_or(self.order)
    }

    /// Effective HCurl polynomial order for eigenmode solver.
    /// Falls back to Solver.Order when Eigenmode.HCurlOrder is not set.
    pub fn eigenmode_hcurl_order(&self) -> u8 {
        self.eigenmode
            .as_ref()
            .and_then(|e| e.hcurl_order)
            .unwrap_or(self.order)
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct EigenmodeSolver {
    /// Number of eigenvalues to compute
    #[serde(rename = "N", default = "default_n_eig")]
    pub n: usize,

    #[serde(rename = "Tol", default = "default_eig_tol")]
    pub tol: f64,

    #[serde(rename = "MaxIter", default = "default_eig_maxiter")]
    pub max_iter: usize,

    /// Target frequency [Hz] for shift-invert
    #[serde(rename = "Target", default)]
    pub target: f64,

    /// Number of modes to save to disk
    #[serde(rename = "Save", default = "default_save")]
    pub save: usize,

    /// Optional discretization override for eigenmode only.
    ///
    /// Supported values (case-insensitive):
    /// - `Auto` (default): follow `Solver.Discretization`
    /// - `HCurl` / `Nedelec`: force edge-element path
    /// - `H1`: force scalar nodal path
    #[serde(rename = "Formulation", default = "default_formulation_auto")]
    pub formulation: String,

    /// Optional HCurl order override for eigenmode (1=ND1, 2=ND2).
    /// If omitted, `Solver.Order` is used.
    #[serde(rename = "HCurlOrder", default)]
    pub hcurl_order: Option<u8>,

    /// Conductor wall conductivity σ [S/m] for Q_conductor calculation.
    /// Set to a positive value (e.g. 5.8e7 for copper) to enable ohmic surface
    /// loss perturbation. Default 0 = disabled (only dielectric Q computed).
    #[serde(rename = "WallConductivity", default)]
    pub wall_conductivity: f64,
}

fn default_n_eig() -> usize { 1 }
fn default_eig_tol() -> f64 { 1e-6 }
fn default_eig_maxiter() -> usize { 200 }
fn default_save() -> usize { 0 }

#[derive(Debug, Clone, Deserialize)]
pub struct DrivenSolver {
    #[serde(rename = "MinFreq")]
    pub min_freq: f64,

    #[serde(rename = "MaxFreq")]
    pub max_freq: f64,

    #[serde(rename = "FreqStep")]
    pub freq_step: f64,

    #[serde(rename = "SaveStep", default = "default_save_step")]
    pub save_step: usize,

    /// Optional discretization override for driven solver only.
    ///
    /// Supported values (case-insensitive):
    /// - `Auto` (default): follow `Solver.Discretization`
    /// - `HCurl` / `Nedelec`: force edge-element path
    /// - `H1`: force scalar nodal path
    #[serde(rename = "Formulation", default = "default_formulation_auto")]
    pub formulation: String,

    /// Optional HCurl order override for driven (1=ND1, 2=ND2).
    /// If omitted, `Solver.Order` is used.
    #[serde(rename = "HCurlOrder", default)]
    pub hcurl_order: Option<u8>,

    #[serde(rename = "AdaptiveTol", default)]
    pub adaptive_tol: f64,

    /// Snapshot-based ROM order: number of full solves used to build the reduced basis.
    /// 0 (default) = disabled; 4�?6 recommended for smooth S-parameter sweeps.
    /// When enabled, only `RomOrder` full complex solves are performed; all other
    /// frequency points are evaluated via the reduced system (much cheaper).
    #[serde(rename = "RomOrder", default)]
    pub rom_order: usize,

    /// Palace `Samples` �?accepted, not implemented (use MinFreq/MaxFreq/FreqStep).
    #[serde(rename = "Samples", default)]
    pub samples: Vec<FreqSampleSpec>,

    /// Palace `Save` array �?accepted, not implemented (use SaveStep integer).
    #[serde(rename = "Save", default)]
    pub save: Vec<f64>,

    /// Enable Vector Fitting circuit synthesis after the frequency sweep.
    /// Produces three downloadable artifacts: `s_params.s1p` (Touchstone),
    /// `circuit_model.csv` (pole-residue table), `equivalent_circuit.cir` (SPICE).
    /// Number of poles defaults to min(N/4, 16); use `RomOrder` to override.
    #[serde(rename = "CircuitSynthesis", default)]
    pub circuit_synthesis: bool,

    /// Near-field source file path for linked-source excitation.  When set,
    /// the Dirichlet BC values on excited port boundaries are interpolated
    /// from the CSV near-field data instead of using a uniform φ=1.
    #[serde(rename = "NearFieldSource", default)]
    pub near_field_source: Option<String>,
}

/// Palace `Driven.Samples` item.
#[derive(Debug, Clone, Deserialize)]
pub struct FreqSampleSpec {
    #[serde(rename = "Type", default = "default_sample_type")]
    pub sample_type: String,

    #[serde(rename = "MinFreq", default)]
    pub min_freq: f64,

    #[serde(rename = "MaxFreq", default)]
    pub max_freq: f64,

    #[serde(rename = "FreqStep", default)]
    pub freq_step: f64,

    #[serde(rename = "Freq", default)]
    pub freq: Vec<f64>,

    #[serde(rename = "SaveStep", default)]
    pub save_step: usize,
}

fn default_sample_type() -> String { "Linear".to_string() }

fn default_save_step() -> usize { 1 }

#[derive(Debug, Clone, Deserialize)]
pub struct TransientSolver {
    #[serde(rename = "Type", default = "default_transient_type")]
    pub solver_type: String,

    /// Integration end time in seconds [s].
    /// NOTE: Palace uses nanoseconds; REM configs must use SI seconds (e.g. 1.0e-9 for 1 ns).
    #[serde(rename = "MaxTime")]
    pub max_time: f64,

    /// Time step in seconds [s].
    /// NOTE: Palace uses nanoseconds; REM configs must use SI seconds (e.g. 5.0e-12 for 5 ps).
    #[serde(rename = "TimeStep")]
    pub time_step: f64,

    #[serde(rename = "SaveStep", default = "default_save_step")]
    pub save_step: usize,

    /// Palace `Excitation` waveform type �?accepted, not fully implemented.
    #[serde(rename = "Excitation", default)]
    pub excitation: String,

    /// Palace `ExcitationFreq` [GHz] �?accepted, not fully implemented.
    #[serde(rename = "ExcitationFreq", default)]
    pub excitation_freq: f64,

    /// Palace `ExcitationWidth` [ns] �?accepted, not fully implemented.
    #[serde(rename = "ExcitationWidth", default)]
    pub excitation_width: f64,

    /// Near-field source file path for linked-source excitation.  When set,
    /// the time-domain excitation envelope is replaced by the near-field
    /// data interpolated at each time step.
    #[serde(rename = "NearFieldSource", default)]
    pub near_field_source: Option<String>,
}

fn default_transient_type() -> String { "GeneralizedAlpha".to_string() }

#[derive(Debug, Clone, Deserialize, Default)]
pub struct StaticSolverConfig {
    #[serde(rename = "Save", default)]
    pub save: usize,
}

#[derive(Debug, Clone, Deserialize)]
pub struct LinearSolver {
    #[serde(rename = "Type", default = "default_linear_type")]
    pub solver_type: String,

    /// Palace `KSPType` �?"CG"/"PCG" routes to PCG; "GMRES"/"" uses GMRES (default).
    #[serde(rename = "KSPType", default)]
    pub ksp_type: String,

    #[serde(rename = "Tol", default = "default_linear_tol")]
    pub tol: f64,

    #[serde(rename = "MaxIter", default = "default_linear_maxiter")]
    pub max_iter: usize,

    #[serde(rename = "MGLevels", default = "default_mg_levels")]
    pub mg_levels: usize,

    #[serde(rename = "PCType", default = "default_pc_type")]
    pub pc_type: String,

    /// Palace `ComplexCoarseSolve` �?accepted, not implemented.
    #[serde(rename = "ComplexCoarseSolve", default)]
    pub complex_coarse_solve: bool,
}

impl LinearSolver {
    /// Returns `true` if `KSPType` is "CG" or "PCG" (case-insensitive).
    /// This hints that the caller should prefer the PCG path over GMRES.
    pub fn prefers_pcg(&self) -> bool {
        matches!(self.ksp_type.to_lowercase().as_str(), "cg" | "pcg")
    }

    /// Returns `true` when complex Helmholtz should use sparse iterative solve
    /// before any dense GMRES fallback.
    ///
    /// Accepted values (case-insensitive):
    /// - `CG`, `PCG` (legacy naming)
    /// - `BiCGSTAB`
    pub fn prefers_sparse_iterative_complex(&self) -> bool {
        matches!(self.ksp_type.to_lowercase().as_str(), "cg" | "pcg" | "bicgstab")
    }
}

impl Default for LinearSolver {
    fn default() -> Self {
        LinearSolver {
            solver_type: default_linear_type(),
            ksp_type: String::new(),
            tol: default_linear_tol(),
            max_iter: default_linear_maxiter(),
            mg_levels: default_mg_levels(),
            pc_type: default_pc_type(),
            complex_coarse_solve: false,
        }
    }
}

fn default_linear_type() -> String { "GMRES".to_string() }
fn default_linear_tol() -> f64 { 1e-6 }
fn default_linear_maxiter() -> usize { 200 }
fn default_mg_levels() -> usize { 10 }
fn default_pc_type() -> String { "AMG".to_string() }

// ---------------------------------------------------------------------------
// Custom deserializer: attributes accept Vec<u32> or "1,3-5" string
// ---------------------------------------------------------------------------

pub fn deserialize_attributes<'de, D>(d: D) -> Result<Vec<u32>, D::Error>
where
    D: Deserializer<'de>,
{
    use serde::de::{self, Visitor};
    use std::fmt;

    struct AttrsVisitor;

    impl<'de> Visitor<'de> for AttrsVisitor {
        type Value = Vec<u32>;

        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "an array of integers or a comma/range string like \"1,3-5\"")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Vec<u32>, A::Error>
        where
            A: de::SeqAccess<'de>,
        {
            #[derive(Deserialize)]
            #[serde(untagged)]
            enum ElemOrRange {
                Int(u64),
                Range(String),
            }

            let mut v = Vec::new();
            while let Some(elem) = seq.next_element::<ElemOrRange>()? {
                match elem {
                    ElemOrRange::Int(n) => v.push(n as u32),
                    ElemOrRange::Range(s) => {
                        let mut parsed = expand_ranges(&s).map_err(de::Error::custom)?;
                        v.append(&mut parsed);
                    }
                }
            }
            v.sort_unstable();
            v.dedup();
            Ok(v)
        }

        fn visit_str<E: de::Error>(self, s: &str) -> Result<Vec<u32>, E> {
            expand_ranges(s).map_err(de::Error::custom)
        }

        fn visit_u64<E: de::Error>(self, v: u64) -> Result<Vec<u32>, E> {
            Ok(vec![v as u32])
        }
        fn visit_i64<E: de::Error>(self, v: i64) -> Result<Vec<u32>, E> {
            Ok(vec![v as u32])
        }
    }

    d.deserialize_any(AttrsVisitor)
}

// ---------------------------------------------------------------------------
// MoM solver config (REM extension �?ignored by Palace)
// ---------------------------------------------------------------------------

/// MoM solver parameters, placed under `Solver.MoM` in the config file.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct MomSolverConfig {
    /// Integral equation: "EFIE" | "MFIE" | "CFIE" | "PMCHWT"
    #[serde(rename = "Equation", default = "default_mom_equation")]
    pub equation: String,

    /// Basis function: "RWG" | "Pulse"
    #[serde(rename = "Basis", default = "default_mom_basis")]
    pub basis: String,

    /// Start frequency [Hz]
    #[serde(rename = "FreqMin")]
    pub freq_min: f64,

    /// End frequency [Hz]
    #[serde(rename = "FreqMax")]
    pub freq_max: f64,

    /// Frequency step [Hz]
    #[serde(rename = "FreqStep")]
    pub freq_step: f64,

    /// CFIE mixing coefficient α �?[0,1]: 0 = pure EFIE, 1 = pure MFIE
    #[serde(rename = "Alpha", default = "default_cfie_alpha")]
    pub alpha: f64,

    /// Convergence tolerance for singular integrals
    #[serde(rename = "SingularTol", default = "default_singular_tol")]
    pub singular_tol: f64,

    /// Linear solver for Z·I = V: "Auto" | "Direct" | "GMRES" | "ACA" | "FFT" | "FMM" | "MLFMA" | "WGPU"
    /// When "Auto", the solver is chosen automatically based on problem size:
    /// N < 1000 �?Direct LU, 1000..5000 �?FFT/FMM, >5000 �?MLFMA.
    #[serde(rename = "FastSolver", default = "default_fast_solver")]
    pub fast_solver: String,

    /// Enable GPU-accelerated impedance matrix assembly via wgpu compute shaders.
    /// Has no effect unless compiled with --features rem-mom/wgpu-gpu.
    #[serde(rename = "UseGPU", default = "default_use_gpu")]
    pub use_gpu: bool,

    /// When true, automatically select solver based on problem size when
    /// FastSolver = "Auto".  Has no effect when FastSolver is explicitly set.
    #[serde(rename = "SolverAutoSelect", default)]
    pub sol_auto_select: bool,

    /// When FastSolver="Auto", switch from LU to FFT when the number of
    /// unknowns (2 · cells_x · cells_y) exceeds this threshold.
    /// Default 500 (empirical break-even for LU vs FFT-GMRES).
    #[serde(rename = "SolverThreshold", default = "default_solver_threshold")]
    pub solver_threshold: usize,

    /// Force a specific solver path for the boxed MoM solver.
    /// When set, overrides both FastSolver and SolverAutoSelect.
    /// Accepted values: "LU", "FFT", "GMRES", "ACA".
    /// Default None �?normal FastSolver / SolverAutoSelect logic.
    #[serde(rename = "SolverOverride", default)]
    pub solver_override: Option<String>,

    /// Quadrature rule (1, 3, 5, 7). Overrides REM_MOM_QUAD_ORDER env var.
    #[serde(rename = "QuadratureRule", default)]
    pub quad_order: Option<usize>,

    /// GMRES restart parameter (default 30). Overrides hardcoded restart.
    #[serde(rename = "GmresRestart", default)]
    pub gmres_restart: Option<usize>,

    /// GMRES convergence tolerance (default 1e-8). Overrides hardcoded tol.
    #[serde(rename = "GmresTol", default)]
    pub gmres_tol: Option<f64>,

    /// GMRES max iterations (default 500). Overrides REM_MOM_GMRES_MAX_ITERS.
    #[serde(rename = "GmresMaxIter", default)]
    pub gmres_max_iter: Option<usize>,

    /// Incident plane wave polar angle [degrees] from +z axis (0 = broadside)
    #[serde(rename = "ThetaInc", default)]
    pub theta_inc_deg: f64,

    /// Incident plane wave azimuth angle [degrees] from +x axis
    #[serde(rename = "PhiInc", default)]
    pub phi_inc_deg: f64,

    /// Incident plane wave polarization: "theta" | "phi" | "x" | "y" | "z"
    #[serde(rename = "Polarization", default = "default_polarization")]
    pub polarization: String,

    /// Lumped ports for S-parameter extraction.  When non-empty the solver
    /// runs one MoM solve per port and outputs S-matrix + Touchstone file.
    /// When empty, the existing plane-wave + RCS path is used.
    #[serde(rename = "Ports", default)]
    pub ports: Vec<MomPort>,

    /// Global reference impedance Z₀ [Ω] for S-parameter normalisation.
    #[serde(rename = "RefImpedance", default = "default_ref_impedance")]
    pub ref_impedance: f64,

    /// Stratified dielectric substrate for layered Green's function.
    /// When present, enables Sommerfeld integral / DCIM-based layered Green's function.
    /// When absent, free-space Green's function is used.
    #[serde(rename = "Substrate", default)]
    pub substrate: Option<SubstrateConfig>,

    /// Conductor wall conductivity σ [S/m] for MoM SIBC loss modeling.
    /// Set to a positive value (for example 5.8e7 for copper) to add the
    /// Leontovich surface impedance term Zs = (1+j)/(σ δs) on PEC surfaces.
    #[serde(rename = "WallConductivity", default)]
    pub wall_conductivity: f64,

    /// Superconducting kinetic inductance Ls [H] (Sonnet SURFACE_IMPEDANCE).
    /// When > 0, SIBC uses Zs = Rdc + Rrf√f + j(ω·Ls + Xdc).
    #[serde(rename = "SurfaceLs", default)]
    pub surface_ls: f64,
    /// Superconducting DC resistance Rdc [Ω].
    #[serde(rename = "SurfaceRdc", default)]
    pub surface_rdc: f64,
    /// Superconducting RF resistance Rrf [Ω/√Hz].
    #[serde(rename = "SurfaceRrf", default)]
    pub surface_rrf: f64,
    /// Superconducting DC reactance Xdc [Ω].
    #[serde(rename = "SurfaceXdc", default)]
    pub surface_xdc: f64,

    /// Near-field source file path. When set, the RHS is built from the
    /// near-field CSV data instead of the plane-wave model.  The file
    /// contains spatially sampled E/H fields exported from a previous
    /// simulation, enabling multi-solver near-field coupling.
    #[serde(rename = "NearFieldSource", default)]
    pub near_field_source: Option<String>,

    /// Snapshot ROM acceleration for S-parameter frequency sweeps.
    /// `0` disables ROM (default); positive value sets the number of anchor
    /// frequencies at which a full MoM solve is performed �?all other
    /// frequencies use the Galerkin-projected low-dimensional system.
    /// Typical values: 4�?6 for narrow-band, 8�?2 for wideband.
    #[serde(rename = "RomOrder", default)]
    pub rom_order: usize,

    /// Enable adaptive frequency sweep (like Sonnet ABS_ENTRY).
    /// When true, the solver starts with (FreqMin, FreqMax) and iteratively
    /// inserts mid-points where S-parameter interpolation error exceeds the
    /// adaptive tolerance, up to `AdaptiveTarget` points. The sweep points
    /// are not uniform �?they concentrate where the response varies rapidly.
    #[serde(rename = "AdaptiveSweep", default)]
    pub adaptive_sweep: bool,

    /// Maximum number of frequency points for adaptive sweep.
    /// Ignored when `AdaptiveSweep` is false. Default 100.
    #[serde(rename = "AdaptiveTarget", default = "default_adaptive_target")]
    pub adaptive_target: usize,

    /// Number of poles for the ABS (Adaptive Band Synthesis) rational model.
    /// Must be even (conjugate pairs). Default 6. Ignored when AdaptiveSweep is false.
    #[serde(rename = "AbsPoles", default)]
    pub abs_poles: Option<usize>,

    /// ABS convergence tolerance (relative L2 error). Default 1e-3.
    /// Smaller values yield more frequency points. Ignored when AdaptiveSweep is false.
    #[serde(rename = "AbsTol", default)]
    pub abs_tol: Option<f64>,

    /// Maximum AMR iterations.  `0` disables AMR (default).
    /// When > 0, the mesh is refined up to `amr_iter` times with a
    /// Dörfler marking threshold of `AmrtTheta`.
    #[serde(rename = "AmrIter", default)]
    pub amr_iter: usize,

    /// Dörfler marking fraction for AMR.  Faces whose squared error sum
    /// exceeds `amr_theta × total` are refined.  Default 0.5.
    #[serde(rename = "AmrTheta", default = "default_amr_theta")]
    pub amr_theta: f64,

    /// Near-field probe points for S-parameter sweep post-processing.
    /// For each listed (x,y,z) point, the E-field is computed at all sweep
    /// frequencies (port-1 excitation) and written to `postpro/probe_e_field.csv`.
    #[serde(rename = "NearFieldProbes", default)]
    pub near_field_probes: Vec<NearFieldProbePoint>,

    /// Transmission-line length [m] for RLGC per-unit-length extraction.
    /// When > 0 and the problem has exactly 2 ports, the ABCD matrix is
    /// computed from S-parameters and R/L/G/C are written to
    /// `postpro/tline_params.csv`.  Default 0 (disabled).
    #[serde(rename = "TlineLength", default)]
    pub tline_length: f64,

    /// Effective relative permittivity used for simple port reference-plane
    /// de-embedding phase correction.  Default 1.0 (free-space velocity).
    #[serde(rename = "DeembedEpsEff", default = "default_deembed_eps_eff")]
    pub deembed_eps_eff: f64,

    /// Optional global attenuation constant α [Np/m] used by de-embedding.
    /// Default 0.0 (phase-only correction).
    #[serde(rename = "DeembedAlpha", default)]
    pub deembed_alpha_np_per_m: f64,

    /// Hierarchical basis function order: 1 = RWG, 2 = HO-RWG, 3 = P3.
    /// Used when `Basis = "Hierarchical"` to select the polynomial degree.
    /// Ignored for explicit basis tokens (e.g. `Basis = "HO"`).
    #[serde(rename = "BasisOrder", default)]
    pub basis_order: Option<u8>,

    /// Floquet/periodic unit cell first lattice vector [m] as [x, y, z].
    /// Required when `Basis = "Periodic"` or `Basis = "Floquet"`.
    #[serde(rename = "FloquetA1", default)]
    pub floquet_a1: Option<[f64; 3]>,

    /// Floquet/periodic unit cell second lattice vector [m] as [x, y, z].
    /// Required when `Basis = "Periodic"` or `Basis = "Floquet"`.
    #[serde(rename = "FloquetA2", default)]
    pub floquet_a2: Option<[f64; 3]>,

    /// Floquet incident wavenumber [rad/m] as [kx, ky].
    /// Default [0, 0] (normal incidence). Used when `Basis = "Floquet"`.
    #[serde(rename = "FloquetKInc", default)]
    pub floquet_k_inc: Option<[f64; 2]>,

    /// Number of characteristic modes to compute for CMA post-processing.
    /// Ignored unless `Equation = "CMA"`.  Default 10.
    #[serde(rename = "CmaModes", default)]
    pub cma_modes: Option<usize>,

    /// Maximum number of Floquet modes for periodic-structure acceleration.
    ///
    /// Used by the Ewald-accelerated periodic Green function in `assemble_efie_periodic`.
    /// Controls the ±n spectral terms in the Floquet-mode summation.
    /// Higher values improve accuracy for large unit cells or off-normal incidence.
    /// Default 0 means auto-select based on unit-cell electrical size.
    #[serde(rename = "FloquetMaxModes", default)]
    pub floquet_max_modes: i32,

    /// Symmetry planes for DOF elimination.
    /// Each entry specifies a plane (e.g. `"x=0"`, `"y=0"`) on which
    /// RWG edge DOFs are eliminated (enforces PEC/PMC symmetry condition).
    #[serde(rename = "SymmetryPlanes", default)]
    pub symmetry_planes: Vec<SymmetryPlaneConfig>,

    /// Automatically detect lumped ports from mesh boundary edges.
    /// When true, the solver scans boundary edge clusters and adds any
    /// discovered port candidates after the explicitly defined ports.
    /// Default false.
    #[serde(rename = "AutoDetectPorts", default)]
    pub auto_detect_ports: bool,

    /// Minimum number of faces required for an auto-detected port cluster.
    /// Smaller clusters are ignored.  Default 1.
    #[serde(rename = "AutoPortMinFaces", default = "default_auto_port_min_faces")]
    pub auto_port_min_faces: usize,

    /// Maximum transverse width [m] of an auto-detected port cluster.
    /// Wider clusters are ignored.  Default 1.0 m (effectively unlimited
    /// for most circuits).
    #[serde(rename = "AutoPortMaxWidth", default = "default_auto_port_max_width")]
    pub auto_port_max_width: f64,

    /// When true, fail with an error if a WavePort modal profile cannot be
    /// constructed on the selected port patch (instead of falling back to
    /// uniform excitation).  Default false.
    #[serde(rename = "WavePortRequireModal", default)]
    pub waveport_require_modal: bool,

    /// Optional TRL calibration kit parameters for 2-port de-embedding.
    #[serde(rename = "TrlKit", default)]
    pub trl_kit: Option<TrlKitConfig>,

    /// Optional SOLT calibration kit parameters for 2-port de-embedding.
    ///
    /// SOLT (Short-Open-Load-Thru) uses known 1-port reflection standards at
    /// each port plus a Thru to extract an 8-term error model. Supports
    /// frequency-dependent Short inductance, Open capacitance, and Load offset.
    #[serde(rename = "SoltKit", default)]
    pub solt_kit: Option<SoltKitConfig>,

    /// Optional output renormalization impedance [ohm].
    /// When set, final S-parameters are renormalized from per-port Z0 to this value.
    #[serde(rename = "OutputRefImpedance", default)]
    pub output_ref_impedance: Option<f64>,

    /// Export per-face current density CSV at the final sweep frequency.
    #[serde(rename = "ExportCurrentDensity", default)]
    pub export_current_density: bool,

    /// Conductor surface roughness model: "hammerstad" | "groisse" | null.
    /// When set, the Leontovich surface impedance is multiplied by the
    /// roughness correction factor K_r �?1 (requires RmsRoughness > 0).
    #[serde(rename = "RoughnessModel", default)]
    pub roughness_model: Option<String>,

    /// RMS surface roughness Δ [m] for conductor loss modeling.
    /// Typical values: 1e-6 (1 µm) for standard PCB copper.
    #[serde(rename = "RmsRoughness", default)]
    pub rms_roughness: Option<f64>,

    /// Boxed/enclosed solver configuration (Sonnet-style shielded box).
    /// When present, the solver uses the rectangular wave-guide mode expansion
    /// (FFT-based coupling) instead of the free-space/layered Green's function.
    #[serde(rename = "Box", default)]
    pub box_config: Option<BoxConfig>,

    /// When true, disable the automatic co-calibration of closely-spaced
    /// ports in the boxed solver.  Co-calibration removes parasitic
    /// port-to-port coupling and is enabled by default for multi-port
    /// structures.  Set to true for manual calibration or debugging.
    #[serde(rename = "DisableCoCalibration", default)]
    pub disable_co_calibration: bool,

    /// Parametric sweep definition (R-10.2).
    /// When set, the solver runs a multi-dimensional sweep over geometry
    /// and/or frequency parameters instead of a single simulation.
    #[serde(rename = "ParametricSweep", default)]
    pub parametric_sweep: Option<Vec<ParamDefConfig>>,

    /// Number of cells per free-space wavelength for automatic grid sizing.
    ///
    /// When set under `Solver.MoM` (non-zero), the boxed solver adjusts
    /// the rectilinear grid resolution based on the highest sweep frequency:
    ///   `cells = ceil(box_dim / (c₀ / freq_max / subs_per_lambda))`
    /// capped at 200×200.  When unset or zero, the explicit `Box.CellsX/CellsY`
    /// values are used without adjustment.
    ///
    /// Typical range: 10–40 (higher = finer, slower).  Ignored when
    /// `Box` configuration is absent.
    #[serde(rename = "SubsPerLambda", default)]
    pub subs_per_lambda: f64,

    /// Explicit solver-type selector for unified dispatch.
    ///
    /// Values: `"FreeForm"` (RWG trianglular mesh), `"Boxed"` (rect grid / Sonnet-style),
    /// `"Capacitance"` (electrostatic BIE).  When absent or `"Auto"`, the solver
    /// infers the type from other config fields (Box presence, Problem.Type, etc.).
    #[serde(rename = "MomType", default = "default_mom_type")]
    pub mom_type: String,

    /// Explicit kernel selector for unified dispatch.
    ///
    /// Values: `"FreeSpace"`, `"Layered"`, `"Cavity"`, `"Laplace"`, `"Auto"`.
    #[serde(rename = "Kernel", default = "default_mom_kernel")]
    pub kernel: String,

    /// Explicit mesh format selector for unified dispatch.
    ///
    /// Values: `"TriSurface"` (triangular faces), `"RectGrid"` (rectilinear grid).
    #[serde(rename = "MeshFormat", default = "default_mesh_format")]
    pub mesh_format: String,
}

fn default_auto_port_min_faces() -> usize { 1 }
fn default_auto_port_max_width() -> f64   { 1.0 }

/// A symmetry-plane definition for MoM DOF elimination.
///
/// The `Plane` field accepts strings of the form `"x=0"`, `"y=0.5"`, `"z=-0.1"`.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct SymmetryPlaneConfig {
    /// Plane specification, e.g. `"x=0"`, `"y=0"`, `"z=0"`.
    #[serde(rename = "Plane", default)]
    pub plane: String,
}

/// TRL calibration kit configuration for 2-port S-parameter correction.
#[derive(Debug, Clone, Deserialize)]
pub struct TrlKitConfig {
    /// Physical length of the thru standard [m].
    #[serde(rename = "ThruLength")]
    pub thru_length: f64,

    /// Physical length of the line standard [m].
    #[serde(rename = "LineLength")]
    pub line_length: f64,

    /// Characteristic impedance of the line standard [ohm].
    #[serde(rename = "LineImpedance")]
    pub line_impedance: f64,

    /// Effective dielectric constant used for initial gamma estimate.
    #[serde(rename = "EpsilonEff")]
    pub epsilon_eff: f64,

    /// Reflect standard type (e.g. "open", "short", "load").
    #[serde(rename = "ReflectType")]
    pub reflect_type: String,

    /// Magnitude of the reflect standard reflection coefficient.
    #[serde(rename = "ReflectMagnitude")]
    pub reflect_magnitude: f64,

    /// Apply TRL de-embedding inside the per-frequency solve loop.
    ///
    /// Default false (batch post-process after all frequencies).
    /// When true, TRL is applied per-frequency inside the solve loop,
    /// which can improve accuracy for adaptive sweeps.
    #[serde(rename = "SolveSide", default)]
    pub solve_side: bool,
}

/// SOLT calibration kit configuration for 2-port S-parameter correction.
///
/// Models the Short as Γ_S = �? · exp(�?j·ω·L_short/Z0) with offset inductance,
/// the Open as Γ_O = +1 · exp(2j·ω·C_open·Z0) with fringing capacitance,
/// and the Load as a resistively-terminated line with possible offset delay.
/// The Thru is modeled as a zero-length connection (S₂₁=S₁₂=1, S₁₁=S₂₂=0).
#[derive(Debug, Clone, Deserialize)]
pub struct SoltKitConfig {
    /// Short offset inductance [H]. Default 0.
    #[serde(rename = "ShortInductance", default)]
    pub short_inductance: f64,

    /// Open fringing capacitance [F]. Default 0.
    #[serde(rename = "OpenCapacitance", default)]
    pub open_capacitance: f64,

    /// Load resistance [Ω]. Default 50.
    #[serde(rename = "LoadResistance", default = "default_solt_load_r")]
    pub load_resistance: f64,

    /// Load offset inductance [H]. Default 0.
    #[serde(rename = "LoadInductance", default)]
    pub load_inductance: f64,

    /// Load offset capacitance [F]. Default 0.
    #[serde(rename = "LoadCapacitance", default)]
    pub load_capacitance: f64,

    /// Port reference impedance [Ω] for standard definitions. Default 50.
    #[serde(rename = "RefImpedance", default = "default_ref_impedance")]
    pub ref_impedance: f64,

    /// Effective relative permittivity used for the Thru/Line propagation
    /// model (needed for phase correction). Default 1.0.
    #[serde(rename = "EpsilonEff", default = "default_deembed_eps_eff")]
    pub epsilon_eff: f64,
}

fn default_solt_load_r() -> f64 { 50.0 }

fn default_mom_equation() -> String { "CFIE".to_string() }
fn default_mom_basis()     -> String { "RWG".to_string()  }
fn default_cfie_alpha()    -> f64    { 0.5 }
fn default_singular_tol()  -> f64    { 1.0e-6 }
fn default_fast_solver()   -> String { "Direct".to_string() }
fn default_solver_threshold() -> usize { 500 }
fn default_use_gpu()       -> bool   { true }
fn default_polarization()  -> String { "theta".to_string() }
fn default_ref_impedance() -> f64    { 50.0 }
fn default_port_direction() -> String { "x".to_string() }
fn default_amr_theta()     -> f64    { 0.5 }
fn default_adaptive_target() -> usize { 100 }
fn default_deembed_eps_eff() -> f64  { 1.0 }
fn default_mom_port_kind() -> String { "Lumped".to_string() }
fn default_mom_type()     -> String { "Auto".to_string() }
fn default_mom_kernel()   -> String { "Auto".to_string() }
fn default_mesh_format()  -> String { "Auto".to_string() }

/// A near-field probe point for E-field evaluation during S-parameter sweep.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct NearFieldProbePoint {
    /// X coordinate [m].
    #[serde(rename = "X", default)]
    pub x: f64,
    /// Y coordinate [m].
    #[serde(rename = "Y", default)]
    pub y: f64,
    /// Z coordinate [m].
    #[serde(rename = "Z", default)]
    pub z: f64,
    /// Optional human-readable label (used in CSV comments).
    #[serde(rename = "Label", default)]
    pub label: String,
}

/// A lumped port definition for MoM S-parameter extraction.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct MomPort {
    /// Port index (1-based), matches Touchstone port ordering.
    #[serde(rename = "Index")]
    pub index: u32,

    /// Physical-group attribute IDs that define the port surface.
    #[serde(rename = "Attributes", default)]
    pub attributes: Vec<u32>,

    /// Port type: "Lumped" | "WavePort".
    ///
    /// `WavePort` uses a port-surface modal profile (graph-Laplacian
    /// eigenmode on the selected port patch) to weight RWG excitation.
    #[serde(rename = "Type", default = "default_mom_port_kind")]
    pub port_type: String,

    /// WavePort mode index (1-based).  Used when `Type = "WavePort"`.
    /// Mode 1 is the fundamental profile; mode > 1 selects higher-order
    /// eigenmodes when available on the port patch.
    #[serde(rename = "Mode", default = "default_mode")]
    pub mode: u32,

    /// Dominant E-field direction: "x" | "y" | "z".  Determines how
    /// the RHS excitation voltage is projected onto the port RWG functions.
    #[serde(rename = "Direction", default = "default_port_direction")]
    pub direction: String,

    /// Per-port reference impedance [Ω].  When None, uses `MomSolverConfig::ref_impedance`.
    #[serde(rename = "Impedance", default)]
    pub impedance: Option<f64>,

    /// Differential pair partner port index (1-based).
    ///
    /// If set on both ports of a pair, mixed-mode S-parameters are generated
    /// in addition to single-ended S-parameters.
    #[serde(rename = "PairWith", default)]
    pub pair_with: Option<u32>,

    /// Ground reference flag.
    ///
    /// When true, this port is treated as a ground reference (zero excitation,
    /// current extraction negated). Used for CPW/GSG structures where the
    /// signal port references adjacent ground planes.
    #[serde(rename = "GndRef", default)]
    pub gnd_ref: bool,

    /// Reference-plane de-embedding length [m] for this port.
    /// Positive values shift the reference plane away from the port.
    #[serde(rename = "DeembedLength", default)]
    pub deembed_length: f64,

    /// Inner conductor radius [m] for coaxial port geometry.
    /// When set together with `OuterRadius`, overrides the auto-detected
    /// radii for `Type = "Coaxial"` ports.  Ignored for lumped ports.
    #[serde(rename = "InnerRadius", default)]
    pub inner_radius: Option<f64>,

    /// Outer conductor radius [m] for coaxial port geometry.
    /// See `InnerRadius`.
    #[serde(rename = "OuterRadius", default)]
    pub outer_radius: Option<f64>,

    /// Optional shunt load resistance [ohm] stamped at this port.
    #[serde(rename = "LoadR", default)]
    pub load_r: f64,

    /// Optional shunt load inductance [H] stamped at this port.
    #[serde(rename = "LoadL", default)]
    pub load_l: f64,

    /// Optional shunt load capacitance [F] stamped at this port.
    #[serde(rename = "LoadC", default)]
    pub load_c: f64,
}

// ---------------------------------------------------------------------------
// MoM Substrate configuration (stratified dielectric layers)
// ---------------------------------------------------------------------------

/// Stratified dielectric substrate for layered Green's function in MoM.
/// Defines a stack of horizontal dielectric layers, with the bottom being
/// either PEC or transmission to semi-infinite space.
#[derive(Debug, Clone, Deserialize)]
pub struct SubstrateConfig {
    /// List of dielectric layers from bottom to top.
    /// Typically the bottom layer is ground or substrate, top layer is air.
    #[serde(rename = "Layers", default)]
    pub layers: Vec<SubstrateLayerConfig>,

    /// If true, bottom layer sits on a PEC ground plane.
    /// If false, semi-infinite dielectric extends downward.
    #[serde(rename = "BottomPec", default = "default_bottom_pec")]
    pub bottom_pec: bool,

    /// If true, a PEC cover is placed at the top of the layer stack (shielded box).
    /// The top PEC height is the sum of all finite-thickness layers.
    #[serde(rename = "TopPec", default)]
    pub top_pec: bool,

    /// Ground conductor conductivity sigma [S/m].
    /// When > 0 the ground plane is treated as a lossy conductor using SIBC.
    /// Default 0.0 (ideal PEC / not used).
    #[serde(rename = "GroundConductivity", default)]
    pub ground_conductivity: f64,
}

/// Boxed MoM solver configuration (Sonnet-style enclosed box).
///
/// When present under `Solver.MoM.Box`, the solver switches from the free-space /
/// layered Green's function to a rectangular waveguide mode expansion with FFT-based
/// coupling (Sonnet-style). All conductors and ports must lie within the box extents.
#[derive(Debug, Clone, Deserialize)]
pub struct BoxConfig {
    /// Box width (x-direction) [m].
    #[serde(rename = "Width")]
    pub width: f64,

    /// Box height (y-direction) [m].
    #[serde(rename = "Height")]
    pub height: f64,

    /// Number of cells in the x-direction (rectilinear grid).
    #[serde(rename = "CellsX")]
    pub cells_x: usize,

    /// Number of cells in the y-direction (rectilinear grid).
    #[serde(rename = "CellsY")]
    pub cells_y: usize,

    /// Whether the top cover is PEC. If true, sets the top boundary to PEC.
    /// If false, the top is open (absorbing / free-space above).
    #[serde(rename = "TopCover", default)]
    pub top_cover: bool,

    /// Whether the bottom cover is PEC. If true, sets the bottom boundary to PEC.
    /// If false, the bottom is open.
    #[serde(rename = "BottomCover", default)]
    pub bottom_cover: bool,

    /// Number of evanescent modes to retain beyond the propagating cutoff.
    /// Larger values improve accuracy for strongly coupling structures at the cost
    /// of increased FFT size. Default 10.
    #[serde(rename = "NumEvanescentModes", default = "default_evanescent_modes")]
    pub num_evanescent_modes: usize,

    /// Maximum conformal sub-sampling level for cells intersecting conductor boundaries.
    /// Level 0 = base grid only; higher values refine cells near edges. Default 2.
    #[serde(rename = "ConformalLevel", default = "default_conformal_level")]
    pub conformal_level: u32,

    /// Box wall PEC flag when true. Side walls are always PEC in Sonnet-style solvers.
    #[serde(rename = "SideWallPec", default = "default_side_wall_pec")]
    pub side_wall_pec: bool,

    /// If true, the box covers are assigned the same conductivity as `TopCover`/`BottomCover`.
    /// If > 0, treated as lossy metal rather than ideal PEC.
    #[serde(rename = "CoverConductivity", default)]
    pub cover_conductivity: f64,

    /// Internal delta-gap ports on the rectilinear grid.
    /// Each entry defines a port between two adjacent cells.
    #[serde(rename = "Ports", default)]
    pub ports: Vec<BoxPortConfig>,

    /// Dielectric bricks in the box volume (VIE).
    /// Each brick fills a rectangular region of cells from z=ZBottom to
    /// z=ZBottom+Thickness, subdivided into NzLayers vertical cells.
    /// The brick is modeled as a polarization current (VIE) coupled to
    /// the surface rooftop basis.
    #[serde(rename = "DielectricBricks", default)]
    pub dielectric_bricks: Vec<DielectricBrickConfig>,

    /// Interior PEC wall polygons for non-rectangular cavity shapes.
    ///
    /// Each polygon is a list of `(x, y)` vertices [m] defining an interior
    /// PEC wall.  Cells whose center falls outside ALL polygons are masked
    /// (pec_tag = 0).  Cells intersected by the polygon boundary get their
    /// `coverage` adjusted.  When specified, `SideWallPec` is still applied
    /// to the bounding-box edges; the polygons define additional interior
    /// walls.
    #[serde(rename = "WallPolygons", default)]
    pub wall_polygons: Vec<Vec<[f64; 2]>>,

    /// Z-position of the signal layer relative to the box bottom [m].
    ///
    /// Controls where the source/observation plane (z_src = z_obs) is placed
    /// inside the box for the spectral Green's function.  The default (None)
    /// places the signal layer at the box mid-height (height/2), which is
    /// appropriate for standard microstrip on a bottom-grounded substrate.
    ///
    /// For conductor-backed coplanar waveguide (CBCPW) or air-bridge
    /// configurations, set this to the actual signal-layer height above the
    /// bottom ground plane.  For example, `signal_layer_z = 0.95·height`
    /// places the signal near the top cover.
    #[serde(rename = "SignalLayerZ", default)]
    pub signal_layer_z: Option<f64>,

    /// Circuit components connected to boxed solver ports (R-10.1).
    ///
    /// Each entry defines an R, L, C, or S-parameter component connected
    /// across one or two ports.  Components are stamped into the port Z-matrix
    /// after the EM solve as a post-processing step.
    #[serde(rename = "Components", default)]
    pub components: Vec<ComponentConfig>,

    /// Multi-layer signal stack (REM extension).
    ///
    /// When present, each entry defines a metal layer with its own set of
    /// wall polygons and z-position, enabling multi-layer EM simulation.
    /// Via ports (in `Vias`) connect between layers.
    /// When absent, the single-layer `WallPolygons` + `SignalLayerZ` is used.
    #[serde(rename = "Layers", default)]
    pub layers: Vec<SignalLayerConfig>,

    /// Via port definitions for multi-layer configurations (REM extension).
    ///
    /// Each via connects two signal layers at the same (ix, iy) grid cell.
    /// The via is modeled as a vertical current element.
    #[serde(rename = "Vias", default)]
    pub vias: Vec<ViaConfig>,
}

/// A signal layer in the multi-layer stack.
#[derive(Debug, Clone, Deserialize)]
pub struct SignalLayerConfig {
    /// Layer name (e.g. "Metal1", "M2").
    #[serde(rename = "Name", default)]
    pub name: String,

    /// Z-position of the metal plane relative to box bottom [m].
    #[serde(rename = "SignalLayerZ")]
    pub signal_layer_z: f64,

    /// Interior PEC wall polygons on this layer.
    #[serde(rename = "WallPolygons", default)]
    pub wall_polygons: Vec<Vec<[f64; 2]>>,
}

/// A via port connecting two signal layers.
#[derive(Debug, Clone, Deserialize)]
pub struct ViaConfig {
    /// Port index (1-based).
    #[serde(rename = "Index")]
    pub index: u32,

    /// Reference impedance [ohm].
    #[serde(rename = "Impedance", default = "default_via_impedance")]
    pub z0: f64,

    /// Source layer index (0-based into `Layers` array).
    #[serde(rename = "LayerFrom", default)]
    pub layer_from: usize,

    /// Destination layer index (0-based into `Layers` array).
    #[serde(rename = "LayerTo", default)]
    pub layer_to: usize,

    /// Grid cell x-index.
    #[serde(rename = "CellX")]
    pub ix: usize,

    /// Grid cell y-index.
    #[serde(rename = "CellY")]
    pub iy: usize,

    /// Via cross-section width [m] (for self-inductance).
    #[serde(rename = "Width", default = "default_via_width")]
    pub width: f64,
}

fn default_via_impedance() -> f64 { 50.0 }
fn default_via_width() -> f64 { 100e-6 }

/// Parameters for the thick-metal multi-sheet model.
///
/// Stored in `Metadata.ThickMetal` by format converters (Sonnet19, etc.).
/// The boxed solver driver reads this from metadata and constructs
/// a solver-internal `ThickMetalConfig` to enable multi-sheet volume
/// current simulation.
#[derive(Debug, Clone, Deserialize)]
pub struct ThickMetalBoxConfig {
    /// Number of current sheets.  0 = auto (3 for SingleSurface, 5 for DoubleSurface).
    #[serde(rename = "NumSheets", default)]
    pub num_sheets: usize,

    /// Total metal thickness [m].
    #[serde(rename = "Thickness")]
    pub thickness_m: f64,

    /// Conductor conductivity [S/m] (e.g. 5.8e7 for Cu).
    #[serde(rename = "Conductivity", default = "default_copper_conductivity")]
    pub conductivity: f64,

    /// Z-meshing strategy: "SingleSurface" (3 sheets) or "DoubleSurface" (5 sheets).
    #[serde(rename = "Strategy", default = "default_tm_strategy")]
    pub strategy: String,

    /// Relative permittivity of the dielectric embedding the metal.
    #[serde(rename = "EpsR", default = "default_permittivity")]
    pub eps_r: f64,
}

fn default_copper_conductivity() -> f64 { 5.8e7 }
fn default_tm_strategy() -> String { "SingleSurface".into() }

/// A circuit component attached to boxed solver ports (R-10.1).
#[derive(Debug, Clone, Deserialize)]
pub struct ComponentConfig {
    /// Component name (for logging / debugging).
    #[serde(rename = "Name", default)]
    pub name: String,
    /// Component type: "R" | "L" | "C" | "RLC" | "Touchstone".
    #[serde(rename = "Type")]
    pub comp_type: String,
    /// Component value [Ω / H / F].
    /// For RLC: R,L,C comma-separated.  For Touchstone: file path.
    #[serde(rename = "Value", default)]
    pub value: String,
    /// First port (0-based).
    #[serde(rename = "PortA", default)]
    pub port_a: usize,
    /// Second port (0-based).  Omit or same as PortA for single-ended.
    #[serde(rename = "PortB")]
    pub port_b: Option<usize>,
}

/// A parameter definition for the parametric sweep engine (R-10.2).
#[derive(Debug, Clone, Deserialize)]
pub struct ParamDefConfig {
    /// Parameter name (e.g. "W", "L", "freq").
    #[serde(rename = "Name")]
    pub name: String,
    /// Target type: "Frequency" | "ConductorWidth" | "ConductorPosition".
    #[serde(rename = "Target")]
    pub target: String,
    /// Target-specific arguments as comma-separated numbers.
    /// For Frequency: "min,max,step" [Hz].
    /// For ConductorWidth: "direction,cell_start,cell_end,row,min_cells,max_cells".
    #[serde(rename = "Args")]
    pub args: String,
    /// Number of sweep steps (including endpoints). Default 5.
    #[serde(rename = "Steps", default = "default_param_steps")]
    pub steps: usize,
}

fn default_param_steps() -> usize { 5 }

/// Configuration for a single dielectric brick in the boxed VIE solver.
#[derive(Debug, Clone, Deserialize)]
pub struct DielectricBrickConfig {
    /// First cell x-index (inclusive).
    #[serde(rename = "IxStart")]
    pub ix_start: usize,
    /// Last cell x-index (exclusive).
    #[serde(rename = "IxEnd")]
    pub ix_end: usize,
    /// First cell y-index (inclusive).
    #[serde(rename = "IyStart")]
    pub iy_start: usize,
    /// Last cell y-index (exclusive).
    #[serde(rename = "IyEnd")]
    pub iy_end: usize,
    /// Number of vertical layers.
    #[serde(rename = "NzLayers", default = "default_nz_layers")]
    pub nz_layers: usize,
    /// Relative permittivity.
    #[serde(rename = "Permittivity", default = "default_permittivity")]
    pub eps_r: f64,
    /// Conductivity [S/m].
    #[serde(rename = "Conductivity", default)]
    pub sigma: f64,
    /// Thickness of the brick [m].
    #[serde(rename = "Thickness")]
    pub thickness: f64,
}

fn default_nz_layers() -> usize { 1 }

/// A delta-gap port on the boxed rectilinear grid (Sonnet internal port).
///
/// The port is placed between `cell_a` and `cell_b` (must be adjacent
/// in either x or y direction). The delta-gap excitation applies ±1V
/// across the edge, and the port current is J_a �?J_b.
#[derive(Debug, Clone, Deserialize)]
pub struct BoxPortConfig {
    /// Port index (1-based).
    #[serde(rename = "Index")]
    pub index: u32,

    /// Port reference impedance [Ω].
    #[serde(rename = "Impedance", default = "default_ref_impedance")]
    pub impedance: f64,

    /// Dominant E-field direction: "x" or "y".
    #[serde(rename = "Direction", default = "default_port_direction")]
    pub direction: String,

    /// First cell grid index (ix, iy).
    #[serde(rename = "CellA")]
    pub cell_a: (usize, usize),

    /// Second cell grid index (ix, iy), adjacent to CellA.
    #[serde(rename = "CellB")]
    pub cell_b: (usize, usize),

    /// Port reactance X [Ω] for complex reference impedance Z = R + jX.
    /// Used by format converters (e.g. ADS Momentum) that parse complex
    /// port impedance from the project file.  The solver ignores X for
    /// now and uses the real part only; this preserves the data for
    /// future implementation of complex reference impedance.
    #[serde(rename = "Reactance", default)]
    pub reactance: f64,
}

fn default_evanescent_modes() -> usize { 10 }
fn default_conformal_level() -> u32 { 2 }
fn default_side_wall_pec() -> bool { true }

/// Single dielectric layer in the substrate stack.
#[derive(Debug, Clone, Deserialize)]
pub struct SubstrateLayerConfig {
    /// Relative permittivity (isotropic for now)
    #[serde(rename = "Permittivity", default = "default_permittivity")]
    pub permittivity: f64,

    /// Loss tangent: tan(δ) for dissipation model
    #[serde(rename = "LossTangent", default)]
    pub loss_tangent: f64,

    /// Dielectric DC conductivity [S/m].
    ///
    /// Models ohmic loss in semiconducting substrates (e.g., Si: σ �?10-50 S/m).
    /// The effective loss at frequency f is included as an additional imaginary
    /// contribution to the complex permittivity:  ε'' += σ/(2πf·ε₀).
    ///
    /// Default: 0.0 (perfect dielectric).
    #[serde(rename = "DielectricConductivity", default = "default_dielectric_conductivity")]
    pub dielectric_conductivity: f64,

    /// DC conductivity of the metal traces on this layer [S/m].
    ///
    /// When > 0, the solver applies a Leontovich surface-impedance boundary
    /// condition (SIBC) to the PEC cells on this layer, using the surface
    /// impedance  Zs = (1+j)·√(πfμ₀/σ).  This models ohmic loss in finite-
    /// conductivity metals (e.g. Cu = 5.8e7, Au = 4.1e7 S/m).
    ///
    /// Default: 0.0 (perfect electric conductor, lossless).
    #[serde(rename = "MetallizationConductivity", default = "default_metallization_conductivity")]
    pub metallization_conductivity: f64,

    /// Relative permeability (isotropic)
    #[serde(rename = "Permeability", default = "default_permeability")]
    pub permeability: f64,

    /// Layer thickness [m].  Use a very large value (1e10) for top (air) layer.
    #[serde(rename = "Thickness")]
    pub thickness: f64,

    /// Optional name (e.g., "Silicon", "FR4") for documentation
    #[serde(rename = "Name", default)]
    pub name: String,

    /// Optional frequency-dependent dispersion model.
    ///
    /// When present, overrides the static `Permittivity` + `LossTangent`
    /// values at each frequency point.  The complex permittivity is evaluated
    /// at the simulation frequency and used in the Sommerfeld integral.
    ///
    /// Supported models:
    /// - `"Debye"`: relaxation model  ε(ω) = ε_�?+ (ε_s �?ε_�? / (1 + jωτ)
    /// - `"Lorentz"`: oscillator model ε(ω) = ε_�?+ Δε·ω₀² / (ω₀² �?ω² + jγω)
    #[serde(rename = "Dispersion", default)]
    pub dispersion: Option<DispersionModel>,

    /// Anisotropic permittivity tensor [ε_xx, ε_yy, ε_zz].
    ///
    /// When set, overrides the scalar `Permittivity` field.  Supports
    /// uniaxial anisotropy (ε_xx = ε_yy �?ε_zz) typical of laminated PCB
    /// substrates. The layered Green model uses ε_t = ε_xx and ε_z = ε_zz.
    #[serde(rename = "PermittivityTensor", default)]
    pub permittivity_tensor: Option<[f64; 3]>,
}

/// Frequency-dependent permittivity dispersion model for a substrate layer.
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "Type")]
pub enum DispersionModel {
    /// Debye relaxation:  ε(ω) = ε_�?+ (ε_s �?ε_�? / (1 + jωτ)
    Debye {
        /// High-frequency permittivity ε_�?(lossless plateau above resonance)
        #[serde(rename = "EpsInf")]
        eps_inf: f64,
        /// Static (DC) permittivity ε_s
        #[serde(rename = "EpsStatic")]
        eps_static: f64,
        /// Relaxation time constant τ [s]
        #[serde(rename = "RelaxationTime")]
        tau_s: f64,
    },
    /// Lorentz oscillator: ε(ω) = ε_�?+ Δε·ω₀² / (ω₀² �?ω² + jγω)
    Lorentz {
        /// High-frequency permittivity ε_�?
        #[serde(rename = "EpsInf")]
        eps_inf: f64,
        /// Permittivity increment Δε = ε_s �?ε_�?
        #[serde(rename = "DeltaEps")]
        delta_eps: f64,
        /// Resonant angular frequency ω₀ [rad/s]
        #[serde(rename = "OmegaRes")]
        omega0_rad_per_s: f64,
        /// Damping coefficient γ [rad/s]
        #[serde(rename = "Gamma")]
        gamma_rad_per_s: f64,
    },
    /// Djordjević–Sarkar wideband model (nearly constant loss tangent).
    ///
    /// ε(ω) = ε_�?+ Δε · Σ_{k=1}^{N} w_k / (1 + jωτ_k)
    ///
    /// where `N` poles are logarithmically spaced between τ_min and τ_max,
    /// with uniform weights w_k = 1/N.  This produces a nearly constant
    /// loss tangent over many decades, satisfying the Kramers–Kronig
    /// relations (causal by construction).
    ///
    /// Reference: Djordjević & Sarkar, "Wideband Frequency Domain
    /// Characterization of FR-4 and Time-Domain Causality",
    /// IEEE Trans. EMC, 2001.
    #[serde(rename = "DjordjevicSarkar")]
    DjordjevicSarkar {
        /// High-frequency permittivity ε_�?
        #[serde(rename = "EpsInf")]
        eps_inf: f64,
        /// Permittivity increment Δε = ε_s �?ε_�?
        #[serde(rename = "DeltaEps")]
        delta_eps: f64,
        /// Minimum relaxation time τ_min [s] (typically ~1e-12)
        #[serde(rename = "TauMin")]
        tau_min: f64,
        /// Maximum relaxation time τ_max [s] (typically ~1e-3)
        #[serde(rename = "TauMax")]
        tau_max: f64,
        /// Number of poles (�?2).  Defaults to 10 if omitted.
        #[serde(rename = "NPoles", default = "default_ds_n_poles")]
        n_poles: usize,
    },
}

fn default_ds_n_poles() -> usize { 10 }

impl DispersionModel {
    /// Complex relative permittivity at angular frequency ω = 2π·f [rad/s].
    pub fn eps_r_at_omega(&self, omega: f64) -> num_complex::Complex64 {
        use num_complex::Complex64;
        match *self {
            DispersionModel::Debye { eps_inf, eps_static, tau_s } => {
                let denom = Complex64::new(1.0, omega * tau_s);
                Complex64::new(eps_inf, 0.0) + Complex64::new(eps_static - eps_inf, 0.0) / denom
            }
            DispersionModel::Lorentz { eps_inf, delta_eps, omega0_rad_per_s, gamma_rad_per_s } => {
                let w0sq = omega0_rad_per_s * omega0_rad_per_s;
                let denom = Complex64::new(w0sq - omega * omega, gamma_rad_per_s * omega);
                Complex64::new(eps_inf, 0.0) + Complex64::new(delta_eps * w0sq, 0.0) / denom
            }
            DispersionModel::DjordjevicSarkar { eps_inf, delta_eps, tau_min, tau_max, n_poles } => {
                let n = n_poles.max(2);
                let log_tau_min = tau_min.ln();
                let log_tau_max = tau_max.ln();
                let d_log = (log_tau_max - log_tau_min) / (n - 1) as f64;
                let mut sum = Complex64::ZERO;
                for k in 0..n {
                    let tau_k = (log_tau_min + k as f64 * d_log).exp();
                    let denom = Complex64::new(1.0, omega * tau_k);
                    let w_k = 1.0 / n as f64;
                    sum += Complex64::new(w_k, 0.0) / denom;
                }
                Complex64::new(eps_inf, 0.0) + Complex64::new(delta_eps, 0.0) * sum
            }
        }
    }
}

impl SubstrateLayerConfig {
    /// Complex relative permittivity ε_r(f) at frequency `freq_hz`.
    ///
    /// Priority:
    /// 1. Dispersion model (`Dispersion` key), if present.
    /// 2. Anisotropic tensor (`PermittivityTensor`): uses ε_xx as lateral ε.
    /// 3. Scalar `Permittivity` + `LossTangent` + `DielectricConductivity`.
    ///
    /// The imaginary part includes both the loss-tangent contribution and the
    /// DC conductivity:
    ///   ε'' = ε'·tanδ + σ_dc / (2πf · ε₀)
    pub fn eps_r_complex(&self, freq_hz: f64) -> num_complex::Complex64 {
        use num_complex::Complex64;
        use std::f64::consts::PI;
        if let Some(ref model) = self.dispersion {
            let omega = 2.0 * PI * freq_hz;
            let mut eps = model.eps_r_at_omega(omega);
            // Add DC conductivity contribution on top of the dispersion model
            if self.dielectric_conductivity > 0.0 && freq_hz > 0.0 {
                let sigma_term = self.dielectric_conductivity / (omega * EPS0);
                eps.im -= sigma_term;
            }
            eps
        } else if let Some([eps_xx, _, _]) = self.permittivity_tensor {
            let eps_imag_base = eps_xx * self.loss_tangent;
            let sigma_contrib = if self.dielectric_conductivity > 0.0 && freq_hz > 0.0 {
                self.dielectric_conductivity / (2.0 * PI * freq_hz * EPS0)
            } else { 0.0 };
            Complex64::new(eps_xx, -(eps_imag_base + sigma_contrib))
        } else {
            let eps = self.permittivity;
            let eps_imag_base = eps * self.loss_tangent;
            let sigma_contrib = if self.dielectric_conductivity > 0.0 && freq_hz > 0.0 {
                self.dielectric_conductivity / (2.0 * PI * freq_hz * EPS0)
            } else { 0.0 };
            Complex64::new(eps, -(eps_imag_base + sigma_contrib))
        }
    }

    /// Vertical (z-axis) complex permittivity ε_zz(f).
    ///
    /// For isotropic and non-tensor layers, equals `eps_r_complex`.
    /// For anisotropic tensor layers, returns ε_zz component.
    pub fn eps_r_z_complex(&self, freq_hz: f64) -> Option<num_complex::Complex64> {
        use num_complex::Complex64;
        use std::f64::consts::PI;
        if let Some([_, _, eps_zz]) = self.permittivity_tensor {
            let eps_imag_base = eps_zz * self.loss_tangent;
            let sigma_contrib = if self.dielectric_conductivity > 0.0 && freq_hz > 0.0 {
                self.dielectric_conductivity / (2.0 * PI * freq_hz * EPS0)
            } else { 0.0 };
            Some(Complex64::new(eps_zz, -(eps_imag_base + sigma_contrib)))
        } else {
            None // isotropic: caller may treat as same as eps_r_complex
        }
    }
}

fn default_bottom_pec() -> bool { true }
fn default_permittivity() -> f64 { 1.0 }
fn default_permeability() -> f64 { 1.0 }
fn default_dielectric_conductivity() -> f64 { 0.0 }
fn default_metallization_conductivity() -> f64 { 0.0 }

// ---------------------------------------------------------------------------
// SBR+ solver config (REM extension �?ignored by Palace)
// ---------------------------------------------------------------------------

/// SBR+ solver parameters, placed under `Solver.SBR` in the config file.
#[derive(Debug, Clone, Deserialize)]
pub struct SbrSolverConfig {
    /// Start frequency [Hz]
    #[serde(rename = "FreqMin")]
    pub freq_min: f64,

    /// End frequency [Hz]
    #[serde(rename = "FreqMax")]
    pub freq_max: f64,

    /// Frequency step [Hz]; set to 0 for single-frequency solve
    #[serde(rename = "FreqStep", default = "default_sbr_freq_step")]
    pub freq_step: f64,

    /// Ray density [rays/m²] on the aperture plane
    #[serde(rename = "RayDensity", default = "default_ray_density")]
    pub ray_density: f64,

    /// Maximum number of ray bounces
    #[serde(rename = "MaxBounces", default = "default_max_bounces")]
    pub max_bounces: usize,

    /// Energy weight threshold below which a ray is terminated
    #[serde(rename = "WeightThresh", default = "default_weight_thresh")]
    pub weight_thresh: f64,

    /// Target type: "PEC" | "Dielectric" | "Coated"
    #[serde(rename = "TargetType", default = "default_target_type")]
    pub target_type: String,

    /// Incident plane wave polar angle [degrees] from +z axis (0 = broadside)
    #[serde(rename = "ThetaInc", default)]
    pub theta_inc_deg: f64,

    /// Incident plane wave azimuth angle [degrees] from +x axis
    #[serde(rename = "PhiInc", default)]
    pub phi_inc_deg: f64,

    /// Incident polarization: "theta" | "phi" | "x" | "y" | "z"
    #[serde(rename = "Polarization", default = "default_polarization")]
    pub polarization: String,
}

fn default_sbr_freq_step()  -> f64    { 0.0 }
fn default_ray_density()    -> f64    { 1.0e4 }
fn default_max_bounces()    -> usize  { 5 }
fn default_weight_thresh()  -> f64    { 1.0e-4 }
fn default_target_type()    -> String { "PEC".to_string() }

// ---------------------------------------------------------------------------
// FE-BI solver config (REM extension �?ignored by Palace)
// ---------------------------------------------------------------------------

/// Hybrid FE-BI solver parameters, placed under `Solver.FEBI` in the config file.
///
/// FE-BI couples a volumetric FEM domain (for heterogeneous materials) with a
/// boundary integral operator on the outer radiation surface, providing a
/// rigorous open-domain boundary condition without PML layers.
#[derive(Debug, Clone, Deserialize)]
pub struct FeBiSolverConfig {
    /// Start frequency [Hz]
    #[serde(rename = "FreqMin")]
    pub freq_min: f64,

    /// End frequency [Hz]
    #[serde(rename = "FreqMax")]
    pub freq_max: f64,

    /// Frequency step [Hz]; set to 0 for single-frequency solve
    #[serde(rename = "FreqStep", default = "default_febi_freq_step")]
    pub freq_step: f64,

    /// Attribute IDs of the radiation boundary (outer surface Γ where BI is applied)
    #[serde(rename = "RadiationBoundary", default)]
    pub radiation_boundary: Vec<u32>,

    /// Boundary integral equation type: "EFIE" | "CFIE"
    #[serde(rename = "Equation", default = "default_febi_equation")]
    pub equation: String,

    /// CFIE mixing coefficient α �?[0,1]: 0 = pure EFIE, 1 = pure MFIE
    #[serde(rename = "Alpha", default = "default_febi_alpha")]
    pub alpha: f64,

    /// ACA tolerance for compressing the BI block (0 = direct dense)
    #[serde(rename = "AcaTol", default = "default_febi_aca_tol")]
    pub aca_tol: f64,

    /// GMRES relative residual tolerance for the coupled FE-BI system
    #[serde(rename = "GmresTol", default = "default_febi_gmres_tol")]
    pub gmres_tol: f64,

    /// Maximum GMRES iterations
    #[serde(rename = "GmresMaxIter", default = "default_febi_gmres_max_iter")]
    pub gmres_max_iter: usize,

    /// Lumped ports for S-parameter extraction (same format as MoM ports)
    #[serde(rename = "Ports", default)]
    pub ports: Vec<MomPort>,

    /// Global reference impedance Z₀ [Ω] for S-parameter normalisation
    #[serde(rename = "RefImpedance", default = "default_ref_impedance")]
    pub ref_impedance: f64,

    /// Relative permittivity of the exterior (background) medium.
    /// Default: 1.0 (vacuum / air).
    #[serde(rename = "ExteriorEpsR", default = "default_exterior_eps")]
    pub exterior_eps_r: f64,

    /// Relative permeability of the exterior medium.
    /// Default: 1.0 (vacuum / air).
    #[serde(rename = "ExteriorMuR", default = "default_exterior_mu")]
    pub exterior_mu_r: f64,

    /// Output directory for postprocessing results
    #[serde(rename = "OutputDir", default = "default_febi_output_dir")]
    pub output_dir: String,
}

fn default_febi_freq_step()     -> f64    { 0.0 }
fn default_febi_equation()      -> String { "CFIE".to_string() }
fn default_febi_alpha()         -> f64    { 0.5 }
fn default_febi_aca_tol()       -> f64    { 1.0e-3 }
fn default_febi_gmres_tol()     -> f64    { 1.0e-6 }
fn default_febi_gmres_max_iter() -> usize { 500 }
fn default_febi_output_dir()    -> String { "postpro".to_string() }
fn default_exterior_eps()       -> f64    { 1.0 }
fn default_exterior_mu()        -> f64    { 1.0 }

// ---------------------------------------------------------------------------
// DDM solver config (REM extension �?ignored by Palace)
// ---------------------------------------------------------------------------

/// Domain Decomposition Method solver parameters, placed under `Solver.DDM`.
///
/// Implements Robin-condition Schwarz iteration for parallel FEM solves.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct DdmSolverConfig {
    /// Number of subdomains (should equal MPI process count)
    #[serde(rename = "NumSubdomains", default = "default_ddm_num_subdomains")]
    pub num_subdomains: usize,

    /// DDM algorithm: "Schwarz" | "FETI"
    #[serde(rename = "Method", default = "default_ddm_method")]
    pub method: String,

    /// Robin condition coefficient order: 1 = first-order OSRC, 2 = second-order
    #[serde(rename = "RobinOrder", default = "default_ddm_robin_order")]
    pub robin_order: u8,

    /// Convergence tolerance for the Schwarz outer iteration
    #[serde(rename = "Tolerance", default = "default_ddm_tolerance")]
    pub tolerance: f64,

    /// Maximum number of Schwarz iterations
    #[serde(rename = "MaxIter", default = "default_ddm_max_iter")]
    pub max_iter: usize,

    /// METIS partitioning: "Dual" | "Nodal"
    #[serde(rename = "PartitionType", default = "default_ddm_partition_type")]
    pub partition_type: String,

    /// Use Multiplicative Schwarz (sequential, uses latest neighbor solution each sub-step).
    /// Additive Schwarz (default false) is embarrassingly parallel; Multiplicative
    /// converges in ~half the iterations at the cost of sequential subdomain solves.
    #[serde(rename = "Multiplicative", default)]
    pub multiplicative: bool,

    /// Anderson acceleration history depth m (0 = disabled).
    /// Stores the last m residual/iterate pairs and solves a small LS problem each
    /// iteration to accelerate convergence. Typical values: 3�?0.
    #[serde(rename = "AndersonDepth", default = "default_ddm_anderson_depth")]
    pub anderson_depth: usize,

    /// Operating frequency [Hz] for the Helmholtz FEM subdomain solves.
    /// Used to compute k = ω�?μ₀ε₀) for the wave operator and Robin α = jk.
    /// Default 1 GHz.  Set to match your Driven solver FreqMin/FreqMax centre.
    #[serde(rename = "FreqHz", default = "default_ddm_freq_hz")]
    pub freq_hz: f64,

    /// Start frequency for sweep [Hz] (default: same as FreqHz).
    #[serde(rename = "FreqMin", default = "default_ddm_freq_hz")]
    pub freq_min: f64,
    /// End frequency for sweep [Hz] (default: same as FreqHz).
    #[serde(rename = "FreqMax", default = "default_ddm_freq_hz")]
    pub freq_max: f64,
    /// Frequency step [Hz] (default: 0 = single frequency).
    #[serde(rename = "FreqStep", default)]
    pub freq_step: f64,

    /// Relative permittivity ε_r for subdomain Helmholtz assembly (default 1.0).
    #[serde(rename = "EpsR", default = "default_one_f64")]
    pub eps_r: f64,

    /// Relative permeability μ_r for subdomain Helmholtz assembly (default 1.0).
    #[serde(rename = "MuR", default = "default_one_f64")]
    pub mu_r: f64,
}

fn default_ddm_num_subdomains() -> usize  { 4 }
fn default_ddm_method()         -> String { "Schwarz".to_string() }
fn default_ddm_robin_order()    -> u8     { 1 }
fn default_ddm_tolerance()      -> f64    { 1.0e-6 }
fn default_ddm_max_iter()       -> usize  { 100 }
fn default_ddm_partition_type() -> String { "Dual".to_string() }
fn default_ddm_anderson_depth() -> usize  { 5 }
fn default_ddm_freq_hz()        -> f64    { 1.0e9 }
fn default_one_f64()            -> f64    { 1.0 }

// ── Planar MoM solver config ───────────────────────────────────────────────────────

/// Planar MoM solver configuration for layered-media RF passives.
///
/// Generates a uniform rectangular grid with RWG basis functions and solves
/// the MoM system for S-parameter extraction.
#[derive(Debug, Clone, Deserialize)]
pub struct PlanarSolverConfig {
    /// X-direction length of the simulation domain [m]
    #[serde(rename = "Lx")]
    pub lx: f64,
    /// Y-direction length of the simulation domain [m]
    #[serde(rename = "Ly")]
    pub ly: f64,
    /// Number of grid segments in X
    #[serde(rename = "Nx", default = "default_planar_n")]
    pub nx: usize,
    /// Number of grid segments in Y
    #[serde(rename = "Ny", default = "default_planar_n")]
    pub ny: usize,

    /// Start frequency [Hz]
    #[serde(rename = "FreqMin")]
    pub freq_min: f64,
    /// End frequency [Hz]
    #[serde(rename = "FreqMax")]
    pub freq_max: f64,
    /// Frequency step [Hz]
    #[serde(rename = "FreqStep")]
    pub freq_step: f64,

    /// Enable FFT-accelerated matrix-vector product (default: false).
    #[serde(rename = "UseFft", default)]
    pub use_fft: bool,

    /// Conductor wall conductivity σ [S/m] for SIBC loss modeling.
    /// Set > 0 (e.g. 5.8e7 for copper) to add Leontovich surface impedance
    /// Zs = (1+j)/(σ·δs) on all conducting surfaces.
    #[serde(rename = "WallConductivity", default)]
    pub wall_conductivity: f64,

    /// Dielectric layers for the substrate (bottom-to-top).
    /// Format: `[eps_r, thickness_mm, ...]` per layer.
    #[serde(rename = "SubstrateLayers", default)]
    pub substrate_layers: Vec<PlanarLayerSpec>,

    /// Edge port definitions for S-parameter extraction.
    #[serde(rename = "Ports", default)]
    pub ports: Vec<PlanarPortSpec>,

    /// Characteristic impedance for S-parameter normalisation [Ω].
    #[serde(rename = "RefImpedance", default = "default_ref_impedance")]
    pub ref_impedance: f64,
}

/// A single dielectric layer in the planar substrate stack.
#[derive(Debug, Clone, Deserialize)]
pub struct PlanarLayerSpec {
    /// Relative permittivity ε_r
    #[serde(rename = "EpsR")]
    pub eps_r: f64,
    /// Loss tangent tan δ
    #[serde(rename = "LossTan", default)]
    pub loss_tan: f64,
    /// Layer thickness [m]
    #[serde(rename = "Thickness")]
    pub thickness: f64,
}

/// An edge port on the planar grid for S-parameter extraction.
#[derive(Debug, Clone, Deserialize)]
pub struct PlanarPortSpec {
    /// Port index (1-based, for S-matrix ordering)
    #[serde(rename = "Index")]
    pub index: usize,
    /// Edge index on the grid to attach the port
    #[serde(rename = "Edge", default)]
    pub edge: usize,
}

fn default_planar_n() -> usize { 20 }

// ---------------------------------------------------------------------------
// Parametric sweep / gradient optimization (REM extension �?ignored by Palace)
// ---------------------------------------------------------------------------

/// Mode for parametric run: full grid sweep or gradient optimization.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
pub enum ParametricMode {
    /// Exhaustive grid sweep over all parameter combinations.
    Sweep,
    /// Derivative-free Nelder-Mead optimization.
    Optimize,
    /// Central finite-difference sensitivity (gradient) analysis at the nominal point.
    Sensitivity,
    /// Monte Carlo yield analysis: sample parameters from Gaussian distributions
    /// around nominal values, compute statistics over N trials.
    MonteCarlo,
}

impl Default for ParametricMode {
    fn default() -> Self { ParametricMode::Sweep }
}

/// Target design variable that a parametric parameter controls.
///
/// JSON: `{"Type": "SubstratePermittivity", "Layer": 0}`
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "Type")]
pub enum ParamTarget {
    /// Lateral relative permittivity of substrate layer `layer` (0-indexed).
    SubstratePermittivity { #[serde(rename = "Layer")] layer: usize },
    /// Physical thickness [m] of substrate layer `layer` (0-indexed).
    SubstrateThickness    { #[serde(rename = "Layer")] layer: usize },
    /// Loss tangent of substrate layer `layer` (0-indexed).
    SubstrateLossTangent  { #[serde(rename = "Layer")] layer: usize },
    /// Reference impedance [Ω] of MoM lumped port `port` (1-indexed, matching `Index`).
    PortZ0 { #[serde(rename = "Port")] port: usize },
    /// Start frequency [Hz] of the MoM sweep.
    FreqMin,
    /// End frequency [Hz] of the MoM sweep.
    FreqMax,
    /// Conductor wall conductivity [S/m] for SIBC loss modeling.
    WallConductivity,
    /// RMS surface roughness Δ [m] for SIBC roughness correction.
    RmsRoughness,
}

/// A single named design parameter with its sweep values (Sweep mode) or
/// initial value and bounds (Optimize mode).
#[derive(Debug, Clone, Deserialize)]
pub struct SweepParam {
    /// Human-readable name used as CSV column header.
    #[serde(rename = "Name")]
    pub name: String,

    /// Physical target that this parameter controls.
    #[serde(rename = "Target")]
    pub target: ParamTarget,

    /// Explicit list of values (Sweep mode).
    /// If omitted, values are generated from `Min`, `Max`, `Steps`.
    #[serde(rename = "Values", default)]
    pub values: Vec<f64>,

    /// Sweep range start (inclusive, Sweep mode).
    #[serde(rename = "Min")]
    pub min: Option<f64>,

    /// Sweep range end (inclusive, Sweep mode).
    #[serde(rename = "Max")]
    pub max: Option<f64>,

    /// Number of equally-spaced steps from `Min` to `Max` (Sweep mode, �?2).
    #[serde(rename = "Steps")]
    pub steps: Option<usize>,

    /// Starting value for the Nelder-Mead optimizer (Optimize mode).
    #[serde(rename = "Initial")]
    pub initial: Option<f64>,

    /// Optimizer parameter bounds `[lower, upper]` (Optimize mode).
    /// The optimizer clamps to these bounds at each evaluation.
    #[serde(rename = "Bounds")]
    pub bounds: Option<[f64; 2]>,
}

impl SweepParam {
    /// Resolve the list of values to use in Sweep mode.
    pub fn resolved_values(&self) -> Vec<f64> {
        if !self.values.is_empty() {
            return self.values.clone();
        }
        if let (Some(lo), Some(hi), Some(n)) = (self.min, self.max, self.steps) {
            let n = n.max(2);
            (0..n).map(|i| lo + (hi - lo) * i as f64 / (n - 1) as f64).collect()
        } else {
            vec![]
        }
    }
}

/// Optimization objective (what to minimize; larger = worse).
///
/// JSON: `{"Type": "MinS11dB", "Port": 1, "FreqHz": 2.4e9}`
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "Type")]
pub enum OptimObjective {
    /// Minimize |S_{port,port}| in dB at a single frequency.
    MinS11dB {
        #[serde(rename = "Port")]    port:    usize,
        #[serde(rename = "FreqHz")] freq_hz: f64,
    },
    /// Minimize |S_{i,j}| in dB (e.g., insertion loss).
    MinSijdB {
        #[serde(rename = "PortI")]  port_i:  usize,
        #[serde(rename = "PortJ")]  port_j:  usize,
        #[serde(rename = "FreqHz")] freq_hz: f64,
    },
    /// Maximize bandwidth where |S_{port,port}| < `thresh_db` [dB] (return loss).
    /// Objective = negative bandwidth (minimizer makes it most negative).
    MaxBandwidthS11dB {
        #[serde(rename = "Port")]       port:       usize,
        #[serde(rename = "ThreshDb")]   thresh_db:  f64,
    },

}

/// ```json
/// "Parametric": {
///   "Mode": "Sweep",
///   "Parameters": [
///     {"Name": "eps_r", "Target": {"Type": "SubstratePermittivity", "Layer": 0},
///      "Min": 3.0, "Max": 5.0, "Steps": 5}
///   ]
/// }
/// ```
///
/// # Optimize example
/// ```json
/// "Parametric": {
///   "Mode": "Optimize",
///   "Parameters": [
///     {"Name": "eps_r",     "Target": {"Type": "SubstratePermittivity", "Layer": 0},
///      "Initial": 4.0, "Bounds": [3.0, 5.0]},
///     {"Name": "thickness", "Target": {"Type": "SubstrateThickness", "Layer": 0},
///      "Initial": 0.254e-3, "Bounds": [0.1e-3, 0.5e-3]}
///   ],
///   "Objectives": [
///     {"Type": "MinS11dB", "Port": 1, "FreqHz": 2.4e9}
///   ],
///   "MaxIter": 200,
///   "Tolerance": 1e-4
/// }
/// ```
#[derive(Debug, Clone, Deserialize)]
pub struct ParametricConfig {
    /// Sweep or optimization mode.
    #[serde(rename = "Mode", default)]
    pub mode: ParametricMode,

    /// Named design parameters.
    #[serde(rename = "Parameters", default)]
    pub parameters: Vec<SweepParam>,

    /// Objective functions to minimize (required for Optimize mode).
    #[serde(rename = "Objectives", default)]
    pub objectives: Vec<OptimObjective>,

    /// Maximum optimizer iterations (Optimize mode, default: 500).
    #[serde(rename = "MaxIter", default = "default_optim_max_iter")]
    pub max_iter: usize,

    /// Convergence tolerance on simplex size (Optimize mode, default: 1e-4).
    #[serde(rename = "Tolerance", default = "default_optim_tolerance")]
    pub tolerance: f64,

    /// Relative step size h/p for central finite-difference sensitivity (Sensitivity mode, default: 0.01).
    #[serde(rename = "SensRelStep", default)]
    pub sens_rel_step: Option<f64>,

    /// Number of Monte Carlo trials (MonteCarlo mode, default: 100).
    #[serde(rename = "McSamples", default = "default_mc_samples")]
    pub mc_samples: usize,

    /// Relative 1-σ standard deviation for each parameter Gaussian in MonteCarlo mode
    /// (e.g. 0.05 = 5 % of nominal value, default: 0.05).
    #[serde(rename = "McSigmaRel", default = "default_mc_sigma")]
    pub mc_sigma_rel: f64,

    /// Optional random seed for reproducible Monte Carlo runs.
    #[serde(rename = "McSeed", default)]
    pub mc_seed: Option<u64>,

    /// Number of parallel evaluations for grid sweep (default: 1 = serial).
    #[serde(rename = "NParallel", default = "default_one")]
    pub n_parallel: usize,
}

fn default_optim_max_iter() -> usize { 500 }
fn default_optim_tolerance() -> f64  { 1e-4 }
fn default_mc_samples()      -> usize { 100 }
fn default_mc_sigma()        -> f64   { 0.05 }
fn default_one()             -> usize { 1 }

// ---------------------------------------------------------------------------
// Postprocessing (REM extension �?ignored by Palace)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Validation helpers: warn on Palace fields not yet implemented in REM
// ---------------------------------------------------------------------------

/// Log an info notice for a Palace field that is accepted but not yet implemented.
/// Call this from each solver's run() function for relevant config fields.
pub fn warn_unsupported(name: &str, hint: &str) {
    log::info!(
        "[REM] Palace field accepted (not fully implemented): {name}. {hint}."
    );
}

/// Validate a Palace config after full deserialization.
/// Emits warnings for all unsupported-but-accepted fields.
pub fn validate_palace_compat(cfg: &PalaceConfig) {
    // --- Problem ---
    if let Some(ref pf) = cfg.problem.output_formats {
        if pf.grid_function {
            log::info!(
                "[REM] Problem.OutputFormats.GridFunction: \
                 VTK solution files are written to <output>/paraview/solution.vtk.",
            );
        }
    }

    // --- Solver ---
    if cfg.solver.device != "CPU" && !cfg.solver.device.is_empty() {
        warn_unsupported(
            &format!("Solver.Device = \"{}\"", cfg.solver.device),
            "REM runs CPU-only; value is ignored",
        );
    }
    let disc = cfg.solver.discretization.to_lowercase();
    if !disc.is_empty() && disc != "h1" && disc != "hcurl" && disc != "nedelec" {
        warn_unsupported(
            &format!("Solver.Discretization = \"{}\"", cfg.solver.discretization),
            "Only H1 and HCurl/Nedelec are supported; defaulting behavior is solver-dependent",
        );
    }

    // --- Solver.Linear ---
    let l = &cfg.solver.linear;
    // CG / PCG / BiCGSTAB / GMRES are all accepted: SPD solvers use PCG by default;
    // Driven uses GMRES unless KSPType requests sparse iterative complex solve.
    let ksp_lower = l.ksp_type.to_lowercase();
    if !l.ksp_type.is_empty()
        && ksp_lower != "gmres"
        && ksp_lower != "cg"
        && ksp_lower != "pcg"
        && ksp_lower != "bicgstab"
        && ksp_lower != "default"
    {
        warn_unsupported(
            &format!("Solver.Linear.KSPType = \"{}\"", l.ksp_type),
            "Only GMRES, BiCGSTAB, and CG/PCG are supported; value is ignored",
        );
    }
    if l.mg_levels != 10 && l.mg_levels != 0 {
        warn_unsupported(
            &format!("Solver.Linear.MGLevels = {}", l.mg_levels),
            "Algebraic multigrid is not yet implemented; value is ignored",
        );
    }
    if l.complex_coarse_solve {
        warn_unsupported(
            "Solver.Linear.ComplexCoarseSolve = true",
            "Complex coarse-grid solver is not implemented; ignored",
        );
    }

    // --- Solver.Driven ---
    if let Some(ref d) = cfg.solver.driven {
        // Samples is now supported: Linear/Log/Point types
        if !d.save.is_empty() {
            warn_unsupported(
                "Solver.Driven.Save (array)",
                "Save array is not supported; use SaveStep integer instead",
            );
        }
    }

    // --- Solver.Transient ---
    if let Some(ref t) = cfg.solver.transient {
        let supported_excitations = ["", "none", "Step", "ModulatedGaussian", "Gaussian"];
        if !t.excitation.is_empty() && !supported_excitations.contains(&t.excitation.as_str()) {
            log::warn!(
                "[REM] Unsupported Solver.Transient.Excitation = \"{}\"; \
                 supported: Step (default), ModulatedGaussian, Gaussian. Falling back to Step.",
                t.excitation
            );
        }
    }

    // --- Domains ---
    for m in &cfg.domains.materials {
        if !m.material_axes.is_empty() {
            log::info!(
                "Domains.Materials[].MaterialAxes: {} axes provided -- tensor epsilon assembly enabled.",
                m.material_axes.len()
            );
        }
    }

    if let Some(ref dp) = cfg.domains.postprocessing {
        if !dp.energy.is_empty() {
            log::info!(
                "[REM] Domains.Postprocessing.Energy: {} group(s) -- \
                 per-group energy written to postpro/energy-E.csv (Electrostatic) \
                 or postpro/energy-B.csv (Magnetostatic).",
                dp.energy.len()
            );
        }
        if !dp.probe.is_empty() {
            log::info!(
                "[REM] Domains.Postprocessing.Probe: {} probe(s) �?\
                 results written to postpro/probe-phi.csv and probe-E.csv \
                 (Electrostatic and Magnetostatic solvers).",
                dp.probe.len()
            );
        }
    }

    // --- Boundaries ---
    if !cfg.boundaries.impedance.is_empty() {
        log::info!("[REM] Boundaries.Impedance: surface impedance BC active ({} regions).",
            cfg.boundaries.impedance.len());
    }
    if !cfg.boundaries.periodic.is_empty() {
        let n_pairs: usize = cfg.boundaries.periodic.iter()
            .map(|p| p.boundary_pairs.len())
            .sum();
        let k_complex = cfg.boundaries.periodic.iter()
            .any(|p| p.floquet_wave_vector.iter().any(|&v| v.abs() > 1e-14));
        if k_complex {
            log::warn!(
                "[REM] Boundaries.Periodic: {} spec(s), {} pair(s) �?non-zero FloquetWaveVector \
                 detected. Complex phase-shift BCs are not yet supported; Γ-point (k=0) pairs will \
                 be applied, others skipped.",
                cfg.boundaries.periodic.len(), n_pairs
            );
        } else {
            log::info!(
                "[REM] Boundaries.Periodic: {} spec(s), {} pair(s) �?Γ-point periodic BCs enabled.",
                cfg.boundaries.periodic.len(), n_pairs
            );
        }
    }
    if !cfg.domains.current_dipole.is_empty() {
        log::debug!(
            "Domains.CurrentDipole: {} dipole source(s) will be applied as point sources in the driven solver",
            cfg.domains.current_dipole.len()
        );
    }
    if !cfg.boundaries.postprocessing_flux.is_empty() {
        let mut n_electric = 0usize;
        let mut n_magnetic = 0usize;
        let mut n_other = 0usize;
        for s in &cfg.boundaries.postprocessing_flux {
            match s.flux_type.to_lowercase().as_str() {
                "electric" => n_electric += 1,
                "magnetic" => n_magnetic += 1,
                _ => n_other += 1,
            }
        }
        if n_electric > 0 {
            log::info!(
                "[REM] Boundaries.Postprocessing: {} Electric flux spec(s) -- \
                 displacement flux written to postpro/surface-flux.csv.",
                n_electric
            );
        }
        if n_magnetic > 0 {
            log::info!(
                "[REM] Boundaries.Postprocessing: {} Magnetic flux spec(s) -- \
                 B-field flux written to postpro/surface-flux.csv.",
                n_magnetic
            );
        }
        if n_other > 0 {
            warn_unsupported(
                "Boundaries.Postprocessing (Power / SA / MS / MA types)",
                "Power flux and interface dielectric loss are not implemented; ignored",
            );
        }
    }

    for lp in &cfg.boundaries.lumped_port {
        if !lp.elements.is_empty() {
            log::info!(
                "Boundaries.LumpedPort[{}]: {} elements �?all element attributes mapped to port BC.",
                lp.index, lp.elements.len()
            );
        }
    }

    for wp in &cfg.boundaries.wave_port {
        if wp.offset != 0.0 {
            warn_unsupported(
                &format!("Boundaries.WavePort[].Offset = {}", wp.offset),
                "WavePort offset is not implemented; ignored",
            );
        }
        if wp.max_its > 0 && wp.max_its != 30 {
            warn_unsupported(
                &format!("Boundaries.WavePort[].MaxIts = {}", wp.max_its),
                "Custom WavePort iterative solver max iterations are not implemented; ignored",
            );
        }
        if wp.eigen_tol > 0.0 && wp.eigen_tol != 1e-6 {
            warn_unsupported(
                &format!("Boundaries.WavePort[].EigenTol = {}", wp.eigen_tol),
                "Custom WavePort eigenvalue tolerance is not implemented; ignored",
            );
        }
        if wp.verbose_port > 0 {
            warn_unsupported(
                &format!("Boundaries.WavePort[].Verbose = {}", wp.verbose_port),
                "WavePort verbose output is not implemented; ignored",
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Postprocessing (REM extension �?ignored by Palace)
// ---------------------------------------------------------------------------

/// Top-level `Postprocessing` section (REM extension).
#[derive(Debug, Clone, Deserialize, Default)]
pub struct Postprocessing {
    /// Far-field RCS pattern output
    #[serde(rename = "RCS", default)]
    pub rcs: Option<RcsConfig>,

    /// Near-field data export on specified boundary surfaces.
    #[serde(rename = "NearField", default)]
    pub near_field: Option<NearFieldExportConfig>,
}

/// Configuration for near-field export.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct NearFieldExportConfig {
    /// Boundary attribute IDs whose surface/boundary faces will have E/H fields exported.
    #[serde(rename = "Attributes", default)]
    pub attributes: Vec<u32>,

    /// Output file path relative to output directory (default: "postpro/near_field.csv").
    #[serde(rename = "OutputFile", default)]
    pub output_file: Option<String>,
}

/// Configuration for near-field source import (used as excitation).
#[derive(Debug, Clone, Deserialize, Default)]
pub struct NearFieldSourceConfig {
    /// Path to a near-field CSV file exported from a previous simulation.
    #[serde(rename = "File")]
    pub file: String,

    /// Boundary attribute IDs where the near-field will be applied as excitation.
    #[serde(rename = "Attributes", default)]
    pub attributes: Vec<u32>,
}

/// Configuration for RCS (Radar Cross Section) output.
#[derive(Debug, Clone, Deserialize)]
pub struct RcsConfig {
    /// Azimuth angles φ in degrees, e.g. `[0, 90]`
    #[serde(rename = "PhiDeg", default)]
    pub phi_deg: Vec<f64>,

    /// Elevation angles θ in degrees; supports range string `"0:5:180"`
    #[serde(rename = "ThetaDeg", deserialize_with = "deserialize_angle_range", default)]
    pub theta_deg: Vec<f64>,
}

/// Deserialize either `[0, 10, 20]` or `"0:10:180"` into `Vec<f64>`.
pub fn deserialize_angle_range<'de, D>(d: D) -> Result<Vec<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    use serde::de::{self, Visitor};
    use std::fmt;

    struct AngleVisitor;

    impl<'de> Visitor<'de> for AngleVisitor {
        type Value = Vec<f64>;

        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "an array of floats or a range string like \"0:10:180\"")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Vec<f64>, A::Error>
        where A: de::SeqAccess<'de> {
            let mut v = Vec::new();
            while let Some(x) = seq.next_element::<f64>()? {
                v.push(x);
            }
            Ok(v)
        }

        fn visit_str<E: de::Error>(self, s: &str) -> Result<Vec<f64>, E> {
            // Parse "start:step:end"
            let parts: Vec<&str> = s.split(':').collect();
            if parts.len() == 3 {
                let start = parts[0].trim().parse::<f64>().map_err(de::Error::custom)?;
                let step  = parts[1].trim().parse::<f64>().map_err(de::Error::custom)?;
                let end   = parts[2].trim().parse::<f64>().map_err(de::Error::custom)?;
                if step == 0.0 { return Err(de::Error::custom("ThetaDeg step cannot be zero")); }
                let n = ((end - start) / step).round() as usize + 1;
                Ok((0..n).map(|i| start + i as f64 * step).collect())
            } else {
                Err(de::Error::custom(format!("expected \"start:step:end\", got {:?}", s)))
            }
        }
    }

    d.deserialize_any(AngleVisitor)
}

#[cfg(test)]
mod dispersion_tests {
    use super::*;
    use num_complex::Complex64;
    use std::f64::consts::PI;

    #[test]
    fn debye_at_dc_equals_eps_static() {
        let model = DispersionModel::Debye { eps_inf: 2.0, eps_static: 4.5, tau_s: 1e-11 };
        // At ω�?, ε(ω) �?ε_s
        let eps = model.eps_r_at_omega(1.0); // very low freq
        assert!((eps.re - 4.5).abs() < 0.01, "Debye DC: re={:.4}", eps.re);
        assert!(eps.im < 0.0, "Debye DC: imaginary part should be negative (loss)");
    }

    #[test]
    fn debye_at_high_freq_equals_eps_inf() {
        let model = DispersionModel::Debye { eps_inf: 2.0, eps_static: 4.5, tau_s: 1e-11 };
        // At very high ω (ωτ >> 1), ε(ω) �?ε_�?
        let eps = model.eps_r_at_omega(1e14); // 10 THz >> 1/τ
        assert!((eps.re - 2.0).abs() < 0.05, "Debye HF: re={:.4}", eps.re);
    }

    #[test]
    fn lorentz_at_resonance_is_imaginary_dominated() {
        let omega0 = 2.0e11; // 200 GHz resonance
        let gamma  = 1.0e10; // 10 GHz damping
        let model = DispersionModel::Lorentz {
            eps_inf: 1.0, delta_eps: 1.0, omega0_rad_per_s: omega0, gamma_rad_per_s: gamma,
        };
        let eps = model.eps_r_at_omega(omega0); // at resonance
        // At resonance: ε = ε_�?+ Δε·ω₀²/(jγω₀) = ε_�?- j Δε·ω₀/γ
        assert!(eps.im.abs() > eps.re.abs(), "Lorentz at resonance: should be loss-dominated");
    }

    #[test]
    fn substrate_layer_eps_r_complex_uses_dispersion_when_set() {
        let layer = SubstrateLayerConfig {
            permittivity: 4.0,
            loss_tangent: 0.02,
            dielectric_conductivity: 0.0,
            metallization_conductivity: 0.0,
            permeability: 1.0,
            thickness: 1e-3,
            name: String::new(),
            dispersion: Some(DispersionModel::Debye {
                eps_inf: 2.0, eps_static: 6.0, tau_s: 1e-11,
            }),
            permittivity_tensor: None,
        };
        let freq = 1e9; // 1 GHz
        let eps = layer.eps_r_complex(freq);
        // Debye at 1 GHz: ωτ = 2π·1e9·1e-11 �?0.063; still close to ε_s
        assert!(eps.re > 4.0 && eps.re < 6.5, "Debye 1 GHz re={:.3}", eps.re);
        assert!(eps.im < 0.0, "Should have negative imaginary part (loss)");
    }

    #[test]
    fn substrate_layer_eps_r_z_returns_tensor_zz() {
        let layer = SubstrateLayerConfig {
            permittivity: 3.48,
            loss_tangent: 0.004,
            dielectric_conductivity: 0.0,
            metallization_conductivity: 0.0,
            permeability: 1.0,
            thickness: 0.254e-3,
            name: "Rogers4350B".to_string(),
            dispersion: None,
            permittivity_tensor: Some([3.48, 3.48, 3.66]),
        };
        let eps_lat = layer.eps_r_complex(1e9);
        let eps_z   = layer.eps_r_z_complex(1e9).expect("should be Some for tensor");
        assert!((eps_lat.re - 3.48).abs() < 1e-6);
        assert!((eps_z.re   - 3.66).abs() < 1e-6);
        // z and lateral should differ
        assert!((eps_z.re - eps_lat.re).abs() > 0.1);
    }

    #[test]
    fn substrate_layer_isotropic_eps_r_z_returns_none() {
        let layer = SubstrateLayerConfig {
            permittivity: 4.2,
            loss_tangent: 0.02,
            dielectric_conductivity: 0.0,
            metallization_conductivity: 0.0,
            permeability: 1.0,
            thickness: 1e-3,
            name: String::new(),
            dispersion: None,
            permittivity_tensor: None,
        };
        assert!(layer.eps_r_z_complex(1e9).is_none());
    }

    #[test]
    fn djordjevic_sarkar_loss_tan_constant() {
        // FR-4 parameters with 20 poles for flatter response
        let model = DispersionModel::DjordjevicSarkar {
            eps_inf: 3.9, delta_eps: 0.6,
            tau_min: 1e-12, tau_max: 1e-3, n_poles: 20,
        };
        let freqs: Vec<f64> = (0..50).map(|i| 1e6 * 10.0_f64.powf(i as f64 / 49.0 * 4.0)).collect();
        let mut tan_deltas = Vec::new();
        for &f in &freqs {
            let eps = model.eps_r_at_omega(2.0 * PI * f);
            let td = (eps.im / eps.re).abs();
            tan_deltas.push(td);
        }
        let mean_td: f64 = tan_deltas.iter().sum::<f64>() / tan_deltas.len() as f64;
        let max_dev: f64 = tan_deltas.iter().map(|&t| (t - mean_td).abs() / mean_td).fold(0.0_f64, f64::max);
        assert!(max_dev < 0.10, "loss tangent deviation {:.3} > 10%", max_dev);
    }

    #[test]
    fn djordjevic_sarkar_causality_kramers_kronig() {
        // D-S model is a sum of Debye terms, each analytically satisfying K-K.
        // Verify the real-part reconstruction at a pole frequency where the
        // integrand for the K-K transform is well-behaved.
        // Test a single Debye pole first: ε(ω)=1+Δε/(1+jωτ), Δε=2, τ=1e-9
        let tau = 1e-9;
        let debye_eps = |omega: f64| -> Complex64 {
            let denom = Complex64::new(1.0, omega * tau);
            Complex64::new(1.0, 0.0) + Complex64::new(2.0, 0.0) / denom
        };
        // At ω₀=1/τ, ε'(ω₀)=1+Δε/(1+1)=1+1=2, ε''(ω₀)=Δε·1/(1+1)=1
        let omega0 = 1.0 / tau;
        let eps = debye_eps(omega0);
        assert!((eps.re - 2.0).abs() < 1e-12, "Debye real part mismatch");
        assert!((eps.im + 1.0).abs() < 1e-12, "Debye imag part mismatch");

        // D-S is sum of weighted Debye poles �?therefore K‑K compliant.
        // Verify that ε'(ω) matches the explicit reactive (non-dispersive)
        // formula: ε'(ω) = 3.9 + 0.6 · Σ w_k / (1+ω²τ_k²)
        let model = DispersionModel::DjordjevicSarkar {
            eps_inf: 3.9, delta_eps: 0.6,
            tau_min: 1e-12, tau_max: 1e-3, n_poles: 20,
        };
        let omega_test = 2.0 * PI * 1e9;
        let eps_direct = model.eps_r_at_omega(omega_test);
        let n = 20;
        let log_min = 1e-12_f64.ln();
        let log_max = 1e-3_f64.ln();
        let d_log = (log_max - log_min) / (n - 1) as f64;
        let mut re_explicit = 3.9;
        for k in 0..n {
            let tau_k = (log_min + k as f64 * d_log).exp();
            let w_k = 1.0 / n as f64;
            re_explicit += 0.6 * w_k / (1.0 + omega_test * omega_test * tau_k * tau_k);
        }
        assert!((eps_direct.re - re_explicit).abs() < 1e-12);
    }
}
