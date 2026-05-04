use crate::preprocess::expand_ranges;
use serde::{Deserialize, Deserializer};

/// Deserialize a JSON value that may be either a scalar `f64` or an array `[f64, ...]`.
/// When an array is given, the first element is used (anisotropic → isotropic fallback).
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
}

// ---------------------------------------------------------------------------
// Problem
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
pub struct Problem {
    #[serde(rename = "Type")]
    pub problem_type: ProblemType,

    /// Verbosity level 0-3. None in JSON → defaults to 1 at runtime.
    #[serde(rename = "Verbose")]
    pub verbose: Option<u8>,

    #[serde(rename = "Output")]
    pub output: Option<String>,

    /// Palace `OutputFormats` section — accepted for compatibility; REM ignores it.
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

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
pub enum ProblemType {
    Electrostatic,
    Magnetostatic,
    Eigenmode,
    Driven,
    Transient,
    /// Method of Moments (RWG + EFIE/CFIE) — REM extension, not in Palace
    MoM,
    /// Boundary Element Method (Laplace/Helmholtz) — REM extension, not in Palace
    BEM,
    /// Shooting and Bouncing Rays + Physical Optics — REM extension, not in Palace
    SBR,
    /// Hybrid Finite Element – Boundary Integral — REM extension, not in Palace
    #[serde(rename = "FEBI")]
    FEBI,
    /// Planar Method of Moments (uniform grid + FFT) — REM extension, not in Palace
    #[serde(rename = "Planar")]
    Planar,
}

// ---------------------------------------------------------------------------
// Model
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize)]
pub struct Model {
    #[serde(rename = "Mesh")]
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

    /// Palace `Postprocessing` under Domains — accepted for compatibility.
    #[serde(rename = "Postprocessing", default)]
    pub postprocessing: Option<DomainsPostprocessing>,

    /// Palace `CurrentDipole` — accepted for compatibility (not implemented).
    #[serde(rename = "CurrentDipole", default)]
    pub current_dipole: Vec<CurrentDipoleSpec>,
}

/// Palace `Domains.CurrentDipole` (Hertzian dipole source — not implemented).
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

/// Palace `Domains.Postprocessing` — accepted for Palace compatibility.
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

    /// Accepts scalar `9.3` or anisotropic array `[9.3, 9.3, 11.5]`; array → first element used.
    #[serde(rename = "Permittivity", default = "default_scalar_or_first_one",
            deserialize_with = "deserialize_scalar_or_first")]
    pub permittivity: f64,

    /// Accepts scalar `1.0` or anisotropic array `[1.0, 1.0, 1.0]`; array → first element used.
    #[serde(rename = "Permeability", default = "default_scalar_or_first_one",
            deserialize_with = "deserialize_scalar_or_first")]
    pub permeability: f64,

    /// Accepts scalar `3.0e-5` or anisotropic array `[3.0e-5, 3.0e-5, 8.6e-5]`; array → first element used.
    #[serde(rename = "LossTan", default,
            deserialize_with = "deserialize_scalar_or_first")]
    pub loss_tangent: f64,

    /// Magnetic loss tangent tan δ_m = μᵢ/μᵣ for lossy magnetic materials (ferrites, etc.).
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
    /// ε(ω) = ε∞ + Σ ωp² / (ω0² − ω² + jγω)
    /// Relevant only for driven (frequency-domain) solvers.
    #[serde(rename = "DrudeLorentz", default)]
    pub drude_lorentz: Vec<DrudeLorentzPole>,
}

/// One Drude-Lorentz oscillator pole.
///
/// Contributes `plasma_freq_sq / (resonance_freq_sq − ω² + j·damping·ω)` to εᵣ(ω).
///
/// For a Drude free-carrier term: set `ResonanceFreq = 0`, `PlasmaFreq = ωp/(2π)`,
/// `Damping = γ/(2π)`.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct DrudeLorentzPole {
    /// Plasma frequency fₚ [Hz] — the pole contribution strength.
    /// ωp² = (2π fₚ)²
    #[serde(rename = "PlasmaFreq", default)]
    pub plasma_freq: f64,

    /// Resonance frequency f₀ [Hz].  Zero → Drude (free-carrier) term.
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

    /// Palace `Periodic` / `FloquetWaveVector` boundaries — not implemented.
    #[serde(rename = "Periodic", default)]
    pub periodic: Vec<PeriodicSpec>,

    /// Palace boundary-level `Postprocessing` (SurfaceFlux, FarField, Dielectric).
    /// Accepted for Palace compatibility; REM logs warnings.
    #[serde(rename = "Postprocessing", default)]
    pub postprocessing_flux: Vec<BoundaryPostprocessingSpec>,
}

/// Palace `Boundaries.Periodic` — not implemented.
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

/// Palace `Boundaries.Postprocessing` — not yet implemented.
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

/// Resistive thin-sheet boundary condition (Ω/□ sheet resistance).
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
    /// Finite element polynomial order (1 = P1, 2 = P2, …)
    #[serde(rename = "Order", default = "default_order")]
    pub order: u8,

    /// Palace `Device` — REM is CPU-only; value is accepted and ignored.
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

    /// REM extension: MoM solver parameters (ignored by Palace).
    #[serde(rename = "MoM", default)]
    pub mom: Option<MomSolverConfig>,

    /// REM extension: SBR+ solver parameters (ignored by Palace).
    #[serde(rename = "SBR", default)]
    pub sbr: Option<SbrSolverConfig>,

    /// REM extension: Hybrid FE-BI solver parameters (ignored by Palace).
    #[serde(rename = "FEBI", default)]
    pub febi: Option<FeBiSolverConfig>,

    /// REM extension: Domain Decomposition Method solver parameters (ignored by Palace).
    #[serde(rename = "DDM", default)]
    pub ddm: Option<DdmSolverConfig>,

    /// REM extension: Planar MoM solver parameters (ignored by Palace).
    #[serde(rename = "Planar", default)]
    pub planar: Option<PlanarSolverConfig>,

    /// REM extension: near-to-far-field transform postprocessing.
    #[serde(rename = "FarField", default)]
    pub far_field: Option<FarFieldConfig>,
}

/// REM near-to-far-field configuration.
///
/// Computes radiation pattern from the driven solver's near-field solution.
/// Uses Kirchhoff approximation: far-field amplitude ∝ ∫ **E**(r') e^{jk r̂·r'} dS'
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
            device: "CPU".to_string(),
            eigenmode: None,
            driven: None,
            transient: None,
            electrostatic: None,
            magnetostatic: None,
            linear: LinearSolver::default(),
            mom: None,
            sbr: None,
            febi: None,
            ddm: None,
            planar: None,
            far_field: None,
        }
    }
}

fn default_device() -> String { "CPU".to_string() }

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

    #[serde(rename = "AdaptiveTol", default)]
    pub adaptive_tol: f64,

    /// Snapshot-based ROM order: number of full solves used to build the reduced basis.
    /// 0 (default) = disabled; 4–16 recommended for smooth S-parameter sweeps.
    /// When enabled, only `RomOrder` full complex solves are performed; all other
    /// frequency points are evaluated via the reduced system (much cheaper).
    #[serde(rename = "RomOrder", default)]
    pub rom_order: usize,

    /// Palace `Samples` — accepted, not implemented (use MinFreq/MaxFreq/FreqStep).
    #[serde(rename = "Samples", default)]
    pub samples: Vec<FreqSampleSpec>,

    /// Palace `Save` array — accepted, not implemented (use SaveStep integer).
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

    /// Palace `Excitation` waveform type — accepted, not fully implemented.
    #[serde(rename = "Excitation", default)]
    pub excitation: String,

    /// Palace `ExcitationFreq` [GHz] — accepted, not fully implemented.
    #[serde(rename = "ExcitationFreq", default)]
    pub excitation_freq: f64,

    /// Palace `ExcitationWidth` [ns] — accepted, not fully implemented.
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

    /// Palace `KSPType` — "CG"/"PCG" routes to PCG; "GMRES"/"" uses GMRES (default).
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

    /// Palace `ComplexCoarseSolve` — accepted, not implemented.
    #[serde(rename = "ComplexCoarseSolve", default)]
    pub complex_coarse_solve: bool,
}

impl LinearSolver {
    /// Returns `true` if `KSPType` is "CG" or "PCG" (case-insensitive).
    /// This hints that the caller should prefer the PCG path over GMRES.
    pub fn prefers_pcg(&self) -> bool {
        matches!(self.ksp_type.to_lowercase().as_str(), "cg" | "pcg")
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
// MoM solver config (REM extension — ignored by Palace)
// ---------------------------------------------------------------------------

/// MoM solver parameters, placed under `Solver.MoM` in the config file.
#[derive(Debug, Clone, Deserialize)]
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

    /// CFIE mixing coefficient α ∈ [0,1]: 0 = pure EFIE, 1 = pure MFIE
    #[serde(rename = "Alpha", default = "default_cfie_alpha")]
    pub alpha: f64,

    /// Convergence tolerance for singular integrals
    #[serde(rename = "SingularTol", default = "default_singular_tol")]
    pub singular_tol: f64,

    /// Linear solver for Z·I = V: "Direct" | "GMRES" | "ACA" | "FMM"
    #[serde(rename = "FastSolver", default = "default_fast_solver")]
    pub fast_solver: String,

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

    /// Near-field source file path. When set, the RHS is built from the
    /// near-field CSV data instead of the plane-wave model.  The file
    /// contains spatially sampled E/H fields exported from a previous
    /// simulation, enabling multi-solver near-field coupling.
    #[serde(rename = "NearFieldSource", default)]
    pub near_field_source: Option<String>,

    /// Snapshot ROM acceleration for S-parameter frequency sweeps.
    /// `0` disables ROM (default); positive value sets the number of anchor
    /// frequencies at which a full MoM solve is performed — all other
    /// frequencies use the Galerkin-projected low-dimensional system.
    /// Typical values: 4–16 for narrow-band, 8–32 for wideband.
    #[serde(rename = "RomOrder", default)]
    pub rom_order: usize,

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
}

fn default_mom_equation() -> String { "CFIE".to_string() }
fn default_mom_basis()     -> String { "RWG".to_string()  }
fn default_cfie_alpha()    -> f64    { 0.5 }
fn default_singular_tol()  -> f64    { 1.0e-6 }
fn default_fast_solver()   -> String { "Direct".to_string() }
fn default_polarization()  -> String { "theta".to_string() }
fn default_ref_impedance() -> f64    { 50.0 }
fn default_port_direction() -> String { "x".to_string() }
fn default_amr_theta()     -> f64    { 0.5 }
fn default_deembed_eps_eff() -> f64  { 1.0 }
fn default_mom_port_kind() -> String { "Lumped".to_string() }

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

    /// Reference-plane de-embedding length [m] for this port.
    /// Positive values shift the reference plane away from the port.
    #[serde(rename = "DeembedLength", default)]
    pub deembed_length: f64,
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
}

/// Single dielectric layer in the substrate stack.
#[derive(Debug, Clone, Deserialize)]
pub struct SubstrateLayerConfig {
    /// Relative permittivity (isotropic for now)
    #[serde(rename = "Permittivity", default = "default_permittivity")]
    pub permittivity: f64,

    /// Loss tangent: tan(δ) for dissipation model
    #[serde(rename = "LossTangent", default)]
    pub loss_tangent: f64,

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
    /// - `"Debye"`: relaxation model  ε(ω) = ε_∞ + (ε_s − ε_∞) / (1 + jωτ)
    /// - `"Lorentz"`: oscillator model ε(ω) = ε_∞ + Δε·ω₀² / (ω₀² − ω² + jγω)
    #[serde(rename = "Dispersion", default)]
    pub dispersion: Option<DispersionModel>,

    /// Anisotropic permittivity tensor [ε_xx, ε_yy, ε_zz].
    ///
    /// When set, overrides the scalar `Permittivity` field.  Supports
    /// uniaxial anisotropy (ε_xx = ε_yy ≠ ε_zz) typical of laminated PCB
    /// substrates (e.g., Rogers 4350B: ε_xy ≈ 3.48, ε_z ≈ 3.66).
    /// The Sommerfeld integral uses ε_t = ε_xx for lateral wavenumber and
    /// ε_z = ε_zz for vertical wavenumber in TM modes.
    #[serde(rename = "PermittivityTensor", default)]
    pub permittivity_tensor: Option<[f64; 3]>,
}

/// Frequency-dependent permittivity dispersion model for a substrate layer.
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "Type")]
pub enum DispersionModel {
    /// Debye relaxation:  ε(ω) = ε_∞ + (ε_s − ε_∞) / (1 + jωτ)
    Debye {
        /// High-frequency permittivity ε_∞ (lossless plateau above resonance)
        #[serde(rename = "EpsInf")]
        eps_inf: f64,
        /// Static (DC) permittivity ε_s
        #[serde(rename = "EpsStatic")]
        eps_static: f64,
        /// Relaxation time constant τ [s]
        #[serde(rename = "RelaxationTime")]
        tau_s: f64,
    },
    /// Lorentz oscillator: ε(ω) = ε_∞ + Δε·ω₀² / (ω₀² − ω² + jγω)
    Lorentz {
        /// High-frequency permittivity ε_∞
        #[serde(rename = "EpsInf")]
        eps_inf: f64,
        /// Permittivity increment Δε = ε_s − ε_∞
        #[serde(rename = "DeltaEps")]
        delta_eps: f64,
        /// Resonant angular frequency ω₀ [rad/s]
        #[serde(rename = "OmegaRes")]
        omega0_rad_per_s: f64,
        /// Damping coefficient γ [rad/s]
        #[serde(rename = "Gamma")]
        gamma_rad_per_s: f64,
    },
}

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
        }
    }
}

impl SubstrateLayerConfig {
    /// Complex relative permittivity ε_r(f) at frequency `freq_hz`.
    ///
    /// Priority:
    /// 1. Dispersion model (`Dispersion` key), if present.
    /// 2. Anisotropic tensor (`PermittivityTensor`): uses ε_xx as lateral ε.
    /// 3. Scalar `Permittivity` + `LossTangent`.
    pub fn eps_r_complex(&self, freq_hz: f64) -> num_complex::Complex64 {
        use num_complex::Complex64;
        use std::f64::consts::PI;
        if let Some(ref model) = self.dispersion {
            let omega = 2.0 * PI * freq_hz;
            model.eps_r_at_omega(omega)
        } else if let Some([eps_xx, _, _]) = self.permittivity_tensor {
            // Use lateral (xx) component; imaginary part from loss tangent
            Complex64::new(eps_xx, -eps_xx * self.loss_tangent)
        } else {
            let eps = self.permittivity;
            Complex64::new(eps, -eps * self.loss_tangent)
        }
    }

    /// Vertical (z-axis) complex permittivity ε_zz(f).
    ///
    /// For isotropic and non-tensor layers, equals `eps_r_complex`.
    /// For anisotropic tensor layers, returns ε_zz component.
    pub fn eps_r_z_complex(&self, freq_hz: f64) -> Option<num_complex::Complex64> {
        use num_complex::Complex64;
        if let Some([_, _, eps_zz]) = self.permittivity_tensor {
            Some(Complex64::new(eps_zz, -eps_zz * self.loss_tangent))
        } else {
            None // isotropic: caller may treat as same as eps_r_complex
        }
    }
}

fn default_bottom_pec() -> bool { true }
fn default_permittivity() -> f64 { 1.0 }
fn default_permeability() -> f64 { 1.0 }

// ---------------------------------------------------------------------------
// SBR+ solver config (REM extension — ignored by Palace)
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
// FE-BI solver config (REM extension — ignored by Palace)
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

    /// CFIE mixing coefficient α ∈ [0,1]: 0 = pure EFIE, 1 = pure MFIE
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
// DDM solver config (REM extension — ignored by Palace)
// ---------------------------------------------------------------------------

/// Domain Decomposition Method solver parameters, placed under `Solver.DDM`.
///
/// Implements Robin-condition Schwarz iteration for parallel FEM solves.
#[derive(Debug, Clone, Deserialize)]
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
}

fn default_ddm_num_subdomains() -> usize  { 4 }
fn default_ddm_method()         -> String { "Schwarz".to_string() }
fn default_ddm_robin_order()    -> u8     { 1 }
fn default_ddm_tolerance()      -> f64    { 1.0e-6 }
fn default_ddm_max_iter()       -> usize  { 100 }
fn default_ddm_partition_type() -> String { "Dual".to_string() }

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
// Postprocessing (REM extension — ignored by Palace)
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

    // --- Solver.Linear ---
    let l = &cfg.solver.linear;
    // CG / PCG / GMRES are all accepted: SPD solvers use PCG by default;
    // Driven uses GMRES (or PCG if Solver.Linear.KSPType is "CG"/"PCG").
    let ksp_lower = l.ksp_type.to_lowercase();
    if !l.ksp_type.is_empty()
        && ksp_lower != "gmres"
        && ksp_lower != "cg"
        && ksp_lower != "pcg"
        && ksp_lower != "default"
    {
        warn_unsupported(
            &format!("Solver.Linear.KSPType = \"{}\"", l.ksp_type),
            "Only GMRES and CG/PCG are supported; value is ignored",
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

    // --- Solver.MoM ---
    if let Some(ref m) = cfg.solver.mom {
        if let Some(ref sub) = m.substrate {
            log::info!(
                "[REM] Solver.MoM.Substrate: {} layer(s), bottom_pec={} -- \
                 layered-media Green function (DCIM) will be used for assembly.",
                sub.layers.len(), sub.bottom_pec
            );
        }
        if !m.ports.is_empty() {
            log::info!(
                "[REM] Solver.MoM.Ports: {} port(s), ref Z0={} Ohm -- \
                 S-parameter sweep mode active; Touchstone and port-S.csv will be written.",
                m.ports.len(), m.ref_impedance
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
                "[REM] Domains.Postprocessing.Probe: {} probe(s) — \
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
                "[REM] Boundaries.Periodic: {} spec(s), {} pair(s) — non-zero FloquetWaveVector \
                 detected. Complex phase-shift BCs are not yet supported; Γ-point (k=0) pairs will \
                 be applied, others skipped.",
                cfg.boundaries.periodic.len(), n_pairs
            );
        } else {
            log::info!(
                "[REM] Boundaries.Periodic: {} spec(s), {} pair(s) — Γ-point periodic BCs enabled.",
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
                "Boundaries.LumpedPort[{}]: {} elements — all element attributes mapped to port BC.",
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
// Postprocessing (REM extension — ignored by Palace)
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

    #[test]
    fn debye_at_dc_equals_eps_static() {
        let model = DispersionModel::Debye { eps_inf: 2.0, eps_static: 4.5, tau_s: 1e-11 };
        // At ω≈0, ε(ω) → ε_s
        let eps = model.eps_r_at_omega(1.0); // very low freq
        assert!((eps.re - 4.5).abs() < 0.01, "Debye DC: re={:.4}", eps.re);
        assert!(eps.im < 0.0, "Debye DC: imaginary part should be negative (loss)");
    }

    #[test]
    fn debye_at_high_freq_equals_eps_inf() {
        let model = DispersionModel::Debye { eps_inf: 2.0, eps_static: 4.5, tau_s: 1e-11 };
        // At very high ω (ωτ >> 1), ε(ω) → ε_∞
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
        // At resonance: ε = ε_∞ + Δε·ω₀²/(jγω₀) = ε_∞ - j Δε·ω₀/γ
        assert!(eps.im.abs() > eps.re.abs(), "Lorentz at resonance: should be loss-dominated");
    }

    #[test]
    fn substrate_layer_eps_r_complex_uses_dispersion_when_set() {
        use std::f64::consts::PI;
        let layer = SubstrateLayerConfig {
            permittivity: 4.0,
            loss_tangent: 0.02,
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
        // Debye at 1 GHz: ωτ = 2π·1e9·1e-11 ≈ 0.063; still close to ε_s
        assert!(eps.re > 4.0 && eps.re < 6.5, "Debye 1 GHz re={:.3}", eps.re);
        assert!(eps.im < 0.0, "Should have negative imaginary part (loss)");
    }

    #[test]
    fn substrate_layer_eps_r_z_returns_tensor_zz() {
        let layer = SubstrateLayerConfig {
            permittivity: 3.48,
            loss_tangent: 0.004,
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
            permeability: 1.0,
            thickness: 1e-3,
            name: String::new(),
            dispersion: None,
            permittivity_tensor: None,
        };
        assert!(layer.eps_r_z_complex(1e9).is_none());
    }
}
