use crate::preprocess::expand_ranges;
use serde::{Deserialize, Deserializer};

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
}

#[derive(Debug, Clone, Deserialize)]
pub struct MaterialSpec {
    /// Physical group attribute IDs. Accepts either `[1,2,3]` or `"1,3-5"`.
    #[serde(rename = "Attributes", deserialize_with = "deserialize_attributes")]
    pub attributes: Vec<u32>,

    #[serde(rename = "Permittivity", default = "default_one")]
    pub permittivity: f64,

    #[serde(rename = "Permeability", default = "default_one")]
    pub permeability: f64,

    #[serde(rename = "LossTan", default)]
    pub loss_tangent: f64,

    #[serde(rename = "Conductivity", default)]
    pub conductivity: f64,
}

fn default_one() -> f64 { 1.0 }

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
}

fn default_order() -> u8 { 1 }

impl Default for SolverConfig {
    fn default() -> Self {
        SolverConfig {
            order: 1,
            eigenmode: None,
            driven: None,
            transient: None,
            electrostatic: None,
            magnetostatic: None,
            linear: LinearSolver::default(),
        }
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
}

fn default_save_step() -> usize { 1 }

#[derive(Debug, Clone, Deserialize)]
pub struct TransientSolver {
    #[serde(rename = "Type", default = "default_transient_type")]
    pub solver_type: String,

    #[serde(rename = "MaxTime")]
    pub max_time: f64,

    #[serde(rename = "TimeStep")]
    pub time_step: f64,

    #[serde(rename = "SaveStep", default = "default_save_step")]
    pub save_step: usize,
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

    #[serde(rename = "Tol", default = "default_linear_tol")]
    pub tol: f64,

    #[serde(rename = "MaxIter", default = "default_linear_maxiter")]
    pub max_iter: usize,

    #[serde(rename = "MGLevels", default = "default_mg_levels")]
    pub mg_levels: usize,

    #[serde(rename = "PCType", default = "default_pc_type")]
    pub pc_type: String,
}

impl Default for LinearSolver {
    fn default() -> Self {
        LinearSolver {
            solver_type: default_linear_type(),
            tol: default_linear_tol(),
            max_iter: default_linear_maxiter(),
            mg_levels: default_mg_levels(),
            pc_type: default_pc_type(),
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
