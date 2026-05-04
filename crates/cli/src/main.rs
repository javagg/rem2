use anyhow::Context;
use clap::{Parser, ValueEnum};
use rem_config::{load_config, ProblemType};
use rem_convert::{convert_project_to_rem, ProjectFormat, Sonnet19Overrides};
use rem_parallel::{NoComm, WorldComm, Comm};
use std::path::PathBuf;

mod output;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum CliProjectFormat {
    Sonnet19,
    Ansys,
    Ads,
}

impl From<CliProjectFormat> for ProjectFormat {
    fn from(value: CliProjectFormat) -> Self {
        match value {
            CliProjectFormat::Sonnet19 => ProjectFormat::Sonnet19,
            CliProjectFormat::Ansys => ProjectFormat::Ansys,
            CliProjectFormat::Ads => ProjectFormat::Ads,
        }
    }
}

#[derive(Parser)]
#[command(name = "rem", about = "Rust Electromagnetic Solver — Palace-compatible")]
struct Args {
    /// Palace-format config file (.json or .yaml)
    config: Option<PathBuf>,

    /// Generic project path for conversion (Sonnet19/Ansys/ADS).
    #[arg(long)]
    project: Option<PathBuf>,

    /// Project format for --project.
    #[arg(long, value_enum)]
    format: Option<CliProjectFormat>,

    /// Sonnet 19 XML project file (.xml) to convert into REM config + .msh
    #[arg(long)]
    sonnet19_xml: Option<PathBuf>,

    /// Output REM config path (JSON). Used with --sonnet19-xml.
    #[arg(long)]
    out_config: Option<PathBuf>,

    /// Output GMSH mesh path (.msh). Used with --sonnet19-xml.
    #[arg(long)]
    out_msh: Option<PathBuf>,

    /// Override frequency start [Hz] in generated MoM config.
    #[arg(long)]
    freq_min: Option<f64>,

    /// Override frequency end [Hz] in generated MoM config.
    #[arg(long)]
    freq_max: Option<f64>,

    /// Override frequency step [Hz] in generated MoM config.
    #[arg(long)]
    freq_step: Option<f64>,

    /// Conversion debug: export pre-meshing geometry STEP into the conversion output directory (Sonnet19 only)
    #[arg(long, default_value_t = false)]
    output_step: bool,

    /// Solver output directory, or conversion output directory in conversion mode
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Increase log verbosity (-v info, -vv debug, -vvv trace)
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    /// JSON API mode: read config from stdin, output results JSON to stdout.
    /// All logging goes to stderr.
    #[arg(long, default_value_t = false)]
    json: bool,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let log_level = match args.verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or(log_level),
    )
    .init();

    // Print banner and system info only in info mode or verbose
    if args.verbose >= 1 {
        output::print_banner();
        output::print_system_info(None, None);
    }

    if let Some(project_path) = args.project.as_ref().or(args.sonnet19_xml.as_ref()) {
        let format = match (args.format, args.sonnet19_xml.as_ref()) {
            (Some(fmt), _) => fmt.into(),
            (None, Some(_)) => ProjectFormat::Sonnet19,
            (None, None) => {
                anyhow::bail!("missing --format with --project; expected one of: sonnet19, ansys, ads")
            }
        };
        let default_out_config = project_path.with_extension("json");
        let default_out_msh = project_path.with_extension("msh");

        let out_config = if let Some(path) = args.out_config.clone() {
            path
        } else if let Some(out_dir) = args.output.clone() {
            let name = default_out_config
                .file_name()
                .map(|n| n.to_owned())
                .unwrap_or_else(|| "converted.json".into());
            out_dir.join(name)
        } else {
            default_out_config
        };

        let out_msh = if let Some(path) = args.out_msh.clone() {
            path
        } else if let Some(out_dir) = args.output.clone() {
            let name = default_out_msh
                .file_name()
                .map(|n| n.to_owned())
                .unwrap_or_else(|| "converted.msh".into());
            out_dir.join(name)
        } else {
            default_out_msh
        };

        let debug_step = if matches!(format, ProjectFormat::Sonnet19) && args.output_step {
            let out_dir = if let Some(dir) = args.output.clone() {
                dir
            } else if let Some(parent) = out_config.parent() {
                parent.to_path_buf()
            } else {
                std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."))
            };
            let stem = out_config
                .file_stem()
                .and_then(|s| s.to_str())
                .filter(|s| !s.is_empty())
                .unwrap_or("converted");
            Some(out_dir.join(format!("{}_geometry.step", stem)))
        } else {
            None
        };

        convert_project_to_rem(
            format,
            project_path,
            &out_config,
            &out_msh,
            Sonnet19Overrides {
                freq_min: args.freq_min,
                freq_max: args.freq_max,
                freq_step: args.freq_step,
                debug_step,
            },
        )?;
        return Ok(());
    }

    // JSON API mode: read config from stdin, output results to stdout
    if args.json {
        let mut input = String::new();
        std::io::Read::read_to_string(&mut std::io::stdin(), &mut input)
            .context("reading config from stdin")?;

        let mut config = rem_config::load_config_from_str(&input, rem_config::ConfigFormat::Json)
            .context("parsing JSON config from stdin")?;

        if let Some(out) = args.output {
            config.problem.output = Some(out.to_string_lossy().into_owned());
        }

        let comm: Box<dyn Comm> = if cfg!(target_arch = "wasm32") {
            Box::new(WorldComm::new())
        } else {
            Box::new(NoComm)
        };

        let start = std::time::Instant::now();
        let result = run_solver(&config, comm.as_ref());

        let elapsed_s = start.elapsed().as_secs_f64();
        let output = serde_json::json!({
            "status": if result.is_ok() { "ok" } else { "error" },
            "problem_type": format!("{:?}", config.problem.problem_type),
            "elapsed_s": elapsed_s,
            "output_dir": config.problem.output_dir(),
            "error": result.as_ref().err().map(|e| e.to_string()),
        });

        // Print JSON result to stdout (logs go to stderr)
        println!("{}", serde_json::to_string_pretty(&output).unwrap());
        return result;
    }

    let config_path = args
        .config
        .as_ref()
        .context("missing config path; use rem <config.json|yaml> or conversion mode: --project <file> --format <sonnet19|ansys|ads>")?;

    let mut config = load_config(config_path)
        .with_context(|| format!("reading config: {}", config_path.display()))?;

    if let Some(out) = args.output {
        config.problem.output = Some(out.to_string_lossy().into_owned());
    }

    let comm = if cfg!(target_arch = "wasm32") {
        Box::new(WorldComm::new()) as Box<dyn Comm>
    } else {
        Box::new(NoComm) as Box<dyn Comm>
    };

    run_solver(&config, comm.as_ref())
}

fn run_solver(config: &rem_config::PalaceConfig, comm: &dyn Comm) -> anyhow::Result<()> {
    match config.problem.problem_type {
        ProblemType::Electrostatic => {
            rem_electrostatic::run(config, comm)?;
        }
        ProblemType::Magnetostatic => {
            rem_magnetostatic::run(config, comm)?;
        }
        ProblemType::Eigenmode => {
            rem_eigenmode::run(config, comm)?;
        }
        ProblemType::Driven => {
            rem_driven::run(config, comm)?;
        }
        ProblemType::Transient => {
            rem_transient::run(config, comm)?;
        }
        ProblemType::MoM => {
            rem_mom::run(config)?;
        }
        ProblemType::BEM => {
            rem_bem::run(config)?;
        }
        ProblemType::Planar => {
            rem_planar::run(config)?;
        }
        ProblemType::SBR => {
            rem_sbr::run(config)?;
        }
        ProblemType::FEBI => {
            rem_febi::run(config)?;
        }
    }

    Ok(())
}
