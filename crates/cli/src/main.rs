use anyhow::Context;
use clap::Parser;
use rem_config::{load_config, ProblemType};
use rem_parallel::{NoComm, WorldComm, Comm};
use std::path::PathBuf;

#[allow(dead_code)]
mod output;

#[derive(Parser)]
#[command(name = "rem", about = "Rust Electromagnetic Solver — Palace-compatible")]
struct Args {
    /// Palace-format config file (.json or .yaml)
    config: Option<PathBuf>,

    /// Solver output directory
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
        .context("missing config path; use rem <config.json|yaml>")?;

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
        ProblemType::BEM => {
            rem_bem::run(config)?;
        }
        ProblemType::MoM
        | ProblemType::Planar
        | ProblemType::SBR
        | ProblemType::FEBI
        | ProblemType::Parametric => {
            anyhow::bail!(
                "{:?} is available only in the private rem-pro workspace; use the pro CLI for this problem type",
                config.problem.problem_type
            );
        }
    }

    Ok(())
}
