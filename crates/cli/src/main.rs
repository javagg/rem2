use anyhow::Context;
use clap::Parser;
use rem_config::{load_config, ProblemType};
use rem_parallel::{NoComm, WorldComm, Comm};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "rem", about = "Rust Electromagnetic Solver — Palace-compatible")]
struct Args {
    /// Palace-format config file (.json or .yaml)
    config: PathBuf,

    /// Override output directory
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Increase log verbosity (-v info, -vv debug, -vvv trace)
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
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

    let mut config = load_config(&args.config)
        .with_context(|| format!("reading config: {}", args.config.display()))?;

    if let Some(out) = args.output {
        config.problem.output = Some(out.to_string_lossy().into_owned());
    }

    let comm = if cfg!(target_arch = "wasm32") {
        Box::new(WorldComm::new()) as Box<dyn Comm>
    } else {
        Box::new(NoComm) as Box<dyn Comm>
    };

    match config.problem.problem_type {
        ProblemType::Electrostatic => {
            rem_electrostatic::run(&config, comm.as_ref() as &dyn Comm)?;
        }
        ProblemType::Magnetostatic => {
            rem_magnetostatic::run(&config, comm.as_ref() as &dyn Comm)?;
        }
        ProblemType::Eigenmode => {
            anyhow::bail!("Eigenmode solver not yet implemented (v0.2)");
        }
        ProblemType::Driven => {
            anyhow::bail!("Driven solver not yet implemented (v0.2)");
        }
        ProblemType::Transient => {
            anyhow::bail!("Transient solver not yet implemented (v1.0)");
        }
    }

    Ok(())
}
