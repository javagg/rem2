//! Terminal UI primitives shared across all solver crates.
//!
//! Provides solver-agnostic output: logo, progress, summary.
//! Solver-specific UI (box model, layer stack) lives in each solver's own
//! `terminal_ui` extension module.

/// REM Pro logo banner.
pub fn ui_banner() {
    eprintln!("");
    eprintln!("╔══════════════════════════════════════════════════════════════════╗");
    eprintln!("║          Rust Electromagnetic Solver (REM)                     ║");
    eprintln!("╚══════════════════════════════════════════════════════════════════╝");
    eprintln!("");
}

/// Start-of-solving separator.
pub fn ui_solving_header() {
    eprintln!("──────────────────────────────────────────────────────────────────");
}

/// Adaptive Band Synthesis header.
pub fn ui_abs_header(target: usize, poles: usize) {
    eprintln!("  Adaptive Band Synthesis (ABS)  |  target={} pts  |  poles={}", target, poles);
}

/// Per-frequency progress line.
pub fn ui_frequency_progress(i: usize, n: usize, freq: f64, elapsed: std::time::Duration) {
    let pad = n.to_string().len();
    eprintln!("  [{:>pad$}/{}]  f = {:.3e} Hz  ({:.2?})", i, n, freq, elapsed, pad = pad);
}

/// Solve complete summary.
pub fn ui_summary(n_solved: usize, n_total: usize, snp_path: &std::path::Path, elapsed: std::time::Duration) {
    eprintln!("┌────────────────────────────────────────────────────────────────┐");
    eprintln!("│  Solve complete                                              │");
    eprintln!("│    Frequencies:  {} / {}                                      │", n_solved, n_total);
    eprintln!("│    Output:       {}                                                │", snp_path.display());
    eprintln!("│    Elapsed:      {:.2?}                                                     │", elapsed);
    eprintln!("└────────────────────────────────────────────────────────────────┘");
}
