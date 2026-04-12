//! Output formatting utilities for Palace-compatible console messages.
//!
//! Provides banner, system info reporting, and solver output formatting
//! following Palace's professional output style.

use std::num::NonZeroUsize;

/// Print the REM banner.
pub fn print_banner() {
    let banner = r#"
  ________  __________   _____  __
 /  _____/_/ ____     \_/     \|  |
|  |   /|  |  |  |  |  \|  Y   \  |
|  |___\|  |  |  |  |  /|  |  /|  |__
 \_____ /|__|  |____/  /|__|  \|____/
      \_/             \/
"#;
    log::info!("{}", banner);
}

/// Print system configuration information.
///
/// # Arguments
/// * `num_threads` - Number of OpenMP threads available
/// * `version` - Optional version string (e.g., git commit hash)
pub fn print_system_info(num_threads: Option<usize>, version: Option<&str>) {
    if let Some(v) = version {
        log::info!("REM version: {}", v);
    }

    let threads = num_threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(NonZeroUsize::get)
            .unwrap_or(1)
    });

    if threads > 1 {
        log::info!("Running with {} threads", threads);
    }
    log::info!("");
}

/// Print a solver section header.
///
/// # Arguments
/// * `solver_name` - Name of the solver (e.g., "Eigenmode", "Driven")
pub fn print_solver_header(solver_name: &str) {
    log::info!("=== {} solver ===", solver_name);
}

/// Print a named section with indentation.
///
/// # Arguments
/// * `section_name` - Name of the section
/// * `content` - Formatted section content (will be indented)
pub fn print_section(section_name: &str, content: &str) {
    log::info!("\n{}:", section_name);
    for line in content.lines() {
        if !line.is_empty() {
            log::info!("  {}", line);
        } else {
            log::info!("");
        }
    }
}

/// Print a parameter with scientific notation if large/small.
///
/// # Arguments
/// * `name` - Parameter name
/// * `value` - Parameter value
/// * `unit` - Optional unit (e.g., "Hz", "m")
pub fn print_parameter(name: &str, value: f64, unit: Option<&str>) {
    if value.abs() < 1e-4 && value != 0.0 || value.abs() >= 1e6 {
        match unit {
            Some(u) => log::info!("  {} = {:.3e} {}", name, value, u),
            None => log::info!("  {} = {:.3e}", name, value),
        }
    } else {
        match unit {
            Some(u) => log::info!("  {} = {} {}", name, value, u),
            None => log::info!("  {} = {}", name, value),
        }
    }
}

/// Print a complete configuration section for a solver.
///
/// Used to report solver-specific parameters before the solve begins.
pub fn print_solver_config(solver_name: &str, params: &[(&str, String)]) {
    log::info!("\nConfiguring {} solver:", solver_name);
    for (name, value) in params {
        log::info!("  {} = {}", name, value);
    }
    log::info!("");
}

/// Print progress/status message with timestamp-like format.
pub fn print_progress(message: &str) {
    log::info!(">> {}", message);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_print_parameter_scientific() {
        // Just verify functions compile and don't panic
        print_parameter("small_val", 1e-6, Some("m"));
        print_parameter("large_val", 1e9, Some("Hz"));
    }

    #[test]
    fn test_print_section() {
        let content = "Line 1\nLine 2\nLine 3";
        print_section("Test Section", content);
    }
}
