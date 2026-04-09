//! Peak memory usage reporting.
//!
//! Platform support:
//!   - Linux:   reads `VmPeak` from `/proc/self/status`
//!   - Windows: uses `GetProcessMemoryInfo` via raw syscall
//!   - WASM:    returns None (JS heap size not reliably available without web_sys)
//!   - Other:   returns None

/// Return the current peak virtual memory usage in bytes, if determinable.
///
/// - Linux:   `VmPeak` from `/proc/self/status` (peak virtual memory since process start)
/// - Windows: `PeakWorkingSetSize` from `GetProcessMemoryInfo`
/// - Other:   returns `None`
pub fn peak_memory_bytes() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        linux_peak_memory()
    }
    #[cfg(target_os = "windows")]
    {
        windows_peak_memory()
    }
    #[cfg(not(any(target_os = "linux", target_os = "windows")))]
    {
        None
    }
}

/// Log peak memory usage at INFO level.
///
/// Call this at the end of a solver run to report memory consumption.
pub fn report_peak_memory(context: &str) {
    match peak_memory_bytes() {
        Some(bytes) => {
            let mib = bytes as f64 / (1024.0 * 1024.0);
            if mib >= 1024.0 {
                log::info!("[REM] {} — peak memory: {:.2} GiB", context, mib / 1024.0);
            } else {
                log::info!("[REM] {} — peak memory: {:.1} MiB", context, mib);
            }
        }
        None => {
            log::debug!("[REM] {} — peak memory: unavailable on this platform", context);
        }
    }
}

// ---------------------------------------------------------------------------
// Linux implementation
// ---------------------------------------------------------------------------

#[cfg(target_os = "linux")]
fn linux_peak_memory() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        // VmPeak: <kB>
        if let Some(rest) = line.strip_prefix("VmPeak:") {
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Windows implementation
// ---------------------------------------------------------------------------

#[cfg(target_os = "windows")]
fn windows_peak_memory() -> Option<u64> {
    // Use GetProcessMemoryInfo from psapi.dll via raw extern.
    // No winapi crate dependency needed.
    #[repr(C)]
    #[allow(non_snake_case)]
    struct ProcessMemoryCounters {
        cb: u32,
        PageFaultCount: u32,
        PeakWorkingSetSize: usize,
        WorkingSetSize: usize,
        QuotaPeakPagedPoolUsage: usize,
        QuotaPagedPoolUsage: usize,
        QuotaPeakNonPagedPoolUsage: usize,
        QuotaNonPagedPoolUsage: usize,
        PagefileUsage: usize,
        PeakPagefileUsage: usize,
    }

    #[link(name = "psapi")]
    extern "system" {
        fn GetCurrentProcess() -> *mut std::ffi::c_void;
        fn GetProcessMemoryInfo(
            Process: *mut std::ffi::c_void,
            ppsmemCounters: *mut ProcessMemoryCounters,
            cb: u32,
        ) -> i32;
    }

    let mut pmc = ProcessMemoryCounters {
        cb: std::mem::size_of::<ProcessMemoryCounters>() as u32,
        PageFaultCount: 0,
        PeakWorkingSetSize: 0,
        WorkingSetSize: 0,
        QuotaPeakPagedPoolUsage: 0,
        QuotaPagedPoolUsage: 0,
        QuotaPeakNonPagedPoolUsage: 0,
        QuotaNonPagedPoolUsage: 0,
        PagefileUsage: 0,
        PeakPagefileUsage: 0,
    };

    unsafe {
        let proc = GetCurrentProcess();
        let ok = GetProcessMemoryInfo(proc, &mut pmc, pmc.cb);
        if ok != 0 {
            Some(pmc.PeakWorkingSetSize as u64)
        } else {
            None
        }
    }
}
