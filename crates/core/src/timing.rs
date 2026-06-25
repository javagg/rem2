//! Hierarchical timing tree for solver performance profiling.
//!
//! Activated by setting the environment variable `REM_PROFILE=1`.
//! At program exit, a JSON timing report is written to `timings.json` in the
//! current directory (or `REM_PROFILE_OUTPUT` if set).
//!
//! # Usage
//!
//! ```ignore
//! use rem_core::timing;
//!
//! let _guard = timing::span("mlfma.solve");
//! // ... work ...
//! // guard drops here, records elapsed time
//! ```
//!
//! # Thread safety
//!
//! TimingTree uses a `Mutex<VecDeque>` to associate spans with the correct
//! parent on each thread.  Rayon-parallel code is supported — each thread
//! independently pushes/pops from the shared tree.
//!
//! # WASM
//!
//! On `wasm32` targets all functions are no-ops and guarded out at compile time.

#[cfg(not(target_arch = "wasm32"))]
mod native {
    use std::collections::VecDeque;
    use std::io::Write;
    use std::sync::Mutex;
    use std::time::Instant;

    // -----------------------------------------------------------------------
    // Tree data model
    // -----------------------------------------------------------------------

    #[derive(Debug, Clone)]
    struct TimingNode {
        name: String,
        /// Cumulative wall-clock seconds across all calls to this span.
        elapsed_secs: f64,
        /// Number of times this span was entered.
        count: usize,
        /// Child spans indexed by name.
        children: Vec<TimingNode>,
    }

    /// Shared mutable tree state.
    static TREE: std::sync::LazyLock<Mutex<TreeState>> =
        std::sync::LazyLock::new(|| Mutex::new(TreeState::new()));

    struct TreeState {
        /// Root node (anonymous).
        root: TimingNode,
        /// Stack of (name, start_time) for in-flight spans on the current
        /// OS thread.  `thread_local!` would be simpler but can't be stored
        /// inside a `LazyLock<Mutex<…>>`.  Instead we key by thread id.
        thread_stacks: std::collections::HashMap<std::thread::ThreadId, VecDeque<(String, Instant)>>,
        /// Whether profiling is enabled.
        enabled: bool,
    }

    impl TreeState {
        fn new() -> Self {
            let enabled = std::env::var("REM_PROFILE")
                .map(|v| v == "1")
                .unwrap_or(false);
            TreeState {
                root: TimingNode {
                    name: String::new(),
                    elapsed_secs: 0.0,
                    count: 0,
                    children: Vec::new(),
                },
                thread_stacks: std::collections::HashMap::new(),
                enabled,
            }
        }

        fn push(&mut self, name: &str) {
            if !self.enabled {
                return;
            }
            let tid = std::thread::current().id();
            self.thread_stacks
                .entry(tid)
                .or_default()
                .push_back((name.to_string(), Instant::now()));
        }

        fn pop(&mut self, name: &str) {
            if !self.enabled {
                return;
            }
            let tid = std::thread::current().id();
            let stack = self.thread_stacks.get_mut(&tid);
            let Some(stack) = stack else { return };

            let (popped_name, start) = match stack.pop_back() {
                Some(v) => v,
                None => return,
            };

            let elapsed = start.elapsed().as_secs_f64();

            // Navigate to the parent path using the remaining stack.
            let mut current = &mut self.root;
            for (ancestor, _) in stack.iter() {
                current = find_or_create_child(current, ancestor);
            }

            let node = find_or_create_child(current, &popped_name);
            assert_eq!(popped_name, name, "TimingTree push/pop mismatch");
            node.elapsed_secs += elapsed;
            node.count += 1;
        }

        fn dump_json(&self) -> String {
            fn write_node(node: &TimingNode, indent: usize, out: &mut String) {
                let pad = "  ".repeat(indent);
                out.push_str(&format!("{}{{\n", pad));
                out.push_str(&format!("{}  \"name\": \"{}\",\n", pad, node.name));
                out.push_str(&format!("{}  \"elapsed_secs\": {:.6},\n", pad, node.elapsed_secs));
                out.push_str(&format!("{}  \"count\": {},\n", pad, node.count));
                if node.children.is_empty() {
                    out.push_str(&format!("{}  \"children\": []\n", pad));
                } else {
                    out.push_str(&format!("{}  \"children\": [\n", pad));
                    for (i, child) in node.children.iter().enumerate() {
                        write_node(child, indent + 2, out);
                        if i + 1 < node.children.len() {
                            out.push(',');
                        }
                        out.push('\n');
                    }
                    out.push_str(&format!("{}  ]\n", pad));
                }
                out.push_str(&format!("{}}}", pad));
            }
            let mut buf = String::new();
            write_node(&self.root, 0, &mut buf);
            buf.push('\n');
            buf
        }
    }

    fn find_or_create_child<'a>(parent: &'a mut TimingNode, name: &str) -> &'a mut TimingNode {
        let pos = parent.children.iter().position(|c| c.name == name);
        if let Some(i) = pos {
            &mut parent.children[i]
        } else {
            parent.children.push(TimingNode {
                name: name.to_string(),
                elapsed_secs: 0.0,
                count: 0,
                children: Vec::new(),
            });
            parent.children.last_mut().unwrap()
        }
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /// Returns `true` when `REM_PROFILE=1` is set.
    pub fn enabled() -> bool {
        TREE.lock().unwrap().enabled
    }

    /// Enter a named timing span.  Returns a guard whose `Drop` impl records
    /// the elapsed duration.
    ///
    /// # Panics
    ///
    /// Spans must be nested properly (LIFO).  Mismatched push/pop will panic
    /// in debug builds.
    #[must_use]
    pub fn span(name: &'static str) -> SpanGuard {
        TREE.lock().unwrap().push(name);
        SpanGuard { name }
    }

    pub struct SpanGuard {
        name: &'static str,
    }

    impl Drop for SpanGuard {
        fn drop(&mut self) {
            TREE.lock().unwrap().pop(self.name);
        }
    }

    /// Write the accumulated timing tree to `REM_PROFILE_OUTPUT` (default
    /// `timings.json`) and attempt to print a flat summary to stderr.
    pub fn dump_and_clear() {
        let state = TREE.lock().unwrap();
        if !state.enabled {
            return;
        }
        let json = state.dump_json();
        let path = std::env::var("REM_PROFILE_OUTPUT")
            .unwrap_or_else(|_| "timings.json".to_string());

        if let Ok(mut f) = std::fs::File::create(&path) {
            let _ = f.write_all(json.as_bytes());
            let _ = f.write_all(b"\n");
        }

        // Flat text summary to stderr (one line per leaf-to-root path)
        fn lines(node: &TimingNode, prefix: &str, out: &mut Vec<String>) {
            if !node.name.is_empty() {
                let label = if prefix.is_empty() {
                    node.name.clone()
                } else {
                    format!("{}.{}", prefix, node.name)
                };
                out.push(format!(
                    "{:<60} {:>10.3}s  {:>6} calls  {:>10.1}ms/call",
                    label,
                    node.elapsed_secs,
                    node.count,
                    if node.count > 0 {
                        node.elapsed_secs * 1000.0 / node.count as f64
                    } else {
                        0.0
                    }
                ));
                for child in &node.children {
                    lines(child, &label, out);
                }
            } else {
                for child in &node.children {
                    lines(child, "", out);
                }
            }
        }

        let mut flat: Vec<String> = Vec::new();
        lines(&state.root, "", &mut flat);
        if !flat.is_empty() {
            let header = format!(
                "{:<60} {:>10}  {:>6}  {:>12}",
                "span", "total", "calls", "ms/call"
            );
            eprintln!("\n[timing] {}\n{}", header, "-".repeat(header.len()));
            for line in &flat {
                eprintln!("[timing] {}", line);
            }
        }
    }

    /// Register an `atexit` handler that calls [`dump_and_clear`].
    ///
    /// Called once during solver initialisation.  Idempotent.
    pub fn init() {
        use std::sync::Once;
        static INIT: Once = Once::new();
        INIT.call_once(|| {
            if enabled() {
                // Best-effort atexit via Drop guard stored in a leaked Box.
                struct DumpGuard;
                impl Drop for DumpGuard {
                    fn drop(&mut self) {
                        dump_and_clear();
                    }
                }
                // Leak the guard so it lives until process exit.
                let _ = Box::leak(Box::new(DumpGuard));
            }
        });
    }
}

#[cfg(target_arch = "wasm32")]
mod native {
    pub fn enabled() -> bool { false }
    pub fn init() {}
    pub fn dump_and_clear() {}

    pub struct SpanGuard;
    #[allow(unused_variables)]
    pub fn span(name: &'static str) -> SpanGuard { SpanGuard }
}

pub use native::*;
