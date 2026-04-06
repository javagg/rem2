mod examples;
mod solver;

use std::collections::BTreeMap;
use std::rc::Rc;

use examples::{ExampleStatus, EXAMPLES};
use solver::SimResult;
use wasm_bindgen_futures::spawn_local;
use web_sys::{HtmlInputElement, HtmlSelectElement};
use yew::prelude::*;

#[cfg(target_arch = "wasm32")]
type MpiJobHandle = Rc<jsmpi::launcher::Job>;

#[cfg(not(target_arch = "wasm32"))]
type MpiJobHandle = Rc<()>;

#[derive(Clone, Copy, PartialEq)]
enum Tab {
    Config,
    Source,
    Mesh,
}

#[derive(Clone, PartialEq)]
struct RankLog {
    rank: u32,
    text: String,
}

const MIN_RANKS: u32 = 1;
const MAX_RANKS: u32 = 7;

fn clamp_ranks(n: u32) -> u32 {
    n.clamp(MIN_RANKS, MAX_RANKS)
}

#[cfg(target_arch = "wasm32")]
fn parse_rank_from_text(text: &str) -> Option<u32> {
    if let Some(rest) = text.strip_prefix("[rank ") {
        return rest
            .split(']')
            .next()
            .and_then(|s| s.trim().parse::<u32>().ok());
    }
    if let Some(rest) = text.strip_prefix("rank ") {
        return rest
            .split(' ')
            .next()
            .and_then(|s| s.trim().parse::<u32>().ok());
    }
    None
}

fn append_line(existing: &str, line: &str) -> String {
    if existing.is_empty() {
        line.to_string()
    } else {
        format!("{}\n{}", existing, line)
    }
}

fn strip_rank_prefix<'a>(text: &'a str, rank: u32) -> &'a str {
    let bracket_prefix = format!("[rank {}] ", rank);
    if let Some(stripped) = text.strip_prefix(&bracket_prefix) {
        return stripped;
    }

    let plain_prefix = format!("rank {} ", rank);
    if let Some(stripped) = text.strip_prefix(&plain_prefix) {
        return stripped;
    }

    text
}

fn is_console_kind(kind: &str) -> bool {
    matches!(kind, "debug" | "info" | "log" | "trace" | "warn" | "error")
}

#[function_component(App)]
fn app() -> Html {
    let selected = use_state(|| "spheres".to_string());
    let running = use_state(|| false);
    let mpi_enabled = use_state(|| false);
    let rank_count = use_state(|| 4u32);
    let mpi_running = use_state(|| false);
    let mpi_status = use_state(|| "Idle".to_string());
    let mpi_job: UseStateHandle<Option<MpiJobHandle>> = use_state(|| None);
    let rank_logs = use_state(BTreeMap::<u32, RankLog>::new);
    let rank_logs_ref = use_mut_ref(BTreeMap::<u32, RankLog>::new);
    let log_text = use_state(|| String::new());
    let result = use_state(|| None::<SimResult>);
    let active_tab = use_state(|| Tab::Config);
    let code_panel_collapsed = use_state(|| true);

    let example = examples::find_example(&selected).unwrap();

    // Example selection handler
    let on_select = {
        let selected = selected.clone();
        let result = result.clone();
        let rank_logs = rank_logs.clone();
        let rank_logs_ref = rank_logs_ref.clone();
        Callback::from(move |e: Event| {
            let el: HtmlSelectElement = e.target_unchecked_into();
            selected.set(el.value());
            result.set(None);
            *rank_logs_ref.borrow_mut() = BTreeMap::new();
            rank_logs.set(BTreeMap::new());
        })
    };

    let on_mpi_toggle = {
        let mpi_enabled = mpi_enabled.clone();
        let mpi_running = mpi_running.clone();
        let mpi_status = mpi_status.clone();
        let mpi_job = mpi_job.clone();
        Callback::from(move |e: Event| {
            let el: HtmlInputElement = e.target_unchecked_into();
            let checked = el.checked();
            if !checked {
                mpi_job.set(None);
                mpi_running.set(false);
                mpi_status.set("Idle".to_string());
            }
            mpi_enabled.set(checked);
        })
    };

    let on_rank_change = {
        let rank_count = rank_count.clone();
        Callback::from(move |e: Event| {
            let el: HtmlInputElement = e.target_unchecked_into();
            if let Ok(v) = el.value().parse::<u32>() {
                rank_count.set(clamp_ranks(v));
            }
        })
    };

    let on_stop = {
        let mpi_job = mpi_job.clone();
        let mpi_running = mpi_running.clone();
        let mpi_status = mpi_status.clone();
        let log_text = log_text.clone();
        Callback::from(move |_: MouseEvent| {
            mpi_job.set(None);
            mpi_running.set(false);
            mpi_status.set("Stopped".to_string());
            log_text.set(append_line(&log_text, "MPI job stopped."));
        })
    };

    // Run simulation
    let on_run = {
        let selected = selected.clone();
        let running = running.clone();
        let mpi_enabled = mpi_enabled.clone();
        let rank_count = rank_count.clone();
        let mpi_running = mpi_running.clone();
        let mpi_status = mpi_status.clone();
        let mpi_job = mpi_job.clone();
        let rank_logs = rank_logs.clone();
        let rank_logs_ref = rank_logs_ref.clone();
        let log_text = log_text.clone();
        let result = result.clone();
        Callback::from(move |_: MouseEvent| {
            let key = (*selected).clone();
            let use_mpi = *mpi_enabled;
            let running = running.clone();
            let mpi_running = mpi_running.clone();
            let mpi_status = mpi_status.clone();
            let mpi_job = mpi_job.clone();
            let rank_logs = rank_logs.clone();
            let rank_logs_ref = rank_logs_ref.clone();
            let log_text = log_text.clone();
            let result = result.clone();

            if use_mpi {
                let ranks = clamp_ranks(*rank_count);
                mpi_job.set(None);

                let mut init_map = BTreeMap::new();
                for r in 0..ranks {
                    init_map.insert(
                        r,
                        RankLog {
                            rank: r,
                            text: "Waiting for output...".to_string(),
                        },
                    );
                }
                *rank_logs_ref.borrow_mut() = init_map.clone();
                rank_logs.set(init_map);
                result.set(None);
                mpi_running.set(true);
                mpi_status.set("Running".to_string());
                log_text.set(format!(
                    "Starting MPI job: example={}, ranks={}\n",
                    key, ranks
                ));

                let module_url = format!("/mpi_rank_probe.js?example={}", key);
                #[cfg(target_arch = "wasm32")]
                {
                    let rank_logs_on_log = rank_logs.clone();
                    let rank_logs_ref_on_log = rank_logs_ref.clone();
                    let log_text_on_log = log_text.clone();
                    let mpi_status_on_state = mpi_status.clone();
                    let mpi_running_on_state = mpi_running.clone();
                    let log_text_on_state = log_text.clone();
                    let mpi_running_on_complete = mpi_running.clone();
                    let mpi_status_on_complete = mpi_status.clone();
                    let log_text_on_complete = log_text.clone();

                    let job = jsmpi::launcher::create_job(
                        "/worker.js",
                        &module_url,
                        ranks,
                        move |kind: String, text: String| {
                            let line = format!("[{}] {}", kind, text);
                            log_text_on_log.set(append_line(&log_text_on_log, &line));

                            if is_console_kind(&kind) {
                                if let Some(rank) = parse_rank_from_text(&text) {
                                    let mut map = rank_logs_ref_on_log.borrow_mut();
                                    let entry = map.entry(rank).or_insert(RankLog {
                                        rank,
                                        text: String::new(),
                                    });
                                    let cleaned = strip_rank_prefix(&text, rank);
                                    entry.text = append_line(&entry.text, cleaned);
                                    rank_logs_on_log.set(map.clone());
                                }
                            }
                        },
                        move |state: String| {
                            mpi_status_on_state.set(state.clone());
                            if state == "finished" {
                                mpi_running_on_state.set(false);
                                log_text_on_state
                                    .set(append_line(&log_text_on_state, "MPI job finished."));
                            }
                        },
                        move |finished, total| {
                            mpi_running_on_complete.set(false);
                            mpi_status_on_complete.set("finished".to_string());
                            log_text_on_complete.set(append_line(
                                &log_text_on_complete,
                                &format!("MPI completion callback: {}/{} ranks done.", finished, total),
                            ));
                        },
                    );
                    mpi_job.set(Some(Rc::new(job)));
                }

                #[cfg(not(target_arch = "wasm32"))]
                {
                    let _ = module_url;
                    mpi_running.set(false);
                    mpi_status.set("Unsupported target".to_string());
                    log_text.set(append_line(
                        &log_text,
                        "MPI launcher is available only on wasm32 target.",
                    ));
                }
                return;
            }

            running.set(true);
            log_text.set(format!("Starting simulation: {}...\n", key));
            result.set(None);

            spawn_local(async move {
                // Yield to let UI update before heavy computation
                gloo::timers::future::sleep(std::time::Duration::from_millis(10)).await;

                match solver::run_example(&key) {
                    Ok(r) => {
                        log_text.set(format!(
                            "Starting simulation: {}...\nSimulation completed.\n",
                            key
                        ));
                        result.set(Some(r));
                    }
                    Err(e) => {
                        log_text.set(format!(
                            "Starting simulation: {}...\nERROR: {}\n",
                            key, e
                        ));
                    }
                }
                running.set(false);
            });
        })
    };

    let on_tab_config = {
        let active_tab = active_tab.clone();
        Callback::from(move |_: MouseEvent| active_tab.set(Tab::Config))
    };

    let on_tab_source = {
        let active_tab = active_tab.clone();
        Callback::from(move |_: MouseEvent| active_tab.set(Tab::Source))
    };

    let on_tab_mesh = {
        let active_tab = active_tab.clone();
        Callback::from(move |_: MouseEvent| active_tab.set(Tab::Mesh))
    };

    let on_toggle_code_panel = {
        let code_panel_collapsed = code_panel_collapsed.clone();
        Callback::from(move |_: MouseEvent| code_panel_collapsed.set(!*code_panel_collapsed))
    };

    let is_unimplemented = example.status == ExampleStatus::Unimplemented;
    let btn_disabled = if *mpi_enabled {
        *mpi_running
    } else {
        *running || is_unimplemented
    };
    let btn_text = if *mpi_enabled && *mpi_running {
        "MPI Running..."
    } else if *running {
        "Running..."
    } else if is_unimplemented {
        "Not Implemented"
    } else {
        "Run Simulation"
    };

    let code_text = match *active_tab {
        Tab::Config => example.config_json,
        Tab::Source => example.source_code,
        Tab::Mesh => "3D mesh preview will be added in a later iteration.\n\nPlanned capabilities:\n- Rotate / pan / zoom\n- Boundary / domain coloring\n- Rank partition overlay in MPI mode\n- Probe point and field value tooltip",
    };

    html! {
        <div id="app">
            <h1>{"REM EM Solver Demo (Yew + WASM)"}</h1>

            <div class="main-layout">
                <div class="controls-panel">
                    <div class="control-group">
                        <label for="example-select">{"Example:"}</label>
                        <select id="example-select"
                            onchange={on_select}
                            disabled={*running || *mpi_running}>
                            { for EXAMPLES.iter().map(|ex| {
                                html! {
                                    <option value={ex.key}
                                        selected={ex.key == selected.as_str()}>
                                        {ex.label}
                                    </option>
                                }
                            })}
                        </select>
                    </div>

                    <div class="control-group mpi-toggle-row">
                        <label for="mpi-mode">{"MPI Mode"}</label>
                        <input
                            id="mpi-mode"
                            type="checkbox"
                            checked={*mpi_enabled}
                            onchange={on_mpi_toggle}
                            disabled={*running || *mpi_running}
                        />
                    </div>

                    if *mpi_enabled {
                        <div class="control-group">
                            <label for="rank-count">{"MPI Ranks (1-7):"}</label>
                            <input
                                id="rank-count"
                                type="number"
                                min="1"
                                max="7"
                                step="1"
                                value={rank_count.to_string()}
                                onchange={on_rank_change}
                                disabled={*mpi_running}
                            />
                        </div>
                        <p class="info-text">{format!("MPI Status: {}", *mpi_status)}</p>
                    }

                    <button type="button" class="run-btn"
                        onclick={on_run}
                        disabled={btn_disabled}>
                        {btn_text}
                    </button>

                    if *mpi_enabled {
                        <button type="button" class="stop-btn"
                            onclick={on_stop}
                            disabled={!*mpi_running}>
                            {"Stop MPI Job"}
                        </button>
                    }

                    if is_unimplemented && !*mpi_enabled {
                        <p class="warning-text">
                            {format!("The {} solver is not yet implemented.", example.problem_type)}
                        </p>
                    }

                    if let Some(r) = &*result {
                        <div class="results-panel">
                            <h3>{"Summary Result:"}</h3>
                            <p><strong>{"Energy: "}</strong>{format!("{:.6} pJ", r.energy * 1e12)}</p>
                            <p><strong>{"Nodes: "}</strong>{r.node_count}</p>
                            if let Some(max_e) = r.max_e {
                                <p><strong>{"Max |E|: "}</strong>{format!("{:.4} V/m", max_e)}</p>
                            }
                            if let Some(max_b) = r.max_b {
                                <p><strong>{"Max |B|: "}</strong>{format!("{:.4} T", max_b)}</p>
                            }
                        </div>
                    }

                    <div class="log-panel">
                        <h3>{"Logs:"}</h3>
                        <pre>{&*log_text}</pre>
                    </div>

                    if *mpi_enabled {
                        <div class="rank-logs-panel">
                            <h3>{"Per-rank Output:"}</h3>
                            <div class="rank-panel-list">
                                { for rank_logs.values().map(|entry| {
                                    let line_count = entry.text.lines().filter(|s| !s.trim().is_empty()).count();
                                    html! {
                                        <section class="rank-output-panel">
                                            <header class="rank-output-header">
                                                <h4>{format!("Rank {}", entry.rank)}</h4>
                                                <span class="rank-output-meta">{format!("{} lines", line_count)}</span>
                                            </header>
                                            <div class="rank-output-body">
                                                <pre>{&entry.text}</pre>
                                            </div>
                                        </section>
                                    }
                                }) }

                                { if rank_logs.is_empty() {
                                    html! {
                                        <section class="rank-output-panel rank-output-empty">
                                            <header class="rank-output-header">
                                                <h4>{"No rank output yet"}</h4>
                                            </header>
                                            <div class="rank-output-body">
                                                <pre>{"Run an MPI job to populate per-rank output panels."}</pre>
                                            </div>
                                        </section>
                                    }
                                } else {
                                    html! {}
                                }}
                            </div>
                        </div>
                    }
                </div>

                <div class="code-panel">
                    <div class="code-panel-header">
                        <h3>{"Config & View"}</h3>
                        <button type="button"
                            class="collapse-btn"
                            onclick={on_toggle_code_panel}>
                            {if *code_panel_collapsed { "Expand" } else { "Collapse" }}
                        </button>
                    </div>

                    if !*code_panel_collapsed {
                        <>
                            <div class="tabs">
                                <button type="button"
                                    class={if *active_tab == Tab::Config { "active" } else { "" }}
                                    onclick={on_tab_config}>
                                    {"Palace Config"}
                                </button>
                                <button type="button"
                                    class={if *active_tab == Tab::Source { "active" } else { "" }}
                                    onclick={on_tab_source}>
                                    {"Test Source (Rust)"}
                                </button>
                                <button type="button"
                                    class={if *active_tab == Tab::Mesh { "active" } else { "" }}
                                    onclick={on_tab_mesh}>
                                    {"Mesh 3D (Soon)"}
                                </button>
                            </div>
                            <div class="code-viewer">
                                <pre>{code_text}</pre>
                            </div>
                        </>
                    } else {
                        <div class="collapsed-hint">
                            {"Config panel is collapsed by default. Expand when you need to inspect JSON, source, or upcoming 3D mesh preview."}
                        </div>
                    }
                </div>
            </div>
        </div>
    }
}

fn main() {
    console_error_panic_hook::set_once();
    console_log::init_with_level(log::Level::Info).ok();
    yew::Renderer::<App>::new().render();
}
