mod examples;
mod solver;

use examples::{ExampleStatus, EXAMPLES};
use solver::SimResult;
use wasm_bindgen_futures::spawn_local;
use web_sys::HtmlSelectElement;
use yew::prelude::*;

#[derive(Clone, Copy, PartialEq)]
enum Tab {
    Config,
    Source,
}

#[function_component(App)]
fn app() -> Html {
    let selected = use_state(|| "spheres".to_string());
    let running = use_state(|| false);
    let log_text = use_state(|| String::new());
    let result = use_state(|| None::<SimResult>);
    let active_tab = use_state(|| Tab::Config);

    let example = examples::find_example(&selected).unwrap();

    // Example selection handler
    let on_select = {
        let selected = selected.clone();
        let result = result.clone();
        Callback::from(move |e: Event| {
            let el: HtmlSelectElement = e.target_unchecked_into();
            selected.set(el.value());
            result.set(None);
        })
    };

    // Run simulation
    let on_run = {
        let selected = selected.clone();
        let running = running.clone();
        let log_text = log_text.clone();
        let result = result.clone();
        Callback::from(move |_: MouseEvent| {
            let key = (*selected).clone();
            let running = running.clone();
            let log_text = log_text.clone();
            let result = result.clone();
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

    let is_unimplemented = example.status == ExampleStatus::Unimplemented;
    let btn_disabled = *running || is_unimplemented;
    let btn_text = if *running {
        "Running..."
    } else if is_unimplemented {
        "Not Implemented"
    } else {
        "Run Simulation"
    };

    let code_text = match *active_tab {
        Tab::Config => example.config_json,
        Tab::Source => example.source_code,
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
                            disabled={*running}>
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

                    <button type="button" class="run-btn"
                        onclick={on_run}
                        disabled={btn_disabled}>
                        {btn_text}
                    </button>

                    if is_unimplemented {
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
                </div>

                <div class="code-panel">
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
                    </div>
                    <div class="code-viewer">
                        <pre>{code_text}</pre>
                    </div>
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
