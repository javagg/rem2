// sim_worker.js — runs rem-wasm simulation off the main thread
// Loaded by sim_bridge.js via new Worker(..., { type: "module" })
//
// Protocol:
//   main → worker:  { type: "run", configJson, meshBytes, wasmJsUrl }
//   worker → main:  { type: "result", value }  (JsValue from run_simulation)
//   worker → main:  { type: "error",  message }
//   worker → main:  { type: "log",    level, text }  (streamed log lines)

// Intercept console output and forward to main thread as log messages.
// This must happen before WASM is loaded so all log::info! calls are captured.
(function patchConsole() {
  const levels = ["log", "info", "warn", "error", "debug"];
  for (const level of levels) {
    const orig = console[level].bind(console);
    console[level] = (...args) => {
      orig(...args);
      try {
        self.postMessage({ type: "log", level, text: args.map(String).join(" ") });
      } catch (_) { /* structured-clone error on exotic objects — ignore */ }
    };
  }
})();

self.addEventListener("message", async (event) => {
  const { type, configJson, meshBytes, wasmJsUrl } = event.data ?? {};
  if (type !== "run") return;

  try {
    // Dynamically import the wasm-bindgen JS glue and initialise the WASM module.
    // Must pass the explicit WASM binary URL so the worker can fetch it
    // (relative URLs resolve against the worker script, not the glue file).
    const wasmBinaryUrl = wasmJsUrl.replace(/\.js$/, "_bg.wasm");
    const wasmModule = await import(wasmJsUrl);

    // wasm-bindgen --target web exports the async init as the default export.
    // Fall back to named 'init' or 'initSync' for other targets.
    const initFn = wasmModule.default ?? wasmModule.init ?? wasmModule.initSync;
    if (typeof initFn !== "function") {
      throw new Error(
        `sim_worker: no init export found in ${wasmJsUrl} ` +
        `(got: ${Object.keys(wasmModule).join(", ")})`
      );
    }
    await initFn(wasmBinaryUrl);
    wasmModule.init_panic_hook();
    wasmModule.init_logger();

    const result = wasmModule.run_simulation(configJson, meshBytes);
    self.postMessage({ type: "result", value: result });
  } catch (err) {
    self.postMessage({
      type: "error",
      message: err instanceof Error ? (err.stack ?? err.message) : String(err),
    });
  }
});
