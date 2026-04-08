// sim_worker.js — runs rem-wasm simulation off the main thread
// Loaded by sim_bridge.js via new Worker(..., { type: "module" })
//
// Protocol:
//   main → worker:  { type: "run", configJson, meshBytes, wasmJsUrl }
//   worker → main:  { type: "result", value }  (JsValue from run_simulation)
//   worker → main:  { type: "error",  message }

self.addEventListener("message", async (event) => {
  const { type, configJson, meshBytes, wasmJsUrl } = event.data ?? {};
  if (type !== "run") return;

  try {
    // Dynamically import the wasm-bindgen JS glue and initialise the WASM module.
    const wasmModule = await import(wasmJsUrl);
    await wasmModule.default();        // calls __wbg_init()
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
