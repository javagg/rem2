// sim_bridge.js — main-thread helper for off-thread simulation
// Exposed as globalThis.remSim

(function () {
  // Find the wasm-bindgen glue JS URL by scanning performance resource entries.
  // By the time runInWorker is called, the glue has already been fetched and
  // its entry is present. The glue filename matches "rem-yew-*.js".
  function getWasmJsUrl() {
    // Prefer the value set by store_wasm_js_url() if it points to the glue.
    const stored = globalThis.__remWasmJsUrl;
    if (stored && !stored.includes("/snippets/")) return stored;

    // Fall back to scanning performance resource entries.
    const entries = performance.getEntriesByType("resource");
    for (const e of entries) {
      if (/rem-yew[^/]*\.js$/.test(e.name) && !e.name.includes("/snippets/")) {
        return e.name;
      }
    }
    return null;
  }

  globalThis.remSim = {
    /**
     * Run a simulation in a dedicated Web Worker.
     * @param {string} configJson  - Palace JSON config
     * @param {Uint8Array} meshBytes - raw mesh bytes
     * @returns {Promise<any>}  - resolves with the JsValue from run_simulation
     */
    runInWorker(configJson, meshBytes) {
      return new Promise((resolve, reject) => {
        const wasmJsUrl = getWasmJsUrl();
        if (!wasmJsUrl) {
          reject(new Error("sim_bridge: WASM JS URL not set (init_logger not called yet?)"));
          return;
        }

        const worker = new Worker(new URL("sim_worker.js", location.href), {
          type: "module",
        });

        worker.addEventListener("message", (event) => {
          const { type, value, message } = event.data ?? {};
          if (type === "result") {
            worker.terminate();
            resolve(value);
          } else if (type === "error") {
            worker.terminate();
            reject(new Error(message ?? "unknown worker error"));
          }
        });

        worker.addEventListener("error", (event) => {
          worker.terminate();
          reject(new Error(event.message ?? "worker uncaught error"));
        });

        // meshBytes is a view into WASM linear memory (non-transferable).
        // Copy it into a fresh ArrayBuffer so we can transfer it zero-copy.
        const meshCopy = meshBytes.slice();   // new Uint8Array with own buffer
        worker.postMessage(
          { type: "run", configJson, meshBytes: meshCopy, wasmJsUrl },
          [meshCopy.buffer],
        );
      });
    },
  };
})();
