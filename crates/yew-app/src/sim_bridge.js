// sim_bridge.js — main-thread helper for off-thread simulation
// Exposed as globalThis.remSim

(function () {
  // The WASM glue module stores its own import.meta.url on globalThis.__remWasmJsUrl
  // via store_wasm_js_url() called from init_logger(). Read it from there.
  function getWasmJsUrl() {
    return globalThis.__remWasmJsUrl ?? null;
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
