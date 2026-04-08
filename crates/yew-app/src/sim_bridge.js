// sim_bridge.js — main-thread helper for off-thread simulation
// Exposed as globalThis.remSim

(function () {
  // Detect the wasm-bindgen JS glue URL from already-loaded <script> tags.
  // Trunk injects it as a module script whose src contains the package name.
  function detectWasmJsUrl() {
    // Look for a <script type="module"> whose src contains "rem-yew"
    const scripts = Array.from(document.querySelectorAll("script[type=module]"));
    for (const s of scripts) {
      if (s.src && s.src.includes("rem-yew")) {
        return s.src;
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
        const wasmJsUrl = detectWasmJsUrl();
        if (!wasmJsUrl) {
          reject(new Error("sim_bridge: could not detect WASM JS URL"));
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

        // Transfer meshBytes buffer to avoid copying (zero-copy transfer)
        const transferable = meshBytes.buffer instanceof SharedArrayBuffer
          ? []
          : [meshBytes.buffer];

        worker.postMessage(
          { type: "run", configJson, meshBytes, wasmJsUrl },
          transferable,
        );
      });
    },
  };
})();
