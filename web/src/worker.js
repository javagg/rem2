import init, { init_panic_hook, init_logger, run_simulation } from '../pkg/rem_wasm.js';

self.onmessage = async (event) => {
    const { rank, size, config, mesh } = event.data;

    let wasmInstance;
    self.jsmpi = {
        Init: () => 0,
        Finalize: () => 0,
        Comm_size: (comm, ptr) => {
            if (wasmInstance) {
                new Int32Array(wasmInstance.memory.buffer)[ptr / 4] = size;
            }
            return 0;
        },
        Comm_rank: (comm, ptr) => {
            if (wasmInstance) {
                new Int32Array(wasmInstance.memory.buffer)[ptr / 4] = rank;
            }
            return 0;
        },
        Barrier: () => {
            console.log(`Worker ${rank} at barrier`);
            return 0;
        },
        Allreduce: (sendptr, recvptr, count, datatype, op) => 0,
        Bcast: (ptr, count, datatype, root) => 0,
        Send: () => 0,
        Recv: () => 0,
    };

    wasmInstance = await init();
    init_panic_hook();
    init_logger();

    try {
        const result = run_simulation(config, mesh);
        self.postMessage({ type: 'RESULT', rank, result });
    } catch (e) {
        self.postMessage({ type: 'ERROR', rank, error: e.toString() });
    }
};
