import init, { init_panic_hook, init_logger, run_simulation } from '../pkg/rem_wasm.js';

self.onmessage = async (event) => {
    const { rank, size, config, mesh } = event.data;

    self.jsmpi = {
        rank,
        size,
        Barrier: () => { console.log(`Worker ${rank} at barrier`); },
        Allreduce: (sendptr, recvptr, count, datatype, op) => {}
    };

    await init();
    init_panic_hook();
    init_logger();

    try {
        const result = run_simulation(config, mesh);
        self.postMessage({ type: 'RESULT', rank, result });
    } catch (e) {
        self.postMessage({ type: 'ERROR', rank, error: e.toString() });
    }
};
