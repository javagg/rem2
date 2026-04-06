function parseExampleFromModuleUrl() {
  try {
    const current = new URL(import.meta.url);
    return current.searchParams.get("example") || "unknown";
  } catch {
    return "unknown";
  }
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

export async function jsmpi_main() {
  const rank = Number(globalThis.__jsmpi_rank ?? -1);
  const size = Number(globalThis.__jsmpi_size ?? -1);
  const example = parseExampleFromModuleUrl();

  console.log(`[rank ${rank}] job started: example=${example}, world_size=${size}`);
  console.info(`[rank ${rank}] phase=init`);

  // Emit multi-phase logs so MPI panel is informative in demo mode.
  await sleep(60 * (rank + 1));
  console.debug(`[rank ${rank}] phase=mesh-load complete`);

  await sleep(40);
  console.info(`[rank ${rank}] phase=assemble complete`);

  await sleep(40);
  console.info(`[rank ${rank}] phase=solve complete`);

  if (typeof globalThis.__jsmpi_barrier_timeout === "function") {
    const ok = globalThis.__jsmpi_barrier_timeout(rank, size, 5000);
    console.log(`[rank ${rank}] phase=barrier ${ok ? "passed" : "timed-out"}`);
  }

  await sleep(30);
  console.info(`[rank ${rank}] phase=postprocess complete`);

  await sleep(30);
  console.log(`[rank ${rank}] finished`);

  if (typeof globalThis.__jsmpi_mark_finished === "function") {
    globalThis.__jsmpi_mark_finished();
  }
}
