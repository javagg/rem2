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

  // Make each rank produce distinct, ordered output for UI verification.
  await sleep(60 * (rank + 1));
  console.log(`[rank ${rank}] local phase complete`);

  if (typeof globalThis.__jsmpi_barrier_timeout === "function") {
    const ok = globalThis.__jsmpi_barrier_timeout(rank, size, 5000);
    console.log(`[rank ${rank}] barrier ${ok ? "passed" : "timed-out"}`);
  }

  await sleep(30);
  console.log(`[rank ${rank}] finished`);

  if (typeof globalThis.__jsmpi_mark_finished === "function") {
    globalThis.__jsmpi_mark_finished();
  }
}
