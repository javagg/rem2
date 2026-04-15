const inbox = [];
const RECV_INDEX = 0;
const BARRIER_INDEX = 1;
const PROTOCOL_VERSION = 1;
let controlSignals = null;

function parseErrorCode(text) {
  const prefix = String(text || "").split(":", 1)[0]?.trim() || "";
  if (/^[A-Z_]+$/.test(prefix) && prefix.includes("_")) {
    return prefix;
  }
  return null;
}

function postStructured(type, payload = {}) {
  globalThis.postMessage({
    type,
    jobId: globalThis.__jsmpi_job_id,
    rank: globalThis.__jsmpi_rank,
    ...payload,
  });
}

function isInteger(value) {
  return Number.isInteger(value);
}

function isValidEnvelope(envelope) {
  if (!envelope || typeof envelope !== "object") {
    return false;
  }

  if (envelope.protocol_version !== PROTOCOL_VERSION) {
    return false;
  }

  if (!isInteger(envelope.src) || !isInteger(envelope.dst) || !isInteger(envelope.tag)) {
    return false;
  }

  if (envelope.src < 0 || envelope.dst < 0) {
    return false;
  }

  if (!Array.isArray(envelope.payload)) {
    return false;
  }

  return true;
}

function matches(envelope, source, tag) {
  const sourceOk = source == null || envelope.src === source;
  const tagOk = tag == null || envelope.tag === tag;
  return sourceOk && tagOk;
}

function takeMatchingEnvelope(source, tag) {
  const index = inbox.findIndex((envelope) => matches(envelope, source, tag));
  if (index === -1) {
    return null;
  }

  return inbox.splice(index, 1)[0] ?? null;
}

function wake(index) {
  if (!controlSignals) {
    return;
  }

  Atomics.add(controlSignals, index, 1);
  Atomics.notify(controlSignals, index, 1);
}

function formatArgs(args) {
  return args
    .map((value) => {
      if (typeof value === "string") {
        return value;
      }
      try {
        return JSON.stringify(value);
      } catch {
        return String(value);
      }
    })
    .join(" ");
}

const originalConsole = {
  debug: console.debug.bind(console),
  info: console.info.bind(console),
  log: console.log.bind(console),
  trace: console.trace.bind(console),
  warn: console.warn.bind(console),
  error: console.error.bind(console),
};

for (const level of ["debug", "info", "log", "trace", "warn", "error"]) {
  console[level] = (...args) => {
    postStructured("jsmpi:log", {
      level,
      text: formatArgs(args),
      phase: "runtime",
      event: "console",
    });
    originalConsole[level](...args);
  };
}

globalThis.addEventListener("error", (event) => {
  postStructured("jsmpi:error", {
    text: event.message,
    errorCode: parseErrorCode(event.message) || "WORKER_RUNTIME_ERROR",
    phase: "runtime",
    event: "error",
  });
});

globalThis.__jsmpi_send = function __jsmpi_send(envelope) {
  if (!isValidEnvelope(envelope)) {
    postStructured("jsmpi:error", {
      text: "PROTO_INVALID_ENVELOPE: outbound envelope validation failed",
      errorCode: "PROTO_INVALID_ENVELOPE",
      phase: "protocol",
      event: "send_validate",
    });
    return;
  }

  postStructured("jsmpi:send", { envelope, phase: "transport", event: "send" });
};

globalThis.__jsmpi_recv = function __jsmpi_recv(source, tag) {
  return takeMatchingEnvelope(source, tag);
};

globalThis.__jsmpi_recv_blocking = function __jsmpi_recv_blocking(source, tag) {
  while (true) {
    const envelope = takeMatchingEnvelope(source, tag);
    if (envelope) {
      return envelope;
    }

    if (!controlSignals) {
      return null;
    }

    const version = Atomics.load(controlSignals, RECV_INDEX);
    const afterCheck = takeMatchingEnvelope(source, tag);
    if (afterCheck) {
      return afterCheck;
    }

    Atomics.wait(controlSignals, RECV_INDEX, version);
  }
};

globalThis.__jsmpi_barrier_timeout = function __jsmpi_barrier_timeout(rank, size, timeoutMs) {
  const version = controlSignals ? Atomics.load(controlSignals, BARRIER_INDEX) : 0;
  postStructured("jsmpi:barrier", { rank, size, phase: "collective", event: "barrier_wait" });

  if (!controlSignals) {
    return true;
  }

  const timeout = typeof timeoutMs === "number" ? timeoutMs : undefined;
  const result = Atomics.wait(controlSignals, BARRIER_INDEX, version, timeout);
  return result !== "timed-out";
};

globalThis.__jsmpi_barrier_begin = function __jsmpi_barrier_begin(rank, size) {
  const version = controlSignals ? Atomics.load(controlSignals, BARRIER_INDEX) : 0;
  postStructured("jsmpi:barrier", { rank, size, phase: "collective", event: "barrier_begin" });
  return version;
};

globalThis.__jsmpi_barrier_ready = function __jsmpi_barrier_ready(version) {
  if (!controlSignals) {
    return true;
  }
  return Atomics.load(controlSignals, BARRIER_INDEX) !== version;
};

globalThis.__jsmpi_mark_finished = function __jsmpi_mark_finished() {
  postStructured("jsmpi:finished", { rank: globalThis.__jsmpi_rank, phase: "lifecycle", event: "finished" });
};

globalThis.addEventListener("message", async (event) => {
  const message = event.data ?? {};

  if (message.type === "jsmpi:init") {
    globalThis.__jsmpi_job_id = message.jobId || "job-unknown";
    globalThis.__jsmpi_rank = message.rank;
    globalThis.__jsmpi_size = message.size;
    globalThis.__jsmpi_collective_timeout_ms = Number.isFinite(message.collectiveTimeoutMs)
      ? message.collectiveTimeoutMs
      : undefined;
    globalThis.__jsmpi_collective_max_retries = Number.isFinite(message.collectiveMaxRetries)
      ? message.collectiveMaxRetries
      : undefined;
    controlSignals = message.controlBuffer ? new Int32Array(message.controlBuffer) : null;
    postStructured("jsmpi:ready", { rank: message.rank, phase: "startup", event: "ready" });

    try {
      if (message.moduleUrl) {
        const mod = await import(message.moduleUrl);
        if (typeof mod.default === "function") {
          await mod.default(message.wasmUrl);
        }
        if (typeof mod.jsmpi_main === "function") {
          const maybeResult = mod.jsmpi_main();
          if (maybeResult && typeof maybeResult.then === "function") {
            await maybeResult;
          }
          return;
        }
      }
      postStructured("jsmpi:finished", { rank: message.rank, phase: "lifecycle", event: "finished" });
    } catch (error) {
      const text = error instanceof Error ? error.stack ?? error.message : String(error);
      postStructured("jsmpi:error", {
        text: error instanceof Error ? error.stack ?? error.message : String(error),
        errorCode: parseErrorCode(text) || "WORKER_RUNTIME_ERROR",
        phase: "runtime",
        event: "module_execute",
      });
    }
    return;
  }

  if (message.type === "jsmpi:deliver" && message.envelope) {
    if (!isValidEnvelope(message.envelope)) {
      postStructured("jsmpi:error", {
        text: "PROTO_INVALID_ENVELOPE: inbound envelope validation failed",
        errorCode: "PROTO_INVALID_ENVELOPE",
        phase: "protocol",
        event: "deliver_validate",
      });
      return;
    }
    inbox.push(message.envelope);
    wake(RECV_INDEX);
    return;
  }

  if (message.type === "jsmpi:barrier-release") {
    wake(BARRIER_INDEX);
  }
});