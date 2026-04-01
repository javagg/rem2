<template>
  <div id="app">
    <h1>REM EM Solver Demo (Vue 3 + Vite)</h1>

    <div class="main-layout">
      <div class="controls-panel">
        <div class="control-group">
          <label for="example-select">Example:</label>
          <select id="example-select" v-model="selectedExample" :disabled="running" title="Select a simulation example">
            <option value="spheres">Spheres (Electrostatic)</option>
            <option value="rings">Rings (Magnetostatic)</option>
          </select>
        </div>

        <div class="control-group">
          <label for="worker-count">MPI Workers:</label>
          <input id="worker-count" type="number" v-model="size" min="1" max="8" :disabled="running" title="Number of MPI workers">
        </div>

        <button type="button" class="run-btn" @click="runSim" :disabled="running">
          {{ running ? 'Running...' : 'Run Simulation' }}
        </button>

        <div class="results-panel" v-if="finalResult">
          <h3>Summary Result:</h3>
          <p><strong>Energy:</strong> {{ (finalResult.energy * 1e12).toFixed(6) }} pJ</p>
          <p><strong>Nodes:</strong> {{ finalResult.phi.length }}</p>
          <div v-if="finalResult.e_field">
            <p><strong>Max |E|:</strong> {{ maxMagnitude(finalResult.e_field).toFixed(4) }} V/m</p>
          </div>
          <div v-if="finalResult.b_field">
            <p><strong>Max |B|:</strong> {{ maxMagnitude(finalResult.b_field).toFixed(4) }} T</p>
          </div>
        </div>

        <div class="log-panel">
          <h3>Logs:</h3>
          <pre>{{ log }}</pre>
        </div>
      </div>

      <div class="code-panel">
        <div class="tabs">
          <button type="button" :class="{ active: activeTab === 'config' }" @click="activeTab = 'config'">Palace Config</button>
          <button type="button" :class="{ active: activeTab === 'source' }" @click="activeTab = 'source'">Test Source (Rust)</button>
        </div>
        <div class="code-viewer">
          <pre><code ref="codeBlock" :class="activeTab === 'config' ? 'language-json' : 'language-rust'">{{ currentCode }}</code></pre>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, watch, onMounted, computed, nextTick } from 'vue';
import init, { get_spheres_mesh, get_rings_mesh } from '../pkg/rem_wasm.js';
import hljs from 'highlight.js';
import 'highlight.js/styles/github-dark.css';

const spheresConfig = {
  Problem: { Type: 'Electrostatic', Output: '.' },
  Model: { Mesh: 'coaxial_2d.msh', L0: 0.001 },
  Domains: { Materials: [{ Attributes: [10], Permittivity: 1.0 }] },
  Boundaries: {
    Ground: { Attributes: [2] },
    Terminal: [{ Index: 1, Attributes: [1] }]
  },
  Solver: { Linear: { Tol: 1e-10, MaxIter: 500 } }
};

const ringsConfig = {
  Problem: { Type: 'Magnetostatic', Output: '.' },
  Model: { Mesh: 'slab_2d.msh', L0: 0.001 },
  Domains: {
    Materials: [
      { Attributes: [10], Permeability: 1000.0 },
      { Attributes: [20], Permeability: 1.0 }
    ]
  },
  Boundaries: {
    Ground: { Attributes: [1] },
    SurfaceCurrent: [{ Index: 1, Attributes: [2], Direction: '+Y' }]
  },
  Solver: { Linear: { Tol: 1e-10, MaxIter: 500 } }
};

const selectedExample = ref('spheres');
const size = ref(1);
const running = ref(false);
const log = ref('');
const results = ref([]);
const finalResult = ref(null);
const activeTab = ref('config');
const codeBlock = ref(null);

const spheresSource = ref('');
const ringsSource = ref('');

const currentCode = computed(() => {
  if (activeTab.value === 'config') {
    return JSON.stringify(selectedExample.value === 'spheres' ? spheresConfig : ringsConfig, null, 2);
  } else {
    return selectedExample.value === 'spheres' ? spheresSource.value : ringsSource.value;
  }
});

watch(currentCode, () => {
  nextTick(() => {
    if (codeBlock.value) {
      hljs.highlightElement(codeBlock.value);
    }
  });
});

onMounted(async () => {
  await init();
  spheresSource.value = `// Simplified view of palace_spheres.rs\n#[test]\nfn solve_spheres() {\n    let mesh = annular_msh(1.0, 4.0, 10, 32, 1, 2, 10);\n    let phi = solve_one(&cfg, &mesh, &dm, Some(1), 1.0, &comm).unwrap();\n    // ... verification\n}`;
  ringsSource.value = `// Simplified view of palace_rings.rs\n#[test]\nfn solve_rings() {\n    let mesh = rect_bimaterial_msh(1.0, 1.0, 20, 20, 1, 2, 10, 20);\n    let az = solve_one(&cfg, &mesh, &dm, Some(1), &comm).unwrap();\n    // ... verification\n}`;

  if (codeBlock.value) {
    hljs.highlightElement(codeBlock.value);
  }
});

function maxMagnitude(vectors) {
  let max = 0;
  for (const v of vectors) {
    const mag = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (mag > max) max = mag;
  }
  return max;
}

async function runSim() {
  running.value = true;
  log.value = `Starting simulation: ${selectedExample.value}...\n`;
  results.value = [];
  finalResult.value = null;

  const workers = [];
  for (let i = 0; i < size.value; i++) {
    const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = (e) => {
      if (e.data.type === 'RESULT') {
        log.value += `Rank ${e.data.rank} finished.\n`;
        results.value.push(e.data.result);
        if (results.value.length === size.value) {
          log.value += 'All ranks completed.\n';
          finalResult.value = results.value[0];
          running.value = false;
        }
      } else if (e.data.type === 'ERROR') {
        log.value += `Rank ${e.data.rank} ERROR: ${e.data.error}\n`;
        running.value = false;
      }
    };
    workers.push(worker);
  }

  const configStr = JSON.stringify(selectedExample.value === 'spheres' ? spheresConfig : ringsConfig);
  const mesh = selectedExample.value === 'spheres' ? get_spheres_mesh() : get_rings_mesh();

  workers.forEach((w, i) => {
    w.postMessage({ rank: i, size: size.value, config: configStr, mesh });
  });
}
</script>

<style>
.main-layout {
  display: flex;
  gap: 20px;
  height: 80vh;
}

.controls-panel {
  flex: 1;
  display: flex;
  flex-direction: column;
  gap: 15px;
  padding: 20px;
  border: 1px solid #ddd;
  border-radius: 8px;
  background: #f8f9fa;
  overflow-y: auto;
}

.code-panel {
  flex: 2;
  display: flex;
  flex-direction: column;
  border: 1px solid #333;
  border-radius: 8px;
  background: #1e1e1e;
  overflow: hidden;
}

.control-group {
  display: flex;
  flex-direction: column;
  gap: 5px;
}

.run-btn {
  padding: 10px;
  background: #007bff;
  color: white;
  border: none;
  border-radius: 4px;
  cursor: pointer;
  font-weight: bold;
}

.run-btn:disabled {
  background: #ccc;
}

.results-panel {
  padding: 10px;
  background: #e7f3ff;
  border: 1px solid #007bff;
  border-radius: 4px;
}

.log-panel pre {
  background: #eee;
  padding: 10px;
  font-size: 12px;
  border-radius: 4px;
  max-height: 200px;
  overflow-y: auto;
}

.tabs {
  display: flex;
  background: #333;
}

.tabs button {
  padding: 10px 20px;
  background: none;
  border: none;
  color: #aaa;
  cursor: pointer;
}

.tabs button.active {
  background: #1e1e1e;
  color: white;
  border-bottom: 2px solid #007bff;
}

.code-viewer {
  flex: 1;
  overflow: auto;
}

.code-viewer pre {
  margin: 0;
  padding: 15px;
}
</style>
