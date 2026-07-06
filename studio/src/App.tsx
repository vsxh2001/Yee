// Yee Studio — S.2 walking skeleton (ADR-0179).
//
// One screen: configure a driven vacuum FDTD job, run it on the in-process
// engine (yee-engine via the Tauri `run_job` command), watch progress
// events stream, and see the probe time series plotted. The plot is a
// dependency-free inline SVG — charting/3D libraries arrive with S.3.

import { Suspense, lazy, useEffect, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import { SliceHeatmap, SpectrumPlot, type Slice } from "./views";

// three.js rides its own chunk, fetched only once a result is on screen.
const FieldSurface3D = lazy(() =>
  import("./FieldSurface3D").then((m) => ({ default: m.FieldSurface3D })),
);

type Backend = "cpu" | "gpu" | "auto";

interface JobResult {
  backend: string;
  dt_s: number;
  probes: number[][];
  slice: Slice | null;
  steps_done: number;
}

interface ProgressEvent {
  step: number;
  total: number;
}

function ProbePlot({ series, dt }: { series: number[]; dt: number }) {
  if (series.length < 2) return null;
  const w = 640;
  const h = 240;
  const pad = 8;
  const max = Math.max(...series.map(Math.abs), 1e-30);
  const points = series
    .map((v, i) => {
      const x = pad + (i / (series.length - 1)) * (w - 2 * pad);
      const y = h / 2 - (v / max) * (h / 2 - pad);
      return `${x.toFixed(1)},${y.toFixed(1)}`;
    })
    .join(" ");
  return (
    <figure className="plot">
      <svg viewBox={`0 0 ${w} ${h}`} role="img" aria-label="Probe time series">
        <line x1={pad} y1={h / 2} x2={w - pad} y2={h / 2} className="axis" />
        <polyline points={points} className="trace" fill="none" />
      </svg>
      <figcaption>
        E-field probe · {series.length} steps · Δt = {(dt * 1e12).toFixed(3)} ps · peak ±
        {max.toExponential(2)}
      </figcaption>
    </figure>
  );
}

export default function App() {
  const [size, setSize] = useState(20);
  const [steps, setSteps] = useState(2000);
  const [backend, setBackend] = useState<Backend>("auto");
  const [progress, setProgress] = useState<ProgressEvent | null>(null);
  const [result, setResult] = useState<JobResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [running, setRunning] = useState(false);
  const unlisten = useRef<UnlistenFn | null>(null);

  useEffect(() => {
    return () => {
      unlisten.current?.();
    };
  }, []);

  async function run() {
    setRunning(true);
    setResult(null);
    setError(null);
    setProgress(null);
    unlisten.current?.();
    unlisten.current = await listen<ProgressEvent>("job://progress", (e) =>
      setProgress(e.payload),
    );
    const c = Math.floor(size / 2);
    try {
      const r = await invoke<JobResult>("run_job", {
        spec: {
          nx: size,
          ny: size,
          nz: size,
          dx_m: 1e-3,
          n_steps: steps,
          boundary: { kind: "pec" },
          sources: [
            {
              kind: "gaussian_ez",
              cell: [c, c, c],
              t0_steps: 12.0,
              sigma_steps: 4.0,
            },
          ],
          ports: [],
          probes: [{ component: "ez", cell: [c + Math.floor(size / 4), c, c] }],
          slice: { component: "ez", k: c },
          backend,
        },
      });
      setResult(r);
    } catch (e) {
      setError(String(e));
    } finally {
      setRunning(false);
      unlisten.current?.();
      unlisten.current = null;
    }
  }

  const pct = progress ? Math.round((100 * progress.step) / progress.total) : 0;

  return (
    <main>
      <h1>Yee Studio</h1>
      <p className="sub">
        GPU/CPU FDTD engine · in-process <code>yee-engine</code> job API
      </p>

      <section className="controls">
        <label>
          Grid
          <input
            type="number"
            min={8}
            max={120}
            value={size}
            disabled={running}
            onChange={(e) => setSize(Number(e.target.value))}
          />
          <span>× {size} × {size} cells</span>
        </label>
        <label>
          Steps
          <input
            type="number"
            min={10}
            max={100000}
            value={steps}
            disabled={running}
            onChange={(e) => setSteps(Number(e.target.value))}
          />
        </label>
        <label>
          Backend
          <select
            value={backend}
            disabled={running}
            onChange={(e) => setBackend(e.target.value as Backend)}
          >
            <option value="auto">auto (GPU → CPU)</option>
            <option value="gpu">gpu</option>
            <option value="cpu">cpu</option>
          </select>
        </label>
        <button onClick={run} disabled={running}>
          {running ? "Running…" : "Run"}
        </button>
      </section>

      {running && (
        <div className="progress" role="progressbar" aria-valuenow={pct}>
          <div style={{ width: `${pct}%` }} />
          <span>
            {progress ? `${progress.step} / ${progress.total}` : "starting…"}
          </span>
        </div>
      )}

      {error && <p className="error">{error}</p>}

      {result && (
        <section>
          <p className="meta">
            ran on <strong>{result.backend}</strong> · {result.steps_done} steps
          </p>
          <ProbePlot series={result.probes[0] ?? []} dt={result.dt_s} />
          <SpectrumPlot series={result.probes[0] ?? []} dt={result.dt_s} />
          {result.slice && <SliceHeatmap slice={result.slice} />}
          {result.slice && (
            <Suspense fallback={<p className="meta">loading 3-D view…</p>}>
              <FieldSurface3D slice={result.slice} />
            </Suspense>
          )}
        </section>
      )}
    </main>
  );
}
