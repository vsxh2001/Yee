// Yee Studio — S.2 walking skeleton (ADR-0179).
//
// One screen: configure a driven vacuum FDTD job, run it on the in-process
// engine (yee-engine via the Tauri `run_job` command), watch progress
// events stream, and see the probe time series plotted. The plot is a
// dependency-free inline SVG — charting/3D libraries arrive with S.3.

import { Suspense, lazy, useEffect, useRef, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, type UnlistenFn } from "@tauri-apps/api/event";
import { SliceHeatmap, SparamPlot, SpectrumPlot, type Slice } from "./views";

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

// R.5 (ADR-0198): the closed-form filter design flow's response.
interface FilterDesignResponse {
  freqs_hz: number[];
  s11_db: number[];
  s21_db: number[];
  s2p: string;
  gerber_copper: string;
  gerber_outline: string;
  line_width_m: number;
  arm_length_m: number;
  gaps_m: number[];
  tap_offset_m: number;
}

function download(name: string, text: string) {
  const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  URL.revokeObjectURL(url);
}

// Spec-entry form + design response + export buttons (R.5 walking
// skeleton): the closed-form design flow — synthesis, dimensions, layout,
// coupling-matrix response, .s2p/Gerber artifacts — is instant; the
// full-wave verify loop stays on the engine job path above.
export function FilterDesignPanel() {
  const [f0Ghz, setF0Ghz] = useState(5.0);
  const [fbw, setFbw] = useState(0.22);
  const [order, setOrder] = useState(3);
  const [rippleDb, setRippleDb] = useState("");
  const [epsR, setEpsR] = useState(4.4);
  const [heightMm, setHeightMm] = useState(0.8);
  const [design, setDesign] = useState<FilterDesignResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  async function runDesign() {
    setBusy(true);
    setError(null);
    setDesign(null);
    try {
      const req: Record<string, unknown> = {
        f0_hz: f0Ghz * 1e9,
        fbw,
        order,
        eps_r: epsR,
        height_m: heightMm * 1e-3,
      };
      const ripple = Number(rippleDb);
      if (rippleDb.trim() !== "" && ripple > 0) req.ripple_db = ripple;
      setDesign(await invoke<FilterDesignResponse>("design_filter", { req }));
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }

  return (
    <section data-testid="filter-design">
      <h2>Filter design</h2>
      <p className="sub">
        hairpin BPF · synthesis → dimensions → layout → design response
      </p>
      <section className="controls">
        <label>
          f₀ (GHz)
          <input
            type="number"
            step={0.1}
            value={f0Ghz}
            disabled={busy}
            onChange={(e) => setF0Ghz(Number(e.target.value))}
          />
        </label>
        <label>
          FBW
          <input
            type="number"
            step={0.01}
            value={fbw}
            disabled={busy}
            onChange={(e) => setFbw(Number(e.target.value))}
          />
        </label>
        <label>
          Order
          <input
            type="number"
            min={2}
            max={9}
            value={order}
            disabled={busy}
            onChange={(e) => setOrder(Number(e.target.value))}
          />
        </label>
        <label>
          Ripple (dB, blank = Butterworth)
          <input
            type="text"
            value={rippleDb}
            disabled={busy}
            onChange={(e) => setRippleDb(e.target.value)}
          />
        </label>
        <label>
          ε_r
          <input
            type="number"
            step={0.1}
            value={epsR}
            disabled={busy}
            onChange={(e) => setEpsR(Number(e.target.value))}
          />
        </label>
        <label>
          h (mm)
          <input
            type="number"
            step={0.1}
            value={heightMm}
            disabled={busy}
            onChange={(e) => setHeightMm(Number(e.target.value))}
          />
        </label>
        <button onClick={runDesign} disabled={busy}>
          {busy ? "Designing…" : "Design"}
        </button>
      </section>

      {error && <p className="error">{error}</p>}

      {design && (
        <section>
          <p className="meta">
            w = {(design.line_width_m * 1e3).toFixed(2)} mm · arm ={" "}
            {(design.arm_length_m * 1e3).toFixed(2)} mm · gaps ={" "}
            {design.gaps_m.map((g) => (g * 1e3).toFixed(2)).join(", ")} mm ·
            tap = {(design.tap_offset_m * 1e3).toFixed(2)} mm
          </p>
          <SparamPlot
            freqsHz={design.freqs_hz}
            s11Db={design.s11_db}
            s21Db={design.s21_db}
          />
          <section className="controls">
            <button onClick={() => download("design.s2p", design.s2p)}>
              Export .s2p
            </button>
            <button
              onClick={() => download("design-F_Cu.gbr", design.gerber_copper)}
            >
              Export copper .gbr
            </button>
            <button
              onClick={() =>
                download("design-Edge_Cuts.gbr", design.gerber_outline)
              }
            >
              Export outline .gbr
            </button>
          </section>
        </section>
      )}
    </section>
  );
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

      <FilterDesignPanel />
    </main>
  );
}
