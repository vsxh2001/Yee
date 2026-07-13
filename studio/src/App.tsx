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

// R.5b (ADR-0199): the full-wave verify flow's response.
interface FilterVerifyResponse {
  freqs_hz: number[];
  measured_s21_db: number[];
  design_s21_db: number[];
  backend: string;
}

interface VerifyProgress {
  phase: string;
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
  const [verify, setVerify] = useState<FilterVerifyResponse | null>(null);
  const [verifyProgress, setVerifyProgress] = useState<VerifyProgress | null>(
    null,
  );
  const [verifying, setVerifying] = useState(false);
  const verifyUnlisten = useRef<UnlistenFn | null>(null);

  useEffect(() => {
    return () => {
      verifyUnlisten.current?.();
    };
  }, []);

  function specRequest(): Record<string, unknown> {
    const req: Record<string, unknown> = {
      f0_hz: f0Ghz * 1e9,
      fbw,
      order,
      eps_r: epsR,
      height_m: heightMm * 1e-3,
    };
    const ripple = Number(rippleDb);
    if (rippleDb.trim() !== "" && ripple > 0) req.ripple_db = ripple;
    return req;
  }

  async function runDesign() {
    setBusy(true);
    setError(null);
    setDesign(null);
    setVerify(null);
    try {
      setDesign(
        await invoke<FilterDesignResponse>("design_filter", {
          req: specRequest(),
        }),
      );
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }

  // R.5b: the full-wave verify loop — two engine solves streamed over
  // verify://progress, minutes of compute at full fidelity.
  async function runVerify() {
    setVerifying(true);
    setError(null);
    setVerify(null);
    setVerifyProgress(null);
    verifyUnlisten.current?.();
    verifyUnlisten.current = await listen<VerifyProgress>(
      "verify://progress",
      (e) => setVerifyProgress(e.payload),
    );
    try {
      setVerify(
        await invoke<FilterVerifyResponse>("verify_filter", {
          req: { design: specRequest() },
        }),
      );
    } catch (e) {
      setError(String(e));
    } finally {
      setVerifying(false);
      verifyUnlisten.current?.();
      verifyUnlisten.current = null;
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
        <button onClick={runDesign} disabled={busy || verifying}>
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
            <button onClick={runVerify} disabled={verifying}>
              {verifying ? "Verifying…" : "Verify (full-wave)"}
            </button>
          </section>
          {verifying && verifyProgress && (
            <p className="meta">
              verify: {verifyProgress.phase} · {verifyProgress.step} /{" "}
              {verifyProgress.total}
            </p>
          )}
          {verify && (
            <>
              <p className="meta">
                full-wave verify ran on <strong>{verify.backend}</strong>
              </p>
              <SparamPlot
                freqsHz={verify.freqs_hz}
                s21Db={verify.measured_s21_db}
                s11Db={verify.design_s21_db}
                labels={["measured |S21|", "designed |S21|"]}
              />
            </>
          )}
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
      <AntennaDesignPanel />
      <ImportPanel />
      <YieldPanel />
    </main>
  );
}

// FS.5c (ADR-0222): the yield-analysis panel — the ADR-0211 Monte-Carlo
// yield machinery (deterministic seeded MC, Wilson 95 % CI) on the
// closed-form patch-resonance testcase the surrogate-yield-001 gate
// certified. Defaults are the ADR-0211 tolerances; same seed ⇒ the same
// numbers every run.
interface YieldResponse {
  yield_frac: number;
  ci95_lo: number;
  ci95_hi: number;
  n_pass: number;
  n_samples: number;
  length_nominal_m: number;
}

export function YieldPanel() {
  const [f0Ghz, setF0Ghz] = useState(2.45);
  const [epsR, setEpsR] = useState(4.4);
  const [sigmaLMm, setSigmaLMm] = useState(0.1);
  const [sigmaEps, setSigmaEps] = useState(0.05);
  const [specMhz, setSpecMhz] = useState(40);
  const [nSamples, setNSamples] = useState(10000);
  const [seed, setSeed] = useState(20260711);
  const [resp, setResp] = useState<YieldResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  async function runYield() {
    setBusy(true);
    setError(null);
    setResp(null);
    try {
      setResp(
        await invoke<YieldResponse>("yield_estimate", {
          req: {
            f0_hz: f0Ghz * 1e9,
            eps_r: epsR,
            sigma_l_m: sigmaLMm * 1e-3,
            sigma_eps_r: sigmaEps,
            spec_halfwidth_hz: specMhz * 1e6,
            n_samples: nSamples,
            seed,
          },
        }),
      );
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }

  return (
    <section data-testid="yield-analysis">
      <h2>Yield analysis</h2>
      <p className="sub">
        patch resonance · Monte-Carlo over fab tolerances · Wilson 95 % CI
      </p>
      <section className="controls">
        <label>
          f₀ (GHz)
          <input
            type="number"
            step={0.05}
            value={f0Ghz}
            disabled={busy}
            onChange={(e) => setF0Ghz(Number(e.target.value))}
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
          σ_L (mm)
          <input
            type="number"
            step={0.01}
            value={sigmaLMm}
            disabled={busy}
            onChange={(e) => setSigmaLMm(Number(e.target.value))}
          />
        </label>
        <label>
          σ_εr
          <input
            type="number"
            step={0.01}
            value={sigmaEps}
            disabled={busy}
            onChange={(e) => setSigmaEps(Number(e.target.value))}
          />
        </label>
        <label>
          Spec ± (MHz)
          <input
            type="number"
            step={1}
            value={specMhz}
            disabled={busy}
            onChange={(e) => setSpecMhz(Number(e.target.value))}
          />
        </label>
        <label>
          Samples
          <input
            type="number"
            min={1}
            step={1000}
            value={nSamples}
            disabled={busy}
            onChange={(e) => setNSamples(Number(e.target.value))}
          />
        </label>
        <label>
          Seed
          <input
            type="number"
            min={0}
            step={1}
            value={seed}
            disabled={busy}
            onChange={(e) => setSeed(Number(e.target.value))}
          />
        </label>
        <button onClick={runYield} disabled={busy}>
          {busy ? "Sampling…" : "Estimate yield"}
        </button>
      </section>

      {error && <p className="error">{error}</p>}

      {resp && (
        <p className="meta" data-testid="yield-result">
          yield <strong>{(resp.yield_frac * 100).toFixed(2)} %</strong> · 95 %
          CI [{(resp.ci95_lo * 100).toFixed(2)},{" "}
          {(resp.ci95_hi * 100).toFixed(2)}] % · {resp.n_pass} /{" "}
          {resp.n_samples} pass · L₀ ={" "}
          {(resp.length_nominal_m * 1e3).toFixed(2)} mm
        </p>
      )}
    </section>
  );
}

// R.5c (ADR-0203): the antenna design + verify panel — Balanis patch dims
// + inset feed + Gerber artifacts (instant), then a one-solve full-wave
// |S11| verify (the A.1 observable under the A.2 open-top boundary).
interface AntennaDesignResponse {
  width_m: number;
  length_m: number;
  eps_eff: number;
  inset_m: number;
  gerber_copper: string;
  gerber_outline: string;
}

interface AntennaVerifyResponse {
  freqs_hz: number[];
  s11_db: number[];
  f_dip_hz: number;
  dip_db: number;
  backend: string;
}

export function AntennaDesignPanel() {
  const [f0Ghz, setF0Ghz] = useState(2.45);
  const [epsR, setEpsR] = useState(4.4);
  const [heightMm, setHeightMm] = useState(1.6);
  const [insetFrac, setInsetFrac] = useState(0.25);
  const [design, setDesign] = useState<AntennaDesignResponse | null>(null);
  const [verify, setVerify] = useState<AntennaVerifyResponse | null>(null);
  const [progress, setProgress] = useState<VerifyProgress | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [verifying, setVerifying] = useState(false);
  const unlisten = useRef<UnlistenFn | null>(null);

  useEffect(() => {
    return () => {
      unlisten.current?.();
    };
  }, []);

  function specRequest(): Record<string, unknown> {
    return {
      f0_hz: f0Ghz * 1e9,
      eps_r: epsR,
      height_m: heightMm * 1e-3,
      inset_frac: insetFrac,
    };
  }

  async function runDesign() {
    setBusy(true);
    setError(null);
    setDesign(null);
    setVerify(null);
    try {
      setDesign(
        await invoke<AntennaDesignResponse>("design_antenna", {
          req: specRequest(),
        }),
      );
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }

  async function runVerify() {
    setVerifying(true);
    setError(null);
    setVerify(null);
    setProgress(null);
    unlisten.current?.();
    unlisten.current = await listen<VerifyProgress>(
      "verify://progress",
      (e) => setProgress(e.payload),
    );
    try {
      setVerify(
        await invoke<AntennaVerifyResponse>("verify_antenna", {
          req: { design: specRequest() },
        }),
      );
    } catch (e) {
      setError(String(e));
    } finally {
      setVerifying(false);
      unlisten.current?.();
      unlisten.current = null;
    }
  }

  return (
    <section data-testid="antenna-design">
      <h2>Antenna design</h2>
      <p className="sub">
        inset-fed patch · Balanis dims → layout → full-wave |S11|
      </p>
      <section className="controls">
        <label>
          f₀ (GHz)
          <input
            type="number"
            step={0.05}
            value={f0Ghz}
            disabled={busy || verifying}
            onChange={(e) => setF0Ghz(Number(e.target.value))}
          />
        </label>
        <label>
          ε_r
          <input
            type="number"
            step={0.1}
            value={epsR}
            disabled={busy || verifying}
            onChange={(e) => setEpsR(Number(e.target.value))}
          />
        </label>
        <label>
          h (mm)
          <input
            type="number"
            step={0.1}
            value={heightMm}
            disabled={busy || verifying}
            onChange={(e) => setHeightMm(Number(e.target.value))}
          />
        </label>
        <label>
          Inset (·L)
          <input
            type="number"
            step={0.01}
            value={insetFrac}
            disabled={busy || verifying}
            onChange={(e) => setInsetFrac(Number(e.target.value))}
          />
        </label>
        <button onClick={runDesign} disabled={busy || verifying}>
          {busy ? "Designing…" : "Design"}
        </button>
      </section>

      {error && <p className="error">{error}</p>}

      {design && (
        <section>
          <p className="meta">
            W = {(design.width_m * 1e3).toFixed(2)} mm · L ={" "}
            {(design.length_m * 1e3).toFixed(2)} mm · ε_eff ={" "}
            {design.eps_eff.toFixed(2)} · inset ={" "}
            {(design.inset_m * 1e3).toFixed(2)} mm
          </p>
          <section className="controls">
            <button
              onClick={() => download("patch-F_Cu.gbr", design.gerber_copper)}
            >
              Export copper .gbr
            </button>
            <button
              onClick={() =>
                download("patch-Edge_Cuts.gbr", design.gerber_outline)
              }
            >
              Export outline .gbr
            </button>
            <button onClick={runVerify} disabled={verifying}>
              {verifying ? "Verifying…" : "Verify (full-wave)"}
            </button>
          </section>
          {verifying && progress && (
            <p className="meta">
              verify: {progress.step} / {progress.total}
            </p>
          )}
          {verify && (
            <>
              <p className="meta">
                dip {(verify.f_dip_hz / 1e9).toFixed(3)} GHz ·{" "}
                {verify.dip_db.toFixed(1)} dB · ran on{" "}
                <strong>{verify.backend}</strong>
              </p>
              <SparamPlot
                freqsHz={verify.freqs_hz}
                s21Db={verify.s11_db}
                labels={["measured |S11|", ""]}
              />
            </>
          )}
        </section>
      )}
    </section>
  );
}

// FS.3.1c (ADR-0209): the board-import panel — paste (or pick) a copper
// Gerber + optional Edge.Cuts outline, supply the stackup and one port
// (Gerber carries neither), and get back the parsed preview plus the
// byte-provable copper echo from `import_gerber`. The echo badge is the
// UI face of gate studio-import-e2e-001: green means what the engine
// understood re-exports byte-identically to what was pasted.
interface ImportResponse {
  trace_count: number;
  bbox_w_m: number;
  bbox_h_m: number;
  svg: string;
  gerber_copper_echo: string;
  outline_m: [number, number][] | null;
  layout_json: string;
}

// Byte-provable losslessness: the echo must equal the input exactly.
// Exported for the vitest gate.
export function echoIsLossless(input: string, echo: string): boolean {
  return input === echo;
}

export function ImportPanel() {
  const [copper, setCopper] = useState("");
  const [outline, setOutline] = useState("");
  const [epsR, setEpsR] = useState(4.4);
  const [heightMm, setHeightMm] = useState(1.6);
  const [lossTangent, setLossTangent] = useState(0.0);
  const [portXMm, setPortXMm] = useState(0.0);
  const [portYMm, setPortYMm] = useState(0.0);
  const [portWMm, setPortWMm] = useState(3.0);
  const [resp, setResp] = useState<ImportResponse | null>(null);
  const [lossless, setLossless] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  async function pickFile(setter: (text: string) => void, file: File | null) {
    if (file) setter(await file.text());
  }

  async function runImport() {
    setBusy(true);
    setError(null);
    setResp(null);
    try {
      const r = await invoke<ImportResponse>("import_gerber", {
        req: {
          copper_gerber: copper,
          outline_gerber: outline.trim() ? outline : null,
          eps_r: epsR,
          height_m: heightMm * 1e-3,
          loss_tangent: lossTangent,
          ports: [
            {
              x_m: portXMm * 1e-3,
              y_m: portYMm * 1e-3,
              width_m: portWMm * 1e-3,
              z0_ohm: 50.0,
            },
          ],
        },
      });
      setResp(r);
      setLossless(echoIsLossless(copper, r.gerber_copper_echo));
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }

  return (
    <section data-testid="board-import">
      <h2>Import board</h2>
      <p className="sub">
        Gerber copper (+ optional Edge.Cuts) → parsed layout · byte-provable
        echo
      </p>
      <section className="controls">
        <label>
          Copper .gbr
          <input
            type="file"
            accept=".gbr,.gtl"
            disabled={busy}
            onChange={(e) => pickFile(setCopper, e.target.files?.[0] ?? null)}
          />
        </label>
        <label>
          Outline .gbr (optional)
          <input
            type="file"
            accept=".gbr,.gm1"
            disabled={busy}
            onChange={(e) => pickFile(setOutline, e.target.files?.[0] ?? null)}
          />
        </label>
      </section>
      <textarea
        aria-label="Copper Gerber text"
        placeholder="…or paste the copper Gerber here"
        value={copper}
        disabled={busy}
        rows={6}
        onChange={(e) => setCopper(e.target.value)}
      />
      <section className="controls">
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
        <label>
          tan δ
          <input
            type="number"
            step={0.005}
            value={lossTangent}
            disabled={busy}
            onChange={(e) => setLossTangent(Number(e.target.value))}
          />
        </label>
        <label>
          Port x (mm)
          <input
            type="number"
            step={0.1}
            value={portXMm}
            disabled={busy}
            onChange={(e) => setPortXMm(Number(e.target.value))}
          />
        </label>
        <label>
          Port y (mm)
          <input
            type="number"
            step={0.1}
            value={portYMm}
            disabled={busy}
            onChange={(e) => setPortYMm(Number(e.target.value))}
          />
        </label>
        <label>
          Port width (mm)
          <input
            type="number"
            step={0.1}
            value={portWMm}
            disabled={busy}
            onChange={(e) => setPortWMm(Number(e.target.value))}
          />
        </label>
        <button onClick={runImport} disabled={busy || !copper.trim()}>
          {busy ? "Importing…" : "Import"}
        </button>
      </section>

      {error && <p className="error">{error}</p>}

      {resp && (
        <section>
          <p className="meta">
            {resp.trace_count} polygons · board{" "}
            {(resp.bbox_w_m * 1e3).toFixed(2)} ×{" "}
            {(resp.bbox_h_m * 1e3).toFixed(2)} mm
            {resp.outline_m ? ` · outline ${resp.outline_m.length} corners` : ""}{" "}
            ·{" "}
            <strong
              data-testid="echo-badge"
              className={lossless ? "ok" : "warn"}
            >
              {lossless
                ? "echo byte-identical (lossless)"
                : "echo differs from input"}
            </strong>
          </p>
          <figure
            className="plot"
            data-testid="import-preview"
            // Trusted content: the SVG is generated by our own layout
            // renderer from the parsed polygons, never from raw input.
            dangerouslySetInnerHTML={{ __html: resp.svg }}
          />
          <section className="controls">
            <button
              onClick={() => download("imported-layout.json", resp.layout_json)}
            >
              Export layout .json
            </button>
            <button
              onClick={() =>
                download("imported-echo-F_Cu.gbr", resp.gerber_copper_echo)
              }
            >
              Export copper echo .gbr
            </button>
          </section>
        </section>
      )}
    </section>
  );
}
