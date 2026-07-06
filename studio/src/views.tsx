// Visualization components (S.3, ADR-0180): frequency-spectrum plot and
// final-field slice heatmap. Both are dependency-free (inline SVG /
// canvas); three.js 3-D volumetric rendering is the S.3 follow-on.

import { useEffect, useRef } from "react";
import { dftMagnitude, divergingColor, peakIndex } from "./analysis";

export function SpectrumPlot({ series, dt }: { series: number[]; dt: number }) {
  if (series.length < 4) return null;
  const spectrum = dftMagnitude(series, dt, 160);
  const w = 640;
  const h = 200;
  const pad = 8;
  const max = Math.max(...spectrum.mags, 1e-30);
  const points = spectrum.mags
    .map((m, i) => {
      const x = pad + (i / (spectrum.mags.length - 1)) * (w - 2 * pad);
      const y = h - pad - (m / max) * (h - 2 * pad);
      return `${x.toFixed(1)},${y.toFixed(1)}`;
    })
    .join(" ");
  const peak = peakIndex(spectrum);
  const peakGhz = spectrum.freqsHz[peak] / 1e9;
  return (
    <figure className="plot" data-testid="spectrum-plot">
      <svg viewBox={`0 0 ${w} ${h}`} role="img" aria-label="Probe spectrum">
        <polyline points={points} className="trace" fill="none" />
      </svg>
      <figcaption>
        |E(f)| single-bin DFT · peak ≈ {peakGhz.toFixed(2)} GHz
      </figcaption>
    </figure>
  );
}

export interface Slice {
  ni: number;
  nj: number;
  data: number[];
}

export function SliceHeatmap({ slice }: { slice: Slice }) {
  const canvas = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const el = canvas.current;
    if (!el) return;
    const ctx = el.getContext("2d");
    if (!ctx) return;
    const { ni, nj, data } = slice;
    const max = Math.max(...data.map(Math.abs), 1e-30);
    const img = ctx.createImageData(nj, ni);
    for (let i = 0; i < ni; i++) {
      for (let j = 0; j < nj; j++) {
        const [r, g, b, a] = divergingColor(data[i * nj + j], max);
        const p = (i * nj + j) * 4;
        img.data[p] = r;
        img.data[p + 1] = g;
        img.data[p + 2] = b;
        img.data[p + 3] = a;
      }
    }
    ctx.putImageData(img, 0, 0);
  }, [slice]);

  return (
    <figure className="plot" data-testid="slice-heatmap">
      <canvas
        ref={canvas}
        width={slice.nj}
        height={slice.ni}
        className="heatmap"
        role="img"
        aria-label="Field slice heatmap"
      />
      <figcaption>
        E_z mid-plane slice · {slice.ni} × {slice.nj} (blue − / red +)
      </figcaption>
    </figure>
  );
}
