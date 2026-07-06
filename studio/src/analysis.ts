// Pure signal-analysis helpers for the studio views (S.3, ADR-0180).
// Kept free of React/DOM so vitest gates them directly.

export interface Spectrum {
  freqsHz: number[];
  mags: number[];
}

/** Single-bin DFT magnitude scan of a real time series sampled at `dt`
 *  seconds, over `nBins` frequencies in (0, fMax]. fMax defaults to the
 *  Nyquist frequency. */
export function dftMagnitude(
  series: number[],
  dt: number,
  nBins = 128,
  fMax?: number,
): Spectrum {
  const n = series.length;
  const nyquist = 1 / (2 * dt);
  const top = fMax ?? nyquist;
  const freqsHz: number[] = [];
  const mags: number[] = [];
  for (let b = 1; b <= nBins; b++) {
    const f = (b / nBins) * top;
    const omega = 2 * Math.PI * f;
    let re = 0;
    let im = 0;
    for (let k = 0; k < n; k++) {
      const phase = omega * k * dt;
      re += series[k] * Math.cos(phase);
      im -= series[k] * Math.sin(phase);
    }
    freqsHz.push(f);
    mags.push(Math.hypot(re, im) / n);
  }
  return { freqsHz, mags };
}

/** Index of the spectrum peak. */
export function peakIndex(spectrum: Spectrum): number {
  let best = 0;
  for (let i = 1; i < spectrum.mags.length; i++) {
    if (spectrum.mags[i] > spectrum.mags[best]) best = i;
  }
  return best;
}

/** Map a value in [-max, +max] to a diverging blue–white–red RGBA pixel. */
export function divergingColor(value: number, max: number): [number, number, number, number] {
  const t = max > 0 ? Math.max(-1, Math.min(1, value / max)) : 0;
  if (t >= 0) {
    // white → red
    const s = 1 - t;
    return [255, Math.round(255 * s), Math.round(255 * s), 255];
  }
  // white → blue
  const s = 1 + t;
  return [Math.round(255 * s), Math.round(255 * s), 255, 255];
}
