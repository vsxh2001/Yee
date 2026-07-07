// S.3 gate (analysis): the DFT scan recovers a known sinusoid's frequency —
// a strong self-contained reference for the studio's spectrum view.

import { describe, expect, it } from "vitest";
import { dftMagnitude, divergingColor, peakIndex } from "./analysis";

describe("dftMagnitude", () => {
  it("peaks at the frequency of a pure sinusoid", () => {
    const dt = 1e-12; // 1 ps sampling → 500 GHz Nyquist
    const f0 = 50e9; // 50 GHz tone
    const n = 512;
    const series = Array.from({ length: n }, (_, k) =>
      Math.sin(2 * Math.PI * f0 * k * dt),
    );
    const spectrum = dftMagnitude(series, dt, 250);
    const peak = spectrum.freqsHz[peakIndex(spectrum)];
    // Bin spacing is 2 GHz here; the peak must land on the nearest bin.
    expect(Math.abs(peak - f0)).toBeLessThan(2e9);
  });

  it("returns the requested number of bins", () => {
    const spectrum = dftMagnitude([0, 1, 0, -1], 1e-12, 32);
    expect(spectrum.freqsHz).toHaveLength(32);
    expect(spectrum.mags).toHaveLength(32);
  });
});

describe("divergingColor", () => {
  it("maps extremes and centre correctly", () => {
    expect(divergingColor(1, 1)).toEqual([255, 0, 0, 255]); // +max → red
    expect(divergingColor(-1, 1)).toEqual([0, 0, 255, 255]); // −max → blue
    expect(divergingColor(0, 1)).toEqual([255, 255, 255, 255]); // 0 → white
  });

  it("clamps out-of-range values", () => {
    expect(divergingColor(5, 1)).toEqual([255, 0, 0, 255]);
    expect(divergingColor(-5, 1)).toEqual([0, 0, 255, 255]);
  });
});
