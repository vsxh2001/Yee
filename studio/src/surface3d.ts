// Pure geometry construction for the 3-D field surface (S.3b, ADR-0181).
// Kept free of three.js/DOM so vitest gates the math directly.

import type { Slice } from "./views";
import { divergingColor } from "./analysis";

export interface SurfaceGeometry {
  /** xyz triples, one vertex per slice sample (row-major over [ni, nj]). */
  positions: Float32Array;
  /** rgb triples per vertex (diverging blue-white-red by field value). */
  colors: Float32Array;
  /** Triangle indices (two triangles per grid cell). */
  indices: Uint32Array;
}

/** Build a height-map surface from a field slice: x/y span [-1, 1] on the
 *  longer axis (aspect preserved), z = value/max · `heightScale`. */
export function buildSurface(slice: Slice, heightScale = 0.35): SurfaceGeometry {
  const { ni, nj, data } = slice;
  const max = Math.max(...data.map(Math.abs), 1e-30);
  const span = Math.max(ni - 1, nj - 1, 1);

  const positions = new Float32Array(ni * nj * 3);
  const colors = new Float32Array(ni * nj * 3);
  for (let i = 0; i < ni; i++) {
    for (let j = 0; j < nj; j++) {
      const v = data[i * nj + j];
      const p = (i * nj + j) * 3;
      positions[p] = (2 * i - (ni - 1)) / span;
      positions[p + 1] = (2 * j - (nj - 1)) / span;
      positions[p + 2] = (v / max) * heightScale;
      const [r, g, b] = divergingColor(v, max);
      colors[p] = r / 255;
      colors[p + 1] = g / 255;
      colors[p + 2] = b / 255;
    }
  }

  const indices = new Uint32Array((ni - 1) * (nj - 1) * 6);
  let t = 0;
  for (let i = 0; i < ni - 1; i++) {
    for (let j = 0; j < nj - 1; j++) {
      const a = i * nj + j;
      const b = a + 1;
      const c = a + nj;
      const d = c + 1;
      indices[t++] = a;
      indices[t++] = b;
      indices[t++] = c;
      indices[t++] = b;
      indices[t++] = d;
      indices[t++] = c;
    }
  }
  return { positions, colors, indices };
}
