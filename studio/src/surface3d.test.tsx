// S.3b gates: the surface geometry is checked against hand-computable
// values, and the 3-D component's WebGL fallback DOM-renders under jsdom
// (which has no WebGL — exactly the fallback path).

import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { FieldSurface3D } from "./FieldSurface3D";
import { buildSurface } from "./surface3d";

describe("buildSurface", () => {
  const slice = {
    ni: 2,
    nj: 3,
    data: [0, 0.5, -1, 1, 0, 0.5],
  };

  it("produces one vertex per sample and two triangles per cell", () => {
    const g = buildSurface(slice, 0.5);
    expect(g.positions).toHaveLength(2 * 3 * 3);
    expect(g.colors).toHaveLength(2 * 3 * 3);
    expect(g.indices).toHaveLength(1 * 2 * 6);
  });

  it("maps extremes to ±heightScale and centres the grid", () => {
    const g = buildSurface(slice, 0.5);
    // max |value| = 1 → z at data[3] (=+1, vertex 3) is +0.5, data[2] is −0.5.
    expect(g.positions[3 * 3 + 2]).toBeCloseTo(0.5);
    expect(g.positions[2 * 3 + 2]).toBeCloseTo(-0.5);
    // x/y are centred: first vertex mirrors the last.
    expect(g.positions[0]).toBeCloseTo(-g.positions[(2 * 3 - 1) * 3]);
    expect(g.positions[1]).toBeCloseTo(-g.positions[(2 * 3 - 1) * 3 + 1]);
  });

  it("indices reference valid vertices", () => {
    const g = buildSurface(slice);
    for (const idx of g.indices) {
      expect(idx).toBeLessThan(2 * 3);
    }
  });
});

describe("FieldSurface3D", () => {
  it("renders the WebGL fallback under jsdom", () => {
    const slice = { ni: 4, nj: 4, data: Array(16).fill(0.1) };
    render(<FieldSurface3D slice={slice} />);
    const fig = screen.getByTestId("field-surface-3d");
    expect(fig.textContent).toContain("3-D surface view unavailable");
  });
});
