// R.5 DOM gates (ADR-0198): the S-parameter response plot renders known
// data, and the filter-design panel renders its spec form.

import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { SparamPlot } from "./views";
import { FilterDesignPanel } from "./App";

describe("SparamPlot", () => {
  it("renders two traces and the band caption for known data", () => {
    const freqsHz = [4.0e9, 4.5e9, 5.0e9, 5.5e9, 6.0e9];
    const s21Db = [-30, -3, 0, -3, -30];
    const s11Db = [-0.5, -6, -25, -6, -0.5];
    render(<SparamPlot freqsHz={freqsHz} s11Db={s11Db} s21Db={s21Db} />);
    const fig = screen.getByTestId("sparam-plot");
    const polylines = fig.querySelectorAll("polyline");
    expect(polylines.length).toBe(2);
    // Every polyline vertex stays inside the viewBox even for values below
    // the -60 dB floor.
    for (const pl of polylines) {
      const ys = (pl.getAttribute("points") ?? "")
        .split(" ")
        .map((p) => Number(p.split(",")[1]));
      for (const y of ys) {
        expect(y).toBeGreaterThanOrEqual(0);
        expect(y).toBeLessThanOrEqual(240);
      }
    }
    expect(fig.textContent).toContain("4.00–6.00 GHz");
  });

  it("renders nothing for a degenerate raster", () => {
    const { container } = render(
      <SparamPlot freqsHz={[5e9]} s11Db={[0]} s21Db={[0]} />,
    );
    expect(container.firstChild).toBeNull();
  });
});

describe("FilterDesignPanel", () => {
  it("renders the spec form with a Design action", () => {
    render(<FilterDesignPanel />);
    const panel = screen.getByTestId("filter-design");
    expect(panel.textContent).toContain("Filter design");
    expect(screen.getByText("Design")).toBeTruthy();
    // The six spec inputs are present.
    expect(panel.querySelectorAll("input").length).toBe(6);
  });
});
