// R.5 DOM gates (ADR-0198): the S-parameter response plot renders known
// data, and the filter-design panel renders its spec form.

import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { SparamPlot } from "./views";
import { AntennaDesignPanel, FilterDesignPanel } from "./App";

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

describe("AntennaDesignPanel", () => {
  it("renders the spec form with Design action", () => {
    render(<AntennaDesignPanel />);
    const panel = screen.getByTestId("antenna-design");
    expect(panel.textContent).toContain("Antenna design");
    expect(panel.querySelectorAll("input").length).toBe(4);
  });
});

describe("SparamPlot single-trace", () => {
  it("renders one trace without the dashed overlay", () => {
    render(
      <SparamPlot
        freqsHz={[2e9, 2.5e9, 3e9]}
        s21Db={[-1, -20, -2]}
        labels={["measured |S11|", ""]}
      />,
    );
    const figs = screen.getAllByTestId("sparam-plot");
    const fig = figs[figs.length - 1];
    expect(fig.querySelectorAll("polyline").length).toBe(1);
  });
});

describe("FilterDesignPanel", () => {
  it("renders the spec form with a Design action", () => {
    render(<FilterDesignPanel />);
    const panel = screen.getByTestId("filter-design");
    expect(panel.textContent).toContain("Filter design");
    expect(panel.querySelector("button")?.textContent).toBe("Design");
    // The six spec inputs are present.
    expect(panel.querySelectorAll("input").length).toBe(6);
  });
});
