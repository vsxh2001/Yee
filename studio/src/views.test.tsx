// S.3 gate (DOM smoke): the spectrum plot and slice heatmap render from
// fixture data — the "DOM-level smoke gates" the roadmap's S.3 row calls
// for, without needing a webview.

import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { SliceHeatmap, SpectrumPlot } from "./views";

describe("SpectrumPlot", () => {
  it("renders an SVG with a peak caption", () => {
    const dt = 1e-12;
    const series = Array.from({ length: 128 }, (_, k) =>
      Math.sin(2 * Math.PI * 50e9 * k * dt),
    );
    render(<SpectrumPlot series={series} dt={dt} />);
    const fig = screen.getByTestId("spectrum-plot");
    expect(fig.querySelector("svg polyline")).toBeTruthy();
    expect(fig.textContent).toContain("peak");
  });

  it("renders nothing for a too-short series", () => {
    const { container } = render(<SpectrumPlot series={[1, 2]} dt={1e-12} />);
    expect(container.querySelector("figure")).toBeNull();
  });
});

describe("SliceHeatmap", () => {
  it("renders a canvas sized to the slice", () => {
    const slice = {
      ni: 5,
      nj: 7,
      data: Array.from({ length: 35 }, (_, i) => Math.sin(i)),
    };
    render(<SliceHeatmap slice={slice} />);
    const fig = screen.getByTestId("slice-heatmap");
    const canvas = fig.querySelector("canvas");
    expect(canvas).toBeTruthy();
    expect(canvas?.getAttribute("width")).toBe("7");
    expect(canvas?.getAttribute("height")).toBe("5");
    expect(fig.textContent).toContain("5 × 7");
  });
});
