// FS.5c DOM gates (ADR-0222), gate studio-yield-dom-001: the yield panel
// renders its ADR-0211-default form, fires the `yield_estimate` command
// with correctly unit-converted arguments, and displays the returned
// yield + Wilson CI. First command-mocking test in the studio suite: the
// tauri `invoke` is replaced via vi.mock, so the gate pins the exact
// request shape the Rust command deserializes.

import { fireEvent, render, screen } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";

const { invokeMock } = vi.hoisted(() => ({ invokeMock: vi.fn() }));
vi.mock("@tauri-apps/api/core", () => ({ invoke: invokeMock }));

import { YieldPanel } from "./App";

// A canned response in the surrogate-yield-001 regime (ADR-0211 measured
// brute-force yield 0.9721 at the gate seed).
const canned = {
  yield_frac: 0.9721,
  ci95_lo: 0.9688,
  ci95_hi: 0.9752,
  n_pass: 9721,
  n_samples: 10000,
  length_nominal_m: 0.037223,
};

beforeEach(() => {
  invokeMock.mockReset();
  invokeMock.mockResolvedValue(canned);
});

describe("YieldPanel", () => {
  it("renders the tolerance form with the ADR-0211 defaults", () => {
    render(<YieldPanel />);
    const panel = screen.getByTestId("yield-analysis");
    expect(panel.textContent).toContain("Yield analysis");
    const inputs = [...panel.querySelectorAll("input[type=number]")].map(
      (i) => (i as HTMLInputElement).value,
    );
    // f0 GHz, eps_r, sigma_L mm, sigma_eps, spec MHz, samples, seed.
    expect(inputs).toEqual([
      "2.45",
      "4.4",
      "0.1",
      "0.05",
      "40",
      "10000",
      "20260711",
    ]);
  });

  it("fires yield_estimate with unit-converted args and shows the result", async () => {
    render(<YieldPanel />);
    // No RTL auto-cleanup here (house pattern, cf. import.test.tsx):
    // take the last render.
    const panels = screen.getAllByTestId("yield-analysis");
    const panel = panels[panels.length - 1];
    const button = [...panel.querySelectorAll("button")].find(
      (b) => b.textContent === "Estimate yield",
    )!;
    fireEvent.click(button);

    expect(invokeMock).toHaveBeenCalledTimes(1);
    expect(invokeMock).toHaveBeenCalledWith("yield_estimate", {
      req: {
        f0_hz: 2.45e9,
        eps_r: 4.4,
        sigma_l_m: 0.1e-3,
        sigma_eps_r: 0.05,
        spec_halfwidth_hz: 40e6,
        n_samples: 10000,
        seed: 20260711,
      },
    });

    const result = await screen.findByTestId("yield-result");
    expect(result.textContent).toContain("97.21 %");
    expect(result.textContent).toContain("[96.88, 97.52] %");
    expect(result.textContent).toContain("9721 / 10000 pass");
    expect(result.textContent).toContain("L₀ = 37.22 mm");
  });

  it("surfaces a command error in the error line", async () => {
    invokeMock.mockRejectedValueOnce("eps_r must be > 1");
    render(<YieldPanel />);
    const panels = screen.getAllByTestId("yield-analysis");
    const panel = panels[panels.length - 1];
    const button = [...panel.querySelectorAll("button")].find(
      (b) => b.textContent === "Estimate yield",
    )!;
    fireEvent.click(button);
    const err = await screen.findByText("eps_r must be > 1");
    expect(err.className).toBe("error");
  });
});
