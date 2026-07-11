// FS.3.1c DOM gates (ADR-0209): the board-import panel renders its form
// (file pickers, paste area, stackup + port fields, gated Import button),
// and the lossless-echo predicate is strict byte equality — the UI face
// of gate studio-import-e2e-001.

import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { ImportPanel, echoIsLossless } from "./App";

describe("echoIsLossless", () => {
  it("is strict byte equality — no normalization", () => {
    const g = "%FSLAX46Y46*%\n%MOMM*%\nG36*\nX0Y0D02*\nG37*\nM02*\n";
    expect(echoIsLossless(g, g)).toBe(true);
    // A single trailing newline difference must NOT count as lossless:
    // "byte-provable" means bytes, not semantics.
    expect(echoIsLossless(g, g + "\n")).toBe(false);
    expect(echoIsLossless(g, g.replace("MOMM", "MOIN"))).toBe(false);
  });
});

describe("ImportPanel", () => {
  it("renders pickers, paste area, stackup + port fields", () => {
    render(<ImportPanel />);
    const panel = screen.getByTestId("board-import");
    expect(panel.textContent).toContain("Import board");
    // 2 file pickers + 6 numeric fields (eps_r, h, tan d, port x/y/w).
    expect(panel.querySelectorAll("input[type=file]").length).toBe(2);
    expect(panel.querySelectorAll("input[type=number]").length).toBe(6);
    expect(panel.querySelector("textarea")).not.toBeNull();
  });

  it("gates the Import action on non-empty copper input", () => {
    render(<ImportPanel />);
    // No RTL auto-cleanup here (house pattern, cf. sparam.test.tsx):
    // take the last render.
    const panels = screen.getAllByTestId("board-import");
    const panel = panels[panels.length - 1];
    const button = [...panel.querySelectorAll("button")].find(
      (b) => b.textContent === "Import",
    )!;
    expect(button.disabled).toBe(true);
    fireEvent.change(panel.querySelector("textarea")!, {
      target: { value: "%MOMM*%\nM02*\n" },
    });
    expect(button.disabled).toBe(false);
  });
});
