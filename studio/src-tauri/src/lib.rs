//! Yee Studio Tauri shell (S.2 walking skeleton, ADR-0179; R.5 design flow,
//! ADR-0198).
//!
//! Two commands: [`run_job`] hands a [`yee_engine::JobSpec`] to the
//! in-process engine, forwards its progress stream to the webview as
//! `job://progress` events, and resolves with the [`yee_engine::JobResult`];
//! [`design_filter`] runs the closed-form filter design flow (synthesis →
//! dimensions → layout → coupling-matrix response → `.s2p`/Gerber export
//! artifacts) — see [`design`]. The same `JobSpec`/`JobEvent` serde protocol
//! rides WebSocket in `yee-server` (S.1) — the frontend is
//! transport-agnostic by construction.

pub mod design;

use tauri::Emitter;
use yee_engine::{JobEvent, JobResult, JobSpec};

use crate::design::{FilterDesignRequest, FilterDesignResponse, design_filter_impl};

/// Run a simulation job on the engine, streaming progress events.
#[tauri::command]
async fn run_job(app: tauri::AppHandle, spec: JobSpec) -> Result<JobResult, String> {
    tauri::async_runtime::spawn_blocking(move || {
        let handle = yee_engine::submit(spec);
        for event in handle.events() {
            match event {
                JobEvent::Progress { step, total } => {
                    let _ = app.emit("job://progress", serde_json::json!({
                        "step": step,
                        "total": total,
                    }));
                }
                JobEvent::Done { result } => return Ok(result),
                JobEvent::Error { message } => return Err(message),
            }
        }
        Err("engine stream ended without a result".into())
    })
    .await
    .map_err(|e| e.to_string())?
}

/// Run the closed-form filter design flow (R.5): spec → synthesis →
/// dimensions → layout → design response + export artifacts.
#[tauri::command]
fn design_filter(req: FilterDesignRequest) -> Result<FilterDesignResponse, String> {
    design_filter_impl(&req)
}

/// Build and run the Tauri application.
pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![run_job, design_filter])
        .run(tauri::generate_context!())
        .expect("error while running yee-studio");
}
