//! Yee Studio Tauri shell (S.2 walking skeleton, ADR-0179).
//!
//! One command: [`run_job`] hands a [`yee_engine::JobSpec`] to the
//! in-process engine, forwards its progress stream to the webview as
//! `job://progress` events, and resolves with the [`yee_engine::JobResult`].
//! The same `JobSpec`/`JobEvent` serde protocol will ride WebSocket in
//! `yee-server` (S.1) — the frontend is transport-agnostic by construction.

use tauri::Emitter;
use yee_engine::{JobEvent, JobResult, JobSpec};

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

/// Build and run the Tauri application.
pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![run_job])
        .run(tauri::generate_context!())
        .expect("error while running yee-studio");
}
