//! WebSocket exposure of the `yee-engine` job protocol (phase S.1,
//! ADR-0180).
//!
//! Routes:
//! - `GET /healthz` — liveness probe, returns `ok`.
//! - `GET /v1/jobs` (WebSocket) — the client sends one JSON
//!   [`yee_engine::JobSpec`] text frame; the server streams every
//!   [`yee_engine::JobEvent`] back as JSON text frames (`progress` …
//!   `done` / `error`) and then closes. Closing the socket early cancels
//!   the job cooperatively.
//!
//! The wire format is exactly the serde protocol the Tauri studio uses
//! in-process — one protocol, many transports (ADR-0179). Serve with
//! [`serve_blocking`] (used by `yee serve`) or embed [`router`] in an
//! existing axum app.

use axum::Router;
use axum::extract::ws::{Message, WebSocket, WebSocketUpgrade};
use axum::response::IntoResponse;
use axum::routing::get;

use yee_engine::{JobEvent, JobSpec};

/// Build the yee-server router (`/healthz`, `/v1/jobs`).
pub fn router() -> Router {
    Router::new()
        .route("/healthz", get(|| async { "ok" }))
        .route("/v1/jobs", get(ws_upgrade))
}

async fn ws_upgrade(ws: WebSocketUpgrade) -> impl IntoResponse {
    ws.on_upgrade(handle_job_socket)
}

/// One job per socket: first text frame is the `JobSpec`; events stream
/// back until `done`/`error`, then the server closes. A client disconnect
/// cancels the running job.
async fn handle_job_socket(mut socket: WebSocket) {
    // ---- receive the spec ----
    let spec: JobSpec = loop {
        match socket.recv().await {
            Some(Ok(Message::Text(text))) => match serde_json::from_str(&text) {
                Ok(spec) => break spec,
                Err(e) => {
                    let _ = send_event(
                        &mut socket,
                        &JobEvent::Error {
                            message: format!("invalid JobSpec: {e}"),
                        },
                    )
                    .await;
                    return;
                }
            },
            // Ignore pings/pongs and binary noise before the spec.
            Some(Ok(Message::Ping(_) | Message::Pong(_) | Message::Binary(_))) => continue,
            _ => return, // closed before sending a spec
        }
    };

    // ---- run the job, streaming events live ----
    // The engine emits on a std::sync channel from its worker thread; a
    // spawn_blocking bridge re-sends each event into a tokio channel so
    // the async task forwards them as they happen (live progress). If the
    // client disconnects mid-run, the job is cancelled cooperatively.
    let handle = yee_engine::submit(spec);
    let canceller = handle.canceller();
    let (tx, mut rx) = tokio::sync::mpsc::unbounded_channel();
    tokio::task::spawn_blocking(move || {
        for event in handle.events() {
            if tx.send(event).is_err() {
                break; // forwarder gone — nothing left to stream to
            }
        }
    });

    loop {
        tokio::select! {
            event = rx.recv() => match event {
                Some(event) => {
                    let terminal = matches!(event, JobEvent::Done { .. } | JobEvent::Error { .. });
                    if send_event(&mut socket, &event).await.is_err() {
                        canceller.cancel();
                        return;
                    }
                    if terminal {
                        let _ = socket.send(Message::Close(None)).await;
                        return;
                    }
                }
                None => {
                    // Engine stream ended without a terminal event.
                    let _ = send_event(&mut socket, &JobEvent::Error {
                        message: "engine stream ended unexpectedly".into(),
                    }).await;
                    return;
                }
            },
            msg = socket.recv() => match msg {
                // Tolerate client chatter (pings etc.) mid-run.
                Some(Ok(Message::Ping(_) | Message::Pong(_) | Message::Text(_) | Message::Binary(_))) => {}
                // Close or transport error: cancel the job and stop.
                _ => {
                    canceller.cancel();
                    return;
                }
            },
        }
    }
}

async fn send_event(socket: &mut WebSocket, event: &JobEvent) -> Result<(), axum::Error> {
    let json = serde_json::to_string(event).expect("JobEvent serializes");
    socket.send(Message::Text(json.into())).await
}

/// Bind `addr` and serve until the process exits, on a fresh multi-thread
/// tokio runtime. This is the `yee serve` entry point — callers without an
/// async context (the CLI) use it directly.
pub fn serve_blocking(addr: std::net::SocketAddr) -> std::io::Result<()> {
    let runtime = tokio::runtime::Runtime::new()?;
    runtime.block_on(async {
        let listener = tokio::net::TcpListener::bind(addr).await?;
        eprintln!("yee-server listening on http://{addr} (WS: /v1/jobs)");
        axum::serve(listener, router()).await
    })
}
