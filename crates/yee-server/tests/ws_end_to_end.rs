//! S.1 gate: end-to-end WebSocket job round trip — bind an ephemeral port,
//! connect with a real WS client, send a `JobSpec`, and assert the streamed
//! `progress` events and the terminal `done` payload (probe series length,
//! non-trivial signal, requested field slice).

use futures_util::{SinkExt, StreamExt};
use tokio_tungstenite::tungstenite::Message;
use yee_engine::{
    BackendChoice, BoundarySpec, JobEvent, JobSpec, ProbeSpec, SliceSpec, SourceSpec,
};

fn cavity_spec() -> JobSpec {
    JobSpec {
        nx: 12,
        ny: 12,
        nz: 12,
        dx_m: 1e-3,
        n_steps: 60,
        boundary: BoundarySpec::Pec,
        sources: vec![SourceSpec::GaussianEz {
            cell: (6, 6, 6),
            t0_steps: 8.0,
            sigma_steps: 3.0,
        }],
        ports: vec![],
        aperture_ports: vec![],
        thin_wires: vec![],
        probes: vec![ProbeSpec {
            component: "ez".into(),
            cell: (8, 6, 6),
        }],
        slice: Some(SliceSpec {
            component: "ez".into(),
            k: 6,
        }),
        ntff: None,
        materials: None,
        dt_s: None,
        spacings: None,
        backend: BackendChoice::Cpu,
    }
}

#[tokio::test]
async fn websocket_job_streams_progress_and_result() {
    // Bind an ephemeral port and serve the router in the background.
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    tokio::spawn(async move {
        axum::serve(listener, yee_server::router()).await.unwrap();
    });

    // Health check.
    let health = tokio::net::TcpStream::connect(addr).await;
    assert!(health.is_ok(), "server not reachable");

    // WS round trip.
    let (mut ws, _) = tokio_tungstenite::connect_async(format!("ws://{addr}/v1/jobs"))
        .await
        .expect("WS connect failed");
    ws.send(Message::Text(
        serde_json::to_string(&cavity_spec()).unwrap().into(),
    ))
    .await
    .unwrap();

    let mut progress_events = 0usize;
    let mut done = None;
    while let Some(frame) = ws.next().await {
        match frame.unwrap() {
            Message::Text(text) => match serde_json::from_str::<JobEvent>(&text).unwrap() {
                JobEvent::Progress { step, total } => {
                    assert!(step <= total && total == 60);
                    progress_events += 1;
                }
                JobEvent::Done { result } => {
                    done = Some(result);
                }
                JobEvent::Error { message } => panic!("job failed: {message}"),
            },
            Message::Close(_) => break,
            _ => {}
        }
    }

    let result = done.expect("no done event");
    assert!(progress_events >= 10, "too few progress events");
    assert_eq!(result.steps_done, 60);
    assert_eq!(result.probes.len(), 1);
    assert_eq!(result.probes[0].len(), 60);
    assert!(result.probes[0].iter().any(|v| *v != 0.0));
    let slice = result.slice.expect("no slice in result");
    assert_eq!((slice.ni, slice.nj), (13, 13));
    assert!(slice.data.iter().any(|v| *v != 0.0));
}

#[tokio::test]
async fn invalid_spec_yields_error_event() {
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    tokio::spawn(async move {
        axum::serve(listener, yee_server::router()).await.unwrap();
    });

    let (mut ws, _) = tokio_tungstenite::connect_async(format!("ws://{addr}/v1/jobs"))
        .await
        .unwrap();
    ws.send(Message::Text("{\"not\": \"a spec\"}".into()))
        .await
        .unwrap();
    let frame = ws.next().await.expect("no reply").unwrap();
    let Message::Text(text) = frame else {
        panic!("expected text frame, got {frame:?}");
    };
    let event: JobEvent = serde_json::from_str(&text).unwrap();
    assert!(
        matches!(event, JobEvent::Error { ref message } if message.contains("invalid JobSpec")),
        "unexpected event: {event:?}"
    );
}
