// FP32 Yee FDTD update kernels (uniform lossless vacuum, PEC box) — E.0.
//
// Six entry points, one per staggered field component, each dispatched over
// its own array extent with in-shader bounds checks. Linearization matches
// the host side (ndarray default C order): idx = (i * dim_j + j) * dim_k + k.
// The E entry points update interior cells only; outer tangential E faces
// are the PEC box and are never written, mirroring the FP64 reference in
// yee-fdtd's update.rs.

struct Params {
    nx: u32,
    ny: u32,
    nz: u32,
    _pad0: u32,
    // dt / (mu0 * mu_r) and dt / (eps0 * eps_r), precomputed in f64 host-side.
    ch: f32,
    ce: f32,
    inv_dx: f32,
    inv_dy: f32,
    inv_dz: f32,
    _pad1: f32,
    _pad2: f32,
    _pad3: f32,
}

@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read_write> ex: array<f32>;
@group(0) @binding(2) var<storage, read_write> ey: array<f32>;
@group(0) @binding(3) var<storage, read_write> ez: array<f32>;
@group(0) @binding(4) var<storage, read_write> hx: array<f32>;
@group(0) @binding(5) var<storage, read_write> hy: array<f32>;
@group(0) @binding(6) var<storage, read_write> hz: array<f32>;

// Per-component linear indices (staggered shapes; see FdtdSpec).
fn iex(i: u32, j: u32, k: u32) -> u32 { // [nx, ny+1, nz+1]
    return (i * (p.ny + 1u) + j) * (p.nz + 1u) + k;
}
fn iey(i: u32, j: u32, k: u32) -> u32 { // [nx+1, ny, nz+1]
    return (i * p.ny + j) * (p.nz + 1u) + k;
}
fn iez(i: u32, j: u32, k: u32) -> u32 { // [nx+1, ny+1, nz]
    return (i * (p.ny + 1u) + j) * p.nz + k;
}
fn ihx(i: u32, j: u32, k: u32) -> u32 { // [nx+1, ny, nz]
    return (i * p.ny + j) * p.nz + k;
}
fn ihy(i: u32, j: u32, k: u32) -> u32 { // [nx, ny+1, nz]
    return (i * (p.ny + 1u) + j) * p.nz + k;
}
fn ihz(i: u32, j: u32, k: u32) -> u32 { // [nx, ny, nz+1]
    return (i * p.ny + j) * (p.nz + 1u) + k;
}

@compute @workgroup_size(4, 4, 4)
fn update_hx(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let k = gid.z;
    if (i > p.nx || j >= p.ny || k >= p.nz) {
        return;
    }
    let dey_dz = (ey[iey(i, j, k + 1u)] - ey[iey(i, j, k)]) * p.inv_dz;
    let dez_dy = (ez[iez(i, j + 1u, k)] - ez[iez(i, j, k)]) * p.inv_dy;
    hx[ihx(i, j, k)] += p.ch * (dey_dz - dez_dy);
}

@compute @workgroup_size(4, 4, 4)
fn update_hy(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let k = gid.z;
    if (i >= p.nx || j > p.ny || k >= p.nz) {
        return;
    }
    let dez_dx = (ez[iez(i + 1u, j, k)] - ez[iez(i, j, k)]) * p.inv_dx;
    let dex_dz = (ex[iex(i, j, k + 1u)] - ex[iex(i, j, k)]) * p.inv_dz;
    hy[ihy(i, j, k)] += p.ch * (dez_dx - dex_dz);
}

@compute @workgroup_size(4, 4, 4)
fn update_hz(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let k = gid.z;
    if (i >= p.nx || j >= p.ny || k > p.nz) {
        return;
    }
    let dex_dy = (ex[iex(i, j + 1u, k)] - ex[iex(i, j, k)]) * p.inv_dy;
    let dey_dx = (ey[iey(i + 1u, j, k)] - ey[iey(i, j, k)]) * p.inv_dx;
    hz[ihz(i, j, k)] += p.ch * (dex_dy - dey_dx);
}

@compute @workgroup_size(4, 4, 4)
fn update_ex(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let k = gid.z;
    // Interior j ∈ [1, ny), k ∈ [1, nz); outer faces are PEC.
    if (i >= p.nx || j == 0u || j >= p.ny || k == 0u || k >= p.nz) {
        return;
    }
    let dhz_dy = (hz[ihz(i, j, k)] - hz[ihz(i, j - 1u, k)]) * p.inv_dy;
    let dhy_dz = (hy[ihy(i, j, k)] - hy[ihy(i, j, k - 1u)]) * p.inv_dz;
    ex[iex(i, j, k)] += p.ce * (dhz_dy - dhy_dz);
}

@compute @workgroup_size(4, 4, 4)
fn update_ey(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let k = gid.z;
    // Interior i ∈ [1, nx), k ∈ [1, nz).
    if (i == 0u || i >= p.nx || j >= p.ny || k == 0u || k >= p.nz) {
        return;
    }
    let dhx_dz = (hx[ihx(i, j, k)] - hx[ihx(i, j, k - 1u)]) * p.inv_dz;
    let dhz_dx = (hz[ihz(i, j, k)] - hz[ihz(i - 1u, j, k)]) * p.inv_dx;
    ey[iey(i, j, k)] += p.ce * (dhx_dz - dhz_dx);
}

@compute @workgroup_size(4, 4, 4)
fn update_ez(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    let k = gid.z;
    // Interior i ∈ [1, nx), j ∈ [1, ny).
    if (i == 0u || i >= p.nx || j == 0u || j >= p.ny || k >= p.nz) {
        return;
    }
    let dhy_dx = (hy[ihy(i, j, k)] - hy[ihy(i - 1u, j, k)]) * p.inv_dx;
    let dhx_dy = (hx[ihx(i, j, k)] - hx[ihx(i, j - 1u, k)]) * p.inv_dy;
    ez[iez(i, j, k)] += p.ce * (dhy_dx - dhx_dy);
}
