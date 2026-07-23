// FP32 Yee FDTD update kernels (E.1: per-cell materials + Roden-Gedney CPML
// + interior PEC masks), fused bulk+CPML, arena-buffer layout.
//
// Nine entry points: six field updates (one per staggered component, each
// dispatched over its own extent with in-shader bounds checks) and three
// PEC-mask clamps for the E components. Linearization matches the host side
// (ndarray default C order): idx = (i * dim_j + j) * dim_k + k.
//
// Arena layout (all offsets derived from nx/ny/nz in the helpers below —
// the host packs buffers in the identical order):
//   fields:   ex | ey | ez | hx | hy | hz
//   coeffs:   ca | cb | ce_cpml | ch          (four [nx+1,ny+1,nz+1] maps)
//   psi:      exy|exz|eyx|eyz|ezx|ezy | hxy|hxz|hyx|hyz|hzx|hzy
//   profiles: b | c | kappa | b_h | c_h | kappa_h   (six npml-vectors)
//   masks:    ex | ey | ez                    (u32 per element)
//
// The bulk update uses the CA/CB form everywhere: the host materializes
// ca = 1 for lossless cells, and `e = 1.0*e + cb*curl` is bit-identical to
// the plain add in IEEE 754. The CPML correction is fused into the same
// kernel — algebraically identical to the reference's separate second pass,
// because it reads the same frozen opposite-family field and each psi cell
// is touched exactly once per step. The E entry points update interior
// cells only; outer tangential E faces are never written (the PEC box is
// enforced host-side by zeroing them at upload).

struct Params {
    nx: u32,
    ny: u32,
    nz: u32,
    npml: u32,
    // Per-face CPML enable (R.3): bit 2*axis + side (side 0 = min face,
    // side 1 = max face).
    faces_mask: u32,
    has_cpml: u32,
    has_mask: u32,
    has_dispersion: u32,
    // ω·Δt for the on-GPU NTFF DFT accumulator (E.5b); 0.0 disables it.
    dft_omega_dt: f32,
}

@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read_write> fields: array<f32>;
@group(0) @binding(2) var<storage, read> coeffs: array<f32>;
@group(0) @binding(3) var<storage, read_write> psi: array<f32>;
@group(0) @binding(4) var<storage, read> profiles: array<f32>;
@group(0) @binding(5) var<storage, read> masks: array<u32>;
// Drive plumbing (E.2). drv_idx: [n_soft, n_ports, n_probes, max_steps],
// then n_soft + n_ports + n_probes field-arena offsets. drv_data:
// [0] step counter (f32; exact to 2^24), [1..] port e_z_prev state, per-port
// alpha and gamma, the per-step soft-amplitude table (max_steps x n_soft),
// the per-step port-EMF table (max_steps x n_ports), and the probe output
// region (max_steps x n_probes). All amplitudes precomputed host-side in
// f64 for the entire run, so a chunk of steps encodes with zero host
// round-trips; the `bump_step` dispatch advances the counter between steps.
@group(0) @binding(6) var<storage, read> drv_idx: array<u32>;
@group(0) @binding(7) var<storage, read_write> drv_data: array<f32>;
// Graded spacings (FS.0b.2, ADR-0214): packed INVERSE primal/dual arrays,
//   inv_sp: inv_xp[nx] | inv_yp[ny] | inv_zp[nz]
//         | inv_xd[nx+1] | inv_yd[ny+1] | inv_zd[nz+1]
// H updates multiply curl-E differences by the inverse PRIMAL cell width at
// the H sample; E updates multiply curl-H differences by the inverse DUAL
// spacing — the FS.0b.0 CPU divisor mapping (cpu.rs), verbatim. The host
// fills this with f64-computed 1/d narrowed once, so the uniform fill is
// bit-equal to the pre-FS.0b.2 scalar Params.inv_* values (compute-020).
// NOTE: 8th storage buffer per stage — the WebGPU default limit exactly;
// the next buffer must be packed into an existing arena.
@group(0) @binding(8) var<storage, read> inv_sp: array<f32>;

// ---- inverse-spacing accessors (p = primal, d = dual) ----
fn inv_xp(i: u32) -> f32 { return inv_sp[i]; }
fn inv_yp(j: u32) -> f32 { return inv_sp[p.nx + j]; }
fn inv_zp(k: u32) -> f32 { return inv_sp[p.nx + p.ny + k]; }
fn inv_xd(i: u32) -> f32 { return inv_sp[p.nx + p.ny + p.nz + i]; }
fn inv_yd(j: u32) -> f32 { return inv_sp[2u * p.nx + p.ny + p.nz + 1u + j]; }
fn inv_zd(k: u32) -> f32 { return inv_sp[2u * p.nx + 2u * p.ny + p.nz + 2u + k]; }

// ---- component lengths ----
fn len_ex() -> u32 { return p.nx * (p.ny + 1u) * (p.nz + 1u); }
fn len_ey() -> u32 { return (p.nx + 1u) * p.ny * (p.nz + 1u); }
fn len_ez() -> u32 { return (p.nx + 1u) * (p.ny + 1u) * p.nz; }
fn len_hx() -> u32 { return (p.nx + 1u) * p.ny * p.nz; }
fn len_hy() -> u32 { return p.nx * (p.ny + 1u) * p.nz; }
fn len_hz() -> u32 { return p.nx * p.ny * (p.nz + 1u); }
fn len_cell() -> u32 { return (p.nx + 1u) * (p.ny + 1u) * (p.nz + 1u); }

// ---- field arena offsets ----
fn off_ex() -> u32 { return 0u; }
fn off_ey() -> u32 { return len_ex(); }
fn off_ez() -> u32 { return off_ey() + len_ey(); }
fn off_hx() -> u32 { return off_ez() + len_ez(); }
fn off_hy() -> u32 { return off_hx() + len_hx(); }
fn off_hz() -> u32 { return off_hy() + len_hy(); }

// ---- per-component linear indices (staggered shapes) ----
fn iex(i: u32, j: u32, k: u32) -> u32 { return (i * (p.ny + 1u) + j) * (p.nz + 1u) + k; }
fn iey(i: u32, j: u32, k: u32) -> u32 { return (i * p.ny + j) * (p.nz + 1u) + k; }
fn iez(i: u32, j: u32, k: u32) -> u32 { return (i * (p.ny + 1u) + j) * p.nz + k; }
fn ihx(i: u32, j: u32, k: u32) -> u32 { return (i * p.ny + j) * p.nz + k; }
fn ihy(i: u32, j: u32, k: u32) -> u32 { return (i * (p.ny + 1u) + j) * p.nz + k; }
fn ihz(i: u32, j: u32, k: u32) -> u32 { return (i * p.ny + j) * (p.nz + 1u) + k; }
fn icell(i: u32, j: u32, k: u32) -> u32 { return (i * (p.ny + 1u) + j) * (p.nz + 1u) + k; }

// ---- field accessors ----
fn ex_at(i: u32, j: u32, k: u32) -> f32 { return fields[off_ex() + iex(i, j, k)]; }
fn ey_at(i: u32, j: u32, k: u32) -> f32 { return fields[off_ey() + iey(i, j, k)]; }
fn ez_at(i: u32, j: u32, k: u32) -> f32 { return fields[off_ez() + iez(i, j, k)]; }
fn hx_at(i: u32, j: u32, k: u32) -> f32 { return fields[off_hx() + ihx(i, j, k)]; }
fn hy_at(i: u32, j: u32, k: u32) -> f32 { return fields[off_hy() + ihy(i, j, k)]; }
fn hz_at(i: u32, j: u32, k: u32) -> f32 { return fields[off_hz() + ihz(i, j, k)]; }

// ---- coefficient accessors (arena order: ca | cb | ce_cpml | ch) ----
fn ca_at(c: u32) -> f32 { return coeffs[c]; }
fn cb_at(c: u32) -> f32 { return coeffs[len_cell() + c]; }
fn ce_cpml_at(c: u32) -> f32 { return coeffs[2u * len_cell() + c]; }
fn ch_at(c: u32) -> f32 { return coeffs[3u * len_cell() + c]; }

// ---- psi arena offsets (E-shaped x6, then H-shaped x6) ----
fn off_psi_exy() -> u32 { return 0u; }
fn off_psi_exz() -> u32 { return len_ex(); }
fn off_psi_eyx() -> u32 { return off_psi_exz() + len_ex(); }
fn off_psi_eyz() -> u32 { return off_psi_eyx() + len_ey(); }
fn off_psi_ezx() -> u32 { return off_psi_eyz() + len_ey(); }
fn off_psi_ezy() -> u32 { return off_psi_ezx() + len_ez(); }
fn off_psi_hxy() -> u32 { return off_psi_ezy() + len_ez(); }
fn off_psi_hxz() -> u32 { return off_psi_hxy() + len_hx(); }
fn off_psi_hyx() -> u32 { return off_psi_hxz() + len_hx(); }
fn off_psi_hyz() -> u32 { return off_psi_hyx() + len_hy(); }
fn off_psi_hzx() -> u32 { return off_psi_hyz() + len_hy(); }
fn off_psi_hzy() -> u32 { return off_psi_hzx() + len_hz(); }

// ---- ADE dispersion extension (E.5c) ----
// Coeff arena gains six per-cell maps after ca|cb|ce_cpml|ch:
//   ce | c0 | c1 | c2 | q | s   (unified ADE form; see dispersive.rs —
//   E' = E + ce*curl + q*(aux1' + s*aux1); the grouping preserves the
//   Lorentz/Debye near-cancellation in f32)
// Psi arena gains six aux maps after the CPML block (or after the 1-element
// dummy when CPML is off): aux1_x|aux1_y|aux1_z|aux2_x|aux2_y|aux2_z.
fn ade_ce(c: u32) -> f32 { return coeffs[4u * len_cell() + c]; }
fn ade_c0(c: u32) -> f32 { return coeffs[5u * len_cell() + c]; }
fn ade_c1(c: u32) -> f32 { return coeffs[6u * len_cell() + c]; }
fn ade_c2(c: u32) -> f32 { return coeffs[7u * len_cell() + c]; }
fn ade_q(c: u32) -> f32 { return coeffs[8u * len_cell() + c]; }
fn ade_s(c: u32) -> f32 { return coeffs[9u * len_cell() + c]; }

fn psi_cpml_len() -> u32 {
    return 2u * (len_ex() + len_ey() + len_ez()) + 2u * (len_hx() + len_hy() + len_hz());
}
fn disp_aux_base() -> u32 {
    if (p.has_cpml != 0u) {
        return psi_cpml_len();
    }
    return 1u;
}
// comp: 0 = x, 1 = y, 2 = z.
fn aux1_off(comp: u32) -> u32 { return disp_aux_base() + comp * len_cell(); }
fn aux2_off(comp: u32) -> u32 { return disp_aux_base() + (3u + comp) * len_cell(); }

// Unified ADE E update: returns the new E and rolls the aux state.
fn ade_update(comp: u32, cc: u32, e_old: f32, curl: f32) -> f32 {
    let a1 = psi[aux1_off(comp) + cc];
    let a2 = psi[aux2_off(comp) + cc];
    let a_new = ade_c0(cc) * a1 + ade_c1(cc) * a2 + ade_c2(cc) * e_old;
    psi[aux2_off(comp) + cc] = a1;
    psi[aux1_off(comp) + cc] = a_new;
    return e_old + ade_ce(cc) * curl + ade_q(cc) * (a_new + ade_s(cc) * a1);
}

// ---- CPML profile accessors (b | c | kappa | b_h | c_h | kappa_h) ----
fn prof_b(d: u32) -> f32 { return profiles[d]; }
fn prof_c(d: u32) -> f32 { return profiles[p.npml + d]; }
fn prof_kappa(d: u32) -> f32 { return profiles[2u * p.npml + d]; }
fn prof_b_h(d: u32) -> f32 { return profiles[3u * p.npml + d]; }
fn prof_c_h(d: u32) -> f32 { return profiles[4u * p.npml + d]; }
fn prof_kappa_h(d: u32) -> f32 { return profiles[5u * p.npml + d]; }

// PML profile depth for absolute index i on an axis of length n, or -1
// outside the PML / on a disabled face (mirrors CpuCpmlState's pml_depth,
// per-face since R.3).
fn pml_depth(axis: u32, i: u32, n: u32) -> i32 {
    if (i < p.npml) {
        if (((p.faces_mask >> (2u * axis)) & 1u) == 0u) {
            return -1;
        }
        return i32(p.npml - 1u - i);
    }
    if (n >= p.npml && i >= n - p.npml) {
        if (((p.faces_mask >> (2u * axis + 1u)) & 1u) == 0u) {
            return -1;
        }
        let depth = i - (n - p.npml);
        if (depth < p.npml) {
            return i32(depth);
        }
    }
    return -1;
}

// ============================= H updates =============================

@compute @workgroup_size(32, 2, 2)
fn update_hx(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    if (i > p.nx || j >= p.ny || k >= p.nz) {
        return;
    }
    let dey_dz = (ey_at(i, j, k + 1u) - ey_at(i, j, k)) * inv_zp(k);
    let dez_dy = (ez_at(i, j + 1u, k) - ez_at(i, j, k)) * inv_yp(j);
    let coeff = ch_at(icell(i, j, k));
    let idx = ihx(i, j, k);
    var h = fields[off_hx() + idx] + coeff * (dey_dz - dez_dy);
    if (p.has_cpml != 0u) {
        let dep_z = pml_depth(2u, k, p.nz);
        if (dep_z >= 0) {
            let d = u32(dep_z);
            let ps = prof_b_h(d) * psi[off_psi_hxz() + idx] + prof_c_h(d) * dey_dz;
            psi[off_psi_hxz() + idx] = ps;
            h += coeff * (ps - (1.0 - 1.0 / prof_kappa_h(d)) * dey_dz);
        }
        let dep_y = pml_depth(1u, j, p.ny);
        if (dep_y >= 0) {
            let d = u32(dep_y);
            let ps = prof_b_h(d) * psi[off_psi_hxy() + idx] + prof_c_h(d) * dez_dy;
            psi[off_psi_hxy() + idx] = ps;
            h -= coeff * (ps - (1.0 - 1.0 / prof_kappa_h(d)) * dez_dy);
        }
    }
    fields[off_hx() + idx] = h;
}

@compute @workgroup_size(32, 2, 2)
fn update_hy(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    if (i >= p.nx || j > p.ny || k >= p.nz) {
        return;
    }
    let dez_dx = (ez_at(i + 1u, j, k) - ez_at(i, j, k)) * inv_xp(i);
    let dex_dz = (ex_at(i, j, k + 1u) - ex_at(i, j, k)) * inv_zp(k);
    let coeff = ch_at(icell(i, j, k));
    let idx = ihy(i, j, k);
    var h = fields[off_hy() + idx] + coeff * (dez_dx - dex_dz);
    if (p.has_cpml != 0u) {
        let dep_x = pml_depth(0u, i, p.nx);
        if (dep_x >= 0) {
            let d = u32(dep_x);
            let ps = prof_b_h(d) * psi[off_psi_hyx() + idx] + prof_c_h(d) * dez_dx;
            psi[off_psi_hyx() + idx] = ps;
            h += coeff * (ps - (1.0 - 1.0 / prof_kappa_h(d)) * dez_dx);
        }
        let dep_z = pml_depth(2u, k, p.nz);
        if (dep_z >= 0) {
            let d = u32(dep_z);
            let ps = prof_b_h(d) * psi[off_psi_hyz() + idx] + prof_c_h(d) * dex_dz;
            psi[off_psi_hyz() + idx] = ps;
            h -= coeff * (ps - (1.0 - 1.0 / prof_kappa_h(d)) * dex_dz);
        }
    }
    fields[off_hy() + idx] = h;
}

@compute @workgroup_size(32, 2, 2)
fn update_hz(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    if (i >= p.nx || j >= p.ny || k > p.nz) {
        return;
    }
    let dex_dy = (ex_at(i, j + 1u, k) - ex_at(i, j, k)) * inv_yp(j);
    let dey_dx = (ey_at(i + 1u, j, k) - ey_at(i, j, k)) * inv_xp(i);
    let coeff = ch_at(icell(i, j, k));
    let idx = ihz(i, j, k);
    var h = fields[off_hz() + idx] + coeff * (dex_dy - dey_dx);
    if (p.has_cpml != 0u) {
        let dep_y = pml_depth(1u, j, p.ny);
        if (dep_y >= 0) {
            let d = u32(dep_y);
            let ps = prof_b_h(d) * psi[off_psi_hzy() + idx] + prof_c_h(d) * dex_dy;
            psi[off_psi_hzy() + idx] = ps;
            h += coeff * (ps - (1.0 - 1.0 / prof_kappa_h(d)) * dex_dy);
        }
        let dep_x = pml_depth(0u, i, p.nx);
        if (dep_x >= 0) {
            let d = u32(dep_x);
            let ps = prof_b_h(d) * psi[off_psi_hzx() + idx] + prof_c_h(d) * dey_dx;
            psi[off_psi_hzx() + idx] = ps;
            h -= coeff * (ps - (1.0 - 1.0 / prof_kappa_h(d)) * dey_dx);
        }
    }
    fields[off_hz() + idx] = h;
}

// ============================= E updates =============================

@compute @workgroup_size(32, 2, 2)
fn update_ex(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    // Interior j ∈ [1, ny), k ∈ [1, nz); outer faces are the PEC box.
    if (i >= p.nx || j == 0u || j >= p.ny || k == 0u || k >= p.nz) {
        return;
    }
    let dhz_dy = (hz_at(i, j, k) - hz_at(i, j - 1u, k)) * inv_yd(j);
    let dhy_dz = (hy_at(i, j, k) - hy_at(i, j, k - 1u)) * inv_zd(k);
    let cellc = icell(i, j, k);
    let idx = iex(i, j, k);
    if (p.has_dispersion != 0u) {
        fields[off_ex() + idx] =
            ade_update(0u, cellc, fields[off_ex() + idx], dhz_dy - dhy_dz);
        return;
    }
    var e = ca_at(cellc) * fields[off_ex() + idx] + cb_at(cellc) * (dhz_dy - dhy_dz);
    if (p.has_cpml != 0u) {
        let ce = ce_cpml_at(cellc);
        let dep_y = pml_depth(1u, j, p.ny + 1u);
        if (dep_y >= 0) {
            let d = u32(dep_y);
            let ps = prof_b(d) * psi[off_psi_exy() + idx] + prof_c(d) * dhz_dy;
            psi[off_psi_exy() + idx] = ps;
            e += ce * (ps - (1.0 - 1.0 / prof_kappa(d)) * dhz_dy);
        }
        let dep_z = pml_depth(2u, k, p.nz + 1u);
        if (dep_z >= 0) {
            let d = u32(dep_z);
            let ps = prof_b(d) * psi[off_psi_exz() + idx] + prof_c(d) * dhy_dz;
            psi[off_psi_exz() + idx] = ps;
            e -= ce * (ps - (1.0 - 1.0 / prof_kappa(d)) * dhy_dz);
        }
    }
    fields[off_ex() + idx] = e;
}

@compute @workgroup_size(32, 2, 2)
fn update_ey(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    // Interior i ∈ [1, nx), k ∈ [1, nz).
    if (i == 0u || i >= p.nx || j >= p.ny || k == 0u || k >= p.nz) {
        return;
    }
    let dhx_dz = (hx_at(i, j, k) - hx_at(i, j, k - 1u)) * inv_zd(k);
    let dhz_dx = (hz_at(i, j, k) - hz_at(i - 1u, j, k)) * inv_xd(i);
    let cellc = icell(i, j, k);
    let idx = iey(i, j, k);
    if (p.has_dispersion != 0u) {
        fields[off_ey() + idx] =
            ade_update(1u, cellc, fields[off_ey() + idx], dhx_dz - dhz_dx);
        return;
    }
    var e = ca_at(cellc) * fields[off_ey() + idx] + cb_at(cellc) * (dhx_dz - dhz_dx);
    if (p.has_cpml != 0u) {
        let ce = ce_cpml_at(cellc);
        let dep_z = pml_depth(2u, k, p.nz + 1u);
        if (dep_z >= 0) {
            let d = u32(dep_z);
            let ps = prof_b(d) * psi[off_psi_eyz() + idx] + prof_c(d) * dhx_dz;
            psi[off_psi_eyz() + idx] = ps;
            e += ce * (ps - (1.0 - 1.0 / prof_kappa(d)) * dhx_dz);
        }
        let dep_x = pml_depth(0u, i, p.nx + 1u);
        if (dep_x >= 0) {
            let d = u32(dep_x);
            let ps = prof_b(d) * psi[off_psi_eyx() + idx] + prof_c(d) * dhz_dx;
            psi[off_psi_eyx() + idx] = ps;
            e -= ce * (ps - (1.0 - 1.0 / prof_kappa(d)) * dhz_dx);
        }
    }
    fields[off_ey() + idx] = e;
}

@compute @workgroup_size(32, 2, 2)
fn update_ez(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    // Interior i ∈ [1, nx), j ∈ [1, ny).
    if (i == 0u || i >= p.nx || j == 0u || j >= p.ny || k >= p.nz) {
        return;
    }
    let dhy_dx = (hy_at(i, j, k) - hy_at(i - 1u, j, k)) * inv_xd(i);
    let dhx_dy = (hx_at(i, j, k) - hx_at(i, j - 1u, k)) * inv_yd(j);
    let cellc = icell(i, j, k);
    let idx = iez(i, j, k);
    if (p.has_dispersion != 0u) {
        fields[off_ez() + idx] =
            ade_update(2u, cellc, fields[off_ez() + idx], dhy_dx - dhx_dy);
        return;
    }
    var e = ca_at(cellc) * fields[off_ez() + idx] + cb_at(cellc) * (dhy_dx - dhx_dy);
    if (p.has_cpml != 0u) {
        let ce = ce_cpml_at(cellc);
        let dep_x = pml_depth(0u, i, p.nx + 1u);
        if (dep_x >= 0) {
            let d = u32(dep_x);
            let ps = prof_b(d) * psi[off_psi_ezx() + idx] + prof_c(d) * dhy_dx;
            psi[off_psi_ezx() + idx] = ps;
            e += ce * (ps - (1.0 - 1.0 / prof_kappa(d)) * dhy_dx);
        }
        let dep_y = pml_depth(1u, j, p.ny + 1u);
        if (dep_y >= 0) {
            let d = u32(dep_y);
            let ps = prof_b(d) * psi[off_psi_ezy() + idx] + prof_c(d) * dhx_dy;
            psi[off_psi_ezy() + idx] = ps;
            e -= ce * (ps - (1.0 - 1.0 / prof_kappa(d)) * dhx_dy);
        }
    }
    fields[off_ez() + idx] = e;
}

// ========================== PEC mask clamps ==========================
// Mask arena order: ex | ey | ez, u32 per element. Dispatched only when
// masks are attached; runs after the E updates so the clamp is the final
// word for the step.

@compute @workgroup_size(32, 2, 2)
fn clamp_ex(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    if (i >= p.nx || j > p.ny || k > p.nz) {
        return;
    }
    let idx = iex(i, j, k);
    if (masks[idx] != 0u) {
        fields[off_ex() + idx] = 0.0;
    }
}

@compute @workgroup_size(32, 2, 2)
fn clamp_ey(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    if (i > p.nx || j >= p.ny || k > p.nz) {
        return;
    }
    let idx = iey(i, j, k);
    if (masks[len_ex() + idx] != 0u) {
        fields[off_ey() + idx] = 0.0;
    }
}

@compute @workgroup_size(32, 2, 2)
fn clamp_ez(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    if (i > p.nx || j > p.ny || k >= p.nz) {
        return;
    }
    let idx = iez(i, j, k);
    if (masks[len_ex() + len_ey() + idx] != 0u) {
        fields[off_ez() + idx] = 0.0;
    }
}

// ============================ drive (E.2) ============================

fn d_nsoft() -> u32 { return drv_idx[0]; }
fn d_nports() -> u32 { return drv_idx[1]; }
fn d_nprobes() -> u32 { return drv_idx[2]; }
fn d_maxsteps() -> u32 { return drv_idx[3]; }
fn soft_field_off(s: u32) -> u32 { return drv_idx[4u + s]; }
fn port_field_off(q: u32) -> u32 { return drv_idx[4u + d_nsoft() + q]; }
fn probe_field_off(q: u32) -> u32 { return drv_idx[4u + d_nsoft() + d_nports() + q]; }

fn dd_step() -> u32 { return u32(drv_data[0]); }
fn dd_state(q: u32) -> u32 { return 1u + q; }
fn dd_alpha(q: u32) -> u32 { return 1u + d_nports() + q; }
fn dd_gamma(q: u32) -> u32 { return 1u + 2u * d_nports() + q; }
fn dd_amp(step: u32, s: u32) -> u32 {
    return 1u + 3u * d_nports() + step * d_nsoft() + s;
}
fn dd_vsrc(step: u32, q: u32) -> u32 {
    return 1u + 3u * d_nports() + d_maxsteps() * d_nsoft() + step * d_nports() + q;
}
fn dd_probe(step: u32, q: u32) -> u32 {
    return 1u + 3u * d_nports() + d_maxsteps() * (d_nsoft() + d_nports())
        + step * d_nprobes() + q;
}

// Soft sources: fields[off] += amp[step]; dispatched between the H and E
// half-steps, matching the reference injection point.
@compute @workgroup_size(64)
fn inject_soft(@builtin(global_invocation_id) gid: vec3<u32>) {
    let s = gid.x;
    if (s >= d_nsoft()) {
        return;
    }
    fields[soft_field_off(s)] += drv_data[dd_amp(dd_step(), s)];
}

// Resistive ports (pure-resistor lumped update, semi-implicit): applied
// after the E half-step + clamps, matching LumpedRlcPort::correct_e.
@compute @workgroup_size(64)
fn apply_ports(@builtin(global_invocation_id) gid: vec3<u32>) {
    let q = gid.x;
    if (q >= d_nports()) {
        return;
    }
    let off = port_field_off(q);
    let e1_star = fields[off];
    let e0 = drv_data[dd_state(q)];
    let v_src = drv_data[dd_vsrc(dd_step(), q)];
    let alpha = drv_data[dd_alpha(q)];
    let gamma = drv_data[dd_gamma(q)];
    let e1 = (e1_star - alpha * e0 + gamma * v_src) / (1.0 + alpha);
    fields[off] = e1;
    drv_data[dd_state(q)] = e1;
}

// Probe recording: one sample per probe per step, written to the output
// region indexed by the step counter.
@compute @workgroup_size(64)
fn record_probes(@builtin(global_invocation_id) gid: vec3<u32>) {
    let q = gid.x;
    if (q >= d_nprobes()) {
        return;
    }
    drv_data[dd_probe(dd_step(), q)] = fields[probe_field_off(q)];
}

// Advance the step counter; dispatched last in each step so every drive
// kernel in the step reads the same index.
@compute @workgroup_size(1)
fn bump_step(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x == 0u) {
        drv_data[0] = drv_data[0] + 1.0;
    }
}

// ======================= aperture ports (R.3) =======================
// Appended after the E.2 regions in both buffers, so every accessor above
// keeps its offset. drv_idx: [n_ap, n_cells x n_ap, cells_start x n_ap,
// flat E_z field-offsets...]; drv_data: [v_prev x n_ap, vcoef x n_ap,
// g x n_ap, back x n_ap, v_src (max_steps x n_ap)].

fn ap_base() -> u32 { return 4u + d_nsoft() + d_nports() + d_nprobes(); }
fn d_nap() -> u32 { return drv_idx[ap_base()]; }
fn ap_ncells(q: u32) -> u32 { return drv_idx[ap_base() + 1u + q]; }
fn ap_cell(q: u32, c: u32) -> u32 {
    let start = drv_idx[ap_base() + 1u + d_nap() + q];
    return drv_idx[ap_base() + 1u + 2u * d_nap() + start + c];
}
fn dd_ap_base() -> u32 {
    return 1u + 3u * d_nports()
        + d_maxsteps() * (d_nsoft() + d_nports() + d_nprobes());
}
fn dd_ap_vprev(q: u32) -> u32 { return dd_ap_base() + q; }
fn dd_ap_vcoef(q: u32) -> u32 { return dd_ap_base() + d_nap() + q; }
fn dd_ap_g(q: u32) -> u32 { return dd_ap_base() + 2u * d_nap() + q; }
fn dd_ap_back(q: u32) -> u32 { return dd_ap_base() + 3u * d_nap() + q; }
fn dd_ap_vsrc(step: u32, q: u32) -> u32 {
    return dd_ap_base() + 4u * d_nap() + step * d_nap() + q;
}

// Aperture ports (LumpedRlcPort::correct_e_aperture, pure-R arm): modal
// V*_T = vcoef * sum(E_z) over the aperture cells (vcoef = dz/n_col),
// semi-implicit midpoint voltage against the cached V_prev, aggregate
// branch current I = (V_mid - V_src) * g with g = 1/(R + beta) precomputed
// host-side (0 for an open port), sheet back-action back*I on every cell,
// then V_prev re-summed from the corrected field — matching the CPU's
// explicit resum. One invocation per port: ports own disjoint cell sets,
// and the serial per-port loop keeps the update consistent without
// cross-invocation synchronization. Dispatched after apply_ports, before
// record_probes (the reference order).
@compute @workgroup_size(64)
fn apply_aperture_ports(@builtin(global_invocation_id) gid: vec3<u32>) {
    let q = gid.x;
    if (q >= d_nap()) {
        return;
    }
    let ncells = ap_ncells(q);
    let vcoef = drv_data[dd_ap_vcoef(q)];
    var v_sum = 0.0;
    for (var c = 0u; c < ncells; c++) {
        v_sum += fields[ap_cell(q, c)];
    }
    let v_star = v_sum * vcoef;
    let v_mid = 0.5 * (v_star + drv_data[dd_ap_vprev(q)]);
    let v_src = drv_data[dd_ap_vsrc(dd_step(), q)];
    let i_branch = (v_mid - v_src) * drv_data[dd_ap_g(q)];
    let back = drv_data[dd_ap_back(q)] * i_branch;
    var v_post = 0.0;
    for (var c = 0u; c < ncells; c++) {
        let idx = ap_cell(q, c);
        let e = fields[idx] - back;
        fields[idx] = e;
        v_post += e;
    }
    drv_data[dd_ap_vprev(q)] = v_post * vcoef;
}


// ======================= NTFF DFT accumulation (E.5b) =======================
// Full-field running DFT at f_probe: after each completed step (fields at
// t = (n+1)·Δt, matching NtffState::sample's timing) accumulate
//   acc_re[i] += F[i]·cos(ω·t),  acc_im[i] −= F[i]·sin(ω·t).
// The phasor pair lives at the end of the psi arena, after the CPML and
// dispersion blocks; the host reads it back once and feeds the reference
// NtffState through two synthetic samples (the accumulation is linear).

fn len_fields() -> u32 {
    return off_hz() + len_hz();
}
fn dft_base() -> u32 {
    if (p.has_dispersion != 0u) {
        return disp_aux_base() + 6u * len_cell();
    }
    return disp_aux_base();
}

@compute @workgroup_size(64)
fn accumulate_dft(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= len_fields()) {
        return;
    }
    let phase = p.dft_omega_dt * (drv_data[0] + 1.0);
    let f = fields[i];
    psi[dft_base() + i] += f * cos(phase);
    psi[dft_base() + len_fields() + i] -= f * sin(phase);
}
