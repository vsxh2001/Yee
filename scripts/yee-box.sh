#!/usr/bin/env bash
#
# yee-box — run a command inside the bounded Yee dev container.
#
# RAM and CPU are capped (via cgroups) so a heavy build or a multi-minute FDTD
# gate can NEVER OOM-kill the host: worst case the container itself is killed
# and the host is untouched. This is what lets the otherwise CI-only heavy
# gates (mom-001, fdtd-coupling-001) run locally.
#
# Usage:
#   scripts/yee-box.sh <command...>
#   YEE_BOX_DIR=worktrees/fdtd-driver scripts/yee-box.sh \
#       cargo test -p yee-voxel --release -- --ignored fdtd_coupling_001 --nocapture
#
# Env knobs:
#   YEE_BOX_MEM   container memory cap          (default 12g; host has ~29g)
#   YEE_BOX_CPUS  container CPU cap             (default 4;   host has 12)
#   YEE_BOX_DIR   directory to mount at /work   (default: repo root; set to a
#                 worktree path to build/test that worktree's checkout)
#
# Build the image once:
#   docker build -t yee-dev:1.92 -f docker/Dockerfile.dev .
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MOUNT_DIR="$(cd "${YEE_BOX_DIR:-$REPO_ROOT}" && pwd)"
MEM="${YEE_BOX_MEM:-12g}"
CPUS="${YEE_BOX_CPUS:-4}"
IMAGE="${YEE_BOX_IMAGE:-yee-dev:1.92}"

if [ "$#" -eq 0 ]; then
    echo "usage: scripts/yee-box.sh <command...>" >&2
    exit 2
fi

# A target/ volume keyed to the mounted dir keeps each checkout's build cache
# separate and off the host filesystem (so the container target never clobbers
# the host's). The cargo registry is shared (read-mostly).
TARGET_VOL="yee-target-$(echo "$MOUNT_DIR" | tr -c 'a-zA-Z0-9' '-' | sed 's/^-*//;s/-*$//')"

exec docker run --rm \
    --memory="$MEM" --memory-swap="$MEM" --cpus="$CPUS" \
    -v "$MOUNT_DIR":/work \
    -v yee-cargo-registry:/usr/local/cargo/registry \
    -v "$TARGET_VOL":/work/target \
    -w /work "$IMAGE" "$@"
