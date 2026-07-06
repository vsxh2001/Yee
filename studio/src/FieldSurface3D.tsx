// 3-D field surface view (S.3b, ADR-0181): the slice rendered as a
// height-mapped, vertex-colored mesh with orbit controls. Falls back to a
// text notice when WebGL is unavailable (e.g. jsdom in the DOM gates).

import { useEffect, useRef, useState } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { buildSurface } from "./surface3d";
import type { Slice } from "./views";

function webglAvailable(): boolean {
  try {
    const canvas = document.createElement("canvas");
    return canvas.getContext("webgl2") != null || canvas.getContext("webgl") != null;
  } catch {
    return false;
  }
}

export function FieldSurface3D({ slice }: { slice: Slice }) {
  const mount = useRef<HTMLDivElement>(null);
  const [supported] = useState(webglAvailable);

  useEffect(() => {
    const el = mount.current;
    if (!el || !supported) return;

    const width = el.clientWidth || 640;
    const height = 320;
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.01, 100);
    camera.position.set(1.4, -1.4, 1.1);
    camera.up.set(0, 0, 1);

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(window.devicePixelRatio);
    el.appendChild(renderer.domElement);

    const { positions, colors, indices } = buildSurface(slice);
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
    geometry.setIndex(new THREE.BufferAttribute(indices, 1));
    geometry.computeVertexNormals();
    const mesh = new THREE.Mesh(
      geometry,
      new THREE.MeshLambertMaterial({
        vertexColors: true,
        side: THREE.DoubleSide,
      }),
    );
    scene.add(mesh);
    scene.add(new THREE.AmbientLight(0xffffff, 0.7));
    const sun = new THREE.DirectionalLight(0xffffff, 1.2);
    sun.position.set(2, -1, 3);
    scene.add(sun);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;

    let alive = true;
    const tick = () => {
      if (!alive) return;
      controls.update();
      renderer.render(scene, camera);
      requestAnimationFrame(tick);
    };
    tick();

    return () => {
      alive = false;
      controls.dispose();
      geometry.dispose();
      renderer.dispose();
      el.removeChild(renderer.domElement);
    };
  }, [slice, supported]);

  return (
    <figure className="plot" data-testid="field-surface-3d">
      {supported ? (
        <div ref={mount} className="surface3d" />
      ) : (
        <p className="surface3d-fallback">
          3-D surface view unavailable (no WebGL context) — see the heatmap above.
        </p>
      )}
      <figcaption>E_z mid-plane as a 3-D surface · drag to orbit</figcaption>
    </figure>
  );
}
