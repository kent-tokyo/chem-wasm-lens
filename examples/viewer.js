import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { loadChem } from 'https://cdn.jsdelivr.net/npm/@kent-tokyo/chem-wasm-lens/dist/chem-wasm-lens.esm.js';

const { MolecularSystem } = await loadChem();

// ── CPK palette & van der Waals radii ────────────────────────────────────────

const CPK = {
  H:0xffffff, C:0xaaaaaa, N:0x4488ff, O:0xff4444,
  S:0xffdd00, P:0xff8800, F:0x44ffcc, Cl:0x44dd44,
  Br:0xaa4400, I:0x880088, Fe:0xe06633, Zn:0x888888,
};
const VDW = { H:0.31, C:0.77, N:0.75, O:0.73, S:1.02, P:1.06 };
const DFLT_COLOR  = 0xff69b4;
const DFLT_RADIUS = 0.80;
const BOND_R      = 0.07;

// ── Presets ──────────────────────────────────────────────────────────────────

const PRESETS = {
  water:   `3\nwater\nO  0.000  0.000  0.119\nH  0.000  0.757 -0.477\nH  0.000 -0.757 -0.477`,
  ethanol: `9\nethanol\nC -1.228  0.000  0.000\nC  0.000  0.831  0.000\nO  1.180  0.083  0.000\nH -1.162  1.019  0.000\nH -1.162 -0.513  0.889\nH -2.178  0.515  0.000\nH  0.000  1.456  0.886\nH  0.000  1.456 -0.886\nH  1.971  0.652  0.000`,
  caffeine:`24\ncaffeine\nN  1.2730  0.9000  0.0000\nC  1.2730 -0.4380  0.0000\nN  0.0000 -1.1260  0.0000\nC -1.2730 -0.4380  0.0000\nN -1.2730  0.9000  0.0000\nC  0.0000  1.5740  0.0000\nC  0.0000 -2.5670  0.0000\nN  2.4180  1.5870  0.0000\nC  3.5110  0.7430  0.0000\nN  3.5110 -0.5970  0.0000\nC  2.4180 -1.2100  0.0000\nO  0.0000  2.7780  0.0000\nO -2.3180 -0.9550  0.0000\nH  0.0000 -2.9300  1.0180\nH  0.0000 -2.9300 -1.0180\nH -0.8820 -2.9300  0.5090\nH  2.4180  2.6630  0.0000\nH  4.4740  1.2590  0.0000\nH  4.4740 -1.1130  0.0000\nH  2.4180 -2.2860  0.0000\nH -1.2730  1.9760  0.0000\nH -1.2730 -1.5140  0.0000\nH  0.8820 -2.9300  0.5090\nH -0.8820 -2.9300 -0.5090`,
  benzene: `12\nbenzene\nC  1.400  0.000  0.000\nC  0.700  1.212  0.000\nC -0.700  1.212  0.000\nC -1.400  0.000  0.000\nC -0.700 -1.212  0.000\nC  0.700 -1.212  0.000\nH  2.480  0.000  0.000\nH  1.240  2.148  0.000\nH -1.240  2.148  0.000\nH -2.480  0.000  0.000\nH -1.240 -2.148  0.000\nH  1.240 -2.148  0.000`,
  helix: (() => {
    const n = 12;
    const lines = [`${n}\npoly-Ala alpha-helix (CA trace)`];
    for (let i = 0; i < n; i++) {
      const a = i * 100 * Math.PI / 180;
      const x = (2.3 * Math.cos(a)).toFixed(3).padStart(7);
      const y = (2.3 * Math.sin(a)).toFixed(3).padStart(7);
      const z = (i * 1.5).toFixed(3).padStart(7);
      lines.push(`C ${x} ${y} ${z}`);
    }
    return lines.join('\n');
  })(),
};

// ── Three.js scene ────────────────────────────────────────────────────────────

const canvas   = document.getElementById('c');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.toneMapping         = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.0;

const scene  = new THREE.Scene();
scene.background = new THREE.Color(0x06060e);

const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
camera.position.set(0, 0, 18);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

scene.add(new THREE.AmbientLight(0xffffff, 0.45));
const sun  = new THREE.DirectionalLight(0xffffff, 1.1);
sun.position.set(8, 12, 10);
scene.add(sun);
const fill = new THREE.DirectionalLight(0x8899ff, 0.35);
fill.position.set(-6, -4, -8);
scene.add(fill);

// ── Bloom post-processing ─────────────────────────────────────────────────────

const composer   = new EffectComposer(renderer);
composer.addPass(new RenderPass(scene, camera));
const bloomPass = new UnrealBloomPass(
  new THREE.Vector2(window.innerWidth, window.innerHeight),
  0.6, 0.4, 0.85,
);
composer.addPass(bloomPass);
composer.addPass(new OutputPass());
let glowOn = true;

// ── Geometry / material helpers ───────────────────────────────────────────────

const geoCache = new Map();
function sphereGeo(r) {
  const k = r.toFixed(3);
  if (!geoCache.has(k)) geoCache.set(k, new THREE.SphereGeometry(r, 24, 16));
  return geoCache.get(k);
}
const cylGeo  = new THREE.CylinderGeometry(BOND_R, BOND_R, 1, 10);
const bondMat = new THREE.MeshStandardMaterial({ color: 0x777799, roughness: 0.6 });

const matCache = new Map();
function elemMat(sym, roughness = 0.35) {
  const k = `${sym}-${roughness}`;
  if (!matCache.has(k)) {
    matCache.set(k, new THREE.MeshStandardMaterial({
      color: CPK[sym] ?? DFLT_COLOR, roughness, metalness: 0.05,
    }));
  }
  return matCache.get(k);
}

// ── Scene state ───────────────────────────────────────────────────────────────

let currentMode = 'bas';
let molGroup    = null;
let atomMeshes  = [];
let lastMol = null, lastSyms = null, lastPos = null, lastN = 0;

function clearScene() {
  if (molGroup) {
    scene.remove(molGroup);
    molGroup.traverse(o => o.geometry?.dispose());
    molGroup = null;
  }
  atomMeshes = [];
}

// ── Molecule build ────────────────────────────────────────────────────────────

function buildScene(mol, syms, pos, n) {
  clearScene();
  molGroup = new THREE.Group();

  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < n; i++) { cx += pos[i*3]; cy += pos[i*3+1]; cz += pos[i*3+2]; }
  cx /= n; cy /= n; cz /= n;

  const mode = currentMode;

  if (mode !== 'stick') {
    const scale = mode === 'cpk' ? 1.0 : 0.55;
    for (let i = 0; i < n; i++) {
      const sym = syms[i];
      const r   = (VDW[sym] ?? DFLT_RADIUS) * scale;
      const m   = new THREE.Mesh(sphereGeo(r), elemMat(sym, mode === 'cpk' ? 0.25 : 0.35));
      m.position.set(pos[i*3]-cx, pos[i*3+1]-cy, pos[i*3+2]-cz);
      m.userData = { atomIndex: i, symbol: sym };
      molGroup.add(m);
      atomMeshes.push(m);
    }
  }

  if (mode === 'bas' || mode === 'stick') {
    const seen = new Set();
    for (let i = 0; i < n; i++) {
      for (const j of mol.get_bonds(i)) {
        const key = i < j ? `${i}-${j}` : `${j}-${i}`;
        if (seen.has(key)) continue;
        seen.add(key);
        const s = new THREE.Vector3(pos[i*3]-cx, pos[i*3+1]-cy, pos[i*3+2]-cz);
        const e = new THREE.Vector3(pos[j*3]-cx, pos[j*3+1]-cy, pos[j*3+2]-cz);
        const len  = s.distanceTo(e);
        const bond = new THREE.Mesh(cylGeo, bondMat);
        bond.scale.set(1, len, 1);
        bond.position.copy(s.clone().add(e).multiplyScalar(0.5));
        bond.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), e.clone().sub(s).normalize());
        molGroup.add(bond);
      }
    }
  }

  if (mode === 'ribbon') {
    const caVec = [];
    for (let i = 0; i < n; i++) {
      if (syms[i] === 'C') {
        caVec.push(new THREE.Vector3(pos[i*3]-cx, pos[i*3+1]-cy, pos[i*3+2]-cz));
      }
    }
    if (caVec.length >= 2) {
      const curve   = new THREE.CatmullRomCurve3(caVec);
      const tubeGeo = new THREE.TubeGeometry(curve, caVec.length * 12, 0.25, 10, false);
      molGroup.add(new THREE.Mesh(tubeGeo, new THREE.MeshStandardMaterial({
        color: 0x5588ff, roughness: 0.3, metalness: 0.1,
      })));
      for (const v of caVec) {
        const m = new THREE.Mesh(sphereGeo(0.22), elemMat('N'));
        m.position.copy(v);
        molGroup.add(m);
      }
    }
  }

  scene.add(molGroup);

  const box    = new THREE.Box3().setFromObject(molGroup);
  const size   = box.getSize(new THREE.Vector3()).length();
  const center = box.getCenter(new THREE.Vector3());
  controls.target.copy(center);
  camera.position.set(center.x, center.y, center.z + size * 1.6);
  controls.update();
}

// ── Loader ────────────────────────────────────────────────────────────────────

function showError(msg) {
  const el = document.getElementById('error');
  el.textContent = msg; el.style.display = 'block';
}
function clearError() { document.getElementById('error').style.display = 'none'; }

function loadMolecule(text) {
  clearError();
  text = text.trim();
  if (!text) return;

  const isPdb = /^(ATOM|HETATM|HEADER|REMARK)/m.test(text);
  const t0 = performance.now();
  let mol;
  try {
    mol = isPdb ? MolecularSystem.from_pdb_string(text) : MolecularSystem.from_xyz_string(text);
  } catch (e) { showError(String(e)); return; }
  const ms = (performance.now() - t0).toFixed(2);

  mol.compute_bonds();

  const n    = mol.atom_count();
  const pos  = mol.get_positions_flat();
  const syms = JSON.parse(mol.get_symbols_json());

  lastMol = mol; lastSyms = syms; lastPos = pos; lastN = n;
  buildScene(mol, syms, pos, n);

  document.getElementById('s-atoms').textContent = n;
  document.getElementById('s-bonds').textContent = mol.bond_count();
  document.getElementById('s-parse').textContent = ms;
}

function rebuildMode() {
  if (lastMol) buildScene(lastMol, lastSyms, lastPos, lastN);
}

// ── Resize ────────────────────────────────────────────────────────────────────

function resize() {
  const wrap = document.getElementById('canvas-wrap');
  const w = wrap.clientWidth, h = wrap.clientHeight;
  renderer.setSize(w, h, false);
  composer.setSize(w, h);
  bloomPass.resolution.set(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
window.addEventListener('resize', resize);
resize();

// ── Render loop ───────────────────────────────────────────────────────────────

renderer.setAnimationLoop(() => {
  controls.update();
  glowOn ? composer.render() : renderer.render(scene, camera);
});

// ── Atom picking ──────────────────────────────────────────────────────────────

const raycaster = new THREE.Raycaster();
const ptr       = new THREE.Vector2();
const tooltip   = document.getElementById('tooltip');

canvas.addEventListener('click', e => {
  if (!atomMeshes.length) return;
  const rect = canvas.getBoundingClientRect();
  ptr.x =  ((e.clientX - rect.left) / rect.width)  * 2 - 1;
  ptr.y = -((e.clientY - rect.top)  / rect.height) * 2 + 1;
  raycaster.setFromCamera(ptr, camera);
  const hits = raycaster.intersectObjects(atomMeshes);
  if (!hits.length) { tooltip.style.display = 'none'; return; }

  const { atomIndex, symbol } = hits[0].object.userData;
  const info = lastMol?.get_atom_info(atomIndex);
  if (!info) return;

  tooltip.textContent = '';
  const lines = [`#${atomIndex}  ${symbol}`];
  if (info.atom_name)                lines.push(`name: ${info.atom_name}`);
  if (info.residue_name)             lines.push(`res:  ${info.residue_name}${info.residue_id}`);
  if (info.chain_id?.trim())         lines.push(`chain: ${info.chain_id}`);
  lines.push(`xyz: ${info.x.toFixed(2)}, ${info.y.toFixed(2)}, ${info.z.toFixed(2)}`);

  for (const line of lines) {
    const div = document.createElement('div');
    div.textContent = line;
    tooltip.appendChild(div);
  }
  tooltip.style.display = 'block';
});

// ── UI wiring ─────────────────────────────────────────────────────────────────

document.getElementById('load-btn').addEventListener('click', () => {
  loadMolecule(document.getElementById('mol-input').value);
});

document.querySelectorAll('.preset-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    const xyz = PRESETS[btn.dataset.preset];
    document.getElementById('mol-input').value = xyz;
    loadMolecule(xyz);
  });
});

document.querySelectorAll('.mode-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.mode-btn').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    currentMode = btn.dataset.mode;
    rebuildMode();
  });
});

const glowToggle = document.getElementById('glow-toggle');
glowToggle.addEventListener('click', () => {
  glowOn = !glowOn;
  glowToggle.classList.toggle('on', glowOn);
});

// ── Initial load ──────────────────────────────────────────────────────────────

loadMolecule(PRESETS.caffeine);
document.getElementById('mol-input').value = PRESETS.caffeine;
