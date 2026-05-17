// fetch_demo.js — loads a real PDB from RCSB and demonstrates chem-wasm-lens

import init, { MolecularSystem } from '../pkg/chem_wasm_lens.js';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

await init();

// ─── Three.js scene setup ────────────────────────────────────────────────────

const canvas = document.getElementById('canvas');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(devicePixelRatio);
renderer.setSize(canvas.width, canvas.height);
renderer.setClearColor(0x1a1a2e);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(45, canvas.width / canvas.height, 0.1, 2000);
camera.position.set(0, 0, 80);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;

scene.add(new THREE.AmbientLight(0xffffff, 0.6));
const dirLight = new THREE.DirectionalLight(0xffffff, 1.0);
dirLight.position.set(10, 20, 15);
scene.add(dirLight);

let molGroup = null;

renderer.setAnimationLoop(() => {
  controls.update();
  renderer.render(scene, camera);
});

// ─── CPK colour table ─────────────────────────────────────────────────────────

const CPK = {
  C: 0x888888, N: 0x4444ff, O: 0xff4444, H: 0xeeeeee,
  S: 0xdddd00, P: 0xff8800, F: 0x44dd44, Cl: 0x44dd44,
  Br: 0xaa4400, I: 0x880088, Fe: 0xcc6600, Zn: 0x8888ff,
};
const DEFAULT_COLOR = 0xbbbbbb;

function atomColor(sym) { return CPK[sym] ?? DEFAULT_COLOR; }

// Sphere geometry cache (indexed by radius string)
const sphereCache = new Map();
function sphereGeo(r) {
  const key = r.toFixed(3);
  if (!sphereCache.has(key)) sphereCache.set(key, new THREE.SphereGeometry(r, 10, 8));
  return sphereCache.get(key);
}
const cylGeo = new THREE.CylinderGeometry(0.12, 0.12, 1, 6);

// ─── Render molecule ──────────────────────────────────────────────────────────

function renderMolecule(mol) {
  if (molGroup) { scene.remove(molGroup); molGroup.traverse(o => { o.geometry?.dispose(); o.material?.dispose(); }); }
  molGroup = new THREE.Group();

  const n = mol.atom_count();
  const pos = mol.get_positions_flat();

  // Centre of mass → shift to origin
  const com = mol.center_of_mass();
  const cx = com[0], cy = com[1], cz = com[2];

  // Atoms
  for (let i = 0; i < n; i++) {
    const sym = mol.get_symbol(i) ?? 'C';
    const r = sym === 'H' ? 0.18 : 0.35;
    const mesh = new THREE.Mesh(sphereGeo(r), new THREE.MeshPhongMaterial({ color: atomColor(sym) }));
    mesh.position.set(pos[i*3] - cx, pos[i*3+1] - cy, pos[i*3+2] - cz);
    molGroup.add(mesh);
  }

  // C-alpha backbone tube (connect consecutive CA atoms in order)
  const caIndices = [];
  for (let i = 0; i < n; i++) {
    if ((mol.get_atom_name(i) ?? '').trim() === 'CA') caIndices.push(i);
  }
  if (caIndices.length >= 2) {
    const pts = caIndices.map(i => new THREE.Vector3(pos[i*3]-cx, pos[i*3+1]-cy, pos[i*3+2]-cz));
    const curve = new THREE.CatmullRomCurve3(pts);
    const tubeGeo = new THREE.TubeGeometry(curve, caIndices.length * 3, 0.15, 6, false);
    molGroup.add(new THREE.Mesh(tubeGeo, new THREE.MeshPhongMaterial({ color: 0x55aaff, opacity: 0.7, transparent: true })));
  }

  scene.add(molGroup);

  // Auto-fit camera
  const bbox = new THREE.Box3().setFromObject(molGroup);
  const size = bbox.getSize(new THREE.Vector3()).length();
  camera.position.set(0, 0, size * 1.2);
  controls.target.set(0, 0, 0);
  controls.update();
}

// ─── Analysis outputs ─────────────────────────────────────────────────────────

function write(id, lines) {
  document.getElementById(id).textContent = lines.join('\n');
}

function countResidues(mol) {
  const seen = new Set();
  for (let i = 0; i < mol.atom_count(); i++) {
    const ch = mol.get_chain_id(i) ?? ' ';
    const rid = mol.get_residue_id(i) ?? 0;
    seen.add(`${ch}:${rid}`);
  }
  return seen.size;
}

function analyseStructure(mol, pdbId) {
  const n = mol.atom_count();
  const residues = countResidues(mol);
  const com = mol.center_of_mass();

  write('out-overview', [
    `PDB ID        : ${pdbId}`,
    `Atoms         : ${n}`,
    `Residues      : ${residues}`,
    `Center of mass: [${com[0].toFixed(2)}, ${com[1].toFixed(2)}, ${com[2].toFixed(2)}] Å`,
  ]);

  // Spatial query around atom 0
  mol.build_spatial_index(4.0);
  const radius = 8.0;
  const nearAtoms = Array.from(mol.get_atoms_within_radius(0, radius));
  const nearRes = mol.get_residues_within_radius(0, radius);
  const atom0name = (mol.get_atom_name(0) ?? '?').trim();
  const atom0res = mol.get_residue_name(0) ?? '?';
  const atom0id = mol.get_residue_id(0) ?? 0;

  write('out-spatial', [
    `Query: atoms within ${radius} Å of atom 0 (${atom0name} of ${atom0res}${atom0id})`,
    ``,
    `  get_atoms_within_radius → ${nearAtoms.length} atoms`,
    `  get_residues_within_radius → ${nearRes.length} residue(s):`,
    ...nearRes.map(r => `    ${r}`),
  ]);

  // Backbone dihedral angles: find N, CA, C, O runs for up to first 5 residues
  const dihedralLines = [
    'Backbone φ/ψ for first 5 residues (requires at least 4 consecutive backbone atoms):',
    '',
  ];
  const backbone = ['N', 'CA', 'C'];
  const byResidue = new Map();
  for (let i = 0; i < n; i++) {
    const name = (mol.get_atom_name(i) ?? '').trim();
    if (!backbone.includes(name)) continue;
    const rid = mol.get_residue_id(i) ?? 0;
    if (!byResidue.has(rid)) byResidue.set(rid, {});
    byResidue.get(rid)[name] = i;
  }
  const rids = [...byResidue.keys()].sort((a, b) => a - b).slice(0, 6);

  for (let k = 1; k < rids.length - 1; k++) {
    const prev = byResidue.get(rids[k-1]);
    const curr = byResidue.get(rids[k]);
    const next = byResidue.get(rids[k+1]);
    const resName = mol.get_residue_name(curr.CA ?? curr.N ?? 0) ?? '???';

    let phi = '—', psi = '—';
    // φ = C(i-1)–N(i)–CA(i)–C(i)
    if (prev?.C != null && curr.N != null && curr.CA != null && curr.C != null) {
      phi = mol.dihedral(prev.C, curr.N, curr.CA, curr.C).toFixed(1) + '°';
    }
    // ψ = N(i)–CA(i)–C(i)–N(i+1)
    if (curr.N != null && curr.CA != null && curr.C != null && next?.N != null) {
      psi = mol.dihedral(curr.N, curr.CA, curr.C, next.N).toFixed(1) + '°';
    }
    dihedralLines.push(`  Residue ${rids[k].toString().padStart(3)} (${resName})  φ = ${phi.padStart(8)}  ψ = ${psi.padStart(8)}`);
  }
  write('out-geometry', dihedralLines);
}

// ─── Fetch PDB from RCSB ──────────────────────────────────────────────────────

const status = document.getElementById('status');

async function loadPdb(pdbId) {
  pdbId = pdbId.toUpperCase().trim();
  if (!pdbId.match(/^[A-Z0-9]{4}$/)) { status.textContent = 'Invalid PDB ID.'; return; }

  status.textContent = `Fetching ${pdbId} from RCSB…`;
  try {
    const t0 = performance.now();
    const resp = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    const text = await resp.text();
    const fetchMs = (performance.now() - t0).toFixed(0);

    const t1 = performance.now();
    const mol = MolecularSystem.from_pdb_string(text);
    const parseMs = (performance.now() - t1).toFixed(0);

    status.textContent = `${pdbId}: ${mol.atom_count()} atoms — fetch ${fetchMs} ms, parse ${parseMs} ms`;

    analyseStructure(mol, pdbId);
    renderMolecule(mol);
  } catch (e) {
    status.textContent = `Error: ${e.message}`;
  }
}

// ─── UI wiring ────────────────────────────────────────────────────────────────

document.getElementById('btn-load').addEventListener('click', () => {
  loadPdb(document.getElementById('pdbid').value);
});
document.getElementById('pdbid').addEventListener('keydown', e => {
  if (e.key === 'Enter') loadPdb(e.target.value);
});

function quickLoad(id) {
  document.getElementById('pdbid').value = id;
  document.querySelectorAll('button.active').forEach(b => b.classList.remove('active'));
  document.getElementById(`btn-${id.toLowerCase()}`).classList.add('active');
  loadPdb(id);
}
document.getElementById('btn-1ubq').addEventListener('click', () => quickLoad('1UBQ'));
document.getElementById('btn-1crn').addEventListener('click', () => quickLoad('1CRN'));
document.getElementById('btn-4hhb').addEventListener('click', () => quickLoad('4HHB'));

// Auto-load 1UBQ on startup
loadPdb('1UBQ');
