// chem-wasm-lens browser demo
// Run after: wasm-pack build --target web
// Serve with: python3 -m http.server 8080

// Use the loadChem convenience wrapper (handles Wasm init automatically).
// When installed from npm: import { loadChem } from '@kent-tokyo/chem-wasm-lens';
import { loadChem } from '../pkg/index.mjs';

const { MolecularSystem } = await loadChem();

// ─── helpers ────────────────────────────────────────────────────────────────

function write(id, lines) {
  document.getElementById(id).textContent = lines.join('\n');
}

// ─── 1. XYZ parsing ─────────────────────────────────────────────────────────

const WATER_XYZ = `3
water
O   0.000  0.000  0.119
H   0.000  0.757 -0.477
H   0.000 -0.757 -0.477
`;

const water = MolecularSystem.from_xyz_string(WATER_XYZ);

write('out-xyz', [
  `atom_count()       = ${water.atom_count()}`,
  `get_symbol(0)      = "${water.get_symbol(0)}"`,
  `get_symbol(1)      = "${water.get_symbol(1)}"`,
  `get_x(0), get_y(0), get_z(0) = ${water.get_x(0)}, ${water.get_y(0)}, ${water.get_z(0)}`,
  `get_symbols_json() = ${water.get_symbols_json()}`,
]);

// ─── 2. Bond detection ───────────────────────────────────────────────────────

const ETHANOL_XYZ = `9
ethanol (C2H5OH)
C  -1.228  0.000  0.000
C   0.000  0.831  0.000
O   1.180  0.083  0.000
H  -1.162  1.019  0.000
H  -1.162 -0.513  0.889
H  -2.178  0.515  0.000
H   0.000  1.456  0.886
H   0.000  1.456 -0.886
H   1.971  0.652  0.000
`;

const ethanol = MolecularSystem.from_xyz_string(ETHANOL_XYZ);

ethanol.compute_bonds();

const bondLines = [`bond_count()  = ${ethanol.bond_count()}  (C-C + C-O + 6 C/O-H)`];
for (let i = 0; i < ethanol.atom_count(); i++) {
  const sym = ethanol.get_symbol(i);
  const neighbors = Array.from(ethanol.get_bonds(i));
  if (neighbors.length > 0) {
    bondLines.push(`  atom ${i} (${sym})  ->  ${neighbors.map(j => `${j}(${ethanol.get_symbol(j)})`).join(', ')}`);
  }
}

write('out-bonds', bondLines);

// ─── 3. Bulk coordinate export ───────────────────────────────────────────────

// get_positions_flat() returns a Float32Array [x0,y0,z0, x1,y1,z1, ...]
const positions = ethanol.get_positions_flat();

write('out-bulk', [
  `get_positions_flat() → Float32Array, length = ${positions.length}  (${ethanol.atom_count()} atoms × 3)`,
  `  atom 0 (C): x=${positions[0].toFixed(3)}, y=${positions[1].toFixed(3)}, z=${positions[2].toFixed(3)}`,
  `  atom 2 (O): x=${positions[6].toFixed(3)}, y=${positions[7].toFixed(3)}, z=${positions[8].toFixed(3)}`,
  `get_symbols_json()   = ${ethanol.get_symbols_json()}`,
]);

// ─── 4. Spatial query ────────────────────────────────────────────────────────

// Build voxel grid (recommended cell size: 3-5 Å)
ethanol.build_spatial_index(3.0);

const center = 0; // C atom
const radius = 1.8; // Å

// get_atoms_within_radius uses the grid automatically after build_spatial_index()
const nearbyIndices = Array.from(ethanol.get_atoms_within_radius(center, radius));

// get_neighbors_info returns a JS array of plain objects
const nearbyInfos = ethanol.get_neighbors_info(center, radius);

const spatialLines = [
  `has_spatial_index() = ${ethanol.has_spatial_index()}`,
  `get_atoms_within_radius(${center}, ${radius}Å) = [${nearbyIndices.join(', ')}]`,
  ``,
  `get_neighbors_info(${center}, ${radius}Å):`,
];
for (const info of nearbyInfos) {
  spatialLines.push(
    `  index=${info.index}  symbol=${info.symbol}  ` +
    `x=${info.x.toFixed(3)}  y=${info.y.toFixed(3)}  z=${info.z.toFixed(3)}`
  );
}

// get_atom_info for the center atom itself
const centerInfo = ethanol.get_atom_info(center);
spatialLines.push(``, `get_atom_info(${center}) = ${JSON.stringify(centerInfo)}`);

write('out-spatial', spatialLines);

// ─── 5. PDB — residue neighbor search ───────────────────────────────────────

// Three-residue protein snippet: GLY(1), ALA(2), VAL(3) with Cα atoms
const PDB_SNIPPET = [
  'ATOM      1  N   GLY A   1       1.885  22.498   3.903  1.00  0.00           N  ',
  'ATOM      2  CA  GLY A   1       2.849  22.100   2.916  1.00  0.00           C  ',
  'ATOM      3  C   GLY A   1       4.218  22.753   3.020  1.00  0.00           C  ',
  'ATOM      4  O   GLY A   1       4.318  23.956   2.839  1.00  0.00           O  ',
  'ATOM      5  N   ALA A   2       5.235  21.946   3.275  1.00  0.00           N  ',
  'ATOM      6  CA  ALA A   2       6.584  22.449   3.412  1.00  0.00           C  ',
  'ATOM      7  C   ALA A   2       7.507  21.517   4.195  1.00  0.00           C  ',
  'ATOM      8  O   ALA A   2       7.162  20.339   4.263  1.00  0.00           O  ',
  'ATOM      9  CB  ALA A   2       7.075  22.654   1.989  1.00  0.00           C  ',
  'ATOM     10  N   VAL A   3       8.778  21.982   4.462  1.00  0.00           N  ',
  'ATOM     11  CA  VAL A   3       9.744  21.117   5.181  1.00  0.00           C  ',
  'ATOM     12  C   VAL A   3      10.886  21.938   5.782  1.00  0.00           C  ',
  'ATOM     13  O   VAL A   3      11.082  23.125   5.537  1.00  0.00           O  ',
  'ATOM     14  CB  VAL A   3       9.194  20.105   6.203  1.00  0.00           C  ',
].join('\n');

const protein = MolecularSystem.from_pdb_string(PDB_SNIPPET);
protein.build_spatial_index(4.0);

const queryAtom = 1; // Cα of GLY(1)
const searchRadius = 5.0;
const residuesNearby = protein.get_residues_within_radius(queryAtom, searchRadius);
const atomsNearby = Array.from(protein.get_atoms_within_radius(queryAtom, searchRadius));

const pdbLines = [
  `Loaded ${protein.atom_count()} atoms from PDB snippet (GLY, ALA, VAL)`,
  ``,
  `Query: atoms within ${searchRadius}Å of atom ${queryAtom}`,
  `  (${protein.get_atom_name(queryAtom)} of ${protein.get_residue_name(queryAtom)}${protein.get_residue_id(queryAtom)})`,
  ``,
  `get_atoms_within_radius   → [${atomsNearby.join(', ')}]  (${atomsNearby.length} atoms)`,
  `get_residues_within_radius → ${JSON.stringify(residuesNearby)}`,
  ``,
  `distance(queryAtom=1, atom=5) = ${protein.distance(1, 5).toFixed(3)} Å  (GLY-Cα to ALA-N)`,
];

write('out-pdb', pdbLines);

// ─── 6. Performance timing ───────────────────────────────────────────────────

function makePdb(n) {
  const lines = [];
  for (let i = 0; i < n; i++) {
    const x = ((i * 0.93) % 50).toFixed(3).padStart(8);
    const y = ((i * 1.71) % 50).toFixed(3).padStart(8);
    const z = ((i * 2.37) % 50).toFixed(3).padStart(8);
    const serial = String((i + 1) % 100000).padStart(5);
    const res    = String((Math.floor(i / 10) + 1) % 10000).padStart(4);
    lines.push(`ATOM  ${serial}  CA  ALA A${res}    ${x}${y}${z}  1.00  0.00           C  `);
  }
  return lines.join('\n');
}

function time(fn) {
  const t0 = performance.now();
  const result = fn();
  return { result, ms: (performance.now() - t0).toFixed(3) };
}

const N = 5_000;
const pdbText = makePdb(N);

// Parse once; build grid on molGrid, leave molLinear without index.
const { result: molGrid, ms: parseMs } = time(() => MolecularSystem.from_pdb_string(pdbText));
const molLinear = MolecularSystem.from_pdb_string(pdbText);

const { ms: indexMs } = time(() => molGrid.build_spatial_index(5.0));
const { result: gridResult,   ms: gridMs   } = time(() => molGrid.get_atoms_within_radius(2500, 5.0));
const { result: linearResult, ms: linearMs } = time(() => molLinear.get_atoms_within_radius(2500, 5.0));

write('out-perf', [
  `synthetic PDB: ${N.toLocaleString()} atoms, 3D-distributed (50 Å cube)`,
  ``,
  `  from_pdb_string()       ${parseMs.padStart(8)} ms`,
  `  build_spatial_index()   ${indexMs.padStart(8)} ms`,
  `  radius_query (grid)     ${gridMs.padStart(8)} ms  (${gridResult.length} atoms within 5 Å)`,
  `  radius_query (linear)   ${linearMs.padStart(8)} ms  (${linearResult.length} atoms within 5 Å)`,
  ``,
  `  grid speedup: ${(parseFloat(linearMs) / parseFloat(gridMs)).toFixed(1)}× vs linear scan`,
]);
