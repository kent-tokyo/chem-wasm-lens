# chem-wasm-lens

[![CI](https://github.com/kent-tokyo/chem-wasm-lens/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chem-wasm-lens/actions)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chem-wasm-lens)](https://www.npmjs.com/package/@kent-tokyo/chem-wasm-lens)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

An ultra-lightweight molecular analysis kernel in **pure Rust**, compiled to **WebAssembly**. Parse PDB / SDF / XYZ / mmCIF files, run spatial queries, detect bonds, compute geometry, and perform fingerprint similarity — all inside a browser Web Worker, with zero C/C++ dependencies.

![chem-wasm-lens 3D viewer — ubiquitin loaded from RCSB](images/screenshot.png)

---

## Install

```bash
npm install @kent-tokyo/chem-wasm-lens
```

```js
// One-liner init — handles Wasm loading automatically
import { loadChem } from '@kent-tokyo/chem-wasm-lens';
const { MolecularSystem } = await loadChem();
```

Or, if you prefer explicit control:

```js
import init, { MolecularSystem } from '@kent-tokyo/chem-wasm-lens/raw';
await init();
```

### CDN (no install)

Drop into any page with a single `<script>` — no npm, no bundler:

```html
<script type="module">
import { loadChem } from
  'https://cdn.jsdelivr.net/npm/@kent-tokyo/chem-wasm-lens/dist/chem-wasm-lens.esm.js';

const { MolecularSystem } = await loadChem();
const mol = MolecularSystem.from_smiles('c1ccccc1');
mol.compute_2d_coords();
const svgEl = new DOMParser()
  .parseFromString(mol.to_svg_string(400, 300), 'image/svg+xml').documentElement;
document.body.appendChild(svgEl);
</script>
```

> Wasm (~411 KB) is inlined in the bundle; no separate `.wasm` file fetch needed.
> Gzip transfer size: ~200 KB.

---

## Quick start (10 lines)

```js
import { loadChem } from '@kent-tokyo/chem-wasm-lens';
const { MolecularSystem } = await loadChem();

// Load a real protein from RCSB
const pdb = await fetch('https://files.rcsb.org/download/1UBQ.pdb').then(r => r.text());
const mol = MolecularSystem.from_pdb_string(pdb);

mol.build_spatial_index(5.0);                       // O(1) avg neighbor lookup
const neighbors = mol.get_neighbors_info(0, 8.0);   // atoms within 8 Å of atom 0
console.log(`${mol.atom_count()} atoms, ${neighbors.length} neighbors near N-terminus`);

// Backbone dihedral angle
console.log(`phi = ${mol.dihedral(0, 1, 2, 3).toFixed(1)}°`);
```

---

## Why chem-wasm-lens?

| | chem-wasm-lens | RDKit.js | Pure JS |
|---|:---:|:---:|:---:|
| Runs in browser | Yes | Yes | Yes |
| Non-blocking (Web Worker) | Yes | Partial | No |
| Bundle size | ~200 KB | ~10 MB+ | varies |
| C/C++ deps | **None** | Yes | None |
| Parse 10k-atom PDB | 1.5 ms | ~20 ms | ~200 ms |
| Radius query (10k atoms) | 1.5 µs | — | ~6 µs |
| XYZ / PDB / SDF | Yes | Partial | varies |
| Spatial index | Yes (voxel grid) | Yes | No |
| Geometry (angle/dihedral) | Yes | Yes | No |
| Kabsch superposition | Yes | Yes | No |
| Scaffold analysis | Yes | No | No |
| SMILES + stereochemistry | Yes | Yes | No |
| SVG 2D rendering | Yes | No | No |
| SMARTS substructure | Yes | Yes | No |
| Fingerprint (ECFP4) + Tanimoto | Yes | Yes | No |
| Lipinski / ADMET descriptors | Yes | Yes | No |

Best for: **structural biology**, **cheminformatics**, **educational tools** — where you need real molecular analysis in the browser without a 10 MB dependency.

---

## Features

### Parsers
| Method | Format | Notes |
|--------|--------|-------|
| `from_xyz_string(s)` | XYZ | Element symbols + 3D coordinates |
| `from_pdb_string(s)` | PDB | ATOM/HETATM, CONECT, chain/residue metadata |
| `from_sdf_string(s)` | SDF V2000 | Explicit bonds from bond block |
| `from_smiles(s)` | SMILES | Organic subset; topology + implicit H |

### Spatial queries
| Method | Returns | Notes |
|--------|---------|-------|
| `build_spatial_index(cell_size)` | `void` | Voxel grid; enables O(1) avg queries |
| `get_atoms_within_radius(center, r)` | `Uint32Array` | Auto-uses grid if built |
| `get_residues_within_radius(center, r)` | `string[]` | `"CHAIN:RESNAME:RESID"` format |
| `get_neighbors_info(center, r)` | `AtomInfo[]` | Full object array |
| `distance(i, j)` | `number` | Euclidean distance (Å) |

### Geometry
| Method | Returns | Notes |
|--------|---------|-------|
| `angle(i, j, k)` | `number` | Bond angle at j (degrees) |
| `dihedral(i, j, k, l)` | `number` | Torsion angle −180..180° |
| `center_of_mass()` | `Float32Array` | Mass-weighted [x, y, z] |
| `rmsd(other)` | `number` | Without superposition |
| `superpose(reference)` | `number` | Kabsch alignment; updates coords in-place |
| `rmsd_aligned(reference)` | `number` | Kabsch RMSD without mutating coords |

### Bond topology
| Method | Returns | Notes |
|--------|---------|-------|
| `compute_bonds()` | `void` | Cordero 2008 covalent radii (18 elements) |
| `get_bonds(atom)` | `Uint32Array` | Neighbor indices |
| `bond_count()` | `number` | Unique bond count |
| `get_fragments()` | `Array<Uint32Array>` | Connected components |
| `has_bonds_computed()` | `boolean` | |
| `murcko_scaffold_indices()` | `Uint32Array` | Ring systems + linker atoms (Murcko framework) |
| `ring_system_count()` | `number` | Connected ring system count |

### 2D structure rendering
| Method | Returns | Notes |
|--------|---------|-------|
| `compute_2d_coords()` | `void` | Ring-aware 2D layout (regular polygon rings + zig-zag chains) |
| `to_svg_string(width, height)` | `string` | SVG markup; H hidden, CPK colours, aromatic inner dashes, wedge/dash stereo bonds |
| `to_svg_string_highlighted(width, height, atoms)` | `string` | SVG with yellow halos on specified atoms (`Uint32Array`) |
| `is_aromatic(index)` | `boolean` | True for SMILES lowercase atoms (c/n/o/s); false for XYZ/PDB/SDF |
| `molecular_formula()` | `string` | Hill-order formula, e.g. `"C9H8O4"` for aspirin |
| `molecular_weight()` | `number` | Sum of atomic masses in Daltons |
| `stereo_center_count()` | `number` | Number of tetrahedral stereo centers parsed from `@`/`@@` |
| `is_stereo_center(index)` | `boolean` | True if atom was annotated with `@` or `@@` in SMILES |
| `tpsa()` | `number` | Topological Polar Surface Area in Å² (Ertl 2000) |
| `h_bond_donors()` | `number` | N/O atoms with at least one H neighbour |
| `h_bond_acceptors()` | `number` | Total N/O atom count |
| `rotatable_bond_count()` | `number` | Non-ring single bonds between heavy atoms (requires `compute_rings()` first) |
| `logp()` | `number` | Estimated logP (Wildman-Crippen atom contributions) |
| `fingerprint_ecfp4()` | `Uint8Array` | 2048-bit Morgan/ECFP4 fingerprint (256 bytes) |
| `tanimoto_similarity(other)` | `number` | Tanimoto coefficient 0.0–1.0 versus another `MolecularSystem` |

### SMARTS substructure search
| Method | Returns | Notes |
|--------|---------|-------|
| `has_substructure(smarts)` | `boolean` | True if molecule contains the SMARTS pattern |
| `find_substructure(smarts)` | `Array<Uint32Array>` | All matches; each `Uint32Array` = one match's atom indices |
| `get_substructure_atoms(smarts)` | `Uint32Array` | Union of all matched atom indices (for highlighting) |

Supported SMARTS primitives: element symbols (`C`, `c`, `N`, `O` …), atomic number (`[#6]`), H count (`[OH]`, `[NH2]`), aromaticity (`a`, `A`), ring membership (`[R]`), negation (`[!#1]`), formal charge (`[+]`, `[-2]`), bond types (`-`, `=`, `#`, `:`, `~`).

### Bulk export
| Method | Returns | Notes |
|--------|---------|-------|
| `get_positions_flat()` | `Float32Array` | Interleaved [x₀,y₀,z₀, x₁,y₁,z₁, ...] |
| `get_symbols_json()` | `string` | `["O","H","H",...]` |
| `get_atom_info(i)` | `AtomInfo \| null` | Full per-atom object |

### Structure quality
| Method | Returns | Notes |
|--------|---------|-------|
| `chain_breaks(ca_cutoff)` | `ChainBreakRow[]` | Sequence gaps or Cα–Cα > cutoff Å |
| `ramachandran_outliers()` | `RamachandranOutlierRow[]` | φ/ψ outside α/β/left-hand helix regions |

---

## Performance (Apple M-series, Wasm)

| Operation | Time |
|-----------|------|
| Parse 10k-atom PDB | 1.54 ms |
| Parse 10k-atom XYZ | 0.73 ms |
| Build spatial index (10k atoms) | 0.52 ms |
| Radius query — grid (10k atoms) | **1.45 µs** |
| Radius query — linear (10k atoms) | 5.7 µs |
| Compute bonds (500 atoms) | 0.25 ms |

Grid is **4× faster** than linear scan for radius queries. For >5k atoms, `build_spatial_index()` pays for itself immediately.

---

## Usage examples

### SVG 2D rendering

```js
import { loadChem } from '@kent-tokyo/chem-wasm-lens';
const { MolecularSystem } = await loadChem();

const mol = MolecularSystem.from_smiles('c1ccc2ccccc2c1');  // naphthalene
mol.compute_2d_coords();                 // ring-aware 2D layout
const svg = mol.to_svg_string(400, 300); // returns an SVG string

// Insert via DOMParser (parses as XML, appends SVGElement directly)
const svgEl = new DOMParser()
  .parseFromString(svg, 'image/svg+xml').documentElement;
document.getElementById('canvas').appendChild(svgEl);

// Or use a blob URL as <img> src
const src = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml' }));
document.querySelector('img#mol').src = src;

// Check aromaticity
console.log(mol.is_aromatic(0));  // true — lowercase c in SMILES
```

Works with SDF too — explicit bond orders (single/double/triple) render correctly:

```js
const mol = MolecularSystem.from_sdf_string(sdfText);
mol.compute_2d_coords();
const svg = mol.to_svg_string(400, 300);
```

### Framework integration (Vite / React / Next.js)

#### Quick start with the official template

```sh
npx degit kent-tokyo/chem-wasm-lens/template my-mol-app
cd my-mol-app
npm install
npm run dev
```

```ts
// vite.config.ts — no special plugin needed
import { defineConfig } from 'vite';
export default defineConfig({
  optimizeDeps: { exclude: ['@kent-tokyo/chem-wasm-lens'] },
});
```

```tsx
// MolViewer.tsx
import { useEffect, useState } from 'react';
import { loadChem } from '@kent-tokyo/chem-wasm-lens';
import type { MolecularSystem } from '@kent-tokyo/chem-wasm-lens';

let MolSystem: typeof MolecularSystem | null = null;

export function MolViewer({ smiles }: { smiles: string }) {
  const [src, setSrc] = useState('');

  useEffect(() => {
    (async () => {
      if (!MolSystem) ({ MolecularSystem: MolSystem } = await loadChem());
      const mol = MolSystem.from_smiles(smiles);
      mol.compute_2d_coords();
      const blob = new Blob([mol.to_svg_string(400, 300)], { type: 'image/svg+xml' });
      setSrc(URL.createObjectURL(blob));
    })();
  }, [smiles]);

  return src ? <img src={src} width={400} height={300} alt={smiles} /> : null;
}
```

### SDF small molecule

```js
const sdf = `aspirin
  -OEChem-

 21 21  0  0  0  0  0  0  0999 V2000
...
M  END
`;

const mol = MolecularSystem.from_sdf_string(sdf);
console.log(mol.atom_count());   // 21
console.log(mol.bond_count());   // 21 — loaded from SDF bond block directly
```

### SMILES topology

```js
const mol = MolecularSystem.from_smiles('CCO');  // ethanol
console.log(mol.atom_count());   // 9 (2C + 1O + 6H)
console.log(mol.bond_count());   // 8
```

### Bond detection and fragment analysis

```js
mol.compute_bonds();                 // Cordero 2008 covalent radii
const frags = mol.get_fragments();   // Array<Uint32Array>
console.log(`${frags.length} fragment(s)`);  // e.g. 2 for a ligand + solvent

for (const frag of frags) {
  console.log([...frag].map(i => mol.get_symbol(i)).join(''));
}
```

### Kabsch superposition

```js
const ref  = MolecularSystem.from_pdb_string(refPdb);
const mobile = MolecularSystem.from_pdb_string(mobilePdb);

// Align mobile onto ref; coords updated in-place
const rmsd = mobile.superpose(ref);
console.log(`RMSD after alignment: ${rmsd.toFixed(3)} Å`);

// Or query without mutating
const rmsdOnly = mobile.rmsd_aligned(ref);
```

### TypeScript

Full TypeScript types ship with the package:

```ts
import { loadChem } from '@kent-tokyo/chem-wasm-lens';
import type { MolecularSystem, AtomInfo } from '@kent-tokyo/chem-wasm-lens';

const { MolecularSystem } = await loadChem();
const mol: MolecularSystem = MolecularSystem.from_pdb_string(pdbText);
const info: AtomInfo | null = mol.get_atom_info(0);
const neighbors: AtomInfo[] = mol.get_neighbors_info(0, 5.0);
```

---

## Build & test

```sh
# Rust toolchain + wasm-pack
rustup target add wasm32-unknown-unknown
cargo install wasm-pack

# Native unit tests (274 tests)
cargo test

# Linting
cargo clippy

# Wasm package (bundler — Vite / Webpack / esbuild)
wasm-pack build --target bundler

# Wasm package (vanilla web — no bundler)
wasm-pack build --target web

# CDN bundle (Wasm inlined as base64 — run after wasm-pack build --target web)
node scripts/build_cdn.mjs
# → pkg/dist/chem-wasm-lens.esm.js  (~600 KB uncompressed)

# Parallel bond detection on native (optional feature)
cargo test --features parallel
cargo bench --features parallel -- bonds

# Browser tests (requires chromedriver)
wasm-pack test --headless --chrome
```

---

## Architecture

```
src/lib.rs
├── ParseError          — enum (EmptyInput, InvalidAtomCount, ...)
├── SpatialGrid         — HashMap-based uniform voxel grid
├── AtomInfo            — serde::Serialize; shape returned by get_atom_info
├── MolecularSystem     — core struct (#[wasm_bindgen])
│   ├── x / y / z: Vec<f32>           — separate flat vecs; cache-efficient
│   ├── symbols: Vec<String>           — element symbols
│   ├── atom_names / residue_* / chain_ids / hetatm_flags  — PDB metadata
│   ├── bonds: Vec<Vec<usize>>         — adjacency list; lazy via compute_bonds()
│   └── spatial_grid: Option<SpatialGrid>
├── parse_xyz()         — pure-Rust, no JsValue; testable with cargo test
├── parse_pdb()         — fixed-width columns; CONECT records → bonds
├── parse_sdf()         — V2000; bond block → adjacency list
├── parse_smiles()      — organic subset; implicit H; topology only
└── impl MolecularSystem
    ├── Parsers:    from_xyz/pdb/sdf/smiles_string()
    ├── Geometry:   distance, angle, dihedral, center_of_mass, rmsd,
    │               superpose, rmsd_aligned
    ├── Topology:   compute_bonds, get_bonds, bond_count, get_fragments
    ├── Spatial:    build_spatial_index, get_atoms_within_radius,
    │               get_residues_within_radius, get_neighbors_info
    ├── 2D layout:  compute_2d_coords, to_svg_string, is_aromatic
    └── Export:     get_positions_flat, get_symbols_json, get_atom_info
```

**Key design decisions**

- Pure Rust parsers (`parse_xyz`, `parse_pdb`, ...) contain no `JsValue` → fully testable without a browser
- Separate `x/y/z` vectors maximise cache locality for vectorized distance loops
- Entire file content passed as `&str` across JS/Wasm boundary — no per-atom round-trips
- Bonds and spatial grid are built lazily; callers control the cost
- `get_atoms_within_radius` auto-falls back to O(N) scan when grid not built

---

## Demo pages

After `wasm-pack build --target web` and `python3 -m http.server 8080`:

- `examples/index.html` — feature walkthrough (XYZ, PDB, bonds, spatial, SDF)
- `examples/viewer.html` — interactive 3D viewer (ball-and-stick, CPK, ribbon)
- `examples/fetch_demo.html` — live RCSB fetch + φ/ψ dihedral analysis
- `examples/smiles_svg.html` — interactive SMILES → SVG renderer

---

## License

Dual-licensed under **MIT** ([LICENSE-MIT](LICENSE-MIT)) or **Apache 2.0** ([LICENSE-APACHE](LICENSE-APACHE)).
