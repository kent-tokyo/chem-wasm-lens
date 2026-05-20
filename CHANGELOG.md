# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.7] - 2026-05-20

### Changed
- Upgrade `quick-xml` 0.37.5 → 0.40.1

### Fixed
- Resolve 4 clippy warnings: `needless_range_loop`, `collapsible_match` (×2), `manual_is_multiple_of`

## [0.2.6] - 2026-05-20

### Added
- **P57** Ring classification: `get_spiro_atoms`, `get_fused_ring_bonds`, `get_bridged_atoms`
- **P56** E/Z stereochemistry API: `get_ez_stereo_bonds`
- **P54** Atom map API: `get_atom_map_numbers`, `set_atom_map_number`, `clear_atom_map_numbers`

## [0.2.5] - 2026-05-19

### Added
- **P55** `add_hs` — explicit hydrogen addition
- **P55** `fingerprint_topological` — topological (Morgan-style) fingerprint
- **P55** `get_stereo_tags` — stereo tag enumeration

## [0.2.4] - 2026-05-19

### Added
- **P54** `remove_hs` — explicit hydrogen removal
- **P54** `get_descriptors` — molecular descriptor bundle (MW, HBD, HBA, TPSA, rotatable bonds, …)
- **P54** `normalize_depiction` — 2D coordinate normalization
- **P54** SDF property block read/write support

## [0.2.3] - 2026-05-19

### Added
- **P53** Ring-size SMARTS matching
- **P53** Atom-pair fingerprint
- **P53** Template-based 2D alignment (`align_to_template`)

## [0.2.2] - 2026-05-19

### Added
- **P50** Molecular descriptors (LogP, ring count, aromatic atoms, …)
- **P51** 3D conformer generation (`generate_3d_coords`)
- **P52** Reaction execution (`Reaction::run`)

## [0.2.1] - 2026-05-18

### Fixed
- Guard optional vectors in `remove_atom` to prevent index-out-of-bounds panic
- Strip explicit H atoms added by `compute_2d_coords` in the editor kernel

## [0.2.0] - 2026-05-18

### Added
- **P42–P49** Structure editor kernel: atom/bond add-remove, valence validation, undo stack, 2D coordinate generation (`compute_2d_coords`)
- Interactive structure editor example page

## [0.1.0] - 2026-05-17

### Added
- Initial release
- XYZ and PDB parsers (`parse_xyz`, `parse_pdb`)
- SDF / SDF-nth parsers (`parse_sdf`, `parse_sdf_nth`)
- SMILES parser (`parse_smiles`)
- CDXML and CML parsers (`parse_cdxml`, `parse_cml`)
- Bond detection and ring detection
- Spatial queries: `get_atoms_within_radius`, `get_residues_within_radius`
- WebAssembly bindings via `wasm-bindgen`
- GitHub Actions CI (Test & Lint + Wasm Build)

[0.2.7]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.2.6...v0.2.7
[0.2.6]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.2.5...v0.2.6
[0.2.5]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.2.4...v0.2.5
[0.2.4]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.2.3...v0.2.4
[0.2.3]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/kent-tokyo/chem-wasm-lens/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/kent-tokyo/chem-wasm-lens/releases/tag/v0.1.0
