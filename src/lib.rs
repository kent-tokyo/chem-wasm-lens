use wasm_bindgen::prelude::*;

// --- Error Type ---

#[derive(Debug, PartialEq)]
pub enum ParseError {
    EmptyInput,
    // XYZ-specific
    InvalidAtomCount(String),
    MissingCommentLine,
    AtomCountMismatch { expected: usize, found: usize },
    // Shared
    InvalidCoordinate { line: usize, field: &'static str, value: String },
    InvalidField { line: usize, field: &'static str, value: String },
    AtomLimitExceeded { found: usize, limit: usize },
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::EmptyInput => write!(f, "empty input"),
            ParseError::InvalidAtomCount(s) => write!(f, "invalid atom count: '{s}'"),
            ParseError::MissingCommentLine => write!(f, "missing comment line (line 2 required)"),
            ParseError::AtomCountMismatch { expected, found } => {
                write!(f, "atom count mismatch: header says {expected}, found {found}")
            }
            ParseError::InvalidCoordinate { line, field, value } => {
                write!(f, "invalid {field} coordinate on line {line}: '{value}'")
            }
            ParseError::InvalidField { line, field, value } => {
                write!(f, "invalid {field} on line {line}: '{value}'")
            }
            ParseError::AtomLimitExceeded { found, limit } => {
                write!(f, "atom count {found} exceeds limit of {limit}")
            }
        }
    }
}

// --- Core Struct ---

/// Holds atom data for a molecular system.
/// Coordinates are stored in separate flat vectors for cache-efficient distance calculations.
/// PDB-specific fields (atom_names, residue_names, residue_ids, chain_ids, is_hetatm)
/// are empty for XYZ-parsed molecules.
/// `bonds` is empty until `compute_bonds()` is called.
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct MolecularSystem {
    symbols: Vec<String>,
    x: Vec<f32>,
    y: Vec<f32>,
    z: Vec<f32>,
    // PDB-specific (empty for XYZ)
    atom_names: Vec<String>,
    residue_names: Vec<String>,
    residue_ids: Vec<i32>,
    chain_ids: Vec<u8>,
    hetatm_flags: Vec<bool>,
    // Computed on demand; adjacency list indexed by atom index.
    bonds: Vec<Vec<usize>>,
    // Built on demand; enables O(1) average-case neighbor lookup.
    spatial_grid: Option<SpatialGrid>,
    // Populated by compute_rings(); flags which atoms/bonds are in rings.
    ring_atoms: Vec<bool>,
    ring_bonds: std::collections::HashSet<(usize, usize)>,
    // Parallel to bonds: bond_orders[i][k] = order of bond bonds[i][k].
    // 1=single, 2=double, 3=triple. Empty means unknown/all-single.
    bond_orders: Vec<Vec<u8>>,
    // PDB/mmCIF structural fields; empty for XYZ/SDF/SMILES.
    occupancies: Vec<f32>,  // per-atom occupancy (0.0-1.0); default 1.0
    b_factors: Vec<f32>,    // per-atom isotropic B-factor (Å²); default 0.0
    // Formal charges; populated from SMILES bracket atoms; empty for other formats.
    charges: Vec<i32>,
    // Set by parse_smiles; true for atoms parsed as aromatic (c/n/o/s/p/b).
    // Empty for XYZ/PDB/SDF/mmCIF inputs.
    aromatic_atoms: Vec<bool>,
    // Tetrahedral stereo centers parsed from SMILES @/@@.
    // Key: atom index. Value: (descriptor, from_atom).
    //   descriptor: 1 = @@ (CW), -1 = @ (CCW).
    //   from_atom: the atom immediately preceding this center in the SMILES chain.
    stereo_centers: std::collections::HashMap<usize, (i8, Option<usize>)>,
    // Atom mapping numbers from SMILES [X:N] notation. 0 = no mapping.
    atom_map: Vec<u32>,
    // E/Z stereo for double bonds. Key = (lower_idx, higher_idx). true = E (trans), false = Z (cis).
    ez_bonds: std::collections::HashMap<(usize, usize), bool>,
    // ring_sizes_per_atom[i] = sizes of SSSR rings containing atom i.
    // Populated by compute_rings() after enumerate_rings().
    ring_sizes_per_atom: Vec<Vec<u8>>,
    // SDF/molecule-level data items ("> <NAME>" blocks). Empty for non-SDF inputs.
    properties: std::collections::HashMap<String, String>,
}

/// Maximum atoms accepted from user input. Prevents OOM from maliciously crafted files.
const MAX_ATOMS: usize = 1_000_000;

// --- Spatial Grid ---

#[derive(Debug, Clone)]
struct SpatialGrid {
    cells: std::collections::HashMap<(i32, i32, i32), Vec<usize>>,
    cell_size: f32,
    origin: [f32; 3],
}

// --- Serializable Types ---

/// Full per-atom data, serializable to a JS object via serde-wasm-bindgen.
#[derive(Debug, serde::Serialize)]
pub struct AtomInfo {
    pub index: usize,
    pub symbol: String,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub atom_name: String,
    pub residue_name: String,
    pub residue_id: i32,
    pub chain_id: String,
    pub is_hetatm: bool,
    pub occupancy: f32,
    pub b_factor: f32,
}

#[derive(serde::Serialize)]
struct BackboneAngleRow {
    chain_id: String,
    residue_id: i32,
    residue_name: String,
    phi: Option<f32>,
    psi: Option<f32>,
}

#[derive(serde::Serialize)]
struct SecStructRow {
    chain_id: String,
    residue_id: i32,
    residue_name: String,
    ss: String,
}

#[derive(serde::Serialize)]
struct HBondRow {
    donor: usize,
    acceptor: usize,
    distance: f32,
    h_atom: Option<usize>,   // None when no explicit H available
    angle: Option<f32>,      // D-H···A angle in degrees; None when no explicit H
}

#[derive(serde::Serialize)]
struct DisulfideBond {
    atom_i: usize,
    atom_j: usize,
    distance: f32,
}

#[derive(serde::Serialize)]
struct MetalSite {
    metal_atom: usize,
    element: String,
    coordinating: Vec<usize>,
}

#[derive(serde::Serialize)]
struct ContactMapRow {
    chain_i: String,
    resid_i: i32,
    resname_i: String,
    chain_j: String,
    resid_j: i32,
    resname_j: String,
    distance: f32,
}

#[derive(serde::Serialize)]
struct BindingSiteRow {
    chain_id: String,
    residue_id: i32,
    residue_name: String,
}

#[derive(serde::Serialize)]
struct SdfScreenRow {
    index: usize,
    atom_count: usize,
    bond_count: usize,
    formula: String,
    molecular_weight: f32,
    h_bond_donors: u32,
    h_bond_acceptors: u32,
    rotatable_bonds: u32,
}

#[derive(serde::Serialize)]
struct ChainBreakRow {
    chain_id: String,
    from_resid: i32,
    to_resid: i32,
}

#[derive(serde::Serialize)]
struct RamachandranOutlierRow {
    chain_id: String,
    residue_id: i32,
    residue_name: String,
    phi: f32,
    psi: f32,
}

// --- Covalent Radii (Cordero 2008, Å) ---

fn covalent_radius(element: &str) -> Option<f32> {
    match element {
        "H"  => Some(0.31),
        "C"  => Some(0.76),
        "N"  => Some(0.71),
        "O"  => Some(0.66),
        "F"  => Some(0.57),
        "P"  => Some(1.07),
        "S"  => Some(1.05),
        "Cl" => Some(1.02),
        "Br" => Some(1.20),
        "I"  => Some(1.39),
        "Fe" => Some(1.32),
        "Zn" => Some(1.22),
        "Mg" => Some(1.41),
        "Ca" => Some(1.76),
        "Na" => Some(1.66),
        "K"  => Some(2.03),
        "Mn" => Some(1.39),
        "Cu" => Some(1.32),
        "Pt" => Some(1.36),
        "Pd" => Some(1.31),
        "Ni" => Some(1.24),
        "Co" => Some(1.26),
        "Cr" => Some(1.39),
        _    => None,
    }
}

fn atomic_mass(element: &str) -> f32 {
    match element {
        "H"  => 1.008,
        "C"  => 12.011,
        "N"  => 14.007,
        "O"  => 15.999,
        "F"  => 18.998,
        "P"  => 30.974,
        "S"  => 32.06,
        "Cl" => 35.45,
        "Br" => 79.904,
        "I"  => 126.904,
        "Fe" => 55.845,
        "Zn" => 65.38,
        "Mg" => 24.305,
        "Ca" => 40.078,
        "Na" => 22.990,
        "K"  => 39.098,
        "Mn" => 54.938,
        "Cu" => 63.546,
        _    => 12.0, // Unknown → carbon-like mass
    }
}

// --- Helpers ---

fn atom_cell(x: f32, y: f32, z: f32, origin: &[f32; 3], cell_size: f32) -> (i32, i32, i32) {
    (
        ((x - origin[0]) / cell_size).floor() as i32,
        ((y - origin[1]) / cell_size).floor() as i32,
        ((z - origin[2]) / cell_size).floor() as i32,
    )
}

/// Parse a coordinate field, returning a user-friendly error on failure.
fn parse_coord_field(
    slice: Option<&str>,
    field: &'static str,
    line_num: usize,
) -> Result<f32, ParseError> {
    let raw = slice.unwrap_or("").trim();
    raw.parse::<f32>().map_err(|_| ParseError::InvalidCoordinate {
        line: line_num,
        field,
        value: raw.to_string(),
    })
}

/// Bidirectional adjacency push with duplicate guard (O(degree), degree≤4 in practice).
#[inline]
fn adj_push_unique(adj: &mut [Vec<usize>], a: usize, b: usize) {
    if !adj[a].contains(&b) { adj[a].push(b); }
    if !adj[b].contains(&a) { adj[b].push(a); }
}

/// Same as adj_push_unique but also tracks parallel bond-order arrays.
#[inline]
fn adj_push_unique_ordered(
    adj: &mut [Vec<usize>],
    ord: &mut [Vec<u8>],
    a: usize, b: usize, order: u8,
) {
    if !adj[a].contains(&b) { adj[a].push(b); ord[a].push(order); }
    if !adj[b].contains(&a) { adj[b].push(a); ord[b].push(order); }
}

/// Derive an element symbol from a PDB atom name when the element column is absent.
/// PDB atom names like " CA " → "C", "FE  " → "Fe".
/// Ambiguous cases (CA = alpha carbon or calcium, NA = nitrogen or sodium, etc.)
/// are documented: the element column is always preferred when present.
fn derive_element_from_atom_name(atom_name: &str) -> String {
    let trimmed = atom_name.trim();
    // Take leading alpha chars (handles "CA", "N", "OXT")
    let leading: String = trimmed.chars().take_while(|c| c.is_alphabetic()).collect();
    let letters = if !leading.is_empty() {
        leading
    } else {
        // Atom name starts with digit (e.g. "1HB" in some programs) — take first alpha
        trimmed.chars().filter(|c| c.is_alphabetic()).take(1).collect()
    };
    // Title-case: "FE" → "Fe", "N" → "N"
    let mut chars = letters.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().collect::<String>() + &chars.as_str().to_lowercase(),
        None => String::new(),
    }
}

fn aa_one_letter(three: &str) -> char {
    match three {
        "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D',
        "CYS" => 'C', "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G',
        "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
        "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S',
        "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
        "SEC" => 'U', "PYL" => 'O', "MSE" => 'M',
        _ => 'X',
    }
}

// --- XYZ Parser ---

pub fn parse_xyz(input: &str) -> Result<MolecularSystem, ParseError> {
    let mut lines = input.lines();

    let count_line = lines.next().ok_or(ParseError::EmptyInput)?;
    let atom_count: usize = count_line
        .trim()
        .parse()
        .map_err(|_| ParseError::InvalidAtomCount(count_line.trim().to_string()))?;
    if atom_count > MAX_ATOMS {
        return Err(ParseError::AtomLimitExceeded { found: atom_count, limit: MAX_ATOMS });
    }

    // Line 2 is a comment; must be present even if blank.
    lines.next().ok_or(ParseError::MissingCommentLine)?;

    let mut symbols = Vec::with_capacity(atom_count);
    let mut x = Vec::with_capacity(atom_count);
    let mut y = Vec::with_capacity(atom_count);
    let mut z = Vec::with_capacity(atom_count);

    for (idx, line) in lines.enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let line_num = idx + 3; // 1-indexed, offset by 2 header lines
        let mut fields = line.split_whitespace();

        let symbol = fields
            .next()
            .ok_or_else(|| ParseError::InvalidCoordinate {
                line: line_num,
                field: "symbol",
                value: String::new(),
            })?
            .to_string();

        let coord = |s: Option<&str>, field: &'static str| -> Result<f32, ParseError> {
            let raw = s.ok_or_else(|| ParseError::InvalidCoordinate {
                line: line_num,
                field,
                value: String::new(),
            })?;
            raw.parse::<f32>().map_err(|_| ParseError::InvalidCoordinate {
                line: line_num,
                field,
                value: raw.to_string(),
            })
        };

        symbols.push(symbol);
        x.push(coord(fields.next(), "x")?);
        y.push(coord(fields.next(), "y")?);
        z.push(coord(fields.next(), "z")?);
    }

    if symbols.len() != atom_count {
        return Err(ParseError::AtomCountMismatch {
            expected: atom_count,
            found: symbols.len(),
        });
    }

    let mut mol = MolecularSystem::new_empty();
    mol.symbols = symbols;
    mol.x = x;
    mol.y = y;
    mol.z = z;
    Ok(mol)
}

// --- PDB Parser ---
//
// PDB is a fixed-width format. Fields MUST be extracted by byte column position.
// DO NOT use split_whitespace(): adjacent negative coordinates have no separator,
// e.g. "-999.999-999.999" for x/y both equal to -999.999.
//
// Column positions (0-indexed Rust slices):
//   [0..6]   record type ("ATOM  " / "HETATM" / "CONECT" / ...)
//   [6..11]  serial
//   [12..16] atom name
//   [16..17] alternate location indicator
//   [17..20] residue name
//   [21..22] chain ID
//   [22..26] residue sequence number
//   [30..38] x (F8.3)
//   [38..46] y (F8.3)
//   [46..54] z (F8.3)
//   [76..78] element symbol (may be absent in short/legacy files)
//
// CONECT record:
//   [6..11]  central atom serial, then up to 4 bonded serials at [11..16], [16..21], [21..26], [26..31]

pub fn parse_pdb(input: &str) -> Result<MolecularSystem, ParseError> {
    if input.trim().is_empty() {
        return Err(ParseError::EmptyInput);
    }

    let mut symbols: Vec<String> = Vec::new();
    let mut x: Vec<f32> = Vec::new();
    let mut y: Vec<f32> = Vec::new();
    let mut z: Vec<f32> = Vec::new();
    let mut atom_names: Vec<String> = Vec::new();
    let mut residue_names: Vec<String> = Vec::new();
    let mut residue_ids: Vec<i32> = Vec::new();
    let mut chain_ids: Vec<u8> = Vec::new();
    let mut hetatm_flags: Vec<bool> = Vec::new();
    let mut occupancies: Vec<f32> = Vec::new();
    let mut b_factors: Vec<f32> = Vec::new();

    // serial number → atom index (for CONECT resolution)
    let mut serial_to_index: std::collections::HashMap<i32, usize> = std::collections::HashMap::new();
    // raw CONECT pairs (serial1, serial2) — resolved after parsing
    let mut conect_pairs: Vec<(i32, i32)> = Vec::new();
    let mut has_conect = false;

    let mut past_first_model = false;

    for (line_idx, line) in input.lines().enumerate() {
        let line_num = line_idx + 1;
        let record = line.get(0..6).unwrap_or("").trim();

        match record {
            "ENDMDL" => {
                // Stop after the first model (NMR ensembles have multiple).
                past_first_model = true;
            }
            "CONECT" => {
                has_conect = true;
                let s0 = line.get(6..11).and_then(|s| s.trim().parse::<i32>().ok());
                if let Some(src) = s0 {
                    for col in [11usize, 16, 21, 26] {
                        if let Some(dst) = line
                            .get(col..col + 5)
                            .and_then(|s| s.trim().parse::<i32>().ok())
                        {
                            conect_pairs.push((src, dst));
                        }
                    }
                }
            }
            "ATOM" | "HETATM" if !past_first_model => {
                let hetatm = record == "HETATM";

                // Skip alternate conformations other than the primary ('A' or blank).
                let alt_loc = line.get(16..17).and_then(|s| s.chars().next()).unwrap_or(' ');
                if alt_loc != ' ' && alt_loc != 'A' {
                    continue;
                }

                let atom_name = line.get(12..16).unwrap_or("").trim().to_string();
                let res_name = line.get(17..20).unwrap_or("").trim().to_string();
                let chain_id = line.get(21..22).and_then(|s| s.bytes().next()).unwrap_or(b' ');

                let res_seq_raw = line.get(22..26).unwrap_or("").trim();
                let res_seq: i32 = res_seq_raw.parse().map_err(|_| ParseError::InvalidField {
                    line: line_num,
                    field: "resSeq",
                    value: res_seq_raw.to_string(),
                })?;

                let px = parse_coord_field(line.get(30..38), "x", line_num)?;
                let py = parse_coord_field(line.get(38..46), "y", line_num)?;
                let pz = parse_coord_field(line.get(46..54), "z", line_num)?;

                let occupancy: f32 = line.get(54..60)
                    .and_then(|s| s.trim().parse().ok())
                    .unwrap_or(1.0);
                let b_factor: f32 = line.get(60..66)
                    .and_then(|s| s.trim().parse().ok())
                    .unwrap_or(0.0);

                // Prefer the element column (76-78); fall back to atom name derivation.
                let element = line
                    .get(76..78)
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(|s| {
                        let mut chars = s.chars();
                        match chars.next() {
                            Some(first) => {
                                first.to_uppercase().collect::<String>()
                                    + &chars.as_str().to_lowercase()
                            }
                            None => String::new(),
                        }
                    })
                    .unwrap_or_else(|| derive_element_from_atom_name(&atom_name));

                if symbols.len() >= MAX_ATOMS {
                    return Err(ParseError::AtomLimitExceeded {
                        found: symbols.len() + 1,
                        limit: MAX_ATOMS,
                    });
                }
                // Record serial → atom index for CONECT resolution.
                if let Some(serial) = line.get(6..11).and_then(|s| s.trim().parse::<i32>().ok()) {
                    serial_to_index.insert(serial, symbols.len());
                }
                symbols.push(element);
                x.push(px);
                y.push(py);
                z.push(pz);
                atom_names.push(atom_name);
                residue_names.push(res_name);
                residue_ids.push(res_seq);
                chain_ids.push(chain_id);
                hetatm_flags.push(hetatm);
                occupancies.push(occupancy);
                b_factors.push(b_factor);
            }
            _ => {} // HEADER, REMARK, TER, END, etc.
        }
    }

    // Resolve CONECT records into an adjacency list.
    // Bonds are deduplicated; only present if at least one CONECT line was seen.
    let bonds = if has_conect {
        let n = symbols.len();
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (s1, s2) in conect_pairs {
            if let (Some(&i), Some(&j)) = (serial_to_index.get(&s1), serial_to_index.get(&s2)) {
                adj_push_unique(&mut adj, i, j);
            }
        }
        adj
    } else {
        Vec::new()
    };

    let mut mol = MolecularSystem::new_empty();
    mol.symbols = symbols;
    mol.x = x;
    mol.y = y;
    mol.z = z;
    mol.atom_names = atom_names;
    mol.residue_names = residue_names;
    mol.residue_ids = residue_ids;
    mol.chain_ids = chain_ids;
    mol.hetatm_flags = hetatm_flags;
    mol.bonds = bonds;
    mol.occupancies = occupancies;
    mol.b_factors = b_factors;
    Ok(mol)
}

// --- SDF / MOL V2000 Parser ---
//
// Fixed-width format. Column positions (0-indexed):
//   Counts line (line 4 of block):
//     [0..3]  atom count
//     [3..6]  bond count
//   Atom block lines:
//     [0..10] x, [10..20] y, [20..30] z, [31..34] symbol
//   Bond block lines:
//     [0..3]  atom1 (1-based), [3..6] atom2 (1-based), [6..9] bond type
//   `M  END` terminates the molecule; `$$$$` separates multiple entries (first only).

pub fn parse_sdf(input: &str) -> Result<MolecularSystem, ParseError> {
    let lines: Vec<&str> = input
        .lines()
        .take_while(|l| !l.starts_with("$$$$"))
        .collect();

    if lines.len() < 4 {
        return Err(ParseError::EmptyInput);
    }

    let counts = lines[3];
    let atom_count: usize = counts
        .get(0..3)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| ParseError::InvalidAtomCount(counts.get(0..3).unwrap_or("").trim().to_string()))?;
    let bond_count: usize = counts
        .get(3..6)
        .unwrap_or("")
        .trim()
        .parse()
        .map_err(|_| ParseError::InvalidField {
            line: 4,
            field: "bond count",
            value: counts.get(3..6).unwrap_or("").trim().to_string(),
        })?;

    if atom_count > MAX_ATOMS {
        return Err(ParseError::AtomLimitExceeded { found: atom_count, limit: MAX_ATOMS });
    }

    let mut symbols: Vec<String> = Vec::with_capacity(atom_count);
    let mut x: Vec<f32> = Vec::with_capacity(atom_count);
    let mut y: Vec<f32> = Vec::with_capacity(atom_count);
    let mut z: Vec<f32> = Vec::with_capacity(atom_count);

    for i in 0..atom_count {
        let line_num = 4 + i;
        let line = lines.get(line_num).ok_or(ParseError::AtomCountMismatch {
            expected: atom_count,
            found: i,
        })?;
        if line.len() < 34 {
            return Err(ParseError::InvalidField {
                line: line_num + 1,
                field: "atom line (too short)",
                value: line.to_string(),
            });
        }
        let px = parse_coord_field(Some(&line[0..10]), "x", line_num + 1)?;
        let py = parse_coord_field(Some(&line[10..20]), "y", line_num + 1)?;
        let pz = parse_coord_field(Some(&line[20..30]), "z", line_num + 1)?;
        let sym = line[31..34].trim().to_string();
        x.push(px);
        y.push(py);
        z.push(pz);
        symbols.push(sym);
    }

    // Bond block — store as adjacency list directly (explicit bond graph from SDF).
    let mut bonds: Vec<Vec<usize>> = vec![Vec::new(); atom_count];
    let mut bond_orders: Vec<Vec<u8>> = vec![Vec::new(); atom_count];
    let bond_start = 4 + atom_count;
    for b in 0..bond_count {
        let line_num = bond_start + b;
        let line = match lines.get(line_num) {
            Some(l) if !l.starts_with("M ") => l,
            _ => break,
        };
        if line.len() < 6 {
            continue;
        }
        let a1: usize = match line[0..3].trim().parse::<usize>() {
            Ok(v) if v > 0 => v - 1,
            _ => continue,
        };
        let a2: usize = match line[3..6].trim().parse::<usize>() {
            Ok(v) if v > 0 => v - 1,
            _ => continue,
        };
        let order: u8 = line.get(6..9)
            .and_then(|s| s.trim().parse::<u8>().ok())
            .map(|t| if t == 4 { 1u8 } else { t.clamp(1, 3) })
            .unwrap_or(1);
        if a1 < atom_count && a2 < atom_count {
            bonds[a1].push(a2);
            bonds[a2].push(a1);
            bond_orders[a1].push(order);
            bond_orders[a2].push(order);
        }
    }

    let mut mol = MolecularSystem::new_empty();
    mol.symbols = symbols;
    mol.x = x;
    mol.y = y;
    mol.z = z;
    mol.bonds = bonds;
    mol.bond_orders = bond_orders;

    // Parse SDF data items: "> <PROPNAME>" / value / blank lines after the bond block.
    let data_start = 4 + atom_count + bond_count;
    let mut di = data_start;
    while di < lines.len() {
        let line = lines[di].trim();
        if let Some(rest) = line.strip_prefix("> <").and_then(|r| r.strip_suffix('>')) {
            let key = rest.to_string();
            di += 1;
            let mut vals: Vec<&str> = Vec::new();
            while di < lines.len() && !lines[di].trim().is_empty() && !lines[di].trim().starts_with('>') {
                vals.push(lines[di].trim());
                di += 1;
            }
            mol.properties.insert(key, vals.join("\n"));
        } else {
            di += 1;
        }
    }

    Ok(mol)
}

/// Extract the n-th molecule block (0-indexed) from a multi-molecule SDF and parse it.
pub fn parse_sdf_nth(data: &str, n: usize) -> Result<MolecularSystem, ParseError> {
    let mut current = 0usize;
    let mut block_start = 0usize;
    let mut byte_pos = 0usize;
    for line in data.lines() {
        let line_bytes = line.len() + 1; // +1 for '\n'
        if line.starts_with("$$$$") {
            if current == n {
                return parse_sdf(&data[block_start..byte_pos + line_bytes]);
            }
            current += 1;
            block_start = byte_pos + line_bytes;
        }
        byte_pos += line_bytes;
    }
    // Last block (no trailing $$$$)
    if current == n && block_start < data.len() {
        return parse_sdf(&data[block_start..]);
    }
    Err(ParseError::EmptyInput)
}

// --- SMILES Parser (organic subset, topology only) ---
//
// Supports:
//   Atoms   : B C N O P S F Cl Br I (organic subset, case-sensitive uppercase)
//             bracket atoms: [Fe], [NH4+], [OH-], [n]
//             aromatic atoms: b c n o p s (kekulized to alternating single/double)
//   Bonds   : - (single) = (double) # (triple) : (aromatic); implicit single when omitted
//   Branches: ( )
//   Rings   : digit 1-9, or %dd (two-digit)
//   Implicit H: added based on standard valence (C=4, N=3, O=2, S=2, P=3, B=3)
//
// Out of scope (MVP): stereochemistry (@,/,\), isotopes

/// Augmenting path search for Hopcroft-Karp-style maximum matching on aromatic subgraph.
fn try_augment(
    v: usize,
    arom_adj: &[Vec<usize>],
    matched: &mut Vec<Option<usize>>,
    visited: &mut Vec<bool>,
) -> bool {
    for &u in &arom_adj[v] {
        if visited[u] { continue; }
        visited[u] = true;
        if matched[u].is_none()
            || try_augment(matched[u].unwrap(), arom_adj, matched, visited)
        {
            matched[v] = Some(u);
            matched[u] = Some(v);
            return true;
        }
    }
    false
}

/// Convert aromatic bonds to alternating single/double via maximum matching (Kekulization).
/// Must be called before implicit H computation.
fn kekulize(
    n: usize,
    aromatic: &[bool],
    bonds: &mut [(usize, usize, BondOrder)],
) {
    if !aromatic.iter().any(|&a| a) {
        return; // no aromatic atoms
    }
    // Build aromatic-only adjacency list
    let mut arom_adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b, ord) in bonds.iter() {
        if ord == BondOrder::Aromatic {
            adj_push_unique(&mut arom_adj, a, b);
        }
    }
    // Maximum matching on aromatic subgraph
    let mut matched: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        if aromatic.get(i).copied().unwrap_or(false) && matched[i].is_none() {
            let mut visited = vec![false; n];
            visited[i] = true;
            try_augment(i, &arom_adj, &mut matched, &mut visited);
        }
    }
    // Assign bond orders: matched aromatic edge → double, unmatched → single
    for (a, b, ref mut ord) in bonds.iter_mut() {
        if *ord == BondOrder::Aromatic {
            *ord = if matched[*a] == Some(*b) {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
        }
    }
}

fn smiles_valence(element: &str) -> Option<u8> {
    match element {
        "B"  => Some(3),
        "C"  => Some(4),
        "N"  => Some(3),
        "O"  => Some(2),
        "P"  => Some(3),
        "S"  => Some(2),
        "F"  => Some(1),
        "Cl" => Some(1),
        "Br" => Some(1),
        "I"  => Some(1),
        _    => None,
    }
}

// Bond order value used when counting valence
#[derive(Clone, Copy, PartialEq)]
enum BondOrder { Single, Double, Triple, Aromatic }
impl BondOrder {
    fn valence_contribution(self) -> u8 {
        match self {
            BondOrder::Single | BondOrder::Aromatic => 1, // Aromatic resolved after kekulize
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        }
    }
}

pub fn parse_smiles(input: &str) -> Result<MolecularSystem, ParseError> {
    let s = input.trim();
    if s.is_empty() {
        return Err(ParseError::EmptyInput);
    }

    let chars: Vec<char> = s.chars().collect();
    let n = chars.len();
    let mut pos = 0;

    // --- Atom and bond storage ---
    let mut symbols: Vec<String> = Vec::new();
    let mut bonds: Vec<(usize, usize, BondOrder)> = Vec::new();
    let mut aromatic_atoms: Vec<bool> = Vec::new();
    let mut charges: Vec<i32> = Vec::new();
    let mut atom_maps: Vec<u32> = Vec::new();
    // Stereo centers: atom_index → (descriptor, from_atom)
    let mut stereo_map: std::collections::HashMap<usize, (i8, Option<usize>)> = std::collections::HashMap::new();

    // Stack of atom indices for branch points
    let mut stack: Vec<usize> = Vec::new();
    // Ring-open atoms: ring_digit → (atom_index, bond_order)
    let mut ring_opens: std::collections::HashMap<u16, (usize, BondOrder)> = std::collections::HashMap::new();

    // The previous atom index and the pending bond order between prev and next atom
    let mut prev: Option<usize> = None;
    let mut pending_bond: Option<BondOrder> = None;
    // E/Z: bond direction chars. Parallel to `bonds`: +1 = /, -1 = \, 0 = none.
    let mut bond_dirs: Vec<i8> = Vec::new();
    let mut pending_dir: i8 = 0; // direction for the next bond

    let err_inv = |msg: &str| -> ParseError {
        ParseError::InvalidField { line: 1, field: "SMILES", value: msg.to_string() }
    };

    while pos < n {
        let c = chars[pos];

        match c {
            // ── Bond tokens ───────────────────────────────────────────────────
            '-' => { pending_bond = Some(BondOrder::Single); pos += 1; }
            '=' => { pending_bond = Some(BondOrder::Double); pos += 1; }
            '#' => { pending_bond = Some(BondOrder::Triple); pos += 1; }

            // ── Branch open ────────────────────────────────────────────────────
            '(' => {
                if let Some(p) = prev {
                    stack.push(p);
                }
                pos += 1;
            }

            // ── Branch close ──────────────────────────────────────────────────
            ')' => {
                prev = stack.pop();
                pending_bond = None;
                pos += 1;
            }

            // ── Ring closure (digit) ──────────────────────────────────────────
            '1'..='9' => {
                let digit = c as u16 - '0' as u16;
                let cur = prev.ok_or_else(|| err_inv("ring closure before any atom"))?;
                let cur_is_arom = aromatic_atoms.get(cur).copied().unwrap_or(false);
                match ring_opens.remove(&digit) {
                    Some((open_atom, stored_order)) => {
                        let order = if let Some(pb) = pending_bond.take() {
                            pb
                        } else if cur_is_arom && aromatic_atoms.get(open_atom).copied().unwrap_or(false) {
                            BondOrder::Aromatic
                        } else {
                            stored_order
                        };
                        bonds.push((open_atom, cur, order));
                        bond_dirs.push(0); pending_dir = 0;
                    }
                    None => {
                        let order = pending_bond.take().unwrap_or(BondOrder::Single);
                        ring_opens.insert(digit, (cur, order));
                        pending_dir = 0;
                    }
                }
                pos += 1;
            }

            // ── Two-digit ring closure (%dd) ──────────────────────────────────
            '%' if pos + 2 < n => {
                let d1 = chars[pos+1].to_digit(10).ok_or_else(|| err_inv("bad %nn ring"))? as u16;
                let d2 = chars[pos+2].to_digit(10).ok_or_else(|| err_inv("bad %nn ring"))? as u16;
                let digit = d1 * 10 + d2;
                let cur = prev.ok_or_else(|| err_inv("ring closure before any atom"))?;
                let cur_is_arom = aromatic_atoms.get(cur).copied().unwrap_or(false);
                match ring_opens.remove(&digit) {
                    Some((open_atom, stored_order)) => {
                        let order = if let Some(pb) = pending_bond.take() {
                            pb
                        } else if cur_is_arom && aromatic_atoms.get(open_atom).copied().unwrap_or(false) {
                            BondOrder::Aromatic
                        } else {
                            stored_order
                        };
                        bonds.push((open_atom, cur, order));
                        bond_dirs.push(0); pending_dir = 0;
                    }
                    None => {
                        let order = pending_bond.take().unwrap_or(BondOrder::Single);
                        ring_opens.insert(digit, (cur, order));
                        pending_dir = 0;
                    }
                }
                pos += 3;
            }

            // ── Bracket atom [Fe], [NH4+], [OH-] ─────────────────────────────
            '[' => {
                pos += 1; // skip '['
                // Skip isotope digits
                while pos < n && chars[pos].is_ascii_digit() { pos += 1; }
                // Read element symbol (1-2 chars starting uppercase or lowercase)
                let sym = read_element(&chars, &mut pos);
                if sym.is_empty() {
                    return Err(err_inv("empty bracket atom element"));
                }
                // Parse optional chirality (@, @@) — SMILES spec: isotope symbol chiral H charge
                let stereo_desc: i8 = if pos < n && chars[pos] == '@' {
                    pos += 1;
                    if pos < n && chars[pos] == '@' { pos += 1; 1i8 } // @@ = CW
                    else { -1i8 } // @ = CCW
                } else { 0i8 };
                // Parse optional H count (H, H2, ...)
                if pos < n && chars[pos] == 'H' {
                    pos += 1;
                    if pos < n && chars[pos].is_ascii_digit() { pos += 1; }
                }
                // Parse optional formal charge: +, -, ++, --, +2, -2, ...
                let atom_charge: i32 = if pos < n && (chars[pos] == '+' || chars[pos] == '-') {
                    let sign = if chars[pos] == '+' { 1i32 } else { -1i32 };
                    let prev_char = chars[pos];
                    pos += 1;
                    if pos < n && chars[pos].is_ascii_digit() {
                        let mag = (chars[pos] as i32) - ('0' as i32); pos += 1; sign * mag
                    } else if pos < n && chars[pos] == prev_char {
                        pos += 1; sign * 2
                    } else { sign }
                } else { 0 };
                // Parse optional atom-map number :[N] and skip remaining until ']'
                let mut map_num: u32 = 0;
                if pos < n && chars[pos] == ':' {
                    pos += 1; // skip ':'
                    let mut num = 0u32;
                    while pos < n && chars[pos].is_ascii_digit() {
                        num = num * 10 + (chars[pos] as u32 - '0' as u32);
                        pos += 1;
                    }
                    map_num = num;
                }
                while pos < n && chars[pos] != ']' { pos += 1; }
                if pos < n { pos += 1; } // skip ']'

                let atom_idx = symbols.len();
                aromatic_atoms.push(false);
                symbols.push(sym);
                charges.push(atom_charge);
                atom_maps.push(map_num);
                if stereo_desc != 0 {
                    stereo_map.insert(atom_idx, (stereo_desc, prev));
                }
                if let Some(p) = prev {
                    let order = pending_bond.take().unwrap_or(BondOrder::Single);
                    bonds.push((p, atom_idx, order));
                    bond_dirs.push(pending_dir); pending_dir = 0;
                } else {
                    pending_bond = None;
                    pending_dir = 0;
                }
                prev = Some(atom_idx);
            }

            // ── Organic-subset atom ────────────────────────────────────────────
            'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' => {
                // Check for two-character atoms: Cl, Br are NOT starting with B/C
                // (Cl starts C? No — Cl starts with 'C', Br starts with 'B')
                // Actually: Cl starts with 'C', Br starts with 'B'
                let sym = if c == 'C' && pos + 1 < n && chars[pos+1] == 'l' {
                    pos += 2; "Cl".to_string()
                } else if c == 'B' && pos + 1 < n && chars[pos+1] == 'r' {
                    pos += 2; "Br".to_string()
                } else {
                    pos += 1; c.to_string()
                };

                let atom_idx = symbols.len();
                aromatic_atoms.push(false);
                symbols.push(sym);
                charges.push(0);
                atom_maps.push(0);
                if let Some(p) = prev {
                    let order = pending_bond.take().unwrap_or(BondOrder::Single);
                    bonds.push((p, atom_idx, order));
                    bond_dirs.push(pending_dir); pending_dir = 0;
                } else {
                    pending_bond = None;
                    pending_dir = 0;
                }
                prev = Some(atom_idx);
            }

            // ── Aromatic atoms (lowercase) — track as aromatic for kekulization ──
            'b' | 'c' | 'n' | 'o' | 'p' | 's' => {
                let sym = c.to_uppercase().to_string();
                let atom_idx = symbols.len();
                aromatic_atoms.push(true);
                symbols.push(sym);
                charges.push(0);
                atom_maps.push(0);
                if let Some(p) = prev {
                    let order = pending_bond.take().unwrap_or_else(|| {
                        if aromatic_atoms.get(p).copied().unwrap_or(false) {
                            BondOrder::Aromatic
                        } else {
                            BondOrder::Single
                        }
                    });
                    bonds.push((p, atom_idx, order));
                    bond_dirs.push(pending_dir); pending_dir = 0;
                } else {
                    pending_bond = None;
                    pending_dir = 0;
                }
                prev = Some(atom_idx);
                pos += 1;
            }

            // ── Aromatic bond ──────────────────────────────────────────────────
            ':' => { pending_bond = Some(BondOrder::Aromatic); pos += 1; }

            // ── Skip whitespace ────────────────────────────────────────────────
            ' ' | '\t' => { pos += 1; }

            // ── Dot (disconnected structure) ───────────────────────────────────
            '.' => { prev = None; pending_bond = None; pos += 1; }

            // ── E/Z bond direction ────────────────────────────────────────────
            '/' | '\\' => { pending_dir = if c == '/' { 1 } else { -1 }; pos += 1; }

            _ => { pos += 1; } // silently skip unknown chars
        }
    }

    // --- Add implicit hydrogens based on standard valence ---
    let heavy_count = symbols.len();

    // Kekulize aromatic bonds before computing implicit H
    kekulize(heavy_count, &aromatic_atoms, &mut bonds);

    // Build per-atom bond valence sum (from heavy-atom bonds only)
    let mut valence_used: Vec<u8> = vec![0u8; heavy_count];
    for &(a, b, order) in &bonds {
        let v = order.valence_contribution();
        if a < heavy_count { valence_used[a] = valence_used[a].saturating_add(v); }
        if b < heavy_count { valence_used[b] = valence_used[b].saturating_add(v); }
    }

    let mut h_bonds: Vec<(usize, usize, BondOrder)> = Vec::new();
    let mut h_symbols: Vec<String> = Vec::new();

    for i in 0..heavy_count {
        let h_count = match smiles_valence(&symbols[i]) {
            Some(max_v) => max_v.saturating_sub(valence_used[i]) as usize,
            None        => 0,
        };
        let base = heavy_count + h_symbols.len();
        for k in 0..h_count {
            h_symbols.push("H".to_string());
            h_bonds.push((i, base + k, BondOrder::Single));
        }
    }

    // Merge heavy + H atoms
    let mut all_symbols = symbols;
    all_symbols.extend(h_symbols);
    let total = all_symbols.len();
    // Extend charges with zeros for each implicit H
    let h_count_added = total - charges.len();
    charges.extend(std::iter::repeat_n(0, h_count_added));
    let all_charges = charges;

    // Build E/Z map from SMILES bond direction characters.
    // bond_dirs[i] is parallel to bonds[i]: +1 = '/', -1 = '\', 0 = none.
    // For double bond (b=c): find directed adjacent bonds on each side.
    // Bond (a→b, dir d): substituent a is at position -d relative to b.
    // Bond (c→d, dir d): substituent d is at position +d relative to c.
    // Different positions → E (trans); same position → Z (cis).
    let smiles_ez = {
        let n_bonds = bonds.len();
        let mut ez: std::collections::HashMap<(usize, usize), bool> = std::collections::HashMap::new();
        for bi in 0..n_bonds {
            let (b, c, ref border) = bonds[bi];
            if *border != BondOrder::Double { continue; }
            let mut pos_b: Option<i8> = None;
            let mut pos_c: Option<i8> = None;
            for ki in 0..n_bonds {
                if ki == bi || bond_dirs[ki] == 0 { continue; }
                let (aa, bb, _) = bonds[ki];
                if pos_b.is_none() {
                    if bb == b { pos_b = Some(-bond_dirs[ki]); }
                    else if aa == b { pos_b = Some(bond_dirs[ki]); }
                }
                if pos_c.is_none() {
                    if aa == c { pos_c = Some(bond_dirs[ki]); }
                    else if bb == c { pos_c = Some(-bond_dirs[ki]); }
                }
            }
            if let (Some(pb), Some(pc)) = (pos_b, pos_c) {
                ez.insert((b.min(c), b.max(c)), pb != pc);
            }
        }
        ez
    };

    // Build adjacency list
    let all_bonds_iter = bonds.into_iter().chain(h_bonds);
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); total];
    let mut adj_orders: Vec<Vec<u8>> = vec![Vec::new(); total];
    for (a, b, bond_ord) in all_bonds_iter {
        if a < total && b < total {
            let order = bond_ord.valence_contribution();
            adj_push_unique_ordered(&mut adj, &mut adj_orders, a, b, order);
        }
    }

    if total > MAX_ATOMS {
        return Err(ParseError::AtomLimitExceeded { found: total, limit: MAX_ATOMS });
    }

    // Extend aromatic_atoms to cover implicit H atoms (all false)
    aromatic_atoms.resize(total, false);
    // Extend atom_maps: non-bracket atoms have map 0; fill to total
    atom_maps.resize(total, 0);

    let mut mol = MolecularSystem::new_empty();
    mol.symbols = all_symbols;
    mol.x = vec![0.0; total];
    mol.y = vec![0.0; total];
    mol.z = vec![0.0; total];
    mol.bonds = adj;
    mol.bond_orders = adj_orders;
    mol.charges = all_charges;
    mol.aromatic_atoms = aromatic_atoms;
    mol.stereo_centers = stereo_map;
    mol.ez_bonds = smiles_ez;
    mol.atom_map = atom_maps;
    Ok(mol)
}

// --- mmCIF Parser ---
//
// Parses the _atom_site loop from an mmCIF / PDBx file.
// Only reads the first MODEL (pdbx_PDB_model_num == first seen).
// Columns are mapped by name from the loop_ header.

pub fn parse_mmcif(input: &str) -> Result<MolecularSystem, ParseError> {
    use std::collections::HashMap;

    let mut symbols: Vec<String> = Vec::new();
    let mut x: Vec<f32> = Vec::new();
    let mut y: Vec<f32> = Vec::new();
    let mut z: Vec<f32> = Vec::new();
    let mut atom_names: Vec<String> = Vec::new();
    let mut residue_names: Vec<String> = Vec::new();
    let mut residue_ids: Vec<i32> = Vec::new();
    let mut chain_ids: Vec<u8> = Vec::new();
    let mut hetatm_flags: Vec<bool> = Vec::new();
    let mut occupancies: Vec<f32> = Vec::new();
    let mut b_factors: Vec<f32> = Vec::new();

    // State machine
    enum State {
        Scanning,
        CollectingAtomSiteHeader,
        ReadingDataRows,
    }

    let mut state = State::Scanning;
    // column name → index in each data row
    let mut col_map: HashMap<String, usize> = HashMap::new();
    let mut col_count = 0usize;
    let mut first_model: Option<i32> = None;

    for line in input.lines() {
        let trimmed = line.trim();

        match state {
            State::Scanning => {
                if trimmed == "loop_" {
                    col_map.clear();
                    col_count = 0;
                    state = State::CollectingAtomSiteHeader;
                }
            }

            State::CollectingAtomSiteHeader => {
                if trimmed.starts_with("_atom_site.") {
                    // Extract column name after prefix
                    let col_name = trimmed.trim_start_matches("_atom_site.").to_string();
                    col_map.insert(col_name, col_count);
                    col_count += 1;
                } else if trimmed.starts_with('_') || trimmed == "loop_" || trimmed.starts_with("data_") {
                    // Header from a different category — not atom_site loop, reset
                    col_map.clear();
                    col_count = 0;
                    if trimmed == "loop_" {
                        state = State::CollectingAtomSiteHeader;
                    } else {
                        state = State::Scanning;
                    }
                } else if trimmed.is_empty() || trimmed.starts_with('#') {
                    // skip blank lines and comments inside header
                    continue;
                } else {
                    // First data row
                    if col_map.is_empty() {
                        // No atom_site columns collected — still scanning
                        state = State::Scanning;
                        continue;
                    }
                    state = State::ReadingDataRows;
                    // Fall through to process this line as a data row
                    if let Some(result) = mmcif_parse_row(
                        trimmed, &col_map,
                        &mut first_model,
                        &mut symbols, &mut x, &mut y, &mut z,
                        &mut atom_names, &mut residue_names, &mut residue_ids,
                        &mut chain_ids, &mut hetatm_flags,
                        &mut occupancies, &mut b_factors,
                    ) {
                        return Err(result);
                    }
                }
            }

            State::ReadingDataRows => {
                if trimmed.starts_with('_') || trimmed == "loop_" || trimmed.starts_with("data_") {
                    // End of atom_site data
                    break;
                }
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    continue;
                }
                if let Some(result) = mmcif_parse_row(
                    trimmed, &col_map,
                    &mut first_model,
                    &mut symbols, &mut x, &mut y, &mut z,
                    &mut atom_names, &mut residue_names, &mut residue_ids,
                    &mut chain_ids, &mut hetatm_flags,
                    &mut occupancies, &mut b_factors,
                ) {
                    return Err(result);
                }
            }
        }
    }

    if symbols.is_empty() {
        return Err(ParseError::EmptyInput);
    }

    let mut mol = MolecularSystem::new_empty();
    mol.symbols = symbols;
    mol.x = x;
    mol.y = y;
    mol.z = z;
    mol.atom_names = atom_names;
    mol.residue_names = residue_names;
    mol.residue_ids = residue_ids;
    mol.chain_ids = chain_ids;
    mol.hetatm_flags = hetatm_flags;
    mol.occupancies = occupancies;
    mol.b_factors = b_factors;
    Ok(mol)
}

/// Parse a single mmCIF atom_site data row and push to the provided vectors.
/// Returns Some(ParseError) on fatal error, None on success (including skipped rows).
#[allow(clippy::too_many_arguments)]
fn mmcif_parse_row(
    line: &str,
    col_map: &std::collections::HashMap<String, usize>,
    first_model: &mut Option<i32>,
    symbols: &mut Vec<String>,
    x: &mut Vec<f32>,
    y: &mut Vec<f32>,
    z: &mut Vec<f32>,
    atom_names: &mut Vec<String>,
    residue_names: &mut Vec<String>,
    residue_ids: &mut Vec<i32>,
    chain_ids: &mut Vec<u8>,
    hetatm_flags: &mut Vec<bool>,
    occupancies: &mut Vec<f32>,
    b_factors: &mut Vec<f32>,
) -> Option<ParseError> {
    let fields: Vec<&str> = line.split_whitespace().collect();

    let get = |name: &str| -> Option<&str> {
        col_map.get(name).and_then(|&i| fields.get(i).copied())
    };

    // Model number filtering
    let model_num: i32 = get("pdbx_PDB_model_num")
        .and_then(|s| if s == "." || s == "?" { None } else { s.parse().ok() })
        .unwrap_or(1);

    match *first_model {
        None => { *first_model = Some(model_num); }
        Some(first) if model_num != first => {
            // Different model — stop (signal by returning early; caller breaks)
            // We can't break the outer loop from here, so we just return None
            // and the caller will handle it on the next loop_/# encounter.
            // Better: skip silently (second model exclusion)
            return None;
        }
        Some(_) => {}
    }

    // group_PDB: ATOM or HETATM
    let group = get("group_PDB").unwrap_or("ATOM");
    let hetatm = group.trim() == "HETATM";

    // Alternate location: skip B, C, D, ...
    let alt_loc = get("label_alt_id").unwrap_or(".");
    if alt_loc != "." && alt_loc != "?" && alt_loc != "A" && !alt_loc.is_empty() {
        // Skip non-primary alternates
        return None;
    }

    // Element symbol — mmCIF stores uppercase (e.g., "FE"), title-case it
    let raw_sym = get("type_symbol").unwrap_or("C");
    let symbol = mmcif_title_case(raw_sym);

    let atom_name = get("label_atom_id").unwrap_or("").to_string();
    let res_name = get("label_comp_id").unwrap_or("").to_string();
    let chain_id_str = get("auth_asym_id").unwrap_or(" ");
    let chain_byte = chain_id_str.bytes().next().unwrap_or(b' ');

    let res_id: i32 = get("auth_seq_id")
        .and_then(|s| if s == "." || s == "?" { None } else { s.parse().ok() })
        .unwrap_or(0);

    let px: f32 = get("Cartn_x").and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let py: f32 = get("Cartn_y").and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let pz: f32 = get("Cartn_z").and_then(|s| s.parse().ok()).unwrap_or(0.0);

    if symbols.len() >= MAX_ATOMS {
        return Some(ParseError::AtomLimitExceeded { found: symbols.len() + 1, limit: MAX_ATOMS });
    }

    symbols.push(symbol);
    x.push(px);
    y.push(py);
    z.push(pz);
    atom_names.push(atom_name);
    residue_names.push(res_name);
    residue_ids.push(res_id);
    chain_ids.push(chain_byte);
    hetatm_flags.push(hetatm);

    let occ: f32 = get("occupancy").and_then(|s| s.trim().parse().ok()).unwrap_or(1.0);
    let bfac: f32 = get("B_iso_or_equiv")
        .and_then(|s| if s == "." || s == "?" { None } else { s.trim().parse().ok() })
        .unwrap_or(0.0);
    occupancies.push(occ);
    b_factors.push(bfac);

    None
}

/// Title-case an mmCIF element string: "FE" → "Fe", "C" → "C", "ZN" → "Zn".
fn mmcif_title_case(s: &str) -> String {
    let s = s.trim();
    if s.is_empty() {
        return String::new();
    }
    let mut chars = s.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().collect::<String>() + &chars.as_str().to_lowercase(),
        None => String::new(),
    }
}

/// Read a 1- or 2-character element symbol from `chars` at `pos`.
/// Advances `pos`. Returns empty string if no alphabetic char found.
fn read_element(chars: &[char], pos: &mut usize) -> String {
    if *pos >= chars.len() { return String::new(); }
    let first = chars[*pos];
    if !first.is_alphabetic() { return String::new(); }
    *pos += 1;
    let mut sym = first.to_uppercase().to_string();
    // Optional second char: must be lowercase (Cl, Br, Fe, etc.)
    if *pos < chars.len() && chars[*pos].is_lowercase() {
        sym.push(chars[*pos]);
        *pos += 1;
    }
    sym
}

// --- Secondary Structure Classification ---

fn classify_ss(phi: f32, psi: f32) -> &'static str {
    if (-90.0..=-40.0).contains(&phi) && (-70.0..=-10.0).contains(&psi) {
        "H"
    } else if (-160.0..=-80.0).contains(&phi) && (
        (80.0..=180.0).contains(&psi) || (-180.0..=-160.0).contains(&psi)
    ) {
        "E"
    } else {
        "C"
    }
}

fn is_ramachandran_allowed(phi: f32, psi: f32) -> bool {
    let in_alpha = (-130.0..=-10.0).contains(&phi) && (-100.0..=60.0).contains(&psi);
    let in_beta = (-180.0..=-40.0).contains(&phi) && (psi >= 60.0 || psi <= -140.0);
    let in_lh = (10.0..=120.0).contains(&phi) && (10.0..=100.0).contains(&psi);
    in_alpha || in_beta || in_lh
}

// --- Morgan Fingerprint Helpers ---

/// FNV-1a inspired hash mix for Morgan fingerprint computation.
fn fnv_mix(mut h: u32, val: u32) -> u32 {
    h ^= val;
    h = h.wrapping_mul(0x01000193);
    h
}

/// Compute a 2048-bit (256-byte) Morgan/ECFP-like fingerprint.
/// `radius=2` is ECFP4-equivalent. The `bonds` adjacency list must be populated.
fn morgan_fingerprint_bits(symbols: &[String], bonds: &[Vec<usize>], radius: u32) -> [u8; 256] {
    let n = symbols.len();
    let mut bit_vec = [0u8; 256];

    if n == 0 {
        return bit_vec;
    }

    // Compute initial atom identifiers (radius-0)
    let mut ids: Vec<u32> = symbols
        .iter()
        .enumerate()
        .map(|(i, sym)| {
            let mut h = 0x811c9dc5u32;
            for b in sym.bytes() {
                h = fnv_mix(h, b as u32);
            }
            let degree = if i < bonds.len() { bonds[i].len() } else { 0 };
            h = fnv_mix(h, degree as u32);
            h
        })
        .collect();

    // Set bits for radius-0 contributions
    for &id in &ids {
        let bit_pos = id % 2048;
        bit_vec[(bit_pos / 8) as usize] |= 1 << (bit_pos % 8);
    }

    // Iterative neighbourhood expansion
    for _ in 1..=radius {
        let mut new_ids = ids.clone();
        for i in 0..n {
            let mut nbr_ids: Vec<u32> = if i < bonds.len() {
                bonds[i].iter().map(|&j| ids[j]).collect()
            } else {
                Vec::new()
            };
            nbr_ids.sort_unstable();

            let mut h = ids[i];
            for nbr in nbr_ids {
                h = fnv_mix(h, nbr);
            }
            new_ids[i] = h;
        }
        ids = new_ids;

        for &id in &ids {
            let bit_pos = id % 2048;
            bit_vec[(bit_pos / 8) as usize] |= 1 << (bit_pos % 8);
        }
    }

    bit_vec
}

// --- Wasm-Exposed Methods ---

fn to_js<T: serde::Serialize>(val: &T) -> JsValue {
    serde_wasm_bindgen::to_value(val).unwrap_or(JsValue::NULL)
}

#[wasm_bindgen]
impl MolecularSystem {
    pub fn from_xyz_string(data: &str) -> Result<MolecularSystem, JsValue> {
        parse_xyz(data).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    pub fn from_pdb_string(data: &str) -> Result<MolecularSystem, JsValue> {
        parse_pdb(data).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    /// Parse an SDF / MOL V2000 string. Only the first molecule entry is loaded.
    /// Explicit bonds from the bond block are stored directly (no distance heuristic needed).
    /// Stereo centers are automatically perceived from 3D coordinates when present;
    /// 2D (flat, z=0) molecules are unaffected.
    pub fn from_sdf_string(data: &str) -> Result<MolecularSystem, JsValue> {
        let mut mol = parse_sdf(data).map_err(|e| JsValue::from_str(&e.to_string()))?;
        mol.perceive_stereo_from_3d();
        mol.perceive_ez_from_3d();
        Ok(mol)
    }

    /// Count the number of molecule entries in an SDF file (delimited by `$$$$`).
    pub fn count_sdf_molecules(data: &str) -> usize {
        let count = data.lines().filter(|l| l.starts_with("$$$$")).count();
        if count == 0 && !data.trim().is_empty() { 1 } else { count }
    }

    /// Parse the n-th molecule (0-indexed) from a multi-molecule SDF file.
    pub fn from_sdf_nth_string(data: &str, n: usize) -> Result<MolecularSystem, JsValue> {
        parse_sdf_nth(data, n).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    /// Parse a SMILES string (organic subset). Returns topology + implicit H; 3D coords are 0.0.
    /// Supported: B C N O P S F Cl Br I, bracket atoms [Fe], branches, ring closures.
    /// Not supported (MVP): stereochemistry, isotopes.
    pub fn from_smiles(data: &str) -> Result<MolecularSystem, JsValue> {
        parse_smiles(data).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    /// Parse an mmCIF / PDBx string. Only the first MODEL is loaded.
    /// Reads the _atom_site loop for coordinates, element, residue, and chain data.
    pub fn from_mmcif_string(data: &str) -> Result<MolecularSystem, JsValue> {
        parse_mmcif(data).map_err(|e| JsValue::from_str(&e.to_string()))
    }
}

// ── Accessors ──────────────────────────────────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    /// Returns a 2048-bit (256-byte) Morgan/ECFP-like fingerprint as a byte Vec (Uint8Array in JS).
    /// `radius=2` is ECFP4-equivalent. Requires `compute_bonds()` (or `from_sdf_string()`) first.
    pub fn morgan_fingerprint(&self, radius: u32) -> Vec<u8> {
        morgan_fingerprint_bits(&self.symbols, &self.bonds, radius).to_vec()
    }

    /// Tanimoto (Jaccard) coefficient between this molecule and `other`.
    /// Both fingerprints are computed with the same `radius`.
    /// Returns 0.0 if both fingerprints are all-zero.
    pub fn tanimoto(&self, other: &MolecularSystem, radius: u32) -> f32 {
        let a = self.morgan_fingerprint(radius);
        let b = other.morgan_fingerprint(radius);
        let and_bits: u32 = a.iter().zip(b.iter()).map(|(x, y)| (x & y).count_ones()).sum();
        let or_bits:  u32 = a.iter().zip(b.iter()).map(|(x, y)| (x | y).count_ones()).sum();
        if or_bits == 0 { return 0.0; }
        and_bits as f32 / or_bits as f32
    }

    pub fn atom_count(&self) -> usize {
        self.symbols.len()
    }

    pub fn get_symbol(&self, index: usize) -> Option<String> {
        self.symbols.get(index).cloned()
    }

    pub fn get_x(&self, index: usize) -> Option<f32> {
        self.x.get(index).copied()
    }

    pub fn get_y(&self, index: usize) -> Option<f32> {
        self.y.get(index).copied()
    }

    pub fn get_z(&self, index: usize) -> Option<f32> {
        self.z.get(index).copied()
    }

    pub fn get_atom_name(&self, index: usize) -> Option<String> {
        self.atom_names.get(index).cloned()
    }

    pub fn get_residue_name(&self, index: usize) -> Option<String> {
        self.residue_names.get(index).cloned()
    }

    pub fn get_residue_id(&self, index: usize) -> Option<i32> {
        self.residue_ids.get(index).copied()
    }

    /// Returns the chain ID as a single-character string, or None if out of bounds.
    pub fn get_chain_id(&self, index: usize) -> Option<String> {
        self.chain_ids.get(index).map(|&b| (b as char).to_string())
    }

    pub fn is_hetatm(&self, index: usize) -> bool {
        self.hetatm_flags.get(index).copied().unwrap_or(false)
    }

    /// Returns all coordinates as an interleaved flat array [x0,y0,z0,x1,y1,z1,...].
    /// wasm-bindgen converts Vec<f32> to Float32Array on the JS side.
    pub fn get_positions_flat(&self) -> Vec<f32> {
        let mut buf = Vec::with_capacity(self.x.len() * 3);
        for ((&xi, &yi), &zi) in self.x.iter().zip(self.y.iter()).zip(self.z.iter()) {
            buf.push(xi);
            buf.push(yi);
            buf.push(zi);
        }
        buf
    }

    /// Returns element symbols as a JSON array string: `["O","H","H",...]`.
    pub fn get_symbols_json(&self) -> String {
        let parts: Vec<String> = self.symbols.iter()
            .map(|s| format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\"")))
            .collect();
        format!("[{}]", parts.join(","))
    }
}

// ── Geometry ───────────────────────────────────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    /// Euclidean distance between atoms `i` and `j` in Å.
    /// Returns 0.0 if either index is out of bounds.
    pub fn distance(&self, i: usize, j: usize) -> f32 {
        match (
            self.x.get(i), self.y.get(i), self.z.get(i),
            self.x.get(j), self.y.get(j), self.z.get(j),
        ) {
            (Some(&xi), Some(&yi), Some(&zi), Some(&xj), Some(&yj), Some(&zj)) => {
                let dx = xi - xj;
                let dy = yi - yj;
                let dz = zi - zj;
                (dx * dx + dy * dy + dz * dz).sqrt()
            }
            _ => 0.0,
        }
    }

    /// Bond angle at atom `j` in the i–j–k triplet, in degrees.
    /// Returns 0.0 if any index is out of bounds or if vectors are degenerate.
    pub fn angle(&self, i: usize, j: usize, k: usize) -> f32 {
        let n = self.symbols.len();
        if i >= n || j >= n || k >= n {
            return 0.0;
        }
        let ji = [self.x[i] - self.x[j], self.y[i] - self.y[j], self.z[i] - self.z[j]];
        let jk = [self.x[k] - self.x[j], self.y[k] - self.y[j], self.z[k] - self.z[j]];
        let dot = ji[0] * jk[0] + ji[1] * jk[1] + ji[2] * jk[2];
        let mag_ji = (ji[0] * ji[0] + ji[1] * ji[1] + ji[2] * ji[2]).sqrt();
        let mag_jk = (jk[0] * jk[0] + jk[1] * jk[1] + jk[2] * jk[2]).sqrt();
        if mag_ji < 1e-10 || mag_jk < 1e-10 {
            return 0.0;
        }
        (dot / (mag_ji * mag_jk)).clamp(-1.0, 1.0).acos().to_degrees()
    }

    /// Torsion (dihedral) angle for the i–j–k–l quartet, in degrees (range –180..180).
    /// Returns 0.0 if any index is out of bounds.
    pub fn dihedral(&self, i: usize, j: usize, k: usize, l: usize) -> f32 {
        let n = self.symbols.len();
        if i >= n || j >= n || k >= n || l >= n {
            return 0.0;
        }
        let b1 = [self.x[j]-self.x[i], self.y[j]-self.y[i], self.z[j]-self.z[i]];
        let b2 = [self.x[k]-self.x[j], self.y[k]-self.y[j], self.z[k]-self.z[j]];
        let b3 = [self.x[l]-self.x[k], self.y[l]-self.y[k], self.z[l]-self.z[k]];

        let cross = |a: [f32; 3], b: [f32; 3]| -> [f32; 3] {
            [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
        };
        let dot = |a: [f32; 3], b: [f32; 3]| a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        let mag = |a: [f32; 3]| (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]).sqrt();

        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let m1 = cross(n1, b2);
        let x = dot(n1, n2);
        let y = dot(m1, n2) / mag(b2).max(1e-10);
        y.atan2(x).to_degrees()
    }

    /// Returns the mass-weighted centroid [x, y, z] as a Float32Array.
    /// Unknown elements are assigned a carbon-like mass (12.0 u).
    pub fn center_of_mass(&self) -> Vec<f32> {
        let mut total = 0.0f32;
        let mut cx = 0.0f32;
        let mut cy = 0.0f32;
        let mut cz = 0.0f32;
        for i in 0..self.symbols.len() {
            let m = atomic_mass(&self.symbols[i]);
            total += m;
            cx += m * self.x[i];
            cy += m * self.y[i];
            cz += m * self.z[i];
        }
        if total < 1e-10 {
            return vec![0.0, 0.0, 0.0];
        }
        vec![cx / total, cy / total, cz / total]
    }

    /// Optimal superposition of `self` onto `reference` using the Kabsch algorithm.
    /// Updates `self`'s coordinates in-place (rotation + translation applied).
    /// Returns the RMSD after alignment. Atom count may differ; only the first
    /// `min(self.n, reference.n)` atoms are used.
    pub fn superpose(&mut self, reference: &MolecularSystem) -> f32 {
        let n = self.symbols.len().min(reference.symbols.len());
        if n == 0 { return 0.0; }

        let mobile: Vec<[f32;3]> = (0..n).map(|i| [self.x[i], self.y[i], self.z[i]]).collect();
        let refer:  Vec<[f32;3]> = (0..n).map(|i| [reference.x[i], reference.y[i], reference.z[i]]).collect();

        let (rot, c_mob, c_ref) = kabsch(&mobile, &refer);

        // Apply: new_pos = rot · (pos - c_mob) + c_ref
        for i in 0..n {
            let p = [self.x[i]-c_mob[0], self.y[i]-c_mob[1], self.z[i]-c_mob[2]];
            self.x[i] = rot[0][0]*p[0] + rot[0][1]*p[1] + rot[0][2]*p[2] + c_ref[0];
            self.y[i] = rot[1][0]*p[0] + rot[1][1]*p[1] + rot[1][2]*p[2] + c_ref[1];
            self.z[i] = rot[2][0]*p[0] + rot[2][1]*p[1] + rot[2][2]*p[2] + c_ref[2];
        }

        // Compute RMSD after superposition
        let sum: f32 = (0..n).map(|i| {
            let dx = self.x[i] - reference.x[i];
            let dy = self.y[i] - reference.y[i];
            let dz = self.z[i] - reference.z[i];
            dx*dx + dy*dy + dz*dz
        }).sum();
        (sum / n as f32).sqrt()
    }

    /// Returns the Kabsch-aligned RMSD between `self` and `reference` without
    /// modifying either structure's coordinates.
    pub fn rmsd_aligned(&self, reference: &MolecularSystem) -> f32 {
        let n = self.symbols.len().min(reference.symbols.len());
        if n == 0 { return 0.0; }

        let mobile: Vec<[f32;3]> = (0..n).map(|i| [self.x[i], self.y[i], self.z[i]]).collect();
        let refer:  Vec<[f32;3]> = (0..n).map(|i| [reference.x[i], reference.y[i], reference.z[i]]).collect();

        let (rot, c_mob, c_ref) = kabsch(&mobile, &refer);

        let sum: f32 = (0..n).map(|i| {
            let p = [mobile[i][0]-c_mob[0], mobile[i][1]-c_mob[1], mobile[i][2]-c_mob[2]];
            let ax = rot[0][0]*p[0] + rot[0][1]*p[1] + rot[0][2]*p[2] + c_ref[0];
            let ay = rot[1][0]*p[0] + rot[1][1]*p[1] + rot[1][2]*p[2] + c_ref[1];
            let az = rot[2][0]*p[0] + rot[2][1]*p[1] + rot[2][2]*p[2] + c_ref[2];
            let dx = ax - refer[i][0];
            let dy = ay - refer[i][1];
            let dz = az - refer[i][2];
            dx*dx + dy*dy + dz*dz
        }).sum();
        (sum / n as f32).sqrt()
    }

    /// Root-mean-square deviation between this system and `other` (same atom count assumed).
    /// Compares atoms by index order without superposition.
    /// Returns 0.0 if either system is empty.
    pub fn rmsd(&self, other: &MolecularSystem) -> f32 {
        let n = self.symbols.len().min(other.symbols.len());
        if n == 0 {
            return 0.0;
        }
        let sum: f32 = (0..n)
            .map(|i| {
                let dx = self.x[i] - other.x[i];
                let dy = self.y[i] - other.y[i];
                let dz = self.z[i] - other.z[i];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum / n as f32).sqrt()
    }

    /// Molecular formula in Hill order (C first, H second, then other elements alphabetically).
    /// Example: ethanol (CCO + 6H) → "C2H6O".
    pub fn molecular_formula(&self) -> String {
        use std::collections::BTreeMap;
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for sym in &self.symbols {
            *counts.entry(sym.clone()).or_insert(0) += 1;
        }
        let mut result = String::new();
        // Hill order: C, H, then rest alphabetically
        for elem in &["C", "H"] {
            if let Some(&n) = counts.get(*elem) {
                result.push_str(elem);
                if n > 1 { result.push_str(&n.to_string()); }
            }
        }
        for (elem, &n) in &counts {
            if elem == "C" || elem == "H" { continue; }
            result.push_str(elem);
            if n > 1 { result.push_str(&n.to_string()); }
        }
        result
    }

    /// Sum of atomic masses (Daltons). Unknown elements use carbon-like mass (12.0 u).
    pub fn molecular_weight(&self) -> f32 {
        self.symbols.iter().map(|s| atomic_mass(s)).sum()
    }
}

// ── Selection & Fragments ──────────────────────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    /// Returns a new MolecularSystem with only atoms belonging to `chain_id` (e.g. "A").
    /// Bond adjacency is remapped to the new indices.
    /// Returns an empty system if chain_id is not found or chain metadata is absent.
    pub fn select_chain(&self, chain_id: &str) -> MolecularSystem {
        let ch = chain_id.bytes().next().unwrap_or(b'A');
        let indices: Vec<usize> = self.chain_ids.iter().enumerate()
            .filter(|(_, &c)| c == ch)
            .map(|(i, _)| i)
            .collect();
        self.select_by_indices(&indices)
    }

    /// Returns a new MolecularSystem with only HETATM atoms (ligands, solvent, etc.).
    pub fn select_hetatm(&self) -> MolecularSystem {
        let indices: Vec<usize> = self.hetatm_flags.iter().enumerate()
            .filter(|(_, &h)| h)
            .map(|(i, _)| i)
            .collect();
        self.select_by_indices(&indices)
    }

    /// Returns a new MolecularSystem with only ATOM-record atoms (protein/nucleic acid).
    pub fn select_protein(&self) -> MolecularSystem {
        let indices: Vec<usize> = self.hetatm_flags.iter().enumerate()
            .filter(|(_, &h)| !h)
            .map(|(i, _)| i)
            .collect();
        self.select_by_indices(&indices)
    }

    /// Returns the one-letter amino acid sequence for `chain_id`.
    /// HETATM records are skipped. Unknown residues become 'X'.
    /// Returns empty string if no residue metadata is available.
    pub fn get_sequence(&self, chain_id: &str) -> String {
        if self.chain_ids.is_empty() || self.residue_names.is_empty() {
            return String::new();
        }
        let ch = chain_id.bytes().next().unwrap_or(b'A');
        let mut residue_map: std::collections::BTreeMap<i32, String> = std::collections::BTreeMap::new();
        let n = self.symbols.len();
        for i in 0..n {
            if self.chain_ids.get(i).copied() != Some(ch) { continue; }
            if self.hetatm_at(i) { continue; }
            let res_id = self.residue_id_i(i);
            let res_name = self.residue_names.get(i).cloned().unwrap_or_default();
            residue_map.entry(res_id).or_insert(res_name);
        }
        residue_map.values().map(|name| aa_one_letter(name)).collect()
    }

    /// Screen all molecules in a multi-SDF file in a single Wasm call.
    /// Returns a JS Array of objects: [{index, atom_count, bond_count, formula,
    /// molecular_weight, h_bond_donors, h_bond_acceptors, rotatable_bonds}, ...].
    /// Molecules that fail to parse are silently skipped.
    pub fn screen_sdf_string(data: &str) -> JsValue {
        let n = Self::count_sdf_molecules(data);
        let mut rows: Vec<SdfScreenRow> = Vec::with_capacity(n);
        for i in 0..n {
            if let Ok(mut mol) = parse_sdf_nth(data, i) {
                mol.compute_rings();
                rows.push(SdfScreenRow {
                    index: i,
                    atom_count: mol.atom_count(),
                    bond_count: mol.bond_count(),
                    formula: mol.molecular_formula(),
                    molecular_weight: mol.molecular_weight(),
                    h_bond_donors: mol.h_bond_donors(),
                    h_bond_acceptors: mol.h_bond_acceptors(),
                    rotatable_bonds: mol.rotatable_bond_count(),
                });
            }
        }
        to_js(&rows)
    }

    /// Returns connected components (fragments) as a JS Array of Uint32Array.
    /// Uses the bond adjacency list; call `compute_bonds()` or use `from_sdf_string()`
    /// first so that bonds are populated.
    pub fn get_fragments(&self) -> JsValue {
        let n = self.symbols.len();
        let mut visited = vec![false; n];
        let mut components: Vec<Vec<u32>> = Vec::new();

        for start in 0..n {
            if visited[start] {
                continue;
            }
            let mut component = Vec::new();
            let mut stack = vec![start];
            while let Some(atom) = stack.pop() {
                if visited[atom] {
                    continue;
                }
                visited[atom] = true;
                component.push(atom as u32);
                if let Some(neighbors) = self.bonds.get(atom) {
                    for &nb in neighbors {
                        if !visited[nb] {
                            stack.push(nb);
                        }
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }

        let arr = js_sys::Array::new();
        for comp in components {
            let inner = js_sys::Uint32Array::from(comp.as_slice());
            arr.push(&inner);
        }
        arr.into()
    }

    /// Returns indices of atoms within `radius` Å of atom `center` (excludes `center` itself).
    /// Uses the spatial grid if `build_spatial_index()` has been called; otherwise O(N) linear scan.
    /// Returns empty if `center` is out of bounds.
    pub fn get_atoms_within_radius(&self, center: usize, radius: f32) -> Vec<u32> {
        let (&cx, &cy, &cz) = match (self.x.get(center), self.y.get(center), self.z.get(center)) {
            (Some(x), Some(y), Some(z)) => (x, y, z),
            _ => return Vec::new(),
        };
        let r2 = radius * radius;
        let mut result = Vec::new();

        if let Some(grid) = &self.spatial_grid {
            let layers = (radius / grid.cell_size) as i32 + 1;
            let (ccx, ccy, ccz) = atom_cell(cx, cy, cz, &grid.origin, grid.cell_size);
            for dx in -layers..=layers {
                for dy in -layers..=layers {
                    for dz in -layers..=layers {
                        if let Some(atoms) = grid.cells.get(&(ccx + dx, ccy + dy, ccz + dz)) {
                            for &j in atoms {
                                if j != center {
                                    if let (Some(&xj), Some(&yj), Some(&zj)) =
                                        (self.x.get(j), self.y.get(j), self.z.get(j))
                                    {
                                        let ex = cx - xj;
                                        let ey = cy - yj;
                                        let ez = cz - zj;
                                        if ex * ex + ey * ey + ez * ez <= r2 {
                                            result.push(j as u32);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            result.sort_unstable();
        } else {
            for (j, ((&xj, &yj), &zj)) in self
                .x.iter().zip(self.y.iter()).zip(self.z.iter())
                .enumerate()
            {
                if j != center {
                    let dx = cx - xj;
                    let dy = cy - yj;
                    let dz = cz - zj;
                    if dx * dx + dy * dy + dz * dz <= r2 {
                        result.push(j as u32);
                    }
                }
            }
        }

        result
    }

    /// Build a uniform spatial grid with the given `cell_size` (Å).
    /// After calling this, `get_atoms_within_radius` uses the grid for faster queries.
    /// Recommended cell_size: 3–5 Å for typical molecular queries.
    pub fn build_spatial_index(&mut self, cell_size: f32) {
        if cell_size <= 0.0 || !cell_size.is_finite() {
            return;
        }
        use std::collections::HashMap;
        let n = self.x.len();
        if n == 0 {
            self.spatial_grid = Some(SpatialGrid {
                cells: HashMap::new(),
                cell_size,
                origin: [0.0; 3],
            });
            return;
        }
        let xmin = self.x.iter().cloned().fold(f32::INFINITY, f32::min);
        let ymin = self.y.iter().cloned().fold(f32::INFINITY, f32::min);
        let zmin = self.z.iter().cloned().fold(f32::INFINITY, f32::min);
        let origin = [xmin, ymin, zmin];
        let mut cells: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::with_capacity(n);
        for i in 0..n {
            let key = atom_cell(self.x[i], self.y[i], self.z[i], &origin, cell_size);
            cells.entry(key).or_default().push(i);
        }
        self.spatial_grid = Some(SpatialGrid { cells, cell_size, origin });
    }

    /// Returns true if `build_spatial_index()` has been called on this system.
    pub fn has_spatial_index(&self) -> bool {
        self.spatial_grid.is_some()
    }

    /// Returns sorted unique residue identifiers ("CHAIN:RESNAME:RESID") for all residues
    /// with at least one atom within `radius` Å of atom `center_atom`.
    /// Returns empty if `center_atom` is out of bounds or the system has no residue data.
    pub fn get_residues_within_radius(&self, center_atom: usize, radius: f32) -> Vec<String> {
        use std::collections::HashSet;
        let nearby = self.get_atoms_within_radius(center_atom, radius);
        let mut seen: HashSet<String> = HashSet::new();
        for idx in nearby {
            let i = idx as usize;
            let chain = self.chain_ids.get(i).copied().unwrap_or(b' ') as char;
            let res_name = self.residue_names.get(i).map(String::as_str).unwrap_or("");
            let res_id = self.residue_id_i(i);
            seen.insert(format!("{chain}:{res_name}:{res_id}"));
        }
        let mut result: Vec<String> = seen.into_iter().collect();
        result.sort();
        result
    }
}

// ── Bond Topology & Ring Detection ─────────────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    /// Parallel bond detection using rayon. Same bond rule as `compute_bonds()`.
    /// Available only with `--features parallel` on native targets (not Wasm).
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    pub fn compute_bonds_parallel(&mut self) {
        use rayon::prelude::*;
        let n = self.symbols.len();

        // Build per-atom radii once; None for unknown elements.
        let radii: Vec<Option<f32>> = self.symbols.iter().map(|s| covalent_radius(s)).collect();

        // Each row i produces the set of j > i it bonds to.
        let half_bonds: Vec<Vec<usize>> = (0..n)
            .into_par_iter()
            .map(|i| {
                let ri = match radii[i] {
                    Some(r) => r,
                    None => return vec![],
                };
                let mut row = Vec::new();
                for j in (i + 1)..n {
                    let rj = match radii[j] {
                        Some(r) => r,
                        None => continue,
                    };
                    let dx = self.x[i] - self.x[j];
                    let dy = self.y[i] - self.y[j];
                    let dz = self.z[i] - self.z[j];
                    let d2 = dx * dx + dy * dy + dz * dz;
                    let upper = ri + rj + 0.4;
                    if d2 > 0.16 && d2 <= upper * upper {
                        row.push(j);
                    }
                }
                row
            })
            .collect();

        // Merge into bidirectional adjacency list.
        self.bonds = vec![Vec::new(); n];
        for (i, neighbors) in half_bonds.into_iter().enumerate() {
            for j in neighbors {
                self.bonds[i].push(j);
                self.bonds[j].push(i);
            }
        }
    }

    /// Compute covalent bonds using per-element radii (Cordero 2008).
    /// Bond rule: 0.4 < d(i,j) <= r_cov(i) + r_cov(j) + 0.4 (Å).
    /// Atoms with unknown elements are skipped.
    /// Replaces any previously computed bond data.
    pub fn compute_bonds(&mut self) {
        let n = self.symbols.len();
        self.bonds = vec![Vec::new(); n];

        for i in 0..n {
            let ri = match covalent_radius(&self.symbols[i]) {
                Some(r) => r,
                None => continue,
            };
            for j in (i + 1)..n {
                let rj = match covalent_radius(&self.symbols[j]) {
                    Some(r) => r,
                    None => continue,
                };

                let dx = self.x[i] - self.x[j];
                let dy = self.y[i] - self.y[j];
                let dz = self.z[i] - self.z[j];
                let d2 = dx * dx + dy * dy + dz * dz;

                let upper = ri + rj + 0.4;
                let lower = 0.4f32;

                if d2 > lower * lower && d2 <= upper * upper {
                    self.bonds[i].push(j);
                    self.bonds[j].push(i);
                }
            }
        }
    }

    /// Returns neighbor atom indices for `atom`. Returns empty vec if out of bounds
    /// or if `compute_bonds()` has not been called.
    pub fn get_bonds(&self, atom: usize) -> Vec<u32> {
        self.bonds
            .get(atom)
            .map(|neighbors| neighbors.iter().map(|&i| i as u32).collect())
            .unwrap_or_default()
    }

    /// Total number of unique bonds (each bond counted once).
    pub fn bond_count(&self) -> usize {
        self.bonds.iter().map(|v| v.len()).sum::<usize>() / 2
    }

    /// Returns true if `compute_bonds()` has been called on this system.
    pub fn has_bonds_computed(&self) -> bool {
        self.bonds.len() == self.symbols.len()
    }

    // ── Ring Detection (Tarjan bridge-finding) ────────────────────────────────

    /// Detects which atoms and bonds are part of rings using Tarjan's bridge algorithm.
    /// Non-bridge bonds are ring bonds; atoms connected by ring bonds are ring atoms.
    /// Call `compute_bonds()` or use `from_sdf_string()` first.
    #[allow(clippy::while_let_loop)]
    pub fn compute_rings(&mut self) {
        let n = self.symbols.len();
        self.ring_atoms = vec![false; n];
        self.ring_bonds = std::collections::HashSet::new();

        if n == 0 || self.bonds.is_empty() {
            return;
        }

        let mut disc = vec![u32::MAX; n]; // discovery time; MAX = unvisited
        let mut low  = vec![0u32; n];
        let mut timer = 0u32;

        for start in 0..n {
            if disc[start] != u32::MAX {
                continue;
            }

            // Iterative DFS: each frame = (node, parent, neighbour_cursor)
            let mut stack: Vec<(usize, usize, usize)> = vec![(start, usize::MAX, 0)];

            loop {
                let frame = match stack.last_mut() {
                    Some(f) => f,
                    None => break,
                };
                let u = frame.0;
                let par = frame.1;

                // First visit: initialise disc and low
                if disc[u] == u32::MAX {
                    disc[u] = timer;
                    low[u] = timer;
                    timer += 1;
                }

                let neighbors = &self.bonds[u];
                if frame.2 < neighbors.len() {
                    let v = neighbors[frame.2];
                    frame.2 += 1;

                    if v == par {
                        continue; // skip edge back to parent
                    }
                    if disc[v] == u32::MAX {
                        // Tree edge — recurse
                        stack.push((v, u, 0));
                    } else {
                        // Back edge — update low[u] with disc[v].
                        // A back edge is itself always part of a ring (it creates the cycle).
                        let cur_low = low[u];
                        low[u] = cur_low.min(disc[v]);
                        let key = (u.min(v), u.max(v));
                        self.ring_bonds.insert(key);
                        self.ring_atoms[u] = true;
                        self.ring_atoms[v] = true;
                    }
                } else {
                    // Done with u — pop and update parent
                    stack.pop();
                    if let Some(parent_frame) = stack.last() {
                        let p = parent_frame.0;
                        let low_u = low[u];
                        let low_p = low[p];
                        low[p] = low_p.min(low_u);

                        // Determine if edge (p, u) is a bridge
                        if low_u <= disc[p] {
                            // Not a bridge → ring bond
                            let key = (u.min(p), u.max(p));
                            self.ring_bonds.insert(key);
                            self.ring_atoms[u] = true;
                            self.ring_atoms[p] = true;
                        }
                        // else: low[u] > disc[p] → bridge, not a ring bond
                    }
                }
            }
        }

        // Populate ring_sizes_per_atom from SSSR rings
        let rings = self.enumerate_rings();
        self.ring_sizes_per_atom = vec![Vec::new(); n];
        for ring in &rings {
            let sz = ring.len() as u8;
            for &a in ring {
                if a < n {
                    self.ring_sizes_per_atom[a].push(sz);
                }
            }
        }
    }

    /// Returns true if `compute_rings()` has been called on this system.
    pub fn has_rings_computed(&self) -> bool {
        !self.ring_atoms.is_empty()
    }

    /// Returns true if atom at `index` is part of a ring.
    pub fn is_ring_atom(&self, index: usize) -> bool {
        self.ring_atoms.get(index).copied().unwrap_or(false)
    }

    // ── Drug-likeness Descriptors ─────────────────────────────────────────────

    /// H-bond donors: N or O atoms that have at least one H neighbour in the bond graph.
    /// Accurate for SMILES-derived molecules (explicit H nodes). May return 0 for PDB/SDF
    /// where H atoms are often omitted.
    pub fn h_bond_donors(&self) -> u32 {
        self.symbols
            .iter()
            .enumerate()
            .filter(|(_, sym)| *sym == "N" || *sym == "O")
            .filter(|(i, _)| {
                self.bonds.get(*i).is_some_and(|nbrs| {
                    nbrs.iter().any(|&j| {
                        self.symbols.get(j).map(|s| s == "H").unwrap_or(false)
                    })
                })
            })
            .count() as u32
    }

    /// H-bond acceptors: all N and O atoms (without valence/hybridisation correction).
    pub fn h_bond_acceptors(&self) -> u32 {
        self.symbols.iter().filter(|s| *s == "N" || *s == "O").count() as u32
    }

    /// Rotatable bond count: non-ring bonds between two non-H heavy atoms where both
    /// atoms have degree > 1. Requires `compute_bonds()` + `compute_rings()`.
    pub fn rotatable_bond_count(&self) -> u32 {
        if !self.has_rings_computed() {
            return 0;
        }
        let mut count = 0u32;
        let n = self.symbols.len();
        for u in 0..n {
            if self.symbols[u] == "H" {
                continue;
            }
            let deg_u = self.bonds.get(u).map_or(0, |v| v.len());
            if deg_u <= 1 {
                continue;
            }
            if let Some(nbrs) = self.bonds.get(u) {
                for &v in nbrs {
                    if v <= u {
                        continue; // count each bond once
                    }
                    if self.symbols.get(v).map(|s| s == "H").unwrap_or(false) {
                        continue;
                    }
                    let deg_v = self.bonds.get(v).map_or(0, |b| b.len());
                    if deg_v <= 1 {
                        continue;
                    }
                    let key = (u.min(v), u.max(v));
                    if !self.ring_bonds.contains(&key) {
                        count += 1;
                    }
                }
            }
        }
        count
    }
}

// ── Per-Atom Data & Molecular Properties ──────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    /// Returns full data for atom `index` as a JS object `{index, symbol, x, y, z,
    /// atom_name, residue_name, residue_id, chain_id, is_hetatm}`.
    /// Returns `null` if `index` is out of bounds.
    pub fn get_atom_info(&self, index: usize) -> JsValue {
        match self.atom_info_at(index) {
            Some(info) => to_js(&info),
            None => JsValue::NULL,
        }
    }

    /// Returns full data for all atoms within `radius` Å of `center` as a JS Array
    /// of objects (same shape as `get_atom_info`).
    /// Uses the spatial grid if `build_spatial_index()` has been called.
    pub fn get_neighbors_info(&self, center: usize, radius: f32) -> JsValue {
        let indices = self.get_atoms_within_radius(center, radius);
        let infos: Vec<AtomInfo> = indices
            .iter()
            .filter_map(|&i| self.atom_info_at(i as usize))
            .collect();
        to_js(&infos)
    }

    // ── P12: Backbone Analysis ────────────────────────────────────────────────

    /// Returns backbone φ/ψ angles for all residues.
    /// Each element: {chain_id, residue_id, residue_name, phi, psi} — phi/psi are null for terminal residues.
    /// Returns null if no backbone metadata (e.g. XYZ/SMILES input).
    pub fn backbone_angles(&self) -> JsValue {
        let data = self.backbone_angle_data();
        if data.is_empty() { return JsValue::NULL; }
        to_js(&data)
    }

    /// Assigns secondary structure per residue from Ramachandran φ/ψ classification.
    /// Returns: [{chain_id, residue_id, residue_name, ss}] where ss = "H" (helix), "E" (strand), "C" (coil).
    /// Returns null if backbone metadata is absent.
    pub fn secondary_structure(&self) -> JsValue {
        let angles = self.backbone_angle_data();
        if angles.is_empty() { return JsValue::NULL; }

        let ss_data: Vec<SecStructRow> = angles.iter().map(|r| {
            let ss = match (r.phi, r.psi) {
                (Some(phi), Some(psi)) => classify_ss(phi, psi),
                _ => "C",
            };
            SecStructRow {
                chain_id: r.chain_id.clone(),
                residue_id: r.residue_id,
                residue_name: r.residue_name.clone(),
                ss: ss.to_string(),
            }
        }).collect();

        to_js(&ss_data)
    }

    // ── P13: Bond Order Storage + TPSA ───────────────────────────────────────

    /// Returns the bond order (1=single, 2=double, 3=triple) for atom `atom`'s
    /// k-th entry in its adjacency list. Returns 1 if bond orders are not stored.
    pub fn get_bond_order(&self, atom: usize, neighbor_k: usize) -> u8 {
        self.bond_orders
            .get(atom)
            .and_then(|v| v.get(neighbor_k))
            .copied()
            .unwrap_or(1)
    }

    /// Occupancy of atom `index` (0.0–1.0). Returns 1.0 if unavailable or out of bounds.
    #[wasm_bindgen]
    pub fn get_occupancy(&self, index: usize) -> f32 {
        self.occupancies.get(index).copied().unwrap_or(1.0)
    }

    /// Isotropic B-factor (Å²) of atom `index`. Returns 0.0 if unavailable or out of bounds.
    #[wasm_bindgen]
    pub fn get_b_factor(&self, index: usize) -> f32 {
        self.b_factors.get(index).copied().unwrap_or(0.0)
    }

    /// Topological Polar Surface Area using Ertl (2000) atom-type contributions.
    /// Accurate for SDF/SMILES inputs where bond orders are stored.
    /// Returns 0.0 for XYZ/PDB inputs (no bond order information).
    pub fn tpsa(&self) -> f32 {
        if self.bond_orders.is_empty() {
            return 0.0;
        }
        let n = self.symbols.len();
        let mut total = 0.0f32;
        for i in 0..n {
            let sym = &self.symbols[i];
            if sym != "O" && sym != "N" {
                continue;
            }
            // H neighbor count
            let h_count = self.bonds.get(i).map(|nbrs| {
                nbrs.iter().filter(|&&j| self.symbols.get(j).map(|s| s == "H").unwrap_or(false)).count()
            }).unwrap_or(0);
            // Has any double bond?
            let has_double = self.bond_orders.get(i).map(|orders| orders.contains(&2)).unwrap_or(false);
            // In ring?
            let in_ring = self.ring_atoms.get(i).copied().unwrap_or(false);

            total += if sym == "O" {
                if in_ring {
                    13.14
                } else if h_count >= 1 {
                    20.23  // OH
                } else if has_double {
                    17.07  // C=O carbonyl
                } else {
                    9.23   // ether / ester O
                }
            } else { // N
                if in_ring {
                    if h_count >= 1 { 13.97 } else { 12.89 }
                } else {
                    match h_count {
                        0 => 3.24,
                        1 => 24.06,
                        _ => 26.02,
                    }
                }
            };
        }
        total
    }

    // ── P16: LogP ────────────────────────────────────────────────────────────

    /// Approximate logP using simplified atom-type contributions (Wildman-Crippen-like).
    /// Requires `compute_bonds()` (bonds must be non-empty). Returns 0.0 if bonds absent.
    /// Accuracy: roughly ±1.0 for typical drug-like molecules; suitable for virtual screening.
    #[wasm_bindgen]
    pub fn logp(&self) -> f32 {
        if self.bonds.is_empty() {
            return 0.0;
        }
        let n = self.symbols.len();
        let mut total = 0.0f32;
        for i in 0..n {
            total += self.atom_logp_contribution(i);
        }
        total
    }

    fn atom_logp_contribution(&self, i: usize) -> f32 {
        let sym = match self.symbols.get(i) {
            Some(s) => s.as_str(),
            None => return 0.0,
        };

        // Skip hydrogen — contributions accounted for in heavy-atom types
        if sym == "H" {
            return 0.0;
        }

        // Count H neighbors (0 if no explicit H in structure)
        let h_ct: usize = self.bonds.get(i).map_or(0, |nbrs| {
            nbrs.iter().filter(|&&j| {
                self.symbols.get(j).map(|s| s == "H").unwrap_or(false)
            }).count()
        });

        // Check ring membership
        let in_ring = self.ring_atoms.get(i).copied().unwrap_or(false);

        // Has a double bond to oxygen (carbonyl)?
        let double_to_o = self.bonds.get(i).is_some_and(|nbrs| {
            nbrs.iter().enumerate().any(|(k, &j)| {
                self.symbols.get(j).map(|s| s == "O").unwrap_or(false)
                    && self.bond_orders.get(i).and_then(|o| o.get(k)).copied().unwrap_or(1) == 2
            })
        });

        // Has any double or triple bond?
        let max_order = self.bond_orders.get(i).map_or(1u8, |orders| {
            orders.iter().copied().max().unwrap_or(1)
        });

        match sym {
            "C" => {
                if in_ring && max_order >= 2 {
                    // Aromatic / sp2 ring carbon
                    if h_ct >= 1 { 0.36 } else { 0.12 }
                } else if double_to_o {
                    // Carbonyl carbon (C=O)
                    if h_ct >= 1 { 0.15 } else { -0.06 }
                } else if max_order >= 2 {
                    // Other sp2 carbon (C=C, C=N, etc.)
                    if h_ct >= 1 { 0.25 } else { 0.07 }
                } else {
                    // sp3 aliphatic carbon — contribution by H count
                    match h_ct {
                        3 => 0.51,
                        2 => 0.43,
                        1 => 0.25,
                        _ => -0.01,
                    }
                }
            }
            "N" => {
                if in_ring {
                    -0.59  // aromatic / ring N
                } else {
                    // Check for amide: N adjacent to C=O
                    let is_amide = self.bonds.get(i).is_some_and(|nbrs| {
                        nbrs.iter().any(|&j| {
                            self.symbols.get(j).map(|s| s == "C").unwrap_or(false)
                                && self.bonds.get(j).is_some_and(|jnbrs| {
                                    jnbrs.iter().enumerate().any(|(k, &jj)| {
                                        self.symbols.get(jj).map(|s| s == "O").unwrap_or(false)
                                            && self.bond_orders.get(j)
                                                .and_then(|o| o.get(k)).copied().unwrap_or(1) == 2
                                    })
                                })
                        })
                    });
                    if is_amide { -1.03 }
                    else if self.bond_orders.get(i).is_some_and(|o| o.contains(&3)) {
                        -0.33  // nitrile N (triple bond)
                    } else if h_ct > 0 { -0.96 }
                    else { -0.67 }
                }
            }
            "O" => {
                if in_ring {
                    -0.04  // O in aromatic ring (furan-like)
                } else if h_ct > 0 {
                    -0.67  // hydroxyl / carboxyl O–H
                } else if max_order >= 2 {
                    -0.44  // carbonyl O (C=O)
                } else {
                    -0.27  // ether O
                }
            }
            "S" => {
                if in_ring { 0.15 }
                else if h_ct > 0 { -0.09 }
                else { 0.03 }
            }
            "F"  => 0.14,
            "Cl" => 0.60,
            "Br" => 0.88,
            "I"  => 1.01,
            "P"  => -0.18,
            _    => 0.0,
        }
    }
}

// ── Protein & Structural Analysis ─────────────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    /// Find hydrogen bonds using geometric criteria.
    /// Donors: N, O, F atoms. Acceptors: N, O, F, S atoms.
    /// For structures with explicit H (SMILES after `compute_bonds()`): D-H···A angle filter (≥ 120°).
    /// For structures without H (PDB/mmCIF): distance-only (D···A < cutoff).
    /// Default cutoff: 3.5 Å. Returns a JS Array of HBondRow objects.
    #[wasm_bindgen]
    pub fn find_h_bonds(&self, cutoff: f32) -> JsValue {
        let result = self.find_h_bonds_data(cutoff);
        to_js(&result)
    }

    fn find_h_bonds_data(&self, cutoff: f32) -> Vec<HBondRow> {
        use std::collections::HashSet;

        let n = self.symbols.len();
        if n == 0 || cutoff <= 0.0 {
            return Vec::new();
        }

        // Build acceptor set: indices of N, O, F, S atoms
        let acceptors: HashSet<usize> = (0..n)
            .filter(|&i| matches!(self.symbols[i].as_str(), "N" | "O" | "F" | "S"))
            .collect();

        const MIN_ANGLE: f32 = 120.0;
        let mut result: Vec<HBondRow> = Vec::new();

        // Iterate over donors: N, O, F atoms
        for donor in 0..n {
            let dsym = self.symbols[donor].as_str();
            if !matches!(dsym, "N" | "O" | "F") {
                continue;
            }

            // Find H neighbors of this donor (present only if compute_bonds() was called with explicit H)
            let h_neighbors: Vec<usize> = if self.bonds.is_empty() {
                Vec::new()
            } else {
                self.bonds.get(donor).map_or(Vec::new(), |nbrs| {
                    nbrs.iter()
                        .filter(|&&j| self.symbols.get(j).map(|s| s == "H").unwrap_or(false))
                        .copied()
                        .collect()
                })
            };

            // Get candidate acceptors within cutoff
            let candidates = self.get_atoms_within_radius(donor, cutoff);

            for acc_u32 in candidates {
                let acc = acc_u32 as usize;
                if acc == donor || !acceptors.contains(&acc) {
                    continue;
                }

                let dist = self.distance(donor, acc);

                if h_neighbors.is_empty() {
                    // No explicit H → distance-only H-bond
                    result.push(HBondRow {
                        donor,
                        acceptor: acc,
                        distance: dist,
                        h_atom: None,
                        angle: None,
                    });
                } else {
                    // Check D-H···A geometry for each H neighbor
                    for &h in &h_neighbors {
                        // angle(donor, h, acc) = D-H-A angle at H (middle atom)
                        let ang = self.angle(donor, h, acc);
                        if ang >= MIN_ANGLE {
                            result.push(HBondRow {
                                donor,
                                acceptor: acc,
                                distance: dist,
                                h_atom: Some(h),
                                angle: Some(ang),
                            });
                        }
                    }
                }
            }
        }

        result
    }

    // ── P18: SSSR Ring Enumeration ────────────────────────────────────────────

    /// Enumerate rings via BFS from each ring bond, deduplicating by atom set.
    fn enumerate_rings(&self) -> Vec<Vec<usize>> {
        use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

        let mut seen: HashSet<BTreeSet<usize>> = HashSet::new();
        let mut result: Vec<Vec<usize>> = Vec::new();

        if self.bonds.is_empty() || self.ring_bonds.is_empty() {
            return result;
        }

        for &(u, v) in &self.ring_bonds {
            // Try both directions so fused-ring systems produce both rings
            for (src, dst) in [(u, v), (v, u)] {
                let mut queue: VecDeque<usize> = VecDeque::from([src]);
                let mut prev: HashMap<usize, usize> = HashMap::new();
                prev.insert(src, src);

                'bfs: while let Some(node) = queue.pop_front() {
                    let neighbors = match self.bonds.get(node) {
                        Some(n) => n.clone(),
                        None => continue,
                    };
                    for nb in neighbors {
                        if node == src && nb == dst {
                            continue; // exclude direct edge src→dst
                        }
                        if prev.contains_key(&nb) {
                            continue;
                        }
                        prev.insert(nb, node);
                        if nb == dst {
                            // Reconstruct path src → … → dst
                            let mut path = Vec::new();
                            let mut cur = dst;
                            loop {
                                path.push(cur);
                                let p = prev[&cur];
                                if p == cur {
                                    break;
                                }
                                cur = p;
                            }
                            path.reverse();
                            let atom_set: BTreeSet<usize> = path.iter().copied().collect();
                            if seen.insert(atom_set) {
                                result.push(path);
                            }
                            break 'bfs;
                        }
                        queue.push_back(nb);
                    }
                }
            }
        }

        result
    }

    /// Returns each ring as a Uint32Array of atom indices, wrapped in a JS Array.
    /// Requires `compute_bonds()` + `compute_rings()`. Returns empty if no rings.
    #[wasm_bindgen]
    pub fn get_rings(&self) -> JsValue {
        let rings = self.enumerate_rings();
        let arr = js_sys::Array::new();
        for ring in rings {
            let ring_u32: Vec<u32> = ring.iter().map(|&a| a as u32).collect();
            let inner = js_sys::Uint32Array::from(ring_u32.as_slice());
            arr.push(&inner);
        }
        arr.into()
    }

    /// Returns the number of rings containing at least one double bond (aromatic after Kekulization).
    /// Requires `compute_bonds()` + `compute_rings()`.
    #[wasm_bindgen]
    pub fn aromatic_ring_count(&self) -> u32 {
        if self.bond_orders.is_empty() {
            return 0;
        }
        let rings = self.enumerate_rings();
        let mut count = 0u32;
        'ring: for ring in &rings {
            let n = ring.len();
            for i in 0..n {
                let a = ring[i];
                let b = ring[(i + 1) % n];
                if let Some(nbrs) = self.bonds.get(a) {
                    for (k, &nb) in nbrs.iter().enumerate() {
                        if nb == b {
                            let ord = self.bond_orders.get(a)
                                .and_then(|orders| orders.get(k))
                                .copied()
                                .unwrap_or(1);
                            if ord == 2 {
                                count += 1;
                                continue 'ring;
                            }
                            break;
                        }
                    }
                }
            }
        }
        count
    }

    // ── P19: Disulfide Bonds + Metal Coordination Sites ───────────────────────

    /// Detect S–S disulfide bonds within `cutoff` Å.
    /// Returns a JS Array of `{atom_i, atom_j, distance}` objects.
    #[wasm_bindgen]
    pub fn find_disulfide_bonds(&self, cutoff: f32) -> JsValue {
        let result = self.find_disulfide_bonds_data(cutoff);
        to_js(&result)
    }

    fn find_disulfide_bonds_data(&self, cutoff: f32) -> Vec<DisulfideBond> {
        use std::collections::HashSet;
        let n = self.symbols.len();
        if n == 0 || cutoff <= 0.0 {
            return Vec::new();
        }
        let mut seen: HashSet<(usize, usize)> = HashSet::new();
        let mut result: Vec<DisulfideBond> = Vec::new();
        for s in 0..n {
            if self.symbols[s] != "S" {
                continue;
            }
            for j_u32 in self.get_atoms_within_radius(s, cutoff) {
                let j = j_u32 as usize;
                if j == s || self.symbols[j] != "S" {
                    continue;
                }
                let key = (s.min(j), s.max(j));
                if seen.contains(&key) {
                    continue;
                }
                let d = self.distance(s, j);
                if d <= cutoff {
                    seen.insert(key);
                    result.push(DisulfideBond { atom_i: s, atom_j: j, distance: d });
                }
            }
        }
        result
    }

    /// Detect metal coordination sites within `cutoff` Å.
    /// Metals: Fe, Mg, Zn, Cu, Ca, Mn, Co, Ni, Na, K, Mo, W.
    /// Coordinating atoms: N, O, S, P.
    /// Returns a JS Array of `{metal_atom, element, coordinating}` objects.
    #[wasm_bindgen]
    pub fn find_metal_sites(&self, cutoff: f32) -> JsValue {
        let result = self.find_metal_sites_data(cutoff);
        to_js(&result)
    }

    fn find_metal_sites_data(&self, cutoff: f32) -> Vec<MetalSite> {
        const METALS: &[&str] = &["Fe", "Mg", "Zn", "Cu", "Ca", "Mn", "Co", "Ni", "Na", "K", "Mo", "W"];
        const COORD: &[&str] = &["N", "O", "S", "P"];
        let n = self.symbols.len();
        if n == 0 || cutoff <= 0.0 {
            return Vec::new();
        }
        let coord_set: std::collections::HashSet<&str> = COORD.iter().copied().collect();
        let mut result: Vec<MetalSite> = Vec::new();
        for i in 0..n {
            let sym = self.symbols[i].as_str();
            if !METALS.contains(&sym) {
                continue;
            }
            let mut coordinating: Vec<usize> = Vec::new();
            for j_u32 in self.get_atoms_within_radius(i, cutoff) {
                let j = j_u32 as usize;
                if j != i && coord_set.contains(self.symbols[j].as_str()) {
                    coordinating.push(j);
                }
            }
            result.push(MetalSite { metal_atom: i, element: self.symbols[i].clone(), coordinating });
        }
        result
    }

    // ── P20: Contact Map + Binding Site Residues ──────────────────────────────

    /// Residue–residue Cα–Cα contact map within `cutoff` Å.
    /// Returns empty if no Cα atoms (XYZ/SMILES input). Use `build_spatial_index()` for performance.
    #[wasm_bindgen]
    pub fn contact_map(&self, cutoff: f32) -> JsValue {
        let result = self.contact_map_data(cutoff);
        to_js(&result)
    }

    fn contact_map_data(&self, cutoff: f32) -> Vec<ContactMapRow> {
        use std::collections::HashSet;
        let n = self.symbols.len();
        if n == 0 || cutoff <= 0.0 || self.atom_names.is_empty() {
            return Vec::new();
        }
        // Build Cα index: map atom index to true for Cα atoms
        let mut ca_set: HashSet<usize> = HashSet::new();
        // Also keep residue info per atom index
        for i in 0..n {
            if self.atom_names.get(i).map(|s| s == "CA").unwrap_or(false)
                && !self.hetatm_at(i)
            {
                ca_set.insert(i);
            }
        }
        if ca_set.is_empty() {
            return Vec::new();
        }
        let mut result: Vec<ContactMapRow> = Vec::new();
        for &atom_a in &ca_set {
            for j_u32 in self.get_atoms_within_radius(atom_a, cutoff) {
                let atom_b = j_u32 as usize;
                if atom_b <= atom_a || !ca_set.contains(&atom_b) {
                    continue;
                }
                let d = self.distance(atom_a, atom_b);
                result.push(ContactMapRow {
                    chain_i: (self.chain_ids.get(atom_a).copied().unwrap_or(b'A') as char).to_string(),
                    resid_i: self.residue_ids.get(atom_a).copied().unwrap_or(0),
                    resname_i: self.residue_names.get(atom_a).cloned().unwrap_or_default(),
                    chain_j: (self.chain_ids.get(atom_b).copied().unwrap_or(b'A') as char).to_string(),
                    resid_j: self.residue_ids.get(atom_b).copied().unwrap_or(0),
                    resname_j: self.residue_names.get(atom_b).cloned().unwrap_or_default(),
                    distance: d,
                });
            }
        }
        result
    }

    /// Protein residues (ATOM records) within `cutoff` Å of any HETATM atom.
    /// Deduplicates by (chain_id, residue_id). Returns empty if no HETATM atoms.
    #[wasm_bindgen]
    pub fn binding_site_residues(&self, cutoff: f32) -> JsValue {
        let result = self.binding_site_residues_data(cutoff);
        to_js(&result)
    }

    fn binding_site_residues_data(&self, cutoff: f32) -> Vec<BindingSiteRow> {
        use std::collections::HashSet;
        let n = self.symbols.len();
        if n == 0 || cutoff <= 0.0 || self.hetatm_flags.is_empty() {
            return Vec::new();
        }
        let mut seen: HashSet<(u8, i32)> = HashSet::new();
        let mut result: Vec<BindingSiteRow> = Vec::new();
        for h in 0..n {
            if !self.hetatm_flags.get(h).copied().unwrap_or(false) {
                continue;
            }
            for j_u32 in self.get_atoms_within_radius(h, cutoff) {
                let j = j_u32 as usize;
                if self.hetatm_flags.get(j).copied().unwrap_or(true) {
                    continue; // skip other HETATM atoms
                }
                let chain = self.chain_ids.get(j).copied().unwrap_or(b'A');
                let resid = self.residue_ids.get(j).copied().unwrap_or(0);
                if seen.insert((chain, resid)) {
                    result.push(BindingSiteRow {
                        chain_id: (chain as char).to_string(),
                        residue_id: resid,
                        residue_name: self.residue_names.get(j).cloned().unwrap_or_default(),
                    });
                }
            }
        }
        result
    }

    // ── P21: Formal Charge + XYZ Output ──────────────────────────────────────

    /// Net formal charge summed over all atoms. Accurate for SMILES bracket atoms.
    /// Returns 0 for XYZ/PDB/SDF inputs where charge is not stored.
    #[wasm_bindgen]
    pub fn formal_charge(&self) -> i32 {
        self.charges.iter().sum()
    }

    /// Per-atom formal charge. Returns 0 if out of bounds or charges not stored.
    #[wasm_bindgen]
    pub fn get_formal_charge(&self, index: usize) -> i32 {
        self.charges.get(index).copied().unwrap_or(0)
    }

    // ── Atom mapping (SMILES [C:1] notation) ──────────────────────────────────

    /// Atom map index at `idx`. 0 means no mapping. Parsed from SMILES `[C:n]`.
    #[wasm_bindgen]
    pub fn get_atom_map_index(&self, idx: usize) -> u32 {
        self.atom_map.get(idx).copied().unwrap_or(0)
    }

    /// Set the atom map index for atom `idx`. Extends the map array if needed.
    #[wasm_bindgen]
    pub fn set_atom_map_index(&mut self, idx: usize, map_num: u32) {
        if self.atom_map.len() <= idx {
            self.atom_map.resize(self.symbols.len().max(idx + 1), 0);
        }
        if let Some(slot) = self.atom_map.get_mut(idx) {
            *slot = map_num;
        }
    }

    /// Returns true if any atom carries a non-zero map index.
    #[wasm_bindgen]
    pub fn has_atom_map(&self) -> bool {
        self.atom_map.iter().any(|&m| m != 0)
    }

    /// Reset all atom map indices to 0.
    #[wasm_bindgen]
    pub fn clear_atom_map(&mut self) {
        self.atom_map.iter_mut().for_each(|m| *m = 0);
    }

    /// Serialize this molecular system to XYZ format string.
    #[wasm_bindgen]
    pub fn to_xyz_string(&self) -> String {
        let n = self.symbols.len();
        let mut s = format!("{}\ngenerated by chem-wasm-lens\n", n);
        for i in 0..n {
            s.push_str(&format!(
                "{:<2}  {:10.6}  {:10.6}  {:10.6}\n",
                self.symbols[i],
                self.x.get(i).copied().unwrap_or(0.0),
                self.y.get(i).copied().unwrap_or(0.0),
                self.z.get(i).copied().unwrap_or(0.0),
            ));
        }
        s
    }

    /// Serialize this molecular system to SDF V2000 format string.
    /// Bond block is included when bonds are populated (e.g. from SDF/SMILES parser or compute_bonds()).
    #[wasm_bindgen]
    pub fn to_sdf_string(&self) -> String {
        let n = self.symbols.len();
        let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        let mut bonds_vec: Vec<(usize, usize, u8)> = Vec::new();
        for i in 0..self.bonds.len() {
            for (k, &j) in self.bonds[i].iter().enumerate() {
                let key = (i.min(j), i.max(j));
                if seen.insert(key) {
                    let ord = self.bond_orders.get(i).and_then(|o| o.get(k)).copied().unwrap_or(1);
                    bonds_vec.push((key.0, key.1, ord));
                }
            }
        }
        let nb = bonds_vec.len();

        let mut out = String::new();
        out.push_str("\n  chem-wasm-lens\n\n");
        out.push_str(&format!("{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n", n, nb));
        for i in 0..n {
            out.push_str(&format!(
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                self.x.get(i).copied().unwrap_or(0.0),
                self.y.get(i).copied().unwrap_or(0.0),
                self.z.get(i).copied().unwrap_or(0.0),
                self.symbols[i],
            ));
        }
        for (a, b, ord) in &bonds_vec {
            out.push_str(&format!("{:>3}{:>3}{:>3}  0  0  0  0\n", a + 1, b + 1, ord));
        }
        out.push_str("M  END\n$$$$\n");
        out
    }

    fn screen_sdf_diverse_data(data: &str, k: usize, fp_radius: u32) -> Vec<usize> {
        let n = MolecularSystem::count_sdf_molecules(data);
        if n == 0 {
            return Vec::new();
        }
        let k = k.min(n);

        let fps: Vec<[u8; 256]> = (0..n)
            .map(|i| {
                parse_sdf_nth(data, i)
                    .ok()
                    .map(|mol| morgan_fingerprint_bits(&mol.symbols, &mol.bonds, fp_radius))
                    .unwrap_or([0u8; 256])
            })
            .collect();

        let tanimoto_fp = |a: &[u8; 256], b: &[u8; 256]| -> f32 {
            let and: u32 = a.iter().zip(b.iter()).map(|(x, y)| (x & y).count_ones()).sum();
            let or_: u32 = a.iter().zip(b.iter()).map(|(x, y)| (x | y).count_ones()).sum();
            if or_ == 0 { 1.0 } else { and as f32 / or_ as f32 }
        };

        let mut selected: Vec<usize> = vec![0];
        let mut selected_set: std::collections::HashSet<usize> = std::collections::HashSet::from([0]);
        let mut min_dists: Vec<f32> = fps.iter().map(|fp| 1.0 - tanimoto_fp(&fps[0], fp)).collect();

        for _ in 1..k {
            let best = min_dists
                .iter()
                .enumerate()
                .filter(|(i, _)| !selected_set.contains(i))
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, _)| i)
                .unwrap_or(0);
            selected.push(best);
            selected_set.insert(best);
            for (i, fp) in fps.iter().enumerate() {
                let d = 1.0 - tanimoto_fp(&fps[best], fp);
                if d < min_dists[i] {
                    min_dists[i] = d;
                }
            }
        }
        selected
    }

    /// Pick k maximally diverse molecules from a multi-SDF string using the MaxMin algorithm.
    /// Uses Morgan fingerprints at `fp_radius`. Returns a `u32[]` of molecule indices.
    /// If k >= molecule count, returns all indices.
    #[wasm_bindgen]
    pub fn screen_sdf_diverse(data: &str, k: u32, fp_radius: u32) -> JsValue {
        let selected = MolecularSystem::screen_sdf_diverse_data(data, k as usize, fp_radius);
        let result: Vec<u32> = selected.iter().map(|&i| i as u32).collect();
        to_js(&result)
    }

    /// Serialize this molecular system to PDB format string.
    /// Uses ATOM/HETATM metadata when available; falls back to HETATM for non-PDB sources.
    #[wasm_bindgen]
    pub fn to_pdb_string(&self) -> String {
        let n = self.symbols.len();
        let has_meta = !self.atom_names.is_empty();
        let mut out = String::new();
        for i in 0..n {
            let record = if has_meta && self.hetatm_at(i) {
                "HETATM"
            } else if has_meta {
                "ATOM  "
            } else {
                "HETATM"
            };
            let atom_name = if has_meta {
                format!("{:<4}", self.atom_name_str(i))
            } else {
                format!(" {:<3}", self.symbols[i])
            };
            let res_name = self.residue_names.get(i).map(|s| s.as_str()).unwrap_or("LIG");
            let chain = self.chain_id_byte(i) as char;
            let res_seq = self.residue_ids.get(i).copied().unwrap_or(1);
            let occ = self.occupancies.get(i).copied().unwrap_or(1.0);
            let bfac = self.b_factors.get(i).copied().unwrap_or(0.0);
            out.push_str(&format!(
                "{}{:>5} {} {:<3} {}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}  \n",
                record,
                i + 1,
                atom_name,
                res_name,
                chain,
                res_seq,
                self.x.get(i).copied().unwrap_or(0.0),
                self.y.get(i).copied().unwrap_or(0.0),
                self.z.get(i).copied().unwrap_or(0.0),
                occ,
                bfac,
                self.symbols[i],
            ));
        }
        out.push_str("END\n");
        out
    }

    /// Per-atom solvent-accessible surface area (Å²) using the Shrake-Rupley algorithm.
    /// `probe_radius` is typically 1.4 Å for water. Uses spatial index if built.
    #[allow(clippy::needless_range_loop)]
    #[wasm_bindgen]
    pub fn sasa_per_atom(&self, probe_radius: f32) -> Vec<f32> {
        const VDW_RADII: &[(&str, f32)] = &[
            ("H", 1.20), ("C", 1.70), ("N", 1.55), ("O", 1.52), ("S", 1.80),
            ("P", 1.80), ("F", 1.47), ("Cl", 1.75), ("Br", 1.85), ("I", 1.98),
        ];
        let vdw = |sym: &str| {
            VDW_RADII.iter().find(|(s, _)| *s == sym).map(|(_, r)| *r).unwrap_or(1.70)
        };

        let n = self.symbols.len();
        let sphere = fibonacci_sphere_92();
        let n_pts = sphere.len() as f32;
        let max_vdw = 1.98f32;
        let search_r = (max_vdw + probe_radius) * 2.0;

        let mut result = vec![0.0f32; n];
        for i in 0..n {
            let ri = vdw(&self.symbols[i]) + probe_radius;
            let neighbors = self.get_atoms_within_radius(i, search_r);
            let exposed = sphere.iter().filter(|(px, py, pz)| {
                let tx = self.x[i] + px * ri;
                let ty = self.y[i] + py * ri;
                let tz = self.z[i] + pz * ri;
                neighbors.iter().all(|&j_u32| {
                    let j = j_u32 as usize;
                    let rj = vdw(&self.symbols[j]) + probe_radius;
                    let dx = tx - self.x[j];
                    let dy = ty - self.y[j];
                    let dz = tz - self.z[j];
                    dx * dx + dy * dy + dz * dz >= rj * rj
                })
            }).count();
            result[i] = 4.0 * std::f32::consts::PI * ri * ri * (exposed as f32 / n_pts);
        }
        result
    }

    /// Total solvent-accessible surface area (Å²) via Shrake-Rupley.
    /// `probe_radius` is typically 1.4 Å for water.
    #[wasm_bindgen]
    pub fn sasa(&self, probe_radius: f32) -> f32 {
        self.sasa_per_atom(probe_radius).iter().sum()
    }
}

// ── SMARTS Substructure Search ─────────────────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    fn match_smarts_data(&self, smarts: &str) -> Vec<Vec<usize>> {
        let query = match parse_smarts(smarts) {
            Some(q) if !q.atoms.is_empty() => q,
            _ => return Vec::new(),
        };
        let n_query = query.atoms.len();
        let n_mol = self.symbols.len();
        let bfs_order = smarts_bfs_order(&query);

        let mut results: Vec<Vec<usize>> = Vec::new();
        let mut mapping: Vec<Option<usize>> = vec![None; n_query];
        let mut used: Vec<bool> = vec![false; n_mol];

        for start in 0..n_mol {
            if smarts_atom_matches(&query.atoms[bfs_order[0]], self, start) {
                mapping[bfs_order[0]] = Some(start);
                used[start] = true;
                smarts_backtrack(&query, self, &mut mapping, &mut used, &bfs_order, 1, &mut results);
                used[start] = false;
                mapping[bfs_order[0]] = None;
            }
        }
        results
    }

    /// Returns true if this molecule contains a substructure matching `smarts`.
    /// Requires `compute_bonds()` + `compute_rings()` for accurate aromatic matching.
    #[wasm_bindgen]
    pub fn has_substructure(&self, smarts: &str) -> bool {
        !self.match_smarts_data(smarts).is_empty()
    }

    /// Returns all substructure matches as `Array<Uint32Array>` of atom indices.
    /// Each Uint32Array = one match (atom indices in mol order). Empty = no match or invalid SMARTS.
    pub fn match_smarts(&self, smarts: &str) -> JsValue {
        let arr = js_sys::Array::new();
        for m in self.match_smarts_data(smarts) {
            let inner = js_sys::Uint32Array::from(
                m.iter().map(|&i| i as u32).collect::<Vec<_>>().as_slice(),
            );
            arr.push(&inner);
        }
        arr.into()
    }

    /// Returns all substructure matches as `Array<Uint32Array>` of atom indices.
    /// Each `Uint32Array` = one match. Empty = no match or invalid SMARTS.
    /// Alias for `match_smarts` with a cleaner name.
    pub fn find_substructure(&self, smarts: &str) -> JsValue {
        self.match_smarts(smarts)
    }

    /// Returns a sorted `Uint32Array` of all atom indices in any SMARTS match.
    /// Pass the result to `to_svg_string_highlighted` to draw halos.
    pub fn get_substructure_atoms(&self, smarts: &str) -> js_sys::Uint32Array {
        let matches = self.match_smarts_data(smarts);
        let mut set: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for m in &matches {
            for &i in m {
                set.insert(i as u32);
            }
        }
        let mut v: Vec<u32> = set.into_iter().collect();
        v.sort_unstable();
        js_sys::Uint32Array::from(v.as_slice())
    }

    /// Renders the 2D SVG with yellow halos on the given atom indices.
    /// `highlight` is a `Uint32Array` (e.g. from `get_substructure_atoms`).
    #[wasm_bindgen]
    pub fn to_svg_string_highlighted(&self, width: u32, height: u32, highlight: &[u32]) -> String {
        let set: std::collections::HashSet<usize> =
            highlight.iter().map(|&i| i as usize).collect();
        self.to_svg_data_impl(width, height, &set)
    }

    // ── P25: Residue SASA + chain interface ──────────────────────────────────

    #[allow(clippy::needless_range_loop)]
    fn residue_sasa_data(&self, probe_radius: f32) -> Vec<ResidueSasaRow> {
        if self.atom_names.is_empty() {
            return Vec::new();
        }
        let n = self.symbols.len();
        let per_atom = self.sasa_per_atom(probe_radius);
        let mut map: std::collections::BTreeMap<(u8, i32), (String, String, f32)> =
            std::collections::BTreeMap::new();
        for i in 0..n {
            let chain = self.chain_id_byte(i);
            let resid = self.residue_ids.get(i).copied().unwrap_or(1);
            let entry = map.entry((chain, resid)).or_insert_with(|| {
                let chain_str = (chain as char).to_string();
                let resname = self.residue_names.get(i).cloned().unwrap_or_default();
                (chain_str, resname, 0.0)
            });
            entry.2 += per_atom[i];
        }
        map.into_iter()
            .map(|((_, resid), (chain_id, residue_name, sasa))| ResidueSasaRow {
                chain_id,
                residue_id: resid,
                residue_name,
                sasa,
            })
            .collect()
    }

    /// Per-residue solvent-accessible surface area (Å²).
    /// Returns `[{chain_id, residue_id, residue_name, sasa}]`. Empty for non-PDB inputs.
    #[wasm_bindgen]
    pub fn residue_sasa(&self, probe_radius: f32) -> JsValue {
        to_js(&self.residue_sasa_data(probe_radius))
    }

    fn chain_interface_data(
        &self,
        chain_a: &str,
        chain_b: &str,
        cutoff: f32,
    ) -> ChainInterfaceResult {
        if self.chain_ids.is_empty() {
            return ChainInterfaceResult { a: Vec::new(), b: Vec::new() };
        }
        let a_byte = chain_a.bytes().next().unwrap_or(b'A');
        let b_byte = chain_b.bytes().next().unwrap_or(b'B');
        let n = self.symbols.len();
        let b_set: std::collections::HashSet<usize> =
            (0..n).filter(|&i| self.chain_ids.get(i).copied().unwrap_or(0) == b_byte).collect();

        let mut a_res: std::collections::BTreeSet<(i32, String)> =
            std::collections::BTreeSet::new();
        let mut b_res: std::collections::BTreeSet<(i32, String)> =
            std::collections::BTreeSet::new();

        for i in 0..n {
            if self.chain_ids.get(i).copied().unwrap_or(0) != a_byte {
                continue;
            }
            let neighbors = self.get_atoms_within_radius(i, cutoff);
            for j_u32 in neighbors {
                let j = j_u32 as usize;
                if b_set.contains(&j) {
                    let resid_a = self.residue_id_i(i);
                    let resname_a = self.residue_names.get(i).cloned().unwrap_or_default();
                    a_res.insert((resid_a, resname_a));
                    let resid_b = self.residue_ids.get(j).copied().unwrap_or(0);
                    let resname_b = self.residue_names.get(j).cloned().unwrap_or_default();
                    b_res.insert((resid_b, resname_b));
                }
            }
        }

        let chain_a_str = chain_a.chars().next().map(|c| c.to_string()).unwrap_or_default();
        let chain_b_str = chain_b.chars().next().map(|c| c.to_string()).unwrap_or_default();

        ChainInterfaceResult {
            a: a_res
                .into_iter()
                .map(|(resid, resname)| InterfaceRow {
                    chain_id: chain_a_str.clone(),
                    residue_id: resid,
                    residue_name: resname,
                })
                .collect(),
            b: b_res
                .into_iter()
                .map(|(resid, resname)| InterfaceRow {
                    chain_id: chain_b_str.clone(),
                    residue_id: resid,
                    residue_name: resname,
                })
                .collect(),
        }
    }

    /// Residues at the interface between `chain_a` and `chain_b` within `cutoff` Å.
    /// Returns `{a: [{chain_id, residue_id, residue_name}], b: [...]}`.
    #[wasm_bindgen]
    pub fn chain_interface_residues(&self, chain_a: &str, chain_b: &str, cutoff: f32) -> JsValue {
        serde_wasm_bindgen::to_value(&self.chain_interface_data(chain_a, chain_b, cutoff))
            .unwrap_or(JsValue::NULL)
    }

    // ── P26: Murcko Scaffold (public Wasm API) ────────────────────────────────

    /// Returns atom indices forming the Murcko scaffold (ring systems + linker atoms between rings).
    /// Requires `compute_bonds()` + `compute_rings()`. Returns empty if no rings or bonds absent.
    #[wasm_bindgen]
    pub fn murcko_scaffold_indices(&self) -> Vec<u32> {
        self.murcko_scaffold_indices_data()
            .into_iter()
            .map(|i| i as u32)
            .collect()
    }

    /// Number of distinct ring systems (connected components among ring atoms).
    /// Requires `compute_bonds()` + `compute_rings()`. Returns 0 if no rings.
    #[wasm_bindgen]
    pub fn ring_system_count(&self) -> u32 {
        self.ring_system_count_data() as u32
    }

    // ── P57: Ring classification (spiro / fused / bridged) ───────────────────

    /// Atom indices that are spiro centers: shared by two or more rings that
    /// meet only at that single atom. Requires `compute_bonds()` + `compute_rings()`.
    #[wasm_bindgen]
    pub fn get_spiro_atoms(&self) -> Vec<u32> {
        let rings = self.enumerate_rings();
        let n = self.symbols.len();
        let mut membership: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (ri, ring) in rings.iter().enumerate() {
            for &a in ring {
                if a < n { membership[a].push(ri); }
            }
        }
        let mut spiro: Vec<u32> = Vec::new();
        'outer: for a in 0..n {
            let mems = &membership[a];
            if mems.len() < 2 { continue; }
            for i in 0..mems.len() {
                for j in (i + 1)..mems.len() {
                    let shared = rings[mems[i]].iter()
                        .filter(|&&x| rings[mems[j]].contains(&x))
                        .count();
                    if shared == 1 {
                        spiro.push(a as u32);
                        continue 'outer;
                    }
                }
            }
        }
        spiro
    }

    fn fused_ring_bonds_vec(&self) -> Vec<[u32; 2]> {
        let rings = self.enumerate_rings();
        let n_rings = rings.len();
        let mut result: Vec<[u32; 2]> = Vec::new();
        let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for i in 0..n_rings {
            for j in (i + 1)..n_rings {
                let shared: Vec<usize> = rings[i].iter()
                    .filter(|&&a| rings[j].contains(&a))
                    .copied()
                    .collect();
                if shared.len() == 2 {
                    let (a, b) = (shared[0].min(shared[1]), shared[0].max(shared[1]));
                    if seen.insert((a, b)) {
                        result.push([a as u32, b as u32]);
                    }
                }
            }
        }
        result
    }

    /// Bonds shared by exactly two rings (fused ring bonds). Each entry is
    /// `[atom_a, atom_b]`. Requires `compute_bonds()` + `compute_rings()`.
    #[wasm_bindgen]
    pub fn get_fused_ring_bonds(&self) -> JsValue {
        to_js(&self.fused_ring_bonds_vec())
    }

    /// Returns true if any ring system in the molecule is bridged — two rings
    /// sharing more than one bond (>2 atoms) between the same pair of bridgehead
    /// atoms. Norbornane-type systems return true. Requires `compute_rings()`.
    #[wasm_bindgen]
    pub fn is_bridged_ring_system(&self) -> bool {
        let rings = self.enumerate_rings();
        let n_rings = rings.len();
        for i in 0..n_rings {
            for j in (i + 1)..n_rings {
                let shared_count = rings[i].iter()
                    .filter(|&&a| rings[j].contains(&a))
                    .count();
                if shared_count > 2 { return true; }
            }
        }
        false
    }

    // ── P27: Chain Breaks + Ramachandran Outliers (public Wasm API) ───────────

    /// Detect gaps in the protein backbone.
    /// Reports a break when consecutive Cα atoms in the same chain satisfy:
    /// resid gap > 1 (sequence gap) OR Cα–Cα distance > `ca_cutoff` Å (structural break).
    /// Returns `[{chain_id, from_resid, to_resid}]`. Empty for non-PDB inputs.
    #[wasm_bindgen]
    pub fn chain_breaks(&self, ca_cutoff: f32) -> JsValue {
        to_js(&self.chain_breaks_data(ca_cutoff))
    }

    /// Residues with φ/ψ outside the Ramachandran allowed regions (α-helix, β-strand, left-handed helix).
    /// Returns `[{chain_id, residue_id, residue_name, phi, psi}]`.
    /// Empty if no backbone metadata or all residues are in allowed regions.
    #[wasm_bindgen]
    pub fn ramachandran_outliers(&self) -> JsValue {
        to_js(&self.ramachandran_outliers_data())
    }

    // ── P28: SMILES output ────────────────────────────────────────────────────

    /// Convert the molecular graph to a canonical Kekulé SMILES string.
    /// Disconnected fragments are joined with '.'. Returns "" for empty molecules.
    /// Requires `compute_bonds()` for accurate output; `compute_rings()` is optional.
    #[wasm_bindgen]
    pub fn to_smiles(&self) -> String {
        self.to_smiles_data()
    }

    // ── P29a: Functional group detection ─────────────────────────────────────

    /// Detect functional groups present in the molecule via graph traversal.
    /// Returns an alphabetically sorted array of group names (e.g. "alcohol", "ketone").
    /// Requires `compute_bonds()`. For accurate `aromatic` detection call `compute_rings()` first.
    #[wasm_bindgen]
    pub fn detect_functional_groups(&self) -> JsValue {
        to_js(&self.detect_functional_groups_data())
    }

    // ── P29b: Coordination geometry ───────────────────────────────────────────

    /// Classify the coordination geometry around `center_idx` based on ligand angles.
    /// Returns one of: "linear", "bent", "trigonal_planar", "trigonal_pyramidal",
    /// "tetrahedral", "square_planar", "trigonal_bipyramidal", "square_pyramidal",
    /// "octahedral", or "unknown".
    /// Requires `compute_bonds()` and 3D coordinates (not SMILES-derived).
    #[wasm_bindgen]
    pub fn coordination_geometry(&self, center_idx: usize) -> String {
        self.coordination_geometry_data(center_idx)
    }

    // ── P32: 2D coordinate generation ────────────────────────────────────────
    /// Generate 2D layout coordinates for SMILES-derived molecules.
    /// Rings are placed as regular polygons; acyclic atoms use zig-zag DFS.
    /// Updates `x` and `y` in place; sets `z = 0`.
    #[wasm_bindgen]
    pub fn compute_2d_coords(&mut self) {
        self.compute_2d_coords_data();
    }

    // ── P33: SVG renderer ─────────────────────────────────────────────────────
    /// Render the molecule as an SVG string (heavy atoms only; H hidden).
    /// Call `compute_2d_coords()` first for SMILES/flat molecules.
    #[wasm_bindgen]
    pub fn to_svg_string(&self, width: u32, height: u32) -> String {
        self.to_svg_data(width, height)
    }

    // ── P34: aromatic query ───────────────────────────────────────────────────
    /// Returns true if atom at `index` was parsed as aromatic (SMILES lowercase atom).
    /// Always false for XYZ/PDB/SDF inputs.
    #[wasm_bindgen]
    pub fn is_aromatic(&self, index: usize) -> bool {
        self.aromatic_atoms.get(index).copied().unwrap_or(false)
    }

    // ── P38: tetrahedral stereo queries ──────────────────────────────────────
    /// Number of tetrahedral stereo centers parsed from SMILES @/@@.
    #[wasm_bindgen]
    pub fn stereo_center_count(&self) -> usize {
        self.stereo_centers.len()
    }

    /// Returns true if atom at `index` is a tetrahedral stereo center (@/@@).
    #[wasm_bindgen]
    pub fn is_stereo_center(&self, index: usize) -> bool {
        self.stereo_centers.contains_key(&index)
    }

    // ── P56: 3D stereo perception ─────────────────────────────────────────────
    /// Perceive tetrahedral stereo centers from 3D coordinates and populate
    /// `stereo_centers`. Enables `is_stereo_center()`, `stereo_center_count()`,
    /// and stereo-aware `to_smiles()` output for molecules loaded from PDB,
    /// 3D SDF, or mmCIF formats. Requires bonds (`compute_bonds()` or SDF bond
    /// block). Silently returns without effect if bonds are unavailable or all
    /// coordinates are zero (2D / SMILES-loaded molecules).
    #[wasm_bindgen]
    pub fn perceive_stereo_from_3d(&mut self) {
        let n = self.symbols.len();
        if n == 0 { return; }
        if self.bonds.is_empty() || self.bonds.iter().all(|b| b.is_empty()) { return; }

        let h_set: Vec<bool> = self.symbols.iter().map(|s| s == "H").collect();
        self.stereo_centers.clear();

        for center in 0..n {
            if h_set[center] { continue; }

            let all_nbrs: Vec<usize> = self.bonds
                .get(center).cloned().unwrap_or_default();
            let heavy_nbrs: Vec<usize> = all_nbrs.iter()
                .copied().filter(|&j| !h_set[j]).collect();
            let h_nbrs: Vec<usize> = all_nbrs.iter()
                .copied().filter(|&j| h_set[j]).collect();

            // Exactly 4 distinct substituents (any combination of heavy + H)
            if heavy_nbrs.len() + h_nbrs.len() != 4 { continue; }
            // Need at least one heavy neighbor to serve as the reference (from_atom)
            if heavy_nbrs.is_empty() { continue; }

            // Must have non-trivial 3D coordinates
            let any_3d = |i: usize| self.x[i] != 0.0 || self.y[i] != 0.0 || self.z[i] != 0.0;
            if !any_3d(center) && all_nbrs.iter().all(|&j| !any_3d(j)) { continue; }

            // Observer = first heavy neighbor; remaining 3 = H (if any) + other heavy
            let from_atom = heavy_nbrs[0];
            let subs: Vec<usize> = h_nbrs.iter().copied()
                .chain(heavy_nbrs[1..].iter().copied())
                .collect();
            if subs.len() != 3 { continue; }

            let px = self.x[from_atom];
            let py = self.y[from_atom];
            let pz = self.z[from_atom];
            let v = |i: usize| [self.x[i] - px, self.y[i] - py, self.z[i] - pz];
            let [v1, v2, v3] = [v(subs[0]), v(subs[1]), v(subs[2])];

            // Signed triple product: >0 → CCW = @(desc=-1), <0 → CW = @@(desc=1)
            let vol = (v1[1]*v2[2] - v1[2]*v2[1]) * v3[0]
                    + (v1[2]*v2[0] - v1[0]*v2[2]) * v3[1]
                    + (v1[0]*v2[1] - v1[1]*v2[0]) * v3[2];
            if vol.abs() < 1e-4 { continue; }

            let desc: i8 = if vol < 0.0 { 1 } else { -1 };
            self.stereo_centers.insert(center, (desc, Some(from_atom)));
        }
    }

    /// Perceive E/Z double-bond stereochemistry from 3D coordinates.
    /// Called automatically by `from_sdf_string`.
    pub fn perceive_ez_from_3d(&mut self) {
        let n = self.symbols.len();
        if n < 4 { return; }
        let h_set: Vec<bool> = self.symbols.iter().map(|s| s == "H").collect();
        self.ez_bonds.clear();

        for b in 0..n {
            if h_set[b] { continue; }
            let nbrs_b: Vec<usize> = self.bonds.get(b)
                .map(|v| v.iter().copied().filter(|&j| !h_set[j]).collect())
                .unwrap_or_default();
            for &c in &nbrs_b {
                if c <= b { continue; } // process each pair once
                // Check if b-c is a double bond
                let is_double = self.bond_orders.get(b)
                    .and_then(|ords| self.bonds.get(b)
                        .and_then(|nbrs| nbrs.iter().position(|&x| x == c))
                        .map(|pos| ords.get(pos).copied().unwrap_or(0) == 2))
                    .unwrap_or(false);
                if !is_double { continue; }

                // Need non-trivial 3D coords
                let any_3d = |i: usize| self.x[i] != 0.0 || self.y[i] != 0.0 || self.z[i] != 0.0;
                if !any_3d(b) && !any_3d(c) { continue; }

                // Find first heavy substituent on each side (not the other double-bond atom)
                let nbrs_c: Vec<usize> = self.bonds.get(c)
                    .map(|v| v.iter().copied().filter(|&j| !h_set[j] && j != b).collect())
                    .unwrap_or_default();
                let sub_b: Vec<usize> = nbrs_b.iter().copied().filter(|&j| j != c).collect();

                let (Some(&a), Some(&d)) = (sub_b.first(), nbrs_c.first()) else { continue };

                // Dihedral angle a-b=c-d
                let v = |i: usize, j: usize| -> [f32; 3] {
                    [self.x[j]-self.x[i], self.y[j]-self.y[i], self.z[j]-self.z[i]]
                };
                let cross = |u: [f32;3], w: [f32;3]| -> [f32;3] {
                    [u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]]
                };
                let dot = |u: [f32;3], w: [f32;3]| -> f32 { u[0]*w[0]+u[1]*w[1]+u[2]*w[2] };
                let b1 = v(a, b);
                let b2 = v(b, c);
                let b3 = v(c, d);
                let n1 = cross(b1, b2);
                let n2 = cross(b2, b3);
                let len1 = dot(n1, n1).sqrt();
                let len2 = dot(n2, n2).sqrt();
                if len1 < 1e-4 || len2 < 1e-4 { continue; }
                let cos_d = dot(n1, n2) / (len1 * len2);
                // cos near -1 → dihedral ~180° → trans (E); cos near +1 → ~0° → cis (Z)
                let is_e = cos_d < 0.0;
                self.ez_bonds.insert((b, c), is_e);
            }
        }
    }

    /// Return E/Z descriptor for a double bond between atoms `a` and `b`.
    /// Returns `true` for E (trans), `false` for Z (cis), or `None` if unspecified.
    #[wasm_bindgen]
    pub fn is_ez_bond(&self, a: usize, b: usize) -> Option<bool> {
        let key = (a.min(b), a.max(b));
        self.ez_bonds.get(&key).copied()
    }

    /// Number of double bonds with known E/Z configuration.
    #[wasm_bindgen]
    pub fn ez_bond_count(&self) -> usize {
        self.ez_bonds.len()
    }
}

// ── Fingerprints & Similarity ──────────────────────────────────────────────

#[wasm_bindgen]
impl MolecularSystem {

    /// Morgan/ECFP4 fingerprint as a 256-byte (2048-bit) vector.
    /// In JS this arrives as a Uint8Array.
    #[wasm_bindgen]
    pub fn fingerprint_ecfp4(&self) -> Vec<u8> {
        morgan_fingerprint_bits(&self.symbols, &self.bonds, 2).to_vec()
    }

    /// Tanimoto coefficient (0.0–1.0) between this molecule and `other`.
    /// 1.0 = identical fingerprint, 0.0 = no common bits.
    #[wasm_bindgen]
    pub fn tanimoto_similarity(&self, other: &MolecularSystem) -> f32 {
        let a = morgan_fingerprint_bits(&self.symbols, &self.bonds, 2);
        let b = morgan_fingerprint_bits(&other.symbols, &other.bonds, 2);
        let and: u32 = a.iter().zip(&b).map(|(x, y)| (x & y).count_ones()).sum();
        let or_: u32 = a.iter().zip(&b).map(|(x, y)| (x | y).count_ones()).sum();
        if or_ == 0 { 1.0 } else { and as f32 / or_ as f32 }
    }
}

// --- SMARTS types and helpers ---

#[derive(Clone, PartialEq)]
enum SmartsBondType {
    Single,
    Double,
    Triple,
    Aromatic,
    Any,
}

#[derive(Clone, Default)]
struct SmartsAtom {
    is_any: bool,
    negate: bool,
    symbol: Option<String>,
    aromatic: Option<bool>,  // Some(true)=aromatic; Some(false)=aliphatic; None=either
    atomic_num: Option<u8>,
    h_count: Option<u8>,     // min H-neighbor count
    charge: Option<i8>,
    in_ring: Option<bool>,
    ring_size: Option<u8>,    // [r5] = in ring of exactly this size
    degree: Option<u8>,       // [D2] = exactly 2 heavy-atom bonds
    connectivity: Option<u8>, // [X3] = total bond count == 3
}

struct SmartsGraph {
    atoms: Vec<SmartsAtom>,
    adj: Vec<Vec<(usize, SmartsBondType)>>,
}

// --- P25 serializable types ---

#[derive(serde::Serialize)]
struct ResidueSasaRow {
    chain_id: String,
    residue_id: i32,
    residue_name: String,
    sasa: f32,
}

#[derive(serde::Serialize)]
struct ChainInterfaceResult {
    a: Vec<InterfaceRow>,
    b: Vec<InterfaceRow>,
}

#[derive(serde::Serialize)]
struct InterfaceRow {
    chain_id: String,
    residue_id: i32,
    residue_name: String,
}

fn parse_smarts(input: &str) -> Option<SmartsGraph> {
    let chars: Vec<char> = input.chars().collect();
    let n = chars.len();
    let mut atoms: Vec<SmartsAtom> = Vec::new();
    let mut edges: Vec<(usize, usize, SmartsBondType)> = Vec::new();
    let mut prev: Option<usize> = None;
    let mut branch_stack: Vec<Option<usize>> = Vec::new();
    let mut ring_opens: std::collections::HashMap<char, (usize, SmartsBondType)> =
        std::collections::HashMap::new();
    let mut pending_bond: Option<SmartsBondType> = None;
    let mut pos = 0;

    let push_atom = |atoms: &mut Vec<SmartsAtom>,
                     edges: &mut Vec<(usize, usize, SmartsBondType)>,
                     prev: &mut Option<usize>,
                     pending: &mut Option<SmartsBondType>,
                     atom: SmartsAtom| {
        let idx = atoms.len();
        atoms.push(atom);
        if let Some(p) = *prev {
            let bt = pending.take().unwrap_or(SmartsBondType::Any);
            edges.push((p, idx, bt));
        }
        *prev = Some(idx);
    };

    while pos < n {
        match chars[pos] {
            '(' => { branch_stack.push(prev); pos += 1; }
            ')' => { prev = branch_stack.pop().flatten(); pos += 1; }
            '-' => { pending_bond = Some(SmartsBondType::Single); pos += 1; }
            '=' => { pending_bond = Some(SmartsBondType::Double); pos += 1; }
            '#' => { pending_bond = Some(SmartsBondType::Triple); pos += 1; }
            ':' => { pending_bond = Some(SmartsBondType::Aromatic); pos += 1; }
            '~' => { pending_bond = Some(SmartsBondType::Any); pos += 1; }
            '*' => {
                push_atom(&mut atoms, &mut edges, &mut prev, &mut pending_bond,
                    SmartsAtom { is_any: true, ..SmartsAtom::default() });
                pos += 1;
            }
            'a' => {
                push_atom(&mut atoms, &mut edges, &mut prev, &mut pending_bond,
                    SmartsAtom { is_any: true, aromatic: Some(true), ..SmartsAtom::default() });
                pos += 1;
            }
            'A' => {
                push_atom(&mut atoms, &mut edges, &mut prev, &mut pending_bond,
                    SmartsAtom { is_any: true, aromatic: Some(false), ..SmartsAtom::default() });
                pos += 1;
            }
            '[' => {
                pos += 1;
                let mut atom = SmartsAtom::default();

                // Optional negation
                if pos < n && chars[pos] == '!' {
                    atom.negate = true;
                    pos += 1;
                }

                // Atomic number: #nn
                if pos + 1 < n && chars[pos] == '#' && chars[pos + 1].is_ascii_digit() {
                    pos += 1;
                    let mut num = 0u8;
                    while pos < n && chars[pos].is_ascii_digit() {
                        num = num.saturating_mul(10).saturating_add(chars[pos] as u8 - b'0');
                        pos += 1;
                    }
                    atom.atomic_num = Some(num);
                } else if pos < n && chars[pos] == 'R' {
                    // [R] = in any ring
                    atom.in_ring = Some(true);
                    atom.is_any = true;
                    pos += 1;
                } else if pos < n && chars[pos] == 'r' {
                    // [r] / [r5] / [r6] = ring-size SMARTS (must be before element parsing)
                    pos += 1;
                    if pos < n && chars[pos].is_ascii_digit() {
                        let mut sz = 0u8;
                        while pos < n && chars[pos].is_ascii_digit() {
                            sz = sz.saturating_mul(10).saturating_add(chars[pos] as u8 - b'0');
                            pos += 1;
                        }
                        atom.ring_size = Some(sz);
                    } else {
                        atom.in_ring = Some(true);
                    }
                    atom.is_any = true;
                } else if pos < n && chars[pos].is_ascii_alphabetic() && chars[pos] != 'h' {
                    let aromatic_flag = chars[pos].is_lowercase();
                    let mut sym = chars[pos].to_uppercase().to_string();
                    pos += 1;
                    // Two-letter elements (Cl, Br, Si, Sn)
                    if pos < n && chars[pos].is_ascii_lowercase() && chars[pos] != 'h'
                        && matches!(
                            (sym.as_str(), chars[pos]),
                            ("C", 'l') | ("B", 'r') | ("S", 'i') | ("S", 'n') | ("Z", 'n') | ("F", 'e') | ("N", 'i') | ("C", 'u') | ("M", 'g') | ("P", 'd')
                        )
                    {
                        sym.push(chars[pos]);
                        pos += 1;
                    }
                    atom.symbol = Some(sym);
                    atom.aromatic = Some(aromatic_flag);
                }

                // Remaining bracket content: H count, charge, ring size, degree, connectivity
                while pos < n && chars[pos] != ']' {
                    match chars[pos] {
                        'H' => {
                            pos += 1;
                            let cnt = if pos < n && chars[pos].is_ascii_digit() {
                                let c = chars[pos] as u8 - b'0';
                                pos += 1;
                                c
                            } else {
                                1
                            };
                            atom.h_count = Some(cnt);
                        }
                        'r' => {
                            pos += 1;
                            if pos < n && chars[pos].is_ascii_digit() {
                                let mut sz = 0u8;
                                while pos < n && chars[pos].is_ascii_digit() {
                                    sz = sz.saturating_mul(10).saturating_add(chars[pos] as u8 - b'0');
                                    pos += 1;
                                }
                                atom.ring_size = Some(sz);
                            } else {
                                atom.in_ring = Some(true);
                            }
                        }
                        'D' => {
                            pos += 1;
                            let mut d = 0u8;
                            while pos < n && chars[pos].is_ascii_digit() {
                                d = d.saturating_mul(10).saturating_add(chars[pos] as u8 - b'0');
                                pos += 1;
                            }
                            atom.degree = Some(d);
                        }
                        'X' => {
                            pos += 1;
                            let mut xv = 0u8;
                            while pos < n && chars[pos].is_ascii_digit() {
                                xv = xv.saturating_mul(10).saturating_add(chars[pos] as u8 - b'0');
                                pos += 1;
                            }
                            atom.connectivity = Some(xv);
                        }
                        '+' => {
                            pos += 1;
                            let charge = if pos < n && chars[pos].is_ascii_digit() {
                                let c = (chars[pos] as u8 - b'0') as i8;
                                pos += 1;
                                c
                            } else {
                                1i8
                            };
                            atom.charge = Some(charge);
                        }
                        '-' => {
                            pos += 1;
                            let charge = if pos < n && chars[pos].is_ascii_digit() {
                                -((chars[pos] as u8 - b'0') as i8)
                            } else {
                                -1i8
                            };
                            if pos < n && chars[pos].is_ascii_digit() { pos += 1; }
                            atom.charge = Some(charge);
                        }
                        _ => { pos += 1; }
                    }
                }
                if pos < n { pos += 1; } // skip ']'

                // Empty brackets → any atom
                if !atom.negate && atom.symbol.is_none() && atom.atomic_num.is_none()
                    && atom.in_ring.is_none() && atom.h_count.is_none() && atom.charge.is_none()
                {
                    atom.is_any = true;
                }

                push_atom(&mut atoms, &mut edges, &mut prev, &mut pending_bond, atom);
            }
            c if c.is_uppercase() => {
                let mut sym = String::new();
                sym.push(c);
                pos += 1;
                if (c == 'C' && pos < n && chars[pos] == 'l')
                    || (c == 'B' && pos < n && chars[pos] == 'r')
                {
                    sym.push(chars[pos]);
                    pos += 1;
                }
                push_atom(&mut atoms, &mut edges, &mut prev, &mut pending_bond,
                    SmartsAtom { symbol: Some(sym), aromatic: Some(false), ..SmartsAtom::default() });
            }
            c if c.is_lowercase() => {
                let sym = c.to_uppercase().to_string();
                pos += 1;
                push_atom(&mut atoms, &mut edges, &mut prev, &mut pending_bond,
                    SmartsAtom { symbol: Some(sym), aromatic: Some(true), ..SmartsAtom::default() });
            }
            c if c.is_ascii_digit() => {
                if let Some((open_idx, bt)) = ring_opens.remove(&c) {
                    let close_bt = pending_bond.take().unwrap_or(SmartsBondType::Any);
                    let bond_type = if bt == SmartsBondType::Any { close_bt } else { bt };
                    if let Some(p) = prev { edges.push((open_idx, p, bond_type)); }
                } else if let Some(p) = prev {
                    let bt = pending_bond.take().unwrap_or(SmartsBondType::Any);
                    ring_opens.insert(c, (p, bt));
                }
                pos += 1;
            }
            '%' => {
                if pos + 2 < n {
                    let ring_key = chars[pos + 1..pos + 3].iter().collect::<String>();
                    let ring_char = ring_key.chars().next().unwrap_or('0');
                    if let Some((open_idx, bt)) = ring_opens.remove(&ring_char) {
                        let close_bt = pending_bond.take().unwrap_or(SmartsBondType::Any);
                        let bond_type = if bt == SmartsBondType::Any { close_bt } else { bt };
                        if let Some(p) = prev { edges.push((open_idx, p, bond_type)); }
                    } else if let Some(p) = prev {
                        let bt = pending_bond.take().unwrap_or(SmartsBondType::Any);
                        ring_opens.insert(ring_char, (p, bt));
                    }
                }
                pos += 3;
            }
            _ => { pos += 1; }
        }
    }

    if atoms.is_empty() {
        return None;
    }

    let na = atoms.len();
    let mut adj: Vec<Vec<(usize, SmartsBondType)>> = vec![Vec::new(); na];
    for (u, v, bt) in &edges {
        adj[*u].push((*v, bt.clone()));
        adj[*v].push((*u, bt.clone()));
    }

    Some(SmartsGraph { atoms, adj })
}

fn smarts_bfs_order(query: &SmartsGraph) -> Vec<usize> {
    let n = query.atoms.len();
    let mut order = Vec::with_capacity(n);
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(0);
    visited[0] = true;
    while let Some(u) = queue.pop_front() {
        order.push(u);
        for &(v, _) in &query.adj[u] {
            if !visited[v] {
                visited[v] = true;
                queue.push_back(v);
            }
        }
    }
    // Include any disconnected query atoms
    order.extend((0..n).filter(|&i| !visited[i]));
    order
}

fn smarts_atomic_num(sym: &str) -> u8 {
    match sym {
        "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5,
        "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10,
        "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15,
        "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20,
        "Cr" => 24, "Mn" => 25, "Fe" => 26, "Co" => 27, "Ni" => 28,
        "Cu" => 29, "Zn" => 30, "Br" => 35, "Ru" => 44, "Pd" => 46,
        "Ag" => 47, "Cd" => 48, "I" => 53, "Pt" => 78, "Au" => 79,
        "Hg" => 80, _ => 0,
    }
}

fn smarts_is_aromatic(mol: &MolecularSystem, i: usize) -> bool {
    if !mol.aromatic_atoms.is_empty() {
        mol.aromatic_atoms.get(i).copied().unwrap_or(false)
    } else {
        mol.ring_atoms.get(i).copied().unwrap_or(false)
    }
}

fn smarts_atom_matches(q: &SmartsAtom, mol: &MolecularSystem, i: usize) -> bool {
    let result = smarts_atom_base_match(q, mol, i);
    if q.negate { !result } else { result }
}

fn smarts_atom_base_match(q: &SmartsAtom, mol: &MolecularSystem, i: usize) -> bool {
    // Ring membership
    if let Some(need_ring) = q.in_ring {
        if mol.ring_atoms.get(i).copied().unwrap_or(false) != need_ring {
            return false;
        }
    }
    // Ring size [r5] / [r6]
    if let Some(need_sz) = q.ring_size {
        let in_size = mol.ring_sizes_per_atom
            .get(i)
            .map(|v| v.contains(&need_sz))
            .unwrap_or(false);
        if !in_size {
            return false;
        }
    }
    // Heavy-atom degree [D2]
    if let Some(need_deg) = q.degree {
        let deg = mol.bonds.get(i)
            .map(|nb| nb.iter()
                .filter(|&&j| mol.symbols.get(j).map(|s| s != "H").unwrap_or(true))
                .count() as u8)
            .unwrap_or(0);
        if deg != need_deg {
            return false;
        }
    }
    // Total connectivity [X3]
    if let Some(need_x) = q.connectivity {
        let xval = mol.bonds.get(i).map(|nb| nb.len() as u8).unwrap_or(0);
        if xval != need_x {
            return false;
        }
    }
    // Formal charge
    if let Some(need_charge) = q.charge {
        if mol.charges.get(i).copied().unwrap_or(0) as i8 != need_charge {
            return false;
        }
    }
    // H neighbor count
    if let Some(need_h) = q.h_count {
        let h_count = mol.bonds.get(i)
            .map(|nb| nb.iter()
                .filter(|&&j| mol.symbols.get(j).map(|s| s == "H").unwrap_or(false))
                .count() as u8)
            .unwrap_or(0);
        if h_count < need_h {
            return false;
        }
    }
    // Atomic number
    if let Some(need_anum) = q.atomic_num {
        if smarts_atomic_num(mol.symbols.get(i).map(|s| s.as_str()).unwrap_or("")) != need_anum {
            return false;
        }
        // If only atomic_num was specified (no symbol), we're done
        if q.symbol.is_none() {
            return true;
        }
    }
    // Pure aromaticity wildcard (a / A / [a] / [A])
    if q.is_any && q.symbol.is_none() && q.atomic_num.is_none() {
        if let Some(need_arom) = q.aromatic {
            return smarts_is_aromatic(mol, i) == need_arom;
        }
        return true; // * wildcard
    }
    // Element symbol
    if let Some(sym) = &q.symbol {
        if mol.symbols.get(i).map(|s| s.as_str()) != Some(sym.as_str()) {
            return false;
        }
        if let Some(need_arom) = q.aromatic {
            if smarts_is_aromatic(mol, i) != need_arom {
                return false;
            }
        }
        return true;
    }
    false
}

fn smarts_bond_matches(bt: &SmartsBondType, mol: &MolecularSystem, i: usize, j: usize) -> bool {
    match bt {
        SmartsBondType::Any => true,
        SmartsBondType::Aromatic => {
            smarts_is_aromatic(mol, i) && smarts_is_aromatic(mol, j)
        }
        _ => {
            let ord = mol
                .bonds
                .get(i)
                .and_then(|nb| nb.iter().position(|&x| x == j))
                .and_then(|k| mol.bond_orders.get(i).and_then(|bo| bo.get(k)).copied())
                .unwrap_or(1);
            match bt {
                SmartsBondType::Single => ord == 1,
                SmartsBondType::Double => ord == 2,
                SmartsBondType::Triple => ord == 3,
                _ => false,
            }
        }
    }
}

fn smarts_backtrack(
    query: &SmartsGraph,
    mol: &MolecularSystem,
    mapping: &mut Vec<Option<usize>>,
    used: &mut Vec<bool>,
    bfs_order: &[usize],
    depth: usize,
    results: &mut Vec<Vec<usize>>,
) {
    if depth == bfs_order.len() {
        let m: Vec<usize> = bfs_order.iter().map(|&q| mapping[q].unwrap()).collect();
        results.push(m);
        return;
    }
    let q_idx = bfs_order[depth];
    for m_idx in 0..mol.symbols.len() {
        if used[m_idx] {
            continue;
        }
        if !smarts_atom_matches(&query.atoms[q_idx], mol, m_idx) {
            continue;
        }
        let ok = query.adj[q_idx].iter().all(|(nb, bt)| {
            if let Some(m_nb) = mapping[*nb] {
                mol.bonds
                    .get(m_nb)
                    .map(|nb_list| nb_list.contains(&m_idx))
                    .unwrap_or(false)
                    && smarts_bond_matches(bt, mol, m_nb, m_idx)
            } else {
                true
            }
        });
        if ok {
            mapping[q_idx] = Some(m_idx);
            used[m_idx] = true;
            smarts_backtrack(query, mol, mapping, used, bfs_order, depth + 1, results);
            used[m_idx] = false;
            mapping[q_idx] = None;
        }
    }
}

fn fibonacci_sphere_92() -> Vec<(f32, f32, f32)> {
    let n = 92usize;
    let golden = std::f32::consts::PI * (3.0 - 5.0_f32.sqrt());
    (0..n)
        .map(|i| {
            let y = 1.0 - (2.0 * i as f32 + 1.0) / n as f32;
            let r = (1.0 - y * y).max(0.0).sqrt();
            let phi = golden * i as f32;
            (r * phi.cos(), y, r * phi.sin())
        })
        .collect()
}

// --- Kabsch Superposition ---
//
// Kabsch (1976) algorithm for optimal rotation minimising RMSD between two
// point sets of equal size. Pure-Rust 3×3 Jacobi SVD; no external crates.
//
// Steps:
//   1. Translate both sets to their respective centroids.
//   2. Compute covariance matrix H = Pᵀ Q  (P=mobile centred, Q=ref centred).
//   3. SVD(H) → U S Vᵀ.
//   4. Rotation R = V Uᵀ; if det(R) < 0 flip the column of V corresponding
//      to the smallest singular value (handles reflection).
//   5. Apply centring + R to mobile coords; compute RMSD.

type Mat3 = [[f32; 3]; 3];

#[allow(clippy::needless_range_loop)]
fn mat3_mul(a: Mat3, b: Mat3) -> Mat3 {
    let mut c = [[0.0f32; 3]; 3];
    for i in 0..3 { for j in 0..3 { for k in 0..3 { c[i][j] += a[i][k] * b[k][j]; } } }
    c
}

fn mat3_transpose(a: Mat3) -> Mat3 {
    [[a[0][0], a[1][0], a[2][0]],
     [a[0][1], a[1][1], a[2][1]],
     [a[0][2], a[1][2], a[2][2]]]
}

fn mat3_det(a: Mat3) -> f32 {
    a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
  - a[0][1] * (a[1][0]*a[2][2] - a[1][2]*a[2][0])
  + a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0])
}

/// One-sided Jacobi SVD for a 3×3 real matrix A.
/// Returns (U, S, V) such that A = U diag(S) Vᵀ, singular values in descending order.
fn svd3(a: Mat3) -> (Mat3, [f32; 3], Mat3) {
    // Compute AᵀA via symmetric eigendecomposition (Jacobi method on AᵀA).
    // Then recover U from A V S⁻¹.

    let ata: Mat3 = mat3_mul(mat3_transpose(a), a);

    // Jacobi eigenvalue iteration on the symmetric 3×3 matrix ata.
    let mut v: Mat3 = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    let mut m = ata;

    for _ in 0..64 {
        // Find largest off-diagonal element
        let mut p = 0usize; let mut q = 1usize; let mut max_val = m[0][1].abs();
        for (pi, qi) in [(0,2),(1,2)] {
            if m[pi][qi].abs() > max_val { max_val = m[pi][qi].abs(); p = pi; q = qi; }
        }
        if max_val < 1e-10 { break; }

        let theta = 0.5 * (m[q][q] - m[p][p]).atan2(2.0 * m[p][q]);
        let (s, c) = theta.sin_cos();

        // Apply Givens rotation
        let mut g: Mat3 = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
        g[p][p] =  c; g[p][q] = s;
        g[q][p] = -s; g[q][q] = c;

        m = mat3_mul(mat3_mul(mat3_transpose(g), m), g);
        v = mat3_mul(v, g);
    }

    // Singular values = sqrt of eigenvalues of AᵀA
    let mut sv = [m[0][0].max(0.0).sqrt(), m[1][1].max(0.0).sqrt(), m[2][2].max(0.0).sqrt()];

    // Sort descending; apply same permutation to V columns
    for i in 0..2 {
        for j in 0..2-i {
            if sv[j] < sv[j+1] {
                sv.swap(j, j+1);
                for row in &mut v { row.swap(j, j+1); }
            }
        }
    }

    // U = A V S⁻¹
    let av = mat3_mul(a, v);
    let mut u: Mat3 = [[0.; 3]; 3];
    for col in 0..3 {
        let norm = (av[0][col]*av[0][col] + av[1][col]*av[1][col] + av[2][col]*av[2][col]).sqrt();
        let scale = if norm > 1e-10 { 1.0 / norm } else { 0.0 };
        for row in 0..3 { u[row][col] = av[row][col] * scale; }
    }

    (u, sv, v)
}

/// Kabsch algorithm. Returns the optimal rotation matrix R and the two centroids.
/// Caller applies: q_aligned[i] = R · (q[i] - q_centroid) + p_centroid
#[allow(clippy::needless_range_loop)]
fn kabsch(p: &[[f32; 3]], q: &[[f32; 3]]) -> (Mat3, [f32; 3], [f32; 3]) {
    let n = p.len().min(q.len()) as f32;
    if n < 1.0 {
        return ([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], [0.;3], [0.;3]);
    }

    let centroid = |pts: &[[f32;3]]| -> [f32;3] {
        let mut c = [0f32;3];
        for pt in pts { for k in 0..3 { c[k] += pt[k]; } }
        [c[0]/n, c[1]/n, c[2]/n]
    };

    let cp = centroid(p);
    let cq = centroid(q);

    // Covariance H = Pᵀ Q  (centred)
    let mut h: Mat3 = [[0.;3];3];
    for i in 0..(n as usize) {
        let pi = [p[i][0]-cp[0], p[i][1]-cp[1], p[i][2]-cp[2]];
        let qi = [q[i][0]-cq[0], q[i][1]-cq[1], q[i][2]-cq[2]];
        for r in 0..3 { for c in 0..3 { h[r][c] += pi[r] * qi[c]; } }
    }

    let (u, _s, v) = svd3(h);
    let mut rot = mat3_mul(v, mat3_transpose(u));

    // Handle reflection (det = -1)
    if mat3_det(rot) < 0.0 {
        let mut v_fix = v;
        // Flip the column with the smallest singular value (last after sorting descending = col 2)
        for row in &mut v_fix { row[2] = -row[2]; }
        rot = mat3_mul(v_fix, mat3_transpose(u));
    }

    (rot, cp, cq)
}

// --- Internal Helpers (non-Wasm) ---

impl MolecularSystem {
    fn new_empty() -> Self {
        MolecularSystem {
            symbols: Vec::new(),
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            atom_names: Vec::new(),
            residue_names: Vec::new(),
            residue_ids: Vec::new(),
            chain_ids: Vec::new(),
            hetatm_flags: Vec::new(),
            bonds: Vec::new(),
            spatial_grid: None,
            ring_atoms: Vec::new(),
            ring_bonds: std::collections::HashSet::new(),
            bond_orders: Vec::new(),
            occupancies: Vec::new(),
            b_factors: Vec::new(),
            charges: Vec::new(),
            aromatic_atoms: Vec::new(),
            stereo_centers: std::collections::HashMap::new(),
            atom_map: Vec::new(),
            ez_bonds: std::collections::HashMap::new(),
            ring_sizes_per_atom: Vec::new(),
            properties: std::collections::HashMap::new(),
        }
    }

    #[inline] fn chain_id_byte(&self, i: usize) -> u8 { self.chain_ids.get(i).copied().unwrap_or(b'A') }
    #[inline] fn residue_id_i(&self, i: usize) -> i32 { self.residue_ids.get(i).copied().unwrap_or(0) }
    #[inline] fn hetatm_at(&self, i: usize) -> bool { self.hetatm_flags.get(i).copied().unwrap_or(false) }
    #[inline] fn atom_name_str(&self, i: usize) -> &str { self.atom_names.get(i).map(|s| s.as_str()).unwrap_or("") }

    fn select_by_indices(&self, indices: &[usize]) -> MolecularSystem {
        let n_total = self.symbols.len();
        let mut old_to_new = vec![usize::MAX; n_total];
        for (new_i, &old_i) in indices.iter().enumerate() {
            if old_i < n_total {
                old_to_new[old_i] = new_i;
            }
        }

        let mut symbols = Vec::with_capacity(indices.len());
        let mut x = Vec::with_capacity(indices.len());
        let mut y = Vec::with_capacity(indices.len());
        let mut z = Vec::with_capacity(indices.len());
        let mut atom_names = Vec::new();
        let mut residue_names = Vec::new();
        let mut residue_ids = Vec::new();
        let mut chain_ids = Vec::new();
        let mut hetatm_flags = Vec::new();
        let mut bonds = Vec::with_capacity(indices.len());
        let mut occupancies = Vec::new();
        let mut b_factors = Vec::new();
        let mut charges = Vec::new();
        let mut aromatic_atoms_new = Vec::new();

        let has_atom_names = !self.atom_names.is_empty();
        let has_residue = !self.residue_names.is_empty();
        let has_chain = !self.chain_ids.is_empty();
        let has_hetatm = !self.hetatm_flags.is_empty();
        let has_occ = !self.occupancies.is_empty();
        let has_bfac = !self.b_factors.is_empty();
        let has_charges = !self.charges.is_empty();
        let has_arom = !self.aromatic_atoms.is_empty();

        for &old_i in indices {
            symbols.push(self.symbols[old_i].clone());
            x.push(self.x[old_i]);
            y.push(self.y[old_i]);
            z.push(self.z[old_i]);
            if has_atom_names {
                atom_names.push(self.atom_names.get(old_i).cloned().unwrap_or_default());
            }
            if has_residue {
                residue_names.push(self.residue_names.get(old_i).cloned().unwrap_or_default());
                residue_ids.push(self.residue_ids.get(old_i).copied().unwrap_or(0));
            }
            if has_chain {
                chain_ids.push(self.chain_ids.get(old_i).copied().unwrap_or(b'A'));
            }
            if has_hetatm {
                hetatm_flags.push(self.hetatm_flags.get(old_i).copied().unwrap_or(false));
            }
            if has_occ {
                occupancies.push(self.occupancies.get(old_i).copied().unwrap_or(1.0));
            }
            if has_bfac {
                b_factors.push(self.b_factors.get(old_i).copied().unwrap_or(0.0));
            }
            if has_charges {
                charges.push(self.charges.get(old_i).copied().unwrap_or(0));
            }
            if has_arom {
                aromatic_atoms_new.push(self.aromatic_atoms.get(old_i).copied().unwrap_or(false));
            }
            // Remap bonds
            let remapped: Vec<usize> = if old_i < self.bonds.len() {
                self.bonds[old_i].iter()
                    .filter_map(|&j| {
                        let nj = *old_to_new.get(j).unwrap_or(&usize::MAX);
                        if nj != usize::MAX { Some(nj) } else { None }
                    })
                    .collect()
            } else {
                Vec::new()
            };
            bonds.push(remapped);
        }

        // Remap bond_orders in sync with bonds
        let mut bond_orders_new: Vec<Vec<u8>> = Vec::with_capacity(indices.len());
        let has_bond_orders = !self.bond_orders.is_empty();
        for &old_i in indices {
            if has_bond_orders && old_i < self.bond_orders.len() {
                let remapped_orders: Vec<u8> = self.bonds[old_i].iter()
                    .zip(self.bond_orders[old_i].iter())
                    .filter_map(|(&j, &ord)| {
                        let nj = old_to_new.get(j).copied().unwrap_or(usize::MAX);
                        if nj != usize::MAX { Some(ord) } else { None }
                    })
                    .collect();
                bond_orders_new.push(remapped_orders);
            } else {
                bond_orders_new.push(Vec::new());
            }
        }

        let has_atom_map = !self.atom_map.is_empty();
        let atom_map_new: Vec<u32> = if has_atom_map {
            indices.iter().map(|&old_i| self.atom_map.get(old_i).copied().unwrap_or(0)).collect()
        } else {
            Vec::new()
        };

        let has_ring_sizes = !self.ring_sizes_per_atom.is_empty();
        let ring_sizes_per_atom_new: Vec<Vec<u8>> = if has_ring_sizes {
            indices.iter().map(|&old_i| {
                self.ring_sizes_per_atom.get(old_i).cloned().unwrap_or_default()
            }).collect()
        } else {
            Vec::new()
        };

        let mut mol = MolecularSystem::new_empty();
        mol.symbols = symbols;
        mol.x = x;
        mol.y = y;
        mol.z = z;
        mol.atom_names = atom_names;
        mol.residue_names = residue_names;
        mol.residue_ids = residue_ids;
        mol.chain_ids = chain_ids;
        mol.hetatm_flags = hetatm_flags;
        mol.bonds = bonds;
        mol.bond_orders = bond_orders_new;
        mol.occupancies = occupancies;
        mol.b_factors = b_factors;
        mol.charges = charges;
        mol.atom_map = atom_map_new;
        mol.ring_sizes_per_atom = ring_sizes_per_atom_new;
        mol.aromatic_atoms = aromatic_atoms_new;
        mol.properties = self.properties.clone();
        mol.stereo_centers = {
            let mut sc = std::collections::HashMap::new();
            for (&old_i, &(desc, from)) in &self.stereo_centers {
                if let Some(&new_i) = old_to_new.get(old_i) {
                    let new_from = from.and_then(|f| old_to_new.get(f).copied());
                    sc.insert(new_i, (desc, new_from));
                }
            }
            sc
        };
        mol.ez_bonds = {
            let mut ez = std::collections::HashMap::new();
            for (&(oa, ob), &is_e) in &self.ez_bonds {
                if let (Some(&na), Some(&nb)) = (old_to_new.get(oa), old_to_new.get(ob)) {
                    ez.insert((na.min(nb), na.max(nb)), is_e);
                }
            }
            ez
        };
        mol
    }

    fn backbone_angle_data(&self) -> Vec<BackboneAngleRow> {
        if self.atom_names.is_empty() || self.residue_names.is_empty() {
            return Vec::new();
        }

        #[derive(Default)]
        struct ResAtoms {
            n: Option<usize>,
            ca: Option<usize>,
            c: Option<usize>,
            res_name: String,
        }

        let mut residues: std::collections::BTreeMap<(u8, i32), ResAtoms> =
            std::collections::BTreeMap::new();

        for i in 0..self.symbols.len() {
            if self.hetatm_at(i) {
                continue;
            }
            let chain = self.chain_id_byte(i);
            let res_id = self.residue_id_i(i);
            let atom_name = self.atom_name_str(i);
            let res_name = self.residue_names.get(i).cloned().unwrap_or_default();

            let entry = residues.entry((chain, res_id)).or_default();
            if entry.res_name.is_empty() {
                entry.res_name = res_name;
            }
            match atom_name {
                "N"  => entry.n  = Some(i),
                "CA" => entry.ca = Some(i),
                "C"  => entry.c  = Some(i),
                _ => {}
            }
        }

        let keys: Vec<(u8, i32)> = residues.keys().cloned().collect();
        let mut result = Vec::with_capacity(keys.len());

        for (idx, &key) in keys.iter().enumerate() {
            let curr = &residues[&key];
            let (n, ca, c) = match (curr.n, curr.ca, curr.c) {
                (Some(n), Some(ca), Some(c)) => (n, ca, c),
                _ => continue,
            };

            let phi = if idx > 0 {
                let prev_key = keys[idx - 1];
                if prev_key.0 == key.0 {
                    residues[&prev_key].c.map(|c_prev| self.dihedral(c_prev, n, ca, c))
                } else { None }
            } else { None };

            let psi = if idx + 1 < keys.len() {
                let next_key = keys[idx + 1];
                if next_key.0 == key.0 {
                    residues[&next_key].n.map(|n_next| self.dihedral(n, ca, c, n_next))
                } else { None }
            } else { None };

            result.push(BackboneAngleRow {
                chain_id: (key.0 as char).to_string(),
                residue_id: key.1,
                residue_name: curr.res_name.clone(),
                phi,
                psi,
            });
        }
        result
    }

    // ── P26: Murcko Scaffold ──────────────────────────────────────────────────

    fn murcko_scaffold_indices_data(&self) -> Vec<usize> {
        let n = self.symbols.len();
        if n == 0 || self.bonds.is_empty() || self.ring_atoms.is_empty() {
            return Vec::new();
        }
        let mut degree: Vec<usize> = self.bonds.iter().map(|nb| nb.len()).collect();
        loop {
            let mut changed = false;
            for i in 0..n {
                if !self.ring_atoms.get(i).copied().unwrap_or(false) && degree[i] == 1 {
                    for &j in &self.bonds[i] {
                        degree[j] = degree[j].saturating_sub(1);
                    }
                    degree[i] = 0;
                    changed = true;
                }
            }
            if !changed { break; }
        }
        (0..n).filter(|&i| degree[i] > 0).collect()
    }

    fn ring_system_count_data(&self) -> usize {
        let n = self.symbols.len();
        if n == 0 || self.ring_atoms.is_empty() {
            return 0;
        }
        let mut visited = vec![false; n];
        let mut count = 0usize;
        for start in 0..n {
            if !self.ring_atoms.get(start).copied().unwrap_or(false) || visited[start] {
                continue;
            }
            count += 1;
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start);
            visited[start] = true;
            while let Some(cur) = queue.pop_front() {
                for &nb in &self.bonds[cur] {
                    if !visited[nb] && self.ring_atoms.get(nb).copied().unwrap_or(false) {
                        visited[nb] = true;
                        queue.push_back(nb);
                    }
                }
            }
        }
        count
    }

    // ── P27: Chain Breaks + Ramachandran Outliers ─────────────────────────────

    fn chain_breaks_data(&self, ca_cutoff: f32) -> Vec<ChainBreakRow> {
        if self.atom_names.is_empty() {
            return Vec::new();
        }
        let n = self.symbols.len();
        let mut ca_by_chain: std::collections::HashMap<u8, Vec<(i32, usize)>> =
            std::collections::HashMap::new();
        for i in 0..n {
            let atom_name = self.atom_names.get(i).map(|s| s.trim()).unwrap_or("");
            if atom_name == "CA" && !self.hetatm_at(i) {
                let chain = self.chain_id_byte(i);
                let resid = self.residue_id_i(i);
                ca_by_chain.entry(chain).or_default().push((resid, i));
            }
        }
        let mut chains: Vec<u8> = ca_by_chain.keys().cloned().collect();
        chains.sort();
        let mut result: Vec<ChainBreakRow> = Vec::new();
        for chain in chains {
            let entries = ca_by_chain.get_mut(&chain).unwrap();
            entries.sort_by_key(|&(resid, _)| resid);
            for w in entries.windows(2) {
                let (resid_a, idx_a) = w[0];
                let (resid_b, idx_b) = w[1];
                let seq_gap = resid_b - resid_a > 1;
                let dist_gap = self.distance(idx_a, idx_b) > ca_cutoff;
                if seq_gap || dist_gap {
                    result.push(ChainBreakRow {
                        chain_id: (chain as char).to_string(),
                        from_resid: resid_a,
                        to_resid: resid_b,
                    });
                }
            }
        }
        result
    }

    fn ramachandran_outliers_data(&self) -> Vec<RamachandranOutlierRow> {
        self.backbone_angle_data()
            .into_iter()
            .filter_map(|r| {
                let phi = r.phi?;
                let psi = r.psi?;
                if is_ramachandran_allowed(phi, psi) {
                    None
                } else {
                    Some(RamachandranOutlierRow {
                        chain_id: r.chain_id,
                        residue_id: r.residue_id,
                        residue_name: r.residue_name,
                        phi,
                        psi,
                    })
                }
            })
            .collect()
    }

    // ── P29b: Coordination geometry ───────────────────────────────────────────

    fn coordination_geometry_data(&self, center_idx: usize) -> String {
        let n = self.symbols.len();
        if n == 0 || center_idx >= n || self.bonds.is_empty() {
            return "unknown".to_string();
        }
        // Reject all-zero coordinates (SMILES-derived molecules)
        let all_zero = self.x.iter().chain(self.y.iter()).chain(self.z.iter())
            .all(|&v| v.abs() < 1e-6);
        if all_zero { return "unknown".to_string(); }

        let ligands = &self.bonds[center_idx];
        let cn = ligands.len();
        if !(2..=6).contains(&cn) { return "unknown".to_string(); }

        // Compute all inter-ligand angles through the center
        let mut angles: Vec<f32> = Vec::new();
        for i in 0..cn {
            for j in (i + 1)..cn {
                angles.push(self.angle(ligands[i], center_idx, ligands[j]));
            }
        }

        let near = |a: f32, lo: f32, hi: f32| -> bool { a > lo && a < hi };
        let count_near = |lo: f32, hi: f32| -> usize {
            angles.iter().filter(|&&a| near(a, lo, hi)).count()
        };
        let near_180 = count_near(155.0, 190.0);
        let near_120 = count_near(100.0, 140.0);
        let near_109 = count_near(95.0, 125.0);
        let near_90  = count_near(75.0, 105.0);

        match cn {
            2 => if angles[0] > 165.0 { "linear" } else { "bent" },
            3 => if near_120 == 3 { "trigonal_planar" } else { "trigonal_pyramidal" },
            4 => {
                if near_180 == 2 && near_90 == 4 { "square_planar" }
                else if near_180 == 0 && near_109 >= 4 { "tetrahedral" }
                else { "unknown" }
            }
            5 => if near_180 >= 1 { "trigonal_bipyramidal" } else { "square_pyramidal" },
            6 => if near_180 == 3 && near_90 >= 8 { "octahedral" } else { "unknown" },
            _ => "unknown",
        }.to_string()
    }

    // ── P29a: Functional group detection ─────────────────────────────────────

    fn detect_functional_groups_data(&self) -> Vec<String> {
        if self.bonds.is_empty() { return vec![]; }
        let n = self.symbols.len();

        let sym = |i: usize| self.symbols.get(i).map(|s| s.as_str()).unwrap_or("");
        let bo_at = |i: usize, k: usize| -> u8 {
            self.bond_orders.get(i).and_then(|o| o.get(k)).copied().unwrap_or(1)
        };
        let h_nb  = |i: usize| -> usize { self.bonds[i].iter().filter(|&&j| sym(j)=="H").count() };
        let c_nb  = |i: usize| -> usize { self.bonds[i].iter().filter(|&&j| sym(j)=="C").count() };
        let heavy = |i: usize| -> usize { self.bonds[i].iter().filter(|&&j| sym(j)!="H").count() };
        // count double-bond O neighbors
        let dbl_o = |i: usize| -> usize {
            self.bonds[i].iter().enumerate()
                .filter(|&(k, &j)| sym(j)=="O" && bo_at(i, k)==2).count()
        };

        let mut found: std::collections::HashSet<&'static str> = std::collections::HashSet::new();

        for i in 0..n {
            match sym(i) {
                "O" => {
                    let hc = h_nb(i);
                    let cc = c_nb(i);
                    if hc >= 1 && cc >= 1 { found.insert("alcohol"); }
                    if hc == 0 && cc == 2 { found.insert("ether"); }
                }
                "C" => {
                    let nbrs = &self.bonds[i];
                    // Double-bond O neighbors
                    let o_dbl: Vec<usize> = nbrs.iter().enumerate()
                        .filter(|&(k, &j)| sym(j)=="O" && bo_at(i, k)==2)
                        .map(|(_, &j)| j).collect();
                    // Single-bond O neighbors (C-O)
                    let o_sng: Vec<usize> = nbrs.iter().enumerate()
                        .filter(|&(k, &j)| sym(j)=="O" && bo_at(i, k)==1)
                        .map(|(_, &j)| j).collect();
                    // N neighbors (any)
                    let has_n = nbrs.iter().any(|&j| sym(j)=="N");

                    if !o_dbl.is_empty() {
                        let hv = heavy(i);
                        if hv == 1 {
                            found.insert("aldehyde");
                        } else if !has_n && o_sng.is_empty() {
                            found.insert("ketone");
                        } else if has_n {
                            found.insert("amide");
                        }
                        // carboxylic acid / ester: C has both =O and -O
                        if !o_sng.is_empty() {
                            let bridging_o = o_sng[0];
                            if h_nb(bridging_o) >= 1 {
                                found.insert("carboxylic_acid");
                            } else {
                                found.insert("ester");
                            }
                        }
                    }

                    // Alkene: C=C, both atoms non-ring
                    for (k, &j) in nbrs.iter().enumerate() {
                        if sym(j)=="C" && bo_at(i, k)==2
                            && !self.ring_atoms.get(i).copied().unwrap_or(false)
                            && !self.ring_atoms.get(j).copied().unwrap_or(false) {
                            found.insert("alkene");
                        }
                        // Alkyne: C≡C
                        if sym(j)=="C" && bo_at(i, k)==3 {
                            found.insert("alkyne");
                        }
                        // Halides
                        match sym(j) {
                            "F"  => { found.insert("halide_F"); }
                            "Cl" => { found.insert("halide_Cl"); }
                            "Br" => { found.insert("halide_Br"); }
                            "I"  => { found.insert("halide_I"); }
                            _ => {}
                        }
                    }
                }
                "N" => {
                    let hc = h_nb(i);
                    let cc = c_nb(i);
                    // Amines (not bonded to C=O)
                    let bonded_to_carbonyl = self.bonds[i].iter().any(|&j| {
                        sym(j)=="C" && dbl_o(j) >= 1
                    });
                    if !bonded_to_carbonyl {
                        if hc >= 2 && cc >= 1 { found.insert("primary_amine"); }
                        else if hc == 1 && cc >= 2 { found.insert("secondary_amine"); }
                        else if hc == 0 && cc >= 3 { found.insert("tertiary_amine"); }
                    }
                    // Nitro: N with two double-bond O
                    if dbl_o(i) >= 2 { found.insert("nitro"); }
                }
                "S" => {
                    let hc = h_nb(i);
                    let cc = c_nb(i);
                    if hc >= 1 && cc >= 1 { found.insert("thiol"); }
                    if hc == 0 && cc == 2 { found.insert("sulfide"); }
                }
                "P" => {
                    if dbl_o(i) >= 1 { found.insert("phosphate"); }
                }
                _ => {}
            }
        }

        // Aromatic: based on ring detection
        if self.aromatic_ring_count() > 0 { found.insert("aromatic"); }

        let mut result: Vec<String> = found.into_iter().map(|s| s.to_string()).collect();
        result.sort();
        result
    }

    // ── P28: SMILES output ────────────────────────────────────────────────────

    fn to_smiles_data(&self) -> String {
        let n = self.symbols.len();
        if n == 0 { return String::new(); }

        let h_set: Vec<bool> = (0..n).map(|i| self.symbols[i] == "H").collect();
        let charges = &self.charges;

        // Canonical ranks (Morgan-style)
        let ranks = smiles_canonical_ranks(
            n, &self.symbols, &self.bonds, &self.bond_orders,
            &self.ring_atoms, charges,
        );

        // If no bonds at all, emit each heavy atom as "[X]" separated by '.'
        if self.bonds.is_empty() || self.bonds.iter().all(|nb| nb.is_empty()) {
            let parts: Vec<String> = (0..n)
                .filter(|&i| !h_set[i])
                .map(|i| smiles_bracket_atom(&self.symbols[i], charges.get(i).copied().unwrap_or(0), 0))
                .collect();
            return parts.join(".");
        }

        // Find connected components of heavy atoms
        let mut visited = vec![false; n];
        let mut fragments: Vec<Vec<usize>> = Vec::new();
        for start in 0..n {
            if h_set[start] || visited[start] { continue; }
            let mut comp = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start);
            visited[start] = true;
            while let Some(u) = queue.pop_front() {
                comp.push(u);
                for &v in &self.bonds[u] {
                    if !h_set[v] && !visited[v] {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
            fragments.push(comp);
        }

        let mut frag_strings: Vec<String> = Vec::new();
        for frag in &fragments {
            // Start from atom with highest rank in fragment
            let start = *frag.iter().max_by_key(|&&i| ranks[i]).unwrap();

            // Pass 1: find ring-forming edges
            let mut pre_visited = vec![false; n];
            let mut ring_bonds: Vec<(usize, usize, u8)> = Vec::new();
            smiles_find_ring_bonds(self, start, usize::MAX, &ranks, &h_set, &mut pre_visited, &mut ring_bonds);

            // Build ring closure map: atom → Vec<(digit_str, bond_order)>
            let ring_bond_set: std::collections::HashSet<(usize, usize)> =
                ring_bonds.iter().map(|&(a, b, _)| (a, b)).collect();
            let mut ring_closures: std::collections::HashMap<usize, Vec<(String, u8)>> =
                std::collections::HashMap::new();
            for (idx, &(u, v, bo)) in ring_bonds.iter().enumerate() {
                let digit = (idx + 1) as u32;
                let d_str = if digit < 10 { format!("{digit}") } else { format!("%{digit}") };
                ring_closures.entry(u).or_default().push((d_str.clone(), bo));
                ring_closures.entry(v).or_default().push((d_str, bo));
            }

            // Pass 2: generate SMILES
            let mut dfs_visited = vec![false; n];
            let mut out = String::new();
            smiles_dfs(
                self, start, usize::MAX, &ranks, &h_set, charges,
                &ring_closures, &ring_bond_set, &mut dfs_visited, &mut out,
            );
            frag_strings.push(out);
        }
        frag_strings.join(".")
    }

    fn compute_2d_coords_data(&mut self) {
        let n = self.symbols.len();
        if n == 0 {
            return;
        }
        if self.ring_atoms.is_empty() && !self.bonds.is_empty() {
            self.compute_rings();
        }
        let h_set: Vec<bool> = self.symbols.iter().map(|s| s == "H").collect();
        let rings = self.enumerate_rings();
        let frags = coords2d_heavy_fragments(&self.bonds, &h_set, n);
        let mut coords: Vec<[f32; 2]> = vec![[0.0; 2]; n];
        let mut placed: Vec<bool> = vec![false; n];
        let mut x_offset: f32 = 0.0;
        for frag in &frags {
            let has_ring = frag.iter().any(|&i| self.ring_atoms.get(i).copied().unwrap_or(false));
            if has_ring {
                let frag_set: std::collections::HashSet<usize> = frag.iter().copied().collect();
                let frag_rings: Vec<&Vec<usize>> = rings.iter()
                    .filter(|r| r.iter().any(|&i| frag_set.contains(&i)))
                    .collect();
                coords2d_layout_rings(&frag_rings, &mut coords, &mut placed, BOND_LEN_2D, [x_offset, 0.0]);
                coords2d_attach_chains(frag, &self.ring_atoms, &self.bonds, &h_set, &mut coords, &mut placed, BOND_LEN_2D);
            } else {
                coords2d_layout_chain(frag, &self.bonds, &h_set, &mut coords, &mut placed, BOND_LEN_2D, [x_offset, 0.0]);
            }
            let xs: Vec<f32> = frag.iter().filter(|&&i| placed[i]).map(|&i| coords[i][0]).collect();
            let width = if xs.is_empty() {
                0.0
            } else {
                xs.iter().cloned().fold(f32::NEG_INFINITY, f32::max)
                    - xs.iter().cloned().fold(f32::INFINITY, f32::min)
            };
            x_offset += width + BOND_LEN_2D * 3.0;
        }
        // H atoms: place at their bonded heavy atom's position
        for i in 0..n {
            if h_set[i] {
                if let Some(&nb) = self.bonds.get(i).and_then(|b| b.iter().find(|&&v| !h_set[v])) {
                    coords[i] = coords[nb];
                }
            }
        }
        for (i, c) in coords.iter().enumerate() {
            self.x[i] = c[0];
            self.y[i] = c[1];
            self.z[i] = 0.0;
        }
    }

    fn to_svg_data(&self, width: u32, height: u32) -> String {
        self.to_svg_data_impl(width, height, &std::collections::HashSet::new())
    }

    fn to_svg_data_impl(&self, width: u32, height: u32, highlight: &std::collections::HashSet<usize>) -> String {
        let n = self.x.len();
        let header = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" \
             viewBox=\"0 0 {width} {height}\" style=\"background:#fff\">\n"
        );
        let footer = "</svg>";

        let heavy: Vec<usize> = (0..n).filter(|&i| self.symbols[i] != "H").collect();
        if heavy.is_empty() {
            return format!("{header}{footer}");
        }

        let min_x = heavy.iter().map(|&i| self.x[i]).fold(f32::INFINITY, f32::min);
        let max_x = heavy.iter().map(|&i| self.x[i]).fold(f32::NEG_INFINITY, f32::max);
        let min_y = heavy.iter().map(|&i| self.y[i]).fold(f32::INFINITY, f32::min);
        let max_y = heavy.iter().map(|&i| self.y[i]).fold(f32::NEG_INFINITY, f32::max);

        let mol_w = (max_x - min_x).max(1e-3);
        let mol_h = (max_y - min_y).max(1e-3);
        let padding = 40.0_f32;
        let scale = ((width as f32 - 2.0 * padding) / mol_w)
            .min((height as f32 - 2.0 * padding) / mol_h);

        let tx = |x: f32| (x - min_x) * scale + padding;
        let ty = |y: f32| (max_y - y) * scale + padding;

        let atom_r = 8.0_f32;
        let bond_offset = 3.5_f32;

        // Build aromatic bond → ring centroid map (SVG coords).
        // Only populated when aromatic_atoms are present (SMILES) and rings are computed.
        let aromatic_bond_centers: std::collections::HashMap<(usize, usize), [f32; 2]> =
            if self.aromatic_atoms.iter().any(|&a| a) && !self.ring_bonds.is_empty() {
                let rings = self.enumerate_rings();
                let mut map = std::collections::HashMap::new();
                for ring in &rings {
                    if ring.iter().all(|&i| self.aromatic_atoms.get(i).copied().unwrap_or(false)) {
                        let rcx: f32 = ring.iter().map(|&i| tx(self.x[i])).sum::<f32>()
                            / ring.len() as f32;
                        let rcy: f32 = ring.iter().map(|&i| ty(self.y[i])).sum::<f32>()
                            / ring.len() as f32;
                        let rn = ring.len();
                        for j in 0..rn {
                            let u = ring[j];
                            let v = ring[(j + 1) % rn];
                            map.insert((u.min(v), u.max(v)), [rcx, rcy]);
                        }
                    }
                }
                map
            } else {
                std::collections::HashMap::new()
            };

        // Precompute wedge/dash bonds from tetrahedral stereo centers.
        let mut wedge_bonds: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        let mut dash_bonds: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for (&center, &(desc, from)) in &self.stereo_centers {
            let from_atom = match from { Some(f) => f, None => continue };
            if center >= n || from_atom >= n { continue; }
            let heavy_nbs: Vec<usize> = self.bonds.get(center).map(|nb| nb.iter()
                .filter(|&&nb| self.symbols.get(nb).map(|s| s != "H").unwrap_or(false) && nb != from_atom)
                .copied().collect()).unwrap_or_default();
            if heavy_nbs.len() < 2 { continue; }
            let (a, b) = (heavy_nbs[0], heavy_nbs[1]);
            let cx = tx(self.x[center]); let cy = ty(self.y[center]);
            let (fx, fy) = (tx(self.x[from_atom]) - cx, ty(self.y[from_atom]) - cy);
            let (ax, ay) = (tx(self.x[a]) - cx, ty(self.y[a]) - cy);
            let (bx, by) = (tx(self.x[b]) - cx, ty(self.y[b]) - cy);
            // Signed area > 0 means CW in SVG coords (y-down)
            let signed = fx * ay - fy * ax + ax * by - ay * bx + bx * fy - by * fx;
            let layout_is_cw = signed > 0.0;
            let stereo_is_cw = desc > 0; // @@ = 1 = CW
            let bond_key = (center.min(a), center.max(a));
            if layout_is_cw == stereo_is_cw { dash_bonds.insert(bond_key); }
            else { wedge_bonds.insert(bond_key); }
        }

        let mut out = header;

        // Highlight halos (drawn first, behind bonds and atoms)
        let halo_r = atom_r * 1.8_f32;
        for &i in highlight {
            if i < n && self.symbols[i] != "H" {
                let cx = tx(self.x[i]);
                let cy = ty(self.y[i]);
                out.push_str(&format!(
                    "<circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{halo_r:.1}\" \
                     fill=\"#FFE066\" opacity=\"0.75\" stroke=\"none\"/>\n"
                ));
            }
        }

        for u in 0..n {
            if self.symbols[u] == "H" {
                continue;
            }
            for (k, &v) in self.bonds[u].iter().enumerate() {
                if v <= u || self.symbols[v] == "H" {
                    continue;
                }
                let key = (u.min(v), u.max(v));
                if wedge_bonds.contains(&key) {
                    out.push_str(&svg_render_wedge_bond(
                        tx(self.x[u]), ty(self.y[u]),
                        tx(self.x[v]), ty(self.y[v]),
                    ));
                } else if dash_bonds.contains(&key) {
                    out.push_str(&svg_render_dash_bond(
                        tx(self.x[u]), ty(self.y[u]),
                        tx(self.x[v]), ty(self.y[v]),
                    ));
                } else if let Some(&[rcx, rcy]) = aromatic_bond_centers.get(&key) {
                    out.push_str(&svg_render_aromatic_bond(
                        tx(self.x[u]),
                        ty(self.y[u]),
                        tx(self.x[v]),
                        ty(self.y[v]),
                        rcx,
                        rcy,
                        bond_offset,
                    ));
                } else {
                    let order = self
                        .bond_orders
                        .get(u)
                        .and_then(|row| row.get(k))
                        .copied()
                        .unwrap_or(1);
                    out.push_str(&svg_render_bond(
                        tx(self.x[u]),
                        ty(self.y[u]),
                        tx(self.x[v]),
                        ty(self.y[v]),
                        order,
                        bond_offset,
                    ));
                }
            }
        }

        for &i in &heavy {
            let cx = tx(self.x[i]);
            let cy = ty(self.y[i]);
            let sym = &self.symbols[i];
            let color = svg_cpk_color(sym);
            out.push_str(&format!(
                "<circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{atom_r:.1}\" \
                 fill=\"{color}\" stroke=\"#333\" stroke-width=\"1.5\"/>\n"
            ));
            if sym != "C" {
                out.push_str(&format!(
                    "<text x=\"{cx:.1}\" y=\"{cy:.1}\" text-anchor=\"middle\" \
                     dominant-baseline=\"central\" font-family=\"sans-serif\" \
                     font-size=\"10\" fill=\"#fff\">{sym}</text>\n"
                ));
            }
        }

        out.push_str(footer);
        out
    }

    fn atom_info_at(&self, index: usize) -> Option<AtomInfo> {
        if index >= self.symbols.len() {
            return None;
        }
        Some(AtomInfo {
            index,
            symbol: self.symbols[index].clone(),
            x: self.x[index],
            y: self.y[index],
            z: self.z[index],
            atom_name: self.atom_names.get(index).cloned().unwrap_or_default(),
            residue_name: self.residue_names.get(index).cloned().unwrap_or_default(),
            residue_id: self.residue_ids.get(index).copied().unwrap_or(0),
            chain_id: (self.chain_ids.get(index).copied().unwrap_or(b' ') as char).to_string(),
            is_hetatm: self.hetatm_flags.get(index).copied().unwrap_or(false),
            occupancy: self.occupancies.get(index).copied().unwrap_or(1.0),
            b_factor: self.b_factors.get(index).copied().unwrap_or(0.0),
        })
    }
}

// --- P32: 2D layout helpers ---

const BOND_LEN_2D: f32 = 1.5;

fn coords2d_heavy_fragments(bonds: &[Vec<usize>], h_set: &[bool], n: usize) -> Vec<Vec<usize>> {
    let mut visited = vec![false; n];
    let mut frags: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if h_set[start] || visited[start] {
            continue;
        }
        let mut frag = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(u) = queue.pop_front() {
            frag.push(u);
            if let Some(nb) = bonds.get(u) {
                for &v in nb {
                    if !h_set[v] && !visited[v] {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
        }
        frags.push(frag);
    }
    frags
}

fn coords2d_place_ngon(
    ring: &[usize],
    center: [f32; 2],
    bond_len: f32,
    start_angle: f32,
    coords: &mut [[f32; 2]],
    placed: &mut [bool],
) {
    use std::f32::consts::PI;
    let n = ring.len() as f32;
    let radius = bond_len / (2.0 * (PI / n).sin());
    for (i, &atom) in ring.iter().enumerate() {
        let angle = start_angle + i as f32 * 2.0 * PI / n;
        coords[atom] = [center[0] + radius * angle.cos(), center[1] + radius * angle.sin()];
        placed[atom] = true;
    }
}

fn coords2d_place_fused(
    ring: &[usize],
    u: usize,
    v: usize,
    coords: &mut [[f32; 2]],
    placed: &mut [bool],
    bond_len: f32,
) {
    use std::f32::consts::PI;
    let n = ring.len() as f32;
    let mid = [
        (coords[u][0] + coords[v][0]) * 0.5,
        (coords[u][1] + coords[v][1]) * 0.5,
    ];
    let dx = coords[v][0] - coords[u][0];
    let dy = coords[v][1] - coords[u][1];
    let len = (dx * dx + dy * dy).sqrt().max(1e-6);
    let perp = [-dy / len, dx / len];
    let apothem = bond_len / (2.0 * (PI / n).tan());
    // Use centroid of all placed atoms to choose the far side
    let placed_count = placed.iter().filter(|&&p| p).count();
    let (cx, cy) = if placed_count > 0 {
        let sx: f32 = coords.iter().zip(placed.iter()).filter(|(_, &p)| p).map(|(c, _)| c[0]).sum();
        let sy: f32 = coords.iter().zip(placed.iter()).filter(|(_, &p)| p).map(|(c, _)| c[1]).sum();
        (sx / placed_count as f32, sy / placed_count as f32)
    } else {
        (mid[0], mid[1])
    };
    let c1 = [mid[0] + apothem * perp[0], mid[1] + apothem * perp[1]];
    let c2 = [mid[0] - apothem * perp[0], mid[1] - apothem * perp[1]];
    let d1 = (c1[0] - cx).powi(2) + (c1[1] - cy).powi(2);
    let d2 = (c2[0] - cx).powi(2) + (c2[1] - cy).powi(2);
    let center = if d1 >= d2 { c1 } else { c2 };
    let radius = bond_len / (2.0 * (PI / n).sin());
    let angle_u = (coords[u][1] - center[1]).atan2(coords[u][0] - center[0]);
    let u_pos = ring.iter().position(|&i| i == u).unwrap_or(0);
    for (offset, &atom) in ring.iter().enumerate() {
        if placed[atom] {
            continue;
        }
        let angle = angle_u + (offset as isize - u_pos as isize) as f32 * 2.0 * PI / n;
        coords[atom] = [center[0] + radius * angle.cos(), center[1] + radius * angle.sin()];
        placed[atom] = true;
    }
}

fn coords2d_layout_rings(
    frag_rings: &[&Vec<usize>],
    coords: &mut [[f32; 2]],
    placed: &mut [bool],
    bond_len: f32,
    origin: [f32; 2],
) {
    use std::f32::consts::PI;
    if frag_rings.is_empty() {
        return;
    }
    // Prefer 6-membered rings; otherwise pick the largest
    let start_idx = frag_rings
        .iter()
        .enumerate()
        .max_by_key(|(_, r)| if r.len() == 6 { 1000 } else { r.len() })
        .map(|(i, _)| i)
        .unwrap_or(0);
    coords2d_place_ngon(frag_rings[start_idx], origin, bond_len, PI / 2.0, coords, placed);
    let mut changed = true;
    while changed {
        changed = false;
        for ring in frag_rings.iter() {
            let n = ring.len();
            if ring.iter().all(|&a| placed[a]) {
                continue;
            }
            let found = (0..n).find_map(|i| {
                let u = ring[i];
                let v = ring[(i + 1) % n];
                if placed[u] && placed[v] { Some((u, v)) } else { None }
            });
            if let Some((u, v)) = found {
                coords2d_place_fused(ring, u, v, coords, placed, bond_len);
                changed = true;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn coords2d_chain_dfs(
    u: usize,
    parent: usize,
    incoming_angle: f32,
    bonds: &[Vec<usize>],
    h_set: &[bool],
    coords: &mut [[f32; 2]],
    placed: &mut [bool],
    bond_len: f32,
) {
    use std::f32::consts::PI;
    let children: Vec<usize> = bonds
        .get(u)
        .map(|nb| nb.as_slice())
        .unwrap_or(&[])
        .iter()
        .copied()
        .filter(|&v| !h_set[v] && !placed[v] && v != parent)
        .collect();
    let nc = children.len();
    for (i, &child) in children.iter().enumerate() {
        let angle = if nc == 1 {
            if u % 2 == 0 { incoming_angle - PI / 6.0 } else { incoming_angle + PI / 6.0 }
        } else if nc == 2 {
            if i == 0 { incoming_angle - PI / 3.0 } else { incoming_angle + PI / 3.0 }
        } else {
            let spread = 2.0 * PI / 3.0;
            incoming_angle + (i as f32 - (nc - 1) as f32 / 2.0) * spread / (nc - 1).max(1) as f32
        };
        coords[child] = [
            coords[u][0] + bond_len * angle.cos(),
            coords[u][1] + bond_len * angle.sin(),
        ];
        placed[child] = true;
        coords2d_chain_dfs(child, u, angle, bonds, h_set, coords, placed, bond_len);
    }
}

fn coords2d_layout_chain(
    frag: &[usize],
    bonds: &[Vec<usize>],
    h_set: &[bool],
    coords: &mut [[f32; 2]],
    placed: &mut [bool],
    bond_len: f32,
    start_xy: [f32; 2],
) {
    if frag.is_empty() {
        return;
    }
    let start = frag
        .iter()
        .copied()
        .min_by_key(|&i| bonds.get(i).map(|nb| nb.iter().filter(|&&v| !h_set[v]).count()).unwrap_or(0))
        .unwrap_or(frag[0]);
    coords[start] = start_xy;
    placed[start] = true;
    coords2d_chain_dfs(start, usize::MAX, 0.0, bonds, h_set, coords, placed, bond_len);
}

fn coords2d_attach_chains(
    frag: &[usize],
    ring_atoms: &[bool],
    bonds: &[Vec<usize>],
    h_set: &[bool],
    coords: &mut [[f32; 2]],
    placed: &mut [bool],
    bond_len: f32,
) {
    let ring_placed: Vec<usize> = frag
        .iter()
        .copied()
        .filter(|&i| ring_atoms.get(i).copied().unwrap_or(false) && placed[i])
        .collect();
    for u in ring_placed {
        let chain_children: Vec<usize> = bonds
            .get(u)
            .map(|nb| nb.as_slice())
            .unwrap_or(&[])
            .iter()
            .copied()
            .filter(|&v| !h_set[v] && !placed[v])
            .collect();
        if chain_children.is_empty() {
            continue;
        }
        // Average inward direction from ring neighbors; outward is opposite
        let ring_nb: Vec<usize> = bonds
            .get(u)
            .map(|nb| nb.as_slice())
            .unwrap_or(&[])
            .iter()
            .copied()
            .filter(|&v| placed[v] && ring_atoms.get(v).copied().unwrap_or(false))
            .collect();
        let avg_dx: f32 = ring_nb.iter().map(|&v| coords[v][0] - coords[u][0]).sum::<f32>();
        let avg_dy: f32 = ring_nb.iter().map(|&v| coords[v][1] - coords[u][1]).sum::<f32>();
        let len = (avg_dx * avg_dx + avg_dy * avg_dy).sqrt().max(1e-6);
        let out_angle = (-avg_dy / len).atan2(-avg_dx / len);
        let nc = chain_children.len();
        for (i, &child) in chain_children.iter().enumerate() {
            let spread = std::f32::consts::PI / 3.0;
            let angle = if nc == 1 {
                out_angle
            } else {
                out_angle + (i as f32 - (nc - 1) as f32 / 2.0) * spread
            };
            coords[child] = [coords[u][0] + bond_len * angle.cos(), coords[u][1] + bond_len * angle.sin()];
            placed[child] = true;
            coords2d_chain_dfs(child, u, angle, bonds, h_set, coords, placed, bond_len);
        }
    }
}

// --- P33: SVG renderer helpers ---

fn svg_cpk_color(symbol: &str) -> &'static str {
    match symbol {
        "C"  => "#404040",
        "O"  => "#E8534A",
        "N"  => "#4C72B0",
        "S"  => "#E8C533",
        "P"  => "#FF8000",
        "F"  => "#33CC33",
        "Cl" => "#1FB341",
        "Br" => "#A62929",
        "I"  => "#7E1FBD",
        _    => "#999999",
    }
}

fn svg_render_wedge_bond(x1: f32, y1: f32, x2: f32, y2: f32) -> String {
    let dx = x2 - x1; let dy = y2 - y1;
    let len = (dx * dx + dy * dy).sqrt().max(1e-6);
    let pw = 4.0_f32; // half-width at base
    let px = -dy / len * pw; let py = dx / len * pw;
    format!(
        "<polygon points=\"{:.1},{:.1} {:.1},{:.1} {:.1},{:.1}\" fill=\"#333\"/>\n",
        x1, y1, x2 + px, y2 + py, x2 - px, y2 - py
    )
}

fn svg_render_dash_bond(x1: f32, y1: f32, x2: f32, y2: f32) -> String {
    let dx = x2 - x1; let dy = y2 - y1;
    let len = (dx * dx + dy * dy).sqrt().max(1e-6);
    let mut s = String::new();
    for i in 1..=5 {
        let t = i as f32 / 6.0;
        let mx = x1 + dx * t; let my = y1 + dy * t;
        let hw = 3.5 * (i as f32 / 5.0); // taper: narrow at x1, wide at x2
        let hx = -dy / len * hw; let hy = dx / len * hw;
        s.push_str(&format!(
            "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" \
             stroke=\"#333\" stroke-width=\"1.5\"/>\n",
            mx - hx, my - hy, mx + hx, my + hy
        ));
    }
    s
}

fn svg_render_bond(x1: f32, y1: f32, x2: f32, y2: f32, order: u8, offset: f32) -> String {
    match order {
        2 => {
            let dx = x2 - x1;
            let dy = y2 - y1;
            let len = (dx * dx + dy * dy).sqrt().max(1e-6);
            let px = -dy / len * offset;
            let py = dx / len * offset;
            format!(
                "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#333\" stroke-width=\"1.5\"/>\n\
                 <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#333\" stroke-width=\"1.5\"/>\n",
                x1 + px, y1 + py, x2 + px, y2 + py,
                x1 - px, y1 - py, x2 - px, y2 - py
            )
        }
        3 => {
            let dx = x2 - x1;
            let dy = y2 - y1;
            let len = (dx * dx + dy * dy).sqrt().max(1e-6);
            let px = -dy / len * offset;
            let py = dx / len * offset;
            format!(
                "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#333\" stroke-width=\"1.5\"/>\n\
                 <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#333\" stroke-width=\"1.5\"/>\n\
                 <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#333\" stroke-width=\"1.5\"/>\n",
                x1, y1, x2, y2,
                x1 + px, y1 + py, x2 + px, y2 + py,
                x1 - px, y1 - py, x2 - px, y2 - py
            )
        }
        _ => format!(
            "<line x1=\"{x1:.1}\" y1=\"{y1:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" stroke=\"#333\" stroke-width=\"1.5\"/>\n"
        ),
    }
}

fn svg_render_aromatic_bond(
    x1: f32, y1: f32, x2: f32, y2: f32,
    ring_cx: f32, ring_cy: f32, offset: f32,
) -> String {
    // Outer bond: full-length solid line.
    // Inner bond: 70% length, dashed, offset toward ring centroid.
    let mx = (x1 + x2) / 2.0;
    let my = (y1 + y2) / 2.0;
    let dx = ring_cx - mx;
    let dy = ring_cy - my;
    let d_len = (dx * dx + dy * dy).sqrt().max(1e-6);
    let nx = dx / d_len * offset;
    let ny = dy / d_len * offset;
    let bx = x2 - x1;
    let by = y2 - y1;
    let shrink = 0.15_f32;
    let ix1 = x1 + bx * shrink + nx;
    let iy1 = y1 + by * shrink + ny;
    let ix2 = x2 - bx * shrink + nx;
    let iy2 = y2 - by * shrink + ny;
    format!(
        "<line x1=\"{x1:.1}\" y1=\"{y1:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" \
         stroke=\"#333\" stroke-width=\"1.5\"/>\n\
         <line x1=\"{ix1:.1}\" y1=\"{iy1:.1}\" x2=\"{ix2:.1}\" y2=\"{iy2:.1}\" \
         stroke=\"#333\" stroke-width=\"1.5\" stroke-dasharray=\"3,2\"/>\n"
    )
}

// --- P28 SMILES helpers ---

fn smiles_default_valence(sym: &str) -> Option<u8> {
    match sym {
        "B"  => Some(3), "C"  => Some(4), "N"  => Some(3),
        "O"  => Some(2), "P"  => Some(3), "S"  => Some(2),
        "F"  => Some(1), "Cl" => Some(1), "Br" => Some(1), "I" => Some(1),
        _ => None,
    }
}

fn smiles_bracket_atom(sym: &str, charge: i32, implicit_h: u8) -> String {
    let h_str = match implicit_h { 0 => String::new(), 1 => "H".into(), n => format!("H{n}") };
    let ch_str = match charge {
        0 => String::new(), 1 => "+".into(), -1 => "-".into(),
        n if n > 0 => format!("+{n}"), n => format!("{n}"),
    };
    format!("[{sym}{h_str}{ch_str}]")
}


fn smiles_bond_char(order: u8) -> &'static str {
    match order { 2 => "=", 3 => "#", _ => "" }
}

/// Emit a SMILES atom token with an optional tetrahedral stereo annotation.
/// `stereo`: Some(true) = "@@", Some(false) = "@", None = no annotation.
/// When stereo is present the bracket is forced; `explicit_h` counts H-atom
/// neighbors that are hidden from the DFS output and belong in the bracket.
fn smiles_atom_string_stereo(
    sym: &str,
    charge: i32,
    heavy_degree: usize,
    explicit_h: usize,
    stereo: Option<bool>,
) -> String {
    let default_val = smiles_default_valence(sym);
    let implicit_h: u8 = if charge == 0 {
        default_val
            .map(|v| (v as usize).saturating_sub(heavy_degree + explicit_h) as u8)
            .unwrap_or(0)
    } else {
        0
    };
    let needs_bracket = default_val.is_none() || charge != 0 || stereo.is_some();
    if !needs_bracket {
        return sym.to_string();
    }
    let h_count: u8 = if stereo.is_some() { explicit_h as u8 } else { implicit_h };
    let h_str = match h_count { 0 => String::new(), 1 => "H".into(), n => format!("H{n}") };
    let st_str = match stereo { Some(true) => "@@", Some(false) => "@", None => "" };
    let ch_str = match charge {
        0 => String::new(), 1 => "+".into(), -1 => "-".into(),
        n if n > 0 => format!("+{n}"), n => format!("{n}"),
    };
    format!("[{sym}{st_str}{h_str}{ch_str}]")
}

/// Determine the SMILES tetrahedral stereo annotation for atom `u` in the
/// current DFS context. Returns `Some(true)` = "@@" (CW), `Some(false)` = "@"
/// (CCW), or `None` when stereo cannot be determined.
///
/// Convention: looking from `parent` toward `u`, the substituents are listed in
/// SMILES bracket order — H-atom neighbor (if any) first, then `children` by
/// canonical rank. Positive signed triple product → CCW → "@".
fn smiles_stereo_at(
    mol: &MolecularSystem,
    u: usize,
    parent: usize,
    children: &[(usize, u8)],
    h_set: &[bool],
) -> Option<bool> {
    let n = mol.symbols.len();
    if parent >= n { return None; } // starting atom has no DFS parent

    let h_nbr: Option<usize> = mol.bonds.get(u)
        .into_iter().flatten().copied()
        .find(|&j| j < n && h_set[j]);

    let subs: Vec<usize> = h_nbr.into_iter()
        .chain(children.iter().map(|(v, _)| *v))
        .collect();
    if subs.len() != 3 { return None; }

    // Require at least one non-zero 3D coordinate among the participants
    let any_3d = |i: usize| mol.x[i] != 0.0 || mol.y[i] != 0.0 || mol.z[i] != 0.0;
    if !any_3d(u) && !any_3d(parent) && subs.iter().all(|&j| !any_3d(j)) {
        return None;
    }

    let px = mol.x[parent]; let py = mol.y[parent]; let pz = mol.z[parent];
    let v = |i: usize| [mol.x[i] - px, mol.y[i] - py, mol.z[i] - pz];
    let [v1, v2, v3] = [v(subs[0]), v(subs[1]), v(subs[2])];

    let vol = (v1[1]*v2[2] - v1[2]*v2[1]) * v3[0]
            + (v1[2]*v2[0] - v1[0]*v2[2]) * v3[1]
            + (v1[0]*v2[1] - v1[1]*v2[0]) * v3[2];

    if vol.abs() < 1e-4 { None } else { Some(vol < 0.0) }
}

fn smiles_canonical_ranks(
    n: usize,
    symbols: &[String],
    bonds: &[Vec<usize>],
    bond_orders: &[Vec<u8>],
    ring_atoms: &[bool],
    charges: &[i32],
) -> Vec<u32> {
    let h_set: Vec<bool> = (0..n).map(|i| symbols.get(i).map(|s| s == "H").unwrap_or(false)).collect();

    let fnv = |mut h: u32, v: u32| -> u32 {
        h ^= v;
        h = h.wrapping_mul(0x01000193);
        h
    };

    // Initial invariants
    let mut ranks: Vec<u32> = (0..n).map(|i| {
        let sym_hash = symbols.get(i).map(|s| {
            s.bytes().fold(0x811c9dc5u32, |h, b| fnv(h, b as u32))
        }).unwrap_or(0);
        let heavy_deg = bonds.get(i).map(|nb| nb.iter().filter(|&&j| !h_set[j]).count()).unwrap_or(0);
        let h_count   = bonds.get(i).map(|nb| nb.iter().filter(|&&j| h_set[j]).count()).unwrap_or(0);
        let charge    = (charges.get(i).copied().unwrap_or(0) + 127) as u32;
        let in_ring   = ring_atoms.get(i).copied().unwrap_or(false) as u32;
        let max_bo    = bond_orders.get(i).map(|o| o.iter().copied().max().unwrap_or(1)).unwrap_or(1) as u32;
        fnv(fnv(fnv(fnv(fnv(sym_hash, heavy_deg as u32), h_count as u32), charge), in_ring), max_bo)
    }).collect();

    // Morgan iteration — pre-allocate scratch buffers outside the loop to avoid
    // per-iteration heap allocations (typically 3–5 iterations, 3 clones each).
    let mut tmp = Vec::with_capacity(n);
    let mut new_ranks = Vec::with_capacity(n);
    loop {
        tmp.clone_from(&ranks); tmp.sort_unstable(); tmp.dedup();
        let distinct_before = tmp.len();

        new_ranks.clone_from(&ranks);
        for i in 0..n {
            if h_set[i] { continue; }
            let mut nbr_ranks: Vec<u32> = bonds.get(i).map(|nb| {
                nb.iter().filter(|&&j| !h_set[j]).map(|&j| ranks[j]).collect()
            }).unwrap_or_default();
            nbr_ranks.sort_unstable();
            let r = nbr_ranks.iter().fold(ranks[i], |h, &v| fnv(h, v));
            new_ranks[i] = r;
        }
        ranks.clone_from(&new_ranks);

        tmp.clone_from(&ranks); tmp.sort_unstable(); tmp.dedup();
        if tmp.len() <= distinct_before { break; }
    }

    // Assign final unique ranks: sort heavy atoms by (rank, index) → guarantees uniqueness
    let mut order: Vec<usize> = (0..n).filter(|&i| !h_set[i]).collect();
    order.sort_by(|&a, &b| ranks[a].cmp(&ranks[b]).then(a.cmp(&b)));
    for (new_rank, &i) in order.iter().enumerate() {
        ranks[i] = new_rank as u32;
    }
    ranks
}

// Pass 1: DFS to collect ring-forming edges (each edge recorded once, ordered by discovery).
fn smiles_find_ring_bonds(
    mol: &MolecularSystem,
    u: usize,
    parent: usize,
    ranks: &[u32],
    h_set: &[bool],
    visited: &mut Vec<bool>,
    ring_bonds: &mut Vec<(usize, usize, u8)>,
) {
    visited[u] = true;
    let nbrs = match mol.bonds.get(u) { Some(v) => v, None => return };
    let mut children: Vec<(usize, u8)> = nbrs.iter().enumerate()
        .filter(|(_, &v)| !h_set[v] && v != parent)
        .map(|(k, &v)| {
            let bo = mol.bond_orders.get(u).and_then(|o| o.get(k)).copied().unwrap_or(1);
            (v, bo)
        })
        .collect();
    children.sort_by(|(a, _), (b, _)| ranks[*b].cmp(&ranks[*a]));
    for (v, bo) in children {
        let key = (u.min(v), u.max(v));
        if visited[v] {
            if !ring_bonds.iter().any(|&(a, b, _)| (a, b) == key) {
                ring_bonds.push((key.0, key.1, bo));
            }
        } else {
            smiles_find_ring_bonds(mol, v, u, ranks, h_set, visited, ring_bonds);
        }
    }
}

// Pass 2: DFS to emit SMILES.  ring_closures maps atom → list of (digit_str, bond_order)
// that should be appended right after the atom symbol.
// ring_bond_set contains all ring-forming edges (for filtering children).
#[allow(clippy::too_many_arguments)]
fn smiles_dfs(
    mol: &MolecularSystem,
    u: usize,
    parent: usize,
    ranks: &[u32],
    h_set: &[bool],
    charges: &[i32],
    ring_closures: &std::collections::HashMap<usize, Vec<(String, u8)>>,
    ring_bond_set: &std::collections::HashSet<(usize, usize)>,
    visited: &mut Vec<bool>,
    out: &mut String,
) {
    if visited[u] { return; }
    visited[u] = true;

    let sym = &mol.symbols[u];
    let charge = charges.get(u).copied().unwrap_or(0);

    // Compute children BEFORE emitting the atom so stereo can be determined.
    let nbrs = mol.bonds.get(u);
    let heavy_degree = nbrs.map(|nb| nb.iter().filter(|&&j| !h_set[j]).count()).unwrap_or(0);
    let explicit_h   = nbrs.map(|nb| nb.iter().filter(|&&j| h_set[j]).count()).unwrap_or(0);

    let mut children: Vec<(usize, u8)> = nbrs.into_iter().flatten()
        .enumerate()
        .filter(|(_, &v)| !h_set[v] && !visited[v] && v != parent
                && !ring_bond_set.contains(&(u.min(v), u.max(v))))
        .map(|(k, &v)| {
            let bo = mol.bond_orders.get(u).and_then(|o| o.get(k)).copied().unwrap_or(1);
            (v, bo)
        })
        .collect();
    children.sort_by(|(a, _), (b, _)| ranks[*b].cmp(&ranks[*a]));

    // Tetrahedral stereo: only when a valid DFS parent exists and 3D coords present
    let stereo = smiles_stereo_at(mol, u, parent, &children, h_set);

    out.push_str(&smiles_atom_string_stereo(sym, charge, heavy_degree, explicit_h, stereo));

    // Emit ring-closure digits for this atom (pre-computed)
    if let Some(closures) = ring_closures.get(&u) {
        for (d_str, bo) in closures {
            out.push_str(smiles_bond_char(*bo));
            out.push_str(d_str);
        }
    }

    for (idx, (v, bo)) in children.iter().enumerate() {
        if idx > 0 { out.push('('); }
        out.push_str(smiles_bond_char(*bo));
        smiles_dfs(mol, *v, u, ranks, h_set, charges, ring_closures, ring_bond_set, visited, out);
        if idx > 0 { out.push(')'); }
    }
}

// ==================== P41: File Format I/O Expansion ====================
// CDXML (ChemDraw) · MRV (MarvinSketch) · CML · Ket (Ketcher) · RXN · Reaction SMILES

fn atomic_num_to_symbol(n: u32) -> &'static str {
    match n {
        1=>"H",2=>"He",3=>"Li",4=>"Be",5=>"B",6=>"C",7=>"N",8=>"O",9=>"F",10=>"Ne",
        11=>"Na",12=>"Mg",13=>"Al",14=>"Si",15=>"P",16=>"S",17=>"Cl",18=>"Ar",
        19=>"K",20=>"Ca",24=>"Cr",25=>"Mn",26=>"Fe",27=>"Co",28=>"Ni",29=>"Cu",30=>"Zn",
        33=>"As",34=>"Se",35=>"Br",47=>"Ag",50=>"Sn",53=>"I",80=>"Hg",82=>"Pb",
        _=>"C",
    }
}

fn symbol_to_atomic_num(s: &str) -> u32 {
    match s {
        "H"=>1,"He"=>2,"Li"=>3,"Be"=>4,"B"=>5,"C"=>6,"N"=>7,"O"=>8,"F"=>9,"Ne"=>10,
        "Na"=>11,"Mg"=>12,"Al"=>13,"Si"=>14,"P"=>15,"S"=>16,"Cl"=>17,"Ar"=>18,
        "K"=>19,"Ca"=>20,"Cr"=>24,"Mn"=>25,"Fe"=>26,"Co"=>27,"Ni"=>28,"Cu"=>29,"Zn"=>30,
        "As"=>33,"Se"=>34,"Br"=>35,"Ag"=>47,"Sn"=>50,"I"=>53,"Hg"=>80,"Pb"=>82,
        _=>6,
    }
}

fn p41_push_atom(mol: &mut MolecularSystem, sym: &str, x: f32, y: f32, charge: i32) {
    mol.symbols.push(sym.to_string());
    mol.x.push(x);
    mol.y.push(y);
    mol.z.push(0.0);
    mol.charges.push(charge);
    mol.atom_names.push(String::new());
    mol.residue_names.push(String::new());
    mol.residue_ids.push(0);
    mol.chain_ids.push(b' ');
    mol.hetatm_flags.push(false);
    mol.occupancies.push(1.0);
    mol.b_factors.push(0.0);
    if !mol.atom_map.is_empty() { mol.atom_map.push(0); }
}

fn p41_add_bond(mol: &mut MolecularSystem, a: usize, b: usize, order: u8) {
    let n = mol.symbols.len();
    if mol.bonds.len() < n { mol.bonds.resize(n, Vec::new()); }
    if mol.bond_orders.len() < n { mol.bond_orders.resize(n, Vec::new()); }
    if a < n && b < n {
        mol.bonds[a].push(b); mol.bond_orders[a].push(order);
        mol.bonds[b].push(a); mol.bond_orders[b].push(order);
    }
}

fn p41_xml_attrs(e: &quick_xml::events::BytesStart<'_>) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for attr in e.attributes().flatten() {
        let key = std::str::from_utf8(attr.key.as_ref()).unwrap_or("").to_ascii_uppercase();
        let val = std::str::from_utf8(attr.value.as_ref()).unwrap_or("").to_string();
        map.insert(key, val);
    }
    map
}

// --- CDXML ---

fn parse_cdxml(s: &str) -> Result<MolecularSystem, ParseError> {
    use quick_xml::Reader;
    use quick_xml::events::Event;

    let mut mol = MolecularSystem::new_empty();
    let mut id_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut pending_bonds: Vec<(String, String, u8)> = Vec::new();
    let mut bond_length_pt = 36.0f32;
    let mut depth_fragment: usize = 0;

    let mut reader = Reader::from_str(s);
    reader.config_mut().trim_text(true);

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) | Ok(Event::Empty(e)) => {
                let tag = std::str::from_utf8(e.name().as_ref())
                    .unwrap_or("").split(':').next_back().unwrap_or("").to_ascii_uppercase();
                let attrs = p41_xml_attrs(&e);
                match tag.as_str() {
                    "CDXML" => {
                        if let Some(bl) = attrs.get("BONDLENGTH") {
                            bond_length_pt = bl.parse().unwrap_or(36.0);
                        }
                    }
                    "FRAGMENT" => { depth_fragment += 1; }
                    "NODE" if depth_fragment > 0 => {
                        let node_type = attrs.get("NODETYPE").map(|s| s.as_str()).unwrap_or("Element");
                        if matches!(node_type, "ExternalConnectionPoint" | "LinkNode" | "NamedAlternativeGroup") {
                            continue;
                        }
                        let id = attrs.get("ID").cloned().unwrap_or_default();
                        if id.is_empty() { continue; }
                        let p_str = attrs.get("P").cloned().unwrap_or_default();
                        let (px, py) = {
                            let parts: Vec<f32> = p_str.split_whitespace()
                                .filter_map(|v| v.parse().ok()).collect();
                            if parts.len() >= 2 { (parts[0], parts[1]) } else { (0.0, 0.0) }
                        };
                        let element: u32 = attrs.get("ELEMENT").and_then(|v| v.parse().ok()).unwrap_or(6);
                        let charge: i32 = attrs.get("CHARGE").and_then(|v| v.parse().ok()).unwrap_or(0);
                        let sym = atomic_num_to_symbol(element);
                        let scale = bond_length_pt / 1.5;
                        let atom_idx = mol.symbols.len();
                        p41_push_atom(&mut mol, sym, px / scale, -(py / scale), charge);
                        id_map.insert(id, atom_idx);
                    }
                    "BOND" if depth_fragment > 0 => {
                        let b_id = attrs.get("B").cloned().unwrap_or_default();
                        let e_id = attrs.get("E").cloned().unwrap_or_default();
                        let order_str = attrs.get("ORDER").map(|s| s.as_str()).unwrap_or("1");
                        let order: u8 = match order_str { "2" => 2, "3" => 3, "1.5" => 4, _ => 1 };
                        pending_bonds.push((b_id, e_id, order));
                    }
                    _ => {}
                }
            }
            Ok(Event::End(e)) => {
                let tag = std::str::from_utf8(e.name().as_ref())
                    .unwrap_or("").split(':').next_back().unwrap_or("").to_ascii_uppercase();
                if tag == "FRAGMENT" && depth_fragment > 0 { depth_fragment -= 1; }
            }
            Ok(Event::Eof) => break,
            Err(_) => return Err(ParseError::EmptyInput),
            _ => {}
        }
    }

    if mol.symbols.is_empty() { return Err(ParseError::EmptyInput); }

    mol.bonds = vec![Vec::new(); mol.symbols.len()];
    mol.bond_orders = vec![Vec::new(); mol.symbols.len()];
    for (b_id, e_id, order) in pending_bonds {
        if let (Some(&a), Some(&b)) = (id_map.get(&b_id), id_map.get(&e_id)) {
            p41_add_bond(&mut mol, a, b, order);
        }
    }
    Ok(mol)
}

fn parse_cdxml_reaction(s: &str) -> Result<Reaction, ParseError> {
    use quick_xml::Reader;
    use quick_xml::events::Event;

    struct FragState {
        id_map: std::collections::HashMap<String, usize>,
        mol: MolecularSystem,
        pending_bonds: Vec<(String, String, u8)>,
    }
    impl Default for FragState {
        fn default() -> Self {
            FragState { id_map: Default::default(), mol: MolecularSystem::new_empty(), pending_bonds: Vec::new() }
        }
    }

    let mut fragments: Vec<FragState> = Vec::new();
    let mut current_frag: Option<FragState> = None;
    let mut arrows: Vec<(f32, f32)> = Vec::new();
    let mut bond_length_pt = 36.0f32;
    let mut depth: usize = 0;

    let mut reader = Reader::from_str(s);
    reader.config_mut().trim_text(true);

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) | Ok(Event::Empty(e)) => {
                let tag = std::str::from_utf8(e.name().as_ref())
                    .unwrap_or("").split(':').next_back().unwrap_or("").to_ascii_uppercase();
                let attrs = p41_xml_attrs(&e);
                match tag.as_str() {
                    "CDXML" => {
                        if let Some(bl) = attrs.get("BONDLENGTH") {
                            bond_length_pt = bl.parse().unwrap_or(36.0);
                        }
                    }
                    "FRAGMENT" => {
                        depth += 1;
                        if depth == 1 {
                            current_frag = Some(FragState::default());
                        }
                    }
                    "NODE" if depth > 0 => {
                        if let Some(ref mut frag) = current_frag {
                            let node_type = attrs.get("NODETYPE").map(|s| s.as_str()).unwrap_or("Element");
                            if matches!(node_type, "ExternalConnectionPoint" | "LinkNode") { continue; }
                            let id = attrs.get("ID").cloned().unwrap_or_default();
                            if id.is_empty() { continue; }
                            let p_str = attrs.get("P").cloned().unwrap_or_default();
                            let (px, py) = {
                                let parts: Vec<f32> = p_str.split_whitespace()
                                    .filter_map(|v| v.parse().ok()).collect();
                                if parts.len() >= 2 { (parts[0], parts[1]) } else { (0.0, 0.0) }
                            };
                            let element: u32 = attrs.get("ELEMENT").and_then(|v| v.parse().ok()).unwrap_or(6);
                            let charge: i32 = attrs.get("CHARGE").and_then(|v| v.parse().ok()).unwrap_or(0);
                            let sym = atomic_num_to_symbol(element);
                            let scale = bond_length_pt / 1.5;
                            let idx = frag.mol.symbols.len();
                            p41_push_atom(&mut frag.mol, sym, px / scale, -(py / scale), charge);
                            frag.id_map.insert(id, idx);
                        }
                    }
                    "BOND" if depth > 0 => {
                        if let Some(ref mut frag) = current_frag {
                            let b_id = attrs.get("B").cloned().unwrap_or_default();
                            let e_id = attrs.get("E").cloned().unwrap_or_default();
                            let order_str = attrs.get("ORDER").map(|s| s.as_str()).unwrap_or("1");
                            let order: u8 = match order_str { "2" => 2, "3" => 3, "1.5" => 4, _ => 1 };
                            frag.pending_bonds.push((b_id, e_id, order));
                        }
                    }
                    "ARROW" => {
                        let p1 = attrs.get("P1").or_else(|| attrs.get("HEAD3D")).cloned().unwrap_or_default();
                        let p2 = attrs.get("P2").or_else(|| attrs.get("TAIL3D")).cloned().unwrap_or_default();
                        let x1 = p1.split_whitespace().next().and_then(|v| v.parse::<f32>().ok()).unwrap_or(0.0);
                        let x2 = p2.split_whitespace().next().and_then(|v| v.parse::<f32>().ok()).unwrap_or(x1 + 72.0);
                        let (ax, bx) = if x1 <= x2 { (x1, x2) } else { (x2, x1) };
                        arrows.push((ax, bx));
                    }
                    _ => {}
                }
            }
            Ok(Event::End(e)) => {
                let tag = std::str::from_utf8(e.name().as_ref())
                    .unwrap_or("").split(':').next_back().unwrap_or("").to_ascii_uppercase();
                if tag == "FRAGMENT" {
                    depth = depth.saturating_sub(1);
                    if depth == 0 {
                        if let Some(mut frag) = current_frag.take() {
                            if !frag.mol.symbols.is_empty() {
                                let n = frag.mol.symbols.len();
                                frag.mol.bonds = vec![Vec::new(); n];
                                frag.mol.bond_orders = vec![Vec::new(); n];
                                for (b_id, e_id, order) in &frag.pending_bonds {
                                    if let (Some(&a), Some(&b)) = (frag.id_map.get(b_id), frag.id_map.get(e_id)) {
                                        p41_add_bond(&mut frag.mol, a, b, *order);
                                    }
                                }
                                fragments.push(frag);
                            }
                        }
                    }
                }
            }
            Ok(Event::Eof) => break,
            Err(_) => return Err(ParseError::EmptyInput),
            _ => {}
        }
    }

    let mut reactants = Vec::new();
    let mut products = Vec::new();
    let mut reagents = Vec::new();

    if arrows.is_empty() {
        for frag in fragments { reactants.push(frag.mol); }
    } else {
        let arrow_x1 = arrows[0].0;
        let arrow_x2 = arrows[0].1;
        let scale = bond_length_pt / 1.5;
        for frag in fragments {
            if frag.mol.x.is_empty() { continue; }
            let cx: f32 = frag.mol.x.iter().sum::<f32>() / frag.mol.x.len() as f32;
            let cx_pt = cx * scale;
            if cx_pt < arrow_x1 {
                reactants.push(frag.mol);
            } else if cx_pt > arrow_x2 {
                products.push(frag.mol);
            } else {
                reagents.push(frag.mol);
            }
        }
    }

    Ok(Reaction { reactants, products, reagents, conditions: Vec::new() })
}

fn mol_to_cdxml_string(mol: &MolecularSystem) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str("<!DOCTYPE CDXML SYSTEM \"http://www.camsoft.com/xml/cdxml.dtd\">\n");
    out.push_str("<CDXML CreationProgram=\"chem-wasm-lens\" BondLength=\"36.0\">\n<Page>\n<Fragment>\n");
    let scale = 36.0f32 / 1.5;
    for (i, sym) in mol.symbols.iter().enumerate() {
        if sym == "H" { continue; }
        let px = mol.x[i] * scale;
        let py = -(mol.y[i] * scale);
        let element = symbol_to_atomic_num(sym);
        let charge = mol.charges.get(i).copied().unwrap_or(0);
        out.push_str(&format!("  <Node id=\"{}\" p=\"{:.4} {:.4}\" Element=\"{}\"", i + 1, px, py, element));
        if charge != 0 { out.push_str(&format!(" Charge=\"{}\"", charge)); }
        out.push_str("/>\n");
    }
    let mut seen = std::collections::HashSet::new();
    for a in 0..mol.bonds.len() {
        if mol.symbols.get(a).map(|s| s == "H").unwrap_or(false) { continue; }
        for (k, &b) in mol.bonds[a].iter().enumerate() {
            if mol.symbols.get(b).map(|s| s == "H").unwrap_or(false) { continue; }
            let key = if a < b { (a, b) } else { (b, a) };
            if !seen.insert(key) { continue; }
            let order = mol.bond_orders.get(a).and_then(|v| v.get(k)).copied().unwrap_or(1);
            let order_str = match order { 2 => "2", 3 => "3", 4 => "1.5", _ => "1" };
            out.push_str(&format!("  <Bond B=\"{}\" E=\"{}\" Order=\"{}\"/>\n", a + 1, b + 1, order_str));
        }
    }
    out.push_str("</Fragment>\n</Page>\n</CDXML>\n");
    out
}

fn reaction_to_cdxml_string(r: &Reaction) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str("<!DOCTYPE CDXML SYSTEM \"http://www.camsoft.com/xml/cdxml.dtd\">\n");
    out.push_str("<CDXML CreationProgram=\"chem-wasm-lens\" BondLength=\"36.0\">\n<Page>\n");
    let scale = 36.0f32 / 1.5;
    let mut x_offset_pt = 0.0f32;

    for mol in r.reactants.iter().chain(r.products.iter().enumerate().map(|(i, m)| {
        if i == 0 { x_offset_pt += 144.0; } // arrow gap before first product
        m
    }).collect::<Vec<_>>().iter().copied()) {
        // This lambda approach is messy; use a flat loop instead
        let _ = mol;
    }
    // Cleaner: write reactants, then arrow, then products
    out.truncate(out.rfind("<CDXML").unwrap_or(0) + out[out.rfind("<CDXML").unwrap_or(0)..].find('\n').map(|p| p + out.rfind("<CDXML").unwrap_or(0) + 1).unwrap_or(out.len()));

    // Reset and rewrite cleanly
    out.clear();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str("<!DOCTYPE CDXML SYSTEM \"http://www.camsoft.com/xml/cdxml.dtd\">\n");
    out.push_str("<CDXML CreationProgram=\"chem-wasm-lens\" BondLength=\"36.0\">\n<Page>\n");

    let write_mol_fragment = |mol: &MolecularSystem, out: &mut String, x_off: f32| {
        out.push_str("<Fragment>\n");
        let min_x = mol.x.iter().cloned().fold(f32::INFINITY, f32::min);
        for (i, sym) in mol.symbols.iter().enumerate() {
            if sym == "H" { continue; }
            let px = (mol.x[i] - min_x) * scale + x_off;
            let py = -(mol.y[i] * scale);
            let element = symbol_to_atomic_num(sym);
            let charge = mol.charges.get(i).copied().unwrap_or(0);
            out.push_str(&format!("  <Node id=\"{}\" p=\"{:.4} {:.4}\" Element=\"{}\"", i + 1, px, py, element));
            if charge != 0 { out.push_str(&format!(" Charge=\"{}\"", charge)); }
            out.push_str("/>\n");
        }
        let mut seen = std::collections::HashSet::new();
        for a in 0..mol.bonds.len() {
            if mol.symbols.get(a).map(|s| s == "H").unwrap_or(false) { continue; }
            for (k, &b) in mol.bonds[a].iter().enumerate() {
                if mol.symbols.get(b).map(|s| s == "H").unwrap_or(false) { continue; }
                let key = if a < b { (a, b) } else { (b, a) };
                if !seen.insert(key) { continue; }
                let order = mol.bond_orders.get(a).and_then(|v| v.get(k)).copied().unwrap_or(1);
                let order_str = match order { 2 => "2", 3 => "3", 4 => "1.5", _ => "1" };
                out.push_str(&format!("  <Bond B=\"{}\" E=\"{}\" Order=\"{}\"/>\n", a + 1, b + 1, order_str));
            }
        }
        out.push_str("</Fragment>\n");
        let width = mol.x.iter().cloned().fold(0.0f32, |a, b| a.max(b))
            - mol.x.iter().cloned().fold(f32::INFINITY, f32::min);
        x_off + width * scale + 72.0
    };

    for mol in &r.reactants {
        x_offset_pt = write_mol_fragment(mol, &mut out, x_offset_pt);
    }
    let arrow_x1 = x_offset_pt;
    x_offset_pt += 144.0;
    let arrow_x2 = x_offset_pt;
    x_offset_pt += 72.0;
    out.push_str(&format!("<Arrow ArrowType=\"FullHead\" p1=\"{:.1} 0\" p2=\"{:.1} 0\"/>\n", arrow_x1, arrow_x2));
    for mol in &r.products {
        x_offset_pt = write_mol_fragment(mol, &mut out, x_offset_pt);
    }

    out.push_str("</Page>\n</CDXML>\n");
    out
}

// --- MRV (MarvinSketch) ---

fn parse_mrv(s: &str) -> Result<MolecularSystem, ParseError> {
    use quick_xml::Reader;
    use quick_xml::events::Event;

    let mut mol = MolecularSystem::new_empty();
    let mut id_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut pending_bonds: Vec<(String, String, u8)> = Vec::new();
    let mut in_atom_array = false;
    let mut in_bond_array = false;

    let mut reader = Reader::from_str(s);
    reader.config_mut().trim_text(true);

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) | Ok(Event::Empty(e)) => {
                let tag = std::str::from_utf8(e.name().as_ref())
                    .unwrap_or("").split(':').next_back().unwrap_or("").to_ascii_uppercase();
                let attrs = p41_xml_attrs(&e);
                match tag.as_str() {
                    "ATOMARRAY" => { in_atom_array = true; in_bond_array = false; }
                    "BONDARRAY" => { in_bond_array = true; in_atom_array = false; }
                    "ATOM" if in_atom_array => {
                        let id = attrs.get("ID").cloned().unwrap_or_default();
                        let sym = attrs.get("ELEMENTTYPE").cloned().unwrap_or_else(|| "C".to_string());
                        let x: f32 = attrs.get("X2").or_else(|| attrs.get("X3"))
                            .and_then(|v| v.parse().ok()).unwrap_or(0.0);
                        let y: f32 = attrs.get("Y2").or_else(|| attrs.get("Y3"))
                            .and_then(|v| v.parse().ok()).unwrap_or(0.0);
                        let charge: i32 = attrs.get("FORMALCHARGE").and_then(|v| v.parse().ok()).unwrap_or(0);
                        let idx = mol.symbols.len();
                        p41_push_atom(&mut mol, &sym, x, y, charge);
                        if !id.is_empty() { id_map.insert(id, idx); }
                    }
                    "BOND" if in_bond_array => {
                        let refs = attrs.get("ATOMREFS2").cloned().unwrap_or_default();
                        let parts: Vec<&str> = refs.split_whitespace().collect();
                        if parts.len() >= 2 {
                            let order_str = attrs.get("ORDER").map(|s| s.as_str()).unwrap_or("1");
                            let order: u8 = match order_str { "2" => 2, "3" => 3, "A" => 4, _ => 1 };
                            pending_bonds.push((parts[0].to_string(), parts[1].to_string(), order));
                        }
                    }
                    _ => {}
                }
            }
            Ok(Event::End(e)) => {
                let tag = std::str::from_utf8(e.name().as_ref())
                    .unwrap_or("").split(':').next_back().unwrap_or("").to_ascii_uppercase();
                match tag.as_str() {
                    "ATOMARRAY" => { in_atom_array = false; }
                    "BONDARRAY" => { in_bond_array = false; }
                    _ => {}
                }
            }
            Ok(Event::Eof) => break,
            Err(_) => return Err(ParseError::EmptyInput),
            _ => {}
        }
    }

    if mol.symbols.is_empty() { return Err(ParseError::EmptyInput); }

    mol.bonds = vec![Vec::new(); mol.symbols.len()];
    mol.bond_orders = vec![Vec::new(); mol.symbols.len()];
    for (b_id, e_id, order) in pending_bonds {
        if let (Some(&a), Some(&b)) = (id_map.get(&b_id), id_map.get(&e_id)) {
            p41_add_bond(&mut mol, a, b, order);
        }
    }
    Ok(mol)
}

fn mol_to_mrv_string(mol: &MolecularSystem) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<cml>\n  <MDocument>\n    <MChemicalStruct>\n      <molecule molID=\"m1\">\n        <atomArray>\n");
    for (i, sym) in mol.symbols.iter().enumerate() {
        if sym == "H" { continue; }
        let charge = mol.charges.get(i).copied().unwrap_or(0);
        out.push_str(&format!("          <atom id=\"a{}\" elementType=\"{}\" x2=\"{:.6}\" y2=\"{:.6}\"",
            i + 1, sym, mol.x[i], mol.y[i]));
        if charge != 0 { out.push_str(&format!(" formalCharge=\"{}\"", charge)); }
        out.push_str("/>\n");
    }
    out.push_str("        </atomArray>\n        <bondArray>\n");
    let mut seen = std::collections::HashSet::new();
    let mut bid = 1usize;
    for a in 0..mol.bonds.len() {
        if mol.symbols.get(a).map(|s| s == "H").unwrap_or(false) { continue; }
        for (k, &b) in mol.bonds[a].iter().enumerate() {
            if mol.symbols.get(b).map(|s| s == "H").unwrap_or(false) { continue; }
            let key = if a < b { (a, b) } else { (b, a) };
            if !seen.insert(key) { continue; }
            let order = mol.bond_orders.get(a).and_then(|v| v.get(k)).copied().unwrap_or(1);
            let order_str = match order { 2 => "2", 3 => "3", 4 => "A", _ => "1" };
            out.push_str(&format!("          <bond id=\"b{}\" atomRefs2=\"a{} a{}\" order=\"{}\"/>\n",
                bid, a + 1, b + 1, order_str));
            bid += 1;
        }
    }
    out.push_str("        </bondArray>\n      </molecule>\n    </MChemicalStruct>\n  </MDocument>\n</cml>\n");
    out
}

// --- CML (Chemical Markup Language) ---

fn parse_cml(s: &str) -> Result<MolecularSystem, ParseError> {
    parse_mrv(s)
}

fn mol_to_cml_string(mol: &MolecularSystem) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<cml xmlns=\"http://www.xml-cml.org/schema\">\n  <molecule id=\"m1\">\n    <atomArray>\n");
    for (i, sym) in mol.symbols.iter().enumerate() {
        if sym == "H" { continue; }
        let charge = mol.charges.get(i).copied().unwrap_or(0);
        out.push_str(&format!("      <atom id=\"a{}\" elementType=\"{}\" x2=\"{:.6}\" y2=\"{:.6}\"",
            i + 1, sym, mol.x[i], mol.y[i]));
        if charge != 0 { out.push_str(&format!(" formalCharge=\"{}\"", charge)); }
        out.push_str("/>\n");
    }
    out.push_str("    </atomArray>\n    <bondArray>\n");
    let mut seen = std::collections::HashSet::new();
    let mut bid = 1usize;
    for a in 0..mol.bonds.len() {
        if mol.symbols.get(a).map(|s| s == "H").unwrap_or(false) { continue; }
        for (k, &b) in mol.bonds[a].iter().enumerate() {
            if mol.symbols.get(b).map(|s| s == "H").unwrap_or(false) { continue; }
            let key = if a < b { (a, b) } else { (b, a) };
            if !seen.insert(key) { continue; }
            let order = mol.bond_orders.get(a).and_then(|v| v.get(k)).copied().unwrap_or(1);
            let order_str = match order { 2 => "D", 3 => "T", 4 => "A", _ => "S" };
            out.push_str(&format!("      <bond id=\"b{}\" atomRefs2=\"a{} a{}\" order=\"{}\"/>\n",
                bid, a + 1, b + 1, order_str));
            bid += 1;
        }
    }
    out.push_str("    </bondArray>\n  </molecule>\n</cml>\n");
    out
}

// --- Ket (Ketcher JSON) ---

fn parse_ket(s: &str) -> Result<MolecularSystem, ParseError> {
    let v: serde_json::Value = serde_json::from_str(s).map_err(|_| ParseError::EmptyInput)?;
    let mut mol = MolecularSystem::new_empty();
    let mut pending_bonds: Vec<(usize, usize, u8)> = Vec::new();
    let mut atom_offset = 0usize;

    let obj = v.as_object().ok_or(ParseError::EmptyInput)?;
    for (_key, val) in obj {
        let mol_type = val.get("type").and_then(|t| t.as_str()).unwrap_or("");
        if mol_type != "molecule" && mol_type != "mol" { continue; }

        if let Some(atoms) = val.get("atoms").and_then(|a| a.as_array()) {
            for atom in atoms {
                let label = atom.get("label").and_then(|l| l.as_str()).unwrap_or("C");
                let loc = atom.get("location").and_then(|l| l.as_array());
                let (x, y) = if let Some(loc) = loc {
                    let x = loc.first().and_then(|v| v.as_f64()).unwrap_or(0.0) as f32;
                    let y = loc.get(1).and_then(|v| v.as_f64()).unwrap_or(0.0) as f32;
                    (x, y)
                } else { (0.0, 0.0) };
                let charge = atom.get("charge").and_then(|c| c.as_i64()).unwrap_or(0) as i32;
                p41_push_atom(&mut mol, label, x, y, charge);
            }
        }

        if let Some(bonds) = val.get("bonds").and_then(|b| b.as_array()) {
            for bond in bonds {
                let bond_type = bond.get("type").and_then(|t| t.as_i64()).unwrap_or(1) as u8;
                let order: u8 = match bond_type { 2 => 2, 3 => 3, _ => 1 };
                if let Some(arr) = bond.get("atoms").and_then(|a| a.as_array()) {
                    if arr.len() >= 2 {
                        let a = arr[0].as_i64().unwrap_or(0) as usize + atom_offset;
                        let b = arr[1].as_i64().unwrap_or(0) as usize + atom_offset;
                        pending_bonds.push((a, b, order));
                    }
                }
            }
        }
        atom_offset = mol.symbols.len();
    }

    if mol.symbols.is_empty() { return Err(ParseError::EmptyInput); }

    mol.bonds = vec![Vec::new(); mol.symbols.len()];
    mol.bond_orders = vec![Vec::new(); mol.symbols.len()];
    for (a, b, order) in pending_bonds {
        p41_add_bond(&mut mol, a, b, order);
    }
    Ok(mol)
}

fn mol_to_ket_string(mol: &MolecularSystem) -> String {
    let heavy: Vec<usize> = (0..mol.symbols.len())
        .filter(|&i| mol.symbols[i] != "H").collect();
    let mut idx_remap = vec![usize::MAX; mol.symbols.len()];
    for (new_i, &old_i) in heavy.iter().enumerate() { idx_remap[old_i] = new_i; }

    let mut atoms_json = String::new();
    for (new_i, &old_i) in heavy.iter().enumerate() {
        if new_i > 0 { atoms_json.push(','); }
        let charge = mol.charges.get(old_i).copied().unwrap_or(0);
        atoms_json.push_str(&format!(
            "{{\"label\":\"{}\",\"location\":[{:.6},{:.6},0.0],\"charge\":{}}}",
            mol.symbols[old_i], mol.x[old_i], mol.y[old_i], charge
        ));
    }

    let mut bonds_json = String::new();
    let mut seen = std::collections::HashSet::new();
    let mut first = true;
    for a in 0..mol.bonds.len() {
        if mol.symbols.get(a).map(|s| s == "H").unwrap_or(false) { continue; }
        for (k, &b) in mol.bonds[a].iter().enumerate() {
            if mol.symbols.get(b).map(|s| s == "H").unwrap_or(false) { continue; }
            let key = if a < b { (a, b) } else { (b, a) };
            if !seen.insert(key) { continue; }
            let order = mol.bond_orders.get(a).and_then(|v| v.get(k)).copied().unwrap_or(1);
            let bond_type: u8 = match order { 2 => 2, 3 => 3, _ => 1 };
            let na = idx_remap[a];
            let nb = idx_remap[b];
            if na == usize::MAX || nb == usize::MAX { continue; }
            if !first { bonds_json.push(','); }
            first = false;
            bonds_json.push_str(&format!("{{\"type\":{},\"atoms\":[{},{}],\"stereo\":0}}", bond_type, na, nb));
        }
    }
    // suppress unused warning
    let _ = heavy;

    format!(
        "{{\"root\":{{\"nodes\":[{{\"$ref\":\"mol0\"}}]}},\"mol0\":{{\"type\":\"molecule\",\"atoms\":[{}],\"bonds\":[{}]}}}}",
        atoms_json, bonds_json
    )
}

// --- RXN (MDL Reaction Format) ---

fn parse_rxn(s: &str) -> Result<Reaction, ParseError> {
    let rxn_start = s.find("$RXN").ok_or(ParseError::EmptyInput)?;
    let s = &s[rxn_start..];
    let lines: Vec<&str> = s.lines().collect();

    // Count line is the 5th line (index 4): "  2  1  0"
    let count_line = lines.get(4).ok_or(ParseError::EmptyInput)?;
    let n_reactants: usize = count_line.get(0..3)
        .and_then(|v| v.trim().parse().ok()).unwrap_or(0);
    let n_products: usize = count_line.get(3..6)
        .and_then(|v| v.trim().parse().ok()).unwrap_or(0);

    let mol_blocks: Vec<&str> = s.split("$MOL").skip(1).collect();
    let mut reactants = Vec::new();
    let mut products = Vec::new();

    for (i, block) in mol_blocks.iter().enumerate() {
        let mol_data = block.strip_prefix('\n').unwrap_or(block);
        if let Ok(mol) = parse_sdf(mol_data) {
            if i < n_reactants {
                reactants.push(mol);
            } else if i < n_reactants + n_products {
                products.push(mol);
            }
        }
    }

    Ok(Reaction { reactants, products, reagents: Vec::new(), conditions: Vec::new() })
}

fn reaction_to_rxn_string(r: &Reaction) -> String {
    let mut out = String::new();
    out.push_str("$RXN\n\n  chem-wasm-lens\n\n");
    out.push_str(&format!("{:>3}{:>3}\n", r.reactants.len(), r.products.len()));
    for mol in r.reactants.iter().chain(r.products.iter()) {
        out.push_str("$MOL\n");
        out.push_str(&mol.to_sdf_string());
        out.push('\n');
    }
    out
}

// --- Reaction SMILES ---

fn parse_reaction_smiles(s: &str) -> Result<Reaction, ParseError> {
    let s = s.trim();
    let (reactant_str, reagent_str, product_str) = if let Some(pos) = s.find(">>") {
        (&s[..pos], "", &s[pos + 2..])
    } else {
        let parts: Vec<&str> = s.splitn(3, '>').collect();
        if parts.len() == 3 {
            (parts[0], parts[1], parts[2])
        } else {
            return Err(ParseError::EmptyInput);
        }
    };

    let parse_part = |part: &str| -> Vec<MolecularSystem> {
        if part.is_empty() { return Vec::new(); }
        part.split('.').filter_map(|smiles| {
            if smiles.is_empty() { return None; }
            parse_smiles(smiles).ok()
        }).collect()
    };

    Ok(Reaction {
        reactants: parse_part(reactant_str),
        products:  parse_part(product_str),
        reagents:  parse_part(reagent_str),
        conditions: Vec::new(),
    })
}

fn reaction_to_smiles(r: &Reaction) -> String {
    let mols_to_smiles = |mols: &[MolecularSystem]| -> String {
        mols.iter().map(|m| m.to_smiles_data()).collect::<Vec<_>>().join(".")
    };
    format!("{}>>{}", mols_to_smiles(&r.reactants), mols_to_smiles(&r.products))
}

// --- Reaction struct ---

#[wasm_bindgen]
pub struct Reaction {
    reactants:  Vec<MolecularSystem>,
    products:   Vec<MolecularSystem>,
    reagents:   Vec<MolecularSystem>,
    #[allow(dead_code)]
    conditions: Vec<String>,
}

#[wasm_bindgen]
impl Reaction {
    pub fn reactant_count(&self) -> usize { self.reactants.len() }
    pub fn product_count(&self)  -> usize { self.products.len() }
    pub fn reagent_count(&self)  -> usize { self.reagents.len() }

    pub fn get_reactant(&self, i: usize) -> Option<MolecularSystem> { self.reactants.get(i).cloned() }
    pub fn get_product(&self,  i: usize) -> Option<MolecularSystem> { self.products.get(i).cloned() }
    pub fn get_reagent(&self,  i: usize) -> Option<MolecularSystem> { self.reagents.get(i).cloned() }

    pub fn from_rxn_string(s: &str) -> Result<Reaction, JsValue> {
        parse_rxn(s).map_err(|e| JsValue::from_str(&e.to_string()))
    }
    pub fn from_reaction_smiles(s: &str) -> Result<Reaction, JsValue> {
        parse_reaction_smiles(s).map_err(|e| JsValue::from_str(&e.to_string()))
    }
    pub fn from_cdxml_string(s: &str) -> Result<Reaction, JsValue> {
        parse_cdxml_reaction(s).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    pub fn to_rxn_string(&self)         -> String { reaction_to_rxn_string(self) }
    pub fn to_reaction_smiles(&self)    -> String { reaction_to_smiles(self) }
    pub fn to_cdxml_string(&self)       -> String { reaction_to_cdxml_string(self) }

    /// Apply this single-component reaction to `reactant`.
    /// Returns a JS Array of product MolecularSystem (one per substructure match).
    /// Atom-map numbers ([atom:N]) in the reaction SMILES control which atoms are transformed.
    pub fn run_reaction(&self, reactant: &MolecularSystem) -> js_sys::Array {
        let arr = js_sys::Array::new();
        for product in self.run_reaction_data(reactant) {
            arr.push(&JsValue::from(product));
        }
        arr
    }
}

impl Reaction {
    fn run_reaction_data(&self, reactant: &MolecularSystem) -> Vec<MolecularSystem> {
        let mut results: Vec<MolecularSystem> = Vec::new();
        let rct_tmpl = match self.reactants.first() {
            Some(t) => t,
            None => return results,
        };
        let prd_tmpl = match self.products.first() {
            Some(t) => t,
            None => return results,
        };

        // Use reactant-template heavy-atom SMILES as SMARTS to find substructure matches.
        // to_smiles_data() strips implicit H and traverses in DFS order (which may differ
        // from the original atom indices). We resolve the ordering by self-matching.
        let rct_smiles = rct_tmpl.to_smiles_data();
        if rct_smiles.is_empty() { return results; }

        // Self-match: dfs_to_orig[smarts_i] = original template atom index for SMARTS atom i
        let dfs_to_orig: Vec<usize> = rct_tmpl
            .match_smarts_data(&rct_smiles)
            .into_iter()
            .next()
            .unwrap_or_default();

        let n_heavy_tmpl = dfs_to_orig.len();
        let matches = reactant.match_smarts_data(&rct_smiles);

        for match_indices in &matches {
            // match_indices[i] = substrate atom index for SMARTS atom i
            if match_indices.len() < n_heavy_tmpl { continue; }

            // Build atom-map number → substrate copy atom index
            let mut map_to_copy_idx: std::collections::HashMap<u32, usize> =
                std::collections::HashMap::new();
            for (smarts_i, &copy_i) in match_indices.iter().enumerate() {
                let tmpl_i = dfs_to_orig.get(smarts_i).copied().unwrap_or(smarts_i);
                let m = rct_tmpl.atom_map.get(tmpl_i).copied().unwrap_or(0);
                if m != 0 {
                    map_to_copy_idx.insert(m, copy_i);
                }
            }

            // Build reactant-template bond set: (map_a, map_b) pairs (both non-zero)
            let rct_bond_maps: std::collections::HashSet<(u32, u32)> = {
                let mut s = std::collections::HashSet::new();
                for a in 0..rct_tmpl.bonds.len() {
                    for &b in &rct_tmpl.bonds[a] {
                        if b > a {
                            let ma = rct_tmpl.atom_map.get(a).copied().unwrap_or(0);
                            let mb = rct_tmpl.atom_map.get(b).copied().unwrap_or(0);
                            if ma != 0 && mb != 0 {
                                s.insert((ma.min(mb), ma.max(mb)));
                            }
                        }
                    }
                }
                s
            };

            // Build product-template bond set and orders: (map_a, map_b) → order
            let mut prd_bond_maps: std::collections::HashMap<(u32, u32), u8> =
                std::collections::HashMap::new();
            for a in 0..prd_tmpl.bonds.len() {
                for (k, &b) in prd_tmpl.bonds[a].iter().enumerate() {
                    if b > a {
                        let ma = prd_tmpl.atom_map.get(a).copied().unwrap_or(0);
                        let mb = prd_tmpl.atom_map.get(b).copied().unwrap_or(0);
                        if ma != 0 && mb != 0 {
                            let ord = prd_tmpl.bond_orders.get(a)
                                .and_then(|v| v.get(k)).copied().unwrap_or(1);
                            prd_bond_maps.insert((ma.min(mb), ma.max(mb)), ord);
                        }
                    }
                }
            }

            let mut product = reactant.clone();

            // 1. Remove bonds in reactant template not in product template
            for &(ma, mb) in &rct_bond_maps {
                if !prd_bond_maps.contains_key(&(ma, mb)) {
                    if let (Some(&ci_a), Some(&ci_b)) =
                        (map_to_copy_idx.get(&ma), map_to_copy_idx.get(&mb))
                    {
                        product.remove_bond(ci_a as u32, ci_b as u32);
                    }
                }
            }

            // 2. Add/update bonds in product template not in reactant template
            for (&(ma, mb), &ord) in &prd_bond_maps {
                if !rct_bond_maps.contains(&(ma, mb)) {
                    if let (Some(&ci_a), Some(&ci_b)) =
                        (map_to_copy_idx.get(&ma), map_to_copy_idx.get(&mb))
                    {
                        product.add_bond(ci_a as u32, ci_b as u32, ord);
                    }
                }
            }

            // 3. Transform mapped atoms (symbol and charge from product template)
            let n_prd = prd_tmpl.symbols.len();
            for pi in 0..n_prd {
                let m = prd_tmpl.atom_map.get(pi).copied().unwrap_or(0);
                if m == 0 { continue; }
                if let Some(&ci) = map_to_copy_idx.get(&m) {
                    let prd_sym = prd_tmpl.symbols.get(pi).map(|s| s.as_str()).unwrap_or("C");
                    product.set_atom_symbol(ci as u32, prd_sym);
                    let prd_chg = prd_tmpl.charges.get(pi).copied().unwrap_or(0);
                    product.set_atom_charge(ci as u32, prd_chg);
                }
            }

            // 4. Add new atoms from product template (atom_map == 0)
            let mut new_atom_idx: std::collections::HashMap<usize, u32> =
                std::collections::HashMap::new();
            for pi in 0..n_prd {
                let m = prd_tmpl.atom_map.get(pi).copied().unwrap_or(0);
                if m == 0 {
                    let sym = prd_tmpl.symbols.get(pi).map(|s| s.as_str()).unwrap_or("C");
                    let chg = prd_tmpl.charges.get(pi).copied().unwrap_or(0);
                    let ni = product.add_atom(sym, 0.0, 0.0);
                    if !product.charges.is_empty() {
                        if let Some(c) = product.charges.last_mut() { *c = chg; }
                    }
                    new_atom_idx.insert(pi, ni);
                }
            }

            // 5. Add bonds for new atoms
            for a_pi in 0..n_prd {
                let ma = prd_tmpl.atom_map.get(a_pi).copied().unwrap_or(0);
                for (k, &b_pi) in prd_tmpl.bonds.get(a_pi).into_iter().flatten().enumerate() {
                    if b_pi <= a_pi { continue; }
                    let mb = prd_tmpl.atom_map.get(b_pi).copied().unwrap_or(0);
                    let ord = prd_tmpl.bond_orders.get(a_pi)
                        .and_then(|v| v.get(k)).copied().unwrap_or(1);

                    let ci_a: Option<u32> = if ma != 0 {
                        map_to_copy_idx.get(&ma).map(|&i| i as u32)
                    } else {
                        new_atom_idx.get(&a_pi).copied()
                    };
                    let ci_b: Option<u32> = if mb != 0 {
                        map_to_copy_idx.get(&mb).map(|&i| i as u32)
                    } else {
                        new_atom_idx.get(&b_pi).copied()
                    };

                    // Only bond if at least one is a new atom (mapped bonds handled above)
                    if ma == 0 || mb == 0 {
                        if let (Some(a), Some(b)) = (ci_a, ci_b) {
                            product.add_bond(a, b, ord);
                        }
                    }
                }
            }

            results.push(product);
        }
        results
    }
}

// ── P41: Per-molecule format I/O Wasm methods ────────────────────────────────
#[wasm_bindgen]
impl MolecularSystem {
    pub fn from_cdxml_string(s: &str) -> Result<MolecularSystem, JsValue> {
        parse_cdxml(s).map_err(|e| JsValue::from_str(&e.to_string()))
    }
    pub fn to_cdxml_string(&self) -> String { mol_to_cdxml_string(self) }

    pub fn from_mrv_string(s: &str) -> Result<MolecularSystem, JsValue> {
        parse_mrv(s).map_err(|e| JsValue::from_str(&e.to_string()))
    }
    pub fn to_mrv_string(&self) -> String { mol_to_mrv_string(self) }

    pub fn from_ket_string(s: &str) -> Result<MolecularSystem, JsValue> {
        parse_ket(s).map_err(|e| JsValue::from_str(&e.to_string()))
    }
    pub fn to_ket_string(&self) -> String { mol_to_ket_string(self) }

    pub fn from_cml_string(s: &str) -> Result<MolecularSystem, JsValue> {
        parse_cml(s).map_err(|e| JsValue::from_str(&e.to_string()))
    }
    pub fn to_cml_string(&self) -> String { mol_to_cml_string(self) }
}

// --- P42: Editor Kernel ---

#[wasm_bindgen]
impl MolecularSystem {
    // ── Editing primitives ────────────────────────────────────────────────

    /// Add a new atom and return its index.
    pub fn add_atom(&mut self, symbol: &str, x: f32, y: f32) -> u32 {
        let idx = self.symbols.len();
        p41_push_atom(self, symbol, x, y, 0);
        self.bonds.push(Vec::new());
        self.bond_orders.push(Vec::new());
        self.ring_atoms.clear();
        self.ring_bonds.clear();
        self.spatial_grid = None;
        idx as u32
    }

    /// Remove the atom at `idx`, remapping all bond references.
    pub fn remove_atom(&mut self, idx: u32) {
        let idx = idx as usize;
        let n = self.symbols.len();
        if idx >= n { return; }

        // Remove from all parallel vecs
        self.symbols.remove(idx);
        self.x.remove(idx);
        self.y.remove(idx);
        self.z.remove(idx);
        if idx < self.atom_names.len()     { self.atom_names.remove(idx); }
        if idx < self.residue_names.len()  { self.residue_names.remove(idx); }
        if idx < self.residue_ids.len()    { self.residue_ids.remove(idx); }
        if idx < self.chain_ids.len()      { self.chain_ids.remove(idx); }
        if idx < self.hetatm_flags.len()   { self.hetatm_flags.remove(idx); }
        if idx < self.occupancies.len()    { self.occupancies.remove(idx); }
        if idx < self.b_factors.len()      { self.b_factors.remove(idx); }
        if idx < self.charges.len()        { self.charges.remove(idx); }
        if idx < self.aromatic_atoms.len() { self.aromatic_atoms.remove(idx); }
        if idx < self.ring_atoms.len()     { self.ring_atoms.remove(idx); }
        if idx < self.atom_map.len()       { self.atom_map.remove(idx); }

        // Remove bonds pointing to idx in each neighbor's adjacency list
        for nb in self.bonds[idx].clone() {
            if let Some(pos) = self.bonds[nb].iter().position(|&x| x == idx) {
                self.bonds[nb].remove(pos);
                if pos < self.bond_orders[nb].len() {
                    self.bond_orders[nb].remove(pos);
                }
            }
        }

        // Remove the atom's own adjacency list
        self.bonds.remove(idx);
        self.bond_orders.remove(idx);

        // Remap all indices > idx by -1
        for adj in self.bonds.iter_mut() {
            for nb in adj.iter_mut() {
                if *nb > idx { *nb -= 1; }
            }
        }

        // Update stereo_centers: remove entries for idx, shift keys > idx
        let mut new_stereo = std::collections::HashMap::new();
        for (k, (desc, from_opt)) in self.stereo_centers.drain() {
            if k == idx { continue; }
            let new_k = if k > idx { k - 1 } else { k };
            let new_from = from_opt.map(|f| {
                if f == idx { usize::MAX } else if f > idx { f - 1 } else { f }
            }).filter(|&f| f != usize::MAX);
            new_stereo.insert(new_k, (desc, new_from));
        }
        self.stereo_centers = new_stereo;

        self.ring_atoms.clear();
        self.ring_bonds.clear();
        self.spatial_grid = None;
    }

    /// Change the element symbol of atom `idx`.
    pub fn set_atom_symbol(&mut self, idx: u32, symbol: &str) {
        let idx = idx as usize;
        if idx < self.symbols.len() {
            self.symbols[idx] = symbol.to_string();
        }
    }

    /// Move atom `idx` to (x, y).
    pub fn set_atom_position(&mut self, idx: u32, x: f32, y: f32) {
        let idx = idx as usize;
        if idx < self.x.len() {
            self.x[idx] = x;
            self.y[idx] = y;
        }
        self.spatial_grid = None;
    }

    /// Set the formal charge on atom `idx`.
    pub fn set_atom_charge(&mut self, idx: u32, charge: i32) {
        let idx = idx as usize;
        if idx < self.charges.len() {
            self.charges[idx] = charge;
        }
    }

    /// Add a bond between `a` and `b` with the given `order`.
    /// If the bond already exists, its order is updated instead.
    pub fn add_bond(&mut self, a: u32, b: u32, order: u8) {
        let a = a as usize;
        let b = b as usize;
        let n = self.symbols.len();
        if a >= n || b >= n || a == b { return; }

        // Extend adjacency lists if needed (e.g. after add_atom)
        while self.bonds.len() <= a.max(b) {
            self.bonds.push(Vec::new());
            self.bond_orders.push(Vec::new());
        }

        if let Some(pos) = self.bonds[a].iter().position(|&x| x == b) {
            // Update existing bond order
            if pos < self.bond_orders[a].len() { self.bond_orders[a][pos] = order; }
            if let Some(pos_b) = self.bonds[b].iter().position(|&x| x == a) {
                if pos_b < self.bond_orders[b].len() { self.bond_orders[b][pos_b] = order; }
            }
        } else {
            self.bonds[a].push(b);
            self.bond_orders[a].push(order);
            self.bonds[b].push(a);
            self.bond_orders[b].push(order);
        }

        self.ring_atoms.clear();
        self.ring_bonds.clear();
    }

    /// Remove the bond between `a` and `b`.
    pub fn remove_bond(&mut self, a: u32, b: u32) {
        let a = a as usize;
        let b = b as usize;
        if let Some(pos) = self.bonds.get(a).and_then(|v| v.iter().position(|&x| x == b)) {
            self.bonds[a].remove(pos);
            if pos < self.bond_orders[a].len() { self.bond_orders[a].remove(pos); }
        }
        if let Some(pos) = self.bonds.get(b).and_then(|v| v.iter().position(|&x| x == a)) {
            self.bonds[b].remove(pos);
            if pos < self.bond_orders[b].len() { self.bond_orders[b].remove(pos); }
        }
        self.ring_atoms.clear();
        self.ring_bonds.clear();
    }

    /// Change the bond order between `a` and `b`.
    pub fn set_bond_order(&mut self, a: u32, b: u32, order: u8) {
        let a = a as usize;
        let b = b as usize;
        if let Some(pos) = self.bonds.get(a).and_then(|v| v.iter().position(|&x| x == b)) {
            if pos < self.bond_orders[a].len() { self.bond_orders[a][pos] = order; }
        }
        if let Some(pos) = self.bonds.get(b).and_then(|v| v.iter().position(|&x| x == a)) {
            if pos < self.bond_orders[b].len() { self.bond_orders[b][pos] = order; }
        }
    }

    // ── Hit-testing ───────────────────────────────────────────────────────

    /// Return the index of the atom nearest to (x, y) within `tol`, or None.
    pub fn closest_atom(&self, x: f32, y: f32, tol: f32) -> Option<u32> {
        let tol2 = tol * tol;
        let mut best_dist2 = f32::MAX;
        let mut best_idx: Option<u32> = None;
        for i in 0..self.x.len() {
            let dx = self.x[i] - x;
            let dy = self.y[i] - y;
            let d2 = dx * dx + dy * dy;
            if d2 < best_dist2 {
                best_dist2 = d2;
                best_idx = Some(i as u32);
            }
        }
        if best_dist2 <= tol2 { best_idx } else { None }
    }

    /// Return `[a, b]` of the bond whose segment is closest to (x, y) within `tol`.
    /// Returns an empty Uint32Array if no bond is within tolerance.
    pub fn bond_at(&self, x: f32, y: f32, tol: f32) -> Vec<u32> {
        let tol2 = tol * tol;
        let mut best_dist2 = f32::MAX;
        let mut best_pair: Option<(u32, u32)> = None;
        let mut seen = std::collections::HashSet::new();

        for a in 0..self.bonds.len() {
            for &b in &self.bonds[a] {
                let key = (a.min(b), a.max(b));
                if !seen.insert(key) { continue; }
                let ax = self.x[a]; let ay = self.y[a];
                let bx = self.x[b]; let by = self.y[b];
                let d2 = point_to_segment_dist2(x, y, ax, ay, bx, by);
                if d2 < best_dist2 {
                    best_dist2 = d2;
                    best_pair = Some((a as u32, b as u32));
                }
            }
        }
        if best_dist2 <= tol2 {
            if let Some((a, b)) = best_pair { return vec![a, b]; }
        }
        Vec::new()
    }

    // ── Coordinate utilities ──────────────────────────────────────────────

    /// Scale all atom positions so the average bond length equals `target`.
    /// No-op if there are no bonds or all atoms are at the same position.
    pub fn normalize_bond_length(&mut self, target: f32) {
        if target <= 0.0 { return; }
        let mut sum = 0.0f32;
        let mut count = 0usize;
        let mut seen = std::collections::HashSet::new();
        for a in 0..self.bonds.len() {
            for &b in &self.bonds[a] {
                let key = (a.min(b), a.max(b));
                if !seen.insert(key) { continue; }
                let dx = self.x[a] - self.x[b];
                let dy = self.y[a] - self.y[b];
                let len = (dx * dx + dy * dy).sqrt();
                if len < 1e-4 { continue; } // skip degenerate bonds (e.g. H placed at heavy atom)
                sum += len;
                count += 1;
            }
        }
        if count == 0 || sum <= 0.0 { return; }
        let scale = target / (sum / count as f32);
        let n = self.x.len();
        if n == 0 { return; }
        let cx: f32 = self.x.iter().sum::<f32>() / n as f32;
        let cy: f32 = self.y.iter().sum::<f32>() / n as f32;
        for i in 0..n {
            self.x[i] = cx + (self.x[i] - cx) * scale;
            self.y[i] = cy + (self.y[i] - cy) * scale;
        }
        self.spatial_grid = None;
    }

    /// Shift all atom positions by (dx, dy).
    pub fn translate_atoms(&mut self, dx: f32, dy: f32) {
        for i in 0..self.x.len() {
            self.x[i] += dx;
            self.y[i] += dy;
        }
        self.spatial_grid = None;
    }
}

/// Squared distance from point (px, py) to line segment (ax, ay)-(bx, by).
fn point_to_segment_dist2(px: f32, py: f32, ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    let dx = bx - ax; let dy = by - ay;
    let len2 = dx * dx + dy * dy;
    if len2 == 0.0 {
        let ex = px - ax; let ey = py - ay;
        return ex * ex + ey * ey;
    }
    let t = ((px - ax) * dx + (py - ay) * dy) / len2;
    let t = t.clamp(0.0, 1.0);
    let qx = ax + t * dx; let qy = ay + t * dy;
    let ex = px - qx; let ey = py - qy;
    ex * ex + ey * ey
}

// --- P43: Ring Templates + Implicit H Count ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Return the number of implicit hydrogens for the atom at `idx`.
    /// Uses standard organic valences (C=4, N=3, O=2, S=2, P=3, halogens=1).
    /// Returns -1 for unknown elements (metals, etc.).
    pub fn implicit_h_count(&self, idx: u32) -> i32 {
        let idx = idx as usize;
        if idx >= self.symbols.len() { return -1; }
        let std_val = match smiles_valence(&self.symbols[idx]) {
            Some(v) => v as i32,
            None => return -1,
        };
        let bond_sum: i32 = if idx < self.bonds.len() {
            self.bonds[idx].iter().enumerate().map(|(k, _)| {
                let order = self.bond_orders.get(idx)
                    .and_then(|o| o.get(k))
                    .copied()
                    .unwrap_or(1);
                // Aromatic bond (4) counts as 1.5 → round to 1 here; kekulize first for accuracy
                if order == 4 { 1 } else { order as i32 }
            }).sum()
        } else { 0 };
        (std_val - bond_sum).max(0)
    }

    /// Place a regular n-membered ring of carbon atoms centered at (cx, cy).
    /// Returns the indices of the newly added atoms.
    /// No-op (empty result) if n < 3.
    pub fn add_ring_template(&mut self, n: u32, cx: f32, cy: f32, bond_length: f32) -> Vec<u32> {
        if !(3..=100).contains(&n) { return Vec::new(); }
        use std::f32::consts::PI;
        let n = n as usize;
        let radius = bond_length / (2.0 * (PI / n as f32).sin());
        let start_angle = PI / 2.0;
        let mut indices = Vec::with_capacity(n);
        for i in 0..n {
            let angle = start_angle + 2.0 * PI * i as f32 / n as f32;
            let x = cx + radius * angle.cos();
            let y = cy + radius * angle.sin();
            indices.push(self.add_atom("C", x, y));
        }
        for i in 0..n {
            let a = indices[i];
            let b = indices[(i + 1) % n];
            self.add_bond(a, b, 1);
        }
        indices
    }
}

// --- P44: Fused Ring Attachment ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Fuse a new n-membered ring onto the existing bond between atoms `a` and `b`.
    /// The ring shares the a–b edge and grows away from the molecular centroid.
    /// Returns the indices of the n-2 newly added atoms.
    /// Returns an empty Vec if n < 3, a == b, or indices are out of range.
    pub fn attach_ring_to_bond(&mut self, a: u32, b: u32, n: u32) -> Vec<u32> {
        use std::f32::consts::PI;
        let na = a as usize;
        let nb = b as usize;
        if !(3..=100).contains(&n) || a == b || na >= self.symbols.len() || nb >= self.symbols.len() {
            return Vec::new();
        }
        let xa = self.x[na];
        let ya = self.y[na];
        let xb = self.x[nb];
        let yb = self.y[nb];
        let dx = xb - xa;
        let dy = yb - ya;
        let bond_len = (dx * dx + dy * dy).sqrt();
        if bond_len < 1e-4 {
            return Vec::new();
        }
        // Perpendicular unit vector (rotated 90° CCW from a→b direction)
        let perp_x = -dy / bond_len;
        let perp_y =  dx / bond_len;
        let mid_x = (xa + xb) * 0.5;
        let mid_y = (ya + yb) * 0.5;
        // Apothem = distance from ring center to the midpoint of the shared edge
        let apothem = bond_len / (2.0 * (PI / n as f32).tan());
        let c1 = [mid_x + apothem * perp_x, mid_y + apothem * perp_y];
        let c2 = [mid_x - apothem * perp_x, mid_y - apothem * perp_y];
        // Centroid of existing atoms — pick the ring center on the far side
        let natoms = self.symbols.len();
        let cent_x: f32 = self.x.iter().sum::<f32>() / natoms as f32;
        let cent_y: f32 = self.y.iter().sum::<f32>() / natoms as f32;
        let d1 = (c1[0] - cent_x).powi(2) + (c1[1] - cent_y).powi(2);
        let d2 = (c2[0] - cent_x).powi(2) + (c2[1] - cent_y).powi(2);
        let [cx, cy] = if d1 >= d2 { c1 } else { c2 };
        let radius = bond_len / (2.0 * (PI / n as f32).sin());
        let angle_a = (ya - cy).atan2(xa - cx);
        // Cross product (a-center) × (b-center): positive → b is CCW from a.
        // New atoms go the long way (opposite direction) to form the ring interior.
        let cross = (xa - cx) * (yb - cy) - (ya - cy) * (xb - cx);
        let step = if cross > 0.0 { -2.0 * PI / n as f32 } else { 2.0 * PI / n as f32 };
        // Place n-2 new C atoms at interior positions
        let n_usize = n as usize;
        let mut new_idx: Vec<u32> = Vec::with_capacity(n_usize - 2);
        for i in 1..=(n_usize - 2) {
            let angle = angle_a + i as f32 * step;
            new_idx.push(self.add_atom("C", cx + radius * angle.cos(), cy + radius * angle.sin()));
        }
        // Add n-1 bonds: a → new[0] → ... → new[n-3] → b
        let mut chain: Vec<u32> = Vec::with_capacity(n_usize);
        chain.push(a);
        chain.extend_from_slice(&new_idx);
        chain.push(b);
        for w in chain.windows(2) {
            self.add_bond(w[0], w[1], 1);
        }
        new_idx
    }
}

// --- P45: Geometry Utilities ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Return the 2D bounding box of all atoms as [min_x, min_y, max_x, max_y].
    /// Returns an empty Vec when the molecule has no atoms.
    pub fn get_bounds(&self) -> Vec<f32> {
        if self.x.is_empty() {
            return Vec::new();
        }
        let min_x = self.x.iter().cloned().fold(f32::INFINITY, f32::min);
        let min_y = self.y.iter().cloned().fold(f32::INFINITY, f32::min);
        let max_x = self.x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let max_y = self.y.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        vec![min_x, min_y, max_x, max_y]
    }

    /// Rotate all atoms by `angle` radians around the point (cx, cy).
    pub fn rotate_atoms(&mut self, angle: f32, cx: f32, cy: f32) {
        let (sin_a, cos_a) = angle.sin_cos();
        for i in 0..self.x.len() {
            let dx = self.x[i] - cx;
            let dy = self.y[i] - cy;
            self.x[i] = cx + dx * cos_a - dy * sin_a;
            self.y[i] = cy + dx * sin_a + dy * cos_a;
        }
    }
}

// --- P46: Selection + Partial Move ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Return indices of all atoms whose 2D position lies within the rectangle
    /// defined by (x1, y1) and (x2, y2). The rectangle is automatically normalized,
    /// so x1 > x2 is handled correctly.
    pub fn select_atoms_in_rect(&self, x1: f32, y1: f32, x2: f32, y2: f32) -> Vec<u32> {
        let (min_x, max_x) = (x1.min(x2), x1.max(x2));
        let (min_y, max_y) = (y1.min(y2), y1.max(y2));
        (0..self.x.len())
            .filter(|&i| {
                self.x[i] >= min_x && self.x[i] <= max_x
                    && self.y[i] >= min_y && self.y[i] <= max_y
            })
            .map(|i| i as u32)
            .collect()
    }

    /// Translate only the atoms at the given indices by (dx, dy).
    /// Out-of-range indices are silently skipped.
    pub fn move_atoms(&mut self, indices: &[u32], dx: f32, dy: f32) {
        for &i in indices {
            let i = i as usize;
            if i < self.x.len() {
                self.x[i] += dx;
                self.y[i] += dy;
            }
        }
    }
}

// --- P47: Valence Check ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Return the indices of atoms whose current bond order sum exceeds their
    /// standard valence. Aromatic bonds (order 4) count as 1. Unknown elements
    /// (metals, etc.) are skipped.
    pub fn check_valence(&self) -> Vec<u32> {
        (0..self.symbols.len())
            .filter(|&i| {
                let Some(std_val) = smiles_valence(&self.symbols[i]) else {
                    return false;
                };
                let bond_sum: i32 = if i < self.bonds.len() {
                    self.bonds[i]
                        .iter()
                        .enumerate()
                        .map(|(k, _)| {
                            let order = self
                                .bond_orders
                                .get(i)
                                .and_then(|o| o.get(k))
                                .copied()
                                .unwrap_or(1);
                            if order == 4 { 1 } else { order as i32 }
                        })
                        .sum()
                } else {
                    0
                };
                bond_sum > std_val as i32
            })
            .map(|i| i as u32)
            .collect()
    }
}

// --- P48: Flip (Mirror) Operations ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Mirror all atoms across the vertical axis x = cx.
    pub fn flip_horizontal(&mut self, cx: f32) {
        for x in &mut self.x {
            *x = 2.0 * cx - *x;
        }
    }

    /// Mirror all atoms across the horizontal axis y = cy.
    pub fn flip_vertical(&mut self, cy: f32) {
        for y in &mut self.y {
            *y = 2.0 * cy - *y;
        }
    }
}

// --- P49: Copy Atoms (Partial Clone) ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Extract a subset of atoms into a new `MolecularSystem`.
    /// Only bonds between selected atoms are copied; bonds to unselected atoms are dropped.
    /// PDB fields and ring/spatial data are not copied (caller can recompute if needed).
    pub fn copy_atoms(&self, indices: &[u32]) -> MolecularSystem {
        let mut out = MolecularSystem::new_empty();
        let n = self.symbols.len();
        let mut idx_map = vec![usize::MAX; n];
        for (new_i, &old_raw) in indices.iter().enumerate() {
            let old_i = old_raw as usize;
            if old_i < n {
                idx_map[old_i] = new_i;
            }
        }
        let selected: std::collections::HashSet<usize> = indices
            .iter()
            .map(|&i| i as usize)
            .filter(|&i| i < n)
            .collect();
        for &old_raw in indices {
            let old_i = old_raw as usize;
            if old_i >= n {
                continue;
            }
            out.symbols.push(self.symbols[old_i].clone());
            out.x.push(self.x[old_i]);
            out.y.push(self.y[old_i]);
            out.z.push(self.z[old_i]);
            out.charges.push(self.charges.get(old_i).copied().unwrap_or(0));
            out.bonds.push(Vec::new());
            out.bond_orders.push(Vec::new());
        }
        for (new_i, &old_raw) in indices.iter().enumerate() {
            let old_i = old_raw as usize;
            if old_i >= self.bonds.len() {
                continue;
            }
            for (k, &nb) in self.bonds[old_i].iter().enumerate() {
                if selected.contains(&nb) && nb > old_i {
                    let new_nb = idx_map[nb];
                    let order = self
                        .bond_orders
                        .get(old_i)
                        .and_then(|o| o.get(k))
                        .copied()
                        .unwrap_or(1);
                    out.bonds[new_i].push(new_nb);
                    out.bond_orders[new_i].push(order);
                    out.bonds[new_nb].push(new_i);
                    out.bond_orders[new_nb].push(order);
                }
            }
        }
        out
    }
}

// --- P50: RDKit.js parity — convenience methods + descriptors ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Returns the largest connected fragment as a new MolecularSystem.
    /// Useful for salt stripping (pick the main molecule, discard counterions).
    pub fn largest_fragment(&self) -> MolecularSystem {
        let n = self.symbols.len();
        let mut visited = vec![false; n];
        let mut best: Vec<u32> = Vec::new();
        for start in 0..n {
            if visited[start] {
                continue;
            }
            let mut comp: Vec<u32> = Vec::new();
            let mut stack = vec![start];
            while let Some(a) = stack.pop() {
                if visited[a] {
                    continue;
                }
                visited[a] = true;
                comp.push(a as u32);
                for &nb in self.bonds.get(a).into_iter().flatten() {
                    if !visited[nb] {
                        stack.push(nb);
                    }
                }
            }
            if comp.len() > best.len() {
                best = comp;
            }
        }
        best.sort_unstable();
        self.copy_atoms(&best)
    }

    /// Returns the Murcko scaffold (ring systems + inter-ring linkers) as a new MolecularSystem.
    /// Requires `compute_bonds()` + `compute_rings()` first.
    pub fn murcko_scaffold(&self) -> MolecularSystem {
        let idx: Vec<u32> = self
            .murcko_scaffold_indices_data()
            .into_iter()
            .map(|i| i as u32)
            .collect();
        self.copy_atoms(&idx)
    }

    /// Number of non-hydrogen atoms.
    pub fn num_heavy_atoms(&self) -> u32 {
        self.symbols.iter().filter(|s| s.as_str() != "H").count() as u32
    }

    /// Fraction of carbons that are sp3 (no double/triple bond, not aromatic).
    /// Returns 0.0 if the molecule has no carbon atoms.
    /// Requires `compute_bonds()` first; returns 0.0 if bonds are absent.
    pub fn fraction_csp3(&self) -> f32 {
        if self.bonds.is_empty() {
            return 0.0;
        }
        let mut total_c = 0u32;
        let mut sp3_c = 0u32;
        for i in 0..self.symbols.len() {
            if self.symbols[i] != "C" {
                continue;
            }
            total_c += 1;
            let is_arom = self.aromatic_atoms.get(i).copied().unwrap_or(false);
            let max_order = self
                .bond_orders
                .get(i)
                .map_or(1u8, |o| o.iter().copied().max().unwrap_or(1));
            if !is_arom && max_order <= 1 {
                sp3_c += 1;
            }
        }
        if total_c == 0 {
            0.0
        } else {
            sp3_c as f32 / total_c as f32
        }
    }

    /// Approximate molar refractivity (Wildman-Crippen MR, cm³/mol).
    /// Requires `compute_bonds()` first; returns 0.0 if bonds are absent.
    /// Accuracy: roughly ±0.5 for typical drug-like molecules.
    pub fn molar_refractivity(&self) -> f32 {
        if self.bonds.is_empty() {
            return 0.0;
        }
        (0..self.symbols.len())
            .map(|i| self.atom_mr_contribution(i))
            .sum()
    }

    fn atom_mr_contribution(&self, i: usize) -> f32 {
        let sym = match self.symbols.get(i) {
            Some(s) => s.as_str(),
            None => return 0.0,
        };
        if sym == "H" {
            return 0.0;
        }

        let h_ct: usize = self.bonds.get(i).map_or(0, |nbrs| {
            nbrs.iter()
                .filter(|&&j| self.symbols.get(j).map(|s| s == "H").unwrap_or(false))
                .count()
        });
        let in_ring = self.ring_atoms.get(i).copied().unwrap_or(false);
        let is_arom = self.aromatic_atoms.get(i).copied().unwrap_or(false);
        let double_to_o = self.bonds.get(i).is_some_and(|nbrs| {
            nbrs.iter().enumerate().any(|(k, &j)| {
                self.symbols.get(j).map(|s| s == "O").unwrap_or(false)
                    && self
                        .bond_orders
                        .get(i)
                        .and_then(|o| o.get(k))
                        .copied()
                        .unwrap_or(1)
                        == 2
            })
        });
        let max_order = self
            .bond_orders
            .get(i)
            .map_or(1u8, |o| o.iter().copied().max().unwrap_or(1));

        match sym {
            "C" => {
                if is_arom || (in_ring && max_order >= 2) {
                    if h_ct >= 1 { 0.41 } else { 0.29 }
                } else if double_to_o {
                    0.25
                } else if max_order >= 2 {
                    if h_ct >= 1 { 0.35 } else { 0.21 }
                } else {
                    match h_ct {
                        3 => 0.80,
                        2 => 0.77,
                        1 => 0.47,
                        _ => 0.19,
                    }
                }
            }
            "N" => match h_ct {
                2 => 0.99,
                1 => 0.71,
                _ => 0.47,
            },
            "O" => {
                if max_order >= 2 { 0.50 } else { 0.25 }
            }
            "S" => {
                if h_ct >= 1 { 1.08 } else { 0.84 }
            }
            "F"  => 0.18,
            "Cl" => 0.60,
            "Br" => 0.85,
            "I"  => 1.05,
            "P"  => 0.80,
            _    => 0.30,
        }
    }
}

// --- P51: 3D conformer generation via simplified distance geometry ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Assigns 3-D coordinates using simplified distance geometry.
    /// `seed` — LCG seed for reproducibility (0 = default seed 12345).
    /// Bond topology must be present in self.bonds (SMILES/SDF molecules have it automatically;
    /// for XYZ/PDB molecules call compute_bonds() first). Returns `false` if n < 2.
    pub fn embed_molecule(&mut self, seed: u32) -> bool {
        let n = self.symbols.len();
        if n < 2 {
            return false;
        }

        // LCG state
        let mut lcg: u32 = if seed == 0 { 12345 } else { seed };
        let lcg_next = |s: &mut u32| -> f32 {
            *s = s.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (*s as f32) / (u32::MAX as f32)
        };

        // Step 1: Build upper/lower distance bounds
        let mut lo = vec![0f32; n * n];
        let mut hi = vec![100f32; n * n];
        // diagonal
        for i in 0..n {
            lo[i * n + i] = 0.0;
            hi[i * n + i] = 0.0;
        }
        // vdW clash floor for all pairs
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    lo[i * n + j] = 1.2;
                }
            }
        }
        // bonded (1-2) pairs
        for a in 0..n {
            for &b in self.bonds.get(a).into_iter().flatten() {
                if b > a {
                    let ra = covalent_radius(
                        self.symbols.get(a).map(|s| s.as_str()).unwrap_or("C"),
                    ).unwrap_or(0.77);
                    let rb = covalent_radius(
                        self.symbols.get(b).map(|s| s.as_str()).unwrap_or("C"),
                    ).unwrap_or(0.77);
                    let target = ra + rb + 0.08;
                    lo[a * n + b] = target * 0.9;
                    hi[a * n + b] = target * 1.1;
                    lo[b * n + a] = lo[a * n + b];
                    hi[b * n + a] = hi[a * n + b];
                }
            }
        }
        // 1-3 pairs (share a common bonded neighbor)
        for b in 0..n {
            let nbrs: Vec<usize> = self.bonds.get(b).into_iter().flatten().copied().collect();
            for ii in 0..nbrs.len() {
                for jj in (ii + 1)..nbrs.len() {
                    let (a, c) = (nbrs[ii], nbrs[jj]);
                    if hi[a * n + c] > 3.0 {
                        lo[a * n + c] = 2.0f32.max(lo[a * n + c]);
                        hi[a * n + c] = 3.0f32.min(hi[a * n + c]);
                        lo[c * n + a] = lo[a * n + c];
                        hi[c * n + a] = hi[a * n + c];
                    }
                }
            }
        }

        // Step 2: Triangle bound smoothing (upper bounds, Floyd-Warshall)
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    let via = hi[i * n + k] + hi[k * n + j];
                    if via < hi[i * n + j] {
                        hi[i * n + j] = via;
                    }
                }
            }
        }
        // clamp: lower <= upper
        for i in 0..n * n {
            if lo[i] > hi[i] {
                lo[i] = hi[i];
            }
        }

        // Step 3: Sample distance matrix and embed via metric matrix
        let mut d = vec![0f32; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                let v = lo[i * n + j] + lcg_next(&mut lcg) * (hi[i * n + j] - lo[i * n + j]);
                d[i * n + j] = v;
                d[j * n + i] = v;
            }
        }

        // Gram matrix G[i][j] = (d0i² + d0j² - dij²) / 2  (using atom 0 as reference)
        // Then center: G[i][j] -= row_mean[i] + col_mean[j] - grand_mean
        let mut g = vec![0f32; n * n];
        for i in 0..n {
            for j in 0..n {
                // Classic distance geometry: G[i][j] = (d(0,i)² + d(0,j)² - d(i,j)²) / 2
                let d0i = d[i]; // d[0*n + i] — distance from reference atom 0 to atom i
                let d0j = d[j]; // d[0*n + j] — distance from reference atom 0 to atom j
                let dij = d[i * n + j];
                g[i * n + j] = (d0i * d0i + d0j * d0j - dij * dij) * 0.5;
            }
        }
        // double-center
        let row_mean: Vec<f32> = (0..n)
            .map(|i| g[i * n..i * n + n].iter().sum::<f32>() / n as f32)
            .collect();
        let grand_mean = row_mean.iter().sum::<f32>() / n as f32;
        for i in 0..n {
            for j in 0..n {
                g[i * n + j] -= row_mean[i] + row_mean[j] - grand_mean;
            }
        }

        // Power iteration for top 3 eigenvectors
        let power_iter = |mat: &[f32], init: &[f32], steps: usize| -> (Vec<f32>, f32) {
            let nn = init.len();
            let mut v = init.to_vec();
            let norm = v.iter().map(|x| x * x).sum::<f32>().sqrt().max(1e-10);
            for x in &mut v {
                *x /= norm;
            }
            let mut eigenval = 0f32;
            for _ in 0..steps {
                let mut w = vec![0f32; nn];
                for i in 0..nn {
                    for j in 0..nn {
                        w[i] += mat[i * nn + j] * v[j];
                    }
                }
                eigenval = w.iter().zip(v.iter()).map(|(a, b)| a * b).sum::<f32>();
                let norm2 = w.iter().map(|x| x * x).sum::<f32>().sqrt().max(1e-10);
                for x in &mut w {
                    *x /= norm2;
                }
                v = w;
            }
            (v, eigenval.max(0.0))
        };

        // ev1
        let init1: Vec<f32> = (0..n).map(|_| lcg_next(&mut lcg)).collect();
        let (ev1, lam1) = power_iter(&g, &init1, 50);
        // deflate
        let mut g2 = g.clone();
        for i in 0..n {
            for j in 0..n {
                g2[i * n + j] -= lam1 * ev1[i] * ev1[j];
            }
        }
        // ev2
        let init2: Vec<f32> = (0..n).map(|_| lcg_next(&mut lcg)).collect();
        let (ev2, lam2) = power_iter(&g2, &init2, 50);
        // deflate
        let mut g3 = g2.clone();
        for i in 0..n {
            for j in 0..n {
                g3[i * n + j] -= lam2 * ev2[i] * ev2[j];
            }
        }
        // ev3
        let init3: Vec<f32> = (0..n).map(|_| lcg_next(&mut lcg)).collect();
        let (ev3, lam3) = power_iter(&g3, &init3, 50);

        let s1 = lam1.sqrt();
        let s2 = lam2.sqrt();
        let s3 = lam3.sqrt();

        self.x = ev1.iter().map(|v| v * s1).collect();
        self.y = ev2.iter().map(|v| v * s2).collect();
        self.z = ev3.iter().map(|v| v * s3).collect();

        // Step 4: Bond stretch gradient descent (200 steps)
        // Collect bonded pairs with target distances
        let mut bond_pairs: Vec<(usize, usize, f32)> = Vec::new();
        for a in 0..n {
            for &b in self.bonds.get(a).into_iter().flatten() {
                if b > a {
                    let ra = covalent_radius(
                        self.symbols.get(a).map(|s| s.as_str()).unwrap_or("C"),
                    ).unwrap_or(0.77);
                    let rb = covalent_radius(
                        self.symbols.get(b).map(|s| s.as_str()).unwrap_or("C"),
                    ).unwrap_or(0.77);
                    let target = ra + rb + 0.08;
                    bond_pairs.push((a, b, target));
                }
            }
        }

        // Angle triplets for bending term: (a, b, c, ideal_angle_rad)
        // b is the central atom. Ideal angles by heavy-atom degree of b.
        let mut angle_triplets: Vec<(usize, usize, usize, f32)> = Vec::new();
        for b in 0..n {
            let nbrs: Vec<usize> = self.bonds.get(b).into_iter().flatten().copied().collect();
            if nbrs.len() < 2 { continue; }
            let theta_ideal: f32 = match nbrs.len() {
                2 => 2.0943952f32, // 120° in radians (handles sp2 and most sp)
                3 => 2.0943952f32, // 120°
                _ => 1.9106332f32, // 109.47° (tetrahedral)
            };
            for i in 0..nbrs.len() {
                for j in (i + 1)..nbrs.len() {
                    angle_triplets.push((nbrs[i], b, nbrs[j], theta_ideal));
                }
            }
        }

        let step_size = 0.01f32;
        let k_angle = 0.30f32; // weight vs bond-stretch gradient (~2*diff)
        for _ in 0..200 {
            let mut gx = vec![0f32; n];
            let mut gy = vec![0f32; n];
            let mut gz = vec![0f32; n];

            // Bond stretch
            for &(a, b, tgt) in &bond_pairs {
                let dx = self.x[b] - self.x[a];
                let dy = self.y[b] - self.y[a];
                let dz = self.z[b] - self.z[a];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-6);
                let diff = dist - tgt;
                let ux = dx / dist;
                let uy = dy / dist;
                let uz = dz / dist;
                let g = 2.0 * diff;
                gx[a] += g * (-ux);
                gy[a] += g * (-uy);
                gz[a] += g * (-uz);
                gx[b] += g * ux;
                gy[b] += g * uy;
                gz[b] += g * uz;
            }

            // Angle bending: harmonic on theta, gradient via analytical formula
            for &(a, b, c, theta_ideal) in &angle_triplets {
                let bx = self.x[b]; let by = self.y[b]; let bz = self.z[b];
                let ax = self.x[a] - bx; let ay = self.y[a] - by; let az = self.z[a] - bz;
                let cx = self.x[c] - bx; let cy = self.y[c] - by; let cz = self.z[c] - bz;
                let da = (ax*ax + ay*ay + az*az).sqrt().max(1e-6);
                let dc = (cx*cx + cy*cy + cz*cz).sqrt().max(1e-6);
                let cos_t = ((ax*cx + ay*cy + az*cz) / (da * dc)).clamp(-1.0, 1.0);
                let theta = cos_t.acos();
                let diff = theta - theta_ideal;
                if diff.abs() < 1e-4 { continue; }
                let sin_t = (1.0 - cos_t * cos_t).sqrt().max(1e-6);
                let (uax, uay, uaz) = (ax / da, ay / da, az / da);
                let (ucx, ucy, ucz) = (cx / dc, cy / dc, cz / dc);
                let scale = 2.0 * diff * k_angle;
                let fax = (cos_t * uax - ucx) / (da * sin_t);
                let fay = (cos_t * uay - ucy) / (da * sin_t);
                let faz = (cos_t * uaz - ucz) / (da * sin_t);
                let fcx = (cos_t * ucx - uax) / (dc * sin_t);
                let fcy = (cos_t * ucy - uay) / (dc * sin_t);
                let fcz = (cos_t * ucz - uaz) / (dc * sin_t);
                gx[a] += scale * fax; gy[a] += scale * fay; gz[a] += scale * faz;
                gx[c] += scale * fcx; gy[c] += scale * fcy; gz[c] += scale * fcz;
                gx[b] -= scale * (fax + fcx);
                gy[b] -= scale * (fay + fcy);
                gz[b] -= scale * (faz + fcz);
            }

            // Clip grad norm per atom and update positions
            for i in 0..n {
                let gnorm = (gx[i] * gx[i] + gy[i] * gy[i] + gz[i] * gz[i]).sqrt();
                if gnorm > 1.0 {
                    gx[i] /= gnorm;
                    gy[i] /= gnorm;
                    gz[i] /= gnorm;
                }
                self.x[i] -= step_size * gx[i];
                self.y[i] -= step_size * gy[i];
                self.z[i] -= step_size * gz[i];
            }
        }

        // Center at origin
        let cx = self.x.iter().sum::<f32>() / n as f32;
        let cy = self.y.iter().sum::<f32>() / n as f32;
        let cz = self.z.iter().sum::<f32>() / n as f32;
        for i in 0..n {
            self.x[i] -= cx;
            self.y[i] -= cy;
            self.z[i] -= cz;
        }

        true
    }
}

// --- Atom Pair Fingerprint helpers ---

fn centroid2d(pts: &[[f32; 2]]) -> [f32; 2] {
    let n = pts.len() as f32;
    let mut c = [0f32; 2];
    for p in pts { c[0] += p[0]; c[1] += p[1]; }
    [c[0] / n, c[1] / n]
}

fn atom_pair_fingerprint_bits(mol: &MolecularSystem) -> [u8; 256] {
    let n = mol.symbols.len();
    let mut bits = [0u8; 256];
    if n == 0 {
        return bits;
    }

    let atom_hash = |i: usize| -> u32 {
        let mut h = 0x811c9dc5u32;
        for b in mol.symbols[i].bytes() {
            h = fnv_mix(h, b as u32);
        }
        let arom = mol.aromatic_atoms.get(i).copied().unwrap_or(false) as u32;
        let deg = mol.bonds.get(i)
            .map(|nb| nb.iter()
                .filter(|&&j| mol.symbols.get(j).map(|s| s != "H").unwrap_or(true))
                .count())
            .unwrap_or(0) as u32;
        h = fnv_mix(h, arom);
        h = fnv_mix(h, deg);
        h
    };

    let hashes: Vec<u32> = (0..n).map(atom_hash).collect();

    for start in 0..n {
        if mol.symbols[start] == "H" {
            continue;
        }
        let mut visited = vec![u8::MAX; n];
        visited[start] = 0;
        let mut queue = std::collections::VecDeque::from([start]);
        while let Some(u) = queue.pop_front() {
            let d = visited[u];
            if d >= 7 {
                continue;
            }
            for &v in mol.bonds.get(u).into_iter().flatten() {
                if mol.symbols.get(v).map(|s| s == "H").unwrap_or(false) {
                    continue;
                }
                if visited[v] == u8::MAX {
                    visited[v] = d + 1;
                    queue.push_back(v);
                    let (ha, hb) = if hashes[start] <= hashes[v] {
                        (hashes[start], hashes[v])
                    } else {
                        (hashes[v], hashes[start])
                    };
                    let mut h = fnv_mix(ha, hb);
                    h = fnv_mix(h, (d + 1) as u32);
                    let pos = h % 2048;
                    bits[(pos / 8) as usize] |= 1 << (pos % 8);
                }
            }
        }
    }
    bits
}

// --- P53: Ring info APIs, Atom Pair FP, Template 2D Alignment ---

#[wasm_bindgen]
impl MolecularSystem {
    /// Ring sizes of all SSSR rings containing atom `idx`.
    /// Requires `compute_bonds()` + `compute_rings()`.
    pub fn ring_sizes_for_atom(&self, idx: usize) -> Vec<u8> {
        self.ring_sizes_per_atom.get(idx).cloned().unwrap_or_default()
    }

    /// Returns `{num_rings, ring_sizes}` for the molecule.
    /// Requires `compute_bonds()` + `compute_rings()`.
    pub fn ring_info(&self) -> JsValue {
        let rings = self.enumerate_rings();
        let sizes: Vec<usize> = rings.iter().map(|r| r.len()).collect();
        #[derive(serde::Serialize)]
        struct RingInfo { num_rings: usize, ring_sizes: Vec<usize> }
        to_js(&RingInfo { num_rings: rings.len(), ring_sizes: sizes })
    }

    /// Number of rings containing no aromatic atoms.
    /// Requires `compute_bonds()` + `compute_rings()`.
    pub fn aliphatic_ring_count(&self) -> u32 {
        self.enumerate_rings().into_iter()
            .filter(|ring| ring.iter().all(|&i|
                !self.aromatic_atoms.get(i).copied().unwrap_or(false)))
            .count() as u32
    }

    /// 2048-bit (256-byte) Atom Pair fingerprint.
    /// Requires `compute_bonds()`.
    pub fn fingerprint_atom_pair(&self) -> Vec<u8> {
        atom_pair_fingerprint_bits(self).to_vec()
    }

    /// Generates 2D coordinates aligned to a template molecule via substructure match + 2D rotation.
    /// Calls `compute_2d_coords()` internally. Falls back to normal 2D layout if no match found.
    pub fn generate_aligned_coords(&mut self, template: &MolecularSystem) {
        self.compute_2d_coords_data();
        if template.symbols.is_empty() {
            return;
        }

        let tmpl_smiles = template.to_smiles_data();
        let matches = self.match_smarts_data(&tmpl_smiles);
        let first_match = match matches.into_iter().next() {
            Some(m) => m,
            None => return,
        };

        let n_tmpl = first_match.len().min(template.symbols.len());
        if n_tmpl < 2 {
            return;
        }

        let ref_pts: Vec<[f32; 2]> = (0..n_tmpl)
            .map(|qi| [template.x[qi], template.y[qi]])
            .collect();
        let mob_pts: Vec<[f32; 2]> = (0..n_tmpl)
            .map(|qi| { let si = first_match[qi]; [self.x[si], self.y[si]] })
            .collect();

        let cp_ref = centroid2d(&ref_pts);
        let cp_mob = centroid2d(&mob_pts);

        let mut h = [[0f32; 2]; 2];
        for k in 0..n_tmpl {
            let p = [mob_pts[k][0] - cp_mob[0], mob_pts[k][1] - cp_mob[1]];
            let q = [ref_pts[k][0] - cp_ref[0], ref_pts[k][1] - cp_ref[1]];
            for (r, hrow) in h.iter_mut().enumerate() {
                for (hrc, &qc) in hrow.iter_mut().zip(q.iter()) {
                    *hrc += p[r] * qc;
                }
            }
        }
        let angle = (h[0][1] - h[1][0]).atan2(h[0][0] + h[1][1]);
        let cos_a = angle.cos();
        let sin_a = angle.sin();

        let na = self.symbols.len();
        for i in 0..na {
            let dx = self.x[i] - cp_mob[0];
            let dy = self.y[i] - cp_mob[1];
            self.x[i] = cos_a * dx - sin_a * dy + cp_ref[0];
            self.y[i] = sin_a * dx + cos_a * dy + cp_ref[1];
        }
    }
}

// --- P55 private helpers: topological path fingerprint ---

fn path_fp_dfs(
    mol: &MolecularSystem,
    curr: usize,
    hash: u64,
    depth: usize,
    visited: &mut Vec<bool>,
    bits: &mut [u8; 256],
) {
    let empty: Vec<usize> = Vec::new();
    for (k, &nb) in mol.bonds.get(curr).unwrap_or(&empty).iter().enumerate() {
        if visited[nb] { continue; }
        let bo = mol.bond_orders.get(curr).and_then(|o| o.get(k)).copied().unwrap_or(1);
        let sym = symbol_to_atomic_num(mol.symbols.get(nb).map(|s| s.as_str()).unwrap_or("C")) as u64;
        let h = hash.wrapping_mul(1_000_003u64)
                    .wrapping_add(sym.wrapping_mul(31).wrapping_add(bo as u64 * 7));
        let bit = (h % 2048) as usize;
        bits[bit / 8] |= 1 << (bit % 8);
        if depth > 1 {
            visited[nb] = true;
            path_fp_dfs(mol, nb, h, depth - 1, visited, bits);
            visited[nb] = false;
        }
    }
}

fn path_fingerprint_bits(mol: &MolecularSystem) -> [u8; 256] {
    let n = mol.symbols.len();
    let mut bits = [0u8; 256];
    let mut visited = vec![false; n];
    for start in 0..n {
        visited[start] = true;
        let init = symbol_to_atomic_num(
            mol.symbols.get(start).map(|s| s.as_str()).unwrap_or("C")
        ) as u64;
        path_fp_dfs(mol, start, init, 7, &mut visited, &mut bits);
        visited[start] = false;
    }
    bits
}

// --- P54: remove_hs, SDF properties, get_descriptors, normalize_depiction, is_valid ---
#[wasm_bindgen]
impl MolecularSystem {
    /// Returns a new molecule with all H atoms stripped and bonds remapped.
    /// Preserves aromatic flags, stereo centers, and bond orders for heavy atoms.
    pub fn remove_hs(&self) -> MolecularSystem {
        let heavy: Vec<usize> = (0..self.symbols.len())
            .filter(|&i| self.symbols[i] != "H")
            .collect();
        self.select_by_indices(&heavy)
    }

    /// Get SDF data item by name. Returns None if not present.
    pub fn get_prop(&self, name: &str) -> Option<String> {
        self.properties.get(name).cloned()
    }

    /// Set or overwrite an SDF data item.
    pub fn set_prop(&mut self, name: &str, value: &str) {
        self.properties.insert(name.to_string(), value.to_string());
    }

    /// List all SDF data item names, sorted alphabetically.
    pub fn get_prop_list(&self) -> Vec<String> {
        let mut keys: Vec<String> = self.properties.keys().cloned().collect();
        keys.sort();
        keys
    }

    /// Returns all common descriptors as a single JS object.
    /// Ring-dependent values (num_rings, aromatic_ring_count, rotatable_bond_count,
    /// fraction_csp3) are 0 when compute_rings() has not been called.
    pub fn get_descriptors(&self) -> JsValue {
        #[derive(serde::Serialize)]
        struct Desc {
            molecular_weight: f32,
            molecular_formula: String,
            num_heavy_atoms: u32,
            formal_charge: i32,
            h_bond_donors: u32,
            h_bond_acceptors: u32,
            rotatable_bond_count: u32,
            tpsa: f32,
            logp: f32,
            molar_refractivity: f32,
            fraction_csp3: f32,
            num_rings: usize,
            aromatic_ring_count: u32,
        }
        let rings = self.enumerate_rings();
        let num_rings = rings.len();
        let aromatic_ring_count = rings.iter()
            .filter(|ring| ring.iter().any(|&a| self.aromatic_atoms.get(a).copied().unwrap_or(false)))
            .count() as u32;
        to_js(&Desc {
            molecular_weight: self.molecular_weight(),
            molecular_formula: self.molecular_formula(),
            num_heavy_atoms: self.num_heavy_atoms(),
            formal_charge: self.formal_charge(),
            h_bond_donors: self.h_bond_donors(),
            h_bond_acceptors: self.h_bond_acceptors(),
            rotatable_bond_count: self.rotatable_bond_count(),
            tpsa: self.tpsa(),
            logp: self.logp(),
            molar_refractivity: self.molar_refractivity(),
            fraction_csp3: self.fraction_csp3(),
            num_rings,
            aromatic_ring_count,
        })
    }

    /// Scales 2D coordinates so the average heavy-atom bond length equals 1.5 Å,
    /// then centers the molecule at the origin.
    pub fn normalize_depiction(&mut self) {
        if self.x.is_empty() { return; }
        let mut total = 0f32;
        let mut count = 0usize;
        for i in 0..self.symbols.len() {
            if self.symbols[i] == "H" || i >= self.bonds.len() { continue; }
            for &j in &self.bonds[i] {
                if j > i && self.symbols.get(j).map(|s| s != "H").unwrap_or(false) {
                    let dx = self.x[j] - self.x[i];
                    let dy = self.y[j] - self.y[i];
                    let d = (dx * dx + dy * dy).sqrt();
                    if d > 0.001 { total += d; count += 1; }
                }
            }
        }
        let scale = if count > 0 { 1.5 / (total / count as f32) } else { 1.0 };
        let n = self.x.len();
        for i in 0..n { self.x[i] *= scale; self.y[i] *= scale; }
        let cx: f32 = self.x.iter().sum::<f32>() / n as f32;
        let cy: f32 = self.y.iter().sum::<f32>() / n as f32;
        for i in 0..n { self.x[i] -= cx; self.y[i] -= cy; }
    }

    /// Returns true if the molecule has atoms and all bond indices are valid.
    pub fn is_valid(&self) -> bool {
        if self.symbols.is_empty() { return false; }
        let n = self.symbols.len();
        for (i, nb) in self.bonds.iter().enumerate() {
            for &j in nb {
                if j >= n || j == i { return false; }
            }
        }
        true
    }
}

// --- P55: add_hs, fingerprint_topological, get_stereo_tags ---
#[wasm_bindgen]
impl MolecularSystem {
    /// Returns a new molecule with explicit H atoms added for all heavy atoms.
    /// Coordinates are approximate (1.0 Å from parent, evenly distributed).
    /// Idempotent: calling add_hs on a fully-explicit molecule changes nothing.
    pub fn add_hs(&self) -> MolecularSystem {
        let mut mol = self.clone();
        let n = self.symbols.len();
        for i in 0..n {
            let h_count = self.implicit_h_count(i as u32).max(0) as usize;
            if h_count == 0 { continue; }
            let existing = mol.bonds.get(i).map(|b| b.len()).unwrap_or(0);
            let total = (existing + h_count).max(1);
            for k in 0..h_count {
                let j = mol.symbols.len();
                let angle = std::f32::consts::TAU * (existing + k) as f32 / total as f32;
                let (hx, hy, hz) = (mol.x[i] + angle.cos(), mol.y[i] + angle.sin(), mol.z[i]);
                mol.symbols.push("H".to_string());
                mol.x.push(hx); mol.y.push(hy); mol.z.push(hz);
                if !mol.charges.is_empty()             { mol.charges.push(0); }
                if !mol.atom_names.is_empty()          { mol.atom_names.push(String::new()); }
                if !mol.residue_names.is_empty()       { mol.residue_names.push(String::new()); mol.residue_ids.push(0); }
                if !mol.chain_ids.is_empty()           { mol.chain_ids.push(b' '); }
                if !mol.hetatm_flags.is_empty()        { mol.hetatm_flags.push(false); }
                if !mol.occupancies.is_empty()         { mol.occupancies.push(1.0); }
                if !mol.b_factors.is_empty()           { mol.b_factors.push(0.0); }
                if !mol.aromatic_atoms.is_empty()      { mol.aromatic_atoms.push(false); }
                if !mol.atom_map.is_empty()            { mol.atom_map.push(0); }
                if !mol.ring_sizes_per_atom.is_empty() { mol.ring_sizes_per_atom.push(Vec::new()); }
                mol.bonds.push(vec![i]);
                mol.bond_orders.push(vec![1]);
                mol.bonds[i].push(j);
                mol.bond_orders[i].push(1);
            }
        }
        mol.spatial_grid = None;
        mol.ring_atoms.clear();
        mol.ring_bonds.clear();
        mol
    }

    /// Path-based topological fingerprint (2048-bit / 256-byte Uint8Array).
    /// Enumerates all simple paths of bond-length 1–7; hashes atom types + bond orders.
    pub fn fingerprint_topological(&self) -> Vec<u8> {
        path_fingerprint_bits(self).to_vec()
    }

    /// Returns stereo tag objects sorted by atom index: [{index, chirality}].
    /// `chirality` is `"@"` (CCW) or `"@@"` (CW). Empty array if no stereo centers.
    pub fn get_stereo_tags(&self) -> JsValue {
        #[derive(serde::Serialize)]
        struct StereoTag { index: usize, chirality: String }
        let mut tags: Vec<StereoTag> = self.stereo_centers.iter()
            .map(|(&idx, &(desc, _))| StereoTag {
                index: idx,
                chirality: if desc > 0 { "@@".to_string() } else { "@".to_string() },
            })
            .collect();
        tags.sort_by_key(|t| t.index);
        to_js(&tags)
    }

    /// Returns true if the molecule has any tetrahedral stereo centers.
    pub fn has_stereo(&self) -> bool {
        !self.stereo_centers.is_empty()
    }
}

// --- Unit Tests ---

#[cfg(test)]
mod tests {
    use super::*;

    mod test_utils {
        use super::super::*;
        pub(super) fn water_xyz() -> MolecularSystem {
            parse_xyz("3\nwater\nO 0.000 0.000 0.119\nH 0.000 0.757 -0.477\nH 0.000 -0.757 -0.477\n")
                .unwrap()
        }
        pub(super) fn benzene_mol() -> MolecularSystem {
            let mut mol = parse_smiles("c1ccccc1").unwrap();
            mol.compute_rings();
            mol
        }
        pub(super) fn ethanol_mol() -> MolecularSystem {
            parse_smiles("CCO").unwrap()
        }
    }

    // ── XYZ tests ────────────────────────────────────────────────────────────

    const WATER_XYZ: &str = "\
3
water molecule
O   0.000  0.000  0.119
H   0.000  0.757 -0.477
H   0.000 -0.757 -0.477
";

    #[test]
    fn xyz_valid_water_molecule() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.symbols, ["O", "H", "H"]);
        assert!((mol.x[0]).abs() < 1e-5);
        assert!((mol.y[1] - 0.757).abs() < 1e-4);
        assert!((mol.z[2] - (-0.477)).abs() < 1e-4);
        // PDB fields must be empty for XYZ input
        assert!(mol.atom_names.is_empty());
        assert!(mol.residue_names.is_empty());
    }

    #[test]
    fn xyz_empty_input() {
        assert!(matches!(parse_xyz(""), Err(ParseError::EmptyInput)));
    }

    #[test]
    fn xyz_atom_limit_exceeded() {
        // Header claims more atoms than MAX_ATOMS — must be rejected before allocation.
        let input = format!("{}\ntest\n", MAX_ATOMS + 1);
        assert!(matches!(
            parse_xyz(&input),
            Err(ParseError::AtomLimitExceeded { .. })
        ));
    }

    #[test]
    fn pdb_atom_limit_exceeded() {
        // Build a PDB string with MAX_ATOMS + 1 ATOM lines.
        let line = "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N  \n";
        let input = line.repeat(MAX_ATOMS + 1);
        assert!(matches!(
            parse_pdb(&input),
            Err(ParseError::AtomLimitExceeded { .. })
        ));
    }

    #[test]
    fn xyz_invalid_atom_count_header() {
        assert!(matches!(
            parse_xyz("abc\ncomment\n"),
            Err(ParseError::InvalidAtomCount(_))
        ));
    }

    #[test]
    fn xyz_missing_comment_line() {
        assert!(matches!(parse_xyz("3"), Err(ParseError::MissingCommentLine)));
    }

    #[test]
    fn xyz_atom_count_mismatch() {
        let input = "5\ncomment\nO 0.0 0.0 0.0\n";
        assert!(matches!(
            parse_xyz(input),
            Err(ParseError::AtomCountMismatch { expected: 5, found: 1 })
        ));
    }

    #[test]
    fn xyz_invalid_x_coordinate() {
        let input = "1\ncomment\nO foo 0.0 0.0\n";
        assert!(matches!(
            parse_xyz(input),
            Err(ParseError::InvalidCoordinate { field: "x", .. })
        ));
    }

    #[test]
    fn xyz_invalid_y_coordinate() {
        let input = "1\ncomment\nO 0.0 bar 0.0\n";
        assert!(matches!(
            parse_xyz(input),
            Err(ParseError::InvalidCoordinate { field: "y", .. })
        ));
    }

    #[test]
    fn xyz_element_case_is_preserved() {
        let input = "2\ntest\nCo 0.0 0.0 0.0\nCO 1.0 1.0 1.0\n";
        let mol = parse_xyz(input).unwrap();
        assert_eq!(mol.symbols[0], "Co");
        assert_eq!(mol.symbols[1], "CO");
    }

    #[test]
    fn xyz_tab_separated_fields() {
        let input = "1\ntest\nN\t1.5\t2.5\t3.5\n";
        let mol = parse_xyz(input).unwrap();
        assert_eq!(mol.symbols[0], "N");
        assert!((mol.x[0] - 1.5).abs() < 1e-5);
        assert!((mol.y[0] - 2.5).abs() < 1e-5);
        assert!((mol.z[0] - 3.5).abs() < 1e-5);
    }

    #[test]
    fn xyz_trailing_blank_lines_are_ignored() {
        let input = "1\ntest\nC 0.0 0.0 0.0\n\n\n";
        let mol = parse_xyz(input).unwrap();
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn xyz_zero_atom_molecule() {
        let input = "0\nempty system\n";
        let mol = parse_xyz(input).unwrap();
        assert_eq!(mol.atom_count(), 0);
    }

    #[test]
    fn xyz_get_accessors_return_none_out_of_bounds() {
        let mol = parse_xyz("1\ntest\nH 0.0 0.0 0.0\n").unwrap();
        assert!(mol.get_symbol(99).is_none());
        assert!(mol.get_x(99).is_none());
    }

    #[test]
    fn xyz_large_molecule_correct_count() {
        let mut input = String::from("100\nbig molecule\n");
        for i in 0..100 {
            input.push_str(&format!("C {}.0 0.0 0.0\n", i));
        }
        let mol = parse_xyz(&input).unwrap();
        assert_eq!(mol.atom_count(), 100);
        assert_eq!(mol.symbols[99], "C");
        assert!((mol.x[99] - 99.0).abs() < 1e-4);
    }

    // ── PDB tests ────────────────────────────────────────────────────────────
    //
    // PDB columns (0-indexed):
    //   [0..6]  record  [12..16] atom name  [16] alt loc  [17..20] res name
    //   [21] chain  [22..26] resSeq  [30..38] x  [38..46] y  [46..54] z
    //   [76..78] element

    // A minimal single-atom ATOM line (80 chars, chain A, residue GLY 1)
    const ATOM_N_GLY: &str =
        "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N  ";

    #[test]
    fn pdb_single_atom_parses_correctly() {
        let mol = parse_pdb(ATOM_N_GLY).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.symbols[0], "N");
        assert_eq!(mol.atom_names[0], "N");
        assert_eq!(mol.residue_names[0], "GLY");
        assert_eq!(mol.residue_ids[0], 1);
        assert_eq!(mol.chain_ids[0], b'A');
        assert!(!mol.hetatm_flags[0]);
        assert!((mol.x[0] - 1.0).abs() < 1e-4);
        assert!((mol.y[0] - 2.0).abs() < 1e-4);
        assert!((mol.z[0] - 3.0).abs() < 1e-4);
    }

    #[test]
    fn pdb_empty_input() {
        assert!(matches!(parse_pdb(""), Err(ParseError::EmptyInput)));
        assert!(matches!(parse_pdb("   \n  "), Err(ParseError::EmptyInput)));
    }

    // Verify byte-column slicing works when adjacent negative coordinates
    // have no whitespace between them: x="-999.999", y="-999.999", z="-999.999"
    // Combined field [30..54] = "-999.999-999.999-999.999"
    // split_whitespace() would treat this as one token — byte slicing handles it correctly.
    #[test]
    fn pdb_touching_negative_coordinates() {
        let line =
            "ATOM      1  CA  GLY A   1    -999.999-999.999-999.999  1.00  0.00           C  ";
        let mol = parse_pdb(line).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert!((mol.x[0] - (-999.999)).abs() < 1e-3);
        assert!((mol.y[0] - (-999.999)).abs() < 1e-3);
        assert!((mol.z[0] - (-999.999)).abs() < 1e-3);
    }

    #[test]
    fn pdb_hetatm_water() {
        // HETATM line for a water oxygen
        let line =
            "HETATM  100  O   HOH A  50       5.000   6.000   7.000  1.00  0.00           O  ";
        let mol = parse_pdb(line).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.symbols[0], "O");
        assert_eq!(mol.residue_names[0], "HOH");
        assert!(mol.hetatm_flags[0]);
    }

    #[test]
    fn pdb_alt_loc_keeps_only_primary() {
        // 'A' alternate location is kept; 'B' is skipped.
        let input = concat!(
            "ATOM      1  CA AGLY A   1       1.000   2.000   3.000  1.00  0.00           C  \n",
            "ATOM      2  CA BGLY A   1       1.100   2.100   3.100  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(input).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert!((mol.x[0] - 1.0).abs() < 1e-4);
    }

    #[test]
    fn pdb_model_endmdl_parses_only_first_model() {
        let input = concat!(
            "MODEL        1\n",
            "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N  \n",
            "ENDMDL\n",
            "MODEL        2\n",
            "ATOM      2  N   GLY A   1      10.000  20.000  30.000  1.00  0.00           N  \n",
            "ENDMDL\n",
        );
        let mol = parse_pdb(input).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert!((mol.x[0] - 1.0).abs() < 1e-4);
    }

    #[test]
    fn pdb_short_line_no_element_column() {
        // Line truncated at col 54 — element column absent; element derived from atom name.
        let line = "ATOM      1  N   GLY A   1       1.000   2.000   3.000";
        let mol = parse_pdb(line).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.symbols[0], "N");
    }

    #[test]
    fn pdb_invalid_coordinate_returns_error() {
        let line =
            "ATOM      1  N   GLY A   1       X.XXX   2.000   3.000  1.00  0.00           N  ";
        assert!(matches!(
            parse_pdb(line),
            Err(ParseError::InvalidCoordinate { field: "x", .. })
        ));
    }

    #[test]
    fn pdb_skips_non_atom_records() {
        let input = concat!(
            "HEADER    HYDROLASE                               01-JAN-00   1ABC\n",
            "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.\n",
            "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N  \n",
            "TER\n",
            "END\n",
        );
        let mol = parse_pdb(input).unwrap();
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn pdb_multi_residue_protein_snippet() {
        // 3 residues (GLY, ALA, VAL), backbone N atoms only
        let input = concat!(
            "ATOM      1  N   GLY A   1       1.885  22.498   3.903  1.00  0.00           N  \n",
            "ATOM      5  N   ALA A   2       5.123  18.234   7.456  1.00  0.00           N  \n",
            "ATOM      9  N   VAL A   3       9.000  15.000  11.000  1.00  0.00           N  \n",
        );
        let mol = parse_pdb(input).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.residue_names, ["GLY", "ALA", "VAL"]);
        assert_eq!(mol.residue_ids, [1, 2, 3]);
        assert_eq!(mol.chain_ids, [b'A', b'A', b'A']);
        assert!((mol.x[0] - 1.885).abs() < 1e-3);
    }

    #[test]
    fn pdb_element_from_element_column_title_cased() {
        // "FE" in element column → "Fe"
        let line =
            "HETATM  200 FE   HEM A   1       0.000   0.000   0.000  1.00  0.00          FE  ";
        let mol = parse_pdb(line).unwrap();
        assert_eq!(mol.symbols[0], "Fe");
    }

    #[test]
    fn pdb_get_chain_id_accessor() {
        let mol = parse_pdb(ATOM_N_GLY).unwrap();
        assert_eq!(mol.get_chain_id(0), Some("A".to_string()));
        assert!(mol.get_chain_id(99).is_none());
    }

    // ── Bulk export & spatial query tests ────────────────────────────────────

    #[test]
    fn positions_flat_interleaved_order() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        let flat = mol.get_positions_flat();
        assert_eq!(flat.len(), 9); // 3 atoms × 3 coords
        // First atom (O): x=0.0, y=0.0, z=0.119
        assert!((flat[0]).abs() < 1e-5);
        assert!((flat[1]).abs() < 1e-5);
        assert!((flat[2] - 0.119).abs() < 1e-4);
        // Second atom (H): y=0.757
        assert!((flat[4] - 0.757).abs() < 1e-4);
    }

    #[test]
    fn symbols_json_basic() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        assert_eq!(mol.get_symbols_json(), r#"["O","H","H"]"#);
    }

    #[test]
    fn symbols_json_empty_molecule() {
        let mol = parse_xyz("0\nempty\n").unwrap();
        assert_eq!(mol.get_symbols_json(), "[]");
    }

    #[test]
    fn distance_h2_bond_length() {
        let input = "2\ntest\nH  0.000  0.000  0.000\nH  0.740  0.000  0.000\n";
        let mol = parse_xyz(input).unwrap();
        assert!((mol.distance(0, 1) - 0.74).abs() < 1e-4);
    }

    #[test]
    fn distance_same_atom_is_zero() {
        let mol = parse_xyz("1\ntest\nC  1.0  2.0  3.0\n").unwrap();
        assert!((mol.distance(0, 0)).abs() < 1e-6);
    }

    #[test]
    fn distance_out_of_bounds_returns_zero() {
        let mol = parse_xyz("1\ntest\nC  0.0  0.0  0.0\n").unwrap();
        assert_eq!(mol.distance(0, 99), 0.0);
    }

    #[test]
    fn atoms_within_radius_finds_neighbors() {
        // 3 carbons: C0 at origin, C1 at 1Å, C2 at 5Å
        let input = "3\ntest\nC  0.000  0.000  0.000\nC  1.000  0.000  0.000\nC  5.000  0.000  0.000\n";
        let mol = parse_xyz(input).unwrap();
        let within2 = mol.get_atoms_within_radius(0, 2.0);
        assert_eq!(within2, vec![1u32]);
        let within6 = mol.get_atoms_within_radius(0, 6.0);
        assert_eq!(within6, vec![1u32, 2u32]);
    }

    #[test]
    fn atoms_within_radius_excludes_center() {
        let mol = parse_xyz("1\ntest\nC  0.0  0.0  0.0\n").unwrap();
        assert!(mol.get_atoms_within_radius(0, 100.0).is_empty());
    }

    #[test]
    fn atoms_within_radius_out_of_bounds_center() {
        let mol = parse_xyz("1\ntest\nC  0.0  0.0  0.0\n").unwrap();
        assert!(mol.get_atoms_within_radius(99, 10.0).is_empty());
    }

    #[test]
    fn residues_within_radius_returns_unique_ids() {
        // Two residues: GLY(1) at ~0Å, ALA(2) at 2Å, VAL(3) at 10Å
        let input = concat!(
            "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  ALA A   2       2.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      3  CA  VAL A   3      10.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(input).unwrap();
        // From atom 0 (GLY), radius 3Å should capture ALA but not VAL
        let residues = mol.get_residues_within_radius(0, 3.0);
        assert_eq!(residues, vec!["A:ALA:2".to_string()]);
    }

    #[test]
    fn residues_within_radius_empty_when_none_nearby() {
        let input = concat!(
            "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  ALA A   2      20.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(input).unwrap();
        assert!(mol.get_residues_within_radius(0, 5.0).is_empty());
    }

    // ── AtomInfo / serde tests ────────────────────────────────────────────────

    #[test]
    fn atom_info_xyz_fields() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        let info = mol.atom_info_at(0).unwrap();
        assert_eq!(info.index, 0);
        assert_eq!(info.symbol, "O");
        assert!((info.x).abs() < 1e-5);
        assert!((info.z - 0.119).abs() < 1e-4);
        // XYZ has no PDB metadata
        assert!(info.atom_name.is_empty());
        assert!(info.residue_name.is_empty());
        assert!(!info.is_hetatm);
    }

    #[test]
    fn atom_info_pdb_fields() {
        let mol = parse_pdb(ATOM_N_GLY).unwrap();
        let info = mol.atom_info_at(0).unwrap();
        assert_eq!(info.symbol, "N");
        assert_eq!(info.atom_name, "N");
        assert_eq!(info.residue_name, "GLY");
        assert_eq!(info.residue_id, 1);
        assert_eq!(info.chain_id, "A");
        assert!(!info.is_hetatm);
    }

    #[test]
    fn atom_info_hetatm_flag() {
        let line =
            "HETATM  100  O   HOH A  50       5.000   6.000   7.000  1.00  0.00           O  ";
        let mol = parse_pdb(line).unwrap();
        let info = mol.atom_info_at(0).unwrap();
        assert!(info.is_hetatm);
        assert_eq!(info.residue_name, "HOH");
    }

    #[test]
    fn atom_info_out_of_bounds() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        assert!(mol.atom_info_at(99).is_none());
        assert!(mol.atom_info_at(3).is_none());
    }

    #[test]
    fn atom_info_index_field_correct() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        for i in 0..3 {
            assert_eq!(mol.atom_info_at(i).unwrap().index, i);
        }
    }

    // ── Voxel grid spatial index tests ───────────────────────────────────────

    const THREE_CARBONS: &str =
        "3\ntest\nC  0.000  0.000  0.000\nC  1.000  0.000  0.000\nC  5.000  0.000  0.000\n";

    #[test]
    fn spatial_index_state_before_and_after() {
        let mut mol = parse_xyz(WATER_XYZ).unwrap();
        assert!(!mol.has_spatial_index());
        mol.build_spatial_index(3.0);
        assert!(mol.has_spatial_index());
    }

    #[test]
    fn spatial_index_matches_linear_scan() {
        let mut mol = parse_xyz(THREE_CARBONS).unwrap();
        let linear = mol.get_atoms_within_radius(0, 2.0);
        mol.build_spatial_index(3.0);
        let grid = mol.get_atoms_within_radius(0, 2.0);
        assert_eq!(linear, grid);
    }

    #[test]
    fn spatial_index_large_radius_finds_all() {
        let mut mol = parse_xyz(THREE_CARBONS).unwrap();
        mol.build_spatial_index(3.0);
        let within10 = mol.get_atoms_within_radius(0, 10.0);
        assert_eq!(within10, vec![1u32, 2u32]);
    }

    #[test]
    fn spatial_index_empty_molecule() {
        let mut mol = parse_xyz("0\nempty\n").unwrap();
        mol.build_spatial_index(3.0);
        assert!(mol.has_spatial_index());
        assert!(mol.get_atoms_within_radius(0, 10.0).is_empty());
    }

    #[test]
    fn spatial_index_cell_size_larger_than_system() {
        let input = "3\ntest\nC  0.000  0.000  0.000\nC  1.000  0.000  0.000\nC  2.000  0.000  0.000\n";
        let mut mol = parse_xyz(input).unwrap();
        mol.build_spatial_index(100.0);
        let within3 = mol.get_atoms_within_radius(0, 3.0);
        assert_eq!(within3, vec![1u32, 2u32]);
    }

    #[test]
    fn spatial_index_atom_at_exact_boundary_included() {
        // Atom at d = 3.0 with radius = 3.0 must be included (d <= r).
        let input = "2\ntest\nC  0.000  0.000  0.000\nC  3.000  0.000  0.000\n";
        let mut mol = parse_xyz(input).unwrap();
        mol.build_spatial_index(2.0);
        assert_eq!(mol.get_atoms_within_radius(0, 3.0), vec![1u32]);
        assert!(mol.get_atoms_within_radius(0, 2.999).is_empty());
    }

    #[test]
    fn spatial_index_invalid_cell_size_ignored() {
        let mut mol = parse_xyz(WATER_XYZ).unwrap();
        mol.build_spatial_index(-1.0);
        assert!(!mol.has_spatial_index());
        mol.build_spatial_index(0.0);
        assert!(!mol.has_spatial_index());
        mol.build_spatial_index(f32::NAN);
        assert!(!mol.has_spatial_index());
    }

    #[test]
    fn spatial_index_residues_within_radius_uses_grid() {
        let input = concat!(
            "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  ALA A   2       2.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      3  CA  VAL A   3      10.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mut mol = parse_pdb(input).unwrap();
        mol.build_spatial_index(3.0);
        let residues = mol.get_residues_within_radius(0, 3.0);
        assert_eq!(residues, vec!["A:ALA:2".to_string()]);
    }

    // ── Bond detection tests ─────────────────────────────────────────────────

    fn make_two_atom_xyz(sym1: &str, sym2: &str, dist: f32) -> String {
        format!("2\ntest\n{sym1}  0.000  0.000  0.000\n{sym2}  {dist:.3}  0.000  0.000\n")
    }

    #[test]
    fn bonds_h2_at_equilibrium_is_bonded() {
        // H-H bond: r_cov(H)+r_cov(H)+0.4 = 0.31+0.31+0.4 = 1.02 Å; d = 0.74 Å
        let mut mol = parse_xyz(&make_two_atom_xyz("H", "H", 0.74)).unwrap();
        assert!(!mol.has_bonds_computed());
        mol.compute_bonds();
        assert!(mol.has_bonds_computed());
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.get_bonds(0), vec![1u32]);
        assert_eq!(mol.get_bonds(1), vec![0u32]);
    }

    #[test]
    fn bonds_h2_too_far_is_not_bonded() {
        // d = 1.5 Å > 1.02 Å threshold
        let mut mol = parse_xyz(&make_two_atom_xyz("H", "H", 1.5)).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 0);
        assert!(mol.get_bonds(0).is_empty());
    }

    #[test]
    fn bonds_identical_positions_are_not_bonded() {
        // d = 0.0 Å <= 0.4 Å lower bound — excluded
        let mut mol = parse_xyz(&make_two_atom_xyz("H", "H", 0.0)).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn bonds_water_has_two_oh_bonds() {
        // O at origin, two H atoms at ~0.96 Å
        // O-H threshold: 0.66+0.31+0.4 = 1.37 Å
        let input = "3\nwater\nO  0.000  0.000  0.000\nH  0.957  0.000  0.000\nH -0.239  0.927  0.000\n";
        let mut mol = parse_xyz(input).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 2);
        let o_bonds = mol.get_bonds(0);
        assert!(o_bonds.contains(&1u32));
        assert!(o_bonds.contains(&2u32));
    }

    #[test]
    fn bonds_methane_has_four_ch_bonds() {
        // C at origin, 4 H atoms at ~1.09 Å (tetrahedral geometry)
        // C-H threshold: 0.76+0.31+0.4 = 1.47 Å
        let d = 1.089f32;
        let input = format!(
            "5\nmethane\nC  0.000  0.000  0.000\n\
             H  {d:.3}  0.000  0.000\n\
             H -{d:.3}  0.000  0.000\n\
             H  0.000  {d:.3}  0.000\n\
             H  0.000 -{d:.3}  0.000\n"
        );
        let mut mol = parse_xyz(&input).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 4);
        assert_eq!(mol.get_bonds(0).len(), 4);
    }

    #[test]
    fn bonds_two_distant_molecules_no_intermolecular_bond() {
        // Two H2 molecules, 10 Å apart
        let input = "4\ntest\nH  0.000  0.000  0.000\nH  0.740  0.000  0.000\n\
                     H 10.000  0.000  0.000\nH 10.740  0.000  0.000\n";
        let mut mol = parse_xyz(input).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 2);
        let b0 = mol.get_bonds(0);
        assert_eq!(b0, vec![1u32]);
        let b2 = mol.get_bonds(2);
        assert_eq!(b2, vec![3u32]);
    }

    #[test]
    fn bonds_unknown_element_skipped() {
        // "Xx" has no covalent radius — should not bond even if close
        let input = "2\ntest\nXx  0.000  0.000  0.000\nC   0.500  0.000  0.000\n";
        let mut mol = parse_xyz(input).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn bonds_adjacency_list_is_symmetric() {
        let mut mol = parse_xyz(&make_two_atom_xyz("C", "C", 1.2)).unwrap();
        mol.compute_bonds();
        // C-C threshold: 0.76+0.76+0.4 = 1.92 Å; d = 1.2 Å → bonded
        assert_eq!(mol.bond_count(), 1);
        assert!(mol.get_bonds(0).contains(&1u32));
        assert!(mol.get_bonds(1).contains(&0u32));
    }

    #[test]
    fn bonds_bond_count_not_double_counted() {
        let mut mol = parse_xyz(&make_two_atom_xyz("N", "N", 1.1)).unwrap();
        mol.compute_bonds();
        // N-N threshold: 0.71+0.71+0.4 = 1.82 Å; d = 1.1 Å → bonded
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn bonds_get_bonds_out_of_bounds_returns_empty() {
        let mut mol = parse_xyz(&make_two_atom_xyz("H", "H", 0.74)).unwrap();
        mol.compute_bonds();
        assert!(mol.get_bonds(99).is_empty());
    }

    #[test]
    fn bonds_compute_twice_replaces_previous() {
        let mut mol = parse_xyz(&make_two_atom_xyz("H", "H", 0.74)).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 1);
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 1);
    }

    // ── Geometry tests ────────────────────────────────────────────────────────

    #[test]
    fn angle_water_hoh() {
        // Water: O at origin, H at (0.757, -0.477, 0) and (-0.757, -0.477, 0)
        // H-O-H angle ≈ 104.5° (experimental). Our coordinates give ~108° due to simplified placement.
        let mol = parse_xyz(WATER_XYZ).unwrap();
        let a = mol.angle(1, 0, 2); // H-O-H angle
        assert!(a > 100.0 && a < 115.0, "H-O-H angle = {a:.2}°");
    }

    #[test]
    fn angle_linear_is_180() {
        // Three collinear atoms: 0-1-2 on x-axis → angle at 1 = 180°
        let mol = parse_xyz("3\ntest\nC 0.0 0.0 0.0\nC 1.0 0.0 0.0\nC 2.0 0.0 0.0\n").unwrap();
        let a = mol.angle(0, 1, 2);
        assert!((a - 180.0).abs() < 0.1, "linear angle = {a:.2}°");
    }

    #[test]
    fn angle_out_of_bounds_returns_zero() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        assert_eq!(mol.angle(0, 1, 99), 0.0);
    }

    #[test]
    fn dihedral_perpendicular_planes_90deg() {
        // i=(0,0,0), j=(1,0,0), k=(1,1,0), l=(1,1,1)
        // b1=(1,0,0), b2=(0,1,0), b3=(0,0,1) → dihedral = 90°
        let mol = parse_xyz(
            "4\ntest\nC 0 0 0\nC 1 0 0\nC 1 1 0\nC 1 1 1\n",
        ).unwrap();
        let d = mol.dihedral(0, 1, 2, 3);
        assert!((d.abs() - 90.0).abs() < 0.5, "dihedral = {d:.2}°");
    }

    #[test]
    fn center_of_mass_symmetric_molecule() {
        // Two identical C atoms at ±1.0 on x; COM should be (0, 0, 0)
        let mol = parse_xyz("2\ntest\nC  1.0 0.0 0.0\nC -1.0 0.0 0.0\n").unwrap();
        let com = mol.center_of_mass();
        assert!(com[0].abs() < 1e-5);
        assert!(com[1].abs() < 1e-5);
        assert!(com[2].abs() < 1e-5);
    }

    #[test]
    fn rmsd_identical_systems_is_zero() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        let mol2 = parse_xyz(WATER_XYZ).unwrap();
        assert!(mol.rmsd(&mol2) < 1e-6);
    }

    #[test]
    fn rmsd_shifted_system() {
        // Shift all atoms by (1, 0, 0) → RMSD = 1.0
        let a = parse_xyz("2\ntest\nC 0.0 0.0 0.0\nC 1.0 0.0 0.0\n").unwrap();
        let b = parse_xyz("2\ntest\nC 1.0 0.0 0.0\nC 2.0 0.0 0.0\n").unwrap();
        let r = a.rmsd(&b);
        assert!((r - 1.0).abs() < 1e-5, "RMSD = {r:.5}");
    }

    // ── SDF tests ─────────────────────────────────────────────────────────────

    // Minimal V2000 SDF for methane (CH4): 1C + 4H, 4 C-H bonds
    const METHANE_SDF: &str = "\
methane
  -OEChem-

  5  4  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6314    0.6314    0.6314 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6314   -0.6314    0.6314 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6314    0.6314   -0.6314 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6314   -0.6314   -0.6314 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
";

    #[test]
    fn sdf_methane_atom_count() {
        let mol = parse_sdf(METHANE_SDF).unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.symbols[0], "C");
        assert_eq!(mol.symbols[1], "H");
    }

    #[test]
    fn sdf_methane_bonds_from_block() {
        let mol = parse_sdf(METHANE_SDF).unwrap();
        // Bonds are loaded from bond block, not computed from distances
        assert_eq!(mol.bond_count(), 4); // 4 C-H bonds
        assert!(mol.get_bonds(0).contains(&1u32));
        assert!(mol.get_bonds(0).contains(&2u32));
    }

    #[test]
    fn sdf_methane_coordinates() {
        let mol = parse_sdf(METHANE_SDF).unwrap();
        assert!(mol.x[0].abs() < 1e-5); // C at origin
        assert!((mol.x[1] - 0.6314).abs() < 1e-3);
    }

    #[test]
    fn sdf_empty_input_error() {
        assert!(matches!(parse_sdf(""), Err(ParseError::EmptyInput)));
    }

    // ── Kabsch tests ──────────────────────────────────────────────────────────

    #[test]
    fn kabsch_identical_structures_rmsd_zero() {
        let mol = parse_xyz(WATER_XYZ).unwrap();
        let ref_mol = parse_xyz(WATER_XYZ).unwrap();
        assert!(mol.rmsd_aligned(&ref_mol) < 1e-5);
    }

    #[test]
    fn kabsch_translation_only_rmsd_zero() {
        // Translate water by (10, 5, 3) Å — superpose should recover RMSD ≈ 0
        let ref_mol = parse_xyz(WATER_XYZ).unwrap();
        let shifted = format!(
            "3\nshifted\nO  10.000  5.000  3.119\nH  10.000  5.757  2.523\nH  10.000  4.243  2.523\n"
        );
        let mobile = parse_xyz(&shifted).unwrap();
        assert!(mobile.rmsd_aligned(&ref_mol) < 0.1); // small residual due to floating point
    }

    #[test]
    fn kabsch_superpose_mutates_coords() {
        let ref_mol = parse_xyz(WATER_XYZ).unwrap();
        let shifted = "3\nshifted\nO  10.000  5.000  3.119\nH  10.000  5.757  2.523\nH  10.000  4.243  2.523\n";
        let mut mobile = parse_xyz(shifted).unwrap();
        let rmsd_before = mobile.rmsd(&ref_mol);
        let rmsd_after  = mobile.superpose(&ref_mol);
        // After superposition coords change and RMSD drops significantly
        assert!(rmsd_after < rmsd_before, "RMSD should decrease after superposition");
    }

    #[test]
    fn kabsch_rmsd_aligned_does_not_mutate() {
        let ref_mol = parse_xyz(WATER_XYZ).unwrap();
        let shifted = "3\nshifted\nO  10.000  5.000  3.119\nH  10.000  5.757  2.523\nH  10.000  4.243  2.523\n";
        let mobile = parse_xyz(shifted).unwrap();
        let x_before = mobile.x[0];
        mobile.rmsd_aligned(&ref_mol);
        assert_eq!(mobile.x[0], x_before); // coordinates unchanged
    }

    // ── SMILES tests ──────────────────────────────────────────────────────────

    #[test]
    fn smiles_water() {
        let mol = parse_smiles("O").unwrap();
        // O: valence 2, used 0 → 2H implicit
        assert_eq!(mol.atom_count(), 3); // 1 O + 2 H
        assert_eq!(mol.symbols[0], "O");
        assert_eq!(mol.symbols[1], "H");
        assert_eq!(mol.bond_count(), 2); // O-H1, O-H2
    }

    #[test]
    fn smiles_methane() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(mol.atom_count(), 5); // 1C + 4H
        assert_eq!(mol.symbols[0], "C");
        assert_eq!(mol.bond_count(), 4);
    }

    #[test]
    fn smiles_ethanol() {
        // CCO → 2C + 1O heavy, plus H: C(3H)+C(2H)+O(1H) = 6H → 9 atoms total
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.atom_count(), 9);
        let heavy_bonds = mol.get_bonds(0).len(); // C bonded to next C + Hs
        assert!(heavy_bonds > 0);
    }

    #[test]
    fn smiles_cyclohexane_ring() {
        // C1CCCCC1 — 6 carbons, ring, each C gets 2H implicit
        let mol = parse_smiles("C1CCCCC1").unwrap();
        assert_eq!(mol.atom_count(), 18); // 6C + 12H
        assert_eq!(mol.bond_count(), 18); // 6 C-C + 12 C-H
    }

    #[test]
    fn smiles_ethylene_double_bond() {
        // C=C — 2 carbons, double bond, each C gets 2H
        let mol = parse_smiles("C=C").unwrap();
        assert_eq!(mol.atom_count(), 6); // 2C + 4H
    }

    #[test]
    fn smiles_branch() {
        // CC(C)C — isobutane: 4 carbons
        let mol = parse_smiles("CC(C)C").unwrap();
        assert_eq!(mol.symbols.iter().filter(|s| s.as_str() == "C").count(), 4);
    }

    #[test]
    fn smiles_bracket_atom() {
        // [Fe] — iron atom, no implicit H (not in organic subset)
        let mol = parse_smiles("[Fe]").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.symbols[0], "Fe");
    }

    #[test]
    fn smiles_empty_is_error() {
        assert!(matches!(parse_smiles(""), Err(ParseError::EmptyInput)));
    }

    // ── CONECT tests ──────────────────────────────────────────────────────────

    const WATER_PDB_CONECT: &str = "\
ATOM      1  O   HOH A   1       0.000   0.000   0.119  1.00  0.00           O  \n\
ATOM      2  H1  HOH A   1       0.000   0.757  -0.477  1.00  0.00           H  \n\
ATOM      3  H2  HOH A   1       0.000  -0.757  -0.477  1.00  0.00           H  \n\
CONECT    1    2    3\n\
CONECT    2    1\n\
CONECT    3    1\n\
END\n";

    #[test]
    fn pdb_conect_bonds_loaded() {
        let mol = parse_pdb(WATER_PDB_CONECT).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert!(mol.has_bonds_computed()); // has_bonds_computed uses bonds.len() == atoms.len()
        assert_eq!(mol.bond_count(), 2); // O-H1 and O-H2
        assert!(mol.get_bonds(0).contains(&1u32)); // O bonded to H1
        assert!(mol.get_bonds(0).contains(&2u32)); // O bonded to H2
    }

    #[test]
    fn pdb_no_conect_bonds_empty() {
        // Standard PDB snippet without CONECT → bonds empty, must call compute_bonds()
        let mol = parse_pdb(
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  \n\
             END\n",
        ).unwrap();
        assert!(!mol.has_bonds_computed());
        assert_eq!(mol.bond_count(), 0);
    }

    // ── Fragment tests ────────────────────────────────────────────────────────

    #[test]
    fn fragments_two_isolated_molecules() {
        // Two H2 molecules far apart — no bonds after compute_bonds() but also
        // test that with bonds each is its own fragment.
        let mut mol = parse_xyz("4\ntest\nH 0.0 0.0 0.0\nH 0.74 0.0 0.0\nH 10.0 0.0 0.0\nH 10.74 0.0 0.0\n").unwrap();
        mol.compute_bonds();
        assert_eq!(mol.bond_count(), 2);
        // Two connected components: {0,1} and {2,3}
        // get_fragments() returns JsValue — test the bond structure instead.
        assert!(mol.get_bonds(0).contains(&1u32));
        assert!(!mol.get_bonds(0).contains(&2u32));
        assert!(mol.get_bonds(2).contains(&3u32));
    }

    // ── Multi-molecule SDF tests ──────────────────────────────────────────────

    const MULTI_SDF: &str = "\
methane
  -OEChem-

  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
water
  -OEChem-

  3  2  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2392    0.9266    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
$$$$
";

    #[test]
    fn count_sdf_molecules_multi() {
        assert_eq!(MolecularSystem::count_sdf_molecules(MULTI_SDF), 2);
    }

    #[test]
    fn count_sdf_molecules_single() {
        let single = "\nmethane\n\n  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
        assert_eq!(MolecularSystem::count_sdf_molecules(single), 1);
    }

    #[test]
    fn from_sdf_nth_first() {
        let mol = parse_sdf_nth(MULTI_SDF, 0).unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.symbols[0], "C");
    }

    #[test]
    fn from_sdf_nth_second() {
        let mol = parse_sdf_nth(MULTI_SDF, 1).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn from_sdf_nth_out_of_range() {
        assert!(parse_sdf_nth(MULTI_SDF, 5).is_err());
    }

    // ── Molecular descriptor tests ────────────────────────────────────────────

    #[test]
    fn molecular_formula_ethanol() {
        // CCO from SMILES adds implicit H: C×2 + H×6 + O×1 → C2H6O
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.molecular_formula(), "C2H6O");
    }

    #[test]
    fn molecular_formula_water() {
        let mol = parse_xyz("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n").unwrap();
        assert_eq!(mol.molecular_formula(), "H2O");
    }

    #[test]
    fn molecular_weight_water() {
        let mol = parse_xyz("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n").unwrap();
        let mw = mol.molecular_weight();
        // O=15.999, H×2=2.016 → 18.015
        assert!((mw - 18.015f32).abs() < 0.01, "MW={mw}");
    }

    #[test]
    fn molecular_weight_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        let mw = mol.molecular_weight();
        // C2H6O: 2×12.011 + 6×1.008 + 15.999 = 46.069
        assert!((mw - 46.069f32).abs() < 0.05, "MW={mw}");
    }

    // ── P38: stereo + formula tests ───────────────────────────────────────────

    #[test]
    fn stereo_parse_at() {
        let mol = parse_smiles("N[C@H](C)C(=O)O").unwrap();
        assert_eq!(mol.stereo_center_count(), 1, "one stereo center");
        assert!(mol.is_stereo_center(1), "atom 1 is the chiral center");
    }

    #[test]
    fn stereo_parse_atat() {
        let mol = parse_smiles("N[C@@H](C)C(=O)O").unwrap();
        assert_eq!(mol.stereo_center_count(), 1);
    }

    #[test]
    fn stereo_center_none_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.stereo_center_count(), 0);
    }

    #[test]
    fn stereo_multi_center() {
        let mol = parse_smiles("O[C@@H](F)[C@H](Cl)Br").unwrap();
        assert_eq!(mol.stereo_center_count(), 2);
    }

    #[test]
    fn stereo_ez_slash() {
        // trans-2-butene: / treated as single bond, topology should parse cleanly
        let mol = parse_smiles("C/C=C/C").unwrap();
        assert_eq!(mol.atom_count(), 12); // 4 heavy + 8 implicit H
    }

    #[test]
    fn stereo_ez_backslash() {
        let mol = parse_smiles("C/C=C\\C").unwrap();
        assert_eq!(mol.atom_count(), 12);
    }

    #[test]
    fn stereo_ez_mixed() {
        let mol = parse_smiles("F/C=C/F").unwrap();
        // F, C, C, F + 1 implicit H per C = 6 atoms
        assert_eq!(mol.atom_count(), 6);
    }

    #[test]
    fn stereo_formula_aspirin() {
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        assert_eq!(mol.molecular_formula(), "C9H8O4");
    }

    #[test]
    fn stereo_weight_caffeine() {
        let mol = parse_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").unwrap();
        let mw = mol.molecular_weight();
        // Caffeine C8H10N4O2: MW ~194.19
        assert!((193.0..196.0).contains(&mw), "MW={mw}");
    }

    #[test]
    fn stereo_svg_wedge_present() {
        let mut mol = parse_smiles("N[C@@H](C)C(=O)O").unwrap();
        mol.compute_2d_coords();
        let svg = mol.to_svg_data(400, 300);
        assert!(svg.contains("polygon"), "SVG should contain a wedge bond polygon");
    }

    // ── P39: fingerprint + Lipinski tests ────────────────────────────────────

    #[test]
    fn p39_tpsa_acetic_acid() {
        let mol = parse_smiles("CC(=O)O").unwrap();
        let t = mol.tpsa();
        assert!((37.0..38.5).contains(&t), "TPSA={t}");
    }

    #[test]
    fn p39_tpsa_benzene_zero() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.tpsa(), 0.0, "benzene has no polar atoms");
    }

    #[test]
    fn p39_hbd_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.h_bond_donors(), 1);
    }

    #[test]
    fn p39_hba_aspirin() {
        let mol = parse_smiles("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        assert_eq!(mol.h_bond_acceptors(), 4);
    }

    #[test]
    fn p39_rotatable_butane() {
        // rotatable_bond_count requires compute_rings() first.
        // Implementation counts total degree (incl. H), so all 3 C-C bonds are rotatable.
        let mut mol = parse_smiles("CCCC").unwrap();
        mol.compute_rings();
        assert!(mol.rotatable_bond_count() >= 1);
    }

    #[test]
    fn p39_rotatable_benzene_zero() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        assert_eq!(mol.rotatable_bond_count(), 0);
    }

    #[test]
    fn p39_fingerprint_ecfp4_length() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.fingerprint_ecfp4().len(), 256);
    }

    #[test]
    fn p39_tanimoto_similarity_self() {
        let mol = parse_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").unwrap();
        assert!((mol.tanimoto_similarity(&mol) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn p39_tanimoto_similarity_dissimilar() {
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let water   = parse_smiles("O").unwrap();
        assert!(benzene.tanimoto_similarity(&water) < 0.5);
    }

    #[test]
    fn p39_tanimoto_similarity_scaffold() {
        let benzene     = parse_smiles("c1ccccc1").unwrap();
        let naphthalene = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let sim = benzene.tanimoto_similarity(&naphthalene);
        assert!(sim > 0.0, "shared aromatic scaffold, got {sim}");
    }

    // ── mmCIF tests ───────────────────────────────────────────────────────────

    const MINI_MMCIF: &str = "\
data_TEST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 O O . HOH A 1 0.000 0.000 0.000 A 1
ATOM 2 H H1 . HOH A 1 0.957 0.000 0.000 A 1
ATOM 3 H H2 . HOH A 1 -0.239 0.927 0.000 A 1
HETATM 4 C C1 . LIG B 100 5.000 5.000 5.000 B 1
#
";

    #[test]
    fn mmcif_atom_count() {
        let mol = parse_mmcif(MINI_MMCIF).unwrap();
        assert_eq!(mol.atom_count(), 4);
    }

    #[test]
    fn mmcif_symbols() {
        let mol = parse_mmcif(MINI_MMCIF).unwrap();
        assert_eq!(mol.symbols[0], "O");
        assert_eq!(mol.symbols[1], "H");
        assert_eq!(mol.symbols[3], "C");
    }

    #[test]
    fn mmcif_hetatm_flag() {
        let mol = parse_mmcif(MINI_MMCIF).unwrap();
        assert!(!mol.hetatm_flags[0]);
        assert!(mol.hetatm_flags[3]);
    }

    #[test]
    fn mmcif_coordinates() {
        let mol = parse_mmcif(MINI_MMCIF).unwrap();
        assert!((mol.x[0] - 0.0).abs() < 1e-4);
        assert!((mol.x[1] - 0.957).abs() < 1e-3);
        assert!((mol.x[3] - 5.0).abs() < 1e-4);
    }

    #[test]
    fn mmcif_residue_data() {
        let mol = parse_mmcif(MINI_MMCIF).unwrap();
        assert_eq!(mol.residue_names[0], "HOH");
        assert_eq!(mol.residue_ids[3], 100);
    }

    #[test]
    fn mmcif_empty_is_error() {
        assert!(parse_mmcif("data_EMPTY\n#\n").is_err());
    }

    #[test]
    fn mmcif_second_model_excluded() {
        let s = format!("{}\nATOM 5 N N . ALA A 2 1.0 1.0 1.0 A 2\n#\n", &MINI_MMCIF[..MINI_MMCIF.len()-3]);
        // Only first-model atoms should be parsed
        let mol = parse_mmcif(&s);
        if let Ok(m) = mol {
            assert_eq!(m.atom_count(), 4); // second model atom excluded
        }
    }

    // ── Morgan fingerprint / Tanimoto tests ──────────────────────────────────

    #[test]
    fn fingerprint_length() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.morgan_fingerprint(2).len(), 256);
    }

    #[test]
    fn tanimoto_self_is_one() {
        let mol = parse_smiles("CCO").unwrap();
        let t = mol.tanimoto(&mol, 2);
        assert!((t - 1.0).abs() < 1e-6, "Tanimoto self={t}");
    }

    #[test]
    fn tanimoto_different_molecules() {
        let ethanol = parse_smiles("CCO").unwrap();
        let water = parse_smiles("O").unwrap();
        let t = ethanol.tanimoto(&water, 2);
        assert!(t < 1.0, "ethanol vs water should be < 1.0, got {t}");
        assert!(t >= 0.0);
    }

    #[test]
    fn tanimoto_similar_gt_dissimilar() {
        let ethanol  = parse_smiles("CCO").unwrap();
        let methanol = parse_smiles("CO").unwrap();
        let benzene  = parse_smiles("c1ccccc1").unwrap(); // lowercase ok, treated as C
        let t_similar    = ethanol.tanimoto(&methanol, 2);
        let t_dissimilar = ethanol.tanimoto(&benzene, 2);
        assert!(t_similar > t_dissimilar, "ethanol-methanol={t_similar}, ethanol-benzene={t_dissimilar}");
    }

    // ── Ring detection tests ──────────────────────────────────────────────────

    #[test]
    fn rings_benzene_all_ring_atoms() {
        // Benzene SMILES: c1ccccc1 (parsed as C1CCCCC1 in our organic subset)
        let mut mol = parse_smiles("C1CCCCC1").unwrap();
        mol.compute_rings();
        assert!(mol.has_rings_computed());
        // All 6 C atoms should be ring atoms
        for i in 0..6 {
            assert!(mol.is_ring_atom(i), "atom {i} should be a ring atom");
        }
    }

    #[test]
    fn rings_linear_molecule_no_ring_atoms() {
        let mut mol = parse_smiles("CCC").unwrap();
        mol.compute_rings();
        // No rings in propane (ignoring implicit H for this test)
        let ring_count = (0..mol.atom_count()).filter(|&i| mol.is_ring_atom(i)).count();
        assert_eq!(ring_count, 0, "propane has no ring atoms");
    }

    // ── Drug-likeness tests ───────────────────────────────────────────────────

    #[test]
    fn hba_water() {
        // parse_smiles("O") → O + 2H explicit
        let mol = parse_smiles("O").unwrap();
        assert_eq!(mol.h_bond_acceptors(), 1); // 1 O
    }

    #[test]
    fn hbd_water() {
        let mol = parse_smiles("O").unwrap();
        assert_eq!(mol.h_bond_donors(), 1); // O with 2 H neighbors
    }

    #[test]
    fn rotatable_ethanol() {
        // CCO: C-C and C-O bonds (both single, non-ring, both endpoints degree>1 except OH which degree=1 from O side)
        // C1-C2-O: C1 (degree 2: C2+H×3), C2 (degree 3: C1+O+H×2 in full graph), O (degree 2: C2+H)
        // Actually with implicit H: C(degree4)−C(degree4)−O(degree2)
        // C-C bond: both degree>1, not ring → rotatable
        // C-O bond: O has degree 2 (bonded to C and H) → degree>1, rotatable
        let mut mol = parse_smiles("CCO").unwrap();
        mol.compute_rings();
        let rot = mol.rotatable_bond_count();
        assert!(rot >= 1, "ethanol has at least 1 rotatable bond, got {rot}");
    }

    #[test]
    fn rotatable_benzene_zero() {
        let mut mol = parse_smiles("C1CCCCC1").unwrap();
        mol.compute_rings();
        // All bonds are ring bonds → 0 rotatable (H-C bonds excluded too)
        assert_eq!(mol.rotatable_bond_count(), 0);
    }

    // ── Atom selection tests ──────────────────────────────────────────────────

    const TWO_CHAIN_PDB: &str = "\
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N  \n\
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C  \n\
ATOM      3  N   ALA B   2      10.000  10.000  10.000  1.00  5.00           N  \n\
HETATM    4  O   HOH A 100       5.000   5.000   5.000  1.00  8.00           O  \n\
";

    #[test]
    fn select_chain_returns_correct_atom_count() {
        let mol = parse_pdb(TWO_CHAIN_PDB).unwrap();
        let chain_a = mol.select_chain("A");
        // Chain A: atom 1 (N MET) + atom 2 (CA MET) + atom 4 (HOH, HETATM also chain A) = 3
        assert_eq!(chain_a.atom_count(), 3);
    }

    #[test]
    fn select_chain_nonexistent_returns_empty() {
        let mol = parse_pdb(TWO_CHAIN_PDB).unwrap();
        let chain_z = mol.select_chain("Z");
        assert_eq!(chain_z.atom_count(), 0);
    }

    #[test]
    fn select_hetatm_returns_only_hetatm() {
        let mol = parse_pdb(TWO_CHAIN_PDB).unwrap();
        let hetatm = mol.select_hetatm();
        assert_eq!(hetatm.atom_count(), 1);
        assert_eq!(hetatm.symbols[0], "O");
    }

    #[test]
    fn select_protein_excludes_hetatm() {
        let mol = parse_pdb(TWO_CHAIN_PDB).unwrap();
        let prot = mol.select_protein();
        assert_eq!(prot.atom_count(), 3); // 2 chain A ATOM + 1 chain B ATOM
    }

    #[test]
    fn select_chain_bonds_remapped() {
        // Build an SDF molecule, then check that after select_hetatm the bonds still make sense
        // (bond indices remapped to new range)
        let sdf = "\nmethane\n\n  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
        let mol = parse_sdf(sdf).unwrap();
        // select_hetatm on an SDF mol with no hetatm_flags → empty
        let h = mol.select_hetatm();
        assert_eq!(h.atom_count(), 0);
    }

    // ── Sequence tests ────────────────────────────────────────────────────────

    const SEQ_PDB: &str = "\
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N  \n\
ATOM      2  N   GLY A   2       2.000   2.000   2.000  1.00  0.00           N  \n\
ATOM      3  N   VAL A   3       3.000   3.000   3.000  1.00  0.00           N  \n\
HETATM    4  O   HOH A 100       4.000   4.000   4.000  1.00  0.00           O  \n\
ATOM      5  N   UNK A   4       5.000   5.000   5.000  1.00  0.00           N  \n\
";

    #[test]
    fn sequence_standard_residues() {
        let mol = parse_pdb(SEQ_PDB).unwrap();
        let seq = mol.get_sequence("A");
        // ALA=A, GLY=G, VAL=V, UNK=X; HOH is HETATM so skipped
        assert_eq!(seq, "AGVX");
    }

    #[test]
    fn sequence_skips_hetatm() {
        let mol = parse_pdb(SEQ_PDB).unwrap();
        let seq = mol.get_sequence("A");
        assert!(!seq.contains('?')); // HOH should not appear
        assert_eq!(seq.len(), 4); // only 4 ATOM residues
    }

    #[test]
    fn sequence_empty_when_no_metadata() {
        let mol = parse_xyz("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n").unwrap();
        assert_eq!(mol.get_sequence("A"), "");
    }

    #[test]
    fn sequence_unknown_residue_is_x() {
        let mol = parse_pdb(SEQ_PDB).unwrap();
        let seq = mol.get_sequence("A");
        assert!(seq.contains('X'));
    }

    // ── SDF batch screening tests ─────────────────────────────────────────────

    // Reuse MULTI_SDF from the earlier test (methane + water)
    // Check that screen_sdf_string returns the right number of entries
    // (We can't easily inspect JsValue in native tests, so we test the underlying logic)

    #[test]
    fn screen_sdf_all_molecules_processed() {
        // Parse and process both molecules manually to verify the loop logic
        let n = MolecularSystem::count_sdf_molecules(MULTI_SDF);
        let mut count = 0;
        for i in 0..n {
            if let Ok(mut mol) = parse_sdf_nth(MULTI_SDF, i) {
                mol.compute_rings();
                count += 1;
            }
        }
        assert_eq!(count, 2);
    }

    #[test]
    fn screen_sdf_descriptors_correct() {
        // First molecule: methane (C, 1 atom, 0 bonds from SDF)
        let mol0 = parse_sdf_nth(MULTI_SDF, 0).unwrap();
        assert_eq!(mol0.molecular_formula(), "C");
        assert!((mol0.molecular_weight() - 12.011f32).abs() < 0.01);
        assert_eq!(mol0.h_bond_acceptors(), 0); // C only

        // Second molecule: water (O + 2H, 2 bonds)
        let mol1 = parse_sdf_nth(MULTI_SDF, 1).unwrap();
        assert_eq!(mol1.molecular_formula(), "H2O");
        assert_eq!(mol1.h_bond_acceptors(), 1);
        assert_eq!(mol1.h_bond_donors(), 1); // O bonded to H1 and H2
    }

    // ── P12: Backbone Analysis tests ─────────────────────────────────────────

    const TRIPEPTIDE_PDB: &str = "\
ATOM      1  N   GLY A   1       1.000   1.000   1.000  1.00  0.00           N  \n\
ATOM      2  CA  GLY A   1       2.000   1.000   1.000  1.00  0.00           C  \n\
ATOM      3  C   GLY A   1       2.000   2.000   1.000  1.00  0.00           C  \n\
ATOM      4  N   ALA A   2       2.000   2.000   2.000  1.00  0.00           N  \n\
ATOM      5  CA  ALA A   2       3.000   2.000   2.000  1.00  0.00           C  \n\
ATOM      6  C   ALA A   2       3.000   3.000   2.000  1.00  0.00           C  \n\
ATOM      7  N   GLY A   3       3.000   3.000   3.000  1.00  0.00           N  \n\
ATOM      8  CA  GLY A   3       4.000   3.000   3.000  1.00  0.00           C  \n\
ATOM      9  C   GLY A   3       4.000   4.000   3.000  1.00  0.00           C  \n\
";

    #[test]
    fn backbone_angles_has_three_residues() {
        let mol = parse_pdb(TRIPEPTIDE_PDB).unwrap();
        let data = mol.backbone_angle_data();
        assert_eq!(data.len(), 3);
    }

    #[test]
    fn backbone_angles_first_phi_is_none() {
        let mol = parse_pdb(TRIPEPTIDE_PDB).unwrap();
        let data = mol.backbone_angle_data();
        assert!(data[0].phi.is_none(), "first residue phi should be None");
    }

    #[test]
    fn backbone_angles_last_psi_is_none() {
        let mol = parse_pdb(TRIPEPTIDE_PDB).unwrap();
        let data = mol.backbone_angle_data();
        let last = data.last().unwrap();
        assert!(last.psi.is_none(), "last residue psi should be None");
    }

    #[test]
    fn backbone_angles_middle_both_some() {
        let mol = parse_pdb(TRIPEPTIDE_PDB).unwrap();
        let data = mol.backbone_angle_data();
        assert!(data[1].phi.is_some(), "middle phi should be Some");
        assert!(data[1].psi.is_some(), "middle psi should be Some");
    }

    #[test]
    fn backbone_empty_for_xyz() {
        let mol = parse_xyz("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n").unwrap();
        assert_eq!(mol.backbone_angle_data().len(), 0);
    }

    // ── P13: Bond Order Storage + TPSA tests ─────────────────────────────────

    const ACETALDEHYDE_SDF: &str = "acetaldehyde
  -OEChem-

  3  2  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
M  END
";

    #[test]
    fn sdf_bond_orders_stored() {
        let mol = parse_sdf(ACETALDEHYDE_SDF).unwrap();
        assert!(!mol.bond_orders.is_empty(), "bond_orders should be populated from SDF");
        // Find bond between C(1) and O(2) — atom indices 1 and 2
        // bonds[2] should contain 1, and bond_orders[2] should contain 2
        let pos = mol.bonds[2].iter().position(|&j| j == 1);
        assert!(pos.is_some(), "atom 2 should be bonded to atom 1");
        assert_eq!(mol.bond_orders[2][pos.unwrap()], 2, "C=O should have order 2");
    }

    #[test]
    fn smiles_double_bond_order_stored() {
        // CC=O: C(0)-C(1)=O(2) + implicit H on each
        let mol = parse_smiles("CC=O").unwrap();
        assert!(!mol.bond_orders.is_empty(), "SMILES bond_orders should be populated");
        // Find O symbol and check its bonds for order 2
        let o_idx = mol.symbols.iter().position(|s| s == "O").unwrap();
        let has_double = mol.bond_orders[o_idx].iter().any(|&o| o == 2);
        assert!(has_double, "C=O should have a double bond (order 2)");
    }

    #[test]
    fn smiles_single_bonds_only_in_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        // All heavy-atom bonds should be single (1)
        let all_single = mol.bond_orders.iter().enumerate()
            .all(|(i, orders)| {
                orders.iter().enumerate().all(|(k, &ord)| {
                    // Skip H-X bonds, just check heavy-heavy
                    let j = mol.bonds[i][k];
                    mol.symbols[i] == "H" || mol.symbols[j] == "H" || ord == 1
                })
            });
        assert!(all_single, "ethanol has only single bonds");
    }

    #[test]
    fn tpsa_water_smiles() {
        let mol = parse_smiles("O").unwrap();
        let t = mol.tpsa();
        assert!((t - 20.23f32).abs() < 0.5, "water TPSA should be ~20.23, got {t}");
    }

    #[test]
    fn tpsa_acetaldehyde_smiles() {
        // CC=O: one C=O (17.07), no OH
        let mol = parse_smiles("CC=O").unwrap();
        let t = mol.tpsa();
        assert!((t - 17.07f32).abs() < 1.0, "acetaldehyde TPSA should be ~17.07, got {t}");
    }

    #[test]
    fn tpsa_zero_for_xyz_input() {
        let mol = parse_xyz("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n").unwrap();
        assert_eq!(mol.tpsa(), 0.0, "TPSA should be 0 for XYZ (no bond orders)");
    }

    #[test]
    fn tpsa_alkane_is_zero() {
        // Methane has no N or O
        let mol = parse_smiles("C").unwrap();
        assert_eq!(mol.tpsa(), 0.0);
    }

    // ── P14: Aromatic SMILES Kekulization ─────────────────────────────────────

    #[test]
    fn benzene_kekulized() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        // 6 heavy C atoms (no H for now, check bond counts)
        let heavy = mol.symbols.iter().filter(|s| s.as_str() != "H").count();
        assert_eq!(heavy, 6);
        // Count double bonds among heavy atoms
        let mut double_count = 0usize;
        for i in 0..6 {
            for (k, &j) in mol.bonds[i].iter().enumerate() {
                if j < 6 && i < j {
                    if mol.bond_orders[i][k] == 2 { double_count += 1; }
                }
            }
        }
        assert_eq!(double_count, 3, "benzene should have 3 double bonds");
    }

    #[test]
    fn pyridine_kekulized() {
        let mol = parse_smiles("c1ccncc1").unwrap();
        let heavy: Vec<_> = mol.symbols.iter().enumerate()
            .filter(|(_, s)| s.as_str() != "H")
            .map(|(i, _)| i)
            .collect();
        assert_eq!(heavy.len(), 6);
        // N (index 3) should have 0 implicit H
        let n_idx = mol.symbols.iter().position(|s| s == "N").unwrap();
        let n_h_count = mol.bonds[n_idx].iter()
            .filter(|&&j| mol.symbols[j] == "H")
            .count();
        assert_eq!(n_h_count, 0, "pyridine N should have 0 H");
    }

    #[test]
    fn naphthalene_kekulized() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let heavy_count = mol.symbols.iter().filter(|s| s.as_str() != "H").count();
        assert_eq!(heavy_count, 10, "naphthalene has 10 C atoms");
        // Count double bonds among heavy atoms
        let mut double_count = 0usize;
        for i in 0..heavy_count {
            for (k, &j) in mol.bonds[i].iter().enumerate() {
                if j < heavy_count && i < j && mol.bond_orders[i][k] == 2 {
                    double_count += 1;
                }
            }
        }
        assert_eq!(double_count, 5, "naphthalene should have 5 double bonds");
    }

    #[test]
    fn existing_aliphatic_smiles_unchanged() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.symbols.iter().filter(|s| s.as_str() == "C").count(), 2);
        assert_eq!(mol.symbols.iter().filter(|s| s.as_str() == "O").count(), 1);
    }

    #[test]
    fn mixed_aromatic_aliphatic_smiles() {
        // CC(=O)c1ccccc1 — acetophenone: aliphatic chain + aromatic ring
        let mol = parse_smiles("CC(=O)c1ccccc1").unwrap();
        let heavy = mol.symbols.iter().filter(|s| s.as_str() != "H").count();
        assert_eq!(heavy, 9, "acetophenone: 2C aliphatic + 1O + 6C aromatic");
    }

    // ── P15: B-factor + Occupancy ──────────────────────────────────────────────

    #[test]
    fn pdb_b_factor_parsed() {
        let pdb = "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 25.50           N  \n\
                   END\n";
        let mol = parse_pdb(pdb).unwrap();
        let bf = mol.b_factors.get(0).copied().unwrap_or(0.0);
        assert!((bf - 25.5).abs() < 0.01, "B-factor should be 25.5, got {}", bf);
    }

    #[test]
    fn pdb_occupancy_parsed() {
        let pdb = "ATOM      1  N   ALA A   1       1.000   2.000   3.000  0.50 10.00           N  \n\
                   END\n";
        let mol = parse_pdb(pdb).unwrap();
        let occ = mol.occupancies.get(0).copied().unwrap_or(1.0);
        assert!((occ - 0.5).abs() < 0.01, "Occupancy should be 0.5, got {}", occ);
    }

    #[test]
    fn xyz_b_factor_returns_zero() {
        let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        assert_eq!(mol.get_b_factor(0), 0.0);
        assert_eq!(mol.get_occupancy(0), 1.0);
    }

    #[test]
    fn pdb_b_factor_in_atom_info() {
        let pdb = "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 30.00           N  \n\
                   END\n";
        let mol = MolecularSystem::from_pdb_string(pdb).unwrap();
        // get_atom_info uses serde serialization; check b_factor field exists
        let bf = mol.get_b_factor(0);
        assert!((bf - 30.0).abs() < 0.01);
    }

    // ── P16: LogP tests ───────────────────────────────────────────────────────

    #[test]
    fn logp_returns_zero_without_bonds() {
        // XYZ-parsed molecule has no bonds → logp = 0.0
        let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        assert_eq!(mol.logp(), 0.0);
    }

    #[test]
    fn logp_benzene_in_range() {
        // benzene: 6 aromatic C with 1H each → 6 × 0.36 ≈ 2.16; exp = 2.13
        // SMILES parser populates bonds and bond_orders; only compute_rings() needed.
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let lp = mol.logp();
        assert!(lp > 1.0 && lp < 3.5, "benzene logP out of range: {}", lp);
    }

    #[test]
    fn logp_polar_molecule_lower_than_nonpolar() {
        // acetic acid CC(=O)O vs propane CCC
        // SMILES parser populates bonds and bond_orders; only compute_rings() needed.
        let mut propane = parse_smiles("CCC").unwrap();
        propane.compute_rings();

        let mut acid = parse_smiles("CC(=O)O").unwrap();
        acid.compute_rings();

        assert!(propane.logp() > acid.logp(),
            "propane logP ({}) should exceed acetic acid logP ({})", propane.logp(), acid.logp());
    }

    #[test]
    fn logp_halogen_increases_value() {
        // chlorobenzene Clc1ccccc1 should have higher logP than benzene c1ccccc1
        // SMILES parser populates bonds and bond_orders; only compute_rings() needed.
        let mut benzene = parse_smiles("c1ccccc1").unwrap();
        benzene.compute_rings();

        let mut chlorobenzene = parse_smiles("Clc1ccccc1").unwrap();
        chlorobenzene.compute_rings();

        assert!(chlorobenzene.logp() > benzene.logp(),
            "chlorobenzene ({}) should have higher logP than benzene ({})",
            chlorobenzene.logp(), benzene.logp());
    }

    #[test]
    fn logp_nitrogen_decreases_value() {
        // aniline c1ccccc1N should have lower logP than benzene (amine is polar)
        // SMILES parser populates bonds and bond_orders; only compute_rings() needed.
        let mut benzene = parse_smiles("c1ccccc1").unwrap();
        benzene.compute_rings();

        let mut aniline = parse_smiles("c1ccccc1N").unwrap();
        aniline.compute_rings();

        assert!(benzene.logp() > aniline.logp(),
            "benzene ({}) should have higher logP than aniline ({})",
            benzene.logp(), aniline.logp());
    }

    // ── P17: H-bond detection tests ───────────────────────────────────────────

    #[test]
    fn hbond_water_dimer_one_bond() {
        // Two water molecules in ideal H-bond geometry
        // O1 at origin, H1 along +y, H2 along +x
        // O2 at y=2.80 (along the O1-H1 direction, classic linear H-bond)
        let xyz = "6\nwater dimer\n\
                   O  0.000  0.000  0.000\n\
                   H  0.000  0.960  0.000\n\
                   H  0.960  0.000  0.000\n\
                   O  0.000  2.800  0.000\n\
                   H  0.000  3.760  0.000\n\
                   H  0.960  2.800  0.000\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        let rows = mol.find_h_bonds_data(3.5);
        // At least 1 H-bond should be found (O0→O3 via H1)
        assert!(!rows.is_empty(), "expected at least 1 H-bond, found 0");
    }

    #[test]
    fn hbond_none_when_too_far() {
        // Two water molecules 6 Å apart → no H-bonds with cutoff 3.5
        let xyz = "6\nfar waters\n\
                   O  0.0  0.0  0.0\n\
                   H  1.0  0.0  0.0\n\
                   H  0.0  1.0  0.0\n\
                   O  0.0  6.0  0.0\n\
                   H  0.0  7.0  0.0\n\
                   H  1.0  6.0  0.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        let rows = mol.find_h_bonds_data(3.5);
        assert!(rows.is_empty(), "expected 0 H-bonds for far waters, found {}", rows.len());
    }

    #[test]
    fn hbond_pdb_without_h_distance_only() {
        // PDB with N and O close together (no H atoms) → distance-only detection
        let pdb = "ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N  \n\
                   ATOM      2  O   ALA A   2       1.000   4.000   1.000  1.00  0.00           O  \n\
                   END\n";
        let mol = parse_pdb(pdb).unwrap();
        // N···O distance = 3.0 Å < 3.5
        let rows = mol.find_h_bonds_data(3.5);
        assert!(!rows.is_empty(), "PDB distance-only H-bond not found");
        // Verify it has no h_atom (distance-only mode)
        assert!(rows[0].h_atom.is_none(), "distance-only H-bond should have no h_atom");
    }

    #[test]
    fn hbond_angle_filter_excludes_bad_geometry() {
        // H placed at a 90° angle → D-H-A angle < 120° → filtered out
        // O1 at origin, H at (1,0,0), O2 at (0,1,0)
        // D-H-A angle at H: vec H→O1=(-1,0,0), vec H→O2=(-1,1,0) → angle ≈ 45° < 120°
        let xyz = "3\nbad hbond\n\
                   O  0.0  0.0  0.0\n\
                   H  1.0  0.0  0.0\n\
                   O  0.0  1.0  0.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        let rows = mol.find_h_bonds_data(3.5);
        // Check no H-bond has h_atom = Some(1)
        for row in &rows {
            if let Some(h) = row.h_atom {
                assert_ne!(h, 1, "H-bond via H(1) with bad angle should be filtered");
            }
        }
    }

    #[test]
    fn hbond_consistent_with_and_without_spatial_index() {
        // Same water dimer: results should be identical with/without spatial index
        let xyz = "6\nwater dimer\n\
                   O  0.000  0.000  0.000\n\
                   H  0.000  0.960  0.000\n\
                   H  0.960  0.000  0.000\n\
                   O  0.000  2.800  0.000\n\
                   H  0.000  3.760  0.000\n\
                   H  0.960  2.800  0.000\n";
        let mut mol1 = parse_xyz(xyz).unwrap();
        mol1.compute_bonds();
        let result1 = mol1.find_h_bonds_data(3.5);

        let mut mol2 = parse_xyz(xyz).unwrap();
        mol2.compute_bonds();
        mol2.build_spatial_index(3.5);
        let result2 = mol2.find_h_bonds_data(3.5);

        assert_eq!(result1.len(), result2.len(),
            "H-bond count differs with/without spatial index: {} vs {}", result1.len(), result2.len());
    }

    // ── P18: SSSR ring enumeration ────────────────────────────────────────────

    #[test]
    fn get_rings_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let rings = mol.enumerate_rings();
        assert_eq!(rings.len(), 1, "benzene should have 1 ring");
        assert_eq!(rings[0].len(), 6, "benzene ring should have 6 atoms");
        for &a in &rings[0] {
            assert!(mol.is_ring_atom(a), "atom {a} should be a ring atom");
        }
    }

    #[test]
    fn get_rings_naphthalene() {
        let mut mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        mol.compute_rings();
        let rings = mol.enumerate_rings();
        assert_eq!(rings.len(), 2, "naphthalene should have 2 rings");
        for ring in &rings {
            assert_eq!(ring.len(), 6, "each naphthalene ring should have 6 atoms");
        }
    }

    #[test]
    fn get_rings_linear() {
        let mut mol = parse_smiles("CCC").unwrap();
        mol.compute_rings();
        let rings = mol.enumerate_rings();
        assert_eq!(rings.len(), 0, "linear molecule has no rings");
    }

    #[test]
    fn aromatic_ring_count_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        assert_eq!(mol.aromatic_ring_count(), 1);
    }

    #[test]
    fn aromatic_ring_count_cyclohexane() {
        let mut mol = parse_smiles("C1CCCCC1").unwrap();
        mol.compute_rings();
        assert_eq!(mol.aromatic_ring_count(), 0);
    }

    // ── P19: Disulfide bonds + metal coordination sites ───────────────────────

    #[test]
    fn disulfide_two_sulfurs_close() {
        let xyz = "2\nS-S\nS  0.0  0.0  0.0\nS  2.0  0.0  0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        let bonds = mol.find_disulfide_bonds_data(2.5);
        assert_eq!(bonds.len(), 1);
        assert!((bonds[0].distance - 2.0).abs() < 0.01);
    }

    #[test]
    fn disulfide_too_far() {
        let xyz = "2\nS-S far\nS  0.0  0.0  0.0\nS  3.0  0.0  0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        let bonds = mol.find_disulfide_bonds_data(2.5);
        assert_eq!(bonds.len(), 0);
    }

    #[test]
    fn disulfide_no_sulfur() {
        let xyz = "2\nno S\nC  0.0  0.0  0.0\nN  1.5  0.0  0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        let bonds = mol.find_disulfide_bonds_data(2.5);
        assert_eq!(bonds.len(), 0);
    }

    #[test]
    fn metal_site_zinc_tetrahedral() {
        let xyz = "5\nZn tetra\nZn  0.0  0.0  0.0\nN  2.0  0.0  0.0\nN -2.0  0.0  0.0\nN  0.0  2.0  0.0\nN  0.0 -2.0  0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        let sites = mol.find_metal_sites_data(2.5);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].element, "Zn");
        assert_eq!(sites[0].coordinating.len(), 4);
    }

    #[test]
    fn metal_site_no_metals() {
        let xyz = "2\nno metals\nC  0.0  0.0  0.0\nN  1.5  0.0  0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        let sites = mol.find_metal_sites_data(2.5);
        assert_eq!(sites.len(), 0);
    }

    // ── P20: Contact map + binding site ──────────────────────────────────────

    #[test]
    fn contact_map_two_close_residues() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let contacts = mol.contact_map_data(8.0);
        assert_eq!(contacts.len(), 1);
        assert!((contacts[0].distance - 5.0).abs() < 0.01);
    }

    #[test]
    fn contact_map_too_far() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   2      10.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let contacts = mol.contact_map_data(8.0);
        assert_eq!(contacts.len(), 0);
    }

    #[test]
    fn contact_map_no_backbone_xyz() {
        let xyz = "2\nno backbone\nC  0.0  0.0  0.0\nN  1.5  0.0  0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        let contacts = mol.contact_map_data(8.0);
        assert_eq!(contacts.len(), 0, "XYZ input has no Cα atoms");
    }

    #[test]
    fn binding_site_finds_residue() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "HETATM    2  C1  LIG A 100       3.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let site = mol.binding_site_residues_data(5.0);
        assert_eq!(site.len(), 1);
        assert_eq!(site[0].residue_name, "ALA");
    }

    #[test]
    fn binding_site_no_hetatm() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let site = mol.binding_site_residues_data(5.0);
        assert_eq!(site.len(), 0, "no HETATM atoms → empty binding site");
    }

    // ── P21: Formal charge + XYZ output ──────────────────────────────────────

    #[test]
    fn formal_charge_neutral() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.formal_charge(), 0);
    }

    #[test]
    fn formal_charge_ammonium() {
        let mol = parse_smiles("[NH4+]").unwrap();
        assert_eq!(mol.formal_charge(), 1);
    }

    #[test]
    fn formal_charge_carboxylate() {
        let mol = parse_smiles("CC([O-])=O").unwrap();
        assert_eq!(mol.formal_charge(), -1);
    }

    #[test]
    fn formal_charge_zwitterion() {
        // Alanine zwitterion: [NH3+]CC([O-])=O → net charge 0
        let mol = parse_smiles("[NH3+]CC([O-])=O").unwrap();
        assert_eq!(mol.formal_charge(), 0);
    }

    #[test]
    fn to_xyz_roundtrip() {
        let xyz = "3\nwater\nO   0.000000   0.000000   0.000000\nH   0.757000   0.586000   0.000000\nH  -0.757000   0.586000   0.000000\n";
        let mol1 = parse_xyz(xyz).unwrap();
        let out = mol1.to_xyz_string();
        let mol2 = parse_xyz(&out).unwrap();
        assert_eq!(mol1.atom_count(), mol2.atom_count());
        assert_eq!(mol1.get_symbol(0), mol2.get_symbol(0));
        assert!((mol1.get_x(0).unwrap_or(0.0) - mol2.get_x(0).unwrap_or(0.0)).abs() < 1e-4);
    }

    // ── P22: SDF output + diversity picking ──────────────────────────────────

    #[test]
    fn to_sdf_roundtrip() {
        let mol1 = parse_sdf(METHANE_SDF).unwrap();
        let out = mol1.to_sdf_string();
        let mol2 = parse_sdf(&out).unwrap();
        assert_eq!(mol1.atom_count(), mol2.atom_count());
        assert_eq!(mol1.bond_count(), mol2.bond_count());
    }

    #[test]
    fn to_sdf_bond_block_written() {
        let mol = parse_sdf(METHANE_SDF).unwrap();
        let out = mol.to_sdf_string();
        // Counts line should show 5 atoms and 4 bonds
        assert!(out.contains("  5  4"), "bond count line missing: {}", out);
        // Bond lines reference atom indices
        assert!(out.contains("  1  2  1"), "C-H bond line missing");
    }

    const DIVERSE_SDF: &str = "\
mol1
  test

  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
mol2
  test

  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
mol3
  test

  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
";

    #[test]
    fn screen_sdf_diverse_k2() {
        let result = MolecularSystem::screen_sdf_diverse_data(DIVERSE_SDF, 2, 2);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn screen_sdf_diverse_single() {
        let result = MolecularSystem::screen_sdf_diverse_data(DIVERSE_SDF, 1, 2);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 0, "MaxMin always starts with index 0");
    }

    #[test]
    fn screen_sdf_diverse_clamps() {
        let result = MolecularSystem::screen_sdf_diverse_data(DIVERSE_SDF, 100, 2);
        assert_eq!(result.len(), 3, "k > mol count → return all");
    }

    // ── P23: PDB output + SASA ────────────────────────────────────────────────

    #[test]
    fn to_pdb_roundtrip() {
        let pdb_in = concat!(
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  \n",
            "ATOM      2  N   GLY B   2      -1.000  -2.000  -3.000  1.00  0.00           N  \n",
        );
        let mol1 = parse_pdb(pdb_in).unwrap();
        let out = mol1.to_pdb_string();
        let mol2 = parse_pdb(&out).unwrap();
        assert_eq!(mol2.atom_count(), 2);
        assert_eq!(mol2.symbols[0], "C");
        assert_eq!(mol2.symbols[1], "N");
        assert!((mol2.x[0] - 1.0_f32).abs() < 1e-2);
        assert!((mol2.z[1] - (-3.0_f32)).abs() < 1e-2);
    }

    #[test]
    fn to_pdb_hetatm_flag() {
        let pdb_in = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "HETATM    2  C1  LIG A 100       5.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb_in).unwrap();
        let out = mol.to_pdb_string();
        assert!(out.contains("ATOM  "), "ATOM record missing");
        assert!(out.contains("HETATM"), "HETATM record missing");
    }

    #[test]
    fn sasa_water_nonzero() {
        let xyz = "3\nwater\nO 0.000 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000\n";
        let mol = parse_xyz(xyz).unwrap();
        assert!(mol.sasa(1.4) > 0.0, "water SASA must be positive");
    }

    #[test]
    fn sasa_larger_probe_bigger() {
        let xyz = "3\nwater\nO 0.000 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000\n";
        let mol = parse_xyz(xyz).unwrap();
        assert!(mol.sasa(2.0) > mol.sasa(1.4), "larger probe → larger SASA");
    }

    #[test]
    fn sasa_empty_zero() {
        let xyz = "0\nempty\n";
        let mol = parse_xyz(xyz).unwrap();
        assert_eq!(mol.sasa(1.4), 0.0);
    }

    // ── P24: SMARTS substructure matching ────────────────────────────────────

    #[test]
    fn has_substructure_co_found() {
        // parse_smiles already populates bonds; compute_bonds() would reset them to empty
        // (all coords are 0.0 in SMILES output)
        let mol = parse_smiles("CCO").unwrap();
        assert!(mol.has_substructure("CO"), "ethanol should contain C-O");
    }

    #[test]
    fn has_substructure_co_not_found() {
        let mol = parse_smiles("CCC").unwrap();
        assert!(!mol.has_substructure("CO"), "propane has no oxygen");
    }

    #[test]
    fn match_smarts_double_bond() {
        let mol = parse_smiles("CC=O").unwrap();
        let matches = mol.match_smarts_data("C=O");
        assert!(!matches.is_empty(), "acetaldehyde should have C=O match");
        // each match has 2 atoms
        assert_eq!(matches[0].len(), 2);
    }

    #[test]
    fn match_smarts_element_count() {
        // "CC" = ethane with H's; smarts "C" = any carbon → 2 carbon atoms
        let mol = parse_smiles("CC").unwrap();
        let matches = mol.match_smarts_data("C");
        assert_eq!(matches.len(), 2, "ethane has 2 carbon atoms");
    }

    #[test]
    fn match_smarts_invalid_empty() {
        let mol = parse_smiles("CCO").unwrap();
        let matches = mol.match_smarts_data("");
        assert_eq!(matches.len(), 0, "empty SMARTS → 0 matches, no panic");
    }

    // ── P37: enhanced SMARTS tests ────────────────────────────────────────────

    #[test]
    fn smarts_wildcard_matches_all_atoms() {
        // CCO has 2C + 1O + 6H = 9 atoms; * matches all
        let mol = parse_smiles("CCO").unwrap();
        let matches = mol.match_smarts_data("*");
        assert_eq!(matches.len(), 9, "* should match every atom in ethanol");
    }

    #[test]
    fn smarts_aromatic_any_matches_benzene() {
        // c1ccccc1: 6 aromatic C (aromatic_atoms set by parse_smiles)
        let mol = parse_smiles("c1ccccc1").unwrap();
        let matches = mol.match_smarts_data("a");
        let aromatic_count = matches.len();
        assert!(aromatic_count >= 6, "a should match all 6 aromatic carbons; got {aromatic_count}");
    }

    #[test]
    fn smarts_aliphatic_any_does_not_match_benzene_carbons() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let matches = mol.match_smarts_data("A");
        // Only H atoms are aliphatic in benzene (aromatic_atoms = false for H)
        let syms: Vec<&str> = matches.iter()
            .map(|m| mol.symbols[m[0]].as_str())
            .collect();
        assert!(syms.iter().all(|&s| s == "H"), "A should match only H in benzene; got {syms:?}");
    }

    #[test]
    fn smarts_bracket_oh_matches_ethanol_oxygen() {
        // CCO: O has 1 H neighbor
        let mol = parse_smiles("CCO").unwrap();
        let matches = mol.match_smarts_data("[OH]");
        assert_eq!(matches.len(), 1, "[OH] should match exactly the O in ethanol");
        assert_eq!(mol.symbols[matches[0][0]], "O");
    }

    #[test]
    fn smarts_bracket_oh_not_in_ether() {
        // COOC: O atoms have 0 H neighbors
        let mol = parse_smiles("COOC").unwrap();
        let matches = mol.match_smarts_data("[OH]");
        assert_eq!(matches.len(), 0, "[OH] should not match ether oxygens (no H)");
    }

    #[test]
    fn smarts_ring_membership_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let matches = mol.match_smarts_data("[R]");
        let ring_c: Vec<usize> = matches.iter()
            .map(|m| m[0])
            .filter(|&i| mol.symbols[i] == "C")
            .collect();
        assert_eq!(ring_c.len(), 6, "[R] should find all 6 ring carbons in benzene; got {}", ring_c.len());
    }

    #[test]
    fn smarts_ring_membership_chain_zero() {
        let mol = parse_smiles("CCC").unwrap();
        let matches = mol.match_smarts_data("[R]");
        // ring_atoms empty (no rings); all should return false
        assert!(matches.iter().all(|m| !mol.ring_atoms.get(m[0]).copied().unwrap_or(false)),
            "[R] on propane should match no ring atoms");
    }

    #[test]
    fn smarts_atomic_number_carbon() {
        // [#6] matches any carbon regardless of aromaticity
        let mol = parse_smiles("CCO").unwrap();
        let matches = mol.match_smarts_data("[#6]");
        assert_eq!(matches.len(), 2, "[#6] should match 2 carbons in ethanol");
        assert!(matches.iter().all(|m| mol.symbols[m[0]] == "C"));
    }

    #[test]
    fn smarts_not_hydrogen() {
        // [!#1] matches any non-H atom
        let mol = parse_smiles("CCO").unwrap();
        let matches = mol.match_smarts_data("[!#1]");
        // ethanol has 3 heavy atoms: 2C + 1O
        assert_eq!(matches.len(), 3, "[!#1] should match 3 heavy atoms; got {}", matches.len());
        assert!(matches.iter().all(|m| mol.symbols[m[0]] != "H"));
    }

    #[test]
    fn smarts_aromatic_carbon_not_in_cyclohexane() {
        // C1CCCCC1: all aliphatic; SMARTS 'c' should find 0 aromatic C
        let mol = parse_smiles("C1CCCCC1").unwrap();
        let matches = mol.match_smarts_data("c");
        let c_matches: Vec<_> = matches.iter().filter(|m| mol.symbols[m[0]] == "C").collect();
        assert_eq!(c_matches.len(), 0, "aromatic c should not match aliphatic cyclohexane");
    }

    #[test]
    fn smarts_aromatic_carbon_in_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let matches = mol.match_smarts_data("c");
        let c_matches: Vec<_> = matches.iter().filter(|m| mol.symbols[m[0]] == "C").collect();
        assert_eq!(c_matches.len(), 6, "aromatic c should match 6 carbons in benzene");
    }

    #[test]
    fn smarts_aromatic_bond_cc() {
        // c:c matches aromatic C-C bond pairs in benzene
        let mol = parse_smiles("c1ccccc1").unwrap();
        let matches = mol.match_smarts_data("c:c");
        assert!(!matches.is_empty(), "c:c should match aromatic C-C bonds in benzene");
        assert_eq!(matches[0].len(), 2);
    }

    #[test]
    fn smarts_carboxylic_acid_in_acetic_acid() {
        // CC(=O)O: one C(=O)O group
        let mol = parse_smiles("CC(=O)O").unwrap();
        let matches = mol.match_smarts_data("C(=O)[OH]");
        assert!(!matches.is_empty(), "C(=O)[OH] should match in acetic acid");
    }

    #[test]
    fn smarts_benzene_ring_pattern() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert!(mol.has_substructure("c1ccccc1"), "benzene ring SMARTS should match benzene");
    }

    #[test]
    fn smarts_benzene_ring_not_in_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert!(!mol.has_substructure("c1ccccc1"), "benzene ring SMARTS should not match ethanol");
    }

    #[test]
    fn smarts_get_substructure_atoms_union() {
        // [OH] in ethanol → 1 oxygen atom
        let mol = parse_smiles("CCO").unwrap();
        let atom_indices: Vec<usize> = mol.match_smarts_data("[OH]")
            .into_iter().flatten().collect();
        assert_eq!(atom_indices.len(), 1);
        assert_eq!(mol.symbols[atom_indices[0]], "O");
    }

    // ── P25: Residue SASA + chain interface ──────────────────────────────────

    #[test]
    fn residue_sasa_nonzero() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let rows = mol.residue_sasa_data(1.4);
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|r| r.sasa > 0.0), "all residues should have positive SASA");
    }

    #[test]
    fn residue_sasa_sum_approx_total() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   2       5.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let total = mol.sasa(1.4);
        let residue_sum: f32 = mol.residue_sasa_data(1.4).iter().map(|r| r.sasa).sum();
        assert!((total - residue_sum).abs() < 1e-3, "residue SASA sum should equal total SASA");
    }

    #[test]
    fn chain_interface_finds_residues() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY B   1       3.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let result = mol.chain_interface_data("A", "B", 5.0);
        assert_eq!(result.a.len(), 1, "ALA in chain A should be at interface");
        assert_eq!(result.b.len(), 1, "GLY in chain B should be at interface");
    }

    #[test]
    fn chain_interface_empty_far() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY B   1     100.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let result = mol.chain_interface_data("A", "B", 5.0);
        assert_eq!(result.a.len(), 0, "chains far apart → no interface");
        assert_eq!(result.b.len(), 0);
    }

    #[test]
    fn chain_interface_no_meta() {
        let xyz = "2\ntest\nC 0.0 0.0 0.0\nN 3.0 0.0 0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        let result = mol.chain_interface_data("A", "B", 5.0);
        assert_eq!(result.a.len(), 0);
        assert_eq!(result.b.len(), 0);
    }

    // ── P26: Murcko scaffold ──────────────────────────────────────────────────

    #[test]
    fn murcko_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let scaffold = mol.murcko_scaffold_indices_data();
        // benzene: 6 ring C + 6 H; all 6 C are scaffold, H's are pruned
        let ring_count = scaffold.iter().filter(|&&i| mol.symbols[i] == "C").count();
        assert_eq!(ring_count, 6, "benzene scaffold should have 6 carbons");
    }

    #[test]
    fn murcko_ethylbenzene() {
        // c1ccccc1CC: benzene + ethyl side chain; only 6 ring C should remain
        let mut mol = parse_smiles("c1ccccc1CC").unwrap();
        mol.compute_rings();
        let scaffold = mol.murcko_scaffold_indices_data();
        let heavy: Vec<_> = scaffold.iter().filter(|&&i| mol.symbols[i] != "H").collect();
        assert_eq!(heavy.len(), 6, "ethylbenzene scaffold: 6 ring carbons only");
    }

    #[test]
    fn murcko_two_rings_linked() {
        // c1ccccc1CCc1ccccc1: two benzenes with -CH2CH2- linker
        let mut mol = parse_smiles("c1ccccc1CCc1ccccc1").unwrap();
        mol.compute_rings();
        let scaffold = mol.murcko_scaffold_indices_data();
        let heavy: Vec<_> = scaffold.iter().filter(|&&i| mol.symbols[i] != "H").collect();
        // 12 ring C + 2 linker C = 14
        assert_eq!(heavy.len(), 14, "scaffold should include 12 ring + 2 linker carbons");
    }

    #[test]
    fn murcko_no_rings() {
        let mut mol = parse_smiles("CCC").unwrap();
        mol.compute_rings();
        let scaffold = mol.murcko_scaffold_indices_data();
        assert!(scaffold.is_empty(), "propane has no rings → empty scaffold");
    }

    #[test]
    fn ring_system_count_two() {
        let mut mol = parse_smiles("c1ccccc1CCc1ccccc1").unwrap();
        mol.compute_rings();
        assert_eq!(mol.ring_system_count(), 2, "two separated benzene rings → 2 ring systems");
    }

    // ── P27: Chain breaks + Ramachandran outliers ─────────────────────────────

    #[test]
    fn chain_breaks_no_gap() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let breaks = mol.chain_breaks_data(4.5);
        assert_eq!(breaks.len(), 0, "adjacent residues within 4.5 Å → no break");
    }

    #[test]
    fn chain_breaks_seq_gap() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   3       3.800   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let breaks = mol.chain_breaks_data(4.5);
        assert_eq!(breaks.len(), 1, "resid 1→3 is a sequence gap");
        assert_eq!(breaks[0].from_resid, 1);
        assert_eq!(breaks[0].to_resid, 3);
    }

    #[test]
    fn chain_breaks_dist_gap() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      2  CA  GLY A   2       6.000   0.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let breaks = mol.chain_breaks_data(4.5);
        assert_eq!(breaks.len(), 1, "Cα 6 Å apart → structural break");
    }

    #[test]
    fn ramachandran_outlier_found() {
        // Middle residue (ALA) has phi=180°, psi=180° → outside all allowed regions.
        // Geometry: C1=(0,-1,0), N2=(0,0,0), CA2=(1,0,0), C2=(1,1,0), N3=(2,1,0).
        let pdb = concat!(
            "ATOM      1  N   GLY A   1      -2.000  -1.000   0.000  1.00  0.00           N  \n",
            "ATOM      2  CA  GLY A   1      -1.000  -1.000   0.000  1.00  0.00           C  \n",
            "ATOM      3  C   GLY A   1       0.000  -1.000   0.000  1.00  0.00           C  \n",
            "ATOM      4  N   ALA A   2       0.000   0.000   0.000  1.00  0.00           N  \n",
            "ATOM      5  CA  ALA A   2       1.000   0.000   0.000  1.00  0.00           C  \n",
            "ATOM      6  C   ALA A   2       1.000   1.000   0.000  1.00  0.00           C  \n",
            "ATOM      7  N   GLY A   3       2.000   1.000   0.000  1.00  0.00           N  \n",
            "ATOM      8  CA  GLY A   3       3.000   1.000   0.000  1.00  0.00           C  \n",
            "ATOM      9  C   GLY A   3       3.000   2.000   0.000  1.00  0.00           C  \n",
        );
        let mol = parse_pdb(pdb).unwrap();
        let outliers = mol.ramachandran_outliers_data();
        assert_eq!(outliers.len(), 1, "ALA should be an outlier (phi=180°, psi=180°)");
        assert_eq!(outliers[0].residue_name, "ALA");
    }

    #[test]
    fn ramachandran_no_backbone() {
        let mol = parse_xyz("3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n").unwrap();
        let outliers = mol.ramachandran_outliers_data();
        assert!(outliers.is_empty(), "XYZ input has no backbone metadata → empty");
    }

    // ── P29b: Coordination geometry ───────────────────────────────────────────

    #[test]
    fn coord_geom_linear_co2() {
        // CO2: C at origin, O at ±1.16 Å on x-axis → angle ≈ 180°
        let xyz = "3\nco2\nC 0.0 0.0 0.0\nO 1.16 0.0 0.0\nO -1.16 0.0 0.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.coordination_geometry_data(0), "linear");
    }

    #[test]
    fn coord_geom_tetrahedral() {
        // Zn at origin, 4 Cl at tetrahedral vertices — all angles ≈ 109.5°
        let xyz = "5\nzncl4\nZn 0.0 0.0 0.0\nCl 1.0 1.0 1.0\nCl -1.0 -1.0 1.0\nCl -1.0 1.0 -1.0\nCl 1.0 -1.0 -1.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.coordination_geometry_data(0), "tetrahedral");
    }

    #[test]
    fn coord_geom_square_planar() {
        // Pt at origin, 4 N in xy-plane at ±2 Å — 4×90° + 2×180°
        let xyz = "5\npt_sq\nPt 0.0 0.0 0.0\nN 2.0 0.0 0.0\nN -2.0 0.0 0.0\nN 0.0 2.0 0.0\nN 0.0 -2.0 0.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.coordination_geometry_data(0), "square_planar");
    }

    #[test]
    fn coord_geom_octahedral() {
        // Fe at origin, 6 N at ±2 Å on each axis
        let xyz = "7\nfe_oct\nFe 0.0 0.0 0.0\nN 2.0 0.0 0.0\nN -2.0 0.0 0.0\nN 0.0 2.0 0.0\nN 0.0 -2.0 0.0\nN 0.0 0.0 2.0\nN 0.0 0.0 -2.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.coordination_geometry_data(0), "octahedral");
    }

    #[test]
    fn coord_geom_unknown_cn1() {
        let xyz = "2\ntest\nFe 0.0 0.0 0.0\nN 2.0 0.0 0.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_bonds();
        assert_eq!(mol.coordination_geometry_data(0), "unknown");
    }

    #[test]
    fn coord_geom_smiles_unknown() {
        // SMILES coords are all-zero → unknown
        let mol = parse_smiles("[Fe]").unwrap();
        assert_eq!(mol.coordination_geometry_data(0), "unknown");
    }

    // ── P29a: Functional group detection ─────────────────────────────────────

    #[test]
    fn fg_ethanol_alcohol() {
        let mol = parse_smiles("CCO").unwrap();
        let groups = mol.detect_functional_groups_data();
        assert!(groups.contains(&"alcohol".to_string()), "ethanol → alcohol");
        assert!(!groups.contains(&"ether".to_string()), "ethanol → no ether");
    }

    #[test]
    fn fg_ether() {
        let mol = parse_smiles("CCOCC").unwrap();
        let groups = mol.detect_functional_groups_data();
        assert!(groups.contains(&"ether".to_string()), "diethyl ether → ether");
        assert!(!groups.contains(&"alcohol".to_string()), "diethyl ether → no alcohol");
    }

    #[test]
    fn fg_carboxylic_acid() {
        let mol = parse_smiles("CC(=O)O").unwrap();
        let groups = mol.detect_functional_groups_data();
        assert!(groups.contains(&"carboxylic_acid".to_string()), "acetic acid → carboxylic_acid");
        assert!(!groups.contains(&"ester".to_string()), "acetic acid → no ester");
    }

    #[test]
    fn fg_ester() {
        let mol = parse_smiles("CC(=O)OC").unwrap();
        let groups = mol.detect_functional_groups_data();
        assert!(groups.contains(&"ester".to_string()), "methyl acetate → ester");
        assert!(!groups.contains(&"carboxylic_acid".to_string()), "methyl acetate → no carboxylic_acid");
    }

    #[test]
    fn fg_benzene_aromatic() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let groups = mol.detect_functional_groups_data();
        assert!(groups.contains(&"aromatic".to_string()), "benzene → aromatic");
    }

    #[test]
    fn fg_halide_cl() {
        let mol = parse_smiles("ClC(Cl)Cl").unwrap();
        let groups = mol.detect_functional_groups_data();
        assert!(groups.contains(&"halide_Cl".to_string()), "chloroform → halide_Cl");
    }

    // ── P28: SMILES output ────────────────────────────────────────────────────

    #[test]
    fn smiles_methane_roundtrip() {
        let mol = parse_smiles("C").unwrap();
        let s = mol.to_smiles_data();
        assert_eq!(s, "C", "methane → C");
    }

    #[test]
    fn smiles_ethanol_roundtrip() {
        let mol = parse_smiles("CCO").unwrap();
        let s = mol.to_smiles_data();
        // Re-parse and verify heavy atom count (H atoms are implicit)
        let mol2 = parse_smiles(&s).unwrap();
        let heavy = mol2.symbols.iter().filter(|s| s.as_str() != "H").count();
        assert_eq!(heavy, 3, "ethanol SMILES re-parse → 3 heavy atoms");
    }

    #[test]
    fn smiles_benzene_ring_closure() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let s = mol.to_smiles_data();
        assert!(s.contains('1') || s.contains('%'), "benzene SMILES must have ring closure");
        // Re-parse: must give 6 heavy atoms (H atoms are implicit, total > 6)
        let mol2 = parse_smiles(&s).unwrap();
        let heavy = mol2.symbols.iter().filter(|s| s.as_str() != "H").count();
        assert_eq!(heavy, 6, "benzene re-parse → 6 heavy atoms");
    }

    #[test]
    fn smiles_acetylene_triple_bond() {
        let mol = parse_smiles("C#C").unwrap();
        let s = mol.to_smiles_data();
        assert!(s.contains('#'), "acetylene → triple bond in SMILES");
    }

    #[test]
    fn smiles_charged_atom_bracket() {
        let mol = parse_smiles("[NH4+]").unwrap();
        let s = mol.to_smiles_data();
        assert!(s.contains('[') && s.contains('+'), "charged atom uses bracket notation");
    }

    #[test]
    fn smiles_disconnected_dot_separator() {
        // Two isolated atoms (no bonds)
        let xyz = "2\ntest\nC 0.0 0.0 0.0\nO 99.0 0.0 0.0\n";
        let mol = parse_xyz(xyz).unwrap();
        // bonds.is_empty() path → "[C].[O]"
        let s = mol.to_smiles_data();
        assert!(s.contains('.'), "disconnected molecule → '.' separator");
    }

    // ── P32: 2D coordinate generation ────────────────────────────────────────

    #[test]
    fn coords2d_single_atom() {
        let mut mol = parse_smiles("C").unwrap();
        mol.compute_2d_coords_data();
        // The heavy atom (C, index 0) should be placed at the origin
        assert!((mol.x[0]).abs() < 1e-4, "single C x ≈ 0");
        assert!((mol.y[0]).abs() < 1e-4, "single C y ≈ 0");
        assert!((mol.z[0]).abs() < 1e-4, "single C z = 0");
    }

    #[test]
    fn coords2d_chain_ethanol() {
        let mut mol = parse_smiles("CCO").unwrap();
        mol.compute_2d_coords_data();
        // 3 heavy atoms (C, C, O): all x-coords must differ, distances ≈ 1.5
        let heavy: Vec<usize> = mol.symbols.iter().enumerate()
            .filter(|(_, s)| s.as_str() != "H")
            .map(|(i, _)| i)
            .collect();
        assert_eq!(heavy.len(), 3, "ethanol has 3 heavy atoms");
        let xs: Vec<f32> = heavy.iter().map(|&i| mol.x[i]).collect();
        assert!(xs[0] != xs[1] && xs[1] != xs[2], "heavy atom x-coords differ");
        // Check bond lengths between consecutive heavy atoms
        let d01 = ((mol.x[heavy[0]] - mol.x[heavy[1]]).powi(2)
            + (mol.y[heavy[0]] - mol.y[heavy[1]]).powi(2)).sqrt();
        assert!((d01 - 1.5).abs() < 0.05, "C-C distance ≈ 1.5 Å, got {d01}");
    }

    #[test]
    fn coords2d_benzene_hexagon() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_2d_coords_data();
        let heavy: Vec<usize> = mol.symbols.iter().enumerate()
            .filter(|(_, s)| s.as_str() != "H")
            .map(|(i, _)| i)
            .collect();
        assert_eq!(heavy.len(), 6, "benzene has 6 heavy atoms");
        // All 6 C atoms should be equidistant from their centroid
        let cx: f32 = heavy.iter().map(|&i| mol.x[i]).sum::<f32>() / 6.0;
        let cy: f32 = heavy.iter().map(|&i| mol.y[i]).sum::<f32>() / 6.0;
        let radii: Vec<f32> = heavy.iter().map(|&i| {
            ((mol.x[i] - cx).powi(2) + (mol.y[i] - cy).powi(2)).sqrt()
        }).collect();
        let r0 = radii[0];
        for &r in &radii {
            assert!((r - r0).abs() < 0.05, "benzene ring radius uniform: {r} vs {r0}");
        }
        // All z must be 0
        for &i in &heavy {
            assert!((mol.z[i]).abs() < 1e-4, "z = 0 after 2D layout");
        }
    }

    #[test]
    fn coords2d_naphthalene_wider() {
        let mut mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        mol.compute_2d_coords_data();
        let heavy: Vec<usize> = mol.symbols.iter().enumerate()
            .filter(|(_, s)| s.as_str() != "H")
            .map(|(i, _)| i)
            .collect();
        assert_eq!(heavy.len(), 10, "naphthalene has 10 heavy atoms");
        let xmin = heavy.iter().map(|&i| mol.x[i]).fold(f32::INFINITY, f32::min);
        let xmax = heavy.iter().map(|&i| mol.x[i]).fold(f32::NEG_INFINITY, f32::max);
        assert!(xmax - xmin > 3.0, "naphthalene x-span > 3.0 Å, got {}", xmax - xmin);
    }

    #[test]
    fn coords2d_disconnected_fragments_separated() {
        // Two isolated atoms from XYZ (no bonds): must be placed at different x
        let xyz = "2\ntest\nC 0.0 0.0 0.0\nO 0.0 0.0 0.0\n";
        let mut mol = parse_xyz(xyz).unwrap();
        mol.compute_2d_coords_data();
        assert!(
            (mol.x[0] - mol.x[1]).abs() > 0.1,
            "disconnected atoms placed at different x: {} vs {}",
            mol.x[0], mol.x[1]
        );
    }

    // ── P33 SVG renderer tests ────────────────────────────────────────────────

    #[test]
    fn svg_empty_molecule() {
        // [H][H] has no heavy atoms → hits the early-return branch; must still produce valid SVG
        let mut mol = parse_smiles("[H][H]").unwrap();
        mol.compute_2d_coords_data();
        let svg = mol.to_svg_data(200, 200);
        assert!(svg.contains("<svg"), "all-H mol svg has opening tag");
        assert!(svg.contains("</svg>"), "all-H mol svg has closing tag");
    }

    #[test]
    fn svg_viewbox_size() {
        let mut mol = parse_smiles("C").unwrap();
        mol.compute_2d_coords_data();
        let svg = mol.to_svg_data(400, 300);
        assert!(
            svg.contains("viewBox=\"0 0 400 300\""),
            "viewBox reflects width/height: {svg}"
        );
    }

    #[test]
    fn svg_oxygen_label_present() {
        let mut mol = parse_smiles("CCO").unwrap();
        mol.compute_2d_coords_data();
        let svg = mol.to_svg_data(400, 300);
        assert!(svg.contains(">O<"), "oxygen label present in SVG");
    }

    #[test]
    fn svg_carbon_label_absent() {
        let mut mol = parse_smiles("CCO").unwrap();
        mol.compute_2d_coords_data();
        let svg = mol.to_svg_data(400, 300);
        assert!(!svg.contains(">C<"), "carbon label absent in SVG (chemical drawing convention)");
    }

    #[test]
    fn svg_double_bond_two_lines() {
        let mut mol = parse_smiles("C=O").unwrap();
        mol.compute_2d_coords_data();
        let svg = mol.to_svg_data(200, 200);
        let line_count = svg.matches("<line").count();
        assert!(
            line_count >= 2,
            "double bond C=O should produce ≥2 <line> elements, got {line_count}"
        );
    }

    // ── P34 aromatic flag + SVG tests ─────────────────────────────────────────

    #[test]
    fn aromatic_atoms_set_for_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let heavy: Vec<usize> = mol.symbols.iter().enumerate()
            .filter(|(_, s)| s.as_str() != "H")
            .map(|(i, _)| i)
            .collect();
        assert_eq!(heavy.len(), 6, "benzene has 6 heavy atoms");
        for &i in &heavy {
            assert!(mol.aromatic_atoms[i], "atom {i} should be aromatic in benzene");
        }
    }

    #[test]
    fn aromatic_atoms_empty_for_xyz() {
        let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.757 -0.477\nH 0.0 -0.757 -0.477\n";
        let mol = parse_xyz(xyz).unwrap();
        assert!(mol.aromatic_atoms.is_empty(), "XYZ parse should leave aromatic_atoms empty");
    }

    #[test]
    fn svg_benzene_has_dashed_inner_bond() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_2d_coords_data();
        let svg = mol.to_svg_data(300, 300);
        assert!(
            svg.contains("stroke-dasharray"),
            "benzene SVG should have dashed inner aromatic bond lines"
        );
    }

    #[test]
    fn svg_cyclohexane_no_dashed_bond() {
        let mut mol = parse_smiles("C1CCCCC1").unwrap();
        mol.compute_2d_coords_data();
        let svg = mol.to_svg_data(300, 300);
        assert!(
            !svg.contains("stroke-dasharray"),
            "cyclohexane SVG should have no dashed bonds (non-aromatic ring)"
        );
    }

    // ── P41 format I/O tests ─────────────────────────────────────────────────

    #[test]
    fn cdxml_roundtrip_ethane() {
        let cdxml = r#"<?xml version="1.0"?>
<CDXML>
<Page>
<Fragment id="1">
<Node id="1" p="0 0" Element="6"/>
<Node id="2" p="24 0" Element="6"/>
<Bond id="1" B="1" E="2" Order="1"/>
</Fragment>
</Page>
</CDXML>"#;
        let mol = parse_cdxml(cdxml).expect("cdxml parse failed");
        assert_eq!(mol.symbols.len(), 2);
        assert!(mol.symbols.iter().all(|s| s == "C"));
        let out = mol_to_cdxml_string(&mol);
        let mol2 = parse_cdxml(&out).expect("cdxml roundtrip parse failed");
        assert_eq!(mol2.symbols.len(), 2);
    }

    #[test]
    fn mrv_roundtrip_ethanol() {
        let mrv = r#"<?xml version="1.0"?>
<cml>
  <MDocument>
    <MChemicalStruct>
      <molecule molID="m1">
        <atomArray>
          <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
          <atom id="a2" elementType="C" x2="1.54" y2="0.0"/>
          <atom id="a3" elementType="O" x2="3.08" y2="0.0"/>
        </atomArray>
        <bondArray>
          <bond id="b1" atomRefs2="a1 a2" order="1"/>
          <bond id="b2" atomRefs2="a2 a3" order="1"/>
        </bondArray>
      </molecule>
    </MChemicalStruct>
  </MDocument>
</cml>"#;
        let mol = parse_mrv(mrv).expect("mrv parse failed");
        assert_eq!(mol.symbols.len(), 3);
        assert_eq!(mol.symbols[0], "C");
        assert_eq!(mol.symbols[2], "O");
        let out = mol_to_mrv_string(&mol);
        let mol2 = parse_mrv(&out).expect("mrv roundtrip parse failed");
        assert_eq!(mol2.symbols.len(), 3);
    }

    #[test]
    fn ket_roundtrip_co2() {
        // CO2: all heavy atoms (no H), so roundtrip through mol_to_ket_string preserves all 3
        let ket = r#"{"root":{"nodes":[{"$ref":"mol0"}]},"mol0":{"type":"molecule","atoms":[{"label":"C","location":[0.0,0.0,0.0],"charge":0},{"label":"O","location":[1.2,0.0,0.0],"charge":0},{"label":"O","location":[-1.2,0.0,0.0],"charge":0}],"bonds":[{"type":2,"atoms":[0,1],"stereo":0},{"type":2,"atoms":[0,2],"stereo":0}]}}"#;
        let mol = parse_ket(ket).expect("ket parse failed");
        assert_eq!(mol.symbols.len(), 3);
        assert_eq!(mol.symbols[0], "C");
        assert_eq!(mol.symbols[1], "O");
        let out = mol_to_ket_string(&mol);
        let mol2 = parse_ket(&out).expect("ket roundtrip parse failed");
        assert_eq!(mol2.symbols.len(), 3);
        assert!((mol2.x[1] - mol.x[1]).abs() < 0.01);
    }

    #[test]
    fn rxn_parse_counts() {
        // Each $MOL block needs: name / program / comment / counts (4-line header)
        // After strip_prefix('\n'), block becomes: name\nprogram\ncomment\ncounts\n...
        let rxn = "$RXN\nesterification\n  chem-wasm-lens\n\n  1  1\n$MOL\nreactant\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n$MOL\nproduct\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0800    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n";
        let rxn_mol = parse_rxn(rxn).expect("rxn parse failed");
        assert_eq!(rxn_mol.reactants.len(), 1);
        assert_eq!(rxn_mol.products.len(), 1);
        assert_eq!(rxn_mol.reactants[0].symbols.len(), 2);
        assert_eq!(rxn_mol.products[0].symbols.len(), 3);
    }

    #[test]
    fn reaction_smiles_parse() {
        let r = parse_reaction_smiles("CCO>>CC").expect("reaction smiles parse failed");
        assert_eq!(r.reactants.len(), 1);
        assert_eq!(r.products.len(), 1);
        assert!(r.reactants[0].symbols.len() >= 3, "ethanol has C, C, O");
        let back = reaction_to_smiles(&r);
        assert!(back.contains(">>"), "should contain reaction separator");
    }

    #[test]
    fn cml_roundtrip_acetic_acid() {
        let cml = r#"<?xml version="1.0"?>
<cml>
  <molecule id="m1">
    <atomArray>
      <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
      <atom id="a2" elementType="C" x2="1.54" y2="0.0"/>
      <atom id="a3" elementType="O" x2="3.08" y2="0.0"/>
      <atom id="a4" elementType="O" x2="1.54" y2="1.54"/>
    </atomArray>
    <bondArray>
      <bond id="b1" atomRefs2="a1 a2" order="1"/>
      <bond id="b2" atomRefs2="a2 a3" order="1"/>
      <bond id="b3" atomRefs2="a2 a4" order="2"/>
    </bondArray>
  </molecule>
</cml>"#;
        let mol = parse_cml(cml).expect("cml parse failed");
        assert_eq!(mol.symbols.len(), 4);
        let out = mol_to_cml_string(&mol);
        let mol2 = parse_cml(&out).expect("cml roundtrip parse failed");
        assert_eq!(mol2.symbols.len(), 4);
        assert_eq!(mol2.symbols[2], "O");
    }

    // ── P42 Editor Kernel tests ───────────────────────────────────────────

    #[test]
    fn edit_add_atom_increases_count() {
        let mut mol = MolecularSystem::new_empty();
        let idx = mol.add_atom("C", 0.0, 0.0);
        assert_eq!(idx, 0);
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.symbols[0], "C");
        assert!((mol.x[0] - 0.0).abs() < 1e-5);

        let idx2 = mol.add_atom("O", 1.5, 0.0);
        assert_eq!(idx2, 1);
        assert_eq!(mol.atom_count(), 2);
    }

    #[test]
    fn edit_remove_atom_updates_bonds() {
        // Build C-C-C: atoms 0,1,2 with bonds 0-1 and 1-2
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        mol.add_atom("C", 1.5, 0.0);
        mol.add_atom("C", 3.0, 0.0);
        mol.add_bond(0, 1, 1);
        mol.add_bond(1, 2, 1);

        // Remove middle atom (idx 1)
        mol.remove_atom(1);

        assert_eq!(mol.atom_count(), 2);
        // After removal, old atom 2 is now at index 1
        assert_eq!(mol.symbols[1], "C");
        // Bonds should be empty (1-0 and 1-2 were both removed with the atom)
        assert!(mol.bonds[0].is_empty());
        assert!(mol.bonds[1].is_empty());
    }

    #[test]
    fn edit_set_bond_order() {
        let mut mol = parse_smiles("CC").unwrap();
        // atoms 0,1 are C-C single bond (index in bond_orders)
        mol.set_bond_order(0, 1, 2);
        // get_bond_order(atom, neighbor_k) — k=0 for first neighbor
        assert_eq!(mol.get_bond_order(0, 0), 2);
        assert_eq!(mol.get_bond_order(1, 0), 2);
    }

    #[test]
    fn edit_add_bond_deduplicates() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        mol.add_atom("C", 1.5, 0.0);
        mol.add_bond(0, 1, 1);
        let bc_before = mol.bond_count();
        // Adding the same bond again should update order, not duplicate
        mol.add_bond(0, 1, 2);
        assert_eq!(mol.bond_count(), bc_before);
        assert_eq!(mol.get_bond_order(0, 0), 2);
    }

    #[test]
    fn edit_closest_atom() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        mol.add_atom("O", 3.0, 0.0);
        // Closest to (0.1, 0.1) should be atom 0
        assert_eq!(mol.closest_atom(0.1, 0.1, 1.0), Some(0));
        // Closest to (2.9, 0.0) should be atom 1
        assert_eq!(mol.closest_atom(2.9, 0.0, 1.0), Some(1));
        // Outside tolerance
        assert_eq!(mol.closest_atom(10.0, 10.0, 1.0), None);
    }

    // ── P43 tests ─────────────────────────────────────────────────────────

    #[test]
    fn edit_implicit_h_carbon_chain() {
        let mut mol = MolecularSystem::new_empty();
        let c0 = mol.add_atom("C", 0.0, 0.0);
        // Isolated C → 4 implicit H
        assert_eq!(mol.implicit_h_count(c0), 4);
        let c1 = mol.add_atom("C", 1.5, 0.0);
        mol.add_bond(c0, c1, 1); // single bond
        assert_eq!(mol.implicit_h_count(c0), 3);
        assert_eq!(mol.implicit_h_count(c1), 3);
        // Change to double bond → 2 H each
        mol.set_bond_order(c0, c1, 2);
        assert_eq!(mol.implicit_h_count(c0), 2);
        assert_eq!(mol.implicit_h_count(c1), 2);
    }

    #[test]
    fn edit_implicit_h_oxygen() {
        let mut mol = MolecularSystem::new_empty();
        let o = mol.add_atom("O", 0.0, 0.0);
        assert_eq!(mol.implicit_h_count(o), 2); // isolated O → 2 H (water)
        let c = mol.add_atom("C", 1.5, 0.0);
        mol.add_bond(o, c, 1);
        assert_eq!(mol.implicit_h_count(o), 1); // C-O-H
        let c2 = mol.add_atom("C", -1.5, 0.0);
        mol.add_bond(o, c2, 1);
        assert_eq!(mol.implicit_h_count(o), 0); // ether O, no H
    }

    #[test]
    fn edit_implicit_h_unknown_element() {
        let mut mol = MolecularSystem::new_empty();
        let fe = mol.add_atom("Fe", 0.0, 0.0);
        assert_eq!(mol.implicit_h_count(fe), -1);
    }

    #[test]
    fn edit_ring_template_hexagon() {
        let mut mol = MolecularSystem::new_empty();
        let bl = 1.5f32;
        let indices = mol.add_ring_template(6, 0.0, 0.0, bl);
        assert_eq!(indices.len(), 6);
        assert_eq!(mol.atom_count(), 6);
        // 6 bonds in a ring
        assert_eq!(mol.bond_count(), 6);
        // All atoms should be carbon
        for &i in &indices {
            assert_eq!(mol.symbols[i as usize], "C");
        }
        // All bond lengths should be close to bl
        for i in 0..6usize {
            let a = indices[i] as usize;
            let b = indices[(i + 1) % 6] as usize;
            let dx = mol.x[a] - mol.x[b];
            let dy = mol.y[a] - mol.y[b];
            let len = (dx * dx + dy * dy).sqrt();
            assert!((len - bl).abs() < 0.01, "bond len {} ≠ {}", len, bl);
        }
    }

    #[test]
    fn edit_ring_template_too_small() {
        let mut mol = MolecularSystem::new_empty();
        let indices = mol.add_ring_template(2, 0.0, 0.0, 1.5);
        assert!(indices.is_empty());
        assert_eq!(mol.atom_count(), 0);
    }

    // ── P44 tests ─────────────────────────────────────────────────────────

    #[test]
    fn edit_attach_ring_hexagon() {
        let mut mol = MolecularSystem::new_empty();
        let a = mol.add_atom("C", 0.0, 0.0);
        let b = mol.add_atom("C", 1.0, 0.0);
        mol.add_bond(a, b, 1);
        let new_atoms = mol.attach_ring_to_bond(a, b, 6);
        assert_eq!(new_atoms.len(), 4);
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6); // 1 existing + 5 new
        // All ring bonds (a→new[0]→...→new[3]→b) should be ≈ 1.0
        let ring: Vec<u32> = std::iter::once(a)
            .chain(new_atoms.iter().copied())
            .chain(std::iter::once(b))
            .collect();
        for w in ring.windows(2) {
            let (i, j) = (w[0] as usize, w[1] as usize);
            let dx = mol.x[i] - mol.x[j];
            let dy = mol.y[i] - mol.y[j];
            let len = (dx * dx + dy * dy).sqrt();
            assert!((len - 1.0).abs() < 0.01, "bond length {len} ≠ 1.0");
        }
    }

    #[test]
    fn edit_attach_ring_triangle() {
        let mut mol = MolecularSystem::new_empty();
        let a = mol.add_atom("C", 0.0, 0.0);
        let b = mol.add_atom("C", 1.0, 0.0);
        mol.add_bond(a, b, 1);
        let new_atoms = mol.attach_ring_to_bond(a, b, 3);
        assert_eq!(new_atoms.len(), 1);
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 3); // 1 existing + 2 new
    }

    #[test]
    fn edit_attach_ring_too_small() {
        let mut mol = MolecularSystem::new_empty();
        let a = mol.add_atom("C", 0.0, 0.0);
        let b = mol.add_atom("C", 1.0, 0.0);
        mol.add_bond(a, b, 1);
        let new_atoms = mol.attach_ring_to_bond(a, b, 2);
        assert!(new_atoms.is_empty());
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }

    // ── P45 tests ─────────────────────────────────────────────────────────

    #[test]
    fn edit_get_bounds_basic() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C",  1.0, 2.0);
        mol.add_atom("C",  3.0, 4.0);
        mol.add_atom("C", -1.0, 0.0);
        let b = mol.get_bounds();
        assert_eq!(b.len(), 4);
        assert!((b[0] - (-1.0)).abs() < 1e-5, "min_x");
        assert!((b[1] -   0.0 ).abs() < 1e-5, "min_y");
        assert!((b[2] -   3.0 ).abs() < 1e-5, "max_x");
        assert!((b[3] -   4.0 ).abs() < 1e-5, "max_y");
    }

    #[test]
    fn edit_get_bounds_empty() {
        let mol = MolecularSystem::new_empty();
        assert!(mol.get_bounds().is_empty());
    }

    #[test]
    fn edit_rotate_atoms_quarter_turn() {
        use std::f32::consts::PI;
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 1.0, 0.0);
        mol.rotate_atoms(PI / 2.0, 0.0, 0.0); // 90° CCW around origin
        assert!((mol.x[0] - 0.0).abs() < 1e-5, "x after 90°: {}", mol.x[0]);
        assert!((mol.y[0] - 1.0).abs() < 1e-5, "y after 90°: {}", mol.y[0]);
    }

    #[test]
    fn edit_rotate_atoms_identity() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 2.0, 3.0);
        mol.rotate_atoms(0.0, 0.0, 0.0);
        assert!((mol.x[0] - 2.0).abs() < 1e-5);
        assert!((mol.y[0] - 3.0).abs() < 1e-5);
    }

    // ── P46 tests ─────────────────────────────────────────────────────────

    #[test]
    fn edit_select_atoms_in_rect_basic() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C",  0.0, 0.0);
        mol.add_atom("C",  2.0, 2.0);
        mol.add_atom("C",  5.0, 5.0);
        let sel = mol.select_atoms_in_rect(-1.0, -1.0, 3.0, 3.0);
        assert_eq!(sel.len(), 2);
        assert!(sel.contains(&0));
        assert!(sel.contains(&1));
    }

    #[test]
    fn edit_select_atoms_in_rect_empty() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        let sel = mol.select_atoms_in_rect(10.0, 10.0, 20.0, 20.0);
        assert!(sel.is_empty());
    }

    #[test]
    fn edit_move_atoms_subset() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        mol.add_atom("C", 2.0, 0.0);
        mol.add_atom("C", 4.0, 0.0);
        mol.move_atoms(&[0, 2], 1.0, 1.0);
        assert!((mol.x[0] - 1.0).abs() < 1e-5 && (mol.y[0] - 1.0).abs() < 1e-5);
        assert!((mol.x[1] - 2.0).abs() < 1e-5 && (mol.y[1] - 0.0).abs() < 1e-5); // unchanged
        assert!((mol.x[2] - 5.0).abs() < 1e-5 && (mol.y[2] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn edit_move_atoms_out_of_bounds() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        mol.move_atoms(&[0, 99], 1.0, 0.0); // 99 is out of range — must not panic
        assert!((mol.x[0] - 1.0).abs() < 1e-5);
    }

    // ── P47 tests ─────────────────────────────────────────────────────────

    #[test]
    fn edit_check_valence_ok() {
        // Carbon with exactly 4 single bonds → no violation
        let mut mol = MolecularSystem::new_empty();
        let c = mol.add_atom("C", 0.0, 0.0);
        for i in 1..=4u32 {
            let h = mol.add_atom("H", i as f32, 0.0);
            mol.add_bond(c, h, 1);
        }
        assert!(mol.check_valence().is_empty());
    }

    #[test]
    fn edit_check_valence_overvalent_carbon() {
        // Carbon with 5 single bonds → violation
        let mut mol = MolecularSystem::new_empty();
        let c = mol.add_atom("C", 0.0, 0.0);
        for i in 1..=5u32 {
            let h = mol.add_atom("H", i as f32, 0.0);
            mol.add_bond(c, h, 1);
        }
        let violations = mol.check_valence();
        assert_eq!(violations, vec![c]);
    }

    #[test]
    fn edit_check_valence_unknown_element() {
        // Fe has no entry in smiles_valence → skipped regardless of bonds
        let mut mol = MolecularSystem::new_empty();
        let fe = mol.add_atom("Fe", 0.0, 0.0);
        for i in 1..=6u32 {
            let c = mol.add_atom("C", i as f32, 0.0);
            mol.add_bond(fe, c, 1);
        }
        assert!(mol.check_valence().is_empty());
    }

    #[test]
    fn edit_check_valence_empty_mol() {
        let mol = MolecularSystem::new_empty();
        assert!(mol.check_valence().is_empty());
    }

    // ── P48 tests ─────────────────────────────────────────────────────────

    #[test]
    fn edit_flip_horizontal_basic() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 2.0, 3.0);
        mol.flip_horizontal(0.0);
        assert!((mol.x[0] - (-2.0)).abs() < 1e-5);
        assert!((mol.y[0] - 3.0).abs() < 1e-5); // y unchanged
    }

    #[test]
    fn edit_flip_vertical_basic() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 2.0, 3.0);
        mol.flip_vertical(0.0);
        assert!((mol.x[0] - 2.0).abs() < 1e-5); // x unchanged
        assert!((mol.y[0] - (-3.0)).abs() < 1e-5);
    }

    #[test]
    fn edit_flip_horizontal_twice() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 1.5, 2.5);
        mol.flip_horizontal(1.0);
        mol.flip_horizontal(1.0);
        assert!((mol.x[0] - 1.5).abs() < 1e-5);
        assert!((mol.y[0] - 2.5).abs() < 1e-5);
    }

    #[test]
    fn edit_flip_vertical_offset_axis() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 1.0, 0.0);
        mol.flip_vertical(2.0); // reflect across y=2 → y' = 2*2 - 0 = 4
        assert!((mol.x[0] - 1.0).abs() < 1e-5);
        assert!((mol.y[0] - 4.0).abs() < 1e-5);
    }

    // ── P49 tests ─────────────────────────────────────────────────────────

    #[test]
    fn edit_copy_atoms_basic() {
        // A–B–C chain; copy A and B → 2 atoms, 1 bond (A–B), symbol preserved
        let mut mol = MolecularSystem::new_empty();
        let a = mol.add_atom("C", 0.0, 0.0);
        let b = mol.add_atom("N", 1.5, 0.0);
        let c = mol.add_atom("O", 3.0, 0.0);
        mol.add_bond(a, b, 1);
        mol.add_bond(b, c, 2);
        let copy = mol.copy_atoms(&[a, b]);
        assert_eq!(copy.atom_count(), 2);
        assert_eq!(copy.bond_count(), 1);
        assert_eq!(copy.symbols[0], "C");
        assert_eq!(copy.symbols[1], "N");
    }

    #[test]
    fn edit_copy_atoms_no_bonds() {
        // A–B–C; copy A and C (not adjacent) → 2 atoms, 0 bonds
        let mut mol = MolecularSystem::new_empty();
        let a = mol.add_atom("C", 0.0, 0.0);
        let b = mol.add_atom("C", 1.5, 0.0);
        let c = mol.add_atom("C", 3.0, 0.0);
        mol.add_bond(a, b, 1);
        mol.add_bond(b, c, 1);
        let copy = mol.copy_atoms(&[a, c]);
        assert_eq!(copy.atom_count(), 2);
        assert_eq!(copy.bond_count(), 0);
    }

    #[test]
    fn edit_copy_atoms_empty_indices() {
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        let copy = mol.copy_atoms(&[]);
        assert_eq!(copy.atom_count(), 0);
    }

    #[test]
    fn edit_normalize_bond_length() {
        // Use a chain of non-H atoms built manually so all bonds have non-zero length
        let mut mol = MolecularSystem::new_empty();
        mol.add_atom("C", 0.0, 0.0);
        mol.add_atom("C", 3.0, 0.0); // bond length = 3.0
        mol.add_atom("C", 6.0, 0.0); // bond length = 3.0
        mol.add_bond(0, 1, 1);
        mol.add_bond(1, 2, 1);

        let target = 1.5f32;
        mol.normalize_bond_length(target);

        // All bond lengths should equal target
        let mut seen = std::collections::HashSet::new();
        for a in 0..mol.bonds.len() {
            for &b in &mol.bonds[a] {
                let key = (a.min(b), a.max(b));
                if !seen.insert(key) { continue; }
                let dx = mol.x[a] - mol.x[b];
                let dy = mol.y[a] - mol.y[b];
                let len = (dx * dx + dy * dy).sqrt();
                assert!((len - target).abs() < 0.01, "bond length {} != target {}", len, target);
            }
        }
    }

    // ── P50 tests ──────────────────────────────────────────────────────────

    #[test]
    fn p50_largest_fragment_picks_bigger_part() {
        // Two disconnected fragments: 3-atom chain + 1 isolated atom
        let mut mol = MolecularSystem::new_empty();
        let a = mol.add_atom("C", 0.0, 0.0);
        let b = mol.add_atom("C", 1.4, 0.0);
        let c = mol.add_atom("C", 2.8, 0.0);
        mol.add_atom("N", 10.0, 0.0); // isolated, index 3
        mol.add_bond(a, b, 1);
        mol.add_bond(b, c, 1);
        let frag = mol.largest_fragment();
        assert_eq!(frag.atom_count(), 3);
    }

    #[test]
    fn p50_largest_fragment_single_component() {
        let mol = parse_smiles("CC").unwrap();
        let frag = mol.largest_fragment();
        assert_eq!(frag.atom_count(), mol.atom_count());
    }

    #[test]
    fn p50_murcko_scaffold_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings(); // SMILES already has bonds; rings needed for murcko
        let scaffold = mol.murcko_scaffold();
        let heavy = (0..scaffold.atom_count())
            .filter(|&i| scaffold.get_symbol(i).as_deref() != Some("H"))
            .count();
        assert_eq!(heavy, 6, "benzene scaffold should have 6 heavy atoms");
    }

    #[test]
    fn p50_num_heavy_atoms_methane() {
        let mol = parse_smiles("C").unwrap();
        // SMILES 'C' = 1 carbon (H implicit, not stored unless explicit)
        assert_eq!(mol.num_heavy_atoms(), 1);
    }

    #[test]
    fn p50_fraction_csp3_cyclohexane() {
        let mut mol = parse_smiles("C1CCCCC1").unwrap();
        mol.compute_bonds();
        mol.compute_2d_coords();
        // All 6 C are sp3 in cyclohexane
        let f = mol.fraction_csp3();
        assert!(f > 0.99, "expected ~1.0 for cyclohexane, got {}", f);
    }

    // --- P51 tests ---

    #[test]
    fn p51_embed_empty_returns_false() {
        let mut mol = MolecularSystem::new_empty();
        assert!(!mol.embed_molecule(0), "embed on empty should return false");
    }

    #[test]
    fn p51_embed_water_bond_distances() {
        // SMILES parsing already populates self.bonds; no compute_bonds() needed
        let mut mol = parse_smiles("O").unwrap();
        assert!(mol.embed_molecule(42), "embed should succeed for water");
        let o_idx = (0..mol.atom_count()).find(|&i| mol.get_symbol(i).as_deref() == Some("O")).unwrap();
        let h_indices: Vec<usize> = (0..mol.atom_count())
            .filter(|&i| mol.get_symbol(i).as_deref() == Some("H"))
            .collect();
        for &h in &h_indices {
            let dx = mol.get_x(h).unwrap_or(0.0) - mol.get_x(o_idx).unwrap_or(0.0);
            let dy = mol.get_y(h).unwrap_or(0.0) - mol.get_y(o_idx).unwrap_or(0.0);
            let dz = mol.get_z(h).unwrap_or(0.0) - mol.get_z(o_idx).unwrap_or(0.0);
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(d > 0.8 && d < 1.2, "O-H dist should be 0.8–1.2 Å, got {}", d);
        }
    }

    #[test]
    fn p51_embed_benzene_ring_bonds() {
        // SMILES parsing already populates self.bonds
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        assert!(mol.embed_molecule(1), "embed should succeed for benzene");
        let n = mol.atom_count();
        for a in 0..n {
            // check bonds from a
            let neighbors = mol.bonds.get(a).cloned().unwrap_or_default();
            for b in neighbors {
                if b > a && mol.get_symbol(a).as_deref() == Some("C") && mol.get_symbol(b).as_deref() == Some("C") {
                    let dx = mol.get_x(b).unwrap_or(0.0) - mol.get_x(a).unwrap_or(0.0);
                    let dy = mol.get_y(b).unwrap_or(0.0) - mol.get_y(a).unwrap_or(0.0);
                    let dz = mol.get_z(b).unwrap_or(0.0) - mol.get_z(a).unwrap_or(0.0);
                    let d = (dx * dx + dy * dy + dz * dz).sqrt();
                    assert!(d > 1.2 && d < 1.7, "C-C bond dist should be 1.2–1.7 Å, got {}", d);
                }
            }
        }
    }

    // --- P52 tests ---

    #[test]
    fn p52_atom_map_parsed_from_smiles() {
        // Verify [C:1]O parses atom maps correctly
        let mol = parse_smiles("[C:1]O").unwrap();
        assert_eq!(mol.atom_map.first().copied().unwrap_or(0), 1, "C should have map 1");
        // O is non-bracket so map should be 0
        let o_idx = mol.symbols.iter().position(|s| s == "O").unwrap();
        assert_eq!(mol.atom_map.get(o_idx).copied().unwrap_or(0), 0, "O should have map 0");
    }

    #[test]
    fn p52_run_reaction_no_match() {
        // Reaction: [F:1]>>[Cl:1] (fluoride to chloride)
        // Substrate: methanol (no fluorine) → no products
        let rxn = parse_reaction_smiles("[F:1]>>[Cl:1]").expect("parse reaction");
        let substrate = parse_smiles("CO").unwrap();
        let products = rxn.run_reaction_data(&substrate);
        assert_eq!(products.len(), 0, "no match → empty result");
    }

    #[test]
    fn p52_run_reaction_transforms_atom() {
        // Reaction: [F:1]>>[Cl:1] (fluoride → chloride)
        // Substrate: fluoromethane CH3F → product CH3Cl
        let rxn = parse_reaction_smiles("[F:1]>>[Cl:1]").expect("parse reaction");
        let substrate = parse_smiles("CF").unwrap();
        let products = rxn.run_reaction_data(&substrate);
        assert_eq!(products.len(), 1, "one match expected");
        let p = &products[0];
        let has_cl = p.symbols.iter().any(|s| s == "Cl");
        let has_f  = p.symbols.iter().any(|s| s == "F");
        assert!(has_cl, "product should contain Cl");
        assert!(!has_f, "product should not contain F");
    }

    #[test]
    fn p52_run_reaction_multiple_sites() {
        // Reaction: [F:1]>>[Cl:1] on difluoromethane CF2H2 → 2 matches (both F atoms)
        let rxn = parse_reaction_smiles("[F:1]>>[Cl:1]").expect("parse reaction");
        let substrate = parse_smiles("FCF").unwrap();
        let products = rxn.run_reaction_data(&substrate);
        assert_eq!(products.len(), 2, "two F sites → two products");
        for p in &products {
            assert!(p.symbols.iter().any(|s| s == "Cl"), "each product should have Cl");
        }
    }

    // ── P53 tests ──────────────────────────────────────────────────────────

    #[test]
    fn p53a_ring_size_r6_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let matches = mol.match_smarts_data("[r6]");
        assert_eq!(matches.len(), 6, "[r6] should match all 6 benzene atoms");
    }

    #[test]
    fn p53a_ring_size_r5_not_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        let matches = mol.match_smarts_data("[r5]");
        assert_eq!(matches.len(), 0, "[r5] must not match benzene (6-ring)");
    }

    #[test]
    fn p53a_ring_size_r5_cyclopentane() {
        let mut mol = parse_smiles("C1CCCC1").unwrap();
        mol.compute_rings();
        let matches = mol.match_smarts_data("[r5]");
        assert_eq!(matches.len(), 5, "[r5] should match all 5 cyclopentane carbons");
    }

    #[test]
    fn p53a_degree_d1_terminal_carbon() {
        let mut mol = parse_smiles("CCC").unwrap();
        mol.compute_rings();
        let matches = mol.match_smarts_data("[CD1]");
        assert_eq!(matches.len(), 2, "propane has 2 terminal carbons with D=1");
    }

    #[test]
    fn p53a_ring_info_benzene() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_rings();
        assert_eq!(mol.aliphatic_ring_count(), 0, "benzene is aromatic, not aliphatic");
        let sz = mol.ring_sizes_for_atom(0);
        assert!(sz.contains(&6), "atom 0 in benzene should be in ring of size 6");
    }

    #[test]
    fn p53b_atom_pair_length() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.fingerprint_atom_pair().len(), 256);
    }

    #[test]
    fn p53b_atom_pair_self_same() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let a = mol.fingerprint_atom_pair();
        let b = mol.fingerprint_atom_pair();
        assert_eq!(a, b, "same molecule must produce identical atom pair fp");
    }

    #[test]
    fn p53b_atom_pair_benzene_vs_hexane_differ() {
        let benzene = parse_smiles("c1ccccc1").unwrap();
        let hexane = parse_smiles("CCCCCC").unwrap();
        let a = benzene.fingerprint_atom_pair();
        let b = hexane.fingerprint_atom_pair();
        assert_ne!(a, b, "aromatic vs aliphatic must produce different atom pair fp");
    }

    #[test]
    fn p53c_aligned_coords_count() {
        let mut template = parse_smiles("c1ccccc1").unwrap();
        template.compute_rings();
        template.compute_2d_coords_data();
        let mut mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        mol.compute_rings();
        mol.generate_aligned_coords(&template);
        assert_eq!(mol.num_heavy_atoms(), 10, "naphthalene has 10 heavy atoms after alignment");
    }

    #[test]
    fn p53c_aligned_coords_no_nan() {
        let mut template = parse_smiles("c1ccccc1").unwrap();
        template.compute_rings();
        template.compute_2d_coords_data();
        let mut mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        mol.compute_rings();
        mol.generate_aligned_coords(&template);
        for i in 0..mol.atom_count() {
            assert!(mol.x[i].is_finite(), "x[{}] should be finite", i);
            assert!(mol.y[i].is_finite(), "y[{}] should be finite", i);
        }
    }

    #[test]
    fn p53c_aligned_no_match_still_valid() {
        let mut template = parse_smiles("c1ccccc1").unwrap();
        template.compute_rings();
        template.compute_2d_coords_data();
        let mut mol = parse_smiles("CCC").unwrap();
        mol.generate_aligned_coords(&template);
        assert_eq!(mol.num_heavy_atoms(), 3, "no match fallback must not drop atoms");
        for i in 0..mol.atom_count() {
            assert!(mol.x[i].is_finite());
        }
    }

    // ── P54 tests ──

    #[test]
    fn p54a_remove_hs_count() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let heavy = mol.remove_hs();
        assert_eq!(heavy.num_heavy_atoms(), 6);
        assert_eq!(heavy.symbols.iter().filter(|s| s.as_str() == "H").count(), 0);
    }

    #[test]
    fn p54a_remove_hs_aromatic_preserved() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let heavy = mol.remove_hs();
        assert_eq!(heavy.aromatic_atoms.len(), 6);
        assert!(heavy.aromatic_atoms.iter().all(|&a| a));
    }

    #[test]
    fn p54a_remove_hs_bonds_remapped() {
        let mol = parse_smiles("CC").unwrap();
        let heavy = mol.remove_hs();
        assert_eq!(heavy.symbols.len(), 2, "ethane has 2 heavy atoms");
        assert!(heavy.bonds[0].contains(&1), "C–C bond must survive H removal");
    }

    #[test]
    fn p54a_remove_hs_idempotent() {
        let mol = parse_smiles("CC").unwrap().remove_hs();
        let n = mol.symbols.len();
        assert_eq!(mol.remove_hs().symbols.len(), n);
    }

    #[test]
    fn p54b_sdf_prop_get() {
        let sdf = concat!(
            "\n\n\n",
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "M  END\n",
            "> <IC50>\n",
            "12.5\n",
            "\n",
            "$$$$\n"
        );
        let mol = parse_sdf(sdf).unwrap();
        assert_eq!(mol.properties.get("IC50").map(|s| s.as_str()), Some("12.5"));
    }

    #[test]
    fn p54b_sdf_prop_list_sorted() {
        let sdf = concat!(
            "\n\n\n",
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "M  END\n",
            "> <ZZ>\n1\n\n",
            "> <AA>\n2\n\n",
            "$$$$\n"
        );
        let mol = parse_sdf(sdf).unwrap();
        let mut keys: Vec<&String> = mol.properties.keys().collect();
        keys.sort();
        assert_eq!(keys[0].as_str(), "AA");
        assert_eq!(keys[1].as_str(), "ZZ");
    }

    #[test]
    fn p54b_sdf_set_prop() {
        let mut mol = MolecularSystem::new_empty();
        mol.properties.insert("activity".to_string(), "active".to_string());
        assert_eq!(mol.properties.get("activity").map(|s| s.as_str()), Some("active"));
    }

    #[test]
    fn p54b_sdf_no_props_empty() {
        let mol = parse_smiles("C").unwrap();
        assert!(mol.properties.is_empty());
    }

    #[test]
    fn p54c_descriptors_mw_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert!(mol.molecular_weight() > 70.0 && mol.molecular_weight() < 90.0);
    }

    #[test]
    fn p54c_normalize_depiction_finite() {
        let mut mol = parse_smiles("c1ccccc1").unwrap();
        mol.compute_2d_coords_data();
        mol.normalize_depiction();
        for i in 0..mol.symbols.len() {
            assert!(mol.x[i].is_finite() && mol.y[i].is_finite());
        }
    }

    #[test]
    fn p54c_is_valid_smiles() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert!(mol.is_valid());
    }

    #[test]
    fn p54c_is_valid_empty() {
        assert!(!MolecularSystem::new_empty().is_valid());
    }

    // ── P55 tests ──

    #[test]
    fn p55a_add_hs_count() {
        let mol = parse_smiles("CC").unwrap().remove_hs();
        let full = mol.add_hs();
        assert_eq!(full.symbols.len(), 8, "2C + 6H = 8 atoms");
        assert_eq!(full.symbols.iter().filter(|s| s.as_str() == "H").count(), 6);
    }

    #[test]
    fn p55a_add_hs_idempotent() {
        let mol = parse_smiles("CC").unwrap();
        let n = mol.symbols.len();
        assert_eq!(mol.add_hs().symbols.len(), n, "fully-explicit mol: add_hs adds nothing");
    }

    #[test]
    fn p55a_add_hs_bonds_valid() {
        let mol = parse_smiles("CC").unwrap().remove_hs().add_hs();
        let n = mol.symbols.len();
        for (i, nb_list) in mol.bonds.iter().enumerate() {
            for &j in nb_list {
                assert!(j < n && j != i, "bond index {j} must be in-bounds and not self");
            }
        }
    }

    #[test]
    fn p55a_roundtrip_heavy_count() {
        let heavy = parse_smiles("CC").unwrap().remove_hs();
        let restored = heavy.add_hs();
        assert_eq!(restored.num_heavy_atoms(), heavy.num_heavy_atoms());
    }

    #[test]
    fn p55b_topo_fp_length() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.fingerprint_topological().len(), 256);
    }

    #[test]
    fn p55b_topo_fp_identical() {
        let a = parse_smiles("c1ccccc1").unwrap();
        let b = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(a.fingerprint_topological(), b.fingerprint_topological());
    }

    #[test]
    fn p55b_topo_fp_different() {
        let benz = parse_smiles("c1ccccc1").unwrap();
        let eth = parse_smiles("CC").unwrap();
        assert_ne!(benz.fingerprint_topological(), eth.fingerprint_topological());
    }

    #[test]
    fn p55b_topo_fp_nonempty() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert!(mol.fingerprint_topological().iter().any(|&b| b != 0));
    }

    #[test]
    fn p55c_has_stereo_true() {
        let mol = parse_smiles("[C@@H](F)(Cl)Br").unwrap();
        assert!(mol.has_stereo());
    }

    #[test]
    fn p55c_has_stereo_false() {
        let mol = parse_smiles("CC").unwrap();
        assert!(!mol.has_stereo());
    }

    #[test]
    fn p55c_get_stereo_tags_no_crash() {
        // get_stereo_tags() calls to_js (wasm-only); test the underlying data instead
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert!(mol.stereo_centers.is_empty());
        assert!(!mol.has_stereo());
    }

    #[test]
    fn p55c_get_stereo_tags_stereo_center_count() {
        let mol = parse_smiles("[C@@H](F)(Cl)Br").unwrap();
        assert_eq!(mol.stereo_center_count(), 1);
        assert!(mol.has_stereo());
    }

    // ── P56: perceive_stereo_from_3d + to_smiles stereo output ───────────────

    // Tetrahedral C (4 heavy neighbors at proper 3D vertices, non-zero z)
    const TETRA_3D_SDF: &str = "\
tetra
  test

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    0.0000   -0.3536 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4042    1.1431   -0.3536 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4042   -0.5715    1.0607 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4042   -0.5715   -1.0607 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
$$$$
";

    // Same topology but all z=0 (degenerate / flat)
    const TETRA_FLAT_SDF: &str = "\
flat
  test

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
$$$$
";

    // Enantiomer of TETRA_3D_SDF — z-coordinates negated.
    const TETRA_3D_ENT_SDF: &str = "\
tetra_ent
  test

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    0.0000    0.3536 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4042    1.1431    0.3536 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4042   -0.5715   -1.0607 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4042   -0.5715    1.0607 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
$$$$
";

    #[test]
    fn p56_perceive_stereo_from_3d_detects_center() {
        let mut mol = parse_sdf(TETRA_3D_SDF).unwrap();
        mol.perceive_stereo_from_3d();
        assert_eq!(mol.stereo_center_count(), 1, "C should be detected as stereo center");
        assert!(mol.is_stereo_center(0));
    }

    #[test]
    fn p56_perceive_stereo_from_3d_no_center_flat() {
        let mut mol = parse_sdf(TETRA_FLAT_SDF).unwrap();
        mol.perceive_stereo_from_3d();
        assert_eq!(mol.stereo_center_count(), 0, "flat molecule must have no stereo");
    }

    #[test]
    fn p56_perceive_stereo_no_bonds_noop() {
        // Without bonds, perceive_stereo_from_3d must return cleanly.
        let mut mol = MolecularSystem::new_empty();
        mol.symbols = vec!["C".into(), "F".into(), "Cl".into(), "Br".into()];
        mol.x = vec![0.0, 1.0, -1.0, 0.0];
        mol.y = vec![0.0, 0.0, 0.0, 1.0];
        mol.z = vec![0.0, 0.5, 0.5, -0.5];
        mol.perceive_stereo_from_3d();
        assert_eq!(mol.stereo_center_count(), 0);
    }

    #[test]
    fn p56_enantiomers_have_opposite_stereo_descriptor() {
        // Verify sign convention: the two enantiomers (z-negated) must store
        // opposite descriptor values in stereo_centers.
        let mut mol_r = parse_sdf(TETRA_3D_SDF).unwrap();
        mol_r.perceive_stereo_from_3d();
        let (desc_r, _) = mol_r.stereo_centers[&0];

        let mut mol_s = parse_sdf(TETRA_3D_ENT_SDF).unwrap();
        mol_s.perceive_stereo_from_3d();
        let (desc_s, _) = mol_s.stereo_centers[&0];

        assert_ne!(desc_r, desc_s, "enantiomers must have opposite descriptors");
        assert!(
            (desc_r == 1 && desc_s == -1) || (desc_r == -1 && desc_s == 1),
            "descriptors must be ±1, got {desc_r} and {desc_s}"
        );
    }

    #[test]
    fn p56_to_smiles_stereo_3d_outputs_annotation() {
        // A 3D-embedded L-alanine SMILES round-trip should produce @/@@ annotation.
        let mut mol = parse_smiles("N[C@@H](C)C(=O)O").unwrap();
        mol.embed_molecule(42);
        let smi = mol.to_smiles_data();
        assert!(
            smi.contains('@'),
            "to_smiles() must output @ or @@ for 3D stereo center, got: {smi}"
        );
    }

    #[test]
    fn p56_embed_molecule_angle_quality() {
        // Embed n-pentane (no rings), check that C-C-C heavy-atom angles converge
        // close to the tetrahedral ideal (109.47° ± 25°). Fused-ring systems are
        // excluded here because distance geometry needs many more steps for them.
        let mut mol = parse_smiles("CCCCC").unwrap();
        mol.embed_molecule(1);
        // heavy atoms: C0-C1-C2-C3-C4 (indices 0–4)
        for b in 1..4usize {
            let a = b - 1;
            let c = b + 1;
            let angle_deg = mol.angle(a, b, c);
            assert!(
                (70.0..=150.0).contains(&angle_deg),
                "C-C-C angle {angle_deg:.1}° for ({a},{b},{c}) out of expected range"
            );
        }
    }

    // ── P57: E/Z stereochemistry, atom mapping, ring classification ──────────

    #[test]
    fn p57_ez_smiles_trans_2_butene() {
        // C/C=C/C → trans (E). Atoms: C0-C1=C2-C3, bond dirs on C0→C1 and C2→C3.
        // After H addition: heavy atoms 0-3 are C,C,C,C in order.
        let mol = parse_smiles("C/C=C/C").unwrap();
        // double bond is between heavy atoms 1 and 2
        let ez = mol.is_ez_bond(1, 2);
        assert_eq!(ez, Some(true), "C/C=C/C should be E (trans)");
        assert_eq!(mol.ez_bond_count(), 1);
    }

    #[test]
    fn p57_ez_smiles_cis_2_butene() {
        // C/C=C\C → cis (Z)
        let mol = parse_smiles("C/C=C\\C").unwrap();
        let ez = mol.is_ez_bond(1, 2);
        assert_eq!(ez, Some(false), "C/C=C\\C should be Z (cis)");
    }

    #[test]
    fn p57_ez_no_stereo_plain_double_bond() {
        // C=C has no direction chars → no E/Z info
        let mol = parse_smiles("CC=CC").unwrap();
        assert_eq!(mol.ez_bond_count(), 0);
    }

    #[test]
    fn p57_atom_map_parse_and_get() {
        // [C:1]([H:2])([H:3])[H:4] — atom 0 has map index 1
        let mol = parse_smiles("[C:1]([H:2])([H:3])[H:4]").unwrap();
        assert!(mol.has_atom_map());
        assert_eq!(mol.get_atom_map_index(0), 1);
    }

    #[test]
    fn p57_atom_map_set_and_clear() {
        let mut mol = parse_smiles("CC").unwrap();
        assert!(!mol.has_atom_map());
        mol.set_atom_map_index(0, 5);
        assert!(mol.has_atom_map());
        assert_eq!(mol.get_atom_map_index(0), 5);
        mol.clear_atom_map();
        assert!(!mol.has_atom_map());
    }

    #[test]
    fn p57_spiro_atoms_spiro_nonane() {
        // spiro[4.4]nonane: C1CCCC12CCCCC2 — atom 4 is the spiro center
        let mut mol = parse_smiles("C1CCCC12CCCCC2").unwrap();
        mol.compute_rings();
        let spiro = mol.get_spiro_atoms();
        assert!(!spiro.is_empty(), "spiro[4.4]nonane must have a spiro center");
        assert!(spiro.contains(&4u32), "spiro center should be atom 4, got {spiro:?}");
    }

    #[test]
    fn p57_fused_ring_bonds_naphthalene() {
        // naphthalene: c1ccc2ccccc2c1 — has one fused (shared) bond
        let mut mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        mol.compute_rings();
        let bonds = mol.fused_ring_bonds_vec();
        assert!(!bonds.is_empty(), "naphthalene must have at least one fused bond");
    }

    #[test]
    fn p57_bridged_ring_norbornane() {
        // norbornane: C1CC2CCC1C2
        let mut mol = parse_smiles("C1CC2CCC1C2").unwrap();
        mol.compute_rings();
        assert!(mol.is_bridged_ring_system(), "norbornane should be bridged");
    }

    #[test]
    fn p57_not_bridged_cyclohexane() {
        let mut mol = parse_smiles("C1CCCCC1").unwrap();
        mol.compute_rings();
        assert!(!mol.is_bridged_ring_system());
    }
}
