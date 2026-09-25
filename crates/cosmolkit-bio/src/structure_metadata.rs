//! Source-state metadata retained from a Gemmi structural model.

use crate::BioTransform;
use std::collections::BTreeMap;

/// Source-defined structure state that is not part of the atom hierarchy.
#[derive(Debug, Clone, PartialEq)]
pub struct BioStructureSourceState {
    pub name: String,
    pub resolution: f64,
    /// PDB serial-to-serial connectivity; duplicate targets encode repeated
    /// connectivity entries and remain distinct in each sorted target list.
    pub conect_map: BTreeMap<i32, Vec<i32>>,
    pub has_d_fraction: bool,
    pub non_ascii_line: i32,
    /// The original single-byte TER status (`0` when no TER record was read).
    pub ter_status: u8,
    pub has_origx: bool,
    pub origx: BioTransform,
    /// Minimal mmCIF tag/value metadata retained by Gemmi's Structure.
    pub info: BTreeMap<String, String>,
    /// Original PDB REMARK records, in source order.
    pub raw_remarks: Vec<String>,
}

impl Default for BioStructureSourceState {
    fn default() -> Self {
        // Gemmi✔️✔️: std::string name;
        // Gemmi✔️✔️: std::map<int, std::vector<int>> conect_map;
        // Gemmi✔️✔️: bool has_d_fraction = false;  // uses Refmac's ccp4_deuterium_fraction
        // Gemmi✔️✔️: int non_ascii_line = 0;  // first PDB line with non-ASCII bytes, or 0
        // Gemmi✔️✔️: char ter_status = '\0';
        // Gemmi✔️✔️: bool has_origx = false;
        // Gemmi✔️✔️: Transform origx;
        // Gemmi✔️✔️: std::map<std::string, std::string> info;
        // Gemmi✔️✔️: std::vector<std::string> raw_remarks;
        // Gemmi✔️✔️: double resolution = 0;
        // Gemmi✔️✔️: struct Transform {
        // Gemmi✔️✔️:   Mat33 mat;
        // Gemmi✔️✔️:   Vec3 vec;
        // Gemmi✔️✔️: };
        // Gemmi✔️✔️: double a[3][3] = { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} };
        // Gemmi✔️✔️: Vec3_() : x(0), y(0), z(0) {}
        // Behavior review: source-defaulted scalars are false/zero/NUL, empty
        // strings match std::string default construction, and `origx` composes
        // the source identity Mat33 with the zero-initialized Vec3. The raw
        // TER byte is kept losslessly without Unicode interpretation.
        // Complexity review: this value initializes fixed-size scalars and
        // arrays in O(1); its empty String allocates no backing storage.
        Self {
            name: String::new(),
            resolution: 0.0,
            conect_map: BTreeMap::new(),
            has_d_fraction: false,
            non_ascii_line: 0,
            ter_status: 0,
            has_origx: false,
            origx: BioTransform::identity(),
            info: BTreeMap::new(),
            raw_remarks: Vec::new(),
        }
    }
}
