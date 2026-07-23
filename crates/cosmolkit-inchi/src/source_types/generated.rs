// Generated from IUPAC InChI v1.07.5 at commit
// 11a87982bb518f57ac013f0b258c283655e1ea1d.
// Input: dev/inchi_type_inventory_wrapper.h plus all 60 production C files.
// Configuration: COMPILE_ANSI_ONLY, TARGET_API_LIB, GCC/Clang LP64.
// Generator declaration baseline: bindgen-cli 0.72.1.
// The checked-in result is pure Rust and contains no FFI declarations.
use super::*;

pub const CURRENT_VER: &[u8; 7] = b"1.07.5\0";
pub const RINCHI_TEST: u32 = 0;
pub const SPRINTF_FLAG: u32 = 2;
pub const APP_DESCRIPTION: &[u8; 47] = b"InChI version 1, Software 1.07.5 (API Library)\0";
pub const TARGET_PLATFORM: &[u8; 6] = b"Linux\0";
pub const BUILD_WITH_ENG_OPTIONS: u32 = 0;
pub const BUILD_WITH_AMI: u32 = 1;
pub const ADD_CMLPP: u32 = 0;
pub const INCHI_VERSION: &[u8; 2] = b"1\0";
pub const INCHI_NAME: &[u8; 6] = b"InChI\0";
pub const INCHI_NAM_VER_DELIM: &[u8; 2] = b"=\0";
pub const INCHI_OPTION_PREFX: u8 = 45u8;
pub const INCHI_PATH_DELIM: u8 = 47u8;
pub const INCHI_ALT_OPT_PREFIX: u8 = 45u8;
pub const INCHI_ACD_LABS_PREFIX: u8 = 45u8;
pub const bRELEASE_VERSION: u32 = 1;
pub const RELEASE_IS_FINAL: u32 = 1;
pub const DISPLAY_DEBUG_DATA_C_POINT: u32 = 0;
pub const DISPLAY_ORIG_AT_NUMBERS: u32 = 1;
pub const DISPLAY_ZZ_AS_STAR: u32 = 1;
pub const FIX_ChCh_STEREO_CANON_BUG: u32 = 1;
pub const ADD_ChCh_STEREO_CANON_CHK: u32 = 0;
pub const FIX_ChCh_CONSTIT_CANON_BUG: u32 = 1;
pub const FIX_EITHER_STEREO_IN_AUX_INFO: u32 = 1;
pub const FIX_NORM_BUG_ADD_ION_PAIR: u32 = 1;
pub const FIX_REM_PROTON_COUNT_BUG: u32 = 1;
pub const FIX_READ_AUX_MEM_LEAK: u32 = 1;
pub const FIX_READ_LONG_LINE_BUG: u32 = 1;
pub const FIX_N_V_METAL_BONDS_GPF: u32 = 1;
pub const BNS_RAD_SEARCH: u32 = 1;
pub const FIX_ODD_THINGS_REM_Plus_BUG: u32 = 0;
pub const FIX_N_MINUS_NORN_BUG: u32 = 0;
pub const FIX_CANCEL_CHARGE_COUNT_BUG: u32 = 0;
pub const FIX_2D_STEREO_BORDER_CASE: u32 = 0;
pub const FIX_REM_ION_PAIRS_Si_BUG: u32 = 0;
pub const FIX_STEREO_SCALING_BUG: u32 = 0;
pub const FIX_EMPTY_LAYER_BUG: u32 = 0;
pub const FIX_EITHER_DB_AS_NONSTEREO: u32 = 0;
pub const FIX_BOND23_IN_TAUT: u32 = 0;
pub const FIX_TACN_POSSIBLE_BUG: u32 = 0;
pub const FIX_KEEP_H_ON_NH_ANION: u32 = 0;
pub const FIX_AVOID_ADP: u32 = 0;
pub const FIX_NUM_TG: u32 = 0;
pub const FIX_CPOINT_BOND_CAP2: u32 = 0;
pub const FIX_ISO_FIXEDH_BUG: u32 = 1;
pub const FIX_ISO_FIXEDH_BUG_READ: u32 = 0;
pub const FIX_DALKE_BUGS: u32 = 1;
pub const FIX_TRANSPOSITION_CHARGE_BUG: u32 = 1;
pub const FIX_I2I_STEREOCONVERSION_BUG: u32 = 1;
pub const FIX_I2I_STEREOCONVERSION_BUG2: u32 = 1;
pub const FIX_I2I_STEREOCONVERSION_BUG3: u32 = 1;
pub const FIX_TERM_H_CHRG_BUG: u32 = 1;
pub const FIX_AROM_RADICAL: u32 = 1;
pub const FIX_STEREOCOUNT_ERR: u32 = 1;
pub const FIX_IMPOSSIBLE_H_ISOTOPE_BUG: u32 = 1;
pub const FIX_OLEAN_SPIRO_CHIRALITY_DETECTION_BUG: u32 = 1;
pub const FIX_RENUM_BUG_FOR_CASE_OF_ACIDIC_OH_AT_P_PLUS: u32 = 1;
pub const CHECK_STRTOL_ATNUMB: u32 = 1;
pub const DISABLE_READ_COMPRESSED_INCHI: u32 = 1;
pub const FIX_GAF_2019_1: u32 = 1;
pub const FIX_GAF_2019_2: u32 = 1;
pub const FIX_GAF_2020_GENERIC: u32 = 1;
pub const FIX_OSS_FUZZ_25604: u32 = 1;
pub const FIX_GAF_2020_25607: u32 = 1;
pub const FIX_GAF_2020_25726: u32 = 1;
pub const FIX_GAF_2020_25741: u32 = 1;
pub const FIX_OSS_FUZZ_30162_30343: u32 = 1;
pub const FIX_OSS_FUZZ_25734_28139: u32 = 1;
pub const FIX_CURE53_ISSUE_OOB_ALREADY_HAVE_THIS_MESSAGE: u32 = 1;
pub const FIX_CURE53_ISSUE_HEAP_BUFFER_OVERFLOW_INCHITOINPATOM: u32 = 1;
pub const FIX_CURE53_ISSUE_NULL_DEREFERENCE_MAKE_A_COPY_OF_T_GROUP_INFO: u32 = 1;
pub const FIX_GAF_2019_3: u32 = 1;
pub const FIX_ONE_LINE_INCHI_INPUT_CONVERSION_ISSUE: u32 = 1;
pub const ALLOW_EMPTY_INCHI_AS_INPUT: u32 = 1;
pub const I2S_MODIFY_OUTPUT: u32 = 0;
pub const FIX_NP_MINUS_BUG: u32 = 1;
pub const FIX_ADJ_RAD: u32 = 0;
pub const SDF_OUTPUT_V2000: u32 = 1;
pub const SDF_OUTPUT_DT: u32 = 1;
pub const CHECK_AROMBOND2ALT: u32 = 1;
pub const READ_INCHI_STRING: u32 = 1;
pub const INCLUDE_NORMALIZATION_ENTRY_POINT: u32 = 0;
pub const KETO_ENOL_TAUT: u32 = 1;
pub const TAUT_15_NON_RING: u32 = 1;
pub const TAUT_PT_22_00: u32 = 1;
pub const TAUT_PT_16_00: u32 = 1;
pub const TAUT_PT_06_00: u32 = 1;
pub const TAUT_PT_39_00: u32 = 1;
pub const TAUT_PT_13_00: u32 = 1;
pub const TAUT_PT_18_00: u32 = 1;
pub const UNDERIVATIZE: u32 = 1;
pub const RING2CHAIN: u32 = 1;
pub const UNDERIVATIZE_REPORT: u32 = 1;
pub const HAL_ACID_H_XCHG: u32 = 1;
pub const CANON_FIXH_TRANS: u32 = 1;
pub const STEREO_WEDGE_ONLY: u32 = 1;
pub const REMOVE_ION_PAIRS_EARLY: u32 = 1;
pub const REMOVE_ION_PAIRS_DISC_STRU: u32 = 1;
pub const REMOVE_ION_PAIRS_FIX_BONDS: u32 = 1;
pub const S_VI_O_PLUS_METAL_FIX_BOND: u32 = 1;
pub const N_V_STEREOBONDS: u32 = 1;
pub const REMOVE_ION_PAIRS_ORIG_STRU: u32 = 0;
pub const DISCONNECT_SALTS: u32 = 1;
pub const TEST_REMOVE_S_ATOMS: u32 = 1;
pub const CHARGED_SALTS_ONLY: u32 = 1;
pub const BNS_PROTECT_FROM_TAUT: u32 = 1;
pub const BNS_MARK_EDGE_2_DISCONNECT: u32 = 1;
pub const REPLACE_ALT_WITH_TAUT: u32 = 1;
pub const MOVE_CHARGES: u32 = 1;
pub const NEUTRALIZE_ENDPOINTS: u32 = 1;
pub const FIX_H_CHECKING_TAUT: u32 = 1;
pub const ALWAYS_ADD_TG_ON_THE_FLY: u32 = 1;
pub const IGNORE_SINGLE_ENDPOINTS: u32 = 1;
pub const INCL_NON_SALT_CANDIDATATES: u32 = 1;
pub const SALT_WITH_PROTONS: u32 = 1;
pub const OPPOSITE_CHARGE_IN_CGROUP: u32 = 1;
pub const MOVE_PPLUS_TO_REMOVE_PROTONS: u32 = 0;
pub const ADD_MOVEABLE_O_PLUS: u32 = 1;
pub const DISCONNECT_METALS: u32 = 1;
pub const RECONNECT_METALS: u32 = 0;
pub const CHECK_METAL_VALENCE: u32 = 0;
pub const bREUSE_INCHI: u32 = 1;
pub const OUTPUT_CONNECTED_METAL_ONLY: u32 = 0;
pub const EMBED_REC_METALS_INCHI: u32 = 1;
pub const bOUTPUT_ONE_STRUCT_TIME: u32 = 1;
pub const INCHI_NUM: u32 = 2;
pub const INCHI_BAS: u32 = 0;
pub const INCHI_REC: u32 = 1;
pub const TAUT_NUM: u32 = 2;
pub const TAUT_NON: u32 = 0;
pub const TAUT_YES: u32 = 1;
pub const TAUT_INI: u32 = 2;
pub const OUT_N1: u32 = 0;
pub const OUT_T1: u32 = 1;
pub const OUT_NT: u32 = 2;
pub const OUT_TN: u32 = 3;
pub const OUT_NN: u32 = 4;
pub const NEW_STEREOCENTER_CHECK: u32 = 1;
pub const MIN_SB_RING_SIZE: u32 = 8;
pub const REMOVE_KNOWN_NONSTEREO: u32 = 1;
pub const REMOVE_CALC_NONSTEREO: u32 = 1;
pub const PROPAGATE_ILL_DEF_STEREO: u32 = 1;
pub const ONLY_DOUBLE_BOND_STEREO: u32 = 0;
pub const ONE_BAD_SB_NEIGHBOR: u32 = 1;
pub const BREAK_ONE_MORE_SC_TIE: u32 = 1;
pub const BREAK_ALSO_NEIGH_TIE: u32 = 0;
pub const BREAK_ALSO_NEIGH_TIE_ROTATE: u32 = 1;
pub const STEREO_CENTER_BONDS_NORM: u32 = 1;
pub const STEREO_CENTER_BOND4_NORM: u32 = 0;
pub const NORMALIZE_INP_COORD: u32 = 0;
pub const CHECK_C2v_S4_SYMM: u32 = 0;
pub const EQL_H_NUM_TOGETHER: u32 = 1;
pub const ABC_CT_NUM_CLOSURES: u32 = 1;
pub const SINGLET_IS_TRIPLET: u32 = 1;
pub const FIND_CANON_NE_EQUITABLE: u32 = 0;
pub const EXTR_KNOWN_USED_TO_REMOVE_PARITY: u32 = 1;
pub const EXTR_CALC_USED_TO_REMOVE_PARITY: u32 = 2;
pub const EXTR_2EQL2CENTER_TO_REMOVE_PARITY: u32 = 4;
pub const EXTR_HAS_ATOM_WITH_DEFINED_PARITY: u32 = 8;
pub const EXTR_REMOVE_PARITY_WARNING: u32 = 16;
pub const EXTR_SALT_WAS_DISCONNECTED: u32 = 32;
pub const EXTR_SALT_PROTON_MOVED: u32 = 64;
pub const EXTR_SALT_PROTON_MOVE_ERR_WARN: u32 = 128;
pub const EXTR_METAL_WAS_DISCONNECTED: u32 = 256;
pub const EXTR_METAL_WAS_NOT_DISCONNECTED: u32 = 512;
pub const EXTR_NON_TRIVIAL_STEREO: u32 = 1024;
pub const EXTR_UNUSUAL_VALENCES: u32 = 2048;
pub const EXTR_HAS_METAL_ATOM: u32 = 4096;
pub const EXTR_TEST_TAUT3_SALTS_DONE: u32 = 8192;
pub const EXTR_CANON_NE_EQUITABLE: u32 = 16384;
pub const EXTR_HAS_PROTON_PN: u32 = 32768;
pub const EXTR_HAS_FEATURE: u32 = 65536;
pub const EXTR_TAUT_TREATMENT_CHARGES: u32 = 131072;
pub const EXTR_TRANSPOSITION_EXAMPLES: u32 = 262144;
pub const EXTR_MASK: u32 = 0;
pub const EXTR_FLAGS: u32 = 0;
pub const TAUT_TROPOLONE_7: u32 = 1;
pub const TAUT_TROPOLONE_5: u32 = 1;
pub const TAUT_4PYRIDINOL_RINGS: u32 = 1;
pub const TAUT_PYRAZOLE_RINGS: u32 = 1;
pub const TAUT_IGNORE_EQL_ENDPOINTS: u32 = 0;
pub const TAUT_RINGS_ATTACH_CHAIN: u32 = 1;
pub const FIND_RING_SYSTEMS: u32 = 1;
pub const FIND_RINS_SYSTEMS_DISTANCES: u32 = 0;
pub const USE_DISTANCES_FOR_RANKING: u32 = 0;
pub const DISPLAY_RING_SYSTEMS: u32 = 0;
pub const TAUT_OTHER: u32 = 1;
pub const APPLY_IMPLICIT_H_DOWN_RULE: u32 = 0;
pub const ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS: u32 = 1;
pub const IGNORE_TGROUP_WITHOUT_H: u32 = 1;
pub const REMOVE_TGROUP_CHARGE: u32 = 0;
pub const INCHI_T_NUM_MOVABLE: u32 = 2;
pub const USE_AUX_RANKING: u32 = 1;
pub const USE_AUX_RANKING_ALL: u32 = 1;
pub const USE_ISO_SORT_KEY_HFIXED: u32 = 0;
pub const REL_RAC_STEREO_IGN_1_SC: u32 = 0;
pub const CMODE_CT: u32 = 1;
pub const CMODE_ISO: u32 = 2;
pub const CMODE_ISO_OUT: u32 = 4;
pub const CMODE_STEREO: u32 = 8;
pub const CMODE_ISO_STEREO: u32 = 16;
pub const CMODE_TAUT: u32 = 32;
pub const CMODE_NOEQ_STEREO: u32 = 64;
pub const CMODE_REDNDNT_STEREO: u32 = 128;
pub const CMODE_NO_ALT_SBONDS: u32 = 256;
pub const CMODE_RELATIVE_STEREO: u32 = 512;
pub const CMODE_RACEMIC_STEREO: u32 = 1024;
pub const CMODE_SC_IGN_ALL_UU: u32 = 2048;
pub const CMODE_SB_IGN_ALL_UU: u32 = 4096;
pub const CANON_MODE_CT: u32 = 1;
pub const CANON_MODE_TAUT: u32 = 33;
pub const CANON_MODE_ISO: u32 = 7;
pub const CANON_MODE_STEREO: u32 = 9;
pub const CANON_MODE_ISO_STEREO: u32 = 23;
pub const CANON_MODE_MASK: u32 = 255;
pub const CT_ATOMID_DONTINCLUDE: u32 = 1;
pub const CT_ATOMID_IS_INITRANK: u32 = 2;
pub const CT_ATOMID_IS_CURRANK: u32 = 3;
pub const CANON_TAUTOMERS: u32 = 1;
pub const HYDROGENS_IN_INIT_RANKS: u32 = 1;
pub const DOUBLE_BOND_NEIGH_LIST: u32 = 0;
pub const INCL_NON_6AROM: u32 = 1;
pub const CT_ATOMID: u32 = 3;
pub const USE_SYMMETRY_TO_ACCELERATE: u32 = 1;
pub const CT_INITVALUE: i32 = -1;
pub const BEST_PARITY: u32 = 1;
pub const WORSE_PARITY: u32 = 2;
pub const ALL_ALT_AS_AROMATIC: u32 = 1;
pub const ANY_ATOM_IN_ALT_CYCLE: u32 = 1;
pub const EXCL_ALL_AROM_BOND_PARITY: u32 = 0;
pub const ADD_6MEMB_AROM_BOND_PARITY: u32 = 1;
pub const MAX_NUM_STEREO_BONDS: u32 = 3;
pub const MAX_NUM_STEREO_BOND_NEIGH: u32 = 3;
pub const MIN_NUM_STEREO_BOND_NEIGH: u32 = 2;
pub const MAX_NUM_STEREO_ATOM_NEIGH: u32 = 4;
pub const STEREO_AT_MARK: u32 = 8;
pub const DRAW_AROM_TAUT: u32 = 1;
pub const TG_FLAG_TEST_TAUT__ATOMS: u32 = 1;
pub const TG_FLAG_DISCONNECT_SALTS: u32 = 2;
pub const TG_FLAG_TEST_TAUT__SALTS: u32 = 4;
pub const TG_FLAG_MOVE_POS_CHARGES: u32 = 8;
pub const TG_FLAG_TEST_TAUT2_SALTS: u32 = 16;
pub const TG_FLAG_ALLOW_NO_NEGTV_O: u32 = 32;
pub const TG_FLAG_MERGE_TAUT_SALTS: u32 = 64;
pub const TG_FLAG_ALL_TAUTOMERIC: u32 = 85;
pub const TG_FLAG_DISCONNECT_COORD: u32 = 128;
pub const TG_FLAG_RECONNECT_COORD: u32 = 256;
pub const TG_FLAG_CHECK_VALENCE_COORD: u32 = 512;
pub const TG_FLAG_MOVE_HPLUS2NEUTR: u32 = 1024;
pub const TG_FLAG_VARIABLE_PROTONS: u32 = 2048;
pub const TG_FLAG_HARD_ADD_REM_PROTONS: u32 = 4096;
pub const TG_FLAG_POINTED_EDGE_STEREO: u32 = 8192;
pub const TG_FLAG_PHOSPHINE_STEREO: u32 = 32768;
pub const TG_FLAG_ARSINE_STEREO: u32 = 65536;
pub const TG_FLAG_H_ALREADY_REMOVED: u32 = 131072;
pub const TG_FLAG_FIX_SP3_BUG: u32 = 262144;
pub const TG_FLAG_KETO_ENOL_TAUT: u32 = 524288;
pub const TG_FLAG_1_5_TAUT: u32 = 1048576;
pub const TG_FLAG_FIX_ISO_FIXEDH_BUG: u32 = 2097152;
pub const TG_FLAG_FIX_TERM_H_CHRG_BUG: u32 = 4194304;
pub const TG_FLAG_PT_22_00: u32 = 8388608;
pub const TG_FLAG_PT_16_00: u32 = 16777216;
pub const TG_FLAG_PT_06_00: u32 = 33554432;
pub const TG_FLAG_PT_39_00: u32 = 67108864;
pub const TG_FLAG_PT_13_00: u32 = 134217728;
pub const TG_FLAG_PT_18_00: u32 = 268435456;
pub const TG_FLAG_MOVE_HPLUS2NEUTR_DONE: u32 = 1;
pub const TG_FLAG_TEST_TAUT__ATOMS_DONE: u32 = 2;
pub const TG_FLAG_DISCONNECT_SALTS_DONE: u32 = 4;
pub const TG_FLAG_TEST_TAUT__SALTS_DONE: u32 = 8;
pub const TG_FLAG_MOVE_POS_CHARGES_DONE: u32 = 16;
pub const TG_FLAG_TEST_TAUT2_SALTS_DONE: u32 = 32;
pub const TG_FLAG_ALLOW_NO_NEGTV_O_DONE: u32 = 64;
pub const TG_FLAG_MERGE_TAUT_SALTS_DONE: u32 = 128;
pub const TG_FLAG_ALL_SALT_DONE: u32 = 168;
pub const TG_FLAG_DISCONNECT_COORD_DONE: u32 = 256;
pub const TG_FLAG_CHECK_VALENCE_COORD_DONE: u32 = 512;
pub const TG_FLAG_MOVE_CHARGE_COORD_DONE: u32 = 1024;
pub const TG_FLAG_FIX_ODD_THINGS_DONE: u32 = 2048;
pub const TG_FLAG_TEST_TAUT3_SALTS_DONE: u32 = 4096;
pub const TG_FLAG_FOUND_SALT_CHARGES_DONE: u32 = 8192;
pub const TG_FLAG_FOUND_ISOTOPIC_H_DONE: u32 = 16384;
pub const TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE: u32 = 32768;
pub const READ_INCHI_OUTPUT_INCHI: u32 = 1;
pub const READ_INCHI_SPLIT_OUTPUT: u32 = 2;
pub const READ_INCHI_KEEP_BALANCE_P: u32 = 4;
pub const READ_INCHI_TO_STRUCTURE: u32 = 8;
pub const INCHI_IOS_TYPE_NONE: u32 = 0;
pub const INCHI_IOS_TYPE_STRING: u32 = 1;
pub const INCHI_IOS_TYPE_FILE: u32 = 2;
pub const TRACE_MEMORY_LEAKS: u32 = 0;
pub const USE_ALLOCA: u32 = 0;
pub const DEBUG_POLYMERS: u32 = 0;
pub const POLYMERS_NO: u32 = 0;
pub const POLYMERS_MODERN: u32 = 1;
pub const POLYMERS_LEGACY: u32 = 2;
pub const POLYMERS_LEGACY_PLUS: u32 = 3;
pub const STEREO_AT_ZZ: u32 = 0;
pub const ALLOW_SUBSTRUCTURE_FILTERING: u32 = 1;
pub const MAX_ATOMS: u32 = 32766;
pub const NORMALLY_ALLOWED_INP_MAX_ATOMS: u32 = 1024;
pub const CHAR_MASK: u32 = 255;
pub const LEN_COORD: u32 = 10;
pub const NUM_COORD: u32 = 3;
pub const MAX_SDF_HEADER: u32 = 64;
pub const MAX_SDF_VALUE: u32 = 255;
pub const ATOM_EL_LEN: u32 = 6;
pub const ATOM_INFO_LEN: u32 = 36;
pub const MAXVAL: u32 = 20;
pub const MAX_STEREO_BONDS: u32 = 3;
pub const NUM_H_ISOTOPES: u32 = 3;
pub const ATW_H: u32 = 1;
pub const MIN_INPUT_BOND_TYPE: u32 = 1;
pub const MAX_INPUT_BOND_TYPE: u32 = 4;
pub const BOND_TYPE_SINGLE: u32 = 1;
pub const BOND_TYPE_DOUBLE: u32 = 2;
pub const BOND_TYPE_TRIPLE: u32 = 3;
pub const BOND_TYPE_ALTERN: u32 = 4;
pub const STEREO_SNGL_UP: u32 = 1;
pub const STEREO_SNGL_EITHER: u32 = 4;
pub const STEREO_SNGL_DOWN: u32 = 6;
pub const STEREO_DBLE_EITHER: u32 = 3;
pub const INPUT_STEREO_SNGL_UP: u32 = 1;
pub const INPUT_STEREO_SNGL_EITHER: u32 = 4;
pub const INPUT_STEREO_SNGL_DOWN: u32 = 6;
pub const INPUT_STEREO_DBLE_EITHER: u32 = 3;
pub const BOND_MARK_PARITY: u32 = 48;
pub const BOND_MARK_HIGHLIGHT: u32 = 64;
pub const BOND_MARK_ODD: u8 = 45u8;
pub const BOND_MARK_EVEN: u8 = 43u8;
pub const BOND_MARK_UNDF: u8 = 63u8;
pub const BOND_MARK_UNKN: u8 = 117u8;
pub const BOND_MARK_ERR: u8 = 42u8;
pub const SALT_DONOR_H: u32 = 1;
pub const SALT_DONOR_Neg: u32 = 2;
pub const SALT_ACCEPTOR: u32 = 4;
pub const SALT_p_DONOR: u32 = 8;
pub const SALT_p_ACCEPTOR: u32 = 16;
pub const SALT_DONOR_ALL: u32 = 27;
pub const SALT_DONOR_Neg2: u32 = 18;
pub const SALT_DONOR_H2: u32 = 9;
pub const SALT_DONOR: u32 = 3;
pub const SALT_SELECTED: u32 = 32;
pub const RADICAL_SINGLET: u32 = 1;
pub const RADICAL_DOUBLET: u32 = 2;
pub const RADICAL_TRIPLET: u32 = 3;
pub const METAL: u32 = 1;
pub const METAL2: u32 = 3;
pub const IS_METAL: u32 = 3;
pub const ZERO_ATW_DIFF: u32 = 127;
pub const TDP_LEN_LBL: u32 = 16;
pub const MAX_NUM_PATHS: u32 = 4;
pub const SD_FMT_END_OF_DATA: &[u8; 5] = b"$$$$\0";
pub const MOL_FMT_INPLINELEN: u32 = 204;
pub const MOL_FMT_MAXLINELEN: u32 = 200;
pub const MOL_FMT_PRESENT: u32 = 1;
pub const MOL_FMT_ABSENT: u32 = 0;
pub const MOL_FMT_QUERY: u32 = 0;
pub const MOL_FMT_CPSS: u32 = 0;
pub const MOL_FMT_REACT: u32 = 0;
pub const MOL_FMT_STRING_DATA: u8 = 83u8;
pub const MOL_FMT_CHAR_INT_DATA: u8 = 67u8;
pub const MOL_FMT_SHORT_INT_DATA: u8 = 78u8;
pub const MOL_FMT_LONG_INT_DATA: u8 = 76u8;
pub const MOL_FMT_DOUBLE_DATA: u8 = 68u8;
pub const MOL_FMT_FLOAT_DATA: u8 = 70u8;
pub const MOL_FMT_JUMP_TO_RIGHT: u8 = 74u8;
pub const MOL_FMT_INT_DATA: u8 = 73u8;
pub const MOL_FMT_MAX_VALUE_LEN: u32 = 32;
pub const MOL_FMT_M_STY_NON: u32 = 0;
pub const MOL_FMT_M_STY_SRU: u32 = 1;
pub const MOL_FMT_M_STY_MON: u32 = 2;
pub const MOL_FMT_M_STY_COP: u32 = 3;
pub const MOL_FMT_M_STY_MOD: u32 = 4;
pub const MOL_FMT_M_STY_CRO: u32 = 5;
pub const MOL_FMT_M_STY_MER: u32 = 6;
pub const MOL_FMT_M_SST_NON: u32 = 0;
pub const MOL_FMT_M_SST_ALT: u32 = 1;
pub const MOL_FMT_M_SST_RAN: u32 = 2;
pub const MOL_FMT_M_SST_BLK: u32 = 3;
pub const MOL_FMT_M_CONN_NON: u32 = 0;
pub const MOL_FMT_M_CONN_HT: u32 = 1;
pub const MOL_FMT_M_CONN_HH: u32 = 2;
pub const MOL_FMT_M_CONN_EU: u32 = 3;
pub const MOL_FMT_V3000_STENON: i32 = -1;
pub const MOL_FMT_V3000_STEABS: u32 = 1;
pub const MOL_FMT_V3000_STEREL: u32 = 2;
pub const MOL_FMT_V3000_STERAC: u32 = 3;
pub const MOL_FMT_V3000_INPLINELEN: u32 = 32004;
pub const MOL_FMT_V3000_MAXLINELEN: u32 = 32000;
pub const MOL_FMT_V3000_MAXFIELDLEN: u32 = 4096;
pub const ISOTOPIC_SHIFT_FLAG: u32 = 10000;
pub const CLOSING_STARRED_SRU_IS_A_MUST: u32 = 1;
pub const ALLOW_CLOSING_SRU_VIA_HIGHER_ORDER_BOND: u32 = 1;
pub const ALLOW_CLOSING_SRU_VIA_DIRADICAL: u32 = 1;
pub const CLOSING_SRU_NOT_APPLICABLE: u32 = 0;
pub const CLOSING_SRU_RING: u32 = 1;
pub const CLOSING_SRU_HIGHER_ORDER_BOND: u32 = 2;
pub const CLOSING_SRU_DIRADICAL: u32 = 3;
pub const CLOSING_SRU_RING_OPENED: u32 = 4;
pub const CLOSING_SRU_MOVED_BRACKETS: u32 = 5;
pub const ATT_NONE: u32 = 0;
pub const ATT_ACIDIC_CO: u32 = 1;
pub const ATT_ACIDIC_S: u32 = 2;
pub const ATT_OO: u32 = 4;
pub const ATT_ZOO: u32 = 8;
pub const ATT_NO: u32 = 16;
pub const ATT_N_O: u32 = 32;
pub const ATT_ATOM_N: u32 = 64;
pub const ATT_ATOM_P: u32 = 128;
pub const ATT_OTHER_NEG_O: u32 = 256;
pub const ATT_OTHER_ZO: u32 = 512;
pub const ATT_OH_MINUS: u32 = 1024;
pub const ATT_O_PLUS: u32 = 2048;
pub const ATT_PROTON: u32 = 4096;
pub const ATT_HalAnion: u32 = 8192;
pub const ATT_HalAcid: u32 = 16384;
pub const ATT_NP_MINUS_V23: u32 = 32768;
pub const AT_FLAG_ISO_H_POINT: u32 = 1;
pub const PERIODIC_NUMBER_H: u32 = 1;
pub const FlagSC_0D: u32 = 1;
pub const FlagSB_0D: u32 = 2;
pub const SB_PARITY_FLAG: u32 = 56;
pub const SB_PARITY_SHFT: u32 = 3;
pub const SB_PARITY_MASK: u32 = 7;
pub const NO_POLYMER: i32 = -1;
pub const POLYMER_REPRESENTATION_SOURCE_BASED: u32 = 1;
pub const POLYMER_REPRESENTATION_STRUCTURE_BASED: u32 = 2;
pub const POLYMER_REPRESENTATION_MIXED: u32 = 3;
pub const POLYMER_REPRESENTATION_UNRECOGNIZED: u32 = 4;
pub const ALLOW_MIXED_SRU_AND_MON: u32 = 1;
pub const POLYMER_STY_NON: u32 = 0;
pub const POLYMER_STY_SRU: u32 = 1;
pub const POLYMER_STY_MON: u32 = 2;
pub const POLYMER_STY_COP: u32 = 3;
pub const POLYMER_STY_MOD: u32 = 4;
pub const POLYMER_STY_CRO: u32 = 5;
pub const POLYMER_STY_MER: u32 = 6;
pub const POLYMER_SST_NON: u32 = 0;
pub const POLYMER_SST_ALT: u32 = 1;
pub const POLYMER_SST_RAN: u32 = 2;
pub const POLYMER_SST_BLK: u32 = 3;
pub const POLYMER_CONN_NON: u32 = 0;
pub const POLYMER_CONN_HT: u32 = 1;
pub const POLYMER_CONN_HH: u32 = 2;
pub const POLYMER_CONN_EU: u32 = 3;
pub const INF_STEREO_ABS: u32 = 1;
pub const INF_STEREO_REL: u32 = 2;
pub const INF_STEREO_RAC: u32 = 4;
pub const INF_STEREO_NORM: u32 = 8;
pub const INF_STEREO_INV: u32 = 16;
pub const INF_STEREO: u32 = 32;
pub const INF_STEREO_ABS_REL_RAC: u32 = 7;
pub const INF_STEREO_NORM_INV: u32 = 24;
pub const MAX_LEN_REMOVED_PROTONS: u32 = 128;
pub const ADD_LEN_STRUCT_FPTRS: u32 = 100;
pub const FLAG_INP_AT_CHIRAL: u32 = 1;
pub const FLAG_INP_AT_NONCHIRAL: u32 = 2;
pub const FLAG_SET_INP_AT_CHIRAL: u32 = 4;
pub const FLAG_SET_INP_AT_NONCHIRAL: u32 = 8;
pub const FLAG_SET_INP_LARGE_MOLS: u32 = 16;
pub const ALPHA_BASE: u32 = 27;
pub const INCHI_BUILD_PLATFORM: &[u8; 13] = b"Linux 64-bit\0";
pub const INCHI_BUILD_DEBUG: &[u8; 1] = b"\0";
pub const INCHI_SRC_REV: &[u8; 1] = b"\0";
pub const NUM_CHEM_ELEMENTS: u32 = 127;
pub const AT_ISO_SORT_KEY_MULT: u32 = 32;
pub const BYTE_BITS: u32 = 8;
pub const BOND_MASK: u32 = 15;
pub const BOND_BITS: u32 = 4;
pub const BOND_SINGLE: u32 = 1;
pub const BOND_DOUBLE: u32 = 2;
pub const BOND_TRIPLE: u32 = 3;
pub const BOND_ALTERN: u32 = 4;
pub const BOND_ALT_123: u32 = 5;
pub const BOND_ALT_13: u32 = 6;
pub const BOND_ALT_23: u32 = 7;
pub const BOND_TAUTOM: u32 = 8;
pub const BOND_ALT12NS: u32 = 9;
pub const BOND_NUMDIF: u32 = 9;
pub const BOND_TYPE_MASK: u32 = 15;
pub const BOND_MARK_ALL: u32 = 240;
pub const BOND_MARK_ALT12: u32 = 16;
pub const BOND_MARK_ALT123: u32 = 32;
pub const BOND_MARK_ALT13: u32 = 48;
pub const BOND_MARK_ALT23: u32 = 64;
pub const BOND_MARK_ALT12NS: u32 = 80;
pub const BOND_MARK_MASK: u32 = 112;
pub const BITS_PARITY: u32 = 7;
pub const MASK_CUMULENE_LEN: u32 = 56;
pub const KNOWN_PARITIES_EQL: u32 = 64;
pub const MAX_CUMULENE_LEN: u32 = 2;
pub const MULT_STEREOBOND: u32 = 8;
pub const AB_PARITY_NONE: u32 = 0;
pub const AB_PARITY_ODD: u32 = 1;
pub const AB_PARITY_EVEN: u32 = 2;
pub const AB_PARITY_UNKN: u32 = 3;
pub const AB_PARITY_UNDF: u32 = 4;
pub const AB_PARITY_IISO: u32 = 5;
pub const AB_PARITY_CALC: u32 = 6;
pub const AB_PARITY_0D: u32 = 8;
pub const AB_INV_PARITY_BITS: u32 = 3;
pub const AB_MAX_KNOWN_PARITY: u32 = 4;
pub const AB_MIN_KNOWN_PARITY: u32 = 1;
pub const AB_MAX_PART_DEFINED_PARITY: u32 = 3;
pub const AB_MIN_PART_DEFINED_PARITY: u32 = 1;
pub const AB_MAX_WELL_DEFINED_PARITY: u32 = 2;
pub const AB_MIN_WELL_DEFINED_PARITY: u32 = 1;
pub const AB_MIN_ILL_DEFINED_PARITY: u32 = 3;
pub const AB_MAX_ILL_DEFINED_PARITY: u32 = 4;
pub const AB_MAX_ANY_PARITY: u32 = 4;
pub const AB_MIN_ANY_PARITY: u32 = 1;
pub const AMBIGUOUS_STEREO: u32 = 1;
pub const AMBIGUOUS_STEREO_ATOM: u32 = 2;
pub const AMBIGUOUS_STEREO_BOND: u32 = 4;
pub const AMBIGUOUS_STEREO_ATOM_ISO: u32 = 8;
pub const AMBIGUOUS_STEREO_BOND_ISO: u32 = 16;
pub const AMBIGUOUS_STEREO_ERROR: u32 = 32;
pub const MIN_DOT_PROD: u32 = 50;
pub const ALWAYS_SET_STEREO_PARITY: u32 = 0;
pub const NO_ISOLATED_NON_6RING_AROM_BOND: u32 = 0;
pub const SAVE_6_AROM_CENTERS: u32 = 0;
pub const REQ_MODE_BASIC: u32 = 1;
pub const REQ_MODE_TAUT: u32 = 2;
pub const REQ_MODE_ISO: u32 = 4;
pub const REQ_MODE_NON_ISO: u32 = 8;
pub const REQ_MODE_STEREO: u32 = 16;
pub const REQ_MODE_ISO_STEREO: u32 = 32;
pub const REQ_MODE_NOEQ_STEREO: u32 = 64;
pub const REQ_MODE_REDNDNT_STEREO: u32 = 128;
pub const REQ_MODE_NO_ALT_SBONDS: u32 = 256;
pub const REQ_MODE_RELATIVE_STEREO: u32 = 512;
pub const REQ_MODE_RACEMIC_STEREO: u32 = 1024;
pub const REQ_MODE_SC_IGN_ALL_UU: u32 = 2048;
pub const REQ_MODE_SB_IGN_ALL_UU: u32 = 4096;
pub const REQ_MODE_CHIR_FLG_STEREO: u32 = 8192;
pub const REQ_MODE_DIFF_UU_STEREO: u32 = 16384;
pub const REQ_MODE_MIN_SB_RING_MASK: u32 = 983040;
pub const REQ_MODE_MIN_SB_RING_SHFT: u32 = 16;
pub const REQ_MODE_DEFAULT: u32 = 31;
pub const WARN_FAILED_STEREO: u32 = 1;
pub const WARN_FAILED_ISOTOPIC: u32 = 2;
pub const WARN_FAILED_ISOTOPIC_STEREO: u32 = 4;
pub const ERR_NO_CANON_RESULTS: u32 = 8;
pub const CMP_COMPONENTS: u32 = 1;
pub const CMP_COMPONENTS_NONISO: u32 = 2;
pub const CMP_COMPONENTS_NONTAUT: u32 = 4;
pub const INCHI_FLAG_ACID_TAUT: u32 = 1;
pub const INCHI_FLAG_REL_STEREO: u32 = 2;
pub const INCHI_FLAG_RAC_STEREO: u32 = 4;
pub const INCHI_FLAG_SC_IGN_ALL_UU: u32 = 8;
pub const INCHI_FLAG_SB_IGN_ALL_UU: u32 = 16;
pub const INCHI_FLAG_SC_IGN_ALL_ISO_UU: u32 = 32;
pub const INCHI_FLAG_SB_IGN_ALL_ISO_UU: u32 = 64;
pub const INCHI_FLAG_HARD_ADD_REM_PROTON: u32 = 128;
pub const INCHI_OUT_NO_AUX_INFO: u32 = 1;
pub const INCHI_OUT_SHORT_AUX_INFO: u32 = 2;
pub const INCHI_OUT_ONLY_AUX_INFO: u32 = 4;
pub const INCHI_OUT_EMBED_REC: u32 = 8;
pub const INCHI_OUT_SDFILE_ONLY: u32 = 16;
pub const INCHI_OUT_XML: u32 = 32;
pub const INCHI_OUT_PLAIN_TEXT: u32 = 64;
pub const INCHI_OUT_PLAIN_TEXT_COMMENTS: u32 = 128;
pub const INCHI_OUT_XML_TEXT_COMMENTS: u32 = 256;
pub const INCHI_OUT_WINCHI_WINDOW: u32 = 512;
pub const INCHI_OUT_TABBED_OUTPUT: u32 = 1024;
pub const INCHI_OUT_SDFILE_ATOMS_DT: u32 = 2048;
pub const INCHI_OUT_SDFILE_SPLIT: u32 = 4096;
pub const INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG: u32 = 8192;
pub const INCHI_OUT_STDINCHI: u32 = 16384;
pub const INCHI_OUT_SAVEOPT: u32 = 32768;
pub const INCHI_OUT_INCHI_GEN_ERROR: u32 = 1;
pub const INCHI_OUT_MISMATCH_AS_ERROR: u32 = 2;
pub const INCHI_OUT_REQ_LARGE_MOL: u32 = 4;
pub const INCHI_OUT_REQ_POLYMERS: u32 = 8;
pub const SAVE_OPT_SLUUD: u32 = 1;
pub const SAVE_OPT_SUU: u32 = 2;
pub const SAVE_OPT_FIXEDH: u32 = 4;
pub const SAVE_OPT_RECMET: u32 = 8;
pub const SAVE_OPT_KET: u32 = 16;
pub const SAVE_OPT_15T: u32 = 32;
pub const SAVE_OPT_PT_22_00: u32 = 64;
pub const SAVE_OPT_PT_16_00: u32 = 128;
pub const SAVE_OPT_PT_06_00: u32 = 256;
pub const SAVE_OPT_PT_39_00: u32 = 512;
pub const SAVE_OPT_PT_13_00: u32 = 1024;
pub const SAVE_OPT_PT_18_00: u32 = 2048;
pub const INCHI_OUT_PRINT_OPTIONS: u32 = 200;
pub const BN_MAX_ALTP: u32 = 16;
pub const BNS_EDGE_FORBIDDEN_MASK: u32 = 1;
pub const BNS_EDGE_FORBIDDEN_TEMP: u32 = 2;
pub const BNS_EDGE_FORBIDDEN_TEST: u32 = 4;
pub const BNS_VERT_TYPE_ATOM: u32 = 1;
pub const BNS_VERT_TYPE_ENDPOINT: u32 = 2;
pub const BNS_VERT_TYPE_TGROUP: u32 = 4;
pub const BNS_VERT_TYPE_C_POINT: u32 = 8;
pub const BNS_VERT_TYPE_C_GROUP: u32 = 16;
pub const BNS_VERT_TYPE_SUPER_TGROUP: u32 = 32;
pub const BNS_VERT_TYPE_TEMP: u32 = 64;
pub const BNS_VERT_TYPE__AUX: u32 = 128;
pub const BNS_VERT_TYPE_C_NEGATIVE: u32 = 256;
pub const BNS_VERT_TYPE_ACID: u32 = 512;
pub const BNS_VERT_TYPE_CARBON_GR: u32 = 1024;
pub const BNS_VERT_TYPE_METAL_GR: u32 = 2048;
pub const BNS_VERT_TYPE_ANY_GROUP: u32 = 52;
pub const BNS_VT_C_POS: u32 = 16;
pub const BNS_VT_C_NEG: u32 = 272;
pub const BNS_VT_C_POS_C: u32 = 1040;
pub const BNS_VT_C_NEG_C: u32 = 1296;
pub const BNS_VT_C_POS_M: u32 = 2064;
pub const BNS_VT_C_NEG_M: u32 = 2320;
pub const BNS_VT_M_GROUP: u32 = 2048;
pub const BNS_VT_C_POS_ALL: u32 = 48;
pub const BNS_VT_C_NEG_ALL: u32 = 304;
pub const BNS_VT_CHRG_STRUCT: u32 = 192;
pub const BNS_VT_YVCONNECTOR: u32 = 128;
pub const BNS_ADD_SUPER_TGROUP: u32 = 1;
pub const NUM_KINDS_OF_GROUPS: u32 = 2;
pub const BNS_ADD_ATOMS: u32 = 2;
pub const BNS_ADD_EDGES: u32 = 1;
pub const Vertex_s: u32 = 0;
pub const Vertex_t: u32 = 1;
pub const NO_VERTEX: i32 = -2;
pub const BLOSSOM_BASE: i32 = -1;
pub const ADD_CAPACITY_RADICAL: u32 = 1;
pub const MAX_BOND_EDGE_CAP: u32 = 2;
pub const AROM_BOND_EDGE_CAP: u32 = 1;
pub const MAX_TGROUP_EDGE_CAP: u32 = 2;
pub const EDGE_FLOW_ST_MASK: u32 = 16383;
pub const EDGE_FLOW_ST_PATH: u32 = 16384;
pub const EDGE_FLOW_MASK: u32 = 16383;
pub const EDGE_FLOW_PATH: u32 = 16384;
pub const MAX_ALT_AATG_ARRAY_LEN: u32 = 127;
pub const AATG_MARK_IN_PATH: u32 = 1;
pub const AATG_MARK_WAS_IN_PATH: u32 = 2;
pub const AATG_MARK_MAIN_TYPE: u32 = 4;
pub const AATG_MARK_OTHER_TYPE: u32 = 8;
pub const ALT_PATH_MODE_TAUTOM: u32 = 1;
pub const ALT_PATH_MODE_CHARGE: u32 = 2;
pub const ALT_PATH_MODE_4_SALT: u32 = 3;
pub const ALT_PATH_MODE_4_SALT2: u32 = 4;
pub const ALT_PATH_MODE_REM2H_CHG: u32 = 5;
pub const ALT_PATH_MODE_ADD2H_CHG: u32 = 6;
pub const ALT_PATH_MODE_REM2H_TST: u32 = 7;
pub const ALT_PATH_MODE_ADD2H_TST: u32 = 8;
pub const ALT_PATH_MODE_REM_PROTON: u32 = 9;
pub const ALT_PATH_MODE_TAUTOM_KET: u32 = 10;
pub const ALT_PATH_MODE_TAUTOM_PT_22_00: u32 = 11;
pub const ALT_PATH_MODE_TAUTOM_PT_16_00: u32 = 12;
pub const ALT_PATH_MODE_TAUTOM_PT_06_00: u32 = 13;
pub const ALT_PATH_MODE_TAUTOM_PT_39_00: u32 = 14;
pub const ALT_PATH_MODE_TAUTOM_PT_13_00: u32 = 15;
pub const ALT_PATH_MODE_TAUTOM_PT_18_00: u32 = 16;
pub const BNS_EF_CHNG_FLOW: u32 = 1;
pub const BNS_EF_RSTR_FLOW: u32 = 2;
pub const BNS_EF_CHNG_RSTR: u32 = 3;
pub const BNS_EF_CHNG_BONDS: u32 = 4;
pub const BNS_EF_ALTR_BONDS: u32 = 8;
pub const BNS_EF_UPD_RAD_ORI: u32 = 16;
pub const BNS_EF_SET_NOSTEREO: u32 = 32;
pub const BNS_EF_UPD_H_CHARGE: u32 = 64;
pub const BNS_EF_SAVE_ALL: u32 = 21;
pub const BNS_EF_ALTR_NS: u32 = 40;
pub const BNS_EF_RAD_SRCH: u32 = 128;
pub const INCHI_STRBUF_INITIAL_SIZE: u32 = 262144;
pub const INCHI_STRBUF_SIZE_INCREMENT: u32 = 262144;
pub const INCHI_STRBUF_SMALLER_INITIAL_SIZE: u32 = 1024;
pub const INCHI_STRBUF_SMALLER_SIZE_INCREMENT: u32 = 4096;
pub const T_NUM_NO_ISOTOPIC: u32 = 2;
pub const T_NUM_ISOTOPIC: u32 = 3;
pub const T_GROUP_HDR_LEN: u32 = 3;
pub const T_GROUP_ISOWT_MULT: u32 = 1024;
pub const TGSO_CURR_ORDER: u32 = 0;
pub const TGSO_SYMM_RANK: u32 = 1;
pub const TGSO_SYMM_IORDER: u32 = 2;
pub const TGSO_SYMM_IRANK: u32 = 3;
pub const TGSO_TOTAL_LEN: u32 = 4;
pub const FLAG_PROTON_NPO_SIMPLE_REMOVED: u32 = 1;
pub const FLAG_PROTON_NP_HARD_REMOVED: u32 = 2;
pub const FLAG_PROTON_AC_SIMPLE_ADDED: u32 = 4;
pub const FLAG_PROTON_AC_SIMPLE_REMOVED: u32 = 8;
pub const FLAG_PROTON_AC_HARD_REMOVED: u32 = 16;
pub const FLAG_PROTON_AC_HARD_ADDED: u32 = 32;
pub const FLAG_PROTON_CHARGE_CANCEL: u32 = 64;
pub const FLAG_PROTON_SINGLE_REMOVED: u32 = 128;
pub const FLAG_NORM_CONSIDER_TAUT: u32 = 255;
pub const FLAG_FORCE_SALT_TAUT: u32 = 50;
pub const CANON_FLAG_NO_H_RECANON: u32 = 1;
pub const CANON_FLAG_NO_TAUT_H_DIFF: u32 = 2;
pub const CANON_FLAG_ISO_ONLY_NON_TAUT_DIFF: u32 = 4;
pub const CANON_FLAG_ISO_TAUT_DIFF: u32 = 8;
pub const CANON_FLAG_ISO_FIXED_H_DIFF: u32 = 16;
pub const _HI: u32 = 1;
pub const _LO: u32 = 0;
pub const NEIGH_LIST_LEN: u32 = 4;
pub const U_LONG_LEN: u32 = 2;
pub const MOL_PART_MASK: i32 = -8;
pub const _IS_OKAY: u32 = 0;
pub const _IS_WARNING: u32 = 1;
pub const _IS_ERROR: u32 = 2;
pub const _IS_FATAL: u32 = 3;
pub const _IS_UNKNOWN: u32 = 4;
pub const _IS_EOF: i32 = -1;
pub const _IS_SKIP: i32 = -2;
pub const CT_ERR_FIRST: i32 = -30000;
pub const CT_OVERFLOW: i32 = -30000;
pub const CT_LEN_MISMATCH: i32 = -30001;
pub const CT_OUT_OF_RAM: i32 = -30002;
pub const CT_RANKING_ERR: i32 = -30003;
pub const CT_ISOCOUNT_ERR: i32 = -30004;
pub const CT_TAUCOUNT_ERR: i32 = -30005;
pub const CT_ISOTAUCOUNT_ERR: i32 = -30006;
pub const CT_MAPCOUNT_ERR: i32 = -30007;
pub const CT_TIMEOUT_ERR: i32 = -30008;
pub const CT_ISO_H_ERR: i32 = -30009;
pub const CT_STEREOCOUNT_ERR: i32 = -30010;
pub const CT_ATOMCOUNT_ERR: i32 = -30011;
pub const CT_STEREOBOND_ERROR: i32 = -30012;
pub const CT_USER_QUIT_ERR: i32 = -30013;
pub const CT_REMOVE_STEREO_ERR: i32 = -30014;
pub const CT_CALC_STEREO_ERR: i32 = -30015;
pub const CT_CANON_ERR: i32 = -30016;
pub const CT_STEREO_CANON_ERR: i32 = -30017;
pub const CT_WRONG_FORMULA: i32 = -30018;
pub const CT_UNKNOWN_ERR: i32 = -30019;
pub const CT_ERR_MIN: i32 = -30019;
pub const CT_ERR_MAX: i32 = -30000;
pub const BNS_ERR: i32 = -9999;
pub const BNS_WRONG_PARMS: i32 = -9999;
pub const BNS_OUT_OF_RAM: i32 = -9998;
pub const BNS_PROGRAM_ERR: i32 = -9997;
pub const BNS_ALTPATH_OVFL: i32 = -9996;
pub const BNS_BOND_ERR: i32 = -9995;
pub const BNS_VERT_NUM_ERR: i32 = -9994;
pub const BNS_VERT_EDGE_OVFL: i32 = -9993;
pub const BNS_SET_ALTP_ERR: i32 = -9992;
pub const BNS_CPOINT_ERR: i32 = -9991;
pub const BNS_CANT_SET_BOND: i32 = -9990;
pub const BNS_CAP_FLOW_ERR: i32 = -9989;
pub const BNS_RADICAL_ERR: i32 = -9988;
pub const BNS_REINIT_ERR: i32 = -9987;
pub const BNS_ALTBOND_ERR: i32 = -9986;
pub const BNS_TIMEOUT: i32 = -9985;
pub const BNS_MAX_ERR_VALUE: i32 = -9980;
pub const INCHI_INP_ERROR_ERR: u32 = 40;
pub const INCHI_INP_ERROR_RET: i32 = -1;
pub const INCHI_INP_FATAL_ERR: u32 = 1;
pub const INCHI_INP_FATAL_RET: u32 = 0;
pub const INCHI_INP_EOF_ERR: u32 = 11;
pub const INCHI_INP_EOF_RET: u32 = 0;
pub const LOG_MASK_WARN: u32 = 1;
pub const LOG_MASK_ERR: u32 = 2;
pub const LOG_MASK_FATAL: u32 = 4;
pub const LOG_MASK_ALL: u32 = 7;
pub const LOG_MASK_NO_WARN: u32 = 6;
pub const USER_ACTION_QUIT: u32 = 1;
pub const STR_ERR_LEN: u32 = 256;
pub const ICR_MAX_ENDP_IN1_ONLY: u32 = 32;
pub const ICR_MAX_ENDP_IN2_ONLY: u32 = 32;
pub const ICR_MAX_DIFF_FIXED_H: u32 = 32;
pub const ICR_MAX_SB_IN1_ONLY: u32 = 32;
pub const ICR_MAX_SB_IN2_ONLY: u32 = 32;
pub const ICR_MAX_SC_IN1_ONLY: u32 = 32;
pub const ICR_MAX_SC_IN2_ONLY: u32 = 32;
pub const ICR_MAX_SB_UNDF: u32 = 32;
pub const ICR_MAX_SC_UNDF: u32 = 32;
pub const EQL_EXISTS: u32 = 1;
pub const EQL_SP3: u32 = 2;
pub const EQL_SP3_INV: u32 = 4;
pub const EQL_SP2: u32 = 8;
pub const EQL_EQU: u32 = 0;
pub const EQL_EQU_TG: u32 = 1;
pub const EQL_EQU_ISO: u32 = 2;
pub const EQL_NUM: u32 = 0;
pub const EQL_NUM_INV: u32 = 1;
pub const EQL_NUM_ISO: u32 = 2;
pub const FLAG_SORT_PRINT_TRANSPOS_BAS: u32 = 1;
pub const FLAG_SORT_PRINT_TRANSPOS_REC: u32 = 2;
pub const FLAG_SORT_PRINT_NO_NFIX_H_BAS: u32 = 4;
pub const FLAG_SORT_PRINT_NO_NFIX_H_REC: u32 = 8;
pub const FLAG_SORT_PRINT_NO_IFIX_H_BAS: u32 = 16;
pub const FLAG_SORT_PRINT_NO_IFIX_H_REC: u32 = 32;
pub const FLAG_SORT_PRINT_ReChI_PREFIX: u32 = 64;
pub const ESC_KEY: u32 = 27;
pub const INCHI_SEGM_BUFLEN: u32 = 524288;
pub const PRINT_INCHI_MAX_TAG_LEN: u32 = 64;
pub const COMP_ORIG_0_MAIN: u32 = 1;
pub const COMP_ORIG_0_RECN: u32 = 2;
pub const COMP_PREP_0_MAIN: u32 = 4;
pub const COMP_PREP_0_RECN: u32 = 8;
pub const COMP_ORIG_1_MAIN: u32 = 16;
pub const COMP_ORIG_1_RECN: u32 = 32;
pub const ITEM_DELIMETER: &[u8; 2] = b",\0";
pub const EXTRA_SPACE: &[u8; 1] = b"\0";
pub const COMMA_EXTRA_SPACE: &[u8; 2] = b",\0";
pub const LEN_EXTRA_SPACE: u32 = 0;
pub const CT_MODE_NO_ORPHANS: u32 = 1;
pub const CT_MODE_ABC_NUMBERS: u32 = 2;
pub const CT_MODE_ATOM_COUNTS: u32 = 4;
pub const CT_MODE_PREDECESSORS: u32 = 8;
pub const CT_MODE_EQL_H_TOGETHER: u32 = 16;
pub const CT_MODE_ABC_NUM_CLOSURES: u32 = 32;
pub const iiSTEREO: u32 = 1;
pub const iiSTEREO_INV: u32 = 2;
pub const iiNUMB: u32 = 4;
pub const iiEQU: u32 = 8;
pub const iitISO: u32 = 16;
pub const iitNONTAUT: u32 = 32;
pub const iiEq2NONTAUT: u32 = 64;
pub const iiEq2ISO: u32 = 128;
pub const iiEq2INV: u32 = 256;
pub const iiEmpty: u32 = 512;
pub const QUEUE_QINT: u32 = 1;
pub const RI_ERR_ALLOC: i32 = -1;
pub const RI_ERR_SYNTAX: i32 = -2;
pub const RI_ERR_PROGR: i32 = -3;
pub const RI_ERR_EOL: i32 = -4;
pub const RI_ERR_EOF: u32 = 0;
pub const RI_ERR_MISMATCH: i32 = -9;
pub const NO_VALUE_INT: u32 = 9999;
pub const NOT_READ_INT: u32 = 9998;
pub const VALUE_OCTET: u32 = 8;
pub const INC_EDGE_LIST_DEFAULT: u32 = 64;
pub const BFS_Q_CLEAR: i32 = -1;
pub const BFS_Q_FREE: i32 = -2;
pub const EXTRACT_STRUCT_NUMBER: u32 = 1;
pub const FIX_STEREO_BOND_ORDER: u32 = 0;
pub const METAL_FREE_CHARGE_VAL: u32 = 1;
pub const ALLOW_METAL_BOND_ZERO: u32 = 1;
pub const INIT_METAL_BOND_ZERO: u32 = 0;
pub const INIT_METAL_BOND_FLOW: u32 = 1;
pub const I2A_FLAG_FIXEDH: u32 = 1;
pub const I2A_FLAG_RECMET: u32 = 2;
pub const ATYPE_H: u32 = 1;
pub const ATYPE_Na: u32 = 2;
pub const ATYPE_Mg: u32 = 3;
pub const ATYPE_B: u32 = 4;
pub const ATYPE_C: u32 = 5;
pub const ATYPE_N: u32 = 6;
pub const ATYPE_O: u32 = 7;
pub const ATYPE_Cl: u32 = 8;
pub const VAL_BASE: u32 = 1;
pub const VAL_MIN_CHARGE: i32 = -1;
pub const VAL_MAX_CHARGE: u32 = 1;
pub const VAL_NUMBER: u32 = 2;
pub const VAL_NEUTR_ORDER: u32 = 3;
pub const VAL_LENGTH: u32 = 4;
pub const VAL_NEGAT_CHARGE: u32 = 0;
pub const VAL_NEUTR_CHARGE: u32 = 1;
pub const VAL_POSIT_CHARGE: u32 = 2;
pub const INI_NUM_TCGROUPS: u32 = 16;
pub const INC_NUM_TCGROUPS: u32 = 16;
pub const EDGE_LIST_CLEAR: i32 = -1;
pub const EDGE_LIST_FREE: i32 = -2;
pub const BOND_MARK_STEREO: u32 = 16;
pub const BOND_TYPE_STEREO: u32 = 17;
pub const RESET_EDGE_FORBIDDEN_MASK: u32 = 0;
pub const TREAT_ATOM_AS_METAL: u32 = 99;
pub const EL_TYPE_O: u32 = 1;
pub const EL_TYPE_S: u32 = 2;
pub const EL_TYPE_N: u32 = 4;
pub const EL_TYPE_P: u32 = 8;
pub const EL_TYPE_C: u32 = 16;
pub const EL_TYPE_X: u32 = 32;
pub const EL_TYPE_MASK: u32 = 63;
pub const EL_TYPE_OSt: u32 = 256;
pub const EL_TYPE_PT: u32 = 512;
pub const MAX_CN_VAL: u32 = 3;
pub const cn_bits_N: u32 = 1;
pub const cn_bits_P: u32 = 2;
pub const cn_bits_M: u32 = 4;
pub const cn_bits_shift: u32 = 3;
pub const MAX_NUM_CN_BITS: u32 = 4;
pub const cn_bits_Me: i32 = -1;
pub const cnListIndexMe: u32 = 17;
pub const MAX_DIFF_FIXH: u32 = 256;
pub const MAX_DIFF_MOBH: u32 = 256;
pub const KEEP_METAL_EDGE_FLOW: u32 = 0;
pub const MOVE_CHARGES_FROM_HETEREO_TO_METAL: u32 = 0;
pub const FIX_ADD_PROTON_FOR_ADP: u32 = 0;
pub const PES_BIT_POINT_EDGE_STEREO: u32 = 1;
pub const PES_BIT_PHOSPHINE_STEREO: u32 = 2;
pub const PES_BIT_ARSINE_STEREO: u32 = 4;
pub const PES_BIT_FIX_SP3_BUG: u32 = 8;
pub const FIX_DOCANON_RETCODE_RESET_BUG: u32 = 1;
pub const ISOTOPIC_SHIFT_MAX: u32 = 100;
pub const NO_ATOM: i32 = -1;
pub const INCHI_STRING_PREFIX: &[u8; 7] = b"InChI=\0";
pub const LEN_INCHI_STRING_PREFIX: u32 = 6;
pub const INCHIKEY_OK: u32 = 0;
pub const INCHIKEY_UNKNOWN_ERROR: u32 = 1;
pub const INCHIKEY_EMPTY_INPUT: u32 = 2;
pub const INCHIKEY_INVALID_INCHI_PREFIX: u32 = 3;
pub const INCHIKEY_NOT_ENOUGH_MEMORY: u32 = 4;
pub const INCHIKEY_INVALID_INCHI: u32 = 20;
pub const INCHIKEY_INVALID_STD_INCHI: u32 = 21;
pub const IXA_ATOM_NATURAL_MASS: u32 = 0;
pub const IXA_EXT_MOLDATA_INVALID: i32 = -1;
pub const IXA_EXT_POLYMER_INVALID: i32 = -1;
pub const IXA_EXT_V3000_INVALID: i32 = -1;
pub const IXA_USES_SMART_ALLOCS: u32 = 1;
pub const LOGGING_ENABLED: u32 = 0;
pub const PERMAXATOMS: u32 = 32767;
pub const INCHI_LINE_LEN: u32 = 262144;
pub const INCHI_LINE_ADD: u32 = 32767;
pub const STB_SPRINTF_MIN: u32 = 512;
pub const INCHI_MAX_NUM_ARG: u32 = 32;
pub const PSTR_BUFFER_SIZE: u32 = 524288;
pub const MOL2INCHI_NO_RAM: u32 = 1001;
pub const MOL2INCHI_BAD_COMMAND_LINE: u32 = 1002;
pub const INCHIMOL_ATOMS_START_SIZE: u32 = 128;
pub const INCHIMOL_BONDS_START_SIZE: u32 = 128;
pub const INCHIMOL_STEREOS_START_SIZE: u32 = 64;
pub const INCHIMOL_MAX_ATOMS: u32 = 32767;
pub const INCHIMOL_POLYMERUNITS_START_SIZE: u32 = 2;
pub type FILE = SourceFile;
#[derive(Clone, Debug, PartialEq)]
pub struct tagOutputString {
    pub pStr: SourceMutPointer<i8>,
    pub nAllocatedLength: i32,
    pub nUsedLength: i32,
    pub nPtr: i32,
}
impl ::std::default::Default for tagOutputString {
    fn default() -> Self {
        Self {
            pStr: ::std::default::Default::default(),
            nAllocatedLength: ::std::default::Default::default(),
            nUsedLength: ::std::default::Default::default(),
            nPtr: ::std::default::Default::default(),
        }
    }
}
pub type INCHI_IOS_STRING = tagOutputString;
#[derive(Clone, Debug, PartialEq)]
pub struct tagOutputStream {
    pub s: INCHI_IOS_STRING,
    pub f: SourceMutPointer<FILE>,
    pub type_: i32,
}
impl ::std::default::Default for tagOutputStream {
    fn default() -> Self {
        Self {
            s: ::std::default::Default::default(),
            f: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
        }
    }
}
pub type INCHI_IOSTREAM = tagOutputStream;
pub type AT_NUMB = u16;
pub type AT_RANK = u16;
pub type NUM_H = i16;
pub type pAT_RANK = SourceMutPointer<AT_RANK>;
pub type ppAT_RANK = SourceMutPointer<pAT_RANK>;
pub type INCHI_MODE = u64;
pub type MOL_COORD = [i8; 32usize];
pub type S_CHAR = i8;
pub type U_CHAR = u8;
pub type S_SHORT = i16;
pub type U_SHORT = u16;
pub type clock_t = i64;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInchiTime {
    pub clockTime: clock_t,
}
impl ::std::default::Default for tagInchiTime {
    fn default() -> Self {
        Self {
            clockTime: ::std::default::Default::default(),
        }
    }
}
pub type inchiTime = tagInchiTime;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHI_CLOCK {
    pub m_MaxPositiveClock: clock_t,
    pub m_MinNegativeClock: clock_t,
    pub m_HalfMaxPositiveClock: clock_t,
    pub m_HalfMinNegativeClock: clock_t,
}
impl ::std::default::Default for tagINCHI_CLOCK {
    fn default() -> Self {
        Self {
            m_MaxPositiveClock: ::std::default::Default::default(),
            m_MinNegativeClock: ::std::default::Default::default(),
            m_HalfMaxPositiveClock: ::std::default::Default::default(),
            m_HalfMinNegativeClock: ::std::default::Default::default(),
        }
    }
}
pub type INCHI_CLOCK = tagINCHI_CLOCK;
pub const tagTblTypes_itBASIC: tagTblTypes = 0;
pub const tagTblTypes_itISOTOPIC: tagTblTypes = 1;
pub const tagTblTypes_itSTEREO: tagTblTypes = 2;
pub const tagTblTypes_TDP_NUM_PAR: tagTblTypes = 3;
pub type tagTblTypes = u32;
pub use self::tagTblTypes as TBL_TYPES;
pub const tagTblLabels_ilSHOWN: tagTblLabels = 0;
pub const tagTblLabels_TDP_NUM_LBL: tagTblLabels = 1;
pub type tagTblLabels = u32;
pub use self::tagTblLabels as TBL_LABELS;
#[derive(Clone, Debug, PartialEq)]
pub struct tagTblDrawPatms {
    pub ReqShownFoundTxt: [[i8; 16usize]; 1usize],
    pub ReqShownFound: [[i8; 3usize]; 1usize],
    pub nOrientation: i32,
    pub bDrawTbl: i32,
}
impl ::std::default::Default for tagTblDrawPatms {
    fn default() -> Self {
        Self {
            ReqShownFoundTxt: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            ReqShownFound: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            nOrientation: ::std::default::Default::default(),
            bDrawTbl: ::std::default::Default::default(),
        }
    }
}
pub type TBL_DRAW_PARMS = tagTblDrawPatms;
#[derive(Clone, Debug, PartialEq)]
pub struct tagDrawParmsSettings {
    pub tdp: SourceMutPointer<TBL_DRAW_PARMS>,
    pub ulDisplTime: u64,
    pub bOrigAtom: i32,
    pub nFontSize: i32,
}
impl ::std::default::Default for tagDrawParmsSettings {
    fn default() -> Self {
        Self {
            tdp: ::std::default::Default::default(),
            ulDisplTime: ::std::default::Default::default(),
            bOrigAtom: ::std::default::Default::default(),
            nFontSize: ::std::default::Default::default(),
        }
    }
}
pub type SET_DRAW_PARMS = tagDrawParmsSettings;
#[derive(Clone, Debug, PartialEq)]
pub struct tagReturnedDrawParms {
    pub bEsc: i32,
}
impl ::std::default::Default for tagReturnedDrawParms {
    fn default() -> Self {
        Self {
            bEsc: ::std::default::Default::default(),
        }
    }
}
pub type RET_DRAW_PARMS = tagReturnedDrawParms;
#[derive(Clone, Debug, PartialEq)]
pub struct tagPersistDrawParms {
    pub rcPict: [i32; 4usize],
}
impl ::std::default::Default for tagPersistDrawParms {
    fn default() -> Self {
        Self {
            rcPict: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type PER_DRAW_PARMS = tagPersistDrawParms;
#[derive(Clone, Debug, PartialEq)]
pub struct tagDrawParms {
    pub sdp: SET_DRAW_PARMS,
    pub rdp: RET_DRAW_PARMS,
    pub pdp: SourceMutPointer<PER_DRAW_PARMS>,
}
impl ::std::default::Default for tagDrawParms {
    fn default() -> Self {
        Self {
            sdp: ::std::default::Default::default(),
            rdp: ::std::default::Default::default(),
            pdp: ::std::default::Default::default(),
        }
    }
}
pub type DRAW_PARMS = tagDrawParms;
pub const tagInputType_INPUT_NONE: tagInputType = 0;
pub const tagInputType_INPUT_MOLFILE: tagInputType = 1;
pub const tagInputType_INPUT_SDFILE: tagInputType = 2;
pub const tagInputType_INPUT_INCHI_XML: tagInputType = 3;
pub const tagInputType_INPUT_INCHI_PLAIN: tagInputType = 4;
pub const tagInputType_INPUT_CMLFILE: tagInputType = 5;
pub const tagInputType_INPUT_INCHI: tagInputType = 6;
pub const tagInputType_INPUT_MAX: tagInputType = 7;
pub type tagInputType = u32;
pub use self::tagInputType as INPUT_TYPE;
pub const tagInChIHashCalc_INCHIHASH_NONE: tagInChIHashCalc = 0;
pub const tagInChIHashCalc_INCHIHASH_KEY: tagInChIHashCalc = 1;
pub const tagInChIHashCalc_INCHIHASH_KEY_XTRA1: tagInChIHashCalc = 2;
pub const tagInChIHashCalc_INCHIHASH_KEY_XTRA2: tagInChIHashCalc = 3;
pub const tagInChIHashCalc_INCHIHASH_KEY_XTRA1_XTRA2: tagInChIHashCalc = 4;
pub type tagInChIHashCalc = u32;
pub use self::tagInChIHashCalc as INCHI_HASH_CALC;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInputParms {
    pub szSdfDataHeader: [i8; 65usize],
    pub pSdfLabel: SourceMutPointer<i8>,
    pub pSdfValue: SourceMutPointer<i8>,
    pub lSdfId: u64,
    pub lMolfileNumber: i64,
    pub path: [SourceConstPointer<i8>; 4usize],
    pub num_paths: i32,
    pub first_struct_number: i64,
    pub last_struct_number: i64,
    pub nInputType: INPUT_TYPE,
    pub nMode: INCHI_MODE,
    pub bAbcNumbers: i32,
    pub bINChIOutputOptions: i32,
    pub bINChIOutputOptions2: i32,
    pub bCtPredecessors: i32,
    pub bDisplayEachComponentINChI: i32,
    pub msec_MaxTime: i64,
    pub msec_LeftTime: i64,
    pub ulDisplTime: i64,
    pub bDisplay: i32,
    pub bDisplayIfRestoreWarnings: i32,
    pub bMergeAllInputStructures: i32,
    pub bSaveWarningStructsAsProblem: i32,
    pub bSaveAllGoodStructsAsProblem: i32,
    pub bGetSdfileId: i32,
    pub bGetMolfileNumber: i32,
    pub bCompareComponents: i32,
    pub bDisplayCompositeResults: i32,
    pub bDoNotAddH: i32,
    pub bNoStructLabels: i32,
    pub bFixNonUniformDraw: i32,
    pub bCalcInChIHash: i32,
    pub bChiralFlag: i32,
    pub bAllowEmptyStructure: i32,
    pub bLargeMolecules: i32,
    pub bLooseTSACheck: i32,
    pub bPolymers: i32,
    pub bFoldPolymerSRU: i32,
    pub bFrameShiftScheme: i32,
    pub bNPZz: i32,
    pub bStereoAtZz: i32,
    pub bMergeHash: i32,
    pub bNoWarnings: i32,
    pub bHideInChI: i32,
    pub bTautFlags: INCHI_MODE,
    pub bTautFlagsDone: INCHI_MODE,
    pub bReadInChIOptions: i32,
    pub bRenumber: i32,
    pub bUnderivatize: i32,
    pub bRing2Chain: i32,
    pub bIgnoreUnchanged: i32,
    pub bFilterSS: i32,
}
impl ::std::default::Default for tagInputParms {
    fn default() -> Self {
        Self {
            szSdfDataHeader: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pSdfLabel: ::std::default::Default::default(),
            pSdfValue: ::std::default::Default::default(),
            lSdfId: ::std::default::Default::default(),
            lMolfileNumber: ::std::default::Default::default(),
            path: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_paths: ::std::default::Default::default(),
            first_struct_number: ::std::default::Default::default(),
            last_struct_number: ::std::default::Default::default(),
            nInputType: ::std::default::Default::default(),
            nMode: ::std::default::Default::default(),
            bAbcNumbers: ::std::default::Default::default(),
            bINChIOutputOptions: ::std::default::Default::default(),
            bINChIOutputOptions2: ::std::default::Default::default(),
            bCtPredecessors: ::std::default::Default::default(),
            bDisplayEachComponentINChI: ::std::default::Default::default(),
            msec_MaxTime: ::std::default::Default::default(),
            msec_LeftTime: ::std::default::Default::default(),
            ulDisplTime: ::std::default::Default::default(),
            bDisplay: ::std::default::Default::default(),
            bDisplayIfRestoreWarnings: ::std::default::Default::default(),
            bMergeAllInputStructures: ::std::default::Default::default(),
            bSaveWarningStructsAsProblem: ::std::default::Default::default(),
            bSaveAllGoodStructsAsProblem: ::std::default::Default::default(),
            bGetSdfileId: ::std::default::Default::default(),
            bGetMolfileNumber: ::std::default::Default::default(),
            bCompareComponents: ::std::default::Default::default(),
            bDisplayCompositeResults: ::std::default::Default::default(),
            bDoNotAddH: ::std::default::Default::default(),
            bNoStructLabels: ::std::default::Default::default(),
            bFixNonUniformDraw: ::std::default::Default::default(),
            bCalcInChIHash: ::std::default::Default::default(),
            bChiralFlag: ::std::default::Default::default(),
            bAllowEmptyStructure: ::std::default::Default::default(),
            bLargeMolecules: ::std::default::Default::default(),
            bLooseTSACheck: ::std::default::Default::default(),
            bPolymers: ::std::default::Default::default(),
            bFoldPolymerSRU: ::std::default::Default::default(),
            bFrameShiftScheme: ::std::default::Default::default(),
            bNPZz: ::std::default::Default::default(),
            bStereoAtZz: ::std::default::Default::default(),
            bMergeHash: ::std::default::Default::default(),
            bNoWarnings: ::std::default::Default::default(),
            bHideInChI: ::std::default::Default::default(),
            bTautFlags: ::std::default::Default::default(),
            bTautFlagsDone: ::std::default::Default::default(),
            bReadInChIOptions: ::std::default::Default::default(),
            bRenumber: ::std::default::Default::default(),
            bUnderivatize: ::std::default::Default::default(),
            bRing2Chain: ::std::default::Default::default(),
            bIgnoreUnchanged: ::std::default::Default::default(),
            bFilterSS: ::std::default::Default::default(),
        }
    }
}
pub type INPUT_PARMS = tagInputParms;
pub const tagFrameShifScheme_FSS_STARS_CYCLED: tagFrameShifScheme = 0;
pub const tagFrameShifScheme_FSS_NONE: tagFrameShifScheme = 1;
pub const tagFrameShifScheme_FSS_STARS_CYCLED_SORTED: tagFrameShifScheme = 2;
pub const tagFrameShifScheme_FSS_STARS_OPENED: tagFrameShifScheme = 3;
pub const tagFrameShifScheme_FSS_STARS_ENDS_OPENED: tagFrameShifScheme = 4;
pub type tagFrameShifScheme = u32;
pub use self::tagFrameShifScheme as FRAME_SHIFT_SCHEME;
#[derive(Clone, Debug, PartialEq)]
pub struct A_NUM_LISTS {
    pub lists: SourceMutPointer<SourceMutPointer<i32>>,
    pub allocated: i32,
    pub used: i32,
    pub increment: i32,
}
impl ::std::default::Default for A_NUM_LISTS {
    fn default() -> Self {
        Self {
            lists: ::std::default::Default::default(),
            allocated: ::std::default::Default::default(),
            used: ::std::default::Default::default(),
            increment: ::std::default::Default::default(),
        }
    }
}
pub type NUM_LISTS = A_NUM_LISTS;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINT_ARRAY {
    pub item: SourceMutPointer<i32>,
    pub allocated: i32,
    pub used: i32,
    pub increment: i32,
}
impl ::std::default::Default for tagINT_ARRAY {
    fn default() -> Self {
        Self {
            item: ::std::default::Default::default(),
            allocated: ::std::default::Default::default(),
            used: ::std::default::Default::default(),
            increment: ::std::default::Default::default(),
        }
    }
}
pub type INT_ARRAY = tagINT_ARRAY;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_SGROUP {
    pub id: i32,
    pub type_: i32,
    pub subtype: i32,
    pub conn: i32,
    pub label: i32,
    pub xbr1: [f64; 4usize],
    pub xbr2: [f64; 4usize],
    pub smt: [i8; 80usize],
    pub alist: INT_ARRAY,
    pub blist: INT_ARRAY,
}
impl ::std::default::Default for A_MOL_FMT_SGROUP {
    fn default() -> Self {
        Self {
            id: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            subtype: ::std::default::Default::default(),
            conn: ::std::default::Default::default(),
            label: ::std::default::Default::default(),
            xbr1: ::std::array::from_fn(|_| ::std::default::Default::default()),
            xbr2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            smt: ::std::array::from_fn(|_| ::std::default::Default::default()),
            alist: ::std::default::Default::default(),
            blist: ::std::default::Default::default(),
        }
    }
}
pub type MOL_FMT_SGROUP = A_MOL_FMT_SGROUP;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_SGROUPS {
    pub group: SourceMutPointer<SourceMutPointer<MOL_FMT_SGROUP>>,
    pub allocated: i32,
    pub used: i32,
    pub increment: i32,
}
impl ::std::default::Default for A_MOL_FMT_SGROUPS {
    fn default() -> Self {
        Self {
            group: ::std::default::Default::default(),
            allocated: ::std::default::Default::default(),
            used: ::std::default::Default::default(),
            increment: ::std::default::Default::default(),
        }
    }
}
pub type MOL_FMT_SGROUPS = A_MOL_FMT_SGROUPS;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_HEADER_BLOCK {
    pub molname: [i8; 201usize],
    pub line2: [i8; 201usize],
    pub user_initls: [i8; 3usize],
    pub prog_name: [i8; 9usize],
    pub month: i8,
    pub day: i8,
    pub year: i8,
    pub hour: i8,
    pub minute: i8,
    pub dim_code: [i8; 3usize],
    pub scaling_factor1: i16,
    pub scaling_factor2: f64,
    pub energy: f64,
    pub internal_regno: i64,
    pub comment: [i8; 81usize],
}
impl ::std::default::Default for A_MOL_FMT_HEADER_BLOCK {
    fn default() -> Self {
        Self {
            molname: ::std::array::from_fn(|_| ::std::default::Default::default()),
            line2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            user_initls: ::std::array::from_fn(|_| ::std::default::Default::default()),
            prog_name: ::std::array::from_fn(|_| ::std::default::Default::default()),
            month: ::std::default::Default::default(),
            day: ::std::default::Default::default(),
            year: ::std::default::Default::default(),
            hour: ::std::default::Default::default(),
            minute: ::std::default::Default::default(),
            dim_code: ::std::array::from_fn(|_| ::std::default::Default::default()),
            scaling_factor1: ::std::default::Default::default(),
            scaling_factor2: ::std::default::Default::default(),
            energy: ::std::default::Default::default(),
            internal_regno: ::std::default::Default::default(),
            comment: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type MOL_FMT_HEADER_BLOCK = A_MOL_FMT_HEADER_BLOCK;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_ATOM {
    pub fx: f64,
    pub fy: f64,
    pub fz: f64,
    pub symbol: [i8; 6usize],
    pub mass_difference: S_CHAR,
    pub charge: S_CHAR,
    pub radical: i8,
    pub stereo_parity: i8,
    pub valence: i8,
    pub my_n_impH: i8,
    pub display_tom: i8,
    pub atom_aliased_flag: i8,
}
impl ::std::default::Default for A_MOL_FMT_ATOM {
    fn default() -> Self {
        Self {
            fx: ::std::default::Default::default(),
            fy: ::std::default::Default::default(),
            fz: ::std::default::Default::default(),
            symbol: ::std::array::from_fn(|_| ::std::default::Default::default()),
            mass_difference: ::std::default::Default::default(),
            charge: ::std::default::Default::default(),
            radical: ::std::default::Default::default(),
            stereo_parity: ::std::default::Default::default(),
            valence: ::std::default::Default::default(),
            my_n_impH: ::std::default::Default::default(),
            display_tom: ::std::default::Default::default(),
            atom_aliased_flag: ::std::default::Default::default(),
        }
    }
}
pub type MOL_FMT_ATOM = A_MOL_FMT_ATOM;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_BOND {
    pub atnum1: i16,
    pub atnum2: i16,
    pub bond_type: i8,
    pub bond_stereo: i8,
}
impl ::std::default::Default for A_MOL_FMT_BOND {
    fn default() -> Self {
        Self {
            atnum1: ::std::default::Default::default(),
            atnum2: ::std::default::Default::default(),
            bond_type: ::std::default::Default::default(),
            bond_stereo: ::std::default::Default::default(),
        }
    }
}
pub type MOL_FMT_BOND = A_MOL_FMT_BOND;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_v3000 {
    pub n_non_star_atoms: i32,
    pub n_star_atoms: i32,
    pub atom_index_orig: SourceMutPointer<i32>,
    pub atom_index_fin: SourceMutPointer<i32>,
    pub n_sgroups: i32,
    pub n_3d_constraints: i32,
    pub n_collections: i32,
    pub n_non_haptic_bonds: i32,
    pub n_haptic_bonds: i32,
    pub haptic_bonds: SourceMutPointer<NUM_LISTS>,
    pub n_steabs: i32,
    pub steabs: SourceMutPointer<NUM_LISTS>,
    pub n_sterel: i32,
    pub sterel: SourceMutPointer<NUM_LISTS>,
    pub n_sterac: i32,
    pub sterac: SourceMutPointer<NUM_LISTS>,
}
impl ::std::default::Default for A_MOL_FMT_v3000 {
    fn default() -> Self {
        Self {
            n_non_star_atoms: ::std::default::Default::default(),
            n_star_atoms: ::std::default::Default::default(),
            atom_index_orig: ::std::default::Default::default(),
            atom_index_fin: ::std::default::Default::default(),
            n_sgroups: ::std::default::Default::default(),
            n_3d_constraints: ::std::default::Default::default(),
            n_collections: ::std::default::Default::default(),
            n_non_haptic_bonds: ::std::default::Default::default(),
            n_haptic_bonds: ::std::default::Default::default(),
            haptic_bonds: ::std::default::Default::default(),
            n_steabs: ::std::default::Default::default(),
            steabs: ::std::default::Default::default(),
            n_sterel: ::std::default::Default::default(),
            sterel: ::std::default::Default::default(),
            n_sterac: ::std::default::Default::default(),
            sterac: ::std::default::Default::default(),
        }
    }
}
pub type MOL_FMT_v3000 = A_MOL_FMT_v3000;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_CTAB {
    pub n_atoms: i32,
    pub n_bonds: i32,
    pub chiral_flag: i8,
    pub n_stext_entries: i16,
    pub n_property_lines: i16,
    pub follow_inchi_1_treating_iso_mass: i16,
    pub version_string: [i8; 7usize],
    pub atoms: SourceMutPointer<MOL_FMT_ATOM>,
    pub bonds: SourceMutPointer<MOL_FMT_BOND>,
    pub coords: SourceMutPointer<MOL_COORD>,
    pub sgroups: MOL_FMT_SGROUPS,
    pub v3000: SourceMutPointer<MOL_FMT_v3000>,
}
impl ::std::default::Default for A_MOL_FMT_CTAB {
    fn default() -> Self {
        Self {
            n_atoms: ::std::default::Default::default(),
            n_bonds: ::std::default::Default::default(),
            chiral_flag: ::std::default::Default::default(),
            n_stext_entries: ::std::default::Default::default(),
            n_property_lines: ::std::default::Default::default(),
            follow_inchi_1_treating_iso_mass: ::std::default::Default::default(),
            version_string: ::std::array::from_fn(|_| ::std::default::Default::default()),
            atoms: ::std::default::Default::default(),
            bonds: ::std::default::Default::default(),
            coords: ::std::default::Default::default(),
            sgroups: ::std::default::Default::default(),
            v3000: ::std::default::Default::default(),
        }
    }
}
pub type MOL_FMT_CTAB = A_MOL_FMT_CTAB;
#[derive(Clone, Debug, PartialEq)]
pub struct A_MOL_FMT_DATA {
    pub hdr: MOL_FMT_HEADER_BLOCK,
    pub ctab: MOL_FMT_CTAB,
}
impl ::std::default::Default for A_MOL_FMT_DATA {
    fn default() -> Self {
        Self {
            hdr: ::std::default::Default::default(),
            ctab: ::std::default::Default::default(),
        }
    }
}
pub type MOL_FMT_DATA = A_MOL_FMT_DATA;
pub type ST_CAP_FLOW = S_SHORT;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInputAtom {
    pub elname: [i8; 6usize],
    pub el_number: U_CHAR,
    pub neighbor: [AT_NUMB; 20usize],
    pub orig_at_number: AT_NUMB,
    pub orig_compt_at_numb: AT_NUMB,
    pub bond_stereo: [S_CHAR; 20usize],
    pub bond_type: [U_CHAR; 20usize],
    pub valence: S_CHAR,
    pub chem_bonds_valence: S_CHAR,
    pub num_H: S_CHAR,
    pub num_iso_H: [S_CHAR; 3usize],
    pub iso_atw_diff: S_CHAR,
    pub charge: S_CHAR,
    pub radical: S_CHAR,
    pub bAmbiguousStereo: S_CHAR,
    pub cFlags: S_CHAR,
    pub at_type: AT_NUMB,
    pub component: AT_NUMB,
    pub endpoint: AT_NUMB,
    pub c_point: AT_NUMB,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub bUsed0DParity: S_CHAR,
    pub p_parity: S_CHAR,
    pub p_orig_at_num: [AT_NUMB; 4usize],
    pub sb_ord: [S_CHAR; 3usize],
    pub sn_ord: [S_CHAR; 3usize],
    pub sb_parity: [S_CHAR; 3usize],
    pub sn_orig_at_num: [AT_NUMB; 3usize],
    pub bCutVertex: S_CHAR,
    pub nRingSystem: AT_NUMB,
    pub nNumAtInRingSystem: AT_NUMB,
    pub nBlockSystem: AT_NUMB,
}
impl ::std::default::Default for tagInputAtom {
    fn default() -> Self {
        Self {
            elname: ::std::array::from_fn(|_| ::std::default::Default::default()),
            el_number: ::std::default::Default::default(),
            neighbor: ::std::array::from_fn(|_| ::std::default::Default::default()),
            orig_at_number: ::std::default::Default::default(),
            orig_compt_at_numb: ::std::default::Default::default(),
            bond_stereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bond_type: ::std::array::from_fn(|_| ::std::default::Default::default()),
            valence: ::std::default::Default::default(),
            chem_bonds_valence: ::std::default::Default::default(),
            num_H: ::std::default::Default::default(),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            iso_atw_diff: ::std::default::Default::default(),
            charge: ::std::default::Default::default(),
            radical: ::std::default::Default::default(),
            bAmbiguousStereo: ::std::default::Default::default(),
            cFlags: ::std::default::Default::default(),
            at_type: ::std::default::Default::default(),
            component: ::std::default::Default::default(),
            endpoint: ::std::default::Default::default(),
            c_point: ::std::default::Default::default(),
            x: ::std::default::Default::default(),
            y: ::std::default::Default::default(),
            z: ::std::default::Default::default(),
            bUsed0DParity: ::std::default::Default::default(),
            p_parity: ::std::default::Default::default(),
            p_orig_at_num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sb_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sn_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sb_parity: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sn_orig_at_num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bCutVertex: ::std::default::Default::default(),
            nRingSystem: ::std::default::Default::default(),
            nNumAtInRingSystem: ::std::default::Default::default(),
            nBlockSystem: ::std::default::Default::default(),
        }
    }
}
pub type inp_ATOM = tagInputAtom;
#[derive(Clone, Debug, PartialEq)]
pub struct OAD_PolymerUnit {
    pub id: i32,
    pub type_: i32,
    pub subtype: i32,
    pub conn: i32,
    pub label: i32,
    pub na: i32,
    pub nb: i32,
    pub cyclizable: i32,
    pub cyclized: i32,
    pub xbr1: [f64; 4usize],
    pub xbr2: [f64; 4usize],
    pub smt: [i8; 80usize],
    pub representation: i32,
    pub cap1: i32,
    pub end_atom1: i32,
    pub end_atom2: i32,
    pub cap2: i32,
    pub cap1_is_undef: i32,
    pub cap2_is_undef: i32,
    pub alist: SourceMutPointer<i32>,
    pub blist: SourceMutPointer<i32>,
    pub maxbkbonds: i32,
    pub nbkbonds: i32,
    pub bkbonds: SourceMutPointer<SourceMutPointer<i32>>,
}
impl ::std::default::Default for OAD_PolymerUnit {
    fn default() -> Self {
        Self {
            id: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            subtype: ::std::default::Default::default(),
            conn: ::std::default::Default::default(),
            label: ::std::default::Default::default(),
            na: ::std::default::Default::default(),
            nb: ::std::default::Default::default(),
            cyclizable: ::std::default::Default::default(),
            cyclized: ::std::default::Default::default(),
            xbr1: ::std::array::from_fn(|_| ::std::default::Default::default()),
            xbr2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            smt: ::std::array::from_fn(|_| ::std::default::Default::default()),
            representation: ::std::default::Default::default(),
            cap1: ::std::default::Default::default(),
            end_atom1: ::std::default::Default::default(),
            end_atom2: ::std::default::Default::default(),
            cap2: ::std::default::Default::default(),
            cap1_is_undef: ::std::default::Default::default(),
            cap2_is_undef: ::std::default::Default::default(),
            alist: ::std::default::Default::default(),
            blist: ::std::default::Default::default(),
            maxbkbonds: ::std::default::Default::default(),
            nbkbonds: ::std::default::Default::default(),
            bkbonds: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct OAD_Polymer {
    pub units: SourceMutPointer<SourceMutPointer<OAD_PolymerUnit>>,
    pub n: i32,
    pub n_pzz: i32,
    pub pzz: SourceMutPointer<i32>,
    pub really_do_frame_shift: i32,
    pub frame_shift_scheme: i32,
    pub treat: i32,
    pub representation: i32,
    pub is_in_reconn: i32,
    pub edit_repeats: i32,
}
impl ::std::default::Default for OAD_Polymer {
    fn default() -> Self {
        Self {
            units: ::std::default::Default::default(),
            n: ::std::default::Default::default(),
            n_pzz: ::std::default::Default::default(),
            pzz: ::std::default::Default::default(),
            really_do_frame_shift: ::std::default::Default::default(),
            frame_shift_scheme: ::std::default::Default::default(),
            treat: ::std::default::Default::default(),
            representation: ::std::default::Default::default(),
            is_in_reconn: ::std::default::Default::default(),
            edit_repeats: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct OAD_AtProps {
    pub erank: i32,
    pub ring_erank: i32,
    pub ring_num: i32,
    pub ring_size: i32,
}
impl ::std::default::Default for OAD_AtProps {
    fn default() -> Self {
        Self {
            erank: ::std::default::Default::default(),
            ring_erank: ::std::default::Default::default(),
            ring_num: ::std::default::Default::default(),
            ring_size: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct OAD_V3000 {
    pub n_non_star_atoms: i32,
    pub n_star_atoms: i32,
    pub atom_index_orig: SourceMutPointer<i32>,
    pub atom_index_fin: SourceMutPointer<i32>,
    pub n_sgroups: i32,
    pub n_3d_constraints: i32,
    pub n_collections: i32,
    pub n_non_haptic_bonds: i32,
    pub n_haptic_bonds: i32,
    pub lists_haptic_bonds: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_steabs: i32,
    pub lists_steabs: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_sterel: i32,
    pub lists_sterel: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_sterac: i32,
    pub lists_sterac: SourceMutPointer<SourceMutPointer<i32>>,
}
impl ::std::default::Default for OAD_V3000 {
    fn default() -> Self {
        Self {
            n_non_star_atoms: ::std::default::Default::default(),
            n_star_atoms: ::std::default::Default::default(),
            atom_index_orig: ::std::default::Default::default(),
            atom_index_fin: ::std::default::Default::default(),
            n_sgroups: ::std::default::Default::default(),
            n_3d_constraints: ::std::default::Default::default(),
            n_collections: ::std::default::Default::default(),
            n_non_haptic_bonds: ::std::default::Default::default(),
            n_haptic_bonds: ::std::default::Default::default(),
            lists_haptic_bonds: ::std::default::Default::default(),
            n_steabs: ::std::default::Default::default(),
            lists_steabs: ::std::default::Default::default(),
            n_sterel: ::std::default::Default::default(),
            lists_sterel: ::std::default::Default::default(),
            n_sterac: ::std::default::Default::default(),
            lists_sterac: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct tagOrigAtom {
    pub at: SourceMutPointer<inp_ATOM>,
    pub num_dimensions: i32,
    pub num_inp_bonds: i32,
    pub num_inp_atoms: i32,
    pub num_components: i32,
    pub bDisconnectSalts: i32,
    pub bDisconnectCoord: i32,
    pub nCurAtLen: SourceMutPointer<AT_NUMB>,
    pub nOldCompNumber: SourceMutPointer<AT_NUMB>,
    pub nNumEquSets: i32,
    pub nEquLabels: SourceMutPointer<AT_NUMB>,
    pub nSortedOrder: SourceMutPointer<AT_NUMB>,
    pub bSavedInINCHI_LIB: [i32; 2usize],
    pub bPreprocessed: [i32; 2usize],
    pub szCoord: SourceMutPointer<MOL_COORD>,
    pub polymer: SourceMutPointer<OAD_Polymer>,
    pub v3000: SourceMutPointer<OAD_V3000>,
    pub valid_polymer: i32,
    pub n_zy: i32,
}
impl ::std::default::Default for tagOrigAtom {
    fn default() -> Self {
        Self {
            at: ::std::default::Default::default(),
            num_dimensions: ::std::default::Default::default(),
            num_inp_bonds: ::std::default::Default::default(),
            num_inp_atoms: ::std::default::Default::default(),
            num_components: ::std::default::Default::default(),
            bDisconnectSalts: ::std::default::Default::default(),
            bDisconnectCoord: ::std::default::Default::default(),
            nCurAtLen: ::std::default::Default::default(),
            nOldCompNumber: ::std::default::Default::default(),
            nNumEquSets: ::std::default::Default::default(),
            nEquLabels: ::std::default::Default::default(),
            nSortedOrder: ::std::default::Default::default(),
            bSavedInINCHI_LIB: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bPreprocessed: ::std::array::from_fn(|_| ::std::default::Default::default()),
            szCoord: ::std::default::Default::default(),
            polymer: ::std::default::Default::default(),
            v3000: ::std::default::Default::default(),
            valid_polymer: ::std::default::Default::default(),
            n_zy: ::std::default::Default::default(),
        }
    }
}
pub type ORIG_ATOM_DATA = tagOrigAtom;
#[derive(Clone, Debug, PartialEq)]
pub struct tagOriginalStruct {
    pub num_atoms: i32,
    pub szAtoms: SourceMutPointer<i8>,
    pub szBonds: SourceMutPointer<i8>,
    pub szCoord: SourceMutPointer<i8>,
    pub polymer: SourceMutPointer<OAD_Polymer>,
    pub v3000: SourceMutPointer<OAD_V3000>,
    pub n_zy: i32,
}
impl ::std::default::Default for tagOriginalStruct {
    fn default() -> Self {
        Self {
            num_atoms: ::std::default::Default::default(),
            szAtoms: ::std::default::Default::default(),
            szBonds: ::std::default::Default::default(),
            szCoord: ::std::default::Default::default(),
            polymer: ::std::default::Default::default(),
            v3000: ::std::default::Default::default(),
            n_zy: ::std::default::Default::default(),
        }
    }
}
pub type ORIG_STRUCT = tagOriginalStruct;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomParmsForDrawing {
    pub at_string: [i8; 36usize],
    pub DrawingLabelLeftShift: i32,
    pub DrawingLabelLength: i32,
    pub nCanonNbr: AT_NUMB,
    pub nCanonEquNbr: AT_NUMB,
    pub nTautGroupCanonNbr: AT_NUMB,
    pub nTautGroupEquNbr: AT_NUMB,
    pub cFlags: S_CHAR,
    pub cHighlightTheAtom: S_CHAR,
    pub cStereoCenterParity: S_CHAR,
    pub cStereoBondParity: [S_CHAR; 3usize],
    pub cStereoBondWarning: [S_CHAR; 3usize],
    pub cStereoBondNumber: [S_CHAR; 3usize],
}
impl ::std::default::Default for tagAtomParmsForDrawing {
    fn default() -> Self {
        Self {
            at_string: ::std::array::from_fn(|_| ::std::default::Default::default()),
            DrawingLabelLeftShift: ::std::default::Default::default(),
            DrawingLabelLength: ::std::default::Default::default(),
            nCanonNbr: ::std::default::Default::default(),
            nCanonEquNbr: ::std::default::Default::default(),
            nTautGroupCanonNbr: ::std::default::Default::default(),
            nTautGroupEquNbr: ::std::default::Default::default(),
            cFlags: ::std::default::Default::default(),
            cHighlightTheAtom: ::std::default::Default::default(),
            cStereoCenterParity: ::std::default::Default::default(),
            cStereoBondParity: ::std::array::from_fn(|_| ::std::default::Default::default()),
            cStereoBondWarning: ::std::array::from_fn(|_| ::std::default::Default::default()),
            cStereoBondNumber: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type inf_ATOM = tagAtomParmsForDrawing;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInfoAtomData {
    pub at: SourceMutPointer<inf_ATOM>,
    pub num_at: i32,
    pub StereoFlags: AT_NUMB,
    pub num_components: AT_NUMB,
    pub pStereoFlags: SourceMutPointer<AT_NUMB>,
    pub nNumRemovedProtons: i32,
    pub num_removed_iso_H: i32,
    pub num_iso_H: [NUM_H; 3usize],
    pub szRemovedProtons: [i8; 128usize],
}
impl ::std::default::Default for tagInfoAtomData {
    fn default() -> Self {
        Self {
            at: ::std::default::Default::default(),
            num_at: ::std::default::Default::default(),
            StereoFlags: ::std::default::Default::default(),
            num_components: ::std::default::Default::default(),
            pStereoFlags: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            num_removed_iso_H: ::std::default::Default::default(),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            szRemovedProtons: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type INF_ATOM_DATA = tagInfoAtomData;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInputAtomData {
    pub at: SourceMutPointer<inp_ATOM>,
    pub at_fixed_bonds: SourceMutPointer<inp_ATOM>,
    pub num_at: i32,
    pub num_removed_H: i32,
    pub num_bonds: i32,
    pub num_isotopic: i32,
    pub bExists: i32,
    pub bDeleted: i32,
    pub bHasIsotopicLayer: i32,
    pub bTautomeric: i32,
    pub bTautPreprocessed: i32,
    pub nNumRemovedProtons: i32,
    pub nNumRemovedProtonsIsotopic: [NUM_H; 3usize],
    pub num_iso_H: [NUM_H; 3usize],
    pub bTautFlags: INCHI_MODE,
    pub bTautFlagsDone: INCHI_MODE,
    pub bNormalizationFlags: INCHI_MODE,
}
impl ::std::default::Default for tagInputAtomData {
    fn default() -> Self {
        Self {
            at: ::std::default::Default::default(),
            at_fixed_bonds: ::std::default::Default::default(),
            num_at: ::std::default::Default::default(),
            num_removed_H: ::std::default::Default::default(),
            num_bonds: ::std::default::Default::default(),
            num_isotopic: ::std::default::Default::default(),
            bExists: ::std::default::Default::default(),
            bDeleted: ::std::default::Default::default(),
            bHasIsotopicLayer: ::std::default::Default::default(),
            bTautomeric: ::std::default::Default::default(),
            bTautPreprocessed: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            nNumRemovedProtonsIsotopic: ::std::array::from_fn(|_| {
                ::std::default::Default::default()
            }),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bTautFlags: ::std::default::Default::default(),
            bTautFlagsDone: ::std::default::Default::default(),
            bNormalizationFlags: ::std::default::Default::default(),
        }
    }
}
pub type INP_ATOM_DATA = tagInputAtomData;
pub type INP_ATOM_DATA2 = [INP_ATOM_DATA; 2usize];
#[derive(Clone, Debug, PartialEq)]
pub struct tagNormCanonFlags {
    pub bTautFlags: [[INCHI_MODE; 2usize]; 2usize],
    pub bTautFlagsDone: [[INCHI_MODE; 2usize]; 2usize],
    pub bNormalizationFlags: [[INCHI_MODE; 2usize]; 2usize],
    pub nCanonFlags: [[i32; 2usize]; 2usize],
}
impl ::std::default::Default for tagNormCanonFlags {
    fn default() -> Self {
        Self {
            bTautFlags: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            bTautFlagsDone: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            bNormalizationFlags: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            nCanonFlags: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
        }
    }
}
pub type NORM_CANON_FLAGS = tagNormCanonFlags;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCompositeAtomData {
    pub at: SourceMutPointer<inp_ATOM>,
    pub num_at: i32,
    pub num_removed_H: i32,
    pub num_bonds: i32,
    pub num_isotopic: i32,
    pub bExists: i32,
    pub bDeleted: i32,
    pub bHasIsotopicLayer: i32,
    pub bTautomeric: i32,
    pub nNumRemovedProtons: i32,
    pub nNumRemovedProtonsIsotopic: [NUM_H; 3usize],
    pub num_iso_H: [NUM_H; 3usize],
    pub nOffsetAtAndH: SourceMutPointer<AT_NUMB>,
    pub num_components: i32,
}
impl ::std::default::Default for tagCompositeAtomData {
    fn default() -> Self {
        Self {
            at: ::std::default::Default::default(),
            num_at: ::std::default::Default::default(),
            num_removed_H: ::std::default::Default::default(),
            num_bonds: ::std::default::Default::default(),
            num_isotopic: ::std::default::Default::default(),
            bExists: ::std::default::Default::default(),
            bDeleted: ::std::default::Default::default(),
            bHasIsotopicLayer: ::std::default::Default::default(),
            bTautomeric: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            nNumRemovedProtonsIsotopic: ::std::array::from_fn(|_| {
                ::std::default::Default::default()
            }),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            nOffsetAtAndH: ::std::default::Default::default(),
            num_components: ::std::default::Default::default(),
        }
    }
}
pub type COMP_ATOM_DATA = tagCompositeAtomData;
pub type INCHI_FPTR = i64;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStructFptrs {
    pub fptr: SourceMutPointer<INCHI_FPTR>,
    pub len_fptr: i32,
    pub cur_fptr: i32,
    pub max_fptr: i32,
}
impl ::std::default::Default for tagStructFptrs {
    fn default() -> Self {
        Self {
            fptr: ::std::default::Default::default(),
            len_fptr: ::std::default::Default::default(),
            cur_fptr: ::std::default::Default::default(),
            max_fptr: ::std::default::Default::default(),
        }
    }
}
pub type STRUCT_FPTRS = tagStructFptrs;
#[derive(Clone, Debug, PartialEq)]
pub struct tagOAD_StructureEdits {
    pub del_atom: SourceMutPointer<INT_ARRAY>,
    pub del_bond: SourceMutPointer<INT_ARRAY>,
    pub new_bond: SourceMutPointer<INT_ARRAY>,
    pub mod_bond: SourceMutPointer<INT_ARRAY>,
    pub mod_coord: SourceMutPointer<INT_ARRAY>,
    pub del_side_chains: i32,
}
impl ::std::default::Default for tagOAD_StructureEdits {
    fn default() -> Self {
        Self {
            del_atom: ::std::default::Default::default(),
            del_bond: ::std::default::Default::default(),
            new_bond: ::std::default::Default::default(),
            mod_bond: ::std::default::Default::default(),
            mod_coord: ::std::default::Default::default(),
            del_side_chains: ::std::default::Default::default(),
        }
    }
}
pub type OAD_StructureEdits = tagOAD_StructureEdits;
#[derive(Clone, Debug, PartialEq)]
pub struct AtData {
    pub element: [i8; 3usize],
    pub maxvalence: i32,
}
impl ::std::default::Default for AtData {
    fn default() -> Self {
        Self {
            element: ::std::array::from_fn(|_| ::std::default::Default::default()),
            maxvalence: ::std::default::Default::default(),
        }
    }
}
pub type AT_ISO_SORT_KEY = i64;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStereoCarb {
    pub at_num: AT_NUMB,
    pub parity: U_CHAR,
}
impl ::std::default::Default for tagStereoCarb {
    fn default() -> Self {
        Self {
            at_num: ::std::default::Default::default(),
            parity: ::std::default::Default::default(),
        }
    }
}
pub type AT_STEREO_CARB = tagStereoCarb;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStereoDble {
    pub at_num1: AT_NUMB,
    pub at_num2: AT_NUMB,
    pub parity: U_CHAR,
}
impl ::std::default::Default for tagStereoDble {
    fn default() -> Self {
        Self {
            at_num1: ::std::default::Default::default(),
            at_num2: ::std::default::Default::default(),
            parity: ::std::default::Default::default(),
        }
    }
}
pub type AT_STEREO_DBLE = tagStereoDble;
#[derive(Clone, Debug, PartialEq)]
pub struct tagIsotopicAtom {
    pub at_num: AT_NUMB,
    pub num_1H: NUM_H,
    pub num_D: NUM_H,
    pub num_T: NUM_H,
    pub iso_atw_diff: NUM_H,
}
impl ::std::default::Default for tagIsotopicAtom {
    fn default() -> Self {
        Self {
            at_num: ::std::default::Default::default(),
            num_1H: ::std::default::Default::default(),
            num_D: ::std::default::Default::default(),
            num_T: ::std::default::Default::default(),
            iso_atw_diff: ::std::default::Default::default(),
        }
    }
}
pub type AT_ISOTOPIC = tagIsotopicAtom;
pub type AT_STEREO = AT_NUMB;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtom {
    pub elname: [i8; 6usize],
    pub neighbor: [AT_NUMB; 20usize],
    pub init_rank: AT_NUMB,
    pub orig_at_number: AT_NUMB,
    pub orig_compt_at_numb: AT_NUMB,
    pub bond_type: [U_CHAR; 20usize],
    pub el_number: U_CHAR,
    pub valence: S_CHAR,
    pub chem_bonds_valence: S_CHAR,
    pub num_H: S_CHAR,
    pub num_iso_H: [S_CHAR; 3usize],
    pub cFlags: S_CHAR,
    pub iso_atw_diff: S_CHAR,
    pub iso_sort_key: AT_ISO_SORT_KEY,
    pub charge: S_CHAR,
    pub radical: S_CHAR,
    pub marked: S_CHAR,
    pub endpoint: AT_NUMB,
    pub stereo_bond_neighbor: [AT_NUMB; 3usize],
    pub stereo_bond_neighbor2: [AT_NUMB; 3usize],
    pub stereo_bond_ord: [S_CHAR; 3usize],
    pub stereo_bond_ord2: [S_CHAR; 3usize],
    pub stereo_bond_z_prod: [S_CHAR; 3usize],
    pub stereo_bond_z_prod2: [S_CHAR; 3usize],
    pub stereo_bond_parity: [S_CHAR; 3usize],
    pub stereo_bond_parity2: [S_CHAR; 3usize],
    pub parity: S_CHAR,
    pub parity2: S_CHAR,
    pub stereo_atom_parity: S_CHAR,
    pub stereo_atom_parity2: S_CHAR,
    pub final_parity: S_CHAR,
    pub final_parity2: S_CHAR,
    pub bAmbiguousStereo: S_CHAR,
    pub bHasStereoOrEquToStereo: S_CHAR,
    pub bHasStereoOrEquToStereo2: S_CHAR,
    pub bCutVertex: S_CHAR,
    pub nRingSystem: AT_NUMB,
    pub nNumAtInRingSystem: AT_NUMB,
    pub nBlockSystem: AT_NUMB,
    pub z_dir: [S_CHAR; 3usize],
}
impl ::std::default::Default for tagAtom {
    fn default() -> Self {
        Self {
            elname: ::std::array::from_fn(|_| ::std::default::Default::default()),
            neighbor: ::std::array::from_fn(|_| ::std::default::Default::default()),
            init_rank: ::std::default::Default::default(),
            orig_at_number: ::std::default::Default::default(),
            orig_compt_at_numb: ::std::default::Default::default(),
            bond_type: ::std::array::from_fn(|_| ::std::default::Default::default()),
            el_number: ::std::default::Default::default(),
            valence: ::std::default::Default::default(),
            chem_bonds_valence: ::std::default::Default::default(),
            num_H: ::std::default::Default::default(),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            cFlags: ::std::default::Default::default(),
            iso_atw_diff: ::std::default::Default::default(),
            iso_sort_key: ::std::default::Default::default(),
            charge: ::std::default::Default::default(),
            radical: ::std::default::Default::default(),
            marked: ::std::default::Default::default(),
            endpoint: ::std::default::Default::default(),
            stereo_bond_neighbor: ::std::array::from_fn(|_| ::std::default::Default::default()),
            stereo_bond_neighbor2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            stereo_bond_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            stereo_bond_ord2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            stereo_bond_z_prod: ::std::array::from_fn(|_| ::std::default::Default::default()),
            stereo_bond_z_prod2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            stereo_bond_parity: ::std::array::from_fn(|_| ::std::default::Default::default()),
            stereo_bond_parity2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            parity: ::std::default::Default::default(),
            parity2: ::std::default::Default::default(),
            stereo_atom_parity: ::std::default::Default::default(),
            stereo_atom_parity2: ::std::default::Default::default(),
            final_parity: ::std::default::Default::default(),
            final_parity2: ::std::default::Default::default(),
            bAmbiguousStereo: ::std::default::Default::default(),
            bHasStereoOrEquToStereo: ::std::default::Default::default(),
            bHasStereoOrEquToStereo2: ::std::default::Default::default(),
            bCutVertex: ::std::default::Default::default(),
            nRingSystem: ::std::default::Default::default(),
            nNumAtInRingSystem: ::std::default::Default::default(),
            nBlockSystem: ::std::default::Default::default(),
            z_dir: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type sp_ATOM = tagAtom;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINChI_IsotopicAtom {
    pub nAtomNumber: AT_NUMB,
    pub nIsoDifference: NUM_H,
    pub nNum_H: NUM_H,
    pub nNum_D: NUM_H,
    pub nNum_T: NUM_H,
}
impl ::std::default::Default for tagINChI_IsotopicAtom {
    fn default() -> Self {
        Self {
            nAtomNumber: ::std::default::Default::default(),
            nIsoDifference: ::std::default::Default::default(),
            nNum_H: ::std::default::Default::default(),
            nNum_D: ::std::default::Default::default(),
            nNum_T: ::std::default::Default::default(),
        }
    }
}
pub type INChI_IsotopicAtom = tagINChI_IsotopicAtom;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINChI_IsotopicTGroup {
    pub nTGroupNumber: AT_NUMB,
    pub nNum_H: AT_NUMB,
    pub nNum_D: AT_NUMB,
    pub nNum_T: AT_NUMB,
}
impl ::std::default::Default for tagINChI_IsotopicTGroup {
    fn default() -> Self {
        Self {
            nTGroupNumber: ::std::default::Default::default(),
            nNum_H: ::std::default::Default::default(),
            nNum_D: ::std::default::Default::default(),
            nNum_T: ::std::default::Default::default(),
        }
    }
}
pub type INChI_IsotopicTGroup = tagINChI_IsotopicTGroup;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINChI_Stereo {
    pub nNumberOfStereoCenters: i32,
    pub nNumber: SourceMutPointer<AT_NUMB>,
    pub t_parity: SourceMutPointer<S_CHAR>,
    pub nNumberInv: SourceMutPointer<AT_NUMB>,
    pub t_parityInv: SourceMutPointer<S_CHAR>,
    pub nCompInv2Abs: i32,
    pub bTrivialInv: i32,
    pub nNumberOfStereoBonds: i32,
    pub nBondAtom1: SourceMutPointer<AT_NUMB>,
    pub nBondAtom2: SourceMutPointer<AT_NUMB>,
    pub b_parity: SourceMutPointer<S_CHAR>,
}
impl ::std::default::Default for tagINChI_Stereo {
    fn default() -> Self {
        Self {
            nNumberOfStereoCenters: ::std::default::Default::default(),
            nNumber: ::std::default::Default::default(),
            t_parity: ::std::default::Default::default(),
            nNumberInv: ::std::default::Default::default(),
            t_parityInv: ::std::default::Default::default(),
            nCompInv2Abs: ::std::default::Default::default(),
            bTrivialInv: ::std::default::Default::default(),
            nNumberOfStereoBonds: ::std::default::Default::default(),
            nBondAtom1: ::std::default::Default::default(),
            nBondAtom2: ::std::default::Default::default(),
            b_parity: ::std::default::Default::default(),
        }
    }
}
pub type INChI_Stereo = tagINChI_Stereo;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINChI {
    pub nErrorCode: i32,
    pub nFlags: INCHI_MODE,
    pub nTotalCharge: i32,
    pub nNumberOfAtoms: i32,
    pub szHillFormula: SourceMutPointer<i8>,
    pub nAtom: SourceMutPointer<U_CHAR>,
    pub lenConnTable: i32,
    pub nConnTable: SourceMutPointer<AT_NUMB>,
    pub lenTautomer: i32,
    pub nTautomer: SourceMutPointer<AT_NUMB>,
    pub nNum_H: SourceMutPointer<S_CHAR>,
    pub nNum_H_fixed: SourceMutPointer<S_CHAR>,
    pub nNumberOfIsotopicAtoms: i32,
    pub IsotopicAtom: SourceMutPointer<INChI_IsotopicAtom>,
    pub nNumberOfIsotopicTGroups: i32,
    pub IsotopicTGroup: SourceMutPointer<INChI_IsotopicTGroup>,
    pub Stereo: SourceMutPointer<INChI_Stereo>,
    pub StereoIsotopic: SourceMutPointer<INChI_Stereo>,
    pub nPossibleLocationsOfIsotopicH: SourceMutPointer<AT_NUMB>,
    pub bDeleted: i32,
    pub nRefCount: i32,
    pub nLink: i32,
}
impl ::std::default::Default for tagINChI {
    fn default() -> Self {
        Self {
            nErrorCode: ::std::default::Default::default(),
            nFlags: ::std::default::Default::default(),
            nTotalCharge: ::std::default::Default::default(),
            nNumberOfAtoms: ::std::default::Default::default(),
            szHillFormula: ::std::default::Default::default(),
            nAtom: ::std::default::Default::default(),
            lenConnTable: ::std::default::Default::default(),
            nConnTable: ::std::default::Default::default(),
            lenTautomer: ::std::default::Default::default(),
            nTautomer: ::std::default::Default::default(),
            nNum_H: ::std::default::Default::default(),
            nNum_H_fixed: ::std::default::Default::default(),
            nNumberOfIsotopicAtoms: ::std::default::Default::default(),
            IsotopicAtom: ::std::default::Default::default(),
            nNumberOfIsotopicTGroups: ::std::default::Default::default(),
            IsotopicTGroup: ::std::default::Default::default(),
            Stereo: ::std::default::Default::default(),
            StereoIsotopic: ::std::default::Default::default(),
            nPossibleLocationsOfIsotopicH: ::std::default::Default::default(),
            bDeleted: ::std::default::Default::default(),
            nRefCount: ::std::default::Default::default(),
            nLink: ::std::default::Default::default(),
        }
    }
}
pub type INChI = tagINChI;
pub type PINChI2 = [SourceMutPointer<INChI>; 2usize];
#[derive(Clone, Debug, PartialEq)]
pub struct tagOrigInfo {
    pub cCharge: S_CHAR,
    pub cRadical: S_CHAR,
    pub cUnusualValence: S_CHAR,
}
impl ::std::default::Default for tagOrigInfo {
    fn default() -> Self {
        Self {
            cCharge: ::std::default::Default::default(),
            cRadical: ::std::default::Default::default(),
            cUnusualValence: ::std::default::Default::default(),
        }
    }
}
pub type ORIG_INFO = tagOrigInfo;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINChI_Aux {
    pub nErrorCode: i32,
    pub nNumberOfAtoms: i32,
    pub nNumberOfTGroups: i32,
    pub bIsIsotopic: i32,
    pub bIsTautomeric: i32,
    pub nOrigAtNosInCanonOrd: SourceMutPointer<AT_NUMB>,
    pub nIsotopicOrigAtNosInCanonOrd: SourceMutPointer<AT_NUMB>,
    pub nOrigAtNosInCanonOrdInv: SourceMutPointer<AT_NUMB>,
    pub nIsotopicOrigAtNosInCanonOrdInv: SourceMutPointer<AT_NUMB>,
    pub nConstitEquNumbers: SourceMutPointer<AT_NUMB>,
    pub nConstitEquTGroupNumbers: SourceMutPointer<AT_NUMB>,
    pub nConstitEquIsotopicNumbers: SourceMutPointer<AT_NUMB>,
    pub nConstitEquIsotopicTGroupNumbers: SourceMutPointer<AT_NUMB>,
    pub nRefCount: i32,
    pub OrigInfo: SourceMutPointer<ORIG_INFO>,
    pub szOrigCoord: SourceMutPointer<MOL_COORD>,
    pub nNumRemovedProtons: NUM_H,
    pub nNumRemovedIsotopicH: [NUM_H; 3usize],
    pub bDeleted: i32,
    pub bTautFlags: INCHI_MODE,
    pub bTautFlagsDone: INCHI_MODE,
    pub bNormalizationFlags: INCHI_MODE,
    pub nCanonFlags: i32,
}
impl ::std::default::Default for tagINChI_Aux {
    fn default() -> Self {
        Self {
            nErrorCode: ::std::default::Default::default(),
            nNumberOfAtoms: ::std::default::Default::default(),
            nNumberOfTGroups: ::std::default::Default::default(),
            bIsIsotopic: ::std::default::Default::default(),
            bIsTautomeric: ::std::default::Default::default(),
            nOrigAtNosInCanonOrd: ::std::default::Default::default(),
            nIsotopicOrigAtNosInCanonOrd: ::std::default::Default::default(),
            nOrigAtNosInCanonOrdInv: ::std::default::Default::default(),
            nIsotopicOrigAtNosInCanonOrdInv: ::std::default::Default::default(),
            nConstitEquNumbers: ::std::default::Default::default(),
            nConstitEquTGroupNumbers: ::std::default::Default::default(),
            nConstitEquIsotopicNumbers: ::std::default::Default::default(),
            nConstitEquIsotopicTGroupNumbers: ::std::default::Default::default(),
            nRefCount: ::std::default::Default::default(),
            OrigInfo: ::std::default::Default::default(),
            szOrigCoord: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            nNumRemovedIsotopicH: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bDeleted: ::std::default::Default::default(),
            bTautFlags: ::std::default::Default::default(),
            bTautFlagsDone: ::std::default::Default::default(),
            bNormalizationFlags: ::std::default::Default::default(),
            nCanonFlags: ::std::default::Default::default(),
        }
    }
}
pub type INChI_Aux = tagINChI_Aux;
pub type PINChI_Aux2 = [SourceMutPointer<INChI_Aux>; 2usize];
#[derive(Clone, Debug, PartialEq)]
pub struct tagINChIforSort {
    pub pINChI: [SourceMutPointer<INChI>; 2usize],
    pub pINChI_Aux: [SourceMutPointer<INChI_Aux>; 2usize],
    pub ord_number: i16,
    pub n1: i16,
    pub n2: i16,
    pub n3: i16,
}
impl ::std::default::Default for tagINChIforSort {
    fn default() -> Self {
        Self {
            pINChI: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pINChI_Aux: ::std::array::from_fn(|_| ::std::default::Default::default()),
            ord_number: ::std::default::Default::default(),
            n1: ::std::default::Default::default(),
            n2: ::std::default::Default::default(),
            n3: ::std::default::Default::default(),
        }
    }
}
pub type INCHI_SORT = tagINChIforSort;
pub type Vertex = i32;
pub type EdgeIndex = i32;
pub type Edge = [i32; 2usize];
pub type BNS_IEDGE = i32;
pub type EdgeFlow = i32;
pub type VertexFlow = i32;
pub const tagAltPathConst_iALTP_MAX_LEN: tagAltPathConst = 0;
pub const tagAltPathConst_iALTP_FLOW: tagAltPathConst = 1;
pub const tagAltPathConst_iALTP_PATH_LEN: tagAltPathConst = 2;
pub const tagAltPathConst_iALTP_START_ATOM: tagAltPathConst = 3;
pub const tagAltPathConst_iALTP_END_ATOM: tagAltPathConst = 4;
pub const tagAltPathConst_iALTP_NEIGHBOR: tagAltPathConst = 5;
pub const tagAltPathConst_iALTP_HDR_LEN: tagAltPathConst = 5;
pub type tagAltPathConst = u32;
pub use self::tagAltPathConst as ALT_CONST;
#[derive(Clone, Debug, PartialEq)]
pub struct BnsEdge {
    pub neighbor1: AT_NUMB,
    pub neighbor12: AT_NUMB,
    pub neigh_ord: [AT_NUMB; 2usize],
    pub cap: EdgeFlow,
    pub cap0: EdgeFlow,
    pub flow: EdgeFlow,
    pub flow0: EdgeFlow,
    pub pass: S_CHAR,
    pub forbidden: S_CHAR,
}
impl ::std::default::Default for BnsEdge {
    fn default() -> Self {
        Self {
            neighbor1: ::std::default::Default::default(),
            neighbor12: ::std::default::Default::default(),
            neigh_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            cap: ::std::default::Default::default(),
            cap0: ::std::default::Default::default(),
            flow: ::std::default::Default::default(),
            flow0: ::std::default::Default::default(),
            pass: ::std::default::Default::default(),
            forbidden: ::std::default::Default::default(),
        }
    }
}
pub type BNS_EDGE = BnsEdge;
#[derive(Clone, Debug, PartialEq)]
pub struct BnsStEdge {
    pub cap: VertexFlow,
    pub cap0: VertexFlow,
    pub flow: VertexFlow,
    pub flow0: VertexFlow,
    pub pass: S_CHAR,
}
impl ::std::default::Default for BnsStEdge {
    fn default() -> Self {
        Self {
            cap: ::std::default::Default::default(),
            cap0: ::std::default::Default::default(),
            flow: ::std::default::Default::default(),
            flow0: ::std::default::Default::default(),
            pass: ::std::default::Default::default(),
        }
    }
}
pub type BNS_ST_EDGE = BnsStEdge;
#[derive(Clone, Debug, PartialEq)]
pub struct BnsVertex {
    pub st_edge: BNS_ST_EDGE,
    pub type_: AT_NUMB,
    pub num_adj_edges: AT_NUMB,
    pub max_adj_edges: AT_NUMB,
    pub iedge: SourceMutPointer<BNS_IEDGE>,
}
impl ::std::default::Default for BnsVertex {
    fn default() -> Self {
        Self {
            st_edge: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            num_adj_edges: ::std::default::Default::default(),
            max_adj_edges: ::std::default::Default::default(),
            iedge: ::std::default::Default::default(),
        }
    }
}
pub type BNS_VERTEX = BnsVertex;
pub type BNS_ALT_PATH = BnsAltPath;
#[derive(Clone, Debug, PartialEq)]
pub struct BalancedNetworkStructure {
    pub num_atoms: i32,
    pub num_added_atoms: i32,
    pub nMaxAddAtoms: i32,
    pub num_c_groups: i32,
    pub num_t_groups: i32,
    pub num_vertices: i32,
    pub num_bonds: i32,
    pub num_edges: i32,
    pub num_iedges: i32,
    pub num_added_edges: i32,
    pub nMaxAddEdges: i32,
    pub max_vertices: i32,
    pub max_edges: i32,
    pub max_iedges: i32,
    pub tot_st_cap: i32,
    pub tot_st_flow: i32,
    pub len_alt_path: i32,
    pub bNotASimplePath: i32,
    pub bChangeFlow: i32,
    pub vert: SourceMutPointer<BNS_VERTEX>,
    pub edge: SourceMutPointer<BNS_EDGE>,
    pub iedge: SourceMutPointer<BNS_IEDGE>,
    pub alt_path: SourceMutPointer<BNS_ALT_PATH>,
    pub altp: [SourceMutPointer<BNS_ALT_PATH>; 16usize],
    pub max_altp: i32,
    pub num_altp: i32,
    pub pbTautFlags: SourceMutPointer<INCHI_MODE>,
    pub pbTautFlagsDone: SourceMutPointer<INCHI_MODE>,
    pub type_TACN: AT_NUMB,
    pub type_T: AT_NUMB,
    pub type_CN: AT_NUMB,
    pub edge_forbidden_mask: S_CHAR,
    pub ic: SourceMutPointer<tagINCHI_CLOCK>,
    pub ulTimeOutTime: SourceMutPointer<tagInchiTime>,
}
impl ::std::default::Default for BalancedNetworkStructure {
    fn default() -> Self {
        Self {
            num_atoms: ::std::default::Default::default(),
            num_added_atoms: ::std::default::Default::default(),
            nMaxAddAtoms: ::std::default::Default::default(),
            num_c_groups: ::std::default::Default::default(),
            num_t_groups: ::std::default::Default::default(),
            num_vertices: ::std::default::Default::default(),
            num_bonds: ::std::default::Default::default(),
            num_edges: ::std::default::Default::default(),
            num_iedges: ::std::default::Default::default(),
            num_added_edges: ::std::default::Default::default(),
            nMaxAddEdges: ::std::default::Default::default(),
            max_vertices: ::std::default::Default::default(),
            max_edges: ::std::default::Default::default(),
            max_iedges: ::std::default::Default::default(),
            tot_st_cap: ::std::default::Default::default(),
            tot_st_flow: ::std::default::Default::default(),
            len_alt_path: ::std::default::Default::default(),
            bNotASimplePath: ::std::default::Default::default(),
            bChangeFlow: ::std::default::Default::default(),
            vert: ::std::default::Default::default(),
            edge: ::std::default::Default::default(),
            iedge: ::std::default::Default::default(),
            alt_path: ::std::default::Default::default(),
            altp: ::std::array::from_fn(|_| ::std::default::Default::default()),
            max_altp: ::std::default::Default::default(),
            num_altp: ::std::default::Default::default(),
            pbTautFlags: ::std::default::Default::default(),
            pbTautFlagsDone: ::std::default::Default::default(),
            type_TACN: ::std::default::Default::default(),
            type_T: ::std::default::Default::default(),
            type_CN: ::std::default::Default::default(),
            edge_forbidden_mask: ::std::default::Default::default(),
            ic: ::std::default::Default::default(),
            ulTimeOutTime: ::std::default::Default::default(),
        }
    }
}
pub type BN_STRUCT = BalancedNetworkStructure;
pub const tagBnsRadSrchMode_RAD_SRCH_NORM: tagBnsRadSrchMode = 0;
pub const tagBnsRadSrchMode_RAD_SRCH_FROM_FICT: tagBnsRadSrchMode = 1;
pub type tagBnsRadSrchMode = u32;
pub use self::tagBnsRadSrchMode as BRS_MODE;
#[derive(Clone, Debug, PartialEq)]
pub struct BalancedNetworkData {
    pub BasePtr: SourceMutPointer<Vertex>,
    pub SwitchEdge: SourceMutPointer<Edge>,
    pub Tree: SourceMutPointer<S_CHAR>,
    pub ScanQ: SourceMutPointer<Vertex>,
    pub QSize: i32,
    pub Pu: SourceMutPointer<Vertex>,
    pub Pv: SourceMutPointer<Vertex>,
    pub max_num_vertices: i32,
    pub max_len_Pu_Pv: i32,
    pub RadEndpoints: SourceMutPointer<Vertex>,
    pub nNumRadEndpoints: i32,
    pub RadEdges: SourceMutPointer<EdgeIndex>,
    pub nNumRadEdges: i32,
    pub nNumRadicals: i32,
    pub bRadSrchMode: BRS_MODE,
}
impl ::std::default::Default for BalancedNetworkData {
    fn default() -> Self {
        Self {
            BasePtr: ::std::default::Default::default(),
            SwitchEdge: ::std::default::Default::default(),
            Tree: ::std::default::Default::default(),
            ScanQ: ::std::default::Default::default(),
            QSize: ::std::default::Default::default(),
            Pu: ::std::default::Default::default(),
            Pv: ::std::default::Default::default(),
            max_num_vertices: ::std::default::Default::default(),
            max_len_Pu_Pv: ::std::default::Default::default(),
            RadEndpoints: ::std::default::Default::default(),
            nNumRadEndpoints: ::std::default::Default::default(),
            RadEdges: ::std::default::Default::default(),
            nNumRadEdges: ::std::default::Default::default(),
            nNumRadicals: ::std::default::Default::default(),
            bRadSrchMode: ::std::default::Default::default(),
        }
    }
}
pub type BN_DATA = BalancedNetworkData;
#[derive(Clone, Debug, PartialEq)]
pub struct BN_AtomsAtTautGroup {
    pub nAllocLen: i32,
    pub nNumFound: i32,
    pub nNumMainAdj2Tgroup: i32,
    pub nNumOthersAdj2Tgroup: i32,
    pub nEndPoint: SourceMutPointer<AT_NUMB>,
    pub nMarkedAtom: SourceMutPointer<S_CHAR>,
    pub nAtTypeTotals: SourceMutPointer<i32>,
    pub t_group_info: SourceMutPointer<tagTautomerGroupsInfo>,
}
impl ::std::default::Default for BN_AtomsAtTautGroup {
    fn default() -> Self {
        Self {
            nAllocLen: ::std::default::Default::default(),
            nNumFound: ::std::default::Default::default(),
            nNumMainAdj2Tgroup: ::std::default::Default::default(),
            nNumOthersAdj2Tgroup: ::std::default::Default::default(),
            nEndPoint: ::std::default::Default::default(),
            nMarkedAtom: ::std::default::Default::default(),
            nAtTypeTotals: ::std::default::Default::default(),
            t_group_info: ::std::default::Default::default(),
        }
    }
}
pub type BN_AATG = BN_AtomsAtTautGroup;
#[derive(Clone, Debug, PartialEq)]
pub struct tagBNS_FLOW_CHANGES {
    pub iedge: BNS_IEDGE,
    pub flow: EdgeFlow,
    pub cap: EdgeFlow,
    pub v1: Vertex,
    pub cap_st1: VertexFlow,
    pub flow_st1: VertexFlow,
    pub v2: Vertex,
    pub cap_st2: VertexFlow,
    pub flow_st2: VertexFlow,
}
impl ::std::default::Default for tagBNS_FLOW_CHANGES {
    fn default() -> Self {
        Self {
            iedge: ::std::default::Default::default(),
            flow: ::std::default::Default::default(),
            cap: ::std::default::Default::default(),
            v1: ::std::default::Default::default(),
            cap_st1: ::std::default::Default::default(),
            flow_st1: ::std::default::Default::default(),
            v2: ::std::default::Default::default(),
            cap_st2: ::std::default::Default::default(),
            flow_st2: ::std::default::Default::default(),
        }
    }
}
pub type BNS_FLOW_CHANGES = tagBNS_FLOW_CHANGES;
pub type bitWord = U_SHORT;
#[derive(Clone, Debug, PartialEq)]
pub struct tagNodeSet {
    pub bitword: SourceMutPointer<SourceMutPointer<bitWord>>,
    pub num_set: i32,
    pub len_set: i32,
}
impl ::std::default::Default for tagNodeSet {
    fn default() -> Self {
        Self {
            bitword: ::std::default::Default::default(),
            num_set: ::std::default::Default::default(),
            len_set: ::std::default::Default::default(),
        }
    }
}
pub type NodeSet = tagNodeSet;
pub type AT_TAUTOMER = AT_NUMB;
pub type T_GROUP_ISOWT = AT_ISO_SORT_KEY;
#[derive(Clone, Debug, PartialEq)]
pub struct tagIsotopicTautomerGroup {
    pub tgroup_num: AT_NUMB,
    pub num: [AT_NUMB; 3usize],
}
impl ::std::default::Default for tagIsotopicTautomerGroup {
    fn default() -> Self {
        Self {
            tgroup_num: ::std::default::Default::default(),
            num: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type AT_ISO_TGROUP = tagIsotopicTautomerGroup;
pub const tagTG_NumDA_TG_Num_dH: tagTG_NumDA = 0;
pub const tagTG_NumDA_TG_Num_dM: tagTG_NumDA = 1;
pub const tagTG_NumDA_TG_Num_aH: tagTG_NumDA = 2;
pub const tagTG_NumDA_TG_Num_aM: tagTG_NumDA = 3;
pub const tagTG_NumDA_TG_Num_dO: tagTG_NumDA = 4;
pub const tagTG_NumDA_TG_Num_aO: tagTG_NumDA = 5;
pub const tagTG_NumDA_TG_NUM_DA: tagTG_NumDA = 6;
pub type tagTG_NumDA = u32;
pub use self::tagTG_NumDA as TGNUMDA;
#[derive(Clone, Debug, PartialEq)]
pub struct tagTautomerGroup {
    pub num: [AT_RANK; 5usize],
    pub num_DA: [AT_RANK; 6usize],
    pub iWeight: T_GROUP_ISOWT,
    pub nGroupNumber: AT_NUMB,
    pub nNumEndpoints: AT_NUMB,
    pub nFirstEndpointAtNoPos: AT_NUMB,
}
impl ::std::default::Default for tagTautomerGroup {
    fn default() -> Self {
        Self {
            num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_DA: ::std::array::from_fn(|_| ::std::default::Default::default()),
            iWeight: ::std::default::Default::default(),
            nGroupNumber: ::std::default::Default::default(),
            nNumEndpoints: ::std::default::Default::default(),
            nFirstEndpointAtNoPos: ::std::default::Default::default(),
        }
    }
}
pub type T_GROUP = tagTautomerGroup;
#[derive(Clone, Debug, PartialEq)]
pub struct tagTautomerNormInfo {
    pub nNumRemovedExplicitH: NUM_H,
    pub nNumRemovedProtons: NUM_H,
    pub nNumRemovedProtonsIsotopic: [NUM_H; 3usize],
    pub bNormalizationFlags: INCHI_MODE,
}
impl ::std::default::Default for tagTautomerNormInfo {
    fn default() -> Self {
        Self {
            nNumRemovedExplicitH: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            nNumRemovedProtonsIsotopic: ::std::array::from_fn(|_| {
                ::std::default::Default::default()
            }),
            bNormalizationFlags: ::std::default::Default::default(),
        }
    }
}
pub type TNI = tagTautomerNormInfo;
#[derive(Clone, Debug, PartialEq)]
pub struct tagTautomerGroupsInfo {
    pub t_group: SourceMutPointer<T_GROUP>,
    pub nEndpointAtomNumber: SourceMutPointer<AT_NUMB>,
    pub tGroupNumber: SourceMutPointer<AT_NUMB>,
    pub nNumEndpoints: i32,
    pub num_t_groups: i32,
    pub max_num_t_groups: i32,
    pub bIgnoreIsotopic: i32,
    pub nIsotopicEndpointAtomNumber: SourceMutPointer<AT_NUMB>,
    pub nNumIsotopicEndpoints: i32,
    pub num_iso_H: [NUM_H; 3usize],
    pub tni: TNI,
    pub bTautFlags: INCHI_MODE,
    pub bTautFlagsDone: INCHI_MODE,
}
impl ::std::default::Default for tagTautomerGroupsInfo {
    fn default() -> Self {
        Self {
            t_group: ::std::default::Default::default(),
            nEndpointAtomNumber: ::std::default::Default::default(),
            tGroupNumber: ::std::default::Default::default(),
            nNumEndpoints: ::std::default::Default::default(),
            num_t_groups: ::std::default::Default::default(),
            max_num_t_groups: ::std::default::Default::default(),
            bIgnoreIsotopic: ::std::default::Default::default(),
            nIsotopicEndpointAtomNumber: ::std::default::Default::default(),
            nNumIsotopicEndpoints: ::std::default::Default::default(),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            tni: ::std::default::Default::default(),
            bTautFlags: ::std::default::Default::default(),
            bTautFlagsDone: ::std::default::Default::default(),
        }
    }
}
pub type T_GROUP_INFO = tagTautomerGroupsInfo;
#[derive(Clone, Debug, PartialEq)]
pub struct tagTautomerEndpoint {
    pub num: [AT_RANK; 5usize],
    pub num_DA: [AT_RANK; 6usize],
    pub nGroupNumber: AT_NUMB,
    pub nEquNumber: AT_NUMB,
    pub nAtomNumber: AT_NUMB,
}
impl ::std::default::Default for tagTautomerEndpoint {
    fn default() -> Self {
        Self {
            num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_DA: ::std::array::from_fn(|_| ::std::default::Default::default()),
            nGroupNumber: ::std::default::Default::default(),
            nEquNumber: ::std::default::Default::default(),
            nAtomNumber: ::std::default::Default::default(),
        }
    }
}
pub type T_ENDPOINT = tagTautomerEndpoint;
#[derive(Clone, Debug, PartialEq)]
pub struct tagTautomerBondLocation {
    pub nAtomNumber: AT_NUMB,
    pub neighbor_index: AT_NUMB,
}
impl ::std::default::Default for tagTautomerBondLocation {
    fn default() -> Self {
        Self {
            nAtomNumber: ::std::default::Default::default(),
            neighbor_index: ::std::default::Default::default(),
        }
    }
}
pub type T_BONDPOS = tagTautomerBondLocation;
#[derive(Clone, Debug, PartialEq)]
pub struct tagEndpointInfo {
    pub cMoveableCharge: S_CHAR,
    pub cNeutralBondsValence: S_CHAR,
    pub cMobile: S_CHAR,
    pub cDonor: S_CHAR,
    pub cAcceptor: S_CHAR,
    pub cKetoEnolCode: S_CHAR,
}
impl ::std::default::Default for tagEndpointInfo {
    fn default() -> Self {
        Self {
            cMoveableCharge: ::std::default::Default::default(),
            cNeutralBondsValence: ::std::default::Default::default(),
            cMobile: ::std::default::Default::default(),
            cDonor: ::std::default::Default::default(),
            cAcceptor: ::std::default::Default::default(),
            cKetoEnolCode: ::std::default::Default::default(),
        }
    }
}
pub type ENDPOINT_INFO = tagEndpointInfo;
#[derive(Clone, Debug, PartialEq)]
pub struct tagChargeCandidate {
    pub atnumber: AT_NUMB,
    pub type_: S_CHAR,
    pub subtype: S_CHAR,
}
impl ::std::default::Default for tagChargeCandidate {
    fn default() -> Self {
        Self {
            atnumber: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            subtype: ::std::default::Default::default(),
        }
    }
}
pub type C_CANDIDATE = tagChargeCandidate;
#[derive(Clone, Debug, PartialEq)]
pub struct tagChargeGroup {
    pub num: [AT_RANK; 2usize],
    pub num_CPoints: AT_RANK,
    pub nGroupNumber: AT_NUMB,
    pub cGroupType: U_CHAR,
}
impl ::std::default::Default for tagChargeGroup {
    fn default() -> Self {
        Self {
            num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_CPoints: ::std::default::Default::default(),
            nGroupNumber: ::std::default::Default::default(),
            cGroupType: ::std::default::Default::default(),
        }
    }
}
pub type C_GROUP = tagChargeGroup;
#[derive(Clone, Debug, PartialEq)]
pub struct tagChargeGroupsInfo {
    pub c_group: SourceMutPointer<C_GROUP>,
    pub num_c_groups: i32,
    pub max_num_c_groups: i32,
    pub c_candidate: SourceMutPointer<C_CANDIDATE>,
    pub max_num_candidates: i32,
    pub num_candidates: i32,
}
impl ::std::default::Default for tagChargeGroupsInfo {
    fn default() -> Self {
        Self {
            c_group: ::std::default::Default::default(),
            num_c_groups: ::std::default::Default::default(),
            max_num_c_groups: ::std::default::Default::default(),
            c_candidate: ::std::default::Default::default(),
            max_num_candidates: ::std::default::Default::default(),
            num_candidates: ::std::default::Default::default(),
        }
    }
}
pub type C_GROUP_INFO = tagChargeGroupsInfo;
#[derive(Clone, Debug, PartialEq)]
pub struct tagSaltChargeCandidate {
    pub atnumber: AT_NUMB,
    pub type_: S_CHAR,
    pub subtype: S_CHAR,
    pub endpoint: AT_NUMB,
}
impl ::std::default::Default for tagSaltChargeCandidate {
    fn default() -> Self {
        Self {
            atnumber: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            subtype: ::std::default::Default::default(),
            endpoint: ::std::default::Default::default(),
        }
    }
}
pub type S_CANDIDATE = tagSaltChargeCandidate;
#[derive(Clone, Debug, PartialEq)]
pub struct tagSaltGroupInfo {
    pub s_candidate: SourceMutPointer<S_CANDIDATE>,
    pub max_num_candidates: i32,
    pub num_candidates: i32,
    pub num_other_candidates: i32,
    pub num_p_only_candidates: i32,
}
impl ::std::default::Default for tagSaltGroupInfo {
    fn default() -> Self {
        Self {
            s_candidate: ::std::default::Default::default(),
            max_num_candidates: ::std::default::Default::default(),
            num_candidates: ::std::default::Default::default(),
            num_other_candidates: ::std::default::Default::default(),
            num_p_only_candidates: ::std::default::Default::default(),
        }
    }
}
pub type S_GROUP_INFO = tagSaltGroupInfo;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomSizes {
    pub nMaxNumStereoAtoms: i32,
    pub nMaxNumStereoBonds: i32,
    pub num_isotopic_atoms: i32,
    pub nLenCT: i32,
    pub nLenBonds: i32,
    pub nLenIsotopic: i32,
    pub nLenCTAtOnly: i32,
    pub nLenLinearCTStereoDble: i32,
    pub nLenLinearCTStereoCarb: i32,
    pub bMayHaveStereo: i32,
    pub bIgnoreIsotopic: i32,
    pub nLenLinearCTTautomer: i32,
    pub nLenLinearCTIsotopicTautomer: i32,
    pub bHasIsotopicTautGroups: i32,
    pub nLenIsotopicEndpoints: i32,
}
impl ::std::default::Default for tagAtomSizes {
    fn default() -> Self {
        Self {
            nMaxNumStereoAtoms: ::std::default::Default::default(),
            nMaxNumStereoBonds: ::std::default::Default::default(),
            num_isotopic_atoms: ::std::default::Default::default(),
            nLenCT: ::std::default::Default::default(),
            nLenBonds: ::std::default::Default::default(),
            nLenIsotopic: ::std::default::Default::default(),
            nLenCTAtOnly: ::std::default::Default::default(),
            nLenLinearCTStereoDble: ::std::default::Default::default(),
            nLenLinearCTStereoCarb: ::std::default::Default::default(),
            bMayHaveStereo: ::std::default::Default::default(),
            bIgnoreIsotopic: ::std::default::Default::default(),
            nLenLinearCTTautomer: ::std::default::Default::default(),
            nLenLinearCTIsotopicTautomer: ::std::default::Default::default(),
            bHasIsotopicTautGroups: ::std::default::Default::default(),
            nLenIsotopicEndpoints: ::std::default::Default::default(),
        }
    }
}
pub type ATOM_SIZES = tagAtomSizes;
#[derive(Clone, Debug, PartialEq)]
pub struct tagDfsPath {
    pub at_no: AT_RANK,
    pub bond_type: U_CHAR,
    pub bond_pos: S_CHAR,
}
impl ::std::default::Default for tagDfsPath {
    fn default() -> Self {
        Self {
            at_no: ::std::default::Default::default(),
            bond_type: ::std::default::Default::default(),
            bond_pos: ::std::default::Default::default(),
        }
    }
}
pub type DFS_PATH = tagDfsPath;
pub type SU_LONG = tagSplitLong;
pub type NEIGH_LIST = SourceMutPointer<AT_RANK>;
#[derive(Clone, Debug, PartialEq)]
pub struct tagEQUIV_INFO {
    pub nNumSets: i32,
    pub nCutVertexAtom: SourceMutPointer<i32>,
    pub nFirstInSet: SourceMutPointer<i32>,
    pub nNumInSet: SourceMutPointer<i32>,
    pub nAtomNo: SourceMutPointer<i32>,
    pub nAddToRank: SourceMutPointer<i32>,
}
impl ::std::default::Default for tagEQUIV_INFO {
    fn default() -> Self {
        Self {
            nNumSets: ::std::default::Default::default(),
            nCutVertexAtom: ::std::default::Default::default(),
            nFirstInSet: ::std::default::Default::default(),
            nNumInSet: ::std::default::Default::default(),
            nAtomNo: ::std::default::Default::default(),
            nAddToRank: ::std::default::Default::default(),
        }
    }
}
pub type EQUIV_INFO = tagEQUIV_INFO;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtData_dch {
    pub element: [i8; 3usize],
    pub valence: i32,
}
impl ::std::default::Default for tagAtData_dch {
    fn default() -> Self {
        Self {
            element: ::std::array::from_fn(|_| ::std::default::Default::default()),
            valence: ::std::default::Default::default(),
        }
    }
}
pub type AT_DATA = tagAtData_dch;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomInvariantBytes {
    pub cNotExactlyHillOrderNumber: S_CHAR,
    pub cNumberOfConnections: S_CHAR,
    pub cAtomicNumber: S_CHAR,
    pub cNumberOfAttachedHydrogens: S_CHAR,
}
impl ::std::default::Default for tagAtomInvariantBytes {
    fn default() -> Self {
        Self {
            cNotExactlyHillOrderNumber: ::std::default::Default::default(),
            cNumberOfConnections: ::std::default::Default::default(),
            cAtomicNumber: ::std::default::Default::default(),
            cNumberOfAttachedHydrogens: ::std::default::Default::default(),
        }
    }
}
pub type ATOM_INVARIANT_BYTES = tagAtomInvariantBytes;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomInvariant {
    pub b: ATOM_INVARIANT_BYTES,
    pub cNum_tautomer: AT_RANK,
    pub cNum_tautomer_num: [AT_RANK; 2usize],
    pub iso_sort_key: AT_ISO_SORT_KEY,
    pub cNum_tautomer_iso: [AT_RANK; 3usize],
}
impl ::std::default::Default for tagAtomInvariant {
    fn default() -> Self {
        Self {
            b: ::std::default::Default::default(),
            cNum_tautomer: ::std::default::Default::default(),
            cNum_tautomer_num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            iso_sort_key: ::std::default::Default::default(),
            cNum_tautomer_iso: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type ATOM_INVARIANT = tagAtomInvariant;
pub const tagAtInvariantIndexes_AT_INV_HILL_ORDER: tagAtInvariantIndexes = 0;
pub const tagAtInvariantIndexes_AT_INV_NUM_CONNECTIONS: tagAtInvariantIndexes = 1;
pub const tagAtInvariantIndexes_AT_INV_NUM_H: tagAtInvariantIndexes = 2;
pub const tagAtInvariantIndexes_AT_INV_NUM_TG_ENDPOINTS: tagAtInvariantIndexes = 3;
pub const tagAtInvariantIndexes_AT_INV_TG_NUMBERS: tagAtInvariantIndexes = 4;
pub const tagAtInvariantIndexes_AT_INV_NUM_H_FIX: tagAtInvariantIndexes = 6;
pub const tagAtInvariantIndexes_AT_INV_BREAK1: tagAtInvariantIndexes = 7;
pub const tagAtInvariantIndexes_AT_INV_TAUT_ISO: tagAtInvariantIndexes = 7;
pub const tagAtInvariantIndexes_AT_INV_LENGTH: tagAtInvariantIndexes = 10;
pub type tagAtInvariantIndexes = u32;
pub use self::tagAtInvariantIndexes as AT_INV_INDEXES;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomInvariant2 {
    pub val: [AT_NUMB; 10usize],
    pub iso_sort_key: AT_ISO_SORT_KEY,
    pub iso_aux_key: S_CHAR,
}
impl ::std::default::Default for tagAtomInvariant2 {
    fn default() -> Self {
        Self {
            val: ::std::array::from_fn(|_| ::std::default::Default::default()),
            iso_sort_key: ::std::default::Default::default(),
            iso_aux_key: ::std::default::Default::default(),
        }
    }
}
pub type ATOM_INVARIANT2 = tagAtomInvariant2;
#[derive(Clone, Debug, PartialEq)]
pub struct tagPartition {
    pub Rank: SourceMutPointer<AT_RANK>,
    pub AtNumber: SourceMutPointer<AT_NUMB>,
}
impl ::std::default::Default for tagPartition {
    fn default() -> Self {
        Self {
            Rank: ::std::default::Default::default(),
            AtNumber: ::std::default::Default::default(),
        }
    }
}
pub type Partition = tagPartition;
#[derive(Clone, Debug, PartialEq)]
pub struct tagFixHOrTautCanonNumbering {
    pub num_at_tg: i32,
    pub num_atoms: i32,
    pub nCanonFlags: i32,
    pub NeighList: SourceMutPointer<NEIGH_LIST>,
    pub LinearCt: SourceMutPointer<AT_RANK>,
    pub nLenLinearCtAtOnly: i32,
    pub nLenLinearCt: i32,
    pub nMaxLenLinearCt: i32,
    pub PartitionCt: Partition,
    pub nSymmRankCt: SourceMutPointer<AT_RANK>,
    pub nNumHOrig: SourceMutPointer<NUM_H>,
    pub nNumH: SourceMutPointer<NUM_H>,
    pub nLenNumH: i32,
    pub nNumHOrigFixH: SourceMutPointer<NUM_H>,
    pub nNumHFixH: SourceMutPointer<NUM_H>,
    pub nLenNumHFixH: i32,
    pub PartitionCtIso: Partition,
    pub nSymmRankCtIso: SourceMutPointer<AT_RANK>,
    pub iso_sort_keys: SourceMutPointer<AT_ISO_SORT_KEY>,
    pub iso_sort_keysOrig: SourceMutPointer<AT_ISO_SORT_KEY>,
    pub len_iso_sort_keys: i32,
    pub iso_exchg_atnos: SourceMutPointer<S_CHAR>,
    pub iso_exchg_atnosOrig: SourceMutPointer<S_CHAR>,
}
impl ::std::default::Default for tagFixHOrTautCanonNumbering {
    fn default() -> Self {
        Self {
            num_at_tg: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            nCanonFlags: ::std::default::Default::default(),
            NeighList: ::std::default::Default::default(),
            LinearCt: ::std::default::Default::default(),
            nLenLinearCtAtOnly: ::std::default::Default::default(),
            nLenLinearCt: ::std::default::Default::default(),
            nMaxLenLinearCt: ::std::default::Default::default(),
            PartitionCt: ::std::default::Default::default(),
            nSymmRankCt: ::std::default::Default::default(),
            nNumHOrig: ::std::default::Default::default(),
            nNumH: ::std::default::Default::default(),
            nLenNumH: ::std::default::Default::default(),
            nNumHOrigFixH: ::std::default::Default::default(),
            nNumHFixH: ::std::default::Default::default(),
            nLenNumHFixH: ::std::default::Default::default(),
            PartitionCtIso: ::std::default::Default::default(),
            nSymmRankCtIso: ::std::default::Default::default(),
            iso_sort_keys: ::std::default::Default::default(),
            iso_sort_keysOrig: ::std::default::Default::default(),
            len_iso_sort_keys: ::std::default::Default::default(),
            iso_exchg_atnos: ::std::default::Default::default(),
            iso_exchg_atnosOrig: ::std::default::Default::default(),
        }
    }
}
pub type FTCN = tagFixHOrTautCanonNumbering;
#[derive(Clone, Debug, PartialEq)]
pub struct tagBaseCanonNumbering {
    pub pRankStack: SourceMutPointer<SourceMutPointer<AT_RANK>>,
    pub nMaxLenRankStack: i32,
    pub num_max: i32,
    pub num_at_tg: i32,
    pub num_atoms: i32,
    pub ulTimeOutTime: SourceMutPointer<tagInchiTime>,
    pub ftcn: [FTCN; 2usize],
}
impl ::std::default::Default for tagBaseCanonNumbering {
    fn default() -> Self {
        Self {
            pRankStack: ::std::default::Default::default(),
            nMaxLenRankStack: ::std::default::Default::default(),
            num_max: ::std::default::Default::default(),
            num_at_tg: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            ulTimeOutTime: ::std::default::Default::default(),
            ftcn: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type BCN = tagBaseCanonNumbering;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCanonStat {
    pub lNumBreakTies: i64,
    pub lNumNeighListIter: i64,
    pub lNumTotCT: i64,
    pub lNumDecreasedCT: i64,
    pub lNumRejectedCT: i64,
    pub lNumEqualCT: i64,
    pub ulTimeOutTime: SourceMutPointer<tagInchiTime>,
    pub lTotalTime: i64,
    pub bFirstCT: i32,
    pub bKeepSymmRank: i32,
    pub bStereoIsBetter: i32,
    pub nCanonFlags: i32,
    pub LinearCT: SourceMutPointer<AT_NUMB>,
    pub LinearCTIsotopic: SourceMutPointer<AT_ISOTOPIC>,
    pub LinearCTIsotopicTautomer: SourceMutPointer<AT_ISO_TGROUP>,
    pub LinearCTStereoDble: SourceMutPointer<AT_STEREO_DBLE>,
    pub LinearCTStereoCarb: SourceMutPointer<AT_STEREO_CARB>,
    pub LinearCTStereoDbleInv: SourceMutPointer<AT_STEREO_DBLE>,
    pub LinearCTStereoCarbInv: SourceMutPointer<AT_STEREO_CARB>,
    pub LinearCTIsotopicStereoDble: SourceMutPointer<AT_STEREO_DBLE>,
    pub LinearCTIsotopicStereoCarb: SourceMutPointer<AT_STEREO_CARB>,
    pub LinearCTIsotopicStereoDbleInv: SourceMutPointer<AT_STEREO_DBLE>,
    pub LinearCTIsotopicStereoCarbInv: SourceMutPointer<AT_STEREO_CARB>,
    pub LinearCTTautomer: SourceMutPointer<AT_TAUTOMER>,
    pub LinearCT2: SourceMutPointer<AT_NUMB>,
    pub nLenLinearCTStereoDble: i32,
    pub nLenLinearCTStereoDbleInv: i32,
    pub nMaxLenLinearCTStereoDble: i32,
    pub bCmpStereo: i32,
    pub nLenLinearCTStereoCarb: i32,
    pub nLenLinearCTStereoCarbInv: i32,
    pub nMaxLenLinearCTStereoCarb: i32,
    pub nLenLinearCTIsotopic: i32,
    pub nMaxLenLinearCTIsotopic: i32,
    pub nLenLinearCTIsotopicTautomer: i32,
    pub nMaxLenLinearCTIsotopicTautomer: i32,
    pub nLenLinearCT: i32,
    pub nLenLinearCT2: i32,
    pub nLenLinearCTAtOnly: i32,
    pub nLenLinearCTAtOnly2: i32,
    pub nMaxLenLinearCT: i32,
    pub nLenLinearCTTautomer: i32,
    pub nMaxLenLinearCTTautomer: i32,
    pub bCmpIsotopicStereo: i32,
    pub nLenLinearCTIsotopicStereoDble: i32,
    pub nLenLinearCTIsotopicStereoDbleInv: i32,
    pub nMaxLenLinearCTIsotopicStereoDble: i32,
    pub nLenLinearCTIsotopicStereoCarb: i32,
    pub nLenLinearCTIsotopicStereoCarbInv: i32,
    pub nMaxLenLinearCTIsotopicStereoCarb: i32,
    pub bRankUsedForStereo: SourceMutPointer<S_CHAR>,
    pub bAtomUsedForStereo: SourceMutPointer<S_CHAR>,
    pub nPrevAtomNumber: SourceMutPointer<AT_RANK>,
    pub nCanonOrd: SourceMutPointer<AT_RANK>,
    pub nSymmRank: SourceMutPointer<AT_RANK>,
    pub nCanonOrdTaut: SourceMutPointer<AT_RANK>,
    pub nSymmRankTaut: SourceMutPointer<AT_RANK>,
    pub nCanonOrdStereo: SourceMutPointer<AT_RANK>,
    pub nCanonOrdStereoInv: SourceMutPointer<AT_RANK>,
    pub nCanonOrdStereoTaut: SourceMutPointer<AT_RANK>,
    pub nSymmRankIsotopic: SourceMutPointer<AT_RANK>,
    pub nCanonOrdIsotopic: SourceMutPointer<AT_RANK>,
    pub nSymmRankIsotopicTaut: SourceMutPointer<AT_RANK>,
    pub nCanonOrdIsotopicTaut: SourceMutPointer<AT_RANK>,
    pub nCanonOrdIsotopicStereo: SourceMutPointer<AT_RANK>,
    pub nCanonOrdIsotopicStereoInv: SourceMutPointer<AT_RANK>,
    pub nCanonOrdIsotopicStereoTaut: SourceMutPointer<AT_RANK>,
    pub nLenCanonOrd: i32,
    pub nLenCanonOrdTaut: i32,
    pub nLenCanonOrdIsotopic: i32,
    pub nLenCanonOrdIsotopicTaut: i32,
    pub nLenCanonOrdStereo: i32,
    pub nLenCanonOrdStereoTaut: i32,
    pub nLenCanonOrdIsotopicStereo: i32,
    pub nLenCanonOrdIsotopicStereoTaut: i32,
    pub bHasIsotopicInTautomerGroups: i32,
    pub t_group_info: SourceMutPointer<T_GROUP_INFO>,
    pub bIgnoreIsotopic: i32,
    pub bDoubleBondSquare: i32,
    pub nMode: INCHI_MODE,
    pub NeighList: SourceMutPointer<NEIGH_LIST>,
    pub pBCN: SourceMutPointer<BCN>,
    pub nNum_H: SourceMutPointer<S_CHAR>,
    pub nNum_H_fixed: SourceMutPointer<S_CHAR>,
    pub nExchgIsoH: SourceMutPointer<S_CHAR>,
}
impl ::std::default::Default for tagCanonStat {
    fn default() -> Self {
        Self {
            lNumBreakTies: ::std::default::Default::default(),
            lNumNeighListIter: ::std::default::Default::default(),
            lNumTotCT: ::std::default::Default::default(),
            lNumDecreasedCT: ::std::default::Default::default(),
            lNumRejectedCT: ::std::default::Default::default(),
            lNumEqualCT: ::std::default::Default::default(),
            ulTimeOutTime: ::std::default::Default::default(),
            lTotalTime: ::std::default::Default::default(),
            bFirstCT: ::std::default::Default::default(),
            bKeepSymmRank: ::std::default::Default::default(),
            bStereoIsBetter: ::std::default::Default::default(),
            nCanonFlags: ::std::default::Default::default(),
            LinearCT: ::std::default::Default::default(),
            LinearCTIsotopic: ::std::default::Default::default(),
            LinearCTIsotopicTautomer: ::std::default::Default::default(),
            LinearCTStereoDble: ::std::default::Default::default(),
            LinearCTStereoCarb: ::std::default::Default::default(),
            LinearCTStereoDbleInv: ::std::default::Default::default(),
            LinearCTStereoCarbInv: ::std::default::Default::default(),
            LinearCTIsotopicStereoDble: ::std::default::Default::default(),
            LinearCTIsotopicStereoCarb: ::std::default::Default::default(),
            LinearCTIsotopicStereoDbleInv: ::std::default::Default::default(),
            LinearCTIsotopicStereoCarbInv: ::std::default::Default::default(),
            LinearCTTautomer: ::std::default::Default::default(),
            LinearCT2: ::std::default::Default::default(),
            nLenLinearCTStereoDble: ::std::default::Default::default(),
            nLenLinearCTStereoDbleInv: ::std::default::Default::default(),
            nMaxLenLinearCTStereoDble: ::std::default::Default::default(),
            bCmpStereo: ::std::default::Default::default(),
            nLenLinearCTStereoCarb: ::std::default::Default::default(),
            nLenLinearCTStereoCarbInv: ::std::default::Default::default(),
            nMaxLenLinearCTStereoCarb: ::std::default::Default::default(),
            nLenLinearCTIsotopic: ::std::default::Default::default(),
            nMaxLenLinearCTIsotopic: ::std::default::Default::default(),
            nLenLinearCTIsotopicTautomer: ::std::default::Default::default(),
            nMaxLenLinearCTIsotopicTautomer: ::std::default::Default::default(),
            nLenLinearCT: ::std::default::Default::default(),
            nLenLinearCT2: ::std::default::Default::default(),
            nLenLinearCTAtOnly: ::std::default::Default::default(),
            nLenLinearCTAtOnly2: ::std::default::Default::default(),
            nMaxLenLinearCT: ::std::default::Default::default(),
            nLenLinearCTTautomer: ::std::default::Default::default(),
            nMaxLenLinearCTTautomer: ::std::default::Default::default(),
            bCmpIsotopicStereo: ::std::default::Default::default(),
            nLenLinearCTIsotopicStereoDble: ::std::default::Default::default(),
            nLenLinearCTIsotopicStereoDbleInv: ::std::default::Default::default(),
            nMaxLenLinearCTIsotopicStereoDble: ::std::default::Default::default(),
            nLenLinearCTIsotopicStereoCarb: ::std::default::Default::default(),
            nLenLinearCTIsotopicStereoCarbInv: ::std::default::Default::default(),
            nMaxLenLinearCTIsotopicStereoCarb: ::std::default::Default::default(),
            bRankUsedForStereo: ::std::default::Default::default(),
            bAtomUsedForStereo: ::std::default::Default::default(),
            nPrevAtomNumber: ::std::default::Default::default(),
            nCanonOrd: ::std::default::Default::default(),
            nSymmRank: ::std::default::Default::default(),
            nCanonOrdTaut: ::std::default::Default::default(),
            nSymmRankTaut: ::std::default::Default::default(),
            nCanonOrdStereo: ::std::default::Default::default(),
            nCanonOrdStereoInv: ::std::default::Default::default(),
            nCanonOrdStereoTaut: ::std::default::Default::default(),
            nSymmRankIsotopic: ::std::default::Default::default(),
            nCanonOrdIsotopic: ::std::default::Default::default(),
            nSymmRankIsotopicTaut: ::std::default::Default::default(),
            nCanonOrdIsotopicTaut: ::std::default::Default::default(),
            nCanonOrdIsotopicStereo: ::std::default::Default::default(),
            nCanonOrdIsotopicStereoInv: ::std::default::Default::default(),
            nCanonOrdIsotopicStereoTaut: ::std::default::Default::default(),
            nLenCanonOrd: ::std::default::Default::default(),
            nLenCanonOrdTaut: ::std::default::Default::default(),
            nLenCanonOrdIsotopic: ::std::default::Default::default(),
            nLenCanonOrdIsotopicTaut: ::std::default::Default::default(),
            nLenCanonOrdStereo: ::std::default::Default::default(),
            nLenCanonOrdStereoTaut: ::std::default::Default::default(),
            nLenCanonOrdIsotopicStereo: ::std::default::Default::default(),
            nLenCanonOrdIsotopicStereoTaut: ::std::default::Default::default(),
            bHasIsotopicInTautomerGroups: ::std::default::Default::default(),
            t_group_info: ::std::default::Default::default(),
            bIgnoreIsotopic: ::std::default::Default::default(),
            bDoubleBondSquare: ::std::default::Default::default(),
            nMode: ::std::default::Default::default(),
            NeighList: ::std::default::Default::default(),
            pBCN: ::std::default::Default::default(),
            nNum_H: ::std::default::Default::default(),
            nNum_H_fixed: ::std::default::Default::default(),
            nExchgIsoH: ::std::default::Default::default(),
        }
    }
}
pub type CANON_STAT = tagCanonStat;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCANON_GLOBALS {
    pub m_pNeighList_RankForSort: SourceConstPointer<NEIGH_LIST>,
    pub m_pAtomInvariant2ForSort: SourceConstPointer<ATOM_INVARIANT2>,
    pub m_pNeighborsForSort: SourceConstPointer<AT_NUMB>,
    pub m_pn_RankForSort: SourceConstPointer<AT_RANK>,
    pub m_nMaxAtNeighRankForSort: AT_RANK,
    pub m_nNumCompNeighborsRanksCountEql: i32,
    pub m_bBit: SourceMutPointer<bitWord>,
    pub m_bBitInitialized: i32,
    pub m_num_bit: i32,
}
impl ::std::default::Default for tagCANON_GLOBALS {
    fn default() -> Self {
        Self {
            m_pNeighList_RankForSort: ::std::default::Default::default(),
            m_pAtomInvariant2ForSort: ::std::default::Default::default(),
            m_pNeighborsForSort: ::std::default::Default::default(),
            m_pn_RankForSort: ::std::default::Default::default(),
            m_nMaxAtNeighRankForSort: ::std::default::Default::default(),
            m_nNumCompNeighborsRanksCountEql: ::std::default::Default::default(),
            m_bBit: ::std::default::Default::default(),
            m_bBitInitialized: ::std::default::Default::default(),
            m_num_bit: ::std::default::Default::default(),
        }
    }
}
pub type CANON_GLOBALS = tagCANON_GLOBALS;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCanonData {
    pub LinearCT: SourceMutPointer<AT_NUMB>,
    pub nMaxLenLinearCT: i32,
    pub nLenLinearCT: i32,
    pub nLenCTAtOnly: i32,
    pub nCanonFlags: i32,
    pub NumH: SourceMutPointer<NUM_H>,
    pub lenNumH: i32,
    pub maxlenNumH: i32,
    pub NumHfixed: SourceMutPointer<NUM_H>,
    pub lenNumHfixed: i32,
    pub maxlenNumHfixed: i32,
    pub iso_sort_key: SourceMutPointer<AT_ISO_SORT_KEY>,
    pub len_iso_sort_key: i32,
    pub maxlen_iso_sort_key: i32,
    pub iso_exchg_atnos: SourceMutPointer<S_CHAR>,
    pub len_iso_exchg_atnos: i32,
    pub maxlen_iso_exchg_atnos: i32,
    pub nAuxRank: SourceMutPointer<AT_RANK>,
    pub ulTimeOutTime: SourceMutPointer<tagInchiTime>,
}
impl ::std::default::Default for tagCanonData {
    fn default() -> Self {
        Self {
            LinearCT: ::std::default::Default::default(),
            nMaxLenLinearCT: ::std::default::Default::default(),
            nLenLinearCT: ::std::default::Default::default(),
            nLenCTAtOnly: ::std::default::Default::default(),
            nCanonFlags: ::std::default::Default::default(),
            NumH: ::std::default::Default::default(),
            lenNumH: ::std::default::Default::default(),
            maxlenNumH: ::std::default::Default::default(),
            NumHfixed: ::std::default::Default::default(),
            lenNumHfixed: ::std::default::Default::default(),
            maxlenNumHfixed: ::std::default::Default::default(),
            iso_sort_key: ::std::default::Default::default(),
            len_iso_sort_key: ::std::default::Default::default(),
            maxlen_iso_sort_key: ::std::default::Default::default(),
            iso_exchg_atnos: ::std::default::Default::default(),
            len_iso_exchg_atnos: ::std::default::Default::default(),
            maxlen_iso_exchg_atnos: ::std::default::Default::default(),
            nAuxRank: ::std::default::Default::default(),
            ulTimeOutTime: ::std::default::Default::default(),
        }
    }
}
pub type CANON_DATA = tagCanonData;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCanonCounts {
    pub lNumBreakTies: i64,
    pub lNumDecreasedCT: i64,
    pub lNumRejectedCT: i64,
    pub lNumEqualCT: i64,
    pub lNumTotCT: i64,
    pub dGroupSize: f64,
    pub lNumGenerators: i64,
    pub lNumStoredIsomorphisms: i64,
}
impl ::std::default::Default for tagCanonCounts {
    fn default() -> Self {
        Self {
            lNumBreakTies: ::std::default::Default::default(),
            lNumDecreasedCT: ::std::default::Default::default(),
            lNumRejectedCT: ::std::default::Default::default(),
            lNumEqualCT: ::std::default::Default::default(),
            lNumTotCT: ::std::default::Default::default(),
            dGroupSize: ::std::default::Default::default(),
            lNumGenerators: ::std::default::Default::default(),
            lNumStoredIsomorphisms: ::std::default::Default::default(),
        }
    }
}
pub type CANON_COUNTS = tagCanonCounts;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCurTree {
    pub tree: SourceMutPointer<AT_NUMB>,
    pub max_len: i32,
    pub cur_len: i32,
    pub incr_len: i32,
}
impl ::std::default::Default for tagCurTree {
    fn default() -> Self {
        Self {
            tree: ::std::default::Default::default(),
            max_len: ::std::default::Default::default(),
            cur_len: ::std::default::Default::default(),
            incr_len: ::std::default::Default::default(),
        }
    }
}
pub type CUR_TREE = tagCurTree;
#[derive(Clone, Debug, PartialEq)]
pub struct tagEquNeigh {
    pub num_to: i32,
    pub to_at: [AT_RANK; 4usize],
    pub from_at: AT_RANK,
    pub rank: AT_RANK,
    pub canon_rank: AT_RANK,
}
impl ::std::default::Default for tagEquNeigh {
    fn default() -> Self {
        Self {
            num_to: ::std::default::Default::default(),
            to_at: ::std::array::from_fn(|_| ::std::default::Default::default()),
            from_at: ::std::default::Default::default(),
            rank: ::std::default::Default::default(),
            canon_rank: ::std::default::Default::default(),
        }
    }
}
pub type EQ_NEIGH = tagEquNeigh;
pub const tagInchiDiffBits_IDIF_PROBLEM: tagInchiDiffBits = 1;
pub const tagInchiDiffBits_IDIF_NUM_AT: tagInchiDiffBits = 1;
pub const tagInchiDiffBits_IDIF_ATOMS: tagInchiDiffBits = 1;
pub const tagInchiDiffBits_IDIF_NUM_EL: tagInchiDiffBits = 1;
pub const tagInchiDiffBits_IDIF_CON_LEN: tagInchiDiffBits = 1;
pub const tagInchiDiffBits_IDIF_CON_TBL: tagInchiDiffBits = 1;
pub const tagInchiDiffBits_IDIF_POSITION_H: tagInchiDiffBits = 2;
pub const tagInchiDiffBits_IDIF_MORE_FH: tagInchiDiffBits = 4;
pub const tagInchiDiffBits_IDIF_LESS_FH: tagInchiDiffBits = 8;
pub const tagInchiDiffBits_IDIF_MORE_H: tagInchiDiffBits = 16;
pub const tagInchiDiffBits_IDIF_LESS_H: tagInchiDiffBits = 32;
pub const tagInchiDiffBits_IDIF_NO_TAUT: tagInchiDiffBits = 64;
pub const tagInchiDiffBits_IDIF_WRONG_TAUT: tagInchiDiffBits = 128;
pub const tagInchiDiffBits_IDIF_SINGLE_TG: tagInchiDiffBits = 256;
pub const tagInchiDiffBits_IDIF_MULTIPLE_TG: tagInchiDiffBits = 512;
pub const tagInchiDiffBits_IDIF_NUM_TG: tagInchiDiffBits = 1024;
pub const tagInchiDiffBits_IDIF_EXTRA_TG_ENDP: tagInchiDiffBits = 2048;
pub const tagInchiDiffBits_IDIF_MISS_TG_ENDP: tagInchiDiffBits = 4096;
pub const tagInchiDiffBits_IDIF_DIFF_TG_ENDP: tagInchiDiffBits = 8192;
pub const tagInchiDiffBits_IDIF_TG: tagInchiDiffBits = 16384;
pub const tagInchiDiffBits_IDIF_NUM_ISO_AT: tagInchiDiffBits = 32768;
pub const tagInchiDiffBits_IDIF_ISO_AT: tagInchiDiffBits = 65536;
pub const tagInchiDiffBits_IDIF_CHARGE: tagInchiDiffBits = 131072;
pub const tagInchiDiffBits_IDIF_REM_PROT: tagInchiDiffBits = 262144;
pub const tagInchiDiffBits_IDIF_REM_ISO_H: tagInchiDiffBits = 524288;
pub const tagInchiDiffBits_IDIF_SC_INV: tagInchiDiffBits = 1048576;
pub const tagInchiDiffBits_IDIF_SC_PARITY: tagInchiDiffBits = 2097152;
pub const tagInchiDiffBits_IDIF_SC_EXTRA_UNDF: tagInchiDiffBits = 4194304;
pub const tagInchiDiffBits_IDIF_SC_EXTRA: tagInchiDiffBits = 8388608;
pub const tagInchiDiffBits_IDIF_SC_MISS_UNDF: tagInchiDiffBits = 16777216;
pub const tagInchiDiffBits_IDIF_SC_MISS: tagInchiDiffBits = 33554432;
pub const tagInchiDiffBits_IDIF_SB_PARITY: tagInchiDiffBits = 67108864;
pub const tagInchiDiffBits_IDIF_SB_EXTRA_UNDF: tagInchiDiffBits = 134217728;
pub const tagInchiDiffBits_IDIF_SB_EXTRA: tagInchiDiffBits = 268435456;
pub const tagInchiDiffBits_IDIF_SB_MISS_UNDF: tagInchiDiffBits = 536870912;
pub const tagInchiDiffBits_IDIF_SB_MISS: tagInchiDiffBits = 1073741824;
pub type tagInchiDiffBits = u32;
pub use self::tagInchiDiffBits as IDIF;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInChICompareResults {
    pub flags: INCHI_MODE,
    pub tot_num_H1: i32,
    pub tot_num_H2: i32,
    pub num_taut_H1: i32,
    pub num_taut_H2: i32,
    pub num_taut_M1: i32,
    pub num_taut_M2: i32,
    pub endp_in1_only: [AT_NUMB; 32usize],
    pub num_endp_in1_only: i32,
    pub endp_in2_only: [AT_NUMB; 32usize],
    pub num_endp_in2_only: i32,
    pub diff_pos_H_at: [AT_NUMB; 32usize],
    pub diff_pos_H_nH: [S_CHAR; 32usize],
    pub num_diff_pos_H: i32,
    pub fixed_H_at1_more: [AT_NUMB; 32usize],
    pub fixed_H_nH1_more: [S_CHAR; 32usize],
    pub num_fixed_H1_more: i32,
    pub fixed_H_at2_more: [AT_NUMB; 32usize],
    pub fixed_H_nH2_more: [S_CHAR; 32usize],
    pub num_fixed_H2_more: i32,
    pub sc_in1_only: [AT_NUMB; 32usize],
    pub num_sc_in1_only: i32,
    pub sc_in2_only: [AT_NUMB; 32usize],
    pub num_sc_in2_only: i32,
    pub sb_in1_only: [AT_NUMB; 32usize],
    pub num_sb_in1_only: i32,
    pub sb_in2_only: [AT_NUMB; 32usize],
    pub num_sb_in2_only: i32,
    pub sb_undef_in1_only: [AT_NUMB; 32usize],
    pub num_sb_undef_in1_only: i32,
    pub sb_undef_in2_only: [AT_NUMB; 32usize],
    pub num_sb_undef_in2_only: i32,
    pub sc_undef_in1_only: [AT_NUMB; 32usize],
    pub num_sc_undef_in1_only: i32,
    pub sc_undef_in2_only: [AT_NUMB; 32usize],
    pub num_sc_undef_in2_only: i32,
}
impl ::std::default::Default for tagInChICompareResults {
    fn default() -> Self {
        Self {
            flags: ::std::default::Default::default(),
            tot_num_H1: ::std::default::Default::default(),
            tot_num_H2: ::std::default::Default::default(),
            num_taut_H1: ::std::default::Default::default(),
            num_taut_H2: ::std::default::Default::default(),
            num_taut_M1: ::std::default::Default::default(),
            num_taut_M2: ::std::default::Default::default(),
            endp_in1_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_endp_in1_only: ::std::default::Default::default(),
            endp_in2_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_endp_in2_only: ::std::default::Default::default(),
            diff_pos_H_at: ::std::array::from_fn(|_| ::std::default::Default::default()),
            diff_pos_H_nH: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_diff_pos_H: ::std::default::Default::default(),
            fixed_H_at1_more: ::std::array::from_fn(|_| ::std::default::Default::default()),
            fixed_H_nH1_more: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_fixed_H1_more: ::std::default::Default::default(),
            fixed_H_at2_more: ::std::array::from_fn(|_| ::std::default::Default::default()),
            fixed_H_nH2_more: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_fixed_H2_more: ::std::default::Default::default(),
            sc_in1_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sc_in1_only: ::std::default::Default::default(),
            sc_in2_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sc_in2_only: ::std::default::Default::default(),
            sb_in1_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sb_in1_only: ::std::default::Default::default(),
            sb_in2_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sb_in2_only: ::std::default::Default::default(),
            sb_undef_in1_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sb_undef_in1_only: ::std::default::Default::default(),
            sb_undef_in2_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sb_undef_in2_only: ::std::default::Default::default(),
            sc_undef_in1_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sc_undef_in1_only: ::std::default::Default::default(),
            sc_undef_in2_only: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_sc_undef_in2_only: ::std::default::Default::default(),
        }
    }
}
pub type ICR = tagInChICompareResults;
pub const tagDiffINChISegments_DIFS_f_FORMULA: tagDiffINChISegments = 0;
pub const tagDiffINChISegments_DIFS_c_CONNECT: tagDiffINChISegments = 1;
pub const tagDiffINChISegments_DIFS_h_H_ATOMS: tagDiffINChISegments = 2;
pub const tagDiffINChISegments_DIFS_q_CHARGE: tagDiffINChISegments = 3;
pub const tagDiffINChISegments_DIFS_p_PROTONS: tagDiffINChISegments = 4;
pub const tagDiffINChISegments_DIFS_b_SBONDS: tagDiffINChISegments = 5;
pub const tagDiffINChISegments_DIFS_t_SATOMS: tagDiffINChISegments = 6;
pub const tagDiffINChISegments_DIFS_m_SP3INV: tagDiffINChISegments = 7;
pub const tagDiffINChISegments_DIFS_s_STYPE: tagDiffINChISegments = 8;
pub const tagDiffINChISegments_DIFS_i_IATOMS: tagDiffINChISegments = 9;
pub const tagDiffINChISegments_DIFS_o_TRANSP: tagDiffINChISegments = 10;
pub const tagDiffINChISegments_DIFS_idf_LENGTH: tagDiffINChISegments = 11;
pub const tagDiffINChISegments_DIFS_LENGTH: tagDiffINChISegments = 11;
pub type tagDiffINChISegments = u32;
pub use self::tagDiffINChISegments as DIF_SEGMENTS;
pub const tagDiffINChILayers_DIFL_M: tagDiffINChILayers = 0;
pub const tagDiffINChILayers_DIFL_MI: tagDiffINChILayers = 1;
pub const tagDiffINChILayers_DIFL_F: tagDiffINChILayers = 2;
pub const tagDiffINChILayers_DIFL_FI: tagDiffINChILayers = 3;
pub const tagDiffINChILayers_DIFL_LENGTH: tagDiffINChILayers = 4;
pub type tagDiffINChILayers = u32;
pub use self::tagDiffINChILayers as DIF_LAYERS;
pub const tagMarkDiff_DIFV_BOTH_EMPTY: tagMarkDiff = 0;
pub const tagMarkDiff_DIFV_EQL2PRECED: tagMarkDiff = 1;
pub const tagMarkDiff_DIFV_NEQ2PRECED: tagMarkDiff = 2;
pub const tagMarkDiff_DIFV_IS_EMPTY: tagMarkDiff = 4;
pub const tagMarkDiff_DIFV_FI_EQ_MI: tagMarkDiff = 8;
pub const tagMarkDiff_DIFV_OUTPUT_EMPTY_T: tagMarkDiff = 4;
pub const tagMarkDiff_DIFV_OUTPUT_EMPTY_F: tagMarkDiff = 11;
pub const tagMarkDiff_DIFV_OUTPUT_OMIT_F: tagMarkDiff = 6;
pub const tagMarkDiff_DIFV_OUTPUT_FILL_T: tagMarkDiff = 11;
pub type tagMarkDiff = u32;
pub use self::tagMarkDiff as DIF_VALUES;
pub const tagINChISegmAction_INCHI_SEGM_OMIT: tagINChISegmAction = 0;
pub const tagINChISegmAction_INCHI_SEGM_FILL: tagINChISegmAction = 1;
pub const tagINChISegmAction_INCHI_SEGM_EMPTY: tagINChISegmAction = 2;
pub type tagINChISegmAction = u32;
pub use self::tagINChISegmAction as INCHI_SEGM_ACTION;
#[derive(Clone, Debug, PartialEq)]
pub struct subgraf_edge {
    pub nbr: i32,
    pub etype: i32,
}
impl ::std::default::Default for subgraf_edge {
    fn default() -> Self {
        Self {
            nbr: ::std::default::Default::default(),
            etype: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct subgraf {
    pub nnodes: i32,
    pub nodes: SourceMutPointer<i32>,
    pub degrees: SourceMutPointer<i32>,
    pub orig2node: SourceMutPointer<i32>,
    pub adj: SourceMutPointer<SourceMutPointer<subgraf_edge>>,
}
impl ::std::default::Default for subgraf {
    fn default() -> Self {
        Self {
            nnodes: ::std::default::Default::default(),
            nodes: ::std::default::Default::default(),
            degrees: ::std::default::Default::default(),
            orig2node: ::std::default::Default::default(),
            adj: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct subgraf_pathfinder {
    pub sg: SourceMutPointer<subgraf>,
    pub start: i32,
    pub end: i32,
    pub maxbonds: i32,
    pub nbonds: i32,
    pub nseen: i32,
    pub seen: SourceMutPointer<i32>,
}
impl ::std::default::Default for subgraf_pathfinder {
    fn default() -> Self {
        Self {
            sg: ::std::default::Default::default(),
            start: ::std::default::Default::default(),
            end: ::std::default::Default::default(),
            maxbonds: ::std::default::Default::default(),
            nbonds: ::std::default::Default::default(),
            nseen: ::std::default::Default::default(),
            seen: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct tagLine {
    pub str_: SourceMutPointer<i8>,
    pub len: i32,
    pub len_alloc: i32,
    pub c: i32,
}
impl ::std::default::Default for tagLine {
    fn default() -> Self {
        Self {
            str_: ::std::default::Default::default(),
            len: ::std::default::Default::default(),
            len_alloc: ::std::default::Default::default(),
            c: ::std::default::Default::default(),
        }
    }
}
pub type SEGM_LINE = tagLine;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStructData {
    pub ulStructTime: u64,
    pub nErrorCode: i32,
    pub nErrorType: i32,
    pub nStructReadError: i32,
    pub pStrErrStruct: [i8; 256usize],
    pub fPtrStart: i64,
    pub fPtrEnd: i64,
    pub bUserQuit: i32,
    pub bUserQuitComponent: i32,
    pub bUserQuitComponentDisplay: i32,
    pub bChiralFlag: i32,
    pub num_taut: [i32; 2usize],
    pub num_non_taut: [i32; 2usize],
    pub bTautFlags: [INCHI_MODE; 2usize],
    pub bTautFlagsDone: [INCHI_MODE; 2usize],
    pub num_components: [i32; 2usize],
}
impl ::std::default::Default for tagStructData {
    fn default() -> Self {
        Self {
            ulStructTime: ::std::default::Default::default(),
            nErrorCode: ::std::default::Default::default(),
            nErrorType: ::std::default::Default::default(),
            nStructReadError: ::std::default::Default::default(),
            pStrErrStruct: ::std::array::from_fn(|_| ::std::default::Default::default()),
            fPtrStart: ::std::default::Default::default(),
            fPtrEnd: ::std::default::Default::default(),
            bUserQuit: ::std::default::Default::default(),
            bUserQuitComponent: ::std::default::Default::default(),
            bUserQuitComponentDisplay: ::std::default::Default::default(),
            bChiralFlag: ::std::default::Default::default(),
            num_taut: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_non_taut: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bTautFlags: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bTautFlagsDone: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_components: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type STRUCT_DATA = tagStructData;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHI_OUT_CTL {
    pub ATOM_MODE: i32,
    pub TAUT_MODE: i32,
    pub pSortPrintINChIFlags: SourceMutPointer<i32>,
    pub bOverflow: i32,
    pub bAlways: i32,
    pub bOutputType: i32,
    pub bOutType: i32,
    pub bPlainTextTags: i32,
    pub bOmitRepetitions: i32,
    pub bUseMulipliers: i32,
    pub bNonTautNonIsoIdentifierNotEmpty: i32,
    pub bNonTautIsoIdentifierNotEmpty: i32,
    pub bSecondNonTautPass: i32,
    pub bTautomericOutputAllowed: i32,
    pub bTautomeric: i32,
    pub bNonTautomeric: i32,
    pub bNonTautIsIdenticalToTaut: i32,
    pub bFhTag: i32,
    pub bRelRac: i32,
    pub bAbcNumbers: i32,
    pub bIsotopic: i32,
    pub bPolymers: i32,
    pub iCurTautMode: i32,
    pub num_components: i32,
    pub nNumRemovedProtons: i32,
    pub nTag: i32,
    pub bTag1: i32,
    pub bTag2: i32,
    pub bTag3: i32,
    pub tot_len: i32,
    pub tot_len2: i32,
    pub nCurINChISegment: i32,
    pub nSegmAction: i32,
    pub num_comp: [i32; 2usize],
    pub num_iso_H: [i32; 3usize],
    pub bAtomEqu: [i32; 2usize],
    pub bTautEqu: [i32; 2usize],
    pub bInvStereo: [i32; 2usize],
    pub bInvStereoOrigNumb: [i32; 2usize],
    pub bRacemicStereo: [i32; 2usize],
    pub bRelativeStereo: [i32; 2usize],
    pub bIsotopicOrigNumb: [i32; 2usize],
    pub bIsotopicAtomEqu: [i32; 2usize],
    pub bIsotopicTautEqu: [i32; 2usize],
    pub bInvIsotopicStereo: [i32; 2usize],
    pub bInvIsotopicStereoOrigNumb: [i32; 2usize],
    pub bIsotopicRacemicStereo: [i32; 2usize],
    pub bIsotopicRelativeStereo: [i32; 2usize],
    pub bIgn_UU_Sp3: [i32; 2usize],
    pub bIgn_UU_Sp2: [i32; 2usize],
    pub bIgn_UU_Sp3_Iso: [i32; 2usize],
    pub bIgn_UU_Sp2_Iso: [i32; 2usize],
    pub bChargesRadVal: [i32; 2usize],
    pub bOrigCoord: [i32; 2usize],
    pub sDifSegs: [[i8; 11usize]; 4usize],
    pub szTag1: [i8; 64usize],
    pub szTag2: [i8; 64usize],
    pub szTag3: [i8; 64usize],
    pub n_pzz: i32,
    pub n_zy: i32,
    pub pINChISortTautAndNonTaut: SourceMutPointer<SourceMutPointer<INCHI_SORT>>,
    pub pINChISort: SourceMutPointer<INCHI_SORT>,
    pub pINChISort2: SourceMutPointer<INCHI_SORT>,
}
impl ::std::default::Default for tagINCHI_OUT_CTL {
    fn default() -> Self {
        Self {
            ATOM_MODE: ::std::default::Default::default(),
            TAUT_MODE: ::std::default::Default::default(),
            pSortPrintINChIFlags: ::std::default::Default::default(),
            bOverflow: ::std::default::Default::default(),
            bAlways: ::std::default::Default::default(),
            bOutputType: ::std::default::Default::default(),
            bOutType: ::std::default::Default::default(),
            bPlainTextTags: ::std::default::Default::default(),
            bOmitRepetitions: ::std::default::Default::default(),
            bUseMulipliers: ::std::default::Default::default(),
            bNonTautNonIsoIdentifierNotEmpty: ::std::default::Default::default(),
            bNonTautIsoIdentifierNotEmpty: ::std::default::Default::default(),
            bSecondNonTautPass: ::std::default::Default::default(),
            bTautomericOutputAllowed: ::std::default::Default::default(),
            bTautomeric: ::std::default::Default::default(),
            bNonTautomeric: ::std::default::Default::default(),
            bNonTautIsIdenticalToTaut: ::std::default::Default::default(),
            bFhTag: ::std::default::Default::default(),
            bRelRac: ::std::default::Default::default(),
            bAbcNumbers: ::std::default::Default::default(),
            bIsotopic: ::std::default::Default::default(),
            bPolymers: ::std::default::Default::default(),
            iCurTautMode: ::std::default::Default::default(),
            num_components: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            nTag: ::std::default::Default::default(),
            bTag1: ::std::default::Default::default(),
            bTag2: ::std::default::Default::default(),
            bTag3: ::std::default::Default::default(),
            tot_len: ::std::default::Default::default(),
            tot_len2: ::std::default::Default::default(),
            nCurINChISegment: ::std::default::Default::default(),
            nSegmAction: ::std::default::Default::default(),
            num_comp: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bAtomEqu: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bTautEqu: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bInvStereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bInvStereoOrigNumb: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bRacemicStereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bRelativeStereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIsotopicOrigNumb: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIsotopicAtomEqu: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIsotopicTautEqu: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bInvIsotopicStereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bInvIsotopicStereoOrigNumb: ::std::array::from_fn(|_| {
                ::std::default::Default::default()
            }),
            bIsotopicRacemicStereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIsotopicRelativeStereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIgn_UU_Sp3: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIgn_UU_Sp2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIgn_UU_Sp3_Iso: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bIgn_UU_Sp2_Iso: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bChargesRadVal: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bOrigCoord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sDifSegs: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            szTag1: ::std::array::from_fn(|_| ::std::default::Default::default()),
            szTag2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            szTag3: ::std::array::from_fn(|_| ::std::default::Default::default()),
            n_pzz: ::std::default::Default::default(),
            n_zy: ::std::default::Default::default(),
            pINChISortTautAndNonTaut: ::std::default::Default::default(),
            pINChISort: ::std::default::Default::default(),
            pINChISort2: ::std::default::Default::default(),
        }
    }
}
pub type INCHI_OUT_CTL = tagINCHI_OUT_CTL;
#[derive(Clone, Debug, PartialEq)]
pub struct tagPOSEContext {
    pub sd: STRUCT_DATA,
    pub ip: INPUT_PARMS,
    pub szTitle: [i8; 575usize],
    pub pINChI2: [SourceMutPointer<PINChI2>; 2usize],
    pub pINChI_Aux2: [SourceMutPointer<PINChI_Aux2>; 2usize],
    pub inp_file: SourceMutPointer<INCHI_IOSTREAM>,
    pub inchi_file: [INCHI_IOSTREAM; 3usize],
    pub log_file: SourceMutPointer<INCHI_IOSTREAM>,
    pub out_file: SourceMutPointer<INCHI_IOSTREAM>,
    pub prb_file: SourceMutPointer<INCHI_IOSTREAM>,
    pub OrigAtData: ORIG_ATOM_DATA,
    pub orig_inp_data: SourceMutPointer<ORIG_ATOM_DATA>,
    pub PrepAtData: [ORIG_ATOM_DATA; 2usize],
    pub prep_inp_data: SourceMutPointer<ORIG_ATOM_DATA>,
    pub num_inp: i64,
    pub temp_string_container: INCHI_IOS_STRING,
    pub strbuf: SourceMutPointer<INCHI_IOS_STRING>,
    pub save_opt_bits: u8,
}
impl ::std::default::Default for tagPOSEContext {
    fn default() -> Self {
        Self {
            sd: ::std::default::Default::default(),
            ip: ::std::default::Default::default(),
            szTitle: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pINChI2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pINChI_Aux2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            inp_file: ::std::default::Default::default(),
            inchi_file: ::std::array::from_fn(|_| ::std::default::Default::default()),
            log_file: ::std::default::Default::default(),
            out_file: ::std::default::Default::default(),
            prb_file: ::std::default::Default::default(),
            OrigAtData: ::std::default::Default::default(),
            orig_inp_data: ::std::default::Default::default(),
            PrepAtData: ::std::array::from_fn(|_| ::std::default::Default::default()),
            prep_inp_data: ::std::default::Default::default(),
            num_inp: ::std::default::Default::default(),
            temp_string_container: ::std::default::Default::default(),
            strbuf: ::std::default::Default::default(),
            save_opt_bits: ::std::default::Default::default(),
        }
    }
}
pub type POSEContext = tagPOSEContext;
pub const MAIN_LOOP_ACTION_DO_NEXT_STEP: MAIN_LOOP_ACTION = 0;
pub const MAIN_LOOP_ACTION_DO_BREAK_MAIN_LOOP: MAIN_LOOP_ACTION = 1;
pub const MAIN_LOOP_ACTION_DO_EXIT_FUNCTION: MAIN_LOOP_ACTION = 2;
pub const MAIN_LOOP_ACTION_DO_CONTINUE_MAIN_LOOP: MAIN_LOOP_ACTION = 3;
pub type MAIN_LOOP_ACTION = u32;
pub type qInt = AT_RANK;
#[derive(Clone, Debug, PartialEq)]
pub struct tagQieue {
    pub Val: SourceMutPointer<qInt>,
    pub nTotLength: i32,
    pub nFirst: i32,
    pub nLength: i32,
}
impl ::std::default::Default for tagQieue {
    fn default() -> Self {
        Self {
            Val: ::std::default::Default::default(),
            nTotLength: ::std::default::Default::default(),
            nFirst: ::std::default::Default::default(),
            nLength: ::std::default::Default::default(),
        }
    }
}
pub type QUEUE = tagQieue;
#[derive(Clone, Debug, PartialEq)]
pub struct tagXYZCoord {
    pub xyz: [f64; 3usize],
}
impl ::std::default::Default for tagXYZCoord {
    fn default() -> Self {
        Self {
            xyz: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type XYZ_COORD = tagXYZCoord;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStructRestoreMode {
    pub bMetalAddFlower: i32,
    pub nMetalMinBondOrder: i32,
    pub nMetalInitEdgeFlow: i32,
    pub nMetalInitBondOrder: i32,
    pub nMetal2EndpointMinBondOrder: i32,
    pub nMetal2EndpointInitBondOrder: i32,
    pub nMetal2EndpointInitEdgeFlow: i32,
    pub nMetalFlowerParam_D: i32,
    pub nMetalMaxCharge_D: i32,
    pub bStereoRemovesMetalFlag: i32,
    pub bFixStereoBonds: i32,
}
impl ::std::default::Default for tagStructRestoreMode {
    fn default() -> Self {
        Self {
            bMetalAddFlower: ::std::default::Default::default(),
            nMetalMinBondOrder: ::std::default::Default::default(),
            nMetalInitEdgeFlow: ::std::default::Default::default(),
            nMetalInitBondOrder: ::std::default::Default::default(),
            nMetal2EndpointMinBondOrder: ::std::default::Default::default(),
            nMetal2EndpointInitBondOrder: ::std::default::Default::default(),
            nMetal2EndpointInitEdgeFlow: ::std::default::Default::default(),
            nMetalFlowerParam_D: ::std::default::Default::default(),
            nMetalMaxCharge_D: ::std::default::Default::default(),
            bStereoRemovesMetalFlag: ::std::default::Default::default(),
            bFixStereoBonds: ::std::default::Default::default(),
        }
    }
}
pub type SRM = tagStructRestoreMode;
#[derive(Clone, Debug, PartialEq)]
pub struct tagReversedInChI {
    pub pINChI: [SourceMutPointer<PINChI2>; 2usize],
    pub pINChI_Aux: [SourceMutPointer<PINChI_Aux2>; 2usize],
    pub num_components: [i32; 2usize],
    pub nRetVal: i32,
}
impl ::std::default::Default for tagReversedInChI {
    fn default() -> Self {
        Self {
            pINChI: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pINChI_Aux: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_components: ::std::array::from_fn(|_| ::std::default::Default::default()),
            nRetVal: ::std::default::Default::default(),
        }
    }
}
pub type REV_INCHI = tagReversedInChI;
#[derive(Clone, Debug, PartialEq)]
pub struct tagBfsQueue {
    pub q: SourceMutPointer<QUEUE>,
    pub nAtomLevel: SourceMutPointer<AT_RANK>,
    pub cSource: SourceMutPointer<S_CHAR>,
    pub num_at: i32,
    pub min_ring_size: AT_RANK,
}
impl ::std::default::Default for tagBfsQueue {
    fn default() -> Self {
        Self {
            q: ::std::default::Default::default(),
            nAtomLevel: ::std::default::Default::default(),
            cSource: ::std::default::Default::default(),
            num_at: ::std::default::Default::default(),
            min_ring_size: ::std::default::Default::default(),
        }
    }
}
pub type BFS_Q = tagBfsQueue;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInpAtomAddParities {
    pub bUsed0DParity: S_CHAR,
    pub p_parity: S_CHAR,
    pub p_orig_at_num: [AT_NUMB; 4usize],
    pub sb_ord: [S_CHAR; 3usize],
    pub sn_ord: [S_CHAR; 3usize],
    pub sb_parity: [S_CHAR; 3usize],
    pub sn_orig_at_num: [AT_NUMB; 3usize],
}
impl ::std::default::Default for tagInpAtomAddParities {
    fn default() -> Self {
        Self {
            bUsed0DParity: ::std::default::Default::default(),
            p_parity: ::std::default::Default::default(),
            p_orig_at_num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sb_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sn_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sb_parity: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sn_orig_at_num: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type inp_ATOM_STEREO = tagInpAtomAddParities;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomIonPrperies {
    pub cDoNotAddH: i8,
    pub cMetal: i8,
    pub cNumBondsToMetal: i8,
    pub cInitFlowToMetal: i8,
    pub cInitValenceToMetal: i8,
    pub cInitOrigValenceToMetal: i8,
    pub cMaxFlowToMetal: i8,
    pub cInitFreeValences: i8,
    pub cInitCharge: S_CHAR,
    pub cNumValenceElectrons: i8,
    pub cPeriodicRowNumber: i8,
    pub cMinRingSize: i8,
    pub cPeriodicNumber: U_CHAR,
    pub cnListIndex: S_CHAR,
    pub nCMinusGroupEdge: i32,
    pub nCPlusGroupEdge: i32,
    pub nMetalGroupEdge: i32,
    pub nTautGroupEdge: i32,
}
impl ::std::default::Default for tagAtomIonPrperies {
    fn default() -> Self {
        Self {
            cDoNotAddH: ::std::default::Default::default(),
            cMetal: ::std::default::Default::default(),
            cNumBondsToMetal: ::std::default::Default::default(),
            cInitFlowToMetal: ::std::default::Default::default(),
            cInitValenceToMetal: ::std::default::Default::default(),
            cInitOrigValenceToMetal: ::std::default::Default::default(),
            cMaxFlowToMetal: ::std::default::Default::default(),
            cInitFreeValences: ::std::default::Default::default(),
            cInitCharge: ::std::default::Default::default(),
            cNumValenceElectrons: ::std::default::Default::default(),
            cPeriodicRowNumber: ::std::default::Default::default(),
            cMinRingSize: ::std::default::Default::default(),
            cPeriodicNumber: ::std::default::Default::default(),
            cnListIndex: ::std::default::Default::default(),
            nCMinusGroupEdge: ::std::default::Default::default(),
            nCPlusGroupEdge: ::std::default::Default::default(),
            nMetalGroupEdge: ::std::default::Default::default(),
            nTautGroupEdge: ::std::default::Default::default(),
        }
    }
}
pub type VAL_AT = tagAtomIonPrperies;
pub const tagTgRestoreFlags_TGRF_MINUS_FIRST: tagTgRestoreFlags = 1;
pub type tagTgRestoreFlags = u32;
pub use self::tagTgRestoreFlags as TGRF;
#[derive(Clone, Debug, PartialEq)]
pub struct tagTCGroup {
    pub type_: i32,
    pub ord_num: i32,
    pub num_edges: i32,
    pub st_cap: i32,
    pub st_flow: i32,
    pub edges_cap: i32,
    pub edges_flow: i32,
    pub nVertexNumber: i32,
    pub nForwardEdge: i32,
    pub nBackwardEdge: i32,
    pub tg_num_H: i16,
    pub tg_num_Minus: i16,
    pub tg_set_Minus: Vertex,
    pub tg_RestoreFlags: i16,
}
impl ::std::default::Default for tagTCGroup {
    fn default() -> Self {
        Self {
            type_: ::std::default::Default::default(),
            ord_num: ::std::default::Default::default(),
            num_edges: ::std::default::Default::default(),
            st_cap: ::std::default::Default::default(),
            st_flow: ::std::default::Default::default(),
            edges_cap: ::std::default::Default::default(),
            edges_flow: ::std::default::Default::default(),
            nVertexNumber: ::std::default::Default::default(),
            nForwardEdge: ::std::default::Default::default(),
            nBackwardEdge: ::std::default::Default::default(),
            tg_num_H: ::std::default::Default::default(),
            tg_num_Minus: ::std::default::Default::default(),
            tg_set_Minus: ::std::default::Default::default(),
            tg_RestoreFlags: ::std::default::Default::default(),
        }
    }
}
pub type TC_GROUP = tagTCGroup;
pub const tagTCGroupTypes_TCG_None: tagTCGroupTypes = -1;
pub const tagTCGroupTypes_TCG_Plus0: tagTCGroupTypes = 0;
pub const tagTCGroupTypes_TCG_Plus1: tagTCGroupTypes = 1;
pub const tagTCGroupTypes_TCG_Minus0: tagTCGroupTypes = 2;
pub const tagTCGroupTypes_TCG_Minus1: tagTCGroupTypes = 3;
pub const tagTCGroupTypes_TCG_Plus_C0: tagTCGroupTypes = 4;
pub const tagTCGroupTypes_TCG_Plus_C1: tagTCGroupTypes = 5;
pub const tagTCGroupTypes_TCG_Minus_C0: tagTCGroupTypes = 6;
pub const tagTCGroupTypes_TCG_Minus_C1: tagTCGroupTypes = 7;
pub const tagTCGroupTypes_TCG_Plus_M0: tagTCGroupTypes = 8;
pub const tagTCGroupTypes_TCG_Plus_M1: tagTCGroupTypes = 9;
pub const tagTCGroupTypes_TCG_Minus_M0: tagTCGroupTypes = 10;
pub const tagTCGroupTypes_TCG_Minus_M1: tagTCGroupTypes = 11;
pub const tagTCGroupTypes_TCG_MeFlower0: tagTCGroupTypes = 12;
pub const tagTCGroupTypes_TCG_MeFlower1: tagTCGroupTypes = 13;
pub const tagTCGroupTypes_TCG_MeFlower2: tagTCGroupTypes = 14;
pub const tagTCGroupTypes_TCG_MeFlower3: tagTCGroupTypes = 15;
pub const tagTCGroupTypes_TCG_Plus: tagTCGroupTypes = 16;
pub const tagTCGroupTypes_TCG_Minus: tagTCGroupTypes = 17;
pub const tagTCGroupTypes_NUM_TCGROUP_TYPES: tagTCGroupTypes = 18;
pub type tagTCGroupTypes = i32;
pub use self::tagTCGroupTypes as TCGR_TYPE;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAllTCGroups {
    pub pTCG: SourceMutPointer<TC_GROUP>,
    pub num_tc_groups: i32,
    pub max_tc_groups: i32,
    pub nGroup: [i32; 18usize],
    pub nVertices: i32,
    pub nEdges: i32,
    pub nAddIedges: i32,
    pub num_atoms: i32,
    pub num_bonds: i32,
    pub num_tgroups: i32,
    pub num_tgroup_edges: i32,
    pub tgroup_charge: i32,
    pub charge_on_atoms: i32,
    pub added_charge: i32,
    pub total_charge: i32,
    pub total_electrons: i32,
    pub total_electrons_metals: i32,
    pub num_metal_atoms: i32,
    pub num_metal_bonds: i32,
    pub nEdge4charge: i32,
    pub nEdgePlus: i32,
    pub nEdgeMinus: i32,
    pub iComponent: i32,
    pub iAtNoOffset: i32,
}
impl ::std::default::Default for tagAllTCGroups {
    fn default() -> Self {
        Self {
            pTCG: ::std::default::Default::default(),
            num_tc_groups: ::std::default::Default::default(),
            max_tc_groups: ::std::default::Default::default(),
            nGroup: ::std::array::from_fn(|_| ::std::default::Default::default()),
            nVertices: ::std::default::Default::default(),
            nEdges: ::std::default::Default::default(),
            nAddIedges: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_bonds: ::std::default::Default::default(),
            num_tgroups: ::std::default::Default::default(),
            num_tgroup_edges: ::std::default::Default::default(),
            tgroup_charge: ::std::default::Default::default(),
            charge_on_atoms: ::std::default::Default::default(),
            added_charge: ::std::default::Default::default(),
            total_charge: ::std::default::Default::default(),
            total_electrons: ::std::default::Default::default(),
            total_electrons_metals: ::std::default::Default::default(),
            num_metal_atoms: ::std::default::Default::default(),
            num_metal_bonds: ::std::default::Default::default(),
            nEdge4charge: ::std::default::Default::default(),
            nEdgePlus: ::std::default::Default::default(),
            nEdgeMinus: ::std::default::Default::default(),
            iComponent: ::std::default::Default::default(),
            iAtNoOffset: ::std::default::Default::default(),
        }
    }
}
pub type ALL_TC_GROUPS = tagAllTCGroups;
#[derive(Clone, Debug, PartialEq)]
pub struct tagEdgeList {
    pub num_alloc: i32,
    pub num_edges: i32,
    pub pnEdges: SourceMutPointer<EdgeIndex>,
}
impl ::std::default::Default for tagEdgeList {
    fn default() -> Self {
        Self {
            num_alloc: ::std::default::Default::default(),
            num_edges: ::std::default::Default::default(),
            pnEdges: ::std::default::Default::default(),
        }
    }
}
pub type EDGE_LIST = tagEdgeList;
#[derive(Clone, Debug, PartialEq)]
pub struct tagChargeValence {
    pub nValence: i32,
    pub nCharge: i32,
    pub nValenceOrderingNumber: i32,
}
impl ::std::default::Default for tagChargeValence {
    fn default() -> Self {
        Self {
            nValence: ::std::default::Default::default(),
            nCharge: ::std::default::Default::default(),
            nValenceOrderingNumber: ::std::default::Default::default(),
        }
    }
}
pub type CHARGE_VAL = tagChargeValence;
#[derive(Clone, Debug, PartialEq)]
pub struct tagChargeChangeCandidate {
    pub iat: Vertex,
    pub num_bonds: i8,
    pub chem_valence: i8,
    pub cMetal: i8,
    pub cNumBondsToMetal: i8,
    pub cNumValenceElectrons: i8,
    pub cPeriodicRowNumber: i8,
    pub cNumChargeStates: i8,
    pub el_number: U_CHAR,
}
impl ::std::default::Default for tagChargeChangeCandidate {
    fn default() -> Self {
        Self {
            iat: ::std::default::Default::default(),
            num_bonds: ::std::default::Default::default(),
            chem_valence: ::std::default::Default::default(),
            cMetal: ::std::default::Default::default(),
            cNumBondsToMetal: ::std::default::Default::default(),
            cNumValenceElectrons: ::std::default::Default::default(),
            cPeriodicRowNumber: ::std::default::Default::default(),
            cNumChargeStates: ::std::default::Default::default(),
            el_number: ::std::default::Default::default(),
        }
    }
}
pub type CC_CAND = tagChargeChangeCandidate;
#[derive(Clone, Debug, PartialEq)]
pub struct tagOneComponentRemovedAndExchangeableH {
    pub nNumRemovedProtons: NUM_H,
    pub nNumRemovedIsotopicH: [NUM_H; 3usize],
}
impl ::std::default::Default for tagOneComponentRemovedAndExchangeableH {
    fn default() -> Self {
        Self {
            nNumRemovedProtons: ::std::default::Default::default(),
            nNumRemovedIsotopicH: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type COMPONENT_REM_PROTONS = tagOneComponentRemovedAndExchangeableH;
#[derive(Clone, Debug, PartialEq)]
pub struct tagRemovedAndExchangeableH {
    pub nNumRemovedProtons: NUM_H,
    pub nNumRemovedIsotopicH: [NUM_H; 3usize],
    pub pNumProtons: SourceMutPointer<COMPONENT_REM_PROTONS>,
}
impl ::std::default::Default for tagRemovedAndExchangeableH {
    fn default() -> Self {
        Self {
            nNumRemovedProtons: ::std::default::Default::default(),
            nNumRemovedIsotopicH: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pNumProtons: ::std::default::Default::default(),
        }
    }
}
pub type REM_PROTONS = tagRemovedAndExchangeableH;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInputInChI {
    pub pInpInChI: [[SourceMutPointer<INChI>; 2usize]; 2usize],
    pub nNumComponents: [[i32; 2usize]; 2usize],
    pub nNumProtons: [[REM_PROTONS; 2usize]; 2usize],
    pub s: [[[i32; 2usize]; 2usize]; 2usize],
    pub num_inp: i64,
    pub atom: SourceMutPointer<inp_ATOM>,
    pub num_atoms: i32,
    pub num_explicit_H: i32,
    pub CompareInchiFlags: [[INCHI_MODE; 2usize]; 2usize],
    pub polymer: SourceMutPointer<OAD_Polymer>,
    pub v3000: SourceMutPointer<OAD_V3000>,
    pub valid_polymer: i32,
}
impl ::std::default::Default for tagInputInChI {
    fn default() -> Self {
        Self {
            pInpInChI: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            nNumComponents: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            nNumProtons: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            s: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| {
                    ::std::array::from_fn(|_| ::std::default::Default::default())
                })
            }),
            num_inp: ::std::default::Default::default(),
            atom: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_explicit_H: ::std::default::Default::default(),
            CompareInchiFlags: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            polymer: ::std::default::Default::default(),
            v3000: ::std::default::Default::default(),
            valid_polymer: ::std::default::Default::default(),
        }
    }
}
pub type InpInChI = tagInputInChI;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStructFromInChI {
    pub at: SourceMutPointer<inp_ATOM>,
    pub st: SourceMutPointer<inp_ATOM_STEREO>,
    pub at2: SourceMutPointer<inp_ATOM>,
    pub ti: T_GROUP_INFO,
    pub endpoint: SourceMutPointer<AT_NUMB>,
    pub fixed_H: SourceMutPointer<S_CHAR>,
    pub pXYZ: SourceMutPointer<XYZ_COORD>,
    pub num_atoms: i32,
    pub num_deleted_H: i32,
    pub nNumRemovedProtonsMobHInChI: i32,
    pub charge: S_CHAR,
    pub bIsotopic: i8,
    pub pBNS: SourceMutPointer<BN_STRUCT>,
    pub pBD: SourceMutPointer<BN_DATA>,
    pub pSrm: SourceConstPointer<SRM>,
    pub bMobileH: i8,
    pub iINCHI: i8,
    pub bFixedHExists: i8,
    pub RevInChI: REV_INCHI,
    pub nRemovedProtonsByNormFromRevrs: i32,
    pub nNumRemovedProtonsByRevrs: i32,
    pub bExtract: i32,
    pub pOneINChI: [SourceMutPointer<INChI>; 2usize],
    pub pOneINChI_Aux: [SourceMutPointer<INChI_Aux>; 2usize],
    pub pOne_norm_data: [SourceMutPointer<INP_ATOM_DATA>; 2usize],
    pub pOne_fixed_H: SourceMutPointer<S_CHAR>,
    pub One_ti: T_GROUP_INFO,
    pub nOneINChI_bMobileH: i32,
    pub nNumRemovedProtons: i32,
    pub nAtno2Canon: [SourceMutPointer<AT_NUMB>; 2usize],
    pub nCanon2Atno: [SourceMutPointer<AT_NUMB>; 2usize],
    pub nError: i32,
    pub iInchiRec: i8,
    pub iMobileH: i8,
    pub bDeleted: i8,
    pub num_inp_actual: i64,
    pub pbfsq: SourceMutPointer<BFS_Q>,
    pub pVA: SourceMutPointer<VAL_AT>,
    pub nLink: i32,
    pub bPostProcessed: i32,
    pub nChargeRevrs: i32,
    pub nChargeInChI: i32,
    pub n_zy: i32,
    pub n_pzz: i32,
}
impl ::std::default::Default for tagStructFromInChI {
    fn default() -> Self {
        Self {
            at: ::std::default::Default::default(),
            st: ::std::default::Default::default(),
            at2: ::std::default::Default::default(),
            ti: ::std::default::Default::default(),
            endpoint: ::std::default::Default::default(),
            fixed_H: ::std::default::Default::default(),
            pXYZ: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_deleted_H: ::std::default::Default::default(),
            nNumRemovedProtonsMobHInChI: ::std::default::Default::default(),
            charge: ::std::default::Default::default(),
            bIsotopic: ::std::default::Default::default(),
            pBNS: ::std::default::Default::default(),
            pBD: ::std::default::Default::default(),
            pSrm: ::std::default::Default::default(),
            bMobileH: ::std::default::Default::default(),
            iINCHI: ::std::default::Default::default(),
            bFixedHExists: ::std::default::Default::default(),
            RevInChI: ::std::default::Default::default(),
            nRemovedProtonsByNormFromRevrs: ::std::default::Default::default(),
            nNumRemovedProtonsByRevrs: ::std::default::Default::default(),
            bExtract: ::std::default::Default::default(),
            pOneINChI: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pOneINChI_Aux: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pOne_norm_data: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pOne_fixed_H: ::std::default::Default::default(),
            One_ti: ::std::default::Default::default(),
            nOneINChI_bMobileH: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            nAtno2Canon: ::std::array::from_fn(|_| ::std::default::Default::default()),
            nCanon2Atno: ::std::array::from_fn(|_| ::std::default::Default::default()),
            nError: ::std::default::Default::default(),
            iInchiRec: ::std::default::Default::default(),
            iMobileH: ::std::default::Default::default(),
            bDeleted: ::std::default::Default::default(),
            num_inp_actual: ::std::default::Default::default(),
            pbfsq: ::std::default::Default::default(),
            pVA: ::std::default::Default::default(),
            nLink: ::std::default::Default::default(),
            bPostProcessed: ::std::default::Default::default(),
            nChargeRevrs: ::std::default::Default::default(),
            nChargeInChI: ::std::default::Default::default(),
            n_zy: ::std::default::Default::default(),
            n_pzz: ::std::default::Default::default(),
        }
    }
}
pub type StrFromINChI = tagStructFromInChI;
#[derive(Clone, Debug, PartialEq)]
pub struct tagVertCapFlow {
    pub type_: S_SHORT,
    pub cap: S_CHAR,
    pub flow: S_CHAR,
    pub valence: S_CHAR,
}
impl ::std::default::Default for tagVertCapFlow {
    fn default() -> Self {
        Self {
            type_: ::std::default::Default::default(),
            cap: ::std::default::Default::default(),
            flow: ::std::default::Default::default(),
            valence: ::std::default::Default::default(),
        }
    }
}
pub type VCF = tagVertCapFlow;
#[derive(Clone, Debug, PartialEq)]
pub struct tagEdgeCapFlow {
    pub neigh: S_SHORT,
    pub cap: S_CHAR,
    pub bForbiddenEdge: S_CHAR,
    pub flow: S_CHAR,
}
impl ::std::default::Default for tagEdgeCapFlow {
    fn default() -> Self {
        Self {
            neigh: ::std::default::Default::default(),
            cap: ::std::default::Default::default(),
            bForbiddenEdge: ::std::default::Default::default(),
            flow: ::std::default::Default::default(),
        }
    }
}
pub type ECF = tagEdgeCapFlow;
#[derive(Clone, Debug, PartialEq)]
pub struct tagChargeNodes {
    pub v: VCF,
    pub e: [ECF; 3usize],
}
impl ::std::default::Default for tagChargeNodes {
    fn default() -> Self {
        Self {
            v: ::std::default::Default::default(),
            e: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type C_NODE = tagChargeNodes;
#[derive(Clone, Debug, PartialEq)]
pub struct tagChargeNodeList {
    pub pCN: SourceConstPointer<C_NODE>,
    pub bits: i32,
    pub nInitialCharge: i32,
    pub len: i32,
}
impl ::std::default::Default for tagChargeNodeList {
    fn default() -> Self {
        Self {
            pCN: ::std::default::Default::default(),
            bits: ::std::default::Default::default(),
            nInitialCharge: ::std::default::Default::default(),
            len: ::std::default::Default::default(),
        }
    }
}
pub type CN_LIST = tagChargeNodeList;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomsCmpTwoFixedH {
    pub endptInChI: AT_NUMB,
    pub endptRevrs: AT_NUMB,
    pub atomNumber: AT_NUMB,
    pub nValElectr: U_CHAR,
    pub nPeriodNum: U_CHAR,
    pub nFixHInChI: S_CHAR,
    pub nFixHRevrs: S_CHAR,
    pub nMobHInChI: S_CHAR,
    pub nMobHRevrs: S_CHAR,
    pub nNumHRevrs: S_CHAR,
    pub nAtChargeRevrs: S_CHAR,
    pub nValue: S_CHAR,
}
impl ::std::default::Default for tagAtomsCmpTwoFixedH {
    fn default() -> Self {
        Self {
            endptInChI: ::std::default::Default::default(),
            endptRevrs: ::std::default::Default::default(),
            atomNumber: ::std::default::Default::default(),
            nValElectr: ::std::default::Default::default(),
            nPeriodNum: ::std::default::Default::default(),
            nFixHInChI: ::std::default::Default::default(),
            nFixHRevrs: ::std::default::Default::default(),
            nMobHInChI: ::std::default::Default::default(),
            nMobHRevrs: ::std::default::Default::default(),
            nNumHRevrs: ::std::default::Default::default(),
            nAtChargeRevrs: ::std::default::Default::default(),
            nValue: ::std::default::Default::default(),
        }
    }
}
pub type CMP2FHATOMS = tagAtomsCmpTwoFixedH;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStructCmpTwoFixedH {
    pub c2at: [CMP2FHATOMS; 256usize],
    pub len_c2at: i16,
    pub nNumRemHInChI: i16,
    pub nNumRemHRevrs: i16,
    pub nNumTgInChI: i16,
    pub nNumTgRevrs: i16,
    pub nNumEndpInChI: i16,
    pub nNumEndpRevrs: i16,
    pub nNumTgDiffMinus: i16,
    pub nNumTgDiffH: i16,
    pub nNumTgMInChI: i16,
    pub nNumTgHInChI: i16,
    pub nNumTgMRevrs: i16,
    pub nNumTgHRevrs: i16,
    pub nChargeFixHInChI: S_CHAR,
    pub nChargeMobHInChI: S_CHAR,
    pub nChargeFixHRevrs: S_CHAR,
    pub nChargeMobHRevrs: S_CHAR,
    pub nChargeFixHRevrsNonMetal: S_CHAR,
    pub nChargeMobHRevrsNonMetal: S_CHAR,
    pub bFixedHLayerExistsRevrs: i8,
    pub bHasDifference: i8,
    pub nNumDiffMobH: U_CHAR,
}
impl ::std::default::Default for tagStructCmpTwoFixedH {
    fn default() -> Self {
        Self {
            c2at: ::std::array::from_fn(|_| ::std::default::Default::default()),
            len_c2at: ::std::default::Default::default(),
            nNumRemHInChI: ::std::default::Default::default(),
            nNumRemHRevrs: ::std::default::Default::default(),
            nNumTgInChI: ::std::default::Default::default(),
            nNumTgRevrs: ::std::default::Default::default(),
            nNumEndpInChI: ::std::default::Default::default(),
            nNumEndpRevrs: ::std::default::Default::default(),
            nNumTgDiffMinus: ::std::default::Default::default(),
            nNumTgDiffH: ::std::default::Default::default(),
            nNumTgMInChI: ::std::default::Default::default(),
            nNumTgHInChI: ::std::default::Default::default(),
            nNumTgMRevrs: ::std::default::Default::default(),
            nNumTgHRevrs: ::std::default::Default::default(),
            nChargeFixHInChI: ::std::default::Default::default(),
            nChargeMobHInChI: ::std::default::Default::default(),
            nChargeFixHRevrs: ::std::default::Default::default(),
            nChargeMobHRevrs: ::std::default::Default::default(),
            nChargeFixHRevrsNonMetal: ::std::default::Default::default(),
            nChargeMobHRevrsNonMetal: ::std::default::Default::default(),
            bFixedHLayerExistsRevrs: ::std::default::Default::default(),
            bHasDifference: ::std::default::Default::default(),
            nNumDiffMobH: ::std::default::Default::default(),
        }
    }
}
pub type CMP2FHINCHI = tagStructCmpTwoFixedH;
#[derive(Clone, Debug, PartialEq)]
pub struct tagAtomsCmpTwoMobileH {
    pub endptInChI: AT_NUMB,
    pub endptRevrs: AT_NUMB,
    pub atomNumber: AT_NUMB,
    pub nValElectr: U_CHAR,
    pub nPeriodNum: U_CHAR,
    pub nMobHInChI: S_CHAR,
    pub nMobHRevrs: S_CHAR,
    pub nNumHRevrs: S_CHAR,
    pub nAtChargeRevrs: S_CHAR,
    pub nValue: S_CHAR,
}
impl ::std::default::Default for tagAtomsCmpTwoMobileH {
    fn default() -> Self {
        Self {
            endptInChI: ::std::default::Default::default(),
            endptRevrs: ::std::default::Default::default(),
            atomNumber: ::std::default::Default::default(),
            nValElectr: ::std::default::Default::default(),
            nPeriodNum: ::std::default::Default::default(),
            nMobHInChI: ::std::default::Default::default(),
            nMobHRevrs: ::std::default::Default::default(),
            nNumHRevrs: ::std::default::Default::default(),
            nAtChargeRevrs: ::std::default::Default::default(),
            nValue: ::std::default::Default::default(),
        }
    }
}
pub type CMP2MHATOMS = tagAtomsCmpTwoMobileH;
#[derive(Clone, Debug, PartialEq)]
pub struct tagStructCmpTwoMobileH {
    pub c2at: [CMP2MHATOMS; 256usize],
    pub len_c2at: i16,
    pub nNumRemHInChI: i16,
    pub nNumRemHRevrs: i16,
    pub nNumTgInChI: i16,
    pub nNumTgRevrs: i16,
    pub nNumEndpInChI: i16,
    pub nNumEndpRevrs: i16,
    pub nNumTgDiffMinus: i16,
    pub nNumTgDiffH: i16,
    pub nNumTgMInChI: i16,
    pub nNumTgHInChI: i16,
    pub nNumTgOInChI: i16,
    pub nNumTgNInChI: i16,
    pub nNumTgMRevrs: i16,
    pub nNumTgHRevrs: i16,
    pub nNumTgORevrs: i16,
    pub nNumTgNRevrs: i16,
    pub nNumTgOMinusRevrs: i16,
    pub nNumTgOHRevrs: i16,
    pub nNumTgDBORevrs: i16,
    pub nNumTgNMinusRevrs: i16,
    pub nNumTgNHMinusRevrs: i16,
    pub nNumTgNHRevrs: i16,
    pub nNumTgNH2Revrs: i16,
    pub nNumTgDBNHRevrs: i16,
    pub nNumTgDBNMinusRevrs: i16,
    pub nNumTgDBNRevrs: i16,
    pub nNumTgOMinusInChI: i16,
    pub nNumTgOHInChI: i16,
    pub nNumTgDBOInChI: i16,
    pub nNumTgNMinusInChI: i16,
    pub nNumTgNHMinusInChI: i16,
    pub nNumTgNHInChI: i16,
    pub nNumTgNH2InChI: i16,
    pub nNumTgDBNHInChI: i16,
    pub nNumTgDBNMinusInChI: i16,
    pub nNumTgDBNInChI: i16,
    pub nChargeMobHInChI: S_CHAR,
    pub nChargeMobHRevrs: S_CHAR,
    pub nChargeMobHRevrsNonMetal: S_CHAR,
    pub bFixedHLayerExistsRevrs: i8,
    pub bHasDifference: i8,
    pub nNumDiffMobH: U_CHAR,
}
impl ::std::default::Default for tagStructCmpTwoMobileH {
    fn default() -> Self {
        Self {
            c2at: ::std::array::from_fn(|_| ::std::default::Default::default()),
            len_c2at: ::std::default::Default::default(),
            nNumRemHInChI: ::std::default::Default::default(),
            nNumRemHRevrs: ::std::default::Default::default(),
            nNumTgInChI: ::std::default::Default::default(),
            nNumTgRevrs: ::std::default::Default::default(),
            nNumEndpInChI: ::std::default::Default::default(),
            nNumEndpRevrs: ::std::default::Default::default(),
            nNumTgDiffMinus: ::std::default::Default::default(),
            nNumTgDiffH: ::std::default::Default::default(),
            nNumTgMInChI: ::std::default::Default::default(),
            nNumTgHInChI: ::std::default::Default::default(),
            nNumTgOInChI: ::std::default::Default::default(),
            nNumTgNInChI: ::std::default::Default::default(),
            nNumTgMRevrs: ::std::default::Default::default(),
            nNumTgHRevrs: ::std::default::Default::default(),
            nNumTgORevrs: ::std::default::Default::default(),
            nNumTgNRevrs: ::std::default::Default::default(),
            nNumTgOMinusRevrs: ::std::default::Default::default(),
            nNumTgOHRevrs: ::std::default::Default::default(),
            nNumTgDBORevrs: ::std::default::Default::default(),
            nNumTgNMinusRevrs: ::std::default::Default::default(),
            nNumTgNHMinusRevrs: ::std::default::Default::default(),
            nNumTgNHRevrs: ::std::default::Default::default(),
            nNumTgNH2Revrs: ::std::default::Default::default(),
            nNumTgDBNHRevrs: ::std::default::Default::default(),
            nNumTgDBNMinusRevrs: ::std::default::Default::default(),
            nNumTgDBNRevrs: ::std::default::Default::default(),
            nNumTgOMinusInChI: ::std::default::Default::default(),
            nNumTgOHInChI: ::std::default::Default::default(),
            nNumTgDBOInChI: ::std::default::Default::default(),
            nNumTgNMinusInChI: ::std::default::Default::default(),
            nNumTgNHMinusInChI: ::std::default::Default::default(),
            nNumTgNHInChI: ::std::default::Default::default(),
            nNumTgNH2InChI: ::std::default::Default::default(),
            nNumTgDBNHInChI: ::std::default::Default::default(),
            nNumTgDBNMinusInChI: ::std::default::Default::default(),
            nNumTgDBNInChI: ::std::default::Default::default(),
            nChargeMobHInChI: ::std::default::Default::default(),
            nChargeMobHRevrs: ::std::default::Default::default(),
            nChargeMobHRevrsNonMetal: ::std::default::Default::default(),
            bFixedHLayerExistsRevrs: ::std::default::Default::default(),
            bHasDifference: ::std::default::Default::default(),
            nNumDiffMobH: ::std::default::Default::default(),
        }
    }
}
pub type CMP2MHINCHI = tagStructCmpTwoMobileH;
pub type UINT32 = u32;
pub type UINT16 = u16;
pub const tagINCHIRadical_INCHI_RADICAL_NONE: tagINCHIRadical = 0;
pub const tagINCHIRadical_INCHI_RADICAL_SINGLET: tagINCHIRadical = 1;
pub const tagINCHIRadical_INCHI_RADICAL_DOUBLET: tagINCHIRadical = 2;
pub const tagINCHIRadical_INCHI_RADICAL_TRIPLET: tagINCHIRadical = 3;
pub type tagINCHIRadical = u32;
pub use self::tagINCHIRadical as inchi_Radical;
pub const tagINCHIBondType_INCHI_BOND_TYPE_NONE: tagINCHIBondType = 0;
pub const tagINCHIBondType_INCHI_BOND_TYPE_SINGLE: tagINCHIBondType = 1;
pub const tagINCHIBondType_INCHI_BOND_TYPE_DOUBLE: tagINCHIBondType = 2;
pub const tagINCHIBondType_INCHI_BOND_TYPE_TRIPLE: tagINCHIBondType = 3;
pub const tagINCHIBondType_INCHI_BOND_TYPE_ALTERN: tagINCHIBondType = 4;
pub type tagINCHIBondType = u32;
pub use self::tagINCHIBondType as inchi_BondType;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_NONE: tagINCHIBondStereo2D = 0;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1UP: tagINCHIBondStereo2D = 1;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1EITHER: tagINCHIBondStereo2D = 4;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_1DOWN: tagINCHIBondStereo2D = 6;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2UP: tagINCHIBondStereo2D = -1;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2EITHER: tagINCHIBondStereo2D = -4;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_SINGLE_2DOWN: tagINCHIBondStereo2D = -6;
pub const tagINCHIBondStereo2D_INCHI_BOND_STEREO_DOUBLE_EITHER: tagINCHIBondStereo2D = 3;
pub type tagINCHIBondStereo2D = i32;
pub use self::tagINCHIBondStereo2D as inchi_BondStereo2D;
pub type AT_NUM = S_SHORT;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInchiAtom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub neighbor: [AT_NUM; 20usize],
    pub bond_type: [S_CHAR; 20usize],
    pub bond_stereo: [S_CHAR; 20usize],
    pub elname: [i8; 6usize],
    pub num_bonds: AT_NUM,
    pub num_iso_H: [S_CHAR; 4usize],
    pub isotopic_mass: AT_NUM,
    pub radical: S_CHAR,
    pub charge: S_CHAR,
}
impl ::std::default::Default for tagInchiAtom {
    fn default() -> Self {
        Self {
            x: ::std::default::Default::default(),
            y: ::std::default::Default::default(),
            z: ::std::default::Default::default(),
            neighbor: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bond_type: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bond_stereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            elname: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_bonds: ::std::default::Default::default(),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            isotopic_mass: ::std::default::Default::default(),
            radical: ::std::default::Default::default(),
            charge: ::std::default::Default::default(),
        }
    }
}
pub type inchi_Atom = tagInchiAtom;
pub const tagINCHIStereoType0D_INCHI_StereoType_None: tagINCHIStereoType0D = 0;
pub const tagINCHIStereoType0D_INCHI_StereoType_DoubleBond: tagINCHIStereoType0D = 1;
pub const tagINCHIStereoType0D_INCHI_StereoType_Tetrahedral: tagINCHIStereoType0D = 2;
pub const tagINCHIStereoType0D_INCHI_StereoType_Allene: tagINCHIStereoType0D = 3;
pub type tagINCHIStereoType0D = u32;
pub use self::tagINCHIStereoType0D as inchi_StereoType0D;
pub const tagINCHIStereoParity0D_INCHI_PARITY_NONE: tagINCHIStereoParity0D = 0;
pub const tagINCHIStereoParity0D_INCHI_PARITY_ODD: tagINCHIStereoParity0D = 1;
pub const tagINCHIStereoParity0D_INCHI_PARITY_EVEN: tagINCHIStereoParity0D = 2;
pub const tagINCHIStereoParity0D_INCHI_PARITY_UNKNOWN: tagINCHIStereoParity0D = 3;
pub const tagINCHIStereoParity0D_INCHI_PARITY_UNDEFINED: tagINCHIStereoParity0D = 4;
pub type tagINCHIStereoParity0D = u32;
pub use self::tagINCHIStereoParity0D as inchi_StereoParity0D;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHIStereo0D {
    pub neighbor: [AT_NUM; 4usize],
    pub central_atom: AT_NUM,
    pub type_: S_CHAR,
    pub parity: S_CHAR,
}
impl ::std::default::Default for tagINCHIStereo0D {
    fn default() -> Self {
        Self {
            neighbor: ::std::array::from_fn(|_| ::std::default::Default::default()),
            central_atom: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            parity: ::std::default::Default::default(),
        }
    }
}
pub type inchi_Stereo0D = tagINCHIStereo0D;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHI_Input {
    pub atom: SourceMutPointer<inchi_Atom>,
    pub stereo0D: SourceMutPointer<inchi_Stereo0D>,
    pub szOptions: SourceMutPointer<i8>,
    pub num_atoms: AT_NUM,
    pub num_stereo0D: AT_NUM,
}
impl ::std::default::Default for tagINCHI_Input {
    fn default() -> Self {
        Self {
            atom: ::std::default::Default::default(),
            stereo0D: ::std::default::Default::default(),
            szOptions: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_stereo0D: ::std::default::Default::default(),
        }
    }
}
pub type inchi_Input = tagINCHI_Input;
#[derive(Clone, Debug, PartialEq)]
pub struct inchi_Input_PolymerUnit {
    pub id: i32,
    pub type_: i32,
    pub subtype: i32,
    pub conn: i32,
    pub label: i32,
    pub na: i32,
    pub nb: i32,
    pub xbr1: [f64; 4usize],
    pub xbr2: [f64; 4usize],
    pub smt: [i8; 80usize],
    pub alist: SourceMutPointer<i32>,
    pub blist: SourceMutPointer<i32>,
}
impl ::std::default::Default for inchi_Input_PolymerUnit {
    fn default() -> Self {
        Self {
            id: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            subtype: ::std::default::Default::default(),
            conn: ::std::default::Default::default(),
            label: ::std::default::Default::default(),
            na: ::std::default::Default::default(),
            nb: ::std::default::Default::default(),
            xbr1: ::std::array::from_fn(|_| ::std::default::Default::default()),
            xbr2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            smt: ::std::array::from_fn(|_| ::std::default::Default::default()),
            alist: ::std::default::Default::default(),
            blist: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct inchi_Input_Polymer {
    pub units: SourceMutPointer<SourceMutPointer<inchi_Input_PolymerUnit>>,
    pub n: i32,
}
impl ::std::default::Default for inchi_Input_Polymer {
    fn default() -> Self {
        Self {
            units: ::std::default::Default::default(),
            n: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct inchi_Input_V3000 {
    pub n_non_star_atoms: i32,
    pub n_star_atoms: i32,
    pub atom_index_orig: SourceMutPointer<i32>,
    pub atom_index_fin: SourceMutPointer<i32>,
    pub n_sgroups: i32,
    pub n_3d_constraints: i32,
    pub n_collections: i32,
    pub n_non_haptic_bonds: i32,
    pub n_haptic_bonds: i32,
    pub lists_haptic_bonds: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_steabs: i32,
    pub lists_steabs: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_sterel: i32,
    pub lists_sterel: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_sterac: i32,
    pub lists_sterac: SourceMutPointer<SourceMutPointer<i32>>,
}
impl ::std::default::Default for inchi_Input_V3000 {
    fn default() -> Self {
        Self {
            n_non_star_atoms: ::std::default::Default::default(),
            n_star_atoms: ::std::default::Default::default(),
            atom_index_orig: ::std::default::Default::default(),
            atom_index_fin: ::std::default::Default::default(),
            n_sgroups: ::std::default::Default::default(),
            n_3d_constraints: ::std::default::Default::default(),
            n_collections: ::std::default::Default::default(),
            n_non_haptic_bonds: ::std::default::Default::default(),
            n_haptic_bonds: ::std::default::Default::default(),
            lists_haptic_bonds: ::std::default::Default::default(),
            n_steabs: ::std::default::Default::default(),
            lists_steabs: ::std::default::Default::default(),
            n_sterel: ::std::default::Default::default(),
            lists_sterel: ::std::default::Default::default(),
            n_sterac: ::std::default::Default::default(),
            lists_sterac: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct inchi_InputEx {
    pub atom: SourceMutPointer<inchi_Atom>,
    pub stereo0D: SourceMutPointer<inchi_Stereo0D>,
    pub szOptions: SourceMutPointer<i8>,
    pub num_atoms: AT_NUM,
    pub num_stereo0D: AT_NUM,
    pub polymer: SourceMutPointer<inchi_Input_Polymer>,
    pub v3000: SourceMutPointer<inchi_Input_V3000>,
}
impl ::std::default::Default for inchi_InputEx {
    fn default() -> Self {
        Self {
            atom: ::std::default::Default::default(),
            stereo0D: ::std::default::Default::default(),
            szOptions: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_stereo0D: ::std::default::Default::default(),
            polymer: ::std::default::Default::default(),
            v3000: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHI_InputINCHI {
    pub szInChI: SourceMutPointer<i8>,
    pub szOptions: SourceMutPointer<i8>,
}
impl ::std::default::Default for tagINCHI_InputINCHI {
    fn default() -> Self {
        Self {
            szInChI: ::std::default::Default::default(),
            szOptions: ::std::default::Default::default(),
        }
    }
}
pub type inchi_InputINCHI = tagINCHI_InputINCHI;
pub type inchi_Output_PolymerUnit = inchi_Input_PolymerUnit;
pub type inchi_Output_Polymer = inchi_Input_Polymer;
pub type inchi_Output_V3000 = inchi_Input_V3000;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHI_Output {
    pub szInChI: SourceMutPointer<i8>,
    pub szAuxInfo: SourceMutPointer<i8>,
    pub szMessage: SourceMutPointer<i8>,
    pub szLog: SourceMutPointer<i8>,
}
impl ::std::default::Default for tagINCHI_Output {
    fn default() -> Self {
        Self {
            szInChI: ::std::default::Default::default(),
            szAuxInfo: ::std::default::Default::default(),
            szMessage: ::std::default::Default::default(),
            szLog: ::std::default::Default::default(),
        }
    }
}
pub type inchi_Output = tagINCHI_Output;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHI_OutputStruct {
    pub atom: SourceMutPointer<inchi_Atom>,
    pub stereo0D: SourceMutPointer<inchi_Stereo0D>,
    pub num_atoms: AT_NUM,
    pub num_stereo0D: AT_NUM,
    pub szMessage: SourceMutPointer<i8>,
    pub szLog: SourceMutPointer<i8>,
    pub WarningFlags: [[u64; 2usize]; 2usize],
}
impl ::std::default::Default for tagINCHI_OutputStruct {
    fn default() -> Self {
        Self {
            atom: ::std::default::Default::default(),
            stereo0D: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_stereo0D: ::std::default::Default::default(),
            szMessage: ::std::default::Default::default(),
            szLog: ::std::default::Default::default(),
            WarningFlags: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
        }
    }
}
pub type inchi_OutputStruct = tagINCHI_OutputStruct;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHI_OutputStructEx {
    pub atom: SourceMutPointer<inchi_Atom>,
    pub stereo0D: SourceMutPointer<inchi_Stereo0D>,
    pub num_atoms: AT_NUM,
    pub num_stereo0D: AT_NUM,
    pub szMessage: SourceMutPointer<i8>,
    pub szLog: SourceMutPointer<i8>,
    pub WarningFlags: [[u64; 2usize]; 2usize],
    pub polymer: SourceMutPointer<inchi_Output_Polymer>,
    pub v3000: SourceMutPointer<inchi_Output_V3000>,
}
impl ::std::default::Default for tagINCHI_OutputStructEx {
    fn default() -> Self {
        Self {
            atom: ::std::default::Default::default(),
            stereo0D: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_stereo0D: ::std::default::Default::default(),
            szMessage: ::std::default::Default::default(),
            szLog: ::std::default::Default::default(),
            WarningFlags: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            polymer: ::std::default::Default::default(),
            v3000: ::std::default::Default::default(),
        }
    }
}
pub type inchi_OutputStructEx = tagINCHI_OutputStructEx;
pub const tagRetValGetINCHI_inchi_Ret_BREAK: tagRetValGetINCHI = -100;
pub const tagRetValGetINCHI_inchi_Ret_SKIP: tagRetValGetINCHI = -2;
pub const tagRetValGetINCHI_inchi_Ret_EOF: tagRetValGetINCHI = -1;
pub const tagRetValGetINCHI_inchi_Ret_OKAY: tagRetValGetINCHI = 0;
pub const tagRetValGetINCHI_inchi_Ret_WARNING: tagRetValGetINCHI = 1;
pub const tagRetValGetINCHI_inchi_Ret_ERROR: tagRetValGetINCHI = 2;
pub const tagRetValGetINCHI_inchi_Ret_FATAL: tagRetValGetINCHI = 3;
pub const tagRetValGetINCHI_inchi_Ret_UNKNOWN: tagRetValGetINCHI = 4;
pub const tagRetValGetINCHI_inchi_Ret_BUSY: tagRetValGetINCHI = 5;
pub type tagRetValGetINCHI = i32;
pub use self::tagRetValGetINCHI as RetValGetINCHI;
pub const tagRetValMOL2INCHI_mol2inchi_Ret_OKAY: tagRetValMOL2INCHI = 0;
pub const tagRetValMOL2INCHI_mol2inchi_Ret_WARNING: tagRetValMOL2INCHI = 1;
pub const tagRetValMOL2INCHI_mol2inchi_Ret_EOF: tagRetValMOL2INCHI = -1;
pub const tagRetValMOL2INCHI_mol2inchi_Ret_ERROR: tagRetValMOL2INCHI = 2;
pub const tagRetValMOL2INCHI_mol2inchi_Ret_ERROR_get: tagRetValMOL2INCHI = 4;
pub const tagRetValMOL2INCHI_mol2inchi_Ret_ERROR_comp: tagRetValMOL2INCHI = 5;
pub type tagRetValMOL2INCHI = i32;
pub use self::tagRetValMOL2INCHI as RetValMol2INCHI;
pub const tagRetValCheckINCHI_INCHI_VALID_STANDARD: tagRetValCheckINCHI = 0;
pub const tagRetValCheckINCHI_INCHI_VALID_NON_STANDARD: tagRetValCheckINCHI = 1;
pub const tagRetValCheckINCHI_INCHI_VALID_BETA: tagRetValCheckINCHI = 2;
pub const tagRetValCheckINCHI_INCHI_INVALID_PREFIX: tagRetValCheckINCHI = 3;
pub const tagRetValCheckINCHI_INCHI_INVALID_VERSION: tagRetValCheckINCHI = 4;
pub const tagRetValCheckINCHI_INCHI_INVALID_LAYOUT: tagRetValCheckINCHI = 5;
pub const tagRetValCheckINCHI_INCHI_FAIL_I2I: tagRetValCheckINCHI = 6;
pub type tagRetValCheckINCHI = u32;
pub use self::tagRetValCheckINCHI as RetValCheckINCHI;
#[derive(Clone, Debug, PartialEq)]
pub struct tagInchiInpData {
    pub pInp: SourceMutPointer<inchi_Input>,
    pub bChiral: i32,
    pub szErrMsg: [i8; 256usize],
}
impl ::std::default::Default for tagInchiInpData {
    fn default() -> Self {
        Self {
            pInp: ::std::default::Default::default(),
            bChiral: ::std::default::Default::default(),
            szErrMsg: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type InchiInpData = tagInchiInpData;
pub const tagRetValGetINCHIKey_INCHIKEY_VALID_STANDARD: tagRetValGetINCHIKey = 0;
pub const tagRetValGetINCHIKey_INCHIKEY_VALID_NON_STANDARD: tagRetValGetINCHIKey = -1;
pub const tagRetValGetINCHIKey_INCHIKEY_INVALID_LENGTH: tagRetValGetINCHIKey = 1;
pub const tagRetValGetINCHIKey_INCHIKEY_INVALID_LAYOUT: tagRetValGetINCHIKey = 2;
pub const tagRetValGetINCHIKey_INCHIKEY_INVALID_VERSION: tagRetValGetINCHIKey = 3;
pub type tagRetValGetINCHIKey = i32;
pub use self::tagRetValGetINCHIKey as RetValCheckINCHIKeyv;
pub type AT_NUMBR = u16;
pub type NUM_HS = i16;
pub type INCHI_MODES = u64;
#[derive(Clone, Debug, PartialEq)]
pub struct tagNormAtom {
    pub elname: [i8; 6usize],
    pub el_number: U_CHAR,
    pub neighbor: [AT_NUMBR; 20usize],
    pub orig_at_number: AT_NUMBR,
    pub orig_compt_at_numb: AT_NUMBR,
    pub bond_stereo: [S_CHAR; 20usize],
    pub bond_type: [U_CHAR; 20usize],
    pub valence: S_CHAR,
    pub chem_bonds_valence: S_CHAR,
    pub num_H: S_CHAR,
    pub num_iso_H: [S_CHAR; 3usize],
    pub iso_atw_diff: S_CHAR,
    pub charge: S_CHAR,
    pub radical: S_CHAR,
    pub bAmbiguousStereo: S_CHAR,
    pub cFlags: S_CHAR,
    pub at_type: AT_NUMBR,
    pub component: AT_NUMBR,
    pub endpoint: AT_NUMBR,
    pub c_point: AT_NUMBR,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub bUsed0DParity: S_CHAR,
    pub p_parity: S_CHAR,
    pub p_orig_at_num: [AT_NUMBR; 4usize],
    pub sb_ord: [S_CHAR; 3usize],
    pub sn_ord: [S_CHAR; 3usize],
    pub sb_parity: [S_CHAR; 3usize],
    pub sn_orig_at_num: [AT_NUMBR; 3usize],
    pub bCutVertex: S_CHAR,
    pub nRingSystem: AT_NUMBR,
    pub nNumAtInRingSystem: AT_NUMBR,
    pub nBlockSystem: AT_NUMBR,
}
impl ::std::default::Default for tagNormAtom {
    fn default() -> Self {
        Self {
            elname: ::std::array::from_fn(|_| ::std::default::Default::default()),
            el_number: ::std::default::Default::default(),
            neighbor: ::std::array::from_fn(|_| ::std::default::Default::default()),
            orig_at_number: ::std::default::Default::default(),
            orig_compt_at_numb: ::std::default::Default::default(),
            bond_stereo: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bond_type: ::std::array::from_fn(|_| ::std::default::Default::default()),
            valence: ::std::default::Default::default(),
            chem_bonds_valence: ::std::default::Default::default(),
            num_H: ::std::default::Default::default(),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            iso_atw_diff: ::std::default::Default::default(),
            charge: ::std::default::Default::default(),
            radical: ::std::default::Default::default(),
            bAmbiguousStereo: ::std::default::Default::default(),
            cFlags: ::std::default::Default::default(),
            at_type: ::std::default::Default::default(),
            component: ::std::default::Default::default(),
            endpoint: ::std::default::Default::default(),
            c_point: ::std::default::Default::default(),
            x: ::std::default::Default::default(),
            y: ::std::default::Default::default(),
            z: ::std::default::Default::default(),
            bUsed0DParity: ::std::default::Default::default(),
            p_parity: ::std::default::Default::default(),
            p_orig_at_num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sb_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sn_ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sb_parity: ::std::array::from_fn(|_| ::std::default::Default::default()),
            sn_orig_at_num: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bCutVertex: ::std::default::Default::default(),
            nRingSystem: ::std::default::Default::default(),
            nNumAtInRingSystem: ::std::default::Default::default(),
            nBlockSystem: ::std::default::Default::default(),
        }
    }
}
pub type NORM_ATOM = tagNormAtom;
#[derive(Clone, Debug, PartialEq)]
pub struct tagNormAtomData {
    pub at: SourceMutPointer<NORM_ATOM>,
    pub at_fixed_bonds: SourceMutPointer<NORM_ATOM>,
    pub num_at: i32,
    pub num_removed_H: i32,
    pub num_bonds: i32,
    pub num_isotopic: i32,
    pub bExists: i32,
    pub bDeleted: i32,
    pub bHasIsotopicLayer: i32,
    pub bTautomeric: i32,
    pub bTautPreprocessed: i32,
    pub nNumRemovedProtons: i32,
    pub nNumRemovedProtonsIsotopic: [NUM_HS; 3usize],
    pub num_iso_H: [NUM_HS; 3usize],
    pub bTautFlags: INCHI_MODES,
    pub bTautFlagsDone: INCHI_MODES,
    pub bNormalizationFlags: INCHI_MODES,
}
impl ::std::default::Default for tagNormAtomData {
    fn default() -> Self {
        Self {
            at: ::std::default::Default::default(),
            at_fixed_bonds: ::std::default::Default::default(),
            num_at: ::std::default::Default::default(),
            num_removed_H: ::std::default::Default::default(),
            num_bonds: ::std::default::Default::default(),
            num_isotopic: ::std::default::Default::default(),
            bExists: ::std::default::Default::default(),
            bDeleted: ::std::default::Default::default(),
            bHasIsotopicLayer: ::std::default::Default::default(),
            bTautomeric: ::std::default::Default::default(),
            bTautPreprocessed: ::std::default::Default::default(),
            nNumRemovedProtons: ::std::default::Default::default(),
            nNumRemovedProtonsIsotopic: ::std::array::from_fn(|_| {
                ::std::default::Default::default()
            }),
            num_iso_H: ::std::array::from_fn(|_| ::std::default::Default::default()),
            bTautFlags: ::std::default::Default::default(),
            bTautFlagsDone: ::std::default::Default::default(),
            bNormalizationFlags: ::std::default::Default::default(),
        }
    }
}
pub type NORM_ATOMS = tagNormAtomData;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHIGEN_DATA {
    pub pStrErrStruct: [i8; 256usize],
    pub num_components: [i32; 2usize],
    pub NormAtomsNontaut: [SourceMutPointer<NORM_ATOMS>; 2usize],
    pub NormAtomsTaut: [SourceMutPointer<NORM_ATOMS>; 2usize],
}
impl ::std::default::Default for tagINCHIGEN_DATA {
    fn default() -> Self {
        Self {
            pStrErrStruct: ::std::array::from_fn(|_| ::std::default::Default::default()),
            num_components: ::std::array::from_fn(|_| ::std::default::Default::default()),
            NormAtomsNontaut: ::std::array::from_fn(|_| ::std::default::Default::default()),
            NormAtomsTaut: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type INCHIGEN_DATA = tagINCHIGEN_DATA;
pub type INCHIGEN_HANDLE = SourceMutPointer<SourceVoid>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_STATUS_HANDLE_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_STATUS_HANDLE_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_STATUS_HANDLE = SourceMutPointer<IXA_STATUS_HANDLE_STRUCT>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_MOL_HANDLE_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_MOL_HANDLE_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_MOL_HANDLE = SourceMutPointer<IXA_MOL_HANDLE_STRUCT>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_INCHIBUILDER_HANDLE_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_INCHIBUILDER_HANDLE_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_INCHIBUILDER_HANDLE = SourceMutPointer<IXA_INCHIBUILDER_HANDLE_STRUCT>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_INCHIKEYBUILDER_HANDLE_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_INCHIKEYBUILDER_HANDLE_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_INCHIKEYBUILDER_HANDLE = SourceMutPointer<IXA_INCHIKEYBUILDER_HANDLE_STRUCT>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_ATOMID_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_ATOMID_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_ATOMID = SourceMutPointer<IXA_ATOMID_STRUCT>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_BONDID_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_BONDID_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_BONDID = SourceMutPointer<IXA_BONDID_STRUCT>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_STEREOID_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_STEREOID_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_STEREOID = SourceMutPointer<IXA_STEREOID_STRUCT>;
#[derive(Clone, Debug, PartialEq)]
pub struct IXA_POLYMERUNITID_STRUCT {
    pub dummy: i32,
}
impl ::std::default::Default for IXA_POLYMERUNITID_STRUCT {
    fn default() -> Self {
        Self {
            dummy: ::std::default::Default::default(),
        }
    }
}
pub type IXA_POLYMERUNITID = SourceMutPointer<IXA_POLYMERUNITID_STRUCT>;
pub const IXA_STATUS_IXA_STATUS_SUCCESS: IXA_STATUS = 0;
pub const IXA_STATUS_IXA_STATUS_WARNING: IXA_STATUS = 1;
pub const IXA_STATUS_IXA_STATUS_ERROR: IXA_STATUS = 2;
pub type IXA_STATUS = u32;
pub const IXA_BOOL_IXA_FALSE: IXA_BOOL = 0;
pub const IXA_BOOL_IXA_TRUE: IXA_BOOL = 1;
pub type IXA_BOOL = u32;
pub const IXA_ATOM_RADICAL_IXA_ATOM_RADICAL_NONE: IXA_ATOM_RADICAL = 0;
pub const IXA_ATOM_RADICAL_IXA_ATOM_RADICAL_SINGLET: IXA_ATOM_RADICAL = 1;
pub const IXA_ATOM_RADICAL_IXA_ATOM_RADICAL_DOUBLET: IXA_ATOM_RADICAL = 2;
pub const IXA_ATOM_RADICAL_IXA_ATOM_RADICAL_TRIPLET: IXA_ATOM_RADICAL = 3;
pub type IXA_ATOM_RADICAL = u32;
pub const IXA_BOND_TYPE_IXA_BOND_TYPE_SINGLE: IXA_BOND_TYPE = 1;
pub const IXA_BOND_TYPE_IXA_BOND_TYPE_DOUBLE: IXA_BOND_TYPE = 2;
pub const IXA_BOND_TYPE_IXA_BOND_TYPE_TRIPLE: IXA_BOND_TYPE = 3;
pub const IXA_BOND_TYPE_IXA_BOND_TYPE_AROMATIC: IXA_BOND_TYPE = 4;
pub type IXA_BOND_TYPE = u32;
pub const IXA_BOND_WEDGE_IXA_BOND_WEDGE_NONE: IXA_BOND_WEDGE = 0;
pub const IXA_BOND_WEDGE_IXA_BOND_WEDGE_UP: IXA_BOND_WEDGE = 1;
pub const IXA_BOND_WEDGE_IXA_BOND_WEDGE_DOWN: IXA_BOND_WEDGE = 2;
pub const IXA_BOND_WEDGE_IXA_BOND_WEDGE_EITHER: IXA_BOND_WEDGE = 3;
pub type IXA_BOND_WEDGE = u32;
pub const IXA_DBLBOND_CONFIG_IXA_DBLBOND_CONFIG_PERCEIVE: IXA_DBLBOND_CONFIG = 0;
pub const IXA_DBLBOND_CONFIG_IXA_DBLBOND_CONFIG_EITHER: IXA_DBLBOND_CONFIG = 1;
pub type IXA_DBLBOND_CONFIG = u32;
pub const IXA_STEREO_TOPOLOGY_IXA_STEREO_TOPOLOGY_INVALID: IXA_STEREO_TOPOLOGY = 0;
pub const IXA_STEREO_TOPOLOGY_IXA_STEREO_TOPOLOGY_TETRAHEDRON: IXA_STEREO_TOPOLOGY = 2;
pub const IXA_STEREO_TOPOLOGY_IXA_STEREO_TOPOLOGY_RECTANGLE: IXA_STEREO_TOPOLOGY = 3;
pub const IXA_STEREO_TOPOLOGY_IXA_STEREO_TOPOLOGY_ANTIRECTANGLE: IXA_STEREO_TOPOLOGY = 4;
pub type IXA_STEREO_TOPOLOGY = u32;
pub const IXA_STEREO_PARITY_IXA_STEREO_PARITY_NONE: IXA_STEREO_PARITY = 0;
pub const IXA_STEREO_PARITY_IXA_STEREO_PARITY_ODD: IXA_STEREO_PARITY = 1;
pub const IXA_STEREO_PARITY_IXA_STEREO_PARITY_EVEN: IXA_STEREO_PARITY = 2;
pub const IXA_STEREO_PARITY_IXA_STEREO_PARITY_UNKNOWN: IXA_STEREO_PARITY = 3;
pub type IXA_STEREO_PARITY = u32;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_NewPsOff: IXA_INCHIBUILDER_OPTION = 0;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_DoNotAddH: IXA_INCHIBUILDER_OPTION = 1;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_SUU: IXA_INCHIBUILDER_OPTION = 2;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_SLUUD: IXA_INCHIBUILDER_OPTION = 3;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_FixedH: IXA_INCHIBUILDER_OPTION = 4;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_RecMet: IXA_INCHIBUILDER_OPTION = 5;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_KET: IXA_INCHIBUILDER_OPTION = 6;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_15T: IXA_INCHIBUILDER_OPTION = 7;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_SaveOpt: IXA_INCHIBUILDER_OPTION = 8;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_AuxNone: IXA_INCHIBUILDER_OPTION = 9;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_WarnOnEmptyStructure:
    IXA_INCHIBUILDER_OPTION = 10;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_LargeMolecules: IXA_INCHIBUILDER_OPTION =
    11;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_Polymers: IXA_INCHIBUILDER_OPTION = 12;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_Polymers105: IXA_INCHIBUILDER_OPTION = 13;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_Polymers105Plus: IXA_INCHIBUILDER_OPTION =
    14;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_FilterSS: IXA_INCHIBUILDER_OPTION = 15;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_InvFilterSS: IXA_INCHIBUILDER_OPTION = 16;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_NPZZ: IXA_INCHIBUILDER_OPTION = 17;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_SATZZ: IXA_INCHIBUILDER_OPTION = 18;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_NoFrameShift: IXA_INCHIBUILDER_OPTION =
    19;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_FoldCRU: IXA_INCHIBUILDER_OPTION = 20;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_NoEdits: IXA_INCHIBUILDER_OPTION = 21;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_LooseTSACheck: IXA_INCHIBUILDER_OPTION =
    22;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_OutErrInChI: IXA_INCHIBUILDER_OPTION = 23;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_NoWarnings: IXA_INCHIBUILDER_OPTION = 24;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_DoDrv: IXA_INCHIBUILDER_OPTION = 25;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_DoDrvReport: IXA_INCHIBUILDER_OPTION = 26;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_DoR2C: IXA_INCHIBUILDER_OPTION = 27;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_DoneOnly: IXA_INCHIBUILDER_OPTION = 28;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_OnlyRecSalt: IXA_INCHIBUILDER_OPTION = 29;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_OnlyExact: IXA_INCHIBUILDER_OPTION = 30;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_OnlyRecMet: IXA_INCHIBUILDER_OPTION = 31;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_PT_22_00: IXA_INCHIBUILDER_OPTION = 32;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_PT_16_00: IXA_INCHIBUILDER_OPTION = 33;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_PT_06_00: IXA_INCHIBUILDER_OPTION = 34;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_PT_39_00: IXA_INCHIBUILDER_OPTION = 35;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_PT_13_00: IXA_INCHIBUILDER_OPTION = 36;
pub const IXA_INCHIBUILDER_OPTION_IXA_INCHIBUILDER_OPTION_PT_18_00: IXA_INCHIBUILDER_OPTION = 37;
pub type IXA_INCHIBUILDER_OPTION = u32;
pub const IXA_INCHIBUILDER_STEREOOPTION_IXA_INCHIBUILDER_STEREOOPTION_SAbs:
    IXA_INCHIBUILDER_STEREOOPTION = 0;
pub const IXA_INCHIBUILDER_STEREOOPTION_IXA_INCHIBUILDER_STEREOOPTION_SNon:
    IXA_INCHIBUILDER_STEREOOPTION = 1;
pub const IXA_INCHIBUILDER_STEREOOPTION_IXA_INCHIBUILDER_STEREOOPTION_SRel:
    IXA_INCHIBUILDER_STEREOOPTION = 2;
pub const IXA_INCHIBUILDER_STEREOOPTION_IXA_INCHIBUILDER_STEREOOPTION_SRac:
    IXA_INCHIBUILDER_STEREOOPTION = 3;
pub const IXA_INCHIBUILDER_STEREOOPTION_IXA_INCHIBUILDER_STEREOOPTION_SUCF:
    IXA_INCHIBUILDER_STEREOOPTION = 4;
pub type IXA_INCHIBUILDER_STEREOOPTION = u32;
pub const tagInchiCompareDiffBits_INCHIDIFF_ZERO: tagInchiCompareDiffBits = 0;
pub const tagInchiCompareDiffBits_INCHIDIFF_PROBLEM: tagInchiCompareDiffBits = 1;
pub const tagInchiCompareDiffBits_INCHIDIFF_NUM_AT: tagInchiCompareDiffBits = 1;
pub const tagInchiCompareDiffBits_INCHIDIFF_ATOMS: tagInchiCompareDiffBits = 1;
pub const tagInchiCompareDiffBits_INCHIDIFF_NUM_EL: tagInchiCompareDiffBits = 1;
pub const tagInchiCompareDiffBits_INCHIDIFF_CON_LEN: tagInchiCompareDiffBits = 1;
pub const tagInchiCompareDiffBits_INCHIDIFF_CON_TBL: tagInchiCompareDiffBits = 1;
pub const tagInchiCompareDiffBits_INCHIDIFF_POSITION_H: tagInchiCompareDiffBits = 2;
pub const tagInchiCompareDiffBits_INCHIDIFF_MORE_FH: tagInchiCompareDiffBits = 4;
pub const tagInchiCompareDiffBits_INCHIDIFF_LESS_FH: tagInchiCompareDiffBits = 4;
pub const tagInchiCompareDiffBits_INCHIDIFF_MORE_H: tagInchiCompareDiffBits = 8;
pub const tagInchiCompareDiffBits_INCHIDIFF_LESS_H: tagInchiCompareDiffBits = 8;
pub const tagInchiCompareDiffBits_INCHIDIFF_NO_TAUT: tagInchiCompareDiffBits = 16;
pub const tagInchiCompareDiffBits_INCHIDIFF_WRONG_TAUT: tagInchiCompareDiffBits = 32;
pub const tagInchiCompareDiffBits_INCHIDIFF_SINGLE_TG: tagInchiCompareDiffBits = 64;
pub const tagInchiCompareDiffBits_INCHIDIFF_MULTIPLE_TG: tagInchiCompareDiffBits = 128;
pub const tagInchiCompareDiffBits_INCHIDIFF_EXTRA_TG_ENDP: tagInchiCompareDiffBits = 256;
pub const tagInchiCompareDiffBits_INCHIDIFF_MISS_TG_ENDP: tagInchiCompareDiffBits = 256;
pub const tagInchiCompareDiffBits_INCHIDIFF_DIFF_TG_ENDP: tagInchiCompareDiffBits = 256;
pub const tagInchiCompareDiffBits_INCHIDIFF_NUM_TG: tagInchiCompareDiffBits = 512;
pub const tagInchiCompareDiffBits_INCHIDIFF_TG: tagInchiCompareDiffBits = 512;
pub const tagInchiCompareDiffBits_INCHIDIFF_NUM_ISO_AT: tagInchiCompareDiffBits = 1024;
pub const tagInchiCompareDiffBits_INCHIDIFF_ISO_AT: tagInchiCompareDiffBits = 1024;
pub const tagInchiCompareDiffBits_INCHIDIFF_REM_ISO_H: tagInchiCompareDiffBits = 2048;
pub const tagInchiCompareDiffBits_INCHIDIFF_MOB_ISO_H: tagInchiCompareDiffBits = 4096;
pub const tagInchiCompareDiffBits_INCHIDIFF_CHARGE: tagInchiCompareDiffBits = 8192;
pub const tagInchiCompareDiffBits_INCHIDIFF_REM_PROT: tagInchiCompareDiffBits = 16384;
pub const tagInchiCompareDiffBits_INCHIDIFF_MOBH_PROTONS: tagInchiCompareDiffBits = 32768;
pub const tagInchiCompareDiffBits_INCHIDIFF_SC_INV: tagInchiCompareDiffBits = 65536;
pub const tagInchiCompareDiffBits_INCHIDIFF_SC_PARITY: tagInchiCompareDiffBits = 131072;
pub const tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA_UNDF: tagInchiCompareDiffBits = 262144;
pub const tagInchiCompareDiffBits_INCHIDIFF_SC_EXTRA: tagInchiCompareDiffBits = 524288;
pub const tagInchiCompareDiffBits_INCHIDIFF_SC_MISS_UNDF: tagInchiCompareDiffBits = 1048576;
pub const tagInchiCompareDiffBits_INCHIDIFF_SC_MISS: tagInchiCompareDiffBits = 2097152;
pub const tagInchiCompareDiffBits_INCHIDIFF_SB_PARITY: tagInchiCompareDiffBits = 4194304;
pub const tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA_UNDF: tagInchiCompareDiffBits = 8388608;
pub const tagInchiCompareDiffBits_INCHIDIFF_SB_EXTRA: tagInchiCompareDiffBits = 16777216;
pub const tagInchiCompareDiffBits_INCHIDIFF_SB_MISS_UNDF: tagInchiCompareDiffBits = 33554432;
pub const tagInchiCompareDiffBits_INCHIDIFF_SB_MISS: tagInchiCompareDiffBits = 67108864;
pub const tagInchiCompareDiffBits_INCHIDIFF_COMP_HLAYER: tagInchiCompareDiffBits = 134217728;
pub const tagInchiCompareDiffBits_INCHIDIFF_COMP_NUMBER: tagInchiCompareDiffBits = 268435456;
pub const tagInchiCompareDiffBits_INCHIDIFF_STR2INCHI_ERR: tagInchiCompareDiffBits = 536870912;
pub type tagInchiCompareDiffBits = u32;
pub use self::tagInchiCompareDiffBits as INCHIDIFF;
pub const tagtagCompareInchiMsgGroupID_IDGRP_ZERO: tagtagCompareInchiMsgGroupID = 0;
pub const tagtagCompareInchiMsgGroupID_IDGRP_ERR: tagtagCompareInchiMsgGroupID = 1;
pub const tagtagCompareInchiMsgGroupID_IDGRP_H: tagtagCompareInchiMsgGroupID = 2;
pub const tagtagCompareInchiMsgGroupID_IDGRP_MOB_GRP: tagtagCompareInchiMsgGroupID = 3;
pub const tagtagCompareInchiMsgGroupID_IDGRP_ISO_AT: tagtagCompareInchiMsgGroupID = 4;
pub const tagtagCompareInchiMsgGroupID_IDGRP_CHARGE: tagtagCompareInchiMsgGroupID = 5;
pub const tagtagCompareInchiMsgGroupID_IDGRP_PROTONS: tagtagCompareInchiMsgGroupID = 6;
pub const tagtagCompareInchiMsgGroupID_IDGRP_ISO_H: tagtagCompareInchiMsgGroupID = 7;
pub const tagtagCompareInchiMsgGroupID_IDGRP_SC: tagtagCompareInchiMsgGroupID = 8;
pub const tagtagCompareInchiMsgGroupID_IDGRP_SB: tagtagCompareInchiMsgGroupID = 9;
pub const tagtagCompareInchiMsgGroupID_IDGRP_HLAYER: tagtagCompareInchiMsgGroupID = 10;
pub const tagtagCompareInchiMsgGroupID_IDGRP_COMP: tagtagCompareInchiMsgGroupID = 11;
pub const tagtagCompareInchiMsgGroupID_IDGRP_CONV_ERR: tagtagCompareInchiMsgGroupID = 12;
pub type tagtagCompareInchiMsgGroupID = u32;
pub use self::tagtagCompareInchiMsgGroupID as CMP_INCHI_MSG_GROUP_ID;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCompareInchiMsg {
    pub nBit: INCHIDIFF,
    pub nGroupID: CMP_INCHI_MSG_GROUP_ID,
    pub szMsg: SourceConstPointer<i8>,
}
impl ::std::default::Default for tagCompareInchiMsg {
    fn default() -> Self {
        Self {
            nBit: ::std::default::Default::default(),
            nGroupID: ::std::default::Default::default(),
            szMsg: ::std::default::Default::default(),
        }
    }
}
pub type CMP_INCHI_MSG = tagCompareInchiMsg;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCompareInchiMsgGroup {
    pub nGroupID: CMP_INCHI_MSG_GROUP_ID,
    pub szGroupName: SourceConstPointer<i8>,
}
impl ::std::default::Default for tagCompareInchiMsgGroup {
    fn default() -> Self {
        Self {
            nGroupID: ::std::default::Default::default(),
            szGroupName: ::std::default::Default::default(),
        }
    }
}
pub type CMP_INCHI_MSG_GROUP = tagCompareInchiMsgGroup;
#[derive(Clone, Debug, PartialEq)]
pub struct ReadINCHI_CtlData {
    pub ulongID: u64,
    pub bTooLongLine: i32,
    pub bHeaderRead: i32,
    pub bErrorMsg: i32,
    pub bRestoreInfo: i32,
}
impl ::std::default::Default for ReadINCHI_CtlData {
    fn default() -> Self {
        Self {
            ulongID: ::std::default::Default::default(),
            bTooLongLine: ::std::default::Default::default(),
            bHeaderRead: ::std::default::Default::default(),
            bErrorMsg: ::std::default::Default::default(),
            bRestoreInfo: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct sha2_context {
    pub total: [u64; 2usize],
    pub state: [u64; 8usize],
    pub buffer: [u8; 64usize],
}
impl ::std::default::Default for sha2_context {
    fn default() -> Self {
        Self {
            total: ::std::array::from_fn(|_| ::std::default::Default::default()),
            state: ::std::array::from_fn(|_| ::std::default::Default::default()),
            buffer: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type va_list = SourceVaList;
pub type STBSP_SPRINTFCB = ::std::option::Option<
    fn(
        buf: SourceConstPointer<i8>,
        user: SourceMutPointer<SourceVoid>,
        len: i32,
    ) -> SourceMutPointer<i8>,
>;
#[derive(Clone, Debug, PartialEq)]
pub struct tagCOMPONENT_TREAT_INFO {
    pub n1: i32,
    pub n2: i32,
    pub num_atoms: i32,
    pub num_at_tg: i32,
    pub num_deleted_H: i32,
    pub num_deleted_H_taut: i32,
    pub nMode: INCHI_MODE,
    pub vt_group_info: T_GROUP_INFO,
    pub vt_group_info_orig: T_GROUP_INFO,
    pub s: [ATOM_SIZES; 2usize],
    pub Bcn: BCN,
    pub bHasIsotopicAtoms: i32,
    pub bMayHaveStereo: i32,
    pub num_taut_at: i32,
    pub bPointedEdgeStereo: i32,
    pub vABParityUnknown: i32,
    pub bLooseTSACheck: i32,
    pub bStereoAtZz: i32,
    pub bTautFlags: INCHI_MODE,
    pub bTautFlagsDone: INCHI_MODE,
    pub nUserMode: INCHI_MODE,
    pub at: [SourceMutPointer<sp_ATOM>; 2usize],
    pub out_at: SourceMutPointer<inp_ATOM>,
    pub fix_isofixedh: i32,
    pub fix_termhchrg: i32,
}
impl ::std::default::Default for tagCOMPONENT_TREAT_INFO {
    fn default() -> Self {
        Self {
            n1: ::std::default::Default::default(),
            n2: ::std::default::Default::default(),
            num_atoms: ::std::default::Default::default(),
            num_at_tg: ::std::default::Default::default(),
            num_deleted_H: ::std::default::Default::default(),
            num_deleted_H_taut: ::std::default::Default::default(),
            nMode: ::std::default::Default::default(),
            vt_group_info: ::std::default::Default::default(),
            vt_group_info_orig: ::std::default::Default::default(),
            s: ::std::array::from_fn(|_| ::std::default::Default::default()),
            Bcn: ::std::default::Default::default(),
            bHasIsotopicAtoms: ::std::default::Default::default(),
            bMayHaveStereo: ::std::default::Default::default(),
            num_taut_at: ::std::default::Default::default(),
            bPointedEdgeStereo: ::std::default::Default::default(),
            vABParityUnknown: ::std::default::Default::default(),
            bLooseTSACheck: ::std::default::Default::default(),
            bStereoAtZz: ::std::default::Default::default(),
            bTautFlags: ::std::default::Default::default(),
            bTautFlagsDone: ::std::default::Default::default(),
            nUserMode: ::std::default::Default::default(),
            at: ::std::array::from_fn(|_| ::std::default::Default::default()),
            out_at: ::std::default::Default::default(),
            fix_isofixedh: ::std::default::Default::default(),
            fix_termhchrg: ::std::default::Default::default(),
        }
    }
}
pub type COMPONENT_TREAT_INFO = tagCOMPONENT_TREAT_INFO;
#[derive(Clone, Debug, PartialEq)]
pub struct tagINCHIGEN_CONTROL {
    pub init_passed: i32,
    pub norm_passed: i32,
    pub canon_passed: i32,
    pub InpParms: INPUT_PARMS,
    pub ulTotalProcessingTime: u64,
    pub szTitle: [i8; 575usize],
    pub strbuf_container: INCHI_IOS_STRING,
    pub num_err: i64,
    pub num_inp: i64,
    pub OrigStruct: ORIG_STRUCT,
    pub OrigInpData: ORIG_ATOM_DATA,
    pub StructData: STRUCT_DATA,
    pub PrepInpData: [ORIG_ATOM_DATA; 2usize],
    pub InpCurAtData: [SourceMutPointer<INP_ATOM_DATA>; 2usize],
    pub InpNormAtData: [SourceMutPointer<INP_ATOM_DATA>; 2usize],
    pub InpNormTautData: [SourceMutPointer<INP_ATOM_DATA>; 2usize],
    pub composite_norm_data: [[COMP_ATOM_DATA; 3usize]; 2usize],
    pub ncFlags: NORM_CANON_FLAGS,
    pub pINChI: [SourceMutPointer<PINChI2>; 2usize],
    pub pINChI_Aux: [SourceMutPointer<PINChI_Aux2>; 2usize],
    pub cti: [SourceMutPointer<COMPONENT_TREAT_INFO>; 2usize],
    pub inchi_file: [INCHI_IOSTREAM; 3usize],
}
impl ::std::default::Default for tagINCHIGEN_CONTROL {
    fn default() -> Self {
        Self {
            init_passed: ::std::default::Default::default(),
            norm_passed: ::std::default::Default::default(),
            canon_passed: ::std::default::Default::default(),
            InpParms: ::std::default::Default::default(),
            ulTotalProcessingTime: ::std::default::Default::default(),
            szTitle: ::std::array::from_fn(|_| ::std::default::Default::default()),
            strbuf_container: ::std::default::Default::default(),
            num_err: ::std::default::Default::default(),
            num_inp: ::std::default::Default::default(),
            OrigStruct: ::std::default::Default::default(),
            OrigInpData: ::std::default::Default::default(),
            StructData: ::std::default::Default::default(),
            PrepInpData: ::std::array::from_fn(|_| ::std::default::Default::default()),
            InpCurAtData: ::std::array::from_fn(|_| ::std::default::Default::default()),
            InpNormAtData: ::std::array::from_fn(|_| ::std::default::Default::default()),
            InpNormTautData: ::std::array::from_fn(|_| ::std::default::Default::default()),
            composite_norm_data: ::std::array::from_fn(|_| {
                ::std::array::from_fn(|_| ::std::default::Default::default())
            }),
            ncFlags: ::std::default::Default::default(),
            pINChI: ::std::array::from_fn(|_| ::std::default::Default::default()),
            pINChI_Aux: ::std::array::from_fn(|_| ::std::default::Default::default()),
            cti: ::std::array::from_fn(|_| ::std::default::Default::default()),
            inchi_file: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
pub type INCHIGEN_CONTROL = tagINCHIGEN_CONTROL;
#[derive(Clone, Debug, PartialEq)]
pub struct INCHIMOL_ATOM {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub atomic_number: i32,
    pub hydrogens: [i32; 4usize],
    pub mass: i32,
    pub radical: IXA_ATOM_RADICAL,
    pub charge: i32,
    pub bond_count: i32,
    pub bonds: [IXA_BONDID; 20usize],
}
impl ::std::default::Default for INCHIMOL_ATOM {
    fn default() -> Self {
        Self {
            x: ::std::default::Default::default(),
            y: ::std::default::Default::default(),
            z: ::std::default::Default::default(),
            atomic_number: ::std::default::Default::default(),
            hydrogens: ::std::array::from_fn(|_| ::std::default::Default::default()),
            mass: ::std::default::Default::default(),
            radical: ::std::default::Default::default(),
            charge: ::std::default::Default::default(),
            bond_count: ::std::default::Default::default(),
            bonds: ::std::array::from_fn(|_| ::std::default::Default::default()),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct INCHIMOL_BOND {
    pub atom1: IXA_ATOMID,
    pub atom2: IXA_ATOMID,
    pub type_: IXA_BOND_TYPE,
    pub config: IXA_DBLBOND_CONFIG,
    pub wedge_from_atom1: IXA_BOND_WEDGE,
    pub wedge_from_atom2: IXA_BOND_WEDGE,
}
impl ::std::default::Default for INCHIMOL_BOND {
    fn default() -> Self {
        Self {
            atom1: ::std::default::Default::default(),
            atom2: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            config: ::std::default::Default::default(),
            wedge_from_atom1: ::std::default::Default::default(),
            wedge_from_atom2: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct INCHIMOL_STEREO {
    pub topology: IXA_STEREO_TOPOLOGY,
    pub vertices: [IXA_ATOMID; 4usize],
    pub central_entity: SourceMutPointer<SourceVoid>,
    pub parity: IXA_STEREO_PARITY,
}
impl ::std::default::Default for INCHIMOL_STEREO {
    fn default() -> Self {
        Self {
            topology: ::std::default::Default::default(),
            vertices: ::std::array::from_fn(|_| ::std::default::Default::default()),
            central_entity: ::std::default::Default::default(),
            parity: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct INCHIMOL_SGROUP {
    pub id: i32,
    pub type_: i32,
    pub subtype: i32,
    pub conn: i32,
    pub label: i32,
    pub na: i32,
    pub nb: i32,
    pub xbr1: [f64; 4usize],
    pub xbr2: [f64; 4usize],
    pub smt: [i8; 80usize],
    pub alist: SourceMutPointer<i32>,
    pub blist: SourceMutPointer<i32>,
}
impl ::std::default::Default for INCHIMOL_SGROUP {
    fn default() -> Self {
        Self {
            id: ::std::default::Default::default(),
            type_: ::std::default::Default::default(),
            subtype: ::std::default::Default::default(),
            conn: ::std::default::Default::default(),
            label: ::std::default::Default::default(),
            na: ::std::default::Default::default(),
            nb: ::std::default::Default::default(),
            xbr1: ::std::array::from_fn(|_| ::std::default::Default::default()),
            xbr2: ::std::array::from_fn(|_| ::std::default::Default::default()),
            smt: ::std::array::from_fn(|_| ::std::default::Default::default()),
            alist: ::std::default::Default::default(),
            blist: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct INCHIMOL_POLYMER {
    pub units: SourceMutPointer<SourceMutPointer<INCHIMOL_SGROUP>>,
    pub n: i32,
}
impl ::std::default::Default for INCHIMOL_POLYMER {
    fn default() -> Self {
        Self {
            units: ::std::default::Default::default(),
            n: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct INCHIMOL_V3000 {
    pub n_non_star_atoms: i32,
    pub n_star_atoms: i32,
    pub atom_index_orig: SourceMutPointer<i32>,
    pub atom_index_fin: SourceMutPointer<i32>,
    pub n_sgroups: i32,
    pub n_3d_constraints: i32,
    pub n_collections: i32,
    pub n_non_haptic_bonds: i32,
    pub n_haptic_bonds: i32,
    pub lists_haptic_bonds: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_steabs: i32,
    pub lists_steabs: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_sterel: i32,
    pub lists_sterel: SourceMutPointer<SourceMutPointer<i32>>,
    pub n_sterac: i32,
    pub lists_sterac: SourceMutPointer<SourceMutPointer<i32>>,
}
impl ::std::default::Default for INCHIMOL_V3000 {
    fn default() -> Self {
        Self {
            n_non_star_atoms: ::std::default::Default::default(),
            n_star_atoms: ::std::default::Default::default(),
            atom_index_orig: ::std::default::Default::default(),
            atom_index_fin: ::std::default::Default::default(),
            n_sgroups: ::std::default::Default::default(),
            n_3d_constraints: ::std::default::Default::default(),
            n_collections: ::std::default::Default::default(),
            n_non_haptic_bonds: ::std::default::Default::default(),
            n_haptic_bonds: ::std::default::Default::default(),
            lists_haptic_bonds: ::std::default::Default::default(),
            n_steabs: ::std::default::Default::default(),
            lists_steabs: ::std::default::Default::default(),
            n_sterel: ::std::default::Default::default(),
            lists_sterel: ::std::default::Default::default(),
            n_sterac: ::std::default::Default::default(),
            lists_sterac: ::std::default::Default::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct INCHIMOL {
    pub atom_count: i32,
    pub atoms: SourceMutPointer<INCHIMOL_ATOM>,
    pub bond_count: i32,
    pub bonds: SourceMutPointer<INCHIMOL_BOND>,
    pub stereo_count: i32,
    pub stereos: SourceMutPointer<INCHIMOL_STEREO>,
    pub polymer: SourceMutPointer<INCHIMOL_POLYMER>,
    pub sgroup_count: i32,
    pub v3000: SourceMutPointer<INCHIMOL_V3000>,
    pub chiral: IXA_BOOL,
    pub reserved_atom_count: i32,
    pub reserved_bond_count: i32,
    pub reserved_stereo_count: i32,
    pub reserved_sgroup_count: i32,
}
impl ::std::default::Default for INCHIMOL {
    fn default() -> Self {
        Self {
            atom_count: ::std::default::Default::default(),
            atoms: ::std::default::Default::default(),
            bond_count: ::std::default::Default::default(),
            bonds: ::std::default::Default::default(),
            stereo_count: ::std::default::Default::default(),
            stereos: ::std::default::Default::default(),
            polymer: ::std::default::Default::default(),
            sgroup_count: ::std::default::Default::default(),
            v3000: ::std::default::Default::default(),
            chiral: ::std::default::Default::default(),
            reserved_atom_count: ::std::default::Default::default(),
            reserved_bond_count: ::std::default::Default::default(),
            reserved_stereo_count: ::std::default::Default::default(),
            reserved_sgroup_count: ::std::default::Default::default(),
        }
    }
}
pub(crate) mod local_ichi_bns {
    use super::*;
    pub const BNS_MARK_ONLY_BLOCKS: u32 = 1;
    pub const ALLOW_ONLY_SIMPLE_ALT_PATH: u32 = 0;
    pub const CHECK_TG_ALT_PATH: u32 = 0;
    pub const FIX_CPOINT_BOND_CAP: u32 = 1;
    pub const RESET_EDGE_FORBIDDEN_MASK: u32 = 1;
    pub const PR_SIMPLE_TYP: u32 = 2240;
    pub const KNOWN_ACIDIC_TYPE: u32 = 31;
    pub const ATTYP_OS: u32 = 543;
    pub const ATTYP_NP: u32 = 192;
    pub const ATTYP_N: u32 = 64;
    pub const ATTYP_P: u32 = 128;
    pub const AR_ANY_OH: u32 = 0;
    pub const AR_SIMPLE_STEPS: u32 = 3;
    pub const AR_SIMPLE_TYP1: u32 = 16415;
    pub const AA_ANY_O_Minus: u32 = 0;
    pub const AA_SIMPLE_STEPS: u32 = 3;
    pub const AA_SIMPLE_TYP1: u32 = 9247;
    pub const AA_SIMPLE_TYP4: u32 = 32768;
    pub const PR_HARD_TYP_POS: u32 = 64;
    pub const PR_HARD_TYP_POSP: u32 = 128;
    pub const PR_HARD_TYP_NEG: u32 = 607;
    pub const PR_HARD_TYP_H: u32 = 607;
    pub const AR_HARD_TYP_POS: u32 = 64;
    pub const AR_HARD_TYP_NEG: u32 = 607;
    pub const AR_HARD_TYP_HA: u32 = 17;
    pub const AR_HARD_TYP_HN: u32 = 607;
    pub const AA_HARD_TYP_POS: u32 = 64;
    pub const AA_HARD_TYP_NEG: u32 = 607;
    pub const AA_HARD_TYP_CO: u32 = 17;
    pub const AA_HARD_TYP_H: u32 = 607;
    pub const BNS_MAX_NUM_FLOW_CHANGES: u32 = 5;
    pub const TREE_NOT_IN_M: u32 = 0;
    pub const TREE_IN_2: u32 = 1;
    pub const TREE_IN_2BLOSS: u32 = 2;
    pub const TREE_IN_1: u32 = 3;
    pub const MAX_NEIGH: u32 = 6;
    pub const ALL_NONMETAL_Z: u32 = 0;
    pub const BNS_CHK_ALTP_NO_ALTPATH: u32 = 0;
    pub const BNS_CHK_ALTP_SAME_TGROUP: u32 = 1;
    pub const BNS_CHK_ALTP_SAME_VERTEX: u32 = 2;
    pub const BNS_CHK_ALTP_SET_SUCCESS: u32 = 4;
    pub const ALT_PATH_TAUTOM: u32 = 1;
    pub const ALT_PATH_CHARGE: u32 = 2;
    pub const ALT_PATH_4_SALT: u32 = 3;
    pub const MAX_NUM_RAD: u32 = 256;
    pub const CHECK_TACN: u32 = 1;
    pub const BT_ALTERN_BOND: u32 = 1;
    pub const BT_OTHER_ALTERN_BOND: u32 = 2;
    pub const BT_ALTERN_NS_BOND: u32 = 4;
    pub const BT_TAUTOM_BOND: u32 = 8;
    pub const BT_ALTERN_UNKN_BOND: u32 = 16;
    pub const BT_IGNORE_BOND: u32 = 0;
    pub const BT_NONSTEREO_MASK: u32 = 12;
    pub const BT_ALT_BOND_MASK: u32 = 3;
    pub const BT_NONTAUT_BOND_MASK: u32 = 7;
    pub const RI_ERR_ALLOC: i32 = -1;
    pub const RI_ERR_SYNTAX: i32 = -2;
    pub const RI_ERR_PROGR: i32 = -3;
    pub const tagAtTypeTotals_ATTOT_NUM_NP_Plus: tagAtTypeTotals = 0;
    pub const tagAtTypeTotals_ATTOT_NUM_NP_Proton: tagAtTypeTotals = 1;
    pub const tagAtTypeTotals_ATTOT_NUM_NP_H: tagAtTypeTotals = 2;
    pub const tagAtTypeTotals_ATTOT_NUM_N_Minus: tagAtTypeTotals = 3;
    pub const tagAtTypeTotals_ATTOT_NUM_NP: tagAtTypeTotals = 4;
    pub const tagAtTypeTotals_ATTOT_NUM_ON: tagAtTypeTotals = 5;
    pub const tagAtTypeTotals_ATTOT_NUM_COH: tagAtTypeTotals = 6;
    pub const tagAtTypeTotals_ATTOT_NUM_CSH: tagAtTypeTotals = 7;
    pub const tagAtTypeTotals_ATTOT_NUM_ZOH: tagAtTypeTotals = 8;
    pub const tagAtTypeTotals_ATTOT_NUM_OOH: tagAtTypeTotals = 9;
    pub const tagAtTypeTotals_ATTOT_NUM_ZOOH: tagAtTypeTotals = 10;
    pub const tagAtTypeTotals_ATTOT_NUM_NOH: tagAtTypeTotals = 11;
    pub const tagAtTypeTotals_ATTOT_NUM_N_OH: tagAtTypeTotals = 12;
    pub const tagAtTypeTotals_ATTOT_NUM_CO: tagAtTypeTotals = 13;
    pub const tagAtTypeTotals_ATTOT_NUM_ZO: tagAtTypeTotals = 14;
    pub const tagAtTypeTotals_ATTOT_NUM_NO: tagAtTypeTotals = 15;
    pub const tagAtTypeTotals_ATTOT_NUM_N_O: tagAtTypeTotals = 16;
    pub const tagAtTypeTotals_ATTOT_NUM_CO_Minus: tagAtTypeTotals = 17;
    pub const tagAtTypeTotals_ATTOT_NUM_CS_Minus: tagAtTypeTotals = 18;
    pub const tagAtTypeTotals_ATTOT_NUM_ZO_Minus: tagAtTypeTotals = 19;
    pub const tagAtTypeTotals_ATTOT_NUM_OO_Minus: tagAtTypeTotals = 20;
    pub const tagAtTypeTotals_ATTOT_NUM_ZOO_Minus: tagAtTypeTotals = 21;
    pub const tagAtTypeTotals_ATTOT_NUM_NO_Minus: tagAtTypeTotals = 22;
    pub const tagAtTypeTotals_ATTOT_NUM_N_O_Minus: tagAtTypeTotals = 23;
    pub const tagAtTypeTotals_ATTOT_NUM_O_Minus: tagAtTypeTotals = 24;
    pub const tagAtTypeTotals_ATTOT_NUM_OH_Plus: tagAtTypeTotals = 25;
    pub const tagAtTypeTotals_ATTOT_NUM_O_Plus: tagAtTypeTotals = 26;
    pub const tagAtTypeTotals_ATTOT_NUM_Proton: tagAtTypeTotals = 27;
    pub const tagAtTypeTotals_ATTOT_NUM_HalAnion: tagAtTypeTotals = 28;
    pub const tagAtTypeTotals_ATTOT_NUM_HalAcid: tagAtTypeTotals = 29;
    pub const tagAtTypeTotals_ATTOT_NUM_Errors: tagAtTypeTotals = 30;
    pub const tagAtTypeTotals_ATTOT_TOT_CHARGE: tagAtTypeTotals = 31;
    pub const tagAtTypeTotals_ATTOT_NUM_CHARGES: tagAtTypeTotals = 32;
    pub const tagAtTypeTotals_ATTOT_ARRAY_LEN: tagAtTypeTotals = 33;
    pub type tagAtTypeTotals = u32;
    pub use self::tagAtTypeTotals as AT_TYPE_TOTALS;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagProtonRemovalMaskAndType {
        pub typePos: i32,
        pub maskPos: i32,
        pub typeNeg: i32,
        pub maskNeg: i32,
        pub typeH: i32,
        pub maskH: i32,
    }
    impl ::std::default::Default for tagProtonRemovalMaskAndType {
        fn default() -> Self {
            Self {
                typePos: ::std::default::Default::default(),
                maskPos: ::std::default::Default::default(),
                typeNeg: ::std::default::Default::default(),
                maskNeg: ::std::default::Default::default(),
                typeH: ::std::default::Default::default(),
                maskH: ::std::default::Default::default(),
            }
        }
    }
    pub type PRMAT = tagProtonRemovalMaskAndType;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagAltPathChanges {
        pub nOldCapsVert: [[VertexFlow; 21usize]; 2usize],
        pub vOldVert: [Vertex; 2usize],
        pub bSetOldCapsVert: [S_CHAR; 2usize],
        pub vNewVertex: [Vertex; 2usize],
        pub bSetNew: [S_CHAR; 2usize],
    }
    impl ::std::default::Default for tagAltPathChanges {
        fn default() -> Self {
            Self {
                nOldCapsVert: ::std::array::from_fn(|_| {
                    ::std::array::from_fn(|_| ::std::default::Default::default())
                }),
                vOldVert: ::std::array::from_fn(|_| ::std::default::Default::default()),
                bSetOldCapsVert: ::std::array::from_fn(|_| ::std::default::Default::default()),
                vNewVertex: ::std::array::from_fn(|_| ::std::default::Default::default()),
                bSetNew: ::std::array::from_fn(|_| ::std::default::Default::default()),
            }
        }
    }
    pub type ALT_PATH_CHANGES = tagAltPathChanges;
}
pub(crate) mod local_ichi_io {
    use super::*;
    pub const INCHI_ADD_STR_LEN: u32 = 32768;
    pub const CONSOLE_LINE_LEN: u32 = 80;
    pub const FORCE_ANSI: u32 = 65536;
    pub const FORCE_UNICODE: u32 = 131072;
}
pub(crate) mod local_ichican2 {
    use super::*;
    pub const MAX_CELLS: u32 = 32766;
    pub const MAX_NODES: u32 = 32766;
    pub const MAX_SET_SIZE: u32 = 32766;
    pub const NORMALLY_ALLOWED_MAX_SET_SIZE: u32 = 2048;
    pub const MAX_LAYERS: u32 = 100;
    pub const SEPARATE_CANON_CALLS: u32 = 0;
    pub const INCHI_CANON_INFINITY: u32 = 32767;
    pub const EMPTY_CT: u32 = 0;
    pub const EMPTY_H_NUMBER: u32 = 32766;
    pub const BASE_H_NUMBER: u32 = 16383;
    pub const ELEM_NAME_LEN: u32 = 2;
    pub type Node = AT_NUMB;
    pub type Graph = NEIGH_LIST;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagUnorderedPartition {
        pub equ2: SourceMutPointer<AT_NUMB>,
    }
    impl ::std::default::Default for tagUnorderedPartition {
        fn default() -> Self {
            Self {
                equ2: ::std::default::Default::default(),
            }
        }
    }
    pub type UnorderedPartition = tagUnorderedPartition;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagCell {
        pub first: i32,
        pub next: i32,
        pub prev: i32,
    }
    impl ::std::default::Default for tagCell {
        fn default() -> Self {
            Self {
                first: ::std::default::Default::default(),
                next: ::std::default::Default::default(),
                prev: ::std::default::Default::default(),
            }
        }
    }
    pub type Cell = tagCell;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagTransposition {
        pub nAtNumb: SourceMutPointer<AT_NUMB>,
    }
    impl ::std::default::Default for tagTransposition {
        fn default() -> Self {
            Self {
                nAtNumb: ::std::default::Default::default(),
            }
        }
    }
    pub type Transposition = tagTransposition;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagCTable {
        pub Ctbl: SourceMutPointer<AT_RANK>,
        pub lenCt: i32,
        pub nLenCTAtOnly: i32,
        pub maxlenCt: i32,
        pub maxPos: i32,
        pub maxVert: i32,
        pub lenPos: i32,
        pub nextAtRank: SourceMutPointer<AT_RANK>,
        pub nextCtblPos: SourceMutPointer<AT_NUMB>,
        pub NumH: SourceMutPointer<NUM_H>,
        pub lenNumH: i32,
        pub maxlenNumH: i32,
        pub NumHfixed: SourceMutPointer<NUM_H>,
        pub iso_sort_key: SourceMutPointer<AT_ISO_SORT_KEY>,
        pub len_iso_sort_key: i32,
        pub maxlen_iso_sort_key: i32,
        pub iso_exchg_atnos: SourceMutPointer<S_CHAR>,
        pub len_iso_exchg_atnos: i32,
        pub maxlen_iso_exchg_atnos: i32,
    }
    impl ::std::default::Default for tagCTable {
        fn default() -> Self {
            Self {
                Ctbl: ::std::default::Default::default(),
                lenCt: ::std::default::Default::default(),
                nLenCTAtOnly: ::std::default::Default::default(),
                maxlenCt: ::std::default::Default::default(),
                maxPos: ::std::default::Default::default(),
                maxVert: ::std::default::Default::default(),
                lenPos: ::std::default::Default::default(),
                nextAtRank: ::std::default::Default::default(),
                nextCtblPos: ::std::default::Default::default(),
                NumH: ::std::default::Default::default(),
                lenNumH: ::std::default::Default::default(),
                maxlenNumH: ::std::default::Default::default(),
                NumHfixed: ::std::default::Default::default(),
                iso_sort_key: ::std::default::Default::default(),
                len_iso_sort_key: ::std::default::Default::default(),
                maxlen_iso_sort_key: ::std::default::Default::default(),
                iso_exchg_atnos: ::std::default::Default::default(),
                len_iso_exchg_atnos: ::std::default::Default::default(),
                maxlen_iso_exchg_atnos: ::std::default::Default::default(),
            }
        }
    }
    pub type ConTable = tagCTable;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagkLeast {
        pub k: i32,
        pub i: i32,
    }
    impl ::std::default::Default for tagkLeast {
        fn default() -> Self {
            Self {
                k: ::std::default::Default::default(),
                i: ::std::default::Default::default(),
            }
        }
    }
    pub type kLeast = tagkLeast;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagNodeValues {
        pub NumH: NUM_H,
        pub iso_sort_key: AT_ISO_SORT_KEY,
        pub NumHfixed: NUM_H,
        pub nAtNumber: AT_NUMB,
    }
    impl ::std::default::Default for tagNodeValues {
        fn default() -> Self {
            Self {
                NumH: ::std::default::Default::default(),
                iso_sort_key: ::std::default::Default::default(),
                NumHfixed: ::std::default::Default::default(),
                nAtNumber: ::std::default::Default::default(),
            }
        }
    }
    pub type NV = tagNodeValues;
}
pub(crate) mod local_ichicano {
    use super::*;
    pub const FullMaxClock: clock_t = -1;
    pub const HalfMaxClock: clock_t = 0;
}
pub(crate) mod local_ichicans {
    use super::*;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagStereoBondNeighbor {
        pub nRank: AT_RANK,
        pub nNeighRank1: AT_RANK,
        pub nNeighRank2: AT_RANK,
        pub num: AT_RANK,
        pub num_any_parity: AT_RANK,
        pub num_defined_parity: AT_RANK,
        pub what2do: AT_RANK,
        pub cumulene_len: U_CHAR,
        pub bond_type: U_CHAR,
    }
    impl ::std::default::Default for tagStereoBondNeighbor {
        fn default() -> Self {
            Self {
                nRank: ::std::default::Default::default(),
                nNeighRank1: ::std::default::Default::default(),
                nNeighRank2: ::std::default::Default::default(),
                num: ::std::default::Default::default(),
                num_any_parity: ::std::default::Default::default(),
                num_defined_parity: ::std::default::Default::default(),
                what2do: ::std::default::Default::default(),
                cumulene_len: ::std::default::Default::default(),
                bond_type: ::std::default::Default::default(),
            }
        }
    }
    pub type STEREO_BOND_NEIGH = tagStereoBondNeighbor;
}
pub(crate) mod local_ichimake {
    use super::*;
    pub const tagSp3StereoTypeTmp_SP3_NONE: tagSp3StereoTypeTmp = 0;
    pub const tagSp3StereoTypeTmp_SP3_ONLY: tagSp3StereoTypeTmp = 1;
    pub const tagSp3StereoTypeTmp_SP3_ABS: tagSp3StereoTypeTmp = 2;
    pub const tagSp3StereoTypeTmp_SP3_REL: tagSp3StereoTypeTmp = 4;
    pub const tagSp3StereoTypeTmp_SP3_RAC: tagSp3StereoTypeTmp = 8;
    pub const tagSp3StereoTypeTmp_SP3_TYPE: tagSp3StereoTypeTmp = 14;
    pub const tagSp3StereoTypeTmp_SP3_ANY: tagSp3StereoTypeTmp = 15;
    pub type tagSp3StereoTypeTmp = u32;
    pub use self::tagSp3StereoTypeTmp as SP3_TYPE_TMP;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagOrderStruct {
        pub m_gDfs4CT_nDfsNumber: SourceMutPointer<AT_NUMB>,
        pub m_gDfs4CT_nNumDescendants: SourceMutPointer<AT_NUMB>,
        pub m_gDfs4CT_nCurrentAtom: i32,
    }
    impl ::std::default::Default for tagOrderStruct {
        fn default() -> Self {
            Self {
                m_gDfs4CT_nDfsNumber: ::std::default::Default::default(),
                m_gDfs4CT_nNumDescendants: ::std::default::Default::default(),
                m_gDfs4CT_nCurrentAtom: ::std::default::Default::default(),
            }
        }
    }
    pub type OrderStruct = tagOrderStruct;
}
pub(crate) mod local_ichimap2 {
    use super::*;
    pub const MAP_MODE_STD: u32 = 0;
    pub const MAP_MODE_C2v: u32 = 1;
    pub const MAP_MODE_C2: u32 = 2;
    pub const MAP_MODE_S4: u32 = 3;
    pub const MAX_OTHER_NEIGH: u32 = 2;
    pub const NEIGH_MODE_RING: u32 = 1;
    pub const NEIGH_MODE_CHAIN: u32 = 2;
    pub const CHECKING_STEREOCENTER: u32 = 1;
    pub const CHECKING_STEREOBOND: u32 = 2;
    pub const COMP_STEREO_SUCCESS: u32 = 1;
    pub const NOT_WELL_DEF_UNKN: u32 = 2;
    pub const NOT_WELL_DEF_UNDF: u32 = 4;
    pub const PARITY_IMPOSSIBLE: u32 = 999;
}
pub(crate) mod local_ichimap4 {
    use super::*;
    pub const SB_DEPTH: u32 = 6;
}
pub(crate) mod local_ichinorm {
    use super::*;
    pub const DERIV_BRIDGE_O: u32 = 1;
    pub const DERIV_BRIDGE_NH: u32 = 2;
    pub const DERIV_AMINE_tN: u32 = 4;
    pub const DERIV_RING_O_OUTSIDE_PRECURSOR: u32 = 8;
    pub const DERIV_RING_NH_OUTSIDE_PRECURSOR: u32 = 16;
    pub const DERIV_X_OXIME: u32 = 32;
    pub const DERIV_UNMARK: u32 = 64;
    pub const DERIV_DUPLIC: u32 = 128;
    pub const DERIV_RO_COX: u32 = 256;
    pub const DERIV_RING_DMOX_DEOX_N: u32 = 512;
    pub const DERIV_RING_DMOX_DEOX_O: u32 = 1024;
    pub const DERIV_RING2_PRRLDD_OUTSIDE_PRECUR: u32 = 2048;
    pub const DERIV_RING2_PPRDN_OUTSIDE_PRECUR: u32 = 4096;
    pub const DERIV_DANSYL: u32 = 8192;
    pub const DERIV_RING_DMOX_DEOX: u32 = 1536;
    pub const DERIV_RING_OUTSIDE_PRECURSOR: u32 = 24;
    pub const DERIV_RING2_OUTSIDE_PRECUR: u32 = 6144;
    pub const DERIV_AT_LEN: u32 = 4;
    pub const DERIV_NOT: u32 = 4096;
    pub const MAX_AT_DERIV: u32 = 13;
    pub const NOT_AT_DERIV: u32 = 99;
    pub const MIN_AT_LEFT_DERIV: u32 = 2;
    pub const NO_ORD_VAL: u32 = 55;
    pub const CFLAG_MARK_BRANCH: u32 = 1;
    pub const CFLAG_MARK_BLOCK: u32 = 2;
    pub const COUNT_ALL_NOT_DERIV: u32 = 1;
    pub const OX_RING_SIZE: u32 = 5;
    pub const PRRLDD_RING_SIZE: u32 = 5;
    pub const PPRDN_RING_SIZE: u32 = 6;
    pub const MIN_PRRLDD_PPRDN_RING_SIZE: u32 = 5;
    pub const MAX_PRRLDD_PPRDN_RING_SIZE: u32 = 6;
    pub const DERIV_RING: u32 = 24;
    pub const UNDERIV_MAX_NUM: u32 = 512;
    pub const UNDERIV_LEN: u32 = 128;
    pub const UNDERIV_LEN2: u32 = 128;
    pub const UNDERIV_LIST_LEN: u32 = 2048;
    pub const UNDERIV_LIST_LEN2: u32 = 2048;
    pub const R2C_EMPTY: u32 = 127;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagAtPair {
        pub at: [AT_NUMB; 2usize],
        pub atno: AT_NUMB,
    }
    impl ::std::default::Default for tagAtPair {
        fn default() -> Self {
            Self {
                at: ::std::array::from_fn(|_| ::std::default::Default::default()),
                atno: ::std::default::Default::default(),
            }
        }
    }
    pub type R2C_ATPAIR = tagAtPair;
    pub const DERIV_UNEXPADABLE: i32 = 15904;
    pub const DERIV_REPL_N_WITH_O: i32 = 544;
    pub const DERIV_REPL_N_WITH_OH: i32 = 6144;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagDerivAttachment {
        pub typ: [i16; 4usize],
        pub ord: [i8; 4usize],
        pub num: [i8; 4usize],
        pub other_atom: AT_NUMB,
    }
    impl ::std::default::Default for tagDerivAttachment {
        fn default() -> Self {
            Self {
                typ: ::std::array::from_fn(|_| ::std::default::Default::default()),
                ord: ::std::array::from_fn(|_| ::std::default::Default::default()),
                num: ::std::array::from_fn(|_| ::std::default::Default::default()),
                other_atom: ::std::default::Default::default(),
            }
        }
    }
    pub type DERIV_AT = tagDerivAttachment;
    pub const tagDerivId_DERIV_ID_Acentonate: tagDerivId = 0;
    pub const tagDerivId_DERIV_ID_Benzlidene: tagDerivId = 1;
    pub const tagDerivId_DERIV_ID_BenzOX: tagDerivId = 2;
    pub const tagDerivId_DERIV_ID_BuBorate: tagDerivId = 3;
    pub const tagDerivId_DERIV_ID_Dansyl: tagDerivId = 4;
    pub const tagDerivId_DERIV_ID_DEOX: tagDerivId = 5;
    pub const tagDerivId_DERIV_ID_DMOX: tagDerivId = 6;
    pub const tagDerivId_DERIV_ID_EtBorate: tagDerivId = 7;
    pub const tagDerivId_DERIV_ID_EtOX: tagDerivId = 8;
    pub const tagDerivId_DERIV_ID_HFB: tagDerivId = 9;
    pub const tagDerivId_DERIV_ID_MeBorate: tagDerivId = 10;
    pub const tagDerivId_DERIV_ID_MOX: tagDerivId = 11;
    pub const tagDerivId_DERIV_ID_PFB: tagDerivId = 12;
    pub const tagDerivId_DERIV_ID_PFP: tagDerivId = 13;
    pub const tagDerivId_DERIV_ID_Piperidine: tagDerivId = 14;
    pub const tagDerivId_DERIV_ID_Pyrrolidide: tagDerivId = 15;
    pub const tagDerivId_DERIV_ID_TBDMS: tagDerivId = 16;
    pub const tagDerivId_DERIV_ID_TFA: tagDerivId = 17;
    pub const tagDerivId_DERIV_ID_TMS: tagDerivId = 18;
    pub const tagDerivId_DERIV_ID_Unknown: tagDerivId = 19;
    pub const tagDerivId_DERIV_ID_Acetate: tagDerivId = 20;
    pub const tagDerivId_DERIV_ID_Benzoate: tagDerivId = 21;
    pub type tagDerivId = u32;
    pub use self::tagDerivId as DerivId;
    pub const tagDerivBit_DERIV_BIT_Acentonate: tagDerivBit = 1;
    pub const tagDerivBit_DERIV_BIT_Benzlidene: tagDerivBit = 2;
    pub const tagDerivBit_DERIV_BIT_BenzOX: tagDerivBit = 4;
    pub const tagDerivBit_DERIV_BIT_BuBorate: tagDerivBit = 8;
    pub const tagDerivBit_DERIV_BIT_Dansyl: tagDerivBit = 16;
    pub const tagDerivBit_DERIV_BIT_DEOX: tagDerivBit = 32;
    pub const tagDerivBit_DERIV_BIT_DMOX: tagDerivBit = 64;
    pub const tagDerivBit_DERIV_BIT_EtBorate: tagDerivBit = 128;
    pub const tagDerivBit_DERIV_BIT_EtOX: tagDerivBit = 256;
    pub const tagDerivBit_DERIV_BIT_HFB: tagDerivBit = 512;
    pub const tagDerivBit_DERIV_BIT_MeBorate: tagDerivBit = 1024;
    pub const tagDerivBit_DERIV_BIT_MOX: tagDerivBit = 2048;
    pub const tagDerivBit_DERIV_BIT_PFB: tagDerivBit = 4096;
    pub const tagDerivBit_DERIV_BIT_PFP: tagDerivBit = 8192;
    pub const tagDerivBit_DERIV_BIT_Piperidine: tagDerivBit = 16384;
    pub const tagDerivBit_DERIV_BIT_Pyrrolidide: tagDerivBit = 32768;
    pub const tagDerivBit_DERIV_BIT_TBDMS: tagDerivBit = 65536;
    pub const tagDerivBit_DERIV_BIT_TFA: tagDerivBit = 131072;
    pub const tagDerivBit_DERIV_BIT_TMS: tagDerivBit = 262144;
    pub const tagDerivBit_DERIV_BIT_Unknown: tagDerivBit = 524288;
    pub const tagDerivBit_DERIV_BIT_Acetate: tagDerivBit = 1048576;
    pub const tagDerivBit_DERIV_BIT_Benzoate: tagDerivBit = 2097152;
    pub const tagDerivBit_DERIV_ID_NUMBER: tagDerivBit = 2097153;
    pub type tagDerivBit = u32;
    pub use self::tagDerivBit as DerivBit;
    pub type BIT_UNDERIV = i32;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagRing2Chain {
        pub type_: i8,
        pub ordW: i8,
        pub ordY: i8,
        pub ordC: i8,
        pub ordCopt: i8,
    }
    impl ::std::default::Default for tagRing2Chain {
        fn default() -> Self {
            Self {
                type_: ::std::default::Default::default(),
                ordW: ::std::default::Default::default(),
                ordY: ::std::default::Default::default(),
                ordC: ::std::default::Default::default(),
                ordCopt: ::std::default::Default::default(),
            }
        }
    }
    pub type R2C_AT = tagRing2Chain;
}
pub(crate) mod local_ichiparm {
    use super::*;
    pub const _POSIX_C_SOURCE: u32 = 200809;
    pub const VER103_DEFAULT_MODE: u32 = 6166;
}
pub(crate) mod local_ichiprt1 {
    use super::*;
    pub const NOT_YET_I2I_FOR_POLYMERS: u32 = 40;
    pub const OUT_NONTAUT: u32 = 4;
    pub const ORIG_STR_BUFLEN: u32 = 142;
    pub const sCompDelim: &[u8; 2] = b";\0";
    pub const sIdenticalValues: &[u8; 2] = b"*\0";
    pub const x_space: &[u8; 19] = b"                  \0";
    pub const x_inchi: &[u8; 6] = b"InChI\0";
    pub const x_inchi_ver: &[u8; 8] = b"version\0";
    pub const x_curr_ver: &[u8; 2] = b"1\0";
    pub const x_structure: &[u8; 10] = b"structure\0";
    pub const x_number: &[u8; 7] = b"number\0";
    pub const x_header: &[u8; 8] = b"id.name\0";
    pub const x_value: &[u8; 9] = b"id.value\0";
    pub const x_empty: &[u8; 1] = b"\0";
    pub const x_type: &[u8; 5] = b"type\0";
    pub const x_message: &[u8; 8] = b"message\0";
    pub const x_text: &[u8; 6] = b"value\0";
    pub const x_ferr: &[u8; 16] = b"fatal (aborted)\0";
    pub const x_err: &[u8; 17] = b"error (no InChI)\0";
    pub const x_warn: &[u8; 8] = b"warning\0";
    pub const x_basic: &[u8; 11] = b"identifier\0";
    pub const x_tautomeric: &[u8; 9] = b"mobile-H\0";
    pub const x_reconnected: &[u8; 12] = b"reconnected\0";
    pub const x_ver: &[u8; 8] = b"version\0";
    pub const x_type_alpha: &[u8; 6] = b"alpha\0";
    pub const x_type_numer: &[u8; 8] = b"numeric\0";
    pub const x_type_predec: &[u8; 4] = b"sct\0";
    pub const x_type_normal: &[u8; 7] = b"normal\0";
    pub const x_type_short: &[u8; 11] = b"compressed\0";
    pub const x_basic_layer: &[u8; 6] = b"basic\0";
    pub const x_aux_basic: &[u8; 26] = b"identifier.auxiliary-info\0";
    pub const x_aux_comm: &[u8; 70] =
        b"!-- This section is NOT a part of the identifier, it is not unique --\0";
    pub const x_ign_uu_sp2: &[u8; 17] = b"omit_undef_dbond\0";
    pub const x_ign_uu_sp3: &[u8; 15] = b"omit_undef_sp3\0";
    pub const x_line_opening: &[u8; 2] = b"<\0";
    pub const x_line_closing: &[u8; 3] = b"</\0";
    pub const x_close_line: &[u8; 2] = b">\0";
    pub const x_abs: &[u8; 2] = b"1\0";
    pub const x_rel: &[u8; 2] = b"2\0";
    pub const x_rac: &[u8; 2] = b"3\0";
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagInchiTag {
        pub szPlainLabel: SourceConstPointer<i8>,
        pub szPlainComment: SourceConstPointer<i8>,
        pub szXmlLabel: SourceConstPointer<i8>,
        pub bAlwaysOutput: i32,
    }
    impl ::std::default::Default for tagInchiTag {
        fn default() -> Self {
            Self {
                szPlainLabel: ::std::default::Default::default(),
                szPlainComment: ::std::default::Default::default(),
                szXmlLabel: ::std::default::Default::default(),
                bAlwaysOutput: ::std::default::Default::default(),
            }
        }
    }
    pub type INCHI_TAG = tagInchiTag;
    pub const tagIdentLblOrd_IL_FIXH_ORD: tagIdentLblOrd = 0;
    pub const tagIdentLblOrd_IL_ISOT_ORD: tagIdentLblOrd = 1;
    pub const tagIdentLblOrd_IL_STER_ORD: tagIdentLblOrd = 2;
    pub const tagIdentLblOrd_IL_VERS_ORD: tagIdentLblOrd = 3;
    pub const tagIdentLblOrd_IL_FML__ORD: tagIdentLblOrd = 4;
    pub const tagIdentLblOrd_IL_CONN_ORD: tagIdentLblOrd = 5;
    pub const tagIdentLblOrd_IL_ALLH_ORD: tagIdentLblOrd = 6;
    pub const tagIdentLblOrd_IL_CHRG_ORD: tagIdentLblOrd = 7;
    pub const tagIdentLblOrd_IL_PROT_ORD: tagIdentLblOrd = 8;
    pub const tagIdentLblOrd_IL_DBND_ORD: tagIdentLblOrd = 9;
    pub const tagIdentLblOrd_IL_SP3S_ORD: tagIdentLblOrd = 10;
    pub const tagIdentLblOrd_IL_INVS_ORD: tagIdentLblOrd = 11;
    pub const tagIdentLblOrd_IL_TYPS_ORD: tagIdentLblOrd = 12;
    pub const tagIdentLblOrd_IL_ATMS_ORD: tagIdentLblOrd = 13;
    pub const tagIdentLblOrd_IL_XCGA_ORD: tagIdentLblOrd = 14;
    pub const tagIdentLblOrd_IL_FMLF_ORD: tagIdentLblOrd = 15;
    pub const tagIdentLblOrd_IL_HFIX_ORD: tagIdentLblOrd = 16;
    pub const tagIdentLblOrd_IL_TRNS_ORD: tagIdentLblOrd = 17;
    pub const tagIdentLblOrd_IL_REC__ORD: tagIdentLblOrd = 18;
    pub const tagIdentLblOrd_IL_MAX_ORD: tagIdentLblOrd = 19;
    pub type tagIdentLblOrd = u32;
    pub use self::tagIdentLblOrd as IDENT_LBL_ORD;
    pub const tagIdentLblBit_IL_FIXH: tagIdentLblBit = 1;
    pub const tagIdentLblBit_IL_ISOT: tagIdentLblBit = 2;
    pub const tagIdentLblBit_IL_STER: tagIdentLblBit = 4;
    pub const tagIdentLblBit_IL_VERS: tagIdentLblBit = 8;
    pub const tagIdentLblBit_IL_FML_: tagIdentLblBit = 16;
    pub const tagIdentLblBit_IL_CONN: tagIdentLblBit = 32;
    pub const tagIdentLblBit_IL_ALLH: tagIdentLblBit = 64;
    pub const tagIdentLblBit_IL_CHRG: tagIdentLblBit = 128;
    pub const tagIdentLblBit_IL_PROT: tagIdentLblBit = 256;
    pub const tagIdentLblBit_IL_DBND: tagIdentLblBit = 512;
    pub const tagIdentLblBit_IL_SP3S: tagIdentLblBit = 1024;
    pub const tagIdentLblBit_IL_INVS: tagIdentLblBit = 2048;
    pub const tagIdentLblBit_IL_TYPS: tagIdentLblBit = 4096;
    pub const tagIdentLblBit_IL_ATMS: tagIdentLblBit = 8192;
    pub const tagIdentLblBit_IL_XCGA: tagIdentLblBit = 16384;
    pub const tagIdentLblBit_IL_FMLF: tagIdentLblBit = 32768;
    pub const tagIdentLblBit_IL_HFIX: tagIdentLblBit = 65536;
    pub const tagIdentLblBit_IL_TRNS: tagIdentLblBit = 131072;
    pub const tagIdentLblBit_IL_REC_: tagIdentLblBit = 262144;
    pub type tagIdentLblBit = u32;
    pub use self::tagIdentLblBit as IDENT_LBL_BIT;
    pub const tagAuxLblOrd_AL_FIXH_ORD: tagAuxLblOrd = 0;
    pub const tagAuxLblOrd_AL_ISOT_ORD: tagAuxLblOrd = 1;
    pub const tagAuxLblOrd_AL_STER_ORD: tagAuxLblOrd = 2;
    pub const tagAuxLblOrd_AL_REVR_ORD: tagAuxLblOrd = 3;
    pub const tagAuxLblOrd_AL_VERS_ORD: tagAuxLblOrd = 4;
    pub const tagAuxLblOrd_AL_NORM_ORD: tagAuxLblOrd = 5;
    pub const tagAuxLblOrd_AL_ANBR_ORD: tagAuxLblOrd = 6;
    pub const tagAuxLblOrd_AL_AEQU_ORD: tagAuxLblOrd = 7;
    pub const tagAuxLblOrd_AL_GEQU_ORD: tagAuxLblOrd = 8;
    pub const tagAuxLblOrd_AL_SP3I_ORD: tagAuxLblOrd = 9;
    pub const tagAuxLblOrd_AL_SP3N_ORD: tagAuxLblOrd = 10;
    pub const tagAuxLblOrd_AL_CRV__ORD: tagAuxLblOrd = 11;
    pub const tagAuxLblOrd_AL_ATMR_ORD: tagAuxLblOrd = 12;
    pub const tagAuxLblOrd_AL_BNDR_ORD: tagAuxLblOrd = 13;
    pub const tagAuxLblOrd_AL_XYZR_ORD: tagAuxLblOrd = 14;
    pub const tagAuxLblOrd_AL_FIXN_ORD: tagAuxLblOrd = 15;
    pub const tagAuxLblOrd_AL_ISON_ORD: tagAuxLblOrd = 16;
    pub const tagAuxLblOrd_AL_REC__ORD: tagAuxLblOrd = 17;
    pub const tagAuxLblOrd_AL_MAX_ORD: tagAuxLblOrd = 18;
    pub type tagAuxLblOrd = u32;
    pub use self::tagAuxLblOrd as AUX_LBL_ORD;
    pub const tagAuxLblBit_AL_FIXH: tagAuxLblBit = 1;
    pub const tagAuxLblBit_AL_ISOT: tagAuxLblBit = 2;
    pub const tagAuxLblBit_AL_STER: tagAuxLblBit = 4;
    pub const tagAuxLblBit_AL_REVR: tagAuxLblBit = 8;
    pub const tagAuxLblBit_AL_VERS: tagAuxLblBit = 16;
    pub const tagAuxLblBit_AL_NORM: tagAuxLblBit = 32;
    pub const tagAuxLblBit_AL_ANBR: tagAuxLblBit = 64;
    pub const tagAuxLblBit_AL_AEQU: tagAuxLblBit = 128;
    pub const tagAuxLblBit_AL_GEQU: tagAuxLblBit = 256;
    pub const tagAuxLblBit_AL_SP3I: tagAuxLblBit = 512;
    pub const tagAuxLblBit_AL_SP3N: tagAuxLblBit = 1024;
    pub const tagAuxLblBit_AL_CRV_: tagAuxLblBit = 2048;
    pub const tagAuxLblBit_AL_ATMR: tagAuxLblBit = 4096;
    pub const tagAuxLblBit_AL_BNDR: tagAuxLblBit = 8192;
    pub const tagAuxLblBit_AL_XYZR: tagAuxLblBit = 16384;
    pub const tagAuxLblBit_AL_FIXN: tagAuxLblBit = 32768;
    pub const tagAuxLblBit_AL_ISON: tagAuxLblBit = 65536;
    pub const tagAuxLblBit_AL_REC_: tagAuxLblBit = 131072;
    pub type tagAuxLblBit = u32;
    pub use self::tagAuxLblBit as AUX_LBL_BIT;
}
pub(crate) mod local_ichiprt2 {
    use super::*;
    pub const INIT_MIN_NUM_H: i32 = -4;
    pub const INIT_MAX_NUM_H: u32 = 16;
    pub const INIT_LEN_NUM_H: u32 = 21;
    pub const ALPHA_MINUS: u8 = 45u8;
    pub const ALPHA_ZERO_VAL: u8 = 46u8;
    pub const ALPHA_ONE: u8 = 97u8;
    pub const ALPHA_ZERO: u8 = 64u8;
    pub const DECIMAL_BASE: u32 = 10;
    pub const DECIMAL_MINUS: u8 = 45u8;
    pub const DECIMAL_ZERO_VAL: u8 = 48u8;
    pub const DECIMAL_ONE: u8 = 49u8;
    pub const DECIMAL_ZERO: u8 = 48u8;
}
pub(crate) mod local_ichiqueu {
    use super::*;
    pub const BOND_WRONG: u32 = 64;
    pub const MAX_DFS_DEPTH: u32 = 16;
    pub const PATH_LEN: u32 = 8;
    pub const NONE: i32 = 65535;
    pub type CHECK_DFS_RING = ::std::option::Option<
        fn(
            pCG: SourceMutPointer<tagCANON_GLOBALS>,
            atom: SourceMutPointer<inp_ATOM>,
            DfsPath: SourceMutPointer<DFS_PATH>,
            nLenDfsPath: i32,
            nStartAtomNeighbor: i32,
            nStartAtomNeighbor2: i32,
            nStartAtomNeighborNeighbor: i32,
            EndPoint: SourceMutPointer<T_ENDPOINT>,
            nMaxNumEndPoint: i32,
            BondPos: SourceMutPointer<T_BONDPOS>,
            nMaxNumBondPos: i32,
            pnNumEndPoint: SourceMutPointer<i32>,
            pnNumBondPos: SourceMutPointer<i32>,
            pBNS: SourceMutPointer<BalancedNetworkStructure>,
            pBD: SourceMutPointer<BalancedNetworkData>,
            num_atoms: i32,
        ) -> i32,
    >;
    pub type CHECK_CENTERPOINT =
        ::std::option::Option<fn(atom: SourceMutPointer<inp_ATOM>, iat: i32) -> i32>;
    pub type CHECK_DFS_PATH = ::std::option::Option<
        fn(
            pCG: SourceMutPointer<tagCANON_GLOBALS>,
            atom: SourceMutPointer<inp_ATOM>,
            DfsPath: SourceMutPointer<DFS_PATH>,
            nLenDfsPath: i32,
            jNxtNeigh: i32,
            nStartAtomNeighbor: i32,
            nStartAtomNeighbor2: i32,
            nStartAtomNeighborNeighbor: i32,
            EndPoint: SourceMutPointer<T_ENDPOINT>,
            nMaxNumEndPoint: i32,
            BondPos: SourceMutPointer<T_BONDPOS>,
            nMaxNumBondPos: i32,
            pnNumEndPoint: SourceMutPointer<i32>,
            pnNumBondPos: SourceMutPointer<i32>,
            pBNS: SourceMutPointer<BalancedNetworkStructure>,
            pBD: SourceMutPointer<BalancedNetworkData>,
            num_atoms: i32,
        ) -> i32,
    >;
    pub type CHECK_DFS_CENTERPOINT = ::std::option::Option<
        fn(
            atom: SourceMutPointer<inp_ATOM>,
            DfsPath: SourceMutPointer<DFS_PATH>,
            nLenDfsPath: i32,
            jNxtNeigh: i32,
            pBNS: SourceMutPointer<BalancedNetworkStructure>,
            pBD: SourceMutPointer<BalancedNetworkData>,
            num_atoms: i32,
        ) -> i32,
    >;
}
pub(crate) mod local_ichiread {
    use super::*;
    pub const SEGM_LINE_ADD: u32 = 128;
    pub const LINKED_BOND_ADD: u32 = 128;
    pub const SEG_END: u8 = 47u8;
    pub const INCHI_TOKEN: &[u8; 6] = b"/\n\r\t\\\0";
    pub const IST_HAPPENED_IN_RECMET: u32 = 100;
    pub const NSTRLEN: u32 = 524288;
    pub const MAX_MSG_LEN: u32 = 512;
    pub const MAX_MSG_BUF_LEN: u32 = 128;
    pub const LAST_AT_LEN: u32 = 256;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagOneLinkedBond {
        pub neigh: AT_NUMB,
        pub prev: AT_NUMB,
    }
    impl ::std::default::Default for tagOneLinkedBond {
        fn default() -> Self {
            Self {
                neigh: ::std::default::Default::default(),
                prev: ::std::default::Default::default(),
            }
        }
    }
    pub type ONE_LINKED_BOND = tagOneLinkedBond;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagLinkedBonds {
        pub pBond: SourceMutPointer<ONE_LINKED_BOND>,
        pub len: i32,
        pub len_alloc: i32,
    }
    impl ::std::default::Default for tagLinkedBonds {
        fn default() -> Self {
            Self {
                pBond: ::std::default::Default::default(),
                len: ::std::default::Default::default(),
                len_alloc: ::std::default::Default::default(),
            }
        }
    }
    pub type LINKED_BONDS = tagLinkedBonds;
    pub const tagModeProtonIsoExchgH_MODE_PIXH_UNDEFINED: tagModeProtonIsoExchgH = 0;
    pub const tagModeProtonIsoExchgH_MODE_PIXH_ADD_TO_FIRST: tagModeProtonIsoExchgH = 1;
    pub const tagModeProtonIsoExchgH_MODE_PIXH_ADD_TO_EACH: tagModeProtonIsoExchgH = 2;
    pub const tagModeProtonIsoExchgH_MODE_PIXH_ADD_A_PIXH_COMPONENT: tagModeProtonIsoExchgH = 3;
    pub const tagModeProtonIsoExchgH_MODE_PIXH_KEEP_TOTALS: tagModeProtonIsoExchgH = 4;
    pub type tagModeProtonIsoExchgH = u32;
    pub use self::tagModeProtonIsoExchgH as MODE_PIXH;
    pub const tagInChI_STATE_IST_MOBILE_H_FORMULA: tagInChI_STATE = 0;
    pub const tagInChI_STATE_IST_MOBILE_H_CONNECTIONS: tagInChI_STATE = 1;
    pub const tagInChI_STATE_IST_MOBILE_H: tagInChI_STATE = 2;
    pub const tagInChI_STATE_IST_MOBILE_H_CHARGE: tagInChI_STATE = 3;
    pub const tagInChI_STATE_IST_MOBILE_H_PROTONS: tagInChI_STATE = 4;
    pub const tagInChI_STATE_IST_MOBILE_H_SP2: tagInChI_STATE = 5;
    pub const tagInChI_STATE_IST_MOBILE_H_SP3: tagInChI_STATE = 6;
    pub const tagInChI_STATE_IST_MOBILE_H_SP3_M: tagInChI_STATE = 7;
    pub const tagInChI_STATE_IST_MOBILE_H_SP3_S: tagInChI_STATE = 8;
    pub const tagInChI_STATE_IST_MOBILE_H_ISO_LAYER_FORK: tagInChI_STATE = 9;
    pub const tagInChI_STATE_IST_MOBILE_H_ISO_ATOMS: tagInChI_STATE = 10;
    pub const tagInChI_STATE_IST_MOBILE_H_ISO_EXCH_H: tagInChI_STATE = 11;
    pub const tagInChI_STATE_IST_MOBILE_H_ISO_SP2: tagInChI_STATE = 12;
    pub const tagInChI_STATE_IST_MOBILE_H_ISO_SP3: tagInChI_STATE = 13;
    pub const tagInChI_STATE_IST_MOBILE_H_ISO_SP3_M: tagInChI_STATE = 14;
    pub const tagInChI_STATE_IST_MOBILE_H_ISO_SP3_S: tagInChI_STATE = 15;
    pub const tagInChI_STATE_IST_FIXED_H_LAYER_FORK: tagInChI_STATE = 16;
    pub const tagInChI_STATE_IST_FIXED_H_FORMULA: tagInChI_STATE = 17;
    pub const tagInChI_STATE_IST_FIXED_H: tagInChI_STATE = 18;
    pub const tagInChI_STATE_IST_FIXED_H_CHARGE: tagInChI_STATE = 19;
    pub const tagInChI_STATE_IST_FIXED_H_SP2: tagInChI_STATE = 20;
    pub const tagInChI_STATE_IST_FIXED_H_SP3: tagInChI_STATE = 21;
    pub const tagInChI_STATE_IST_FIXED_H_SP3_M: tagInChI_STATE = 22;
    pub const tagInChI_STATE_IST_FIXED_H_SP3_S: tagInChI_STATE = 23;
    pub const tagInChI_STATE_IST_FIXED_H_PERMUTATION: tagInChI_STATE = 24;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_LAYER_FORK: tagInChI_STATE = 25;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_ATOMS: tagInChI_STATE = 26;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_LAYER: tagInChI_STATE = 27;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_SP2: tagInChI_STATE = 28;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_SP3: tagInChI_STATE = 29;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_SP3_M: tagInChI_STATE = 30;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_SP3_S: tagInChI_STATE = 31;
    pub const tagInChI_STATE_IST_FIXED_H_ISO_PERMUTATION: tagInChI_STATE = 32;
    pub const tagInChI_STATE_IST_RECONNECTED_LAYER_FORK: tagInChI_STATE = 33;
    pub const tagInChI_STATE_IST_RECONNECTED_FORMULA: tagInChI_STATE = 34;
    pub const tagInChI_STATE_IST_MATERIAL_BALANCE_ERROR: tagInChI_STATE = 35;
    pub const tagInChI_STATE_IST_MOBILE_H_POLYMER: tagInChI_STATE = 36;
    pub const tagInChI_STATE_IST_END: tagInChI_STATE = -1;
    pub type tagInChI_STATE = i32;
    pub use self::tagInChI_STATE as INCHI_STATE;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagInchiReadErrMsg {
        pub stat: i32,
        pub msg: SourceConstPointer<i8>,
    }
    impl ::std::default::Default for tagInchiReadErrMsg {
        fn default() -> Self {
            Self {
                stat: ::std::default::Default::default(),
                msg: ::std::default::Default::default(),
            }
        }
    }
    pub type INCHI_READ_ERR_MSG = tagInchiReadErrMsg;
    pub const tagCopySegmentType_CPY_SP2: tagCopySegmentType = 0;
    pub const tagCopySegmentType_CPY_SP3: tagCopySegmentType = 1;
    pub const tagCopySegmentType_CPY_SP3_M: tagCopySegmentType = 2;
    pub const tagCopySegmentType_CPY_SP3_S: tagCopySegmentType = 3;
    pub const tagCopySegmentType_CPY_ISO_AT: tagCopySegmentType = 4;
    pub type tagCopySegmentType = u32;
    pub use self::tagCopySegmentType as COPY_SEG_TYPE;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagNumElem {
        pub num: i32,
    }
    impl ::std::default::Default for tagNumElem {
        fn default() -> Self {
            Self {
                num: ::std::default::Default::default(),
            }
        }
    }
    pub type NUM_ELEM = tagNumElem;
    pub const tagAuxInfoState_AST_VERSION: tagAuxInfoState = 0;
    pub const tagAuxInfoState_AST_MOBILE_H_NUMBERS: tagAuxInfoState = 1;
    pub const tagAuxInfoState_AST_MOBILE_H_ATOM_EQ: tagAuxInfoState = 2;
    pub const tagAuxInfoState_AST_MOBILE_H_GROUP_EQ: tagAuxInfoState = 3;
    pub const tagAuxInfoState_AST_MOBILE_H_SP3_INV: tagAuxInfoState = 4;
    pub const tagAuxInfoState_AST_MOBILE_H_SP3_INV_NUMBERS: tagAuxInfoState = 5;
    pub const tagAuxInfoState_AST_MOBILE_H_ISO_LAYER_FORK: tagAuxInfoState = 6;
    pub const tagAuxInfoState_AST_MOBILE_H_ISO_NUMBERS: tagAuxInfoState = 7;
    pub const tagAuxInfoState_AST_MOBILE_H_ISO_ATOM_EQ: tagAuxInfoState = 8;
    pub const tagAuxInfoState_AST_MOBILE_H_ISO_GROUP_EQ: tagAuxInfoState = 9;
    pub const tagAuxInfoState_AST_MOBILE_H_ISO_SP3_INV: tagAuxInfoState = 10;
    pub const tagAuxInfoState_AST_MOBILE_H_ISO_SP3_INV_NUMBERS: tagAuxInfoState = 11;
    pub const tagAuxInfoState_AST_FIXED_H_LAYER_FORK: tagAuxInfoState = 12;
    pub const tagAuxInfoState_AST_FIXED_H_NUMBERS: tagAuxInfoState = 13;
    pub const tagAuxInfoState_AST_FIXED_H_ATOM_EQ: tagAuxInfoState = 14;
    pub const tagAuxInfoState_AST_FIXED_H_SP3_INV: tagAuxInfoState = 15;
    pub const tagAuxInfoState_AST_FIXED_H_SP3_INV_NUMBERS: tagAuxInfoState = 16;
    pub const tagAuxInfoState_AST_FIXED_H_ISO_LAYER_FORK: tagAuxInfoState = 17;
    pub const tagAuxInfoState_AST_FIXED_H_ISO_NUMBERS: tagAuxInfoState = 18;
    pub const tagAuxInfoState_AST_FIXED_H_ISO_ATOM_EQ: tagAuxInfoState = 19;
    pub const tagAuxInfoState_AST_FIXED_H_ISO_SP3_INV: tagAuxInfoState = 20;
    pub const tagAuxInfoState_AST_FIXED_H_ISO_SP3_INV_NUMBERS: tagAuxInfoState = 21;
    pub const tagAuxInfoState_AST_REVERSE_INFO_CRV: tagAuxInfoState = 22;
    pub const tagAuxInfoState_AST_REVERSE_INFO_ATOMS: tagAuxInfoState = 23;
    pub const tagAuxInfoState_AST_REVERSE_INFO_BONDS: tagAuxInfoState = 24;
    pub const tagAuxInfoState_AST_REVERSE_INFO_XYZ: tagAuxInfoState = 25;
    pub const tagAuxInfoState_AST_RECONNECTED_LAYER_FORK: tagAuxInfoState = 26;
    pub const tagAuxInfoState_AST_RECONNECTED_LAYER_NUMBERS: tagAuxInfoState = 27;
    pub type tagAuxInfoState = u32;
    pub use self::tagAuxInfoState as AUX_INFO_STATE;
}
pub(crate) mod local_ichirvr1 {
    use super::*;
    pub const MIN_ATOM_CHARGE: i32 = -2;
    pub const MAX_ATOM_CHARGE: u32 = 2;
    pub const NEUTRAL_STATE: u32 = 2;
    pub const NUM_ATOM_CHARGES: u32 = 5;
    pub const MAX_NUM_VALENCES: u32 = 5;
    pub const NUM_VF: u32 = 3;
    pub const VF_USED_IN: u32 = 1;
    pub const VF_USED_OUT: u32 = 2;
    pub const VF_USED_ALL: u32 = 3;
    pub const cnListNumEl: i32 = 18;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagVertexFlow {
        pub type_: i32,
        pub v: Vertex,
        pub e_In: EdgeIndex,
        pub e_Out: EdgeIndex,
        pub delta_In: EdgeFlow,
        pub delta_Out: EdgeFlow,
        pub bUsed: Vertex,
    }
    impl ::std::default::Default for tagVertexFlow {
        fn default() -> Self {
            Self {
                type_: ::std::default::Default::default(),
                v: ::std::default::Default::default(),
                e_In: ::std::default::Default::default(),
                e_Out: ::std::default::Default::default(),
                delta_In: ::std::default::Default::default(),
                delta_Out: ::std::default::Default::default(),
                bUsed: ::std::default::Default::default(),
            }
        }
    }
    pub type VF = tagVertexFlow;
}
pub(crate) mod local_ichirvr2 {
    use super::*;
    pub const FIX_BOND_ADD_ALLOC: u32 = 128;
    pub const INC_EDGE_LIST: u32 = 16;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagMobileHGroups {
        pub group_number: AT_NUMB,
        pub atom_number: AT_NUMB,
        pub atom_type_pVA: AT_NUMB,
        pub ineigh: S_CHAR,
        pub bond_type: S_CHAR,
        pub forbidden: S_CHAR,
        pub endpoint_valence: S_CHAR,
        pub num_bonds: S_CHAR,
        pub bonds_valence: S_CHAR,
        pub num_bonds_non_metal: S_CHAR,
        pub bonds_valence_non_metal: S_CHAR,
    }
    impl ::std::default::Default for tagMobileHGroups {
        fn default() -> Self {
            Self {
                group_number: ::std::default::Default::default(),
                atom_number: ::std::default::Default::default(),
                atom_type_pVA: ::std::default::Default::default(),
                ineigh: ::std::default::Default::default(),
                bond_type: ::std::default::Default::default(),
                forbidden: ::std::default::Default::default(),
                endpoint_valence: ::std::default::Default::default(),
                num_bonds: ::std::default::Default::default(),
                bonds_valence: ::std::default::Default::default(),
                num_bonds_non_metal: ::std::default::Default::default(),
                bonds_valence_non_metal: ::std::default::Default::default(),
            }
        }
    }
    pub type MOBILE_GR = tagMobileHGroups;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagMobileGroupList {
        pub group_number: AT_NUMB,
        pub num: AT_NUMB,
    }
    impl ::std::default::Default for tagMobileGroupList {
        fn default() -> Self {
            Self {
                group_number: ::std::default::Default::default(),
                num: ::std::default::Default::default(),
            }
        }
    }
    pub type MGROUPS = tagMobileGroupList;
}
pub(crate) mod local_ichirvr3 {
    use super::*;
    pub const INC_ADD_EDGE: u32 = 64;
    pub const fNumRPosChgH: u32 = 0;
    pub const fNumRPosChgU: u32 = 1;
    pub const fNumRNegChgO: u32 = 2;
    pub const fNumRNegChgN: u32 = 3;
    pub const fNumRNeutrlH: u32 = 4;
    pub const fNumNPosChgH: u32 = 5;
    pub const fNumNPosChgU: u32 = 6;
    pub const fNumNNegChgO: u32 = 7;
    pub const fNumNNegChgN: u32 = 8;
    pub const fNumNNeutrlH: u32 = 9;
    pub const fNumAllChgT: u32 = 10;
    pub const CHG_SET_NOOH: u32 = 0;
    pub const CHG_SET_WRONG_TAUT: u32 = 1;
    pub const CHG_SET_TAUT: u32 = 2;
    pub const CHG_LAST_SET: u32 = 2;
    pub const CHG_SET_O_FIXED: u32 = 3;
    pub const CHG_SET_NUM: u32 = 4;
    pub const CHG_SET_MISSED_TAUT: u32 = 0;
    pub const CHG_SET_OTHER_TAUT_O: u32 = 1;
    pub const CHG_SET_OTHER_TAUT_N: u32 = 2;
    pub const CHG_SET_NN: u32 = 3;
    pub const CHG_SET_AVOID: u32 = 4;
    pub const CHG_SET_MISSED_TAUT_1: u32 = 0;
    pub const CHG_SET_MISSED_TAUT_ALL: u32 = 1;
    pub const CHG_SET_OTHER_TAUT_1: u32 = 2;
    pub const CHG_SET_OTHER_TAUT_ALL: u32 = 3;
    pub const CHG_SET_NO_IN_NO2M2: u32 = 4;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagTgDiffHChgFH {
        pub itg: i16,
        pub nNumHInchi: i16,
        pub nNumHRevrs: i16,
        pub nNumHNorml: i16,
        pub nNumMInchi: i16,
        pub nNumMRevrs: i16,
        pub nNumMNorml: i16,
        pub nNumPRevrs: i16,
        pub nNumPNorml: i16,
        pub n: [i16; 10usize],
        pub i: [i16; 10usize],
    }
    impl ::std::default::Default for tagTgDiffHChgFH {
        fn default() -> Self {
            Self {
                itg: ::std::default::Default::default(),
                nNumHInchi: ::std::default::Default::default(),
                nNumHRevrs: ::std::default::Default::default(),
                nNumHNorml: ::std::default::Default::default(),
                nNumMInchi: ::std::default::Default::default(),
                nNumMRevrs: ::std::default::Default::default(),
                nNumMNorml: ::std::default::Default::default(),
                nNumPRevrs: ::std::default::Default::default(),
                nNumPNorml: ::std::default::Default::default(),
                n: ::std::array::from_fn(|_| ::std::default::Default::default()),
                i: ::std::array::from_fn(|_| ::std::default::Default::default()),
            }
        }
    }
    pub type TgDiffHChgFH = tagTgDiffHChgFH;
}
pub(crate) mod local_ichirvr4 {
    use super::*;
    pub const MAX_NUM_CARBON_CHARGE_EDGES: u32 = 2;
}
pub(crate) mod local_ichirvr5 {
    use super::*;
    pub const INC_ADD_EDGE: u32 = 64;
    pub const CHG_SET_WRONG_TAUT_N: u32 = 0;
    pub const CHG_SET_WRONG_TAUT_O: u32 = 1;
    pub const CHG_SET_WRONG_TAUT_ALL: u32 = 2;
    pub const CHG_LAST_SET: u32 = 2;
    pub const CHG_SET_O_FIXED: u32 = 3;
    pub const CHG_SET_NUM: u32 = 4;
}
pub(crate) mod local_ichirvr6 {
    use super::*;
    pub const INC_ADD_EDGE: u32 = 64;
}
pub(crate) mod local_ichister {
    use super::*;
    pub const ZTYPE_DOWN: i32 = -1;
    pub const ZTYPE_NONE: u32 = 0;
    pub const ZTYPE_UP: u32 = 1;
    pub const ZTYPE_3D: u32 = 3;
    pub const ZTYPE_EITHER: u32 = 9999;
    pub const ARR_DIM: u32 = 3;
    pub const MIN_ANGLE: f64 = 0.1;
    pub const MIN_SINE: f64 = 0.03;
    pub const MIN_ANGLE_RELAXED: f64 = 0.001;
    pub const MIN_SINE_RELAXED: f64 = 0.001;
    pub const MIN_ANGLE_DBOND: f64 = 0.087156;
    pub const MIN_SINE_OUTSIDE: f64 = 0.06;
    pub const MIN_SINE_SQUARE: f64 = 0.125;
    pub const MIN_SINE_EDGE: f64 = 0.167;
    pub const MIN_LEN_STRAIGHT: f64 = 1.9;
    pub const MAX_SINE: f64 = 0.7071067811865476;
    pub const MIN_BOND_LEN: f64 = 0.000001;
    pub const ZERO_LENGTH: f64 = 0.000001;
    pub const BOND_PARITY_UNDEFINED: u32 = 64;
    pub const MPY_SINE: f64 = 1.0;
    pub const MAX_EDGE_RATIO: f64 = 2.5;
    pub const T2D_OKAY: u32 = 1;
    pub const T2D_WARN: u32 = 2;
    pub const T2D_UNDF: u32 = 4;
    pub const ZERO_ANGLE: f64 = 0.000001;
    pub const PHOSPHINE_STEREO: u32 = 19;
    pub const ARSINE_STEREO: u32 = 20;
    pub const AB_NEGATIVE: u32 = 16;
    pub const AB_UNKNOWN: u32 = 32;
    pub const ADD_EXPLICIT_HYDROGEN_NEIGH: u32 = 1;
    pub const ADD_EXPLICIT_LONE_PAIR_NEIGH: u32 = 2;
}
pub(crate) mod local_ichitaut {
    use super::*;
    pub const C_SUBTYPE_CHARGED: u32 = 0;
    pub const C_SUBTYPE_p_DONOR: u32 = 1;
    pub const C_SUBTYPE_p_ACCEPT: u32 = 2;
    pub const C_SUBTYPE_H_ACCEPT: u32 = 4;
    pub const C_SUBTYPE_H_DONOR: u32 = 8;
    pub const C_SUBTYPE_NEUTRAL: u32 = 16;
    pub const MAX_STACK_ARRAY_LEN: u32 = 127;
    pub const MAX_TGROUP_ARRAY_LEN: u32 = 127;
    pub const C_SUBTYPE_CHARGED_NON_TAUT: u32 = 0;
    pub const C_SUBTYPE_CHARGED_p_DONOR: u32 = 1;
    pub const C_SUBTYPE_CHARGED_H_ACCEPT: u32 = 4;
    pub const C_SUBTYPE_CHARGED_H_ACCEPT_p_DONOR: u32 = 5;
    pub const C_SUBTYPE_CHARGED_H_DONOR: u32 = 9;
    pub const C_SUBTYPE_NEUTRAL_NON_TAUT: u32 = 16;
    pub const C_SUBTYPE_NEUTRAL_H_ACCEPT: u32 = 20;
    pub const C_SUBTYPE_NEUTRAL_H_ACCEPT_p_ACCEPT: u32 = 22;
    pub const C_SUBTYPE_NEUTRAL_H_DONOR: u32 = 24;
    pub const ALT_PATH_FOUND: u32 = 32767;
    pub const NO_ENDPOINT: u32 = 32768;
    pub const DISABLE_CANDIDATE: u32 = 10;
    pub const ACCEPTOR_PAIR: u32 = 1;
    pub const DONOR_PAIR: u32 = 2;
    pub const MAX_LOCAL_TGNUM: u32 = 0;
    pub const MAX_ALT_PATH_LEN: u32 = 8;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagChargeType {
        pub elname: [i8; 3usize],
        pub charge: S_CHAR,
        pub neutral_valence: S_CHAR,
        pub neutral_bonds_valence: S_CHAR,
        pub cChangeValence: S_CHAR,
        pub cChargeType: S_CHAR,
        pub num_bonds: S_CHAR,
    }
    impl ::std::default::Default for tagChargeType {
        fn default() -> Self {
            Self {
                elname: ::std::array::from_fn(|_| ::std::default::Default::default()),
                charge: ::std::default::Default::default(),
                neutral_valence: ::std::default::Default::default(),
                neutral_bonds_valence: ::std::default::Default::default(),
                cChangeValence: ::std::default::Default::default(),
                cChargeType: ::std::default::Default::default(),
                num_bonds: ::std::default::Default::default(),
            }
        }
    }
    pub type CHARGE_TYPE = tagChargeType;
}
pub(crate) mod local_ikey_dll {
    use super::*;
    pub const INCHIKEY_DEBUG: u32 = 0;
    pub const INCHIKEY_FLAG_OK: u32 = 0;
    pub const INCHIKEY_NOT_VALID_FLAG: u32 = 1;
    pub const MINOUTLENGTH: _bindgen_ty_2 = 256;
    pub type _bindgen_ty_2 = u32;
}
pub(crate) mod local_mol_fmt2 {
    use super::*;
    pub const MIN_STDATA_X_COORD: f64 = 0.0;
    pub const MAX_STDATA_X_COORD: f64 = 256.0;
    pub const MIN_STDATA_Y_COORD: f64 = 0.0;
    pub const MAX_STDATA_Y_COORD: f64 = 256.0;
    pub const MIN_STDATA_Z_COORD: f64 = 0.0;
    pub const MAX_STDATA_Z_COORD: f64 = 256.0;
    pub const MAX_STDATA_AVE_BOND_LENGTH: f64 = 20.0;
    pub const MIN_STDATA_AVE_BOND_LENGTH: f64 = 10.0;
}
pub(crate) mod local_mol_fmt4 {
    use super::*;
    pub const sdf_data_hdr_name: &[u8; 5] = b"NAME\0";
    pub const sdf_data_hdr_comm: &[u8; 8] = b"COMMENT\0";
    pub const SDF_START: _bindgen_ty_3 = 0;
    pub const SDF_DATA_HEADER: _bindgen_ty_3 = 1;
    pub const SDF_DATA_HEADER_NAME: _bindgen_ty_3 = 2;
    pub const SDF_DATA_HEADER_COMMENT: _bindgen_ty_3 = 3;
    pub const SDF_DATA_HEADER_CAS: _bindgen_ty_3 = 4;
    pub const SDF_DATA_HEADER_USER: _bindgen_ty_3 = 5;
    pub const SDF_DATA_LINE: _bindgen_ty_3 = 6;
    pub const SD_FMT_END_OF_DATA_ITEM: _bindgen_ty_3 = 7;
    pub const SDF_EMPTY_LINE: _bindgen_ty_3 = 8;
    pub const SD_FMT_END_OF_DATA_BLOCK: _bindgen_ty_3 = 9;
    pub type _bindgen_ty_3 = u32;
}
pub(crate) mod local_readinch {
    use super::*;
    pub const ISOLATED_ATOM: u32 = 15;
    pub const MAX_CHAIN_LEN: u32 = 20;
}
pub(crate) mod local_runichi {
    use super::*;
    pub const gsMissing: &[u8; 11] = b"is missing\0";
    pub const gsEmpty: &[u8; 1] = b"\0";
    pub const gsSpace: &[u8; 2] = b" \0";
    pub const gsEqual: &[u8; 2] = b"=\0";
}
pub(crate) mod local_runichi2 {
    use super::*;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagModSCenterInfo {
        pub num: i32,
        pub valence: i32,
        pub n_stereo: i32,
        pub nbr: [i32; 20usize],
        pub new_nbr: [i32; 20usize],
    }
    impl ::std::default::Default for tagModSCenterInfo {
        fn default() -> Self {
            Self {
                num: ::std::default::Default::default(),
                valence: ::std::default::Default::default(),
                n_stereo: ::std::default::Default::default(),
                nbr: ::std::array::from_fn(|_| ::std::default::Default::default()),
                new_nbr: ::std::array::from_fn(|_| ::std::default::Default::default()),
            }
        }
    }
    pub type ModSCenterInfo = tagModSCenterInfo;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagDiylFrag {
        pub na: i32,
        pub nb: i32,
        pub end1: i32,
        pub end2: i32,
        pub alist: SourceMutPointer<i32>,
        pub xclist: SourceMutPointer<i32>,
        pub sig: INCHI_IOS_STRING,
    }
    impl ::std::default::Default for tagDiylFrag {
        fn default() -> Self {
            Self {
                na: ::std::default::Default::default(),
                nb: ::std::default::Default::default(),
                end1: ::std::default::Default::default(),
                end2: ::std::default::Default::default(),
                alist: ::std::default::Default::default(),
                xclist: ::std::default::Default::default(),
                sig: ::std::default::Default::default(),
            }
        }
    }
    pub type DiylFrag = tagDiylFrag;
}
pub(crate) mod local_sha2 {
    use super::*;
    pub const _CRT_SECURE_NO_DEPRECATE: u32 = 1;
}
pub(crate) mod local_strutil {
    use super::*;
    pub const FIRST_NEIGHB2: u32 = 4;
    pub const FIRST_CENTER2: u32 = 5;
    pub const MAX_NEIGH: u32 = 6;
    pub const NUM_SEGM: u32 = 20;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagTreeAtom {
        pub neighbor: [AT_NUMB; 20usize],
        pub valence: S_CHAR,
        pub nRingSystem: AT_NUMB,
        pub nBlockSystem: AT_NUMB,
        pub bCutVertex: S_CHAR,
    }
    impl ::std::default::Default for tagTreeAtom {
        fn default() -> Self {
            Self {
                neighbor: ::std::array::from_fn(|_| ::std::default::Default::default()),
                valence: ::std::default::Default::default(),
                nRingSystem: ::std::default::Default::default(),
                nBlockSystem: ::std::default::Default::default(),
                bCutVertex: ::std::default::Default::default(),
            }
        }
    }
    pub type tre_ATOM = tagTreeAtom;
    pub const tagIonAtomType_IAT_H: tagIonAtomType = 0;
    pub const tagIonAtomType_IAT_C: tagIonAtomType = 1;
    pub const tagIonAtomType_IAT_N: tagIonAtomType = 2;
    pub const tagIonAtomType_IAT_P: tagIonAtomType = 3;
    pub const tagIonAtomType_IAT_O: tagIonAtomType = 4;
    pub const tagIonAtomType_IAT_S: tagIonAtomType = 5;
    pub const tagIonAtomType_IAT_Se: tagIonAtomType = 6;
    pub const tagIonAtomType_IAT_Te: tagIonAtomType = 7;
    pub const tagIonAtomType_IAT_F: tagIonAtomType = 8;
    pub const tagIonAtomType_IAT_Cl: tagIonAtomType = 9;
    pub const tagIonAtomType_IAT_Br: tagIonAtomType = 10;
    pub const tagIonAtomType_IAT_I: tagIonAtomType = 11;
    pub const tagIonAtomType_IAT_MAX: tagIonAtomType = 12;
    pub type tagIonAtomType = u32;
    pub use self::tagIonAtomType as ION_ATOM_TYPE;
}
pub(crate) mod local_util {
    use super::*;
    pub const MIN_ATOM_CHARGE: i32 = -2;
    pub const MAX_ATOM_CHARGE: u32 = 2;
    pub const NEUTRAL_STATE: u32 = 2;
    pub const NUM_ATOM_CHARGES: u32 = 5;
    pub const MAX_NUM_VALENCES: u32 = 5;
    pub const MAXQ: u32 = 16;
    #[derive(Clone, Debug, PartialEq)]
    pub struct tagElData {
        pub szElName: SourceConstPointer<i8>,
        pub nAtMass: i32,
        pub nNormAtMass: i32,
        pub dAtMass: f64,
        pub nType: i32,
        pub nElNegPauling10: i32,
        pub bSkipAddingH: i32,
        pub cValence: [[S_CHAR; 5usize]; 5usize],
    }
    impl ::std::default::Default for tagElData {
        fn default() -> Self {
            Self {
                szElName: ::std::default::Default::default(),
                nAtMass: ::std::default::Default::default(),
                nNormAtMass: ::std::default::Default::default(),
                dAtMass: ::std::default::Default::default(),
                nType: ::std::default::Default::default(),
                nElNegPauling10: ::std::default::Default::default(),
                bSkipAddingH: ::std::default::Default::default(),
                cValence: ::std::array::from_fn(|_| {
                    ::std::array::from_fn(|_| ::std::default::Default::default())
                }),
            }
        }
    }
    pub type ELDATA = tagElData;
    pub const ERR_ELEM: i32 = 255;
    pub const nElDataLen: i32 = 122;
}
pub(crate) mod local_inchi_dll {
    use super::*;
    pub const REPEAT_ALL: u32 = 0;
    pub const MAX_MSG_LEN: u32 = 512;
    pub const bInterrupted: i32 = 0;
}
pub(crate) mod local_inchi_dll_b {
    use super::*;
    pub const ISOLATED_ATOM: i32 = -15;
    pub const MAX_CHAIN_LEN: u32 = 20;
}
pub(crate) mod local_inchi_dll_main {
    use super::*;
    pub const dummy_inchi_dll_main: i32 = 0;
}
pub(crate) mod local_ixa_builder {
    use super::*;
    pub const OPTION_PREFIX: &[u8; 2] = b"-\0";
    #[derive(Clone, Debug, PartialEq)]
    pub struct BuilderMolecule {
        pub num_atoms: AT_NUM,
        pub atom: SourceMutPointer<inchi_Atom>,
        pub num_stereo0D: AT_NUM,
        pub stereo0D: SourceMutPointer<inchi_Stereo0D>,
        pub chiral: i32,
        pub polymer: SourceMutPointer<inchi_Input_Polymer>,
        pub v3000: SourceMutPointer<inchi_Input_V3000>,
    }
    impl ::std::default::Default for BuilderMolecule {
        fn default() -> Self {
            Self {
                num_atoms: ::std::default::Default::default(),
                atom: ::std::default::Default::default(),
                num_stereo0D: ::std::default::Default::default(),
                stereo0D: ::std::default::Default::default(),
                chiral: ::std::default::Default::default(),
                polymer: ::std::default::Default::default(),
                v3000: ::std::default::Default::default(),
            }
        }
    }
    #[derive(Clone, Debug, PartialEq)]
    pub struct INCHIBUILDER {
        pub molecule: BuilderMolecule,
        pub option_WMnumber: i64,
        pub option_stereo: IXA_INCHIBUILDER_STEREOOPTION,
        pub option_NEWPSOFF: IXA_BOOL,
        pub option_DoNotAddH: IXA_BOOL,
        pub option_SUU: IXA_BOOL,
        pub option_SLUUD: IXA_BOOL,
        pub option_FixedH: IXA_BOOL,
        pub option_RecMet: IXA_BOOL,
        pub option_KET: IXA_BOOL,
        pub option_15T: IXA_BOOL,
        pub option_PT_22_00: IXA_BOOL,
        pub option_PT_16_00: IXA_BOOL,
        pub option_PT_06_00: IXA_BOOL,
        pub option_PT_39_00: IXA_BOOL,
        pub option_PT_13_00: IXA_BOOL,
        pub option_PT_18_00: IXA_BOOL,
        pub option_AuxNone: IXA_BOOL,
        pub option_WarnOnEmptyStructure: IXA_BOOL,
        pub option_LargeMolecules: IXA_BOOL,
        pub option_Polymers: IXA_BOOL,
        pub option_Polymers105: IXA_BOOL,
        pub option_Polymers105Plus: IXA_BOOL,
        pub option_FilterSS: IXA_BOOL,
        pub option_InvFilterSS: IXA_BOOL,
        pub option_NPZz: IXA_BOOL,
        pub option_SAtZz: IXA_BOOL,
        pub option_NoFrameShift: IXA_BOOL,
        pub option_FoldCRU: IXA_BOOL,
        pub option_NoEdits: IXA_BOOL,
        pub option_LooseTSACheck: IXA_BOOL,
        pub option_NoWarnings: IXA_BOOL,
        pub option_OutErrInChI: IXA_BOOL,
        pub option_SaveOpt: IXA_BOOL,
        pub option_DoDrv: IXA_BOOL,
        pub option_DoDrvReport: IXA_BOOL,
        pub option_DoR2C: IXA_BOOL,
        pub option_DoneOnly: IXA_BOOL,
        pub option_OnlyRecSalt: IXA_BOOL,
        pub option_OnlyExact: IXA_BOOL,
        pub option_OnlyRecMet: IXA_BOOL,
        pub inchi: SourceMutPointer<i8>,
        pub auxinfo: SourceMutPointer<i8>,
        pub log: SourceMutPointer<i8>,
        pub dirty: IXA_BOOL,
    }
    impl ::std::default::Default for INCHIBUILDER {
        fn default() -> Self {
            Self {
                molecule: ::std::default::Default::default(),
                option_WMnumber: ::std::default::Default::default(),
                option_stereo: ::std::default::Default::default(),
                option_NEWPSOFF: ::std::default::Default::default(),
                option_DoNotAddH: ::std::default::Default::default(),
                option_SUU: ::std::default::Default::default(),
                option_SLUUD: ::std::default::Default::default(),
                option_FixedH: ::std::default::Default::default(),
                option_RecMet: ::std::default::Default::default(),
                option_KET: ::std::default::Default::default(),
                option_15T: ::std::default::Default::default(),
                option_PT_22_00: ::std::default::Default::default(),
                option_PT_16_00: ::std::default::Default::default(),
                option_PT_06_00: ::std::default::Default::default(),
                option_PT_39_00: ::std::default::Default::default(),
                option_PT_13_00: ::std::default::Default::default(),
                option_PT_18_00: ::std::default::Default::default(),
                option_AuxNone: ::std::default::Default::default(),
                option_WarnOnEmptyStructure: ::std::default::Default::default(),
                option_LargeMolecules: ::std::default::Default::default(),
                option_Polymers: ::std::default::Default::default(),
                option_Polymers105: ::std::default::Default::default(),
                option_Polymers105Plus: ::std::default::Default::default(),
                option_FilterSS: ::std::default::Default::default(),
                option_InvFilterSS: ::std::default::Default::default(),
                option_NPZz: ::std::default::Default::default(),
                option_SAtZz: ::std::default::Default::default(),
                option_NoFrameShift: ::std::default::Default::default(),
                option_FoldCRU: ::std::default::Default::default(),
                option_NoEdits: ::std::default::Default::default(),
                option_LooseTSACheck: ::std::default::Default::default(),
                option_NoWarnings: ::std::default::Default::default(),
                option_OutErrInChI: ::std::default::Default::default(),
                option_SaveOpt: ::std::default::Default::default(),
                option_DoDrv: ::std::default::Default::default(),
                option_DoDrvReport: ::std::default::Default::default(),
                option_DoR2C: ::std::default::Default::default(),
                option_DoneOnly: ::std::default::Default::default(),
                option_OnlyRecSalt: ::std::default::Default::default(),
                option_OnlyExact: ::std::default::Default::default(),
                option_OnlyRecMet: ::std::default::Default::default(),
                inchi: ::std::default::Default::default(),
                auxinfo: ::std::default::Default::default(),
                log: ::std::default::Default::default(),
                dirty: ::std::default::Default::default(),
            }
        }
    }
}
pub(crate) mod local_ixa_inchikey_builder {
    use super::*;
    #[derive(Clone, Debug, PartialEq)]
    pub struct INCHIKEYBUILDER {
        pub inchi: SourceMutPointer<i8>,
        pub inchi_key: [i8; 256usize],
        pub dirty: i32,
    }
    impl ::std::default::Default for INCHIKEYBUILDER {
        fn default() -> Self {
            Self {
                inchi: ::std::default::Default::default(),
                inchi_key: ::std::array::from_fn(|_| ::std::default::Default::default()),
                dirty: ::std::default::Default::default(),
            }
        }
    }
}
pub(crate) mod local_ixa_status {
    use super::*;
    pub const MAXERR: u32 = 50;
    pub const MAX_MSGLENGTH: u32 = 1024;
    #[derive(Clone, Debug, PartialEq)]
    pub struct StatusItem {
        pub severity: IXA_STATUS,
        pub message: SourceMutPointer<i8>,
    }
    impl ::std::default::Default for StatusItem {
        fn default() -> Self {
            Self {
                severity: ::std::default::Default::default(),
                message: ::std::default::Default::default(),
            }
        }
    }
    #[derive(Clone, Debug, PartialEq)]
    pub struct StatusBlockTag {
        pub data: [i8; 1024usize],
        pub next: SourceMutPointer<StatusBlockTag>,
    }
    impl ::std::default::Default for StatusBlockTag {
        fn default() -> Self {
            Self {
                data: ::std::array::from_fn(|_| ::std::default::Default::default()),
                next: ::std::default::Default::default(),
            }
        }
    }
    pub type StatusBlock = StatusBlockTag;
    #[derive(Clone, Debug, PartialEq)]
    pub struct INCHISTATUS {
        pub item_count: i32,
        pub items: [StatusItem; 50usize],
        pub first_block: StatusBlock,
        pub current_block: SourceMutPointer<StatusBlock>,
        pub current_position: SourceMutPointer<i8>,
    }
    impl ::std::default::Default for INCHISTATUS {
        fn default() -> Self {
            Self {
                item_count: ::std::default::Default::default(),
                items: ::std::array::from_fn(|_| ::std::default::Default::default()),
                first_block: ::std::default::Default::default(),
                current_block: ::std::default::Default::default(),
                current_position: ::std::default::Default::default(),
            }
        }
    }
}
