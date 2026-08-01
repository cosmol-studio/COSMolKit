#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/../../.." && pwd)
rdkit_root="$repo_root/third_party/rdkit"
source_file="$rdkit_root/External/INCHI-API/inchi.cpp"
header_file="$rdkit_root/External/INCHI-API/inchi.h"
build_root="$repo_root/target/rdkit-inchi-cpp-oracle"
assign_bond_dirs_fragment="$build_root/assign_bond_dirs.inc"
find_alternating_bonds_fragment="$build_root/find_alternating_bonds.inc"
neighboring_si_fragment="$build_root/neighboring_si.inc"
valence4n_cleanup1_fragment="$build_root/valence4n_cleanup1.inc"
valence4n_cleanup2_fragment="$build_root/valence4n_cleanup2.inc"
valence5n_cleanup1_fragment="$build_root/valence5n_cleanup1.inc"
valence5n_cleanup2_fragment="$build_root/valence5n_cleanup2.inc"
valence5n_cleanup3_fragment="$build_root/valence5n_cleanup3.inc"
valence5n_cleanup4_fragment="$build_root/valence5n_cleanup4.inc"
valence5n_cleanup5_fragment="$build_root/valence5n_cleanup5.inc"
valence5n_cleanup6_fragment="$build_root/valence5n_cleanup6.inc"
valence5n_cleanup7_fragment="$build_root/valence5n_cleanup7.inc"
valence5n_cleanup8_fragment="$build_root/valence5n_cleanup8.inc"
valence5n_cleanup9_fragment="$build_root/valence5n_cleanup9.inc"
valence5n_cleanupa_fragment="$build_root/valence5n_cleanupa.inc"
valence5n_cleanupb_fragment="$build_root/valence5n_cleanupb.inc"
valence7s_cleanup1_fragment="$build_root/valence7s_cleanup1.inc"
valence7s_cleanup2_fragment="$build_root/valence7s_cleanup2.inc"
valence7s_cleanup3_fragment="$build_root/valence7s_cleanup3.inc"
valence8s_cleanup1_fragment="$build_root/valence8s_cleanup1.inc"
valence8cl_cleanup1_fragment="$build_root/valence8cl_cleanup1.inc"
valence5cl_cleanup1_fragment="$build_root/valence5cl_cleanup1.inc"
valence3cl_cleanup1_fragment="$build_root/valence3cl_cleanup1.inc"
clean_up_fragment="$build_root/clean_up.inc"
inchi_to_mol_fragment="$build_root/inchi_to_mol.inc"
fix_option_symbol_fragment="$build_root/fix_option_symbol.inc"
r_clean_up_fragment="$build_root/r_clean_up.inc"
mol_to_inchi_fragment="$build_root/mol_to_inchi.inc"
inchi_to_inchi_key_fragment="$build_root/inchi_to_inchi_key.inc"
mol_to_inchi_key_fragment="$build_root/mol_to_inchi_key.inc"
oracle_binary="$build_root/rdkit-inchi-cpp-oracle"

check_sha256() {
    expected=$1
    path=$2
    actual=$(sha256sum "$path" | awk '{print $1}')
    if [ "$actual" != "$expected" ]; then
        printf '%s\n' "pinned RDKit source provenance mismatch: $path" >&2
        printf '%s\n' "expected: $expected" >&2
        printf '%s\n' "actual:   $actual" >&2
        exit 65
    fi
}

check_sha256 104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f "$source_file"
check_sha256 27b4eaef1714869c42dbf2998807018a03389b4f9ce40438e843248ebfc3614e "$header_file"
check_sha256 15f873f32dac816103b77a28e6f45e651deffca89637fc8b869f51bf4f1c35e1 "$rdkit_root/Code/GraphMol/Bond.h"
check_sha256 e7356f312f5db00fba6a6afe876ffd84ef80928a90a2734b00a46102882fdddb "$rdkit_root/Code/GraphMol/ROMol.cpp"
check_sha256 ffc389fcbe388012b641e59ab72e24c4da28a1ab65a3f7aef23f471add54356a "$rdkit_root/Code/RDGeneral/Invariant.h"
check_sha256 4a41840314dbb5c5f836140fd87ef6d71f0f932659e7e01ab5be6984dac1d647 "$rdkit_root/CMakeLists.txt"
check_sha256 bfa51d29a6b18c4bfa373b2803ee79499b7ffbfa40d5838efb273bba2078961e "$rdkit_root/Code/GraphMol/Atom.cpp"
check_sha256 a9c2829c484c98df285fa503870351d1ac8d947cfe46ad38f3cb6664f1c9ede2 "$rdkit_root/Code/GraphMol/Substruct/SubstructMatch.cpp"
check_sha256 9c0a6d4c56f8fcbc8b4c3144ffd2893547780ebba5602858b2c5415beaaca1cc "$rdkit_root/Code/GraphMol/Substruct/SubstructMatch.h"
check_sha256 76dbf84efc3b628f206de39b213373d55031e2762b116026c2fd05581af95c5f "$rdkit_root/Code/GraphMol/Substruct/SubstructUtils.cpp"
check_sha256 270b0457fbee294a9a1a1d7ada856538d4faec96fd782ffff9d2c5bb9ec0ca59 "$rdkit_root/Code/GraphMol/SmilesParse/SmilesParse.cpp"

mkdir -p "$build_root"
exec 9>"$build_root/build.lock"
flock 9

sed -n '89,168p' "$source_file" > "$assign_bond_dirs_fragment"
check_sha256 582b6396ac3c3952878f88b6b5072f0dfd061b16448630b7df02d5a723151e7a "$assign_bond_dirs_fragment"
sed -n '207,304p' "$source_file" > "$find_alternating_bonds_fragment"
check_sha256 ef0a05d3a018d6928d66af0d5238dea5bb22742e233f1fe8229f62578fe7d003 "$find_alternating_bonds_fragment"
sed -n '306,321p' "$source_file" > "$neighboring_si_fragment"
check_sha256 a5d70ac89adb467c691ba0a8a16878e951833d598c7f7d15a1a9d241f0c17374 "$neighboring_si_fragment"
sed -n '324,373p' "$source_file" > "$valence4n_cleanup1_fragment"
check_sha256 11132300ed0a447a05c38e8eed0e710762514b9180aad984f12acc697255d87d "$valence4n_cleanup1_fragment"
sed -n '376,390p' "$source_file" > "$valence4n_cleanup2_fragment"
check_sha256 c2fe8b4125473236e4fa1c638d93c23365c1c7df416f17537259dcc24c4d7587 "$valence4n_cleanup2_fragment"
sed -n '393,414p' "$source_file" > "$valence5n_cleanup1_fragment"
check_sha256 58b4273fc255b5ed398a0d81b30069e58d68478d083b499f41a77ec508250f06 "$valence5n_cleanup1_fragment"
sed -n '417,440p' "$source_file" > "$valence5n_cleanup2_fragment"
check_sha256 07654fddc42017f1aa76c39ec11013301f92e460b199218c6a2f3bf6bcfb4484 "$valence5n_cleanup2_fragment"
sed -n '443,473p' "$source_file" > "$valence5n_cleanup3_fragment"
check_sha256 bf11af3caae7ff7d1dd7b381defd449e6cb9b193e9dab328a35a6ddd5d73e25d "$valence5n_cleanup3_fragment"
sed -n '477,542p' "$source_file" > "$valence5n_cleanup4_fragment"
check_sha256 701347673195ee11dec11ba3013f6eb10746d62e0a879ce63d0de98ec0642746 "$valence5n_cleanup4_fragment"
sed -n '544,621p' "$source_file" > "$valence5n_cleanup5_fragment"
check_sha256 572b129f3625cbffe00811214403e660abdcf283c582ff90287011755698b73e "$valence5n_cleanup5_fragment"
sed -n '625,675p' "$source_file" > "$valence5n_cleanup6_fragment"
check_sha256 33f731b27422f52bbe7b9304df66f3126a325a3458a577bcf9b29c01d4b9f47d "$valence5n_cleanup6_fragment"
sed -n '679,746p' "$source_file" > "$valence5n_cleanup7_fragment"
check_sha256 52035cacfe7ced70b44df5ac007a41f101a968fbee2a5338abf5ff81a4efa562 "$valence5n_cleanup7_fragment"
sed -n '750,800p' "$source_file" > "$valence5n_cleanup8_fragment"
check_sha256 9874b6974fcf056f0a7313671fcb32d818b00e9d66c7e11934a5919569440c65 "$valence5n_cleanup8_fragment"
sed -n '804,852p' "$source_file" > "$valence5n_cleanup9_fragment"
check_sha256 c2b7f87dfdd510b485982edfc6a4f0327ca042e86b5e35eee7484e2193096ef0 "$valence5n_cleanup9_fragment"
sed -n '855,913p' "$source_file" > "$valence5n_cleanupa_fragment"
check_sha256 2d4518bea1a8abb8d5736f1c784fa7e30cb033950dcd6ba5da13846bb46d6392 "$valence5n_cleanupa_fragment"
sed -n '916,938p' "$source_file" > "$valence5n_cleanupb_fragment"
check_sha256 9d097bf179047ec55e31d8f45ffe573f449e3d5828833e5430d2faf355a5defc "$valence5n_cleanupb_fragment"
sed -n '941,986p' "$source_file" > "$valence7s_cleanup1_fragment"
check_sha256 a03efbad1fb09fa47513380a625fefccda34d23e3d4f5fd0b2f2b2cf3f8d076e "$valence7s_cleanup1_fragment"
sed -n '989,1018p' "$source_file" > "$valence7s_cleanup2_fragment"
check_sha256 d8e8ca7577874896efc2f1df36738d0bee34a8dcabb35493be0095017a2e9618 "$valence7s_cleanup2_fragment"
sed -n '1021,1041p' "$source_file" > "$valence7s_cleanup3_fragment"
check_sha256 b61785dfecb61f3b567183f96ec1290ef56f576436bfcd07a808907ebbf4b4ed "$valence7s_cleanup3_fragment"
sed -n '1044,1074p' "$source_file" > "$valence8s_cleanup1_fragment"
check_sha256 e8bda5f9c53e83b3f5a0cce7f26861384606e709e5ef66e7adedac5c0ebcfd78 "$valence8s_cleanup1_fragment"
sed -n '1079,1111p' "$source_file" > "$valence8cl_cleanup1_fragment"
check_sha256 8e97a538eab77c78c2fb02be47d47cf6be989761cc280cf4857b9500885b1f3a "$valence8cl_cleanup1_fragment"
sed -n '1114,1131p' "$source_file" > "$valence5cl_cleanup1_fragment"
check_sha256 102c9d4a05c61706c9b0de30d61fc725fc394f80d532810d0b94b0103d6ff250 "$valence5cl_cleanup1_fragment"
sed -n '1134,1149p' "$source_file" > "$valence3cl_cleanup1_fragment"
check_sha256 078efac5f5acece4f6323b171f6aab2b1795873b55a505e56fb44aa54fa7d982 "$valence3cl_cleanup1_fragment"
sed -n '1151,1251p' "$source_file" > "$clean_up_fragment"
check_sha256 25458d16f3e2888aed0a60ba3b1c8f0a4e435a912215af7da93fbc491b34b371 "$clean_up_fragment"
sed -n '1254,1672p' "$source_file" > "$inchi_to_mol_fragment"
check_sha256 3765efd5f4f8a855a2c2231c1d62a81e7d5880da76804f47f400ed8f5a6464a1 "$inchi_to_mol_fragment"
sed -n '1674,1691p' "$source_file" > "$fix_option_symbol_fragment"
check_sha256 5770d6ab655210bd67e5351ef3ebe39aadc0c8b8e96c1b6e2ec0935acb8844d9 "$fix_option_symbol_fragment"
sed -n '1694,1745p' "$source_file" > "$r_clean_up_fragment"
check_sha256 81165fe80fd36e56b94f44fc7f062b0971bdd86d96fec1ebbdd947c954c64858 "$r_clean_up_fragment"
sed -n '1747,2104p' "$source_file" > "$mol_to_inchi_fragment"
check_sha256 000d4eff4353e06aaa802057674fa24f91726872b0c0f9bb593d7a1937721ac2 "$mol_to_inchi_fragment"
sed -n '2145,2175p' "$source_file" > "$inchi_to_inchi_key_fragment"
check_sha256 b2d19205873e8444cef0ce660bbcf36bbf35cc528e77802b40b9b28ec15e67ae "$inchi_to_inchi_key_fragment"
sed -n '107,110p' "$header_file" > "$mol_to_inchi_key_fragment"
check_sha256 f6cc64a50170bb62e6c7a97325815f1399a6d3c74a9bc1b4a8bf69da90ccae51 "$mol_to_inchi_key_fragment"

c++ -std=c++17 -O2 -Wall -Wextra -Werror -Wno-unused-variable \
    -I"$repo_root/third_party/InChI/INCHI-1-SRC/INCHI_BASE/src" \
    -DASSIGN_BOND_DIRS_SOURCE='"'"$assign_bond_dirs_fragment"'"' \
    -DFIND_ALTERNATING_BONDS_SOURCE='"'"$find_alternating_bonds_fragment"'"' \
    -DNEIGHBORING_SI_SOURCE='"'"$neighboring_si_fragment"'"' \
    -DVALENCE4N_CLEANUP1_SOURCE='"'"$valence4n_cleanup1_fragment"'"' \
    -DVALENCE4N_CLEANUP2_SOURCE='"'"$valence4n_cleanup2_fragment"'"' \
    -DVALENCE5N_CLEANUP1_SOURCE='"'"$valence5n_cleanup1_fragment"'"' \
    -DVALENCE5N_CLEANUP2_SOURCE='"'"$valence5n_cleanup2_fragment"'"' \
    -DVALENCE5N_CLEANUP3_SOURCE='"'"$valence5n_cleanup3_fragment"'"' \
    -DVALENCE5N_CLEANUP4_SOURCE='"'"$valence5n_cleanup4_fragment"'"' \
    -DVALENCE5N_CLEANUP5_SOURCE='"'"$valence5n_cleanup5_fragment"'"' \
    -DVALENCE5N_CLEANUP6_SOURCE='"'"$valence5n_cleanup6_fragment"'"' \
    -DVALENCE5N_CLEANUP7_SOURCE='"'"$valence5n_cleanup7_fragment"'"' \
    -DVALENCE5N_CLEANUP8_SOURCE='"'"$valence5n_cleanup8_fragment"'"' \
    -DVALENCE5N_CLEANUP9_SOURCE='"'"$valence5n_cleanup9_fragment"'"' \
    -DVALENCE5N_CLEANUPA_SOURCE='"'"$valence5n_cleanupa_fragment"'"' \
    -DVALENCE5N_CLEANUPB_SOURCE='"'"$valence5n_cleanupb_fragment"'"' \
    -DVALENCE7S_CLEANUP1_SOURCE='"'"$valence7s_cleanup1_fragment"'"' \
    -DVALENCE7S_CLEANUP2_SOURCE='"'"$valence7s_cleanup2_fragment"'"' \
    -DVALENCE7S_CLEANUP3_SOURCE='"'"$valence7s_cleanup3_fragment"'"' \
    -DVALENCE8S_CLEANUP1_SOURCE='"'"$valence8s_cleanup1_fragment"'"' \
    -DVALENCE8CL_CLEANUP1_SOURCE='"'"$valence8cl_cleanup1_fragment"'"' \
    -DVALENCE5CL_CLEANUP1_SOURCE='"'"$valence5cl_cleanup1_fragment"'"' \
    -DVALENCE3CL_CLEANUP1_SOURCE='"'"$valence3cl_cleanup1_fragment"'"' \
    -DCLEAN_UP_SOURCE='"'"$clean_up_fragment"'"' \
    -DINCHI_TO_MOL_SOURCE='"'"$inchi_to_mol_fragment"'"' \
    -DFIX_OPTION_SYMBOL_SOURCE='"'"$fix_option_symbol_fragment"'"' \
    -DR_CLEAN_UP_SOURCE='"'"$r_clean_up_fragment"'"' \
    -DMOL_TO_INCHI_SOURCE='"'"$mol_to_inchi_fragment"'"' \
    -DINCHI_TO_INCHI_KEY_SOURCE='"'"$inchi_to_inchi_key_fragment"'"' \
    -DMOL_TO_INCHI_KEY_SOURCE='"'"$mol_to_inchi_key_fragment"'"' \
    "$script_dir/oracle.cpp" -o "$oracle_binary"

exec "$oracle_binary" "$@"
