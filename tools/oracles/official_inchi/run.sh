#!/bin/sh
set -eu

script_dir=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH='' cd -- "$script_dir/../../.." && pwd)
source_tree="$repo_root/third_party/InChI"
source_dir="$source_tree/INCHI-1-SRC/INCHI_API/libinchi/src"
include_dir="$source_tree/INCHI-1-SRC/INCHI_BASE/src"
build_root="$repo_root/target/inchi-official-c-oracle"
library_build="$build_root/libinchi"
oracle_binary="$build_root/official-c-oracle"
oracle_object="$build_root/oracle.o"
sha2_process_oracle_object="$build_root/sha2-process-oracle.o"
object_response="$build_root/libinchi-objects.rsp"
official_ichirvr4_source="$include_dir/ichirvr4.c"
instrumented_ichirvr4_source="$build_root/ichirvr4.oracle.c"
restored_ichirvr4_source="$build_root/ichirvr4.restored.c"
instrumented_ichirvr4_object="$build_root/ichirvr4.oracle.o"
instrumented_hooks_header="$build_root/ichirvr4.oracle-hooks.h"
instrumented_patch="$script_dir/ichirvr4_normalize_hooks.patch"
expected_tree_sha256=4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd

actual_tree_sha256=$(
    cd "$repo_root"
    find third_party/InChI -type f \
        ! -path 'third_party/InChI/.git' \
        ! -path 'third_party/InChI/.git/*' \
        -print0 |
        LC_ALL=C sort -z |
        xargs -0 sha256sum |
        sha256sum |
        awk '{print $1}'
)
if [ "$actual_tree_sha256" != "$expected_tree_sha256" ]; then
    printf '%s\n' "official InChI source provenance mismatch" >&2
    printf '%s\n' "expected: $expected_tree_sha256" >&2
    printf '%s\n' "actual:   $actual_tree_sha256" >&2
    exit 65
fi

mkdir -p "$build_root"
exec 9>"$build_root/build.lock"
flock 9

cmake -S "$source_dir" -B "$library_build" -DCMAKE_BUILD_TYPE=Release 1>&2
cmake --build "$library_build" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-2}" 1>&2
official_ichirvr4_objects=$(find "$library_build/CMakeFiles/libinchi.dir" \
    -type f -name 'ichirvr4.c.o' -print)
official_ichirvr4_object_count=$(printf '%s\n' "$official_ichirvr4_objects" |
    awk 'NF { count++ } END { print count + 0 }')
if [ "$official_ichirvr4_object_count" -ne 1 ]; then
    printf '%s\n' "expected exactly one official ichirvr4.c.o, found $official_ichirvr4_object_count" >&2
    exit 66
fi

cp "$official_ichirvr4_source" "$instrumented_ichirvr4_source"
patch --batch --fuzz=0 "$instrumented_ichirvr4_source" < "$instrumented_patch" 1>&2
fill_hook_count=$(LC_ALL=C sed -n \
    's/[^o]*oracle_test_FillOutExtraFixedHDataRestr/HOOK/p' \
    "$instrumented_ichirvr4_source" | wc -l)
fix_less_hook_count=$(LC_ALL=C sed -n \
    's/[^o]*oracle_test_FixLessHydrogenInFormula/HOOK/p' \
    "$instrumented_ichirvr4_source" | wc -l)
fix_more_hook_count=$(LC_ALL=C sed -n \
    's/[^o]*oracle_test_FixMoreHydrogenInFormula/HOOK/p' \
    "$instrumented_ichirvr4_source" | wc -l)
fix_extra_hook_count=$(LC_ALL=C sed -n \
    's/[^o]*oracle_test_FixRemoveExtraTautEndpoints/HOOK/p' \
    "$instrumented_ichirvr4_source" | wc -l)
if [ "$fill_hook_count" -ne 4 ] || [ "$fix_less_hook_count" -ne 1 ] || \
   [ "$fix_more_hook_count" -ne 1 ] || [ "$fix_extra_hook_count" -ne 1 ]; then
    printf '%s\n' "unexpected NormalizeAndCompare hook count" >&2
    exit 67
fi
sed \
    -e 's/oracle_test_FillOutExtraFixedHDataRestr/FillOutExtraFixedHDataRestr/g' \
    -e 's/oracle_test_FixLessHydrogenInFormula/FixLessHydrogenInFormula/g' \
    -e 's/oracle_test_FixMoreHydrogenInFormula/FixMoreHydrogenInFormula/g' \
    -e 's/oracle_test_FixRemoveExtraTautEndpoints/FixRemoveExtraTautEndpoints/g' \
    "$instrumented_ichirvr4_source" > "$restored_ichirvr4_source"
if ! cmp -s "$official_ichirvr4_source" "$restored_ichirvr4_source"; then
    printf '%s\n' "instrumented ichirvr4.c differs beyond verified hook tokens" >&2
    exit 68
fi
printf '%s\n' \
    'int oracle_test_FillOutExtraFixedHDataRestr();' \
    'int oracle_test_FixLessHydrogenInFormula();' \
    'int oracle_test_FixMoreHydrogenInFormula();' \
    'int oracle_test_FixRemoveExtraTautEndpoints();' \
    > "$instrumented_hooks_header"
cc -std=c11 -O3 -DNDEBUG -fPIC -g -O1 -c -fno-strict-aliasing -Wno-all \
    -DCOMPILE_ANSI_ONLY -DTARGET_API_LIB -Dlibinchi_EXPORTS \
    -I"$include_dir" -I"$source_dir" -I"$source_dir/ixa" -I"$library_build" \
    -include "$instrumented_hooks_header" \
    -c "$instrumented_ichirvr4_source" -o "$instrumented_ichirvr4_object"
find "$library_build/CMakeFiles/libinchi.dir" -type f -name '*.o' -print |
    awk -v excluded="$official_ichirvr4_objects" '$0 != excluded' |
    LC_ALL=C sort > "$object_response"
cc -std=c11 -O2 -fno-strict-aliasing -Wno-all \
    -DCOMPILE_ANSI_ONLY -DTARGET_API_LIB \
    -I"$include_dir" -I"$source_dir" -I"$source_dir/ixa" \
    -c "$script_dir/oracle.c" -o "$oracle_object"
cc -std=c11 -O2 -fno-strict-aliasing -Wno-all \
    -DCOMPILE_ANSI_ONLY -DTARGET_API_LIB \
    -I"$include_dir" -I"$source_dir" -I"$source_dir/ixa" \
    -c "$script_dir/sha2_process_oracle.c" -o "$sha2_process_oracle_object"
cc "$oracle_object" "$sha2_process_oracle_object" \
    "$instrumented_ichirvr4_object" @"$object_response" -lm \
    -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc -Wl,--wrap=free \
    -Wl,--wrap=MakeOneInChIOutOfStrFromINChI2 \
    -Wl,--wrap=CompareReversedINChI2 \
    -Wl,--wrap=FixFixedHRestoredStructure \
    -Wl,--wrap=FixMobileHRestoredStructure \
    -Wl,--wrap=FixRestoredStructureStereo \
    -Wl,--wrap=inchi_strbuf_init -Wl,--wrap=MergeZzInHillFormula \
    -Wl,--wrap=inchi_strbuf_close \
    -Wl,--wrap=Free_INChI -Wl,--wrap=Free_INChI_Aux \
    -Wl,--wrap=FreeInpAtomData -Wl,--wrap=free_t_group_info \
    -o "$oracle_binary"

exec "$oracle_binary" "$@"
