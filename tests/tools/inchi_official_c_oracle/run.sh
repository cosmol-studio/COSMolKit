#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/../../.." && pwd)
source_tree="$repo_root/third_party/InChI"
source_dir="$source_tree/INCHI-1-SRC/INCHI_API/libinchi/src"
include_dir="$source_tree/INCHI-1-SRC/INCHI_BASE/src"
build_root="$repo_root/target/inchi-official-c-oracle"
library_build="$build_root/libinchi"
library_dir="$library_build/lib"
oracle_binary="$build_root/official-c-oracle"
oracle_object="$build_root/oracle.o"
object_response="$build_root/libinchi-objects.rsp"
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
find "$library_build/CMakeFiles/libinchi.dir" -type f -name '*.o' -print |
    LC_ALL=C sort > "$object_response"
cc -std=c11 -O2 -fno-strict-aliasing -Wno-all \
    -DCOMPILE_ANSI_ONLY -DTARGET_API_LIB \
    -I"$include_dir" -I"$source_dir" -I"$source_dir/ixa" \
    -c "$script_dir/oracle.c" -o "$oracle_object"
cc "$oracle_object" @"$object_response" -lm \
    -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=free \
    -o "$oracle_binary"

exec "$oracle_binary" "$@"
