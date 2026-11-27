#!/usr/bin/env bash
# Build a test executable from a clean RDKit checkout; never change that checkout.
set -euo pipefail
if [[ $# != 4 || ( $1 != master && $1 != exact-ipc ) ]]; then
    echo 'usage: bash build.sh {master|exact-ipc} RDKIT_SOURCE BUILD_DIR ORACLE_ROOT' >&2
    exit 2
fi
here=$(cd -- "$(dirname -- "$0")" && pwd)
kind=$1
source_dir=$(git -C "$(realpath "$2")" rev-parse --show-toplevel)
build_root=$(realpath -m "$3")
out=$(realpath -m "$4/$kind")
for dest in "$build_root" "$out"; do
    if [[ $dest == "$source_dir" || $dest == "$source_dir"/* ]]; then
        echo 'Build and output directories must be outside the supplied RDKit checkout.' >&2
        exit 1
    fi
done
revision=$(jq -er --arg kind "$kind" '.[$kind]' "$here/../../../testdata/reference/descriptor_foundations.json")
head=$(git -C "$source_dir" rev-parse HEAD)
dirty=$(git --no-optional-locks -C "$source_dir" status --porcelain)
if [[ $head != "$revision" || -n $dirty ]]; then
    echo "The $kind oracle requires a clean RDKit checkout at $revision." >&2
    exit 1
fi

# RDKit configures headers and unpacks dependencies in its source tree. Fetch
# only the pinned commit into a private repository, never the caller's ignored files.
mkdir -p "$build_root"
work=$(mktemp -d "$build_root/rdkit-XXXXXXXX")
stage=
trap 'rm -rf -- "$work"; if [[ -n $stage ]]; then rm -rf -- "$stage"; fi' EXIT
git init --quiet "$work/source"
git -C "$work/source" fetch --quiet --no-tags --depth=1 "$source_dir" "$revision"
git -C "$work/source" checkout --quiet --detach FETCH_HEAD
build_dir=$work/build
cmake_args=(-DRDK_INSTALL_INTREE=OFF -DRDK_INSTALL_STATIC_LIBS=ON)
for feature in CPP_TESTS PYTHON_WRAPPERS BOOST_PYTHON_WRAPPERS NANOBIND_WRAPPERS SWIG_WRAPPERS INCHI_SUPPORT AVALON_SUPPORT PGSQL CONTRIB DESCRIPTORS3D FREETYPE_SUPPORT CAIRO_SUPPORT QT_SUPPORT CHEMDRAW_SUPPORT COORDGEN_SUPPORT MAEPARSER_SUPPORT MOLINTERCHANGE_SUPPORT PUBCHEMSHAPE_SUPPORT SLN_SUPPORT FREESASA_SUPPORT YAEHMOP_SUPPORT STRUCTCHECKER_SUPPORT XYZ2MOL_SUPPORT COMPRESSED_SUPPLIERS; do
    cmake_args+=("-DRDK_BUILD_$feature=OFF")
done
cmake -S "$work/source" -B "$build_dir" -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PROJECT_INCLUDE="$here/rdkit.cmake" -DDESCRIPTOR_REFERENCE="$kind" \
    -DRDK_INSTALL_COMIC_FONTS=OFF "${cmake_args[@]}"
cmake --build "$build_dir" --target descriptor_oracle --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-2}"
exe=$build_dir/bin/descriptor_oracle
libs=$(ldd "$exe")
if [[ $libs == *libRDKit* || $libs == *'not found'* ]]; then
    echo 'The oracle must contain its RDKit implementation and resolve its system libraries.' >&2
    exit 1
fi
mkdir -p "$out"
stage=$(mktemp -d "$out/build-XXXXXX")
cp "$exe" "$stage/oracle"
sha=$(sha256sum "$stage/oracle" | cut -d' ' -f1)
"$stage/oracle" --build-info | jq -e --arg sha "$sha" '{binary_sha256: $sha, build_info: .}' > "$stage/identity.json"
if [[ $(jq -er '.build_info.revision' "$stage/identity.json") != "$revision" ]]; then
    echo 'The installed executable does not report the pinned RDKit revision.' >&2
    exit 1
fi
mv "$stage/oracle" "$out/$sha"
mv "$stage/identity.json" "$out/identity.json"
