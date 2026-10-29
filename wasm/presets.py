"""Published npm distributions, composed from Rust facade features."""
import re

PRESETS = {
    "core": ["core"],
    "core-search": ["core", "search"],
    "core-fingerprints": ["core", "fingerprints"],
    "core-analysis": ["core", "search", "fingerprints", "descriptors"],
    "core-reaction": ["core", "reaction"],
    "core-depict": ["core", "depict"],
    "core-3d": ["core", "conformer"],
    "core-bio": ["core", "bio"],
    "core-inchi": ["core", "inchi"],
    "full": ["full"],
}


def npm_release(version, preset):
    if not re.fullmatch(r"\d+\.\d+\.\d+(?:-(?:rc|dev)\.\d+)?", version):
        raise ValueError(f"Unsupported release version: {version}")
    if preset not in PRESETS:
        raise ValueError(f"Unknown WASM preset: {preset}")
    name = "@cosmol-studio/cosmolkit"
    if preset != "full":
        name += f"-{preset}"
    return name, version, version.split("-", 1)[1].split(".", 1)[0] if "-" in version else "latest"
