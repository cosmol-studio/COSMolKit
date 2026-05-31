"""Public Python API example: read an XYZ block.

Usage:
    .venv/bin/python python/examples/read_xyz.py
    .venv/bin/python python/examples/read_xyz.py path/to/molecule.xyz
"""

from __future__ import annotations

import sys
from pathlib import Path

from cosmolkit import Molecule


DEMO_XYZ = """\
3
water
O 0.000 0.000 0.000
H 0.758 0.000 0.504
H -0.758 0.000 0.504
"""


def main() -> None:
    xyz = Path(sys.argv[1]).read_text(encoding="utf-8") if len(sys.argv) > 1 else DEMO_XYZ
    mol = Molecule.from_xyz_block(xyz)

    print("atoms:", mol.num_atoms())
    print("bonds:", mol.num_bonds())
    print("3d conformers:", mol.num_conformers())
    print("coords:")
    print(mol.coords_3d())


if __name__ == "__main__":
    main()
