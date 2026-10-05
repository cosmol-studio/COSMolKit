"""Local native-oracle framing regression; does not execute an oracle."""
import importlib.util
from pathlib import Path

import pytest


def test_native_coordinate_transport_preserves_padded_end_and_escaped_bytes():
    path = Path(__file__).resolve().parents[2] / "tools/oracles/gemmi/bio_pdb_values.py"
    spec = importlib.util.spec_from_file_location("bio_reference_transport", path)
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    assert adapter.decode_output(b"81\tEND" + b" " * 77 + b"\\n\n") == "END" + " " * 77 + "\n"
    assert adapter.decode_output(b"5\t\\xc3\\xa9\\\\\\t\\n\n") == "é\\\t\n"
    for invalid in (b"2\tA\n", b"1\t\\q\n", b"1\tA\nextra\n"):
        with pytest.raises(ValueError):
            adapter.decode_output(invalid)
