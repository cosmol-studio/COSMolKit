"""P5 delivery proposal; independent p1 review and ROOT ruling are pending."""
import json
from pathlib import Path
import pytest
import cosmolkit as ck

ROOT = Path(__file__).resolve().parents[2]
PROFILE = json.loads((ROOT / "testdata/bio/gemmi_mmcif_writer_profile.json").read_text())

@pytest.mark.parametrize("case", PROFILE["cases"], ids=lambda case: case["case_id"])
def test_pinned_native_full_default_bytes(case):
    structure = ck.BioStructure.read(ROOT / case["input"])
    expected = (ROOT / "testdata/bio/expected/gemmi/bio_mmcif_writer" / case["output"]).read_bytes()
    assert structure.to_mmcif().encode() == expected
    assert structure.to_mmcif_with_params(ck.BioMmcifWriteParams()).encode() == expected
