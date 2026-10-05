"""P5 native switch coverage proposal, pending independent review."""
import json
import os
import hashlib
from pathlib import Path
import pytest
import cosmolkit as ck

ROOT=Path(__file__).resolve().parents[2]
MANIFEST=Path(os.environ.get("COSMOLKIT_BIO_SWITCH_CASES", ROOT/"target/bio157-native/cases-proposal.json"))
CASES=json.loads(MANIFEST.read_text())

@pytest.mark.parametrize("case",CASES,ids=lambda c:f"{c['case_id']}-{int(c['all_groups'])}-{c['flag']}")
def test_native_single_switch_bytes(case):
    structure=ck.BioStructure.read(ROOT/case["input"])
    params=ck.BioMmcifWriteParams(all_groups=case["all_groups"],**{case["flag"]:case["value"]})
    expected=(MANIFEST.parent/case["output"]).read_bytes()
    assert len(expected)==case["bytes"]
    assert hashlib.sha256(expected).hexdigest()==case["sha256"]
    assert structure.to_mmcif_with_params(params).encode()==expected
