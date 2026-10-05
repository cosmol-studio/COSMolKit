"""Adversarial controls for the native5000 delivery proposal comparators."""
import ast
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[2]
def comparator(family):
    path=ROOT/f"python/tests/test_canonical_{family}_native_5000.py"
    tree=ast.parse(path.read_text())
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="metadata_record")
    scope={};exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),"exec"),scope)
    return scope["metadata_record"]
class Output:
    def __init__(self,field=None,wrong_id=7):self.field=field;self.wrong_id=wrong_id
    def atom_counts(self):return [1]
    def atom_to_bits(self):return [[self.wrong_id if self.field=="atom_to_bits" else 7]]
    def bit_info_map(self):return {self.wrong_id if self.field=="bit_info_map" else 7:[(0,1)]}
    def atoms_per_bit(self):return {self.wrong_id if self.field=="atoms_per_bit" else 7:[(0,1)]}
EXPECTED={"atom_counts":[1],"atom_to_bits":[[7]],"bit_info_map":{"7":[[0,1]]},"atoms_per_bit":{"7":[[0,1]]}}
@pytest.mark.parametrize("family",["morgan","atom_pair"])
def test_source_value_positive_control(family):
    assert comparator(family)(Output())==EXPECTED
@pytest.mark.parametrize("family",["morgan","atom_pair"])
@pytest.mark.parametrize("field",["atom_to_bits","bit_info_map","atoms_per_bit"])
@pytest.mark.parametrize("wrong_id",[8,4294967303])
def test_comparison_rejects_ordinary_and_high32_wrong_ids(family,field,wrong_id):
    actual=comparator(family)(Output(field,wrong_id))
    assert actual!=EXPECTED
    assert actual[field]!=EXPECTED[field]
