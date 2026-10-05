"""Fixed source metadata and native value transport, pin 351f8f378f8ad6bbd517980c38896e66bf907af8.

119 literal rows match accepted facade literals and pinned atomic_data.cpp
SHA256 7f9cee6e430b60d303a0a7fa9e33c45e5c529ee204f86afeab6a20f68b6b0631.
Tests do not read chemistry source or derive expectations from the owner.
"""

import ast
import inspect
import struct
from pathlib import Path
from typing import cast

import cosmolkit as ck
import pytest

Row = tuple[str, str, int, int, int, list[int], float, float, str, str]
ROWS: list[Row] = [
    ('DUMMY', '*', 0, 0, 0, [-1], 0.0, 0.0, '0000000000000000', '0000000000000000'),
    ('H', 'H', 1, 1, 1, [1], 0.33, 1.008, '3fd51eb851eb851f', '3ff020c49ba5e354'),
    ('HE', 'He', 2, 1, 2, [0], 0.7, 4.003, '3fe6666666666666', '401003126e978d50'),
    ('LI', 'Li', 3, 2, 1, [1, -1], 1.23, 6.941, '3ff3ae147ae147ae', '401bc395810624dd'),
    ('BE', 'Be', 4, 2, 2, [2], 0.9, 9.012, '3feccccccccccccd', '40220624dd2f1aa0'),
    ('B', 'B', 5, 2, 3, [3], 0.82, 10.812, '3fea3d70a3d70a3d', '40259fbe76c8b439'),
    ('C', 'C', 6, 2, 4, [4], 0.77, 12.011, '3fe8a3d70a3d70a4', '402805a1cac08312'),
    ('N', 'N', 7, 2, 5, [3], 0.7, 14.007, '3fe6666666666666', '402c0395810624dd'),
    ('O', 'O', 8, 2, 6, [2], 0.66, 15.999, '3fe51eb851eb851f', '402fff7ced916873'),
    ('F', 'F', 9, 2, 7, [1], 0.611, 18.998, '3fe38d4fdf3b645a', '4032ff7ced916873'),
    ('NE', 'Ne', 10, 2, 8, [0], 0.7, 20.18, '3fe6666666666666', '40342e147ae147ae'),
    ('NA', 'Na', 11, 3, 1, [1, -1], 1.54, 22.99, '3ff8a3d70a3d70a4', '4036fd70a3d70a3d'),
    ('MG', 'Mg', 12, 3, 2, [2, -1], 1.36, 24.305, '3ff5c28f5c28f5c3', '40384e147ae147ae'),
    ('AL', 'Al', 13, 3, 3, [3], 1.18, 26.982, '3ff2e147ae147ae1', '403afb645a1cac08'),
    ('SI', 'Si', 14, 3, 4, [4], 0.937, 28.086, '3fedfbe76c8b4396', '403c1604189374bc'),
    ('P', 'P', 15, 3, 5, [3, 5], 0.89, 30.974, '3fec7ae147ae147b', '403ef95810624dd3'),
    ('S', 'S', 16, 3, 6, [2, 4, 6], 1.04, 32.067, '3ff0a3d70a3d70a4', '4040089374bc6a7f'),
    ('CL', 'Cl', 17, 3, 7, [1], 0.997, 35.453, '3fefe76c8b439581', '4041b9fbe76c8b44'),
    ('AR', 'Ar', 18, 3, 8, [0], 1.74, 39.948, '3ffbd70a3d70a3d7', '4043f95810624dd3'),
    ('K', 'K', 19, 4, 1, [1, -1], 2.03, 39.098, '40003d70a3d70a3d', '40438c8b43958106'),
    ('CA', 'Ca', 20, 4, 2, [2, -1], 1.74, 40.078, '3ffbd70a3d70a3d7', '404409fbe76c8b44'),
    ('SC', 'Sc', 21, 4, 3, [-1], 1.44, 44.956, '3ff70a3d70a3d70a', '40467a5e353f7cee'),
    ('TI', 'Ti', 22, 4, 4, [-1], 1.32, 47.867, '3ff51eb851eb851f', '4047eef9db22d0e5'),
    ('V', 'V', 23, 4, 5, [-1], 1.22, 50.944, '3ff3851eb851eb85', '404978d4fdf3b646'),
    ('CR', 'Cr', 24, 4, 6, [-1], 1.18, 51.996, '3ff2e147ae147ae1', '4049ff7ced916873'),
    ('MN', 'Mn', 25, 4, 7, [-1], 1.17, 54.938, '3ff2b851eb851eb8', '404b7810624dd2f2'),
    ('FE', 'Fe', 26, 4, 8, [-1], 1.17, 55.845, '3ff2b851eb851eb8', '404bec28f5c28f5c'),
    ('CO', 'Co', 27, 4, 9, [-1], 1.16, 58.933, '3ff28f5c28f5c28f', '404d776c8b439581'),
    ('NI', 'Ni', 28, 4, 10, [-1], 1.15, 58.693, '3ff2666666666666', '404d58b439581062'),
    ('CU', 'Cu', 29, 4, 11, [-1], 1.17, 63.546, '3ff2b851eb851eb8', '404fc5e353f7ced9'),
    ('ZN', 'Zn', 30, 4, 2, [-1], 1.25, 65.39, '3ff4000000000000', '405058f5c28f5c29'),
    ('GA', 'Ga', 31, 4, 3, [3], 1.26, 69.723, '3ff428f5c28f5c29', '40516e45a1cac083'),
    ('GE', 'Ge', 32, 4, 4, [4], 1.188, 72.61, '3ff3020c49ba5e35', '4052270a3d70a3d7'),
    ('AS', 'As', 33, 4, 5, [3, 5], 1.2, 74.922, '3ff3333333333333', '4052bb020c49ba5e'),
    ('SE', 'Se', 34, 4, 6, [2, 4, 6], 1.17, 78.96, '3ff2b851eb851eb8', '4053bd70a3d70a3d'),
    ('BR', 'Br', 35, 4, 7, [1], 1.167, 79.904, '3ff2ac083126e979', '4053f9db22d0e560'),
    ('KR', 'Kr', 36, 4, 8, [0], 1.91, 83.8, '3ffe8f5c28f5c28f', '4054f33333333333'),
    ('RB', 'Rb', 37, 5, 1, [1, -1], 2.16, 85.468, '400147ae147ae148', '40555df3b645a1cb'),
    ('SR', 'Sr', 38, 5, 2, [2, -1], 1.91, 87.62, '3ffe8f5c28f5c28f', '4055e7ae147ae148'),
    ('Y', 'Y', 39, 5, 3, [-1], 1.62, 88.906, '3ff9eb851eb851ec', '405639fbe76c8b44'),
    ('ZR', 'Zr', 40, 5, 4, [-1], 1.45, 91.224, '3ff7333333333333', '4056ce5604189375'),
    ('NB', 'Nb', 41, 5, 5, [-1], 1.34, 92.906, '3ff570a3d70a3d71', '405739fbe76c8b44'),
    ('MO', 'Mo', 42, 5, 6, [-1], 1.3, 95.94, '3ff4cccccccccccd', '4057fc28f5c28f5c'),
    ('TC', 'Tc', 43, 5, 7, [-1], 1.27, 98.0, '3ff451eb851eb852', '4058800000000000'),
    ('RU', 'Ru', 44, 5, 8, [-1], 1.25, 101.07, '3ff4000000000000', '4059447ae147ae14'),
    ('RH', 'Rh', 45, 5, 9, [-1], 1.25, 102.906, '3ff4000000000000', '4059b9fbe76c8b44'),
    ('PD', 'Pd', 46, 5, 10, [-1], 1.28, 106.42, '3ff47ae147ae147b', '405a9ae147ae147b'),
    ('AG', 'Ag', 47, 5, 11, [-1], 1.34, 107.868, '3ff570a3d70a3d71', '405af78d4fdf3b64'),
    ('CD', 'Cd', 48, 5, 2, [-1], 1.48, 112.412, '3ff7ae147ae147ae', '405c1a5e353f7cee'),
    ('IN', 'In', 49, 5, 3, [3], 1.44, 114.818, '3ff70a3d70a3d70a', '405cb45a1cac0831'),
    ('SN', 'Sn', 50, 5, 4, [2, 4], 1.385, 118.711, '3ff628f5c28f5c29', '405dad810624dd2f'),
    ('SB', 'Sb', 51, 5, 5, [3, 5], 1.4, 121.76, '3ff6666666666666', '405e70a3d70a3d71'),
    ('TE', 'Te', 52, 5, 6, [2, 4, 6], 1.378, 127.6, '3ff60c49ba5e353f', '405fe66666666666'),
    ('I', 'I', 53, 5, 7, [1, 3, 5], 1.387, 126.904, '3ff63126e978d4fe', '405fb9db22d0e560'),
    ('XE', 'Xe', 54, 5, 8, [0, 2, 4, 6], 1.98, 131.29, '3fffae147ae147ae', '40606947ae147ae1'),
    ('CS', 'Cs', 55, 6, 1, [1], 2.35, 132.905, '4002cccccccccccd', '40609cf5c28f5c29'),
    ('BA', 'Ba', 56, 6, 2, [2, -1], 1.98, 137.328, '3fffae147ae147ae', '40612a7ef9db22d1'),
    ('LA', 'La', 57, 6, 3, [-1], 1.69, 138.906, '3ffb0a3d70a3d70a', '40615cfdf3b645a2'),
    ('CE', 'Ce', 58, 6, 4, [-1], 1.83, 140.116, '3ffd47ae147ae148', '406183b645a1cac1'),
    ('PR', 'Pr', 59, 6, 3, [-1], 1.82, 140.908, '3ffd1eb851eb851f', '40619d0e56041893'),
    ('ND', 'Nd', 60, 6, 4, [-1], 1.81, 144.24, '3ffcf5c28f5c28f6', '406207ae147ae148'),
    ('PM', 'Pm', 61, 6, 5, [-1], 1.8, 145.0, '3ffccccccccccccd', '4062200000000000'),
    ('SM', 'Sm', 62, 6, 6, [-1], 1.8, 150.36, '3ffccccccccccccd', '4062cb851eb851ec'),
    ('EU', 'Eu', 63, 6, 7, [-1], 1.99, 151.964, '3fffd70a3d70a3d7', '4062fed916872b02'),
    ('GD', 'Gd', 64, 6, 8, [-1], 1.79, 157.25, '3ffca3d70a3d70a4', '4063a80000000000'),
    ('TB', 'Tb', 65, 6, 9, [-1], 1.76, 158.925, '3ffc28f5c28f5c29', '4063dd999999999a'),
    ('DY', 'Dy', 66, 6, 10, [-1], 1.75, 162.5, '3ffc000000000000', '4064500000000000'),
    ('HO', 'Ho', 67, 6, 11, [-1], 1.74, 164.93, '3ffbd70a3d70a3d7', '40649dc28f5c28f6'),
    ('ER', 'Er', 68, 6, 12, [-1], 1.73, 167.26, '3ffbae147ae147ae', '4064e851eb851eb8'),
    ('TM', 'Tm', 69, 6, 13, [-1], 1.72, 168.934, '3ffb851eb851eb85', '40651de353f7ced9'),
    ('YB', 'Yb', 70, 6, 14, [-1], 1.94, 173.04, '3fff0a3d70a3d70a', '4065a147ae147ae1'),
    ('LU', 'Lu', 71, 6, 15, [-1], 1.72, 174.967, '3ffb851eb851eb85', '4065def1a9fbe76d'),
    ('HF', 'Hf', 72, 6, 4, [-1], 1.44, 178.49, '3ff70a3d70a3d70a', '40664fae147ae148'),
    ('TA', 'Ta', 73, 6, 5, [-1], 1.34, 180.948, '3ff570a3d70a3d71', '40669e5604189375'),
    ('W', 'W', 74, 6, 6, [-1], 1.3, 183.84, '3ff4cccccccccccd', '4066fae147ae147b'),
    ('RE', 'Re', 75, 6, 7, [-1], 1.28, 186.207, '3ff47ae147ae147b', '4067469fbe76c8b4'),
    ('OS', 'Os', 76, 6, 8, [-1], 1.26, 190.23, '3ff428f5c28f5c29', '4067c75c28f5c28f'),
    ('IR', 'Ir', 77, 6, 9, [-1], 1.27, 192.217, '3ff451eb851eb852', '406806f1a9fbe76d'),
    ('PT', 'Pt', 78, 6, 10, [-1], 1.3, 195.078, '3ff4cccccccccccd', '4068627ef9db22d1'),
    ('AU', 'Au', 79, 6, 11, [-1], 1.34, 196.967, '3ff570a3d70a3d71', '40689ef1a9fbe76d'),
    ('HG', 'Hg', 80, 6, 2, [-1], 1.49, 200.59, '3ff7d70a3d70a3d7', '406912e147ae147b'),
    ('TL', 'Tl', 81, 6, 3, [-1], 1.48, 204.383, '3ff7ae147ae147ae', '40698c4189374bc7'),
    ('PB', 'Pb', 82, 6, 4, [2, 4], 1.48, 207.2, '3ff7ae147ae147ae', '4069e66666666666'),
    ('BI', 'Bi', 83, 6, 5, [3, 5], 1.45, 208.98, '3ff7333333333333', '406a1f5c28f5c28f'),
    ('PO', 'Po', 84, 6, 6, [2, 4, 6], 1.46, 209.0, '3ff75c28f5c28f5c', '406a200000000000'),
    ('AT', 'At', 85, 6, 7, [1, 3, 5], 1.45, 210.0, '3ff7333333333333', '406a400000000000'),
    ('RN', 'Rn', 86, 6, 8, [0], 2.4, 222.0, '4003333333333333', '406bc00000000000'),
    ('FR', 'Fr', 87, 7, 1, [1], 2.0, 223.0, '4000000000000000', '406be00000000000'),
    ('RA', 'Ra', 88, 7, 2, [2, -1], 1.9, 226.0, '3ffe666666666666', '406c400000000000'),
    ('AC', 'Ac', 89, 7, 3, [-1], 1.88, 227.0, '3ffe147ae147ae14', '406c600000000000'),
    ('TH', 'Th', 90, 7, 4, [-1], 1.79, 232.038, '3ffca3d70a3d70a4', '406d01374bc6a7f0'),
    ('PA', 'Pa', 91, 7, 3, [-1], 1.61, 231.036, '3ff9c28f5c28f5c3', '406ce126e978d4fe'),
    ('U', 'U', 92, 7, 4, [-1], 1.58, 238.029, '3ff947ae147ae148', '406dc0ed916872b0'),
    ('NP', 'Np', 93, 7, 5, [-1], 1.55, 237.0, '3ff8cccccccccccd', '406da00000000000'),
    ('PU', 'Pu', 94, 7, 6, [-1], 1.53, 244.0, '3ff87ae147ae147b', '406e800000000000'),
    ('AM', 'Am', 95, 7, 7, [-1], 1.07, 243.0, '3ff11eb851eb851f', '406e600000000000'),
    ('CM', 'Cm', 96, 7, 8, [-1], 0.0, 247.0, '0000000000000000', '406ee00000000000'),
    ('BK', 'Bk', 97, 7, 9, [-1], 0.0, 247.0, '0000000000000000', '406ee00000000000'),
    ('CF', 'Cf', 98, 7, 10, [-1], 0.0, 251.0, '0000000000000000', '406f600000000000'),
    ('ES', 'Es', 99, 7, 11, [-1], 0.0, 252.0, '0000000000000000', '406f800000000000'),
    ('FM', 'Fm', 100, 7, 12, [-1], 0.0, 257.0, '0000000000000000', '4070100000000000'),
    ('MD', 'Md', 101, 7, 13, [-1], 0.0, 258.0, '0000000000000000', '4070200000000000'),
    ('NO', 'No', 102, 7, 14, [-1], 0.0, 259.0, '0000000000000000', '4070300000000000'),
    ('LR', 'Lr', 103, 7, 15, [-1], 0.0, 262.0, '0000000000000000', '4070600000000000'),
    ('RF', 'Rf', 104, 7, 2, [-1], 0.0, 267.0, '0000000000000000', '4070b00000000000'),
    ('DB', 'Db', 105, 7, 2, [-1], 0.0, 268.0, '0000000000000000', '4070c00000000000'),
    ('SG', 'Sg', 106, 7, 2, [-1], 0.0, 269.0, '0000000000000000', '4070d00000000000'),
    ('BH', 'Bh', 107, 7, 2, [-1], 0.0, 270.0, '0000000000000000', '4070e00000000000'),
    ('HS', 'Hs', 108, 7, 2, [-1], 0.0, 269.0, '0000000000000000', '4070d00000000000'),
    ('MT', 'Mt', 109, 7, 2, [-1], 0.0, 278.0, '0000000000000000', '4071600000000000'),
    ('DS', 'Ds', 110, 7, 2, [-1], 0.0, 281.0, '0000000000000000', '4071900000000000'),
    ('RG', 'Rg', 111, 7, 2, [-1], 0.0, 281.0, '0000000000000000', '4071900000000000'),
    ('CN', 'Cn', 112, 7, 2, [-1], 0.0, 285.0, '0000000000000000', '4071d00000000000'),
    ('NH', 'Nh', 113, 7, 2, [-1], 0.0, 284.0, '0000000000000000', '4071c00000000000'),
    ('FL', 'Fl', 114, 7, 2, [-1], 0.0, 289.0, '0000000000000000', '4072100000000000'),
    ('MC', 'Mc', 115, 7, 2, [-1], 0.0, 288.0, '0000000000000000', '4072000000000000'),
    ('LV', 'Lv', 116, 7, 2, [-1], 0.0, 293.0, '0000000000000000', '4072500000000000'),
    ('TS', 'Ts', 117, 7, 2, [-1], 0.0, 292.0, '0000000000000000', '4072400000000000'),
    ('OG', 'Og', 118, 7, 2, [-1], 0.0, 294.0, '0000000000000000', '4072600000000000'),
]

INFO_READS = ["element", "symbol", "atomic_number", "period", "outer_electrons", "valences", "rb0", "atomic_weight"]
ELEMENT_READS = ["atomic_number", "symbol"]
STUB = Path(__file__).resolve().parents[1] / "cosmolkit.pyi"


def invoke(value: object, *args: object, **kwargs: object) -> object:
    assert callable(value)
    return value(*args, **kwargs)


def class_nodes() -> dict[str, ast.ClassDef]:
    return {n.name: n for n in ast.parse(STUB.read_text()).body if isinstance(n, ast.ClassDef)}


def annotation(node: ast.expr | None) -> str:
    assert node is not None
    return ast.unparse(node)


class TestElementMetadata:
    @pytest.mark.parametrize("row", ROWS, ids=[f"e{r[2]:03d}-{r[0]}" for r in ROWS])
    def test_all_eight_fields(self, row: Row):
        name, symbol, number, period, outer, valences, rb0, weight, rb0_bits, weight_bits = row
        constant = cast(ck.Element, getattr(ck.Element, name))
        from_number = ck.Element.from_atomic_number(number)
        from_symbol = ck.Element.from_symbol(symbol)
        assert from_number is not None and from_symbol is not None
        assert type(constant) is ck.Element
        assert constant == from_number == from_symbol
        assert not constant != from_number
        assert hash(constant) == hash(from_number) == hash(from_symbol)
        assert len({constant, from_number, from_symbol}) == 1
        assert str(constant) == symbol
        assert repr(constant) == f"Element(symbol='{symbol}', atomic_number={number})"
        assert ck.Element.from_atomic_number(atomic_number=number) == constant
        assert ck.Element.from_symbol(symbol=symbol) == constant
        info = ck.element_info(element=constant)
        assert type(info) is ck.ElementInfo
        assert repr(info) == f"ElementInfo(symbol='{symbol}', atomic_number={number}, period={period})"
        for _ in range(3):
            assert constant.atomic_number() == number and type(constant.atomic_number()) is int
            assert constant.symbol() == symbol and type(constant.symbol()) is str
            assert info.element() == constant and type(info.element()) is ck.Element
            assert info.symbol() == symbol and type(info.symbol()) is str
            assert info.atomic_number() == number and type(info.atomic_number()) is int
            assert info.period() == period and type(info.period()) is int
            assert info.outer_electrons() == outer and type(info.outer_electrons()) is int
            actual = info.valences()
            assert type(actual) is list and actual == valences
            assert all(type(v) is int for v in actual)
            assert info.rb0() == rb0 and type(info.rb0()) is float
            assert info.atomic_weight() == weight and type(info.atomic_weight()) is float
            assert struct.pack(">d", info.rb0()).hex() == rb0_bits
            assert struct.pack(">d", info.atomic_weight()).hex() == weight_bits
        first = info.valences()
        second = info.valences()
        assert first is not second
        first.append(2**40)
        assert second == valences and info.valences() == valences
        assert ck.element_info(constant).valences() == valences


class TestElementValue:
    @pytest.mark.parametrize("number", range(119, 256), ids=[f"n{n}" for n in range(119, 256)])
    def test_in_u8_outside_element_domain(self, number: int):
        assert ck.Element.from_atomic_number(number) is None

    @pytest.mark.parametrize("number", [-1, 256, 2**64, -(2**64)], ids=["negative", "256", "huge_positive", "huge_negative"])
    def test_number_overflow(self, number: int):
        with pytest.raises(OverflowError):
            _ = ck.Element.from_atomic_number(number)

    @pytest.mark.parametrize("value", [6.0, "6", None, [], {}], ids=["float", "string", "none", "list", "dict"])
    def test_number_type_error(self, value: object):
        with pytest.raises(TypeError):
            _ = invoke(ck.Element.from_atomic_number, value)

    @pytest.mark.parametrize("value, expected", [(False, ck.Element.DUMMY), (True, ck.Element.H)], ids=["false", "true"])
    def test_bool_number(self, value: bool, expected: ck.Element):
        assert ck.Element.from_atomic_number(value) == expected

    @pytest.mark.parametrize("symbol, expected", [("Uut", ck.Element.NH), ("Uup", ck.Element.MC)], ids=["Uut", "Uup"])
    def test_source_symbol_alias(self, symbol: str, expected: ck.Element):
        assert ck.Element.from_symbol(symbol) == expected

    @pytest.mark.parametrize("symbol", ["", "c", "cl", "C ", " C", "6", "Uuo", "unknown"], ids=["empty", "lower_carbon", "lower_chlorine", "trailing_space", "leading_space", "numeric", "Uuo", "unknown"])
    def test_absent_symbols(self, symbol: str):
        assert ck.Element.from_symbol(symbol) is None

    @pytest.mark.parametrize("value", [6, None, b"C", []], ids=["int", "none", "bytes", "list"])
    def test_symbol_type_error(self, value: object):
        with pytest.raises(TypeError):
            _ = invoke(ck.Element.from_symbol, value)

    @pytest.mark.parametrize("value", [0, 6, True, "C", None, ck.element_info(ck.Element.C), object()], ids=["zero", "carbon_number", "bool", "symbol", "none", "metadata_result", "object"])
    def test_module_requires_element(self, value: object):
        with pytest.raises(TypeError):
            _ = invoke(ck.element_info, value)

    @pytest.mark.parametrize("cls", [ck.Element, ck.ElementInfo], ids=["Element", "ElementInfo"])
    def test_no_public_constructor(self, cls: object):
        with pytest.raises(TypeError):
            _ = invoke(cls)

    @pytest.mark.parametrize("name, error", [(n, k) for n in ["from_atomic_number", "from_symbol"] for k in ["missing", "extra", "wrong_keyword"]], ids=[f"{n}-{k}" for n in ["from_atomic_number", "from_symbol"] for k in ["missing", "extra", "wrong_keyword"]])
    def test_static_argument_error(self, name: str, error: str):
        method = cast(object, getattr(ck.Element, name))
        value: object = 6 if name == "from_atomic_number" else "C"
        with pytest.raises(TypeError):
            if error == "missing":
                _ = invoke(method)
            elif error == "extra":
                _ = invoke(method, value, value)
            else:
                _ = invoke(method, wrong=value)

    @pytest.mark.parametrize("error", ["missing", "extra", "wrong_keyword"])
    def test_module_argument_error(self, error: str):
        with pytest.raises(TypeError):
            if error == "missing":
                _ = invoke(ck.element_info)
            elif error == "extra":
                _ = invoke(ck.element_info, ck.Element.C, ck.Element.C)
            else:
                _ = invoke(ck.element_info, wrong=ck.Element.C)

    @pytest.mark.parametrize("cls, name, kind", [(c, n, k) for c, ns in [("Element", ELEMENT_READS), ("ElementInfo", INFO_READS)] for n in ns for k in ["positional", "keyword"]], ids=[f"{c}-{n}-{k}" for c, ns in [("Element", ELEMENT_READS), ("ElementInfo", INFO_READS)] for n in ns for k in ["positional", "keyword"]])
    def test_read_argument_error(self, cls: str, name: str, kind: str):
        receiver: object = ck.Element.C if cls == "Element" else ck.element_info(ck.Element.C)
        method = cast(object, getattr(receiver, name))
        with pytest.raises(TypeError):
            if kind == "positional":
                _ = invoke(method, 1)
            else:
                _ = invoke(method, extra=1)

    @pytest.mark.parametrize("label", [*['info_' + n for n in INFO_READS], "info_extra", "element_extra"])
    def test_frozen(self, label: str):
        value: object = ck.Element.C if label == "element_extra" else ck.element_info(ck.Element.C)
        name = label.removeprefix("element_").removeprefix("info_")
        with pytest.raises(AttributeError):
            setattr(value, name, 1)
        if isinstance(value, ck.ElementInfo):
            assert value.atomic_number() == 6 and value.valences() == [4]
        else:
            assert ck.Element.C.atomic_number() == 6

    @pytest.mark.parametrize("value", [6, "C", None, object()], ids=["int", "string", "none", "object"])
    def test_unrelated_equality(self, value: object):
        assert ck.Element.C.__eq__(value) is NotImplemented
        assert ck.Element.C.__ne__(value) is NotImplemented
        assert (ck.Element.C == value) is False
        assert (ck.Element.C != value) is True

    def test_canonical_namespace_only(self):
        for name in ["ELEMENT_MAP", "element_from_symbol", "get_element_info"]:
            assert not hasattr(ck, name)
        for name in ["UUT", "UUP", "UUO"]:
            assert not hasattr(ck.Element, name)
        assert not issubclass(ck.Element, int)
        assert ck.Element.C != ck.Element.N
        assert ck.Element.C.__eq__(ck.Element.from_symbol("C")) is True
        assert ck.Element.C.__ne__(ck.Element.from_symbol("C")) is False

    def test_exact_119_constants(self):
        expected = {r[0] for r in ROWS}
        native = {n for n in vars(ck.Element) if n.isupper()}
        assert len(native) == 119 and native == expected
        fields = {n.target.id: annotation(n.annotation) for n in class_nodes()["Element"].body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)}
        assert set(fields) == expected
        assert set(fields.values()) == {"typing.ClassVar[Element]"}
        for name, _, number, *_ in ROWS:
            assert cast(ck.Element, getattr(ck.Element, name)).atomic_number() == number

    @pytest.mark.parametrize("name", ["Element", "ElementInfo"])
    def test_native_stub_methods_and_types(self, name: str):
        cls = class_nodes()[name]
        methods = {n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)}
        expected = ({"from_atomic_number": "typing.Optional[Element]", "from_symbol": "typing.Optional[Element]", "atomic_number": "builtins.int", "symbol": "builtins.str", "__str__": "builtins.str", "__repr__": "builtins.str", "__hash__": "builtins.int", "__eq__": "builtins.bool | types.NotImplementedType", "__ne__": "builtins.bool | types.NotImplementedType"}
                    if name == "Element" else {"element": "Element", "symbol": "builtins.str", "atomic_number": "builtins.int", "period": "builtins.int", "outer_electrons": "builtins.int", "valences": "builtins.list[builtins.int]", "rb0": "builtins.float", "atomic_weight": "builtins.float", "__repr__": "builtins.str"})
        assert set(methods) == set(expected)
        assert {n: annotation(m.returns) for n, m in methods.items()} == expected
        runtime = ck.Element if name == "Element" else ck.ElementInfo
        for method in expected:
            assert callable(cast(object, getattr(runtime, method)))
            assert "property" not in [ast.unparse(d) for d in methods[method].decorator_list]
        assert [ast.unparse(d) for d in cls.decorator_list] == ["typing.final"]
        for method in (ELEMENT_READS if name == "Element" else INFO_READS):
            assert [a.arg for a in methods[method].args.args] == ["self"]
            assert not methods[method].args.defaults

    def test_required_signatures(self):
        for method, names in [(ck.Element.from_atomic_number, ["atomic_number"]), (ck.Element.from_symbol, ["symbol"]), (ck.element_info, ["element"])]:
            sig = inspect.signature(method)
            assert list(sig.parameters) == names
            assert all(cast(object, p.default) is inspect.Parameter.empty for p in sig.parameters.values())
        nodes = class_nodes()["Element"]
        methods = {n.name: n for n in nodes.body if isinstance(n, ast.FunctionDef)}
        for method, param, kind in [("from_atomic_number", "atomic_number", "builtins.int"), ("from_symbol", "symbol", "builtins.str")]:
            m = methods[method]
            assert [a.arg for a in m.args.args] == [param]
            assert annotation(m.args.args[0].annotation) == kind
            assert not m.args.defaults and not m.args.kw_defaults
            assert [ast.unparse(d) for d in m.decorator_list] == ["staticmethod"]
        for method in ["__eq__", "__ne__"]:
            assert [a.arg for a in methods[method].args.posonlyargs] == ["self", "value"]
            assert annotation(methods[method].args.posonlyargs[1].annotation) == "builtins.object"
        fn = next(n for n in ast.parse(STUB.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == "element_info")
        assert [a.arg for a in fn.args.args] == ["element"]
        assert annotation(fn.args.args[0].annotation) == "Element" and annotation(fn.returns) == "ElementInfo"
        assert not fn.args.defaults

    @pytest.mark.parametrize("number", [6, 119, 256, -1], ids=["carbon", "invalid119", "overflow256", "negative"])
    def test_index_conversion(self, number: int):
        class Index:
            def __index__(self) -> int:
                return number
        if number < 0 or number > 255:
            with pytest.raises(OverflowError):
                _ = invoke(ck.Element.from_atomic_number, Index())
        else:
            assert invoke(ck.Element.from_atomic_number, Index()) == (ck.Element.C if number == 6 else None)

    @pytest.mark.parametrize("case", ["RuntimeError", "non_integer_return"])
    def test_index_exception(self, case: str):
        class BrokenIndex:
            def __index__(self) -> object:
                if case == "RuntimeError":
                    raise RuntimeError("index sentinel")
                return "6"
        with pytest.raises(RuntimeError if case == "RuntimeError" else TypeError) as caught:
            _ = invoke(ck.Element.from_atomic_number, BrokenIndex())
        if case == "RuntimeError":
            assert str(caught.value) == "index sentinel"

    @pytest.mark.parametrize("case", ["int_subclass", "str_subclass"])
    def test_native_input_subclasses(self, case: str):
        class Number(int):
            pass
        class Symbol(str):
            pass
        if case == "int_subclass":
            assert ck.Element.from_atomic_number(Number(6)) == ck.Element.C
        else:
            assert ck.Element.from_symbol(Symbol("C")) == ck.Element.C
