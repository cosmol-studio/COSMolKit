//! Frozen field-wise ResidueInfo serialization products (BIO-IDENTITY1-36).
//!
//! Four source-literal rows x two repetitions = eight actual
//! serializations with exact field names/values and exactly 7 keys each.
//! The f32 weight JSON cast reference is frozen from the SOURCE literal,
//! not from actual serialized output. A custom recording Serializer proves
//! serialize_struct name/len, field call ORDER, and first-field error
//! propagation. Fresh whole-value and weight-bit checkpoints are captured
//! immediately before/after EVERY call, including error calls. No generic
//! NaN/float claim is made from these four finite fixtures; no production
//! helper exists for test convenience.

use cosmolkit_bio::{ResidueCode, ResidueInfoKind};
use serde::Serialize;
use serde::ser::SerializeStruct as _;
use serde::ser::Serializer as _;

struct FrozenRow {
    label: &'static str,
    idx: usize,
    code: ResidueCode,
    name: &'static str,
    kind: ResidueInfoKind,
    linking: u8,
    letter: char,
    hydrogens: u8,
    weight: f32,
    weight_json: serde_json::Value,
}

fn frozen_rows() -> [FrozenRow; 4] {
    [
        FrozenRow {
            label: "MSE17",
            idx: 17,
            code: ResidueCode::MSE,
            name: "MSE",
            kind: ResidueInfoKind::Aa,
            linking: 1,
            letter: 'm',
            hydrogens: 11,
            weight: 196.106f32,
            weight_json: serde_json::to_value(196.106f32).unwrap(),
        },
        FrozenRow {
            label: "HOH154",
            idx: 154,
            code: ResidueCode::HOH,
            name: "HOH",
            kind: ResidueInfoKind::Hoh,
            linking: 0,
            letter: ' ',
            hydrogens: 2,
            weight: 18.0153f32,
            weight_json: serde_json::to_value(18.0153f32).unwrap(),
        },
        FrozenRow {
            label: "UNK25",
            idx: 25,
            code: ResidueCode::UNK,
            name: "UNK",
            kind: ResidueInfoKind::Aa,
            linking: 1,
            letter: 'X',
            hydrogens: 9,
            weight: 103.120f32,
            weight_json: serde_json::to_value(103.120f32).unwrap(),
        },
        FrozenRow {
            label: "UNKNOWN367",
            idx: 367,
            code: ResidueCode::UNKNOWN,
            name: "",
            kind: ResidueInfoKind::Unknown,
            linking: 0,
            letter: ' ',
            hydrogens: 0,
            weight: 0.0f32,
            weight_json: serde_json::to_value(0.0f32).unwrap(),
        },
    ]
}

#[test]
fn info_serializer_eight_rows_exact_fields_and_order() {
    let rows = frozen_rows();
    let mut calls = 0usize;
    let mut discrepancies: Vec<String> = Vec::new();
    for _repetition in 0..2 {
        for row in &rows {
            let info = cosmolkit_bio::residue_info(row.idx);
            if info.code != row.code {
                discrepancies.push(format!("{}: code prerequisite", row.label));
                continue;
            }
            if info.name != row.name {
                discrepancies.push(format!("{}: name prerequisite", row.label));
                continue;
            }
            if info.kind != row.kind {
                discrepancies.push(format!("{}: kind prerequisite", row.label));
                continue;
            }
            if info.linking_type != row.linking {
                discrepancies.push(format!("{}: linking prerequisite", row.label));
                continue;
            }
            if info.one_letter_code != row.letter {
                discrepancies.push(format!("{}: letter prerequisite", row.label));
                continue;
            }
            if info.hydrogen_count != row.hydrogens {
                discrepancies.push(format!("{}: hydrogen prerequisite", row.label));
                continue;
            }
            if info.weight.to_bits() != row.weight.to_bits() {
                discrepancies.push(format!("{}: weight-bit prerequisite", row.label));
                continue;
            }
            let before = info;
            let before_weight_bits = before.weight.to_bits();
            let value = serde_json::to_value(&info);
            calls += 1;
            if info.weight.to_bits() != before_weight_bits {
                discrepancies.push(format!("{}: weight bits mutated", row.label));
            }
            if info != before {
                discrepancies.push(format!("{}: whole value mutated", row.label));
            }
            let value = value.expect("serialize ResidueInfo");
            let object = match value {
                serde_json::Value::Object(map) => map,
                other => {
                    discrepancies.push(format!("{}: not an object: {other:?}", row.label));
                    continue;
                }
            };
            if object.len() != 7 {
                discrepancies.push(format!("{}: {} keys != exactly 7", row.label, object.len()));
            }
            // serde_json's object map is a SORTED map, so its key order is
            // alphabetical and does not observe serialization call order;
            // the exact seven-field CALL order is proven by the recording
            // Serializer test below, not by this JSON key iteration.
            let code_string = serde_json::to_value(&row.code).unwrap();
            let kind_string = serde_json::to_value(&row.kind).unwrap();
            let checks = [
                ("code", object.get("code"), code_string),
                ("name", object.get("name"), serde_json::json!(row.name)),
                ("kind", object.get("kind"), kind_string),
                (
                    "linking_type",
                    object.get("linking_type"),
                    serde_json::json!(row.linking),
                ),
                (
                    "one_letter_code",
                    object.get("one_letter_code"),
                    serde_json::json!(row.letter),
                ),
                (
                    "hydrogen_count",
                    object.get("hydrogen_count"),
                    serde_json::json!(row.hydrogens),
                ),
                ("weight", object.get("weight"), row.weight_json.clone()),
            ];
            for (field, actual, expected) in checks {
                if actual != Some(&expected) {
                    discrepancies.push(format!(
                        "{}#{}: {actual:?} != source-frozen {expected:?}",
                        row.label, field
                    ));
                }
            }
        }
    }
    assert_eq!(calls, 8, "exact eight actual serializations");
    assert!(
        discrepancies.is_empty(),
        "collected discrepancies after all eight calls: {discrepancies:?}"
    );
}

thread_local! {
    static RECORD: std::cell::RefCell<Option<(Option<&'static str>, Option<usize>, Vec<&'static str>)>> =
        const { std::cell::RefCell::new(None) };
}

#[derive(Debug)]
struct RecError(String);

impl std::fmt::Display for RecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for RecError {}

impl serde::ser::Error for RecError {
    fn custom<T: std::fmt::Display>(msg: T) -> Self {
        RecError(msg.to_string())
    }
}

type Impossible = serde::ser::Impossible<(), RecError>;

struct RecSerializer {
    fail_at: Option<&'static str>,
}

struct RecStruct {
    fail_at: Option<&'static str>,
    fields: Vec<&'static str>,
    struct_name: Option<&'static str>,
    struct_len: Option<usize>,
}

impl serde::ser::SerializeStruct for RecStruct {
    type Ok = ();
    type Error = RecError;

    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), Self::Error> {
        if self.fail_at == Some(key) {
            self.fields.push(key);
            RECORD.with(|cell| {
                *cell.borrow_mut() = Some((
                    self.struct_name,
                    self.struct_len,
                    std::mem::take(&mut self.fields),
                ));
            });
            return Err(RecError(format!("injected at {key}")));
        }
        let _ = serde_json::to_value(value).map_err(|e| RecError(e.to_string()))?;
        self.fields.push(key);
        Ok(())
    }

    fn end(mut self) -> Result<(), Self::Error> {
        RECORD.with(|cell| {
            *cell.borrow_mut() = Some((
                self.struct_name,
                self.struct_len,
                std::mem::take(&mut self.fields),
            ));
        });
        Ok(())
    }
}

macro_rules! leaf {
    () => {
        Err(RecError("unsupported leaf".to_string()))
    };
}

impl serde::ser::Serializer for RecSerializer {
    type Ok = ();
    type Error = RecError;
    type SerializeSeq = Impossible;
    type SerializeTuple = Impossible;
    type SerializeTupleStruct = Impossible;
    type SerializeTupleVariant = Impossible;
    type SerializeMap = Impossible;
    type SerializeStruct = RecStruct;
    type SerializeStructVariant = Impossible;

    fn serialize_struct(
        self,
        name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeStruct, Self::Error> {
        Ok(RecStruct {
            fail_at: self.fail_at,
            fields: Vec::new(),
            struct_name: Some(name),
            struct_len: Some(len),
        })
    }

    fn serialize_bool(self, _v: bool) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_i8(self, _v: i8) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_i16(self, _v: i16) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_i32(self, _v: i32) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_i64(self, _v: i64) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_u8(self, _v: u8) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_u16(self, _v: u16) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_u32(self, _v: u32) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_u64(self, _v: u64) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_f32(self, _v: f32) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_f64(self, _v: f64) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_char(self, _v: char) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_str(self, _v: &str) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_bytes(self, _v: &[u8]) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_none(self) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_some<T: ?Sized + Serialize>(self, _v: &T) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_unit(self) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_unit_struct(self, _name: &'static str) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_unit_variant(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
    ) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_newtype_struct<T: ?Sized + Serialize>(
        self,
        _name: &'static str,
        _v: &T,
    ) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_newtype_variant<T: ?Sized + Serialize>(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _v: &T,
    ) -> Result<Self::Ok, Self::Error> {
        leaf!()
    }
    fn serialize_seq(self, _len: Option<usize>) -> Result<Self::SerializeSeq, Self::Error> {
        leaf!()
    }
    fn serialize_tuple(self, _len: usize) -> Result<Self::SerializeTuple, Self::Error> {
        leaf!()
    }
    fn serialize_tuple_struct(
        self,
        _name: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleStruct, Self::Error> {
        leaf!()
    }
    fn serialize_tuple_variant(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeTupleVariant, Self::Error> {
        leaf!()
    }
    fn serialize_map(self, _len: Option<usize>) -> Result<Self::SerializeMap, Self::Error> {
        leaf!()
    }
    fn serialize_struct_variant(
        self,
        _name: &'static str,
        _index: u32,
        _variant: &'static str,
        _len: usize,
    ) -> Result<Self::SerializeStructVariant, Self::Error> {
        leaf!()
    }
}

#[test]
fn info_serializer_records_struct_and_first_error_propagation() {
    // TWO original literal cases, each with its OWN fresh full-info +
    // weight-bit capture immediately BEFORE serialize, the Result stored,
    // and state compared immediately AFTER the call BEFORE expect/unwrap.
    let cases: [(&str, Option<&'static str>); 2] = [("success", None), ("error", Some("name"))];
    let mut invocations = 0usize;
    for (case, fail_at) in cases {
        let info = cosmolkit_bio::residue_info(17);
        assert_eq!(info.code, ResidueCode::MSE);
        assert_eq!(info.name, "MSE");
        let before = info;
        let before_bits = info.weight.to_bits();

        let serializer = RecSerializer { fail_at };
        let result = info.serialize(serializer);
        invocations += 1;
        // State comparison AFTER the actual call, BEFORE expect/unwrap.
        if info.weight.to_bits() != before_bits {
            panic!("{case}: serialize mutated weight bits");
        }
        if info != before {
            panic!("{case}: serialize mutated the whole value");
        }
        match fail_at {
            None => {
                result.expect("recorded serialization");
                let (struct_name, struct_len, fields) = RECORD
                    .with(|cell| cell.borrow_mut().take())
                    .expect("recording");
                assert_eq!(struct_name, Some("ResidueInfo"));
                assert_eq!(struct_len, Some(7));
                assert_eq!(
                    fields,
                    [
                        "code",
                        "name",
                        "kind",
                        "linking_type",
                        "one_letter_code",
                        "hydrogen_count",
                        "weight"
                    ],
                    "exact field call order"
                );
            }
            Some(_) => {
                let error = result.expect_err("injected first-field error");
                assert!(error.0.contains("injected at name"));
                let (_name, _len, error_fields) = RECORD
                    .with(|cell| cell.borrow_mut().take())
                    .expect("recording on error");
                assert_eq!(
                    error_fields,
                    ["code", "name"],
                    "code + name visited, name failed, no later fields"
                );
                assert_eq!(info.code, ResidueCode::MSE);
            }
        }
    }
    assert_eq!(invocations, 2, "count 2 after actual invocations");
}
