#![allow(dead_code)]

use std::fmt::Debug;

use crate::source_types::{SourceConstPointer, SourceHeap, SourceMutPointer};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SourcePortFixture<Input, Expected> {
    pub(crate) case_id: &'static str,
    pub(crate) input: Input,
    pub(crate) expected: Expected,
}

impl<Input, Expected> SourcePortFixture<Input, Expected> {
    pub(crate) const fn new(case_id: &'static str, input: Input, expected: Expected) -> Self {
        Self {
            case_id,
            input,
            expected,
        }
    }
}

pub(crate) fn assert_source_port_fixtures<Input, Expected>(
    fixtures: &[SourcePortFixture<Input, Expected>],
    mut operation: impl FnMut(Input) -> Expected,
) where
    Input: Clone,
    Expected: Debug + PartialEq,
{
    for fixture in fixtures {
        let actual = operation(fixture.input.clone());
        assert_eq!(
            actual, fixture.expected,
            "source-port fixture failed: {}",
            fixture.case_id
        );
    }
}

pub(crate) fn source_i8_bytes(bytes: &[u8]) -> Vec<i8> {
    bytes.iter().map(|byte| *byte as i8).collect()
}

pub(crate) fn source_i8_c_string(bytes_without_nul: &[u8]) -> Vec<i8> {
    assert!(
        !bytes_without_nul.contains(&0),
        "source C-string fixture contains an interior NUL"
    );
    let mut bytes = source_i8_bytes(bytes_without_nul);
    bytes.push(0);
    bytes
}

pub(crate) fn allocate_source_fixture<T: 'static>(heap: &mut SourceHeap, values: Vec<T>) -> SourceMutPointer<T> {
    heap.allocate(values)
        .expect("source fixture allocation identifier exhausted")
}

pub(crate) fn assert_source_slice_eq<T>(
    case_id: &str,
    heap: &SourceHeap,
    pointer: SourceConstPointer<T>,
    expected: &[T],
) where
    T: Debug + PartialEq + 'static,
{
    let actual = heap
        .slice(pointer)
        .unwrap_or_else(|error| panic!("source-port fixture {case_id} is unreadable: {error:?}"));
    assert_eq!(actual, expected, "source-port fixture failed: {case_id}");
}

pub(crate) fn assert_f32_bits_eq(case_id: &str, actual: f32, expected: f32) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "source-port fixture failed: {case_id}"
    );
}

pub(crate) fn assert_f64_bits_eq(case_id: &str, actual: f64, expected: f64) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "source-port fixture failed: {case_id}"
    );
}
