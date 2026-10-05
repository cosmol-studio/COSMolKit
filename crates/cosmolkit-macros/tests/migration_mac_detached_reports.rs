#![allow(dead_code)]

#[path = "../src/status.rs"]
mod status;
extern crate proc_macro2 as proc_macro;
#[path = "../src/declaration.rs"]
mod declaration;
#[path = "../src/matrices.rs"]
mod matrices;
#[path = "../src/projection.rs"]
mod projection;
#[path = "../src/wrappers.rs"]
mod wrappers;

fn operation(name: &str, extra: &str) -> String {
    format!(
        r#"
        op {name}(amount: usize, mode: u8) {{
            method: {name}, impl_fn: crate::operations::body,
            kind: weak, access: {{ read: [], write: [] }},
            derived_effects: {{ recompute: [], preserve: [], invalidate: [], operation_defined: [] }},
            cip_state: preserve, feature: crate::FEATURE,
            parity: not_applicable, invariant_profile: "detached-report",
            {extra}
        }}
    "#
    )
}

#[test]
fn detached_reports_reject_ambiguous_or_pending_shapes() {
    for (fields, expected) in [
        (
            "report_result_type: crate::CanonicalResult,",
            "requires report_type",
        ),
        (
            "report_result_type: crate::CanonicalResult, inplace_result_type: usize, inplace: true,",
            "requires report_type",
        ),
        (
            "report_type: usize, report_result_type: crate::CanonicalResult, assemble_fn: crate::assemble,",
            "cannot be combined",
        ),
        (
            "report_type: usize, inplace_result_type: usize, inplace: true,",
            "mutually exclusive",
        ),
        (
            "report_type: usize, result_type: crate::PendingReport,",
            "cannot be combined",
        ),
        (
            "inplace_result_type: usize, assemble_fn: crate::assemble, inplace: true,",
            "cannot be combined",
        ),
        (
            "report_type: usize, output: multiple,",
            "require single-output",
        ),
        ("inplace_result_type: usize,", "requires inplace: true"),
        (
            "result_type: crate::PendingReport, inplace: true,",
            "pending result_type cannot generate an in-place wrapper",
        ),
    ] {
        let error = syn::parse_str::<declaration::MoleculeRegistry>(&operation("invalid", fields))
            .err()
            .expect("invalid report declaration accepted")
            .to_string();
        assert!(error.contains(expected), "{fields}: {error}");
    }
}

#[test]
fn generated_reports_execute_once_and_survive_only_successful_commit() {
    let defaults = "inplace: true, default_method: report_default, default_args: [0],";
    let source = operation("report", &format!("report_type: crate::Report, {defaults}"))
        + &operation(
            "position",
            "inplace_result_type: crate::Report, inplace: true, default_method: position_default, default_args: [0],",
        );
    let source = source
        + &operation(
            "canonical",
            "report_type: crate::Report, report_result_type: crate::CanonicalResult, inplace: true, default_method: canonical_default, default_args: [0],",
        );
    let registry = syn::parse_str::<declaration::MoleculeRegistry>(&source).unwrap();
    let generated = wrappers::expand_molecule_wrappers(&registry).unwrap();
    let fixture = r#"
        use std::{cell::RefCell, rc::Rc};
        type FunctionStatus = ();
        struct Spec { status: FunctionStatus }
        const REPORT_SPEC: Spec = Spec { status: () };
        const POSITION_SPEC: Spec = Spec { status: () };
        const CANONICAL_SPEC: Spec = Spec { status: () };
        mod ops { #[derive(Debug, PartialEq)] pub enum OperationError { Body, Finish, New } }
        use ops::OperationError;
        #[derive(Clone)] struct Molecule { value: usize, events: Rc<RefCell<Vec<&'static str>>> }
        struct Report { value: usize, events: Rc<RefCell<Vec<&'static str>>> }
        impl Drop for Report { fn drop(&mut self) { self.events.borrow_mut().push("drop-report"); } }
        struct CanonicalResult { molecule: Molecule, report: Report }
        impl From<(Molecule, Report)> for CanonicalResult {
            fn from((molecule, report): (Molecule, Report)) -> Self {
                molecule.events.borrow_mut().push("convert");
                Self { molecule, report }
            }
        }
        struct OpParts<'a> { target: Option<&'a mut Molecule>, value: usize, events: Rc<RefCell<Vec<&'static str>>>, fail: bool }
        impl<'a> OpParts<'a> {
            fn new(mol: &Molecule, _: &Spec) -> Result<Self, OperationError> {
                mol.events.borrow_mut().push("new");
                if mol.value == 999 { return Err(OperationError::New); }
                Ok(Self { target: None, value: mol.value, events: mol.events.clone(), fail: false })
            }
            fn new_in_place(mol: &'a mut Molecule, _: &Spec) -> Result<Self, OperationError> {
                mol.events.borrow_mut().push("new-inplace");
                if mol.value == 999 { return Err(OperationError::New); }
                Ok(Self { value: mol.value, events: mol.events.clone(), target: Some(mol), fail: false })
            }
            fn finish(self) -> Result<Molecule, OperationError> {
                self.events.borrow_mut().push("finish");
                if self.fail { return Err(OperationError::Finish); }
                Ok(Molecule { value: self.value, events: self.events })
            }
            fn finish_in_place(self) -> Result<(), OperationError> {
                self.events.borrow_mut().push("finish-inplace");
                if self.fail { return Err(OperationError::Finish); }
                self.target.unwrap().value = self.value;
                Ok(())
            }
            fn abort_in_place(self) { self.events.borrow_mut().push("abort"); }
        }
        mod operations {
            use super::*;
            pub fn body(parts: &mut OpParts<'_>, amount: usize, mode: u8) -> Result<Report, OperationError> {
                parts.events.borrow_mut().push("body");
                parts.value += amount;
                if mode == 1 { return Err(OperationError::Body); }
                parts.fail = mode == 2;
                Ok(Report { value: parts.value, events: parts.events.clone() })
            }
        }
        fn molecule() -> Molecule { Molecule { value: 10, events: Rc::default() } }
        fn main() {
            // Both short methods forward leading arguments and default mode.
            let mut mol = molecule();
            let (out, report) = mol.report_default(3).unwrap();
            assert_eq!((mol.value, out.value, report.value), (10, 13, 13));
            assert_eq!(&*mol.events.borrow(), &["new", "body", "finish"]);
            drop(report);
            mol.events.borrow_mut().clear();
            let report = mol.report_default_(4).unwrap();
            assert_eq!((mol.value, report.value), (14, 14));
            assert_eq!(&*mol.events.borrow(), &["new-inplace", "body", "finish-inplace"]);
            drop(report);
            let mut mol = molecule();
            let out: Molecule = mol.position_default(5).unwrap();
            assert_eq!((mol.value, out.value), (10, 15));
            assert_eq!(&*mol.events.borrow(), &["new", "body", "finish", "drop-report"]);
            mol.events.borrow_mut().clear();
            let report = mol.position_default_(7).unwrap();
            assert_eq!((mol.value, report.value), (17, 17));
            assert_eq!(&*mol.events.borrow(), &["new-inplace", "body", "finish-inplace"]);
            drop(report);
            for mode in [1, 2] {
                for inplace in [false, true] {
                    for scalar in [false, true] {
                        let mut mol = molecule();
                        let error = if inplace {
                            if scalar { mol.position_(3, mode).err() } else { mol.report_(3, mode).err() }
                        } else if scalar { mol.position(3, mode).err() } else { mol.report(3, mode).err() };
                        assert_eq!(error, Some(if mode == 1 { OperationError::Body } else { OperationError::Finish }));
                        assert_eq!(mol.value, 10);
                        let mut expected = vec![if inplace { "new-inplace" } else { "new" }, "body"];
                        if mode == 1 && inplace { expected.push("abort"); }
                        if mode == 2 { expected.extend([if inplace { "finish-inplace" } else { "finish" }, "drop-report"]); }
                        assert_eq!(*mol.events.borrow(), expected);
                    }
                }
            }
            let mut mol = molecule();
            let canonical: CanonicalResult = mol.canonical_default(6).unwrap();
            assert_eq!((mol.value, canonical.molecule.value, canonical.report.value), (10, 16, 16));
            assert_eq!(&*mol.events.borrow(), &["new", "body", "finish", "convert"]);
            drop(canonical);
            mol.events.borrow_mut().clear();
            let canonical: CanonicalResult = mol.canonical_default_(7).unwrap();
            assert_eq!((mol.value, canonical.molecule.value, canonical.report.value), (17, 17, 17));
            assert!(Rc::ptr_eq(&mol.events, &canonical.molecule.events));
            assert_eq!(&*mol.events.borrow(), &["new-inplace", "body", "finish-inplace", "convert"]);
            drop(canonical);
            for mode in [1, 2] {
                for inplace in [false, true] {
                    let mut mol = molecule();
                    let error = if inplace { mol.canonical_(3, mode).err() } else { mol.canonical(3, mode).err() };
                    assert!(error == Some(if mode == 1 { OperationError::Body } else { OperationError::Finish }));
                    assert_eq!(mol.value, 10);
                    assert!(!mol.events.borrow().contains(&"convert"));
                }
            }
            for inplace in [false, true] {
                let mut mol = molecule(); mol.value = 999;
                let error = if inplace { mol.report_(1, 0).err() } else { mol.report(1, 0).err() };
                assert_eq!(error, Some(OperationError::New));
                assert_eq!(mol.value, 999);
                assert_eq!(mol.events.borrow().len(), 1);
            }
        }
    "#;
    let temp =
        std::env::temp_dir().join(format!("cosmolkit-detached-reports-{}", std::process::id()));
    std::fs::create_dir_all(&temp).unwrap();
    let file = temp.join("reports.rs");
    std::fs::write(&file, format!("{fixture}\n{generated}")).unwrap();
    let binary = temp.join("reports");
    let compiled = std::process::Command::new("rustc")
        .args(["--edition=2024", "-o"])
        .arg(&binary)
        .arg(&file)
        .output()
        .unwrap();
    assert!(
        compiled.status.success(),
        "{}",
        String::from_utf8_lossy(&compiled.stderr)
    );
    let executed = std::process::Command::new(&binary).output().unwrap();
    assert!(
        executed.status.success(),
        "{}",
        String::from_utf8_lossy(&executed.stderr)
    );
    std::fs::remove_dir_all(temp).unwrap();
}

#[test]
fn detached_value_shapes_are_declared_in_the_generated_operation_matrix() {
    for (fields, expected) in [
        (
            "report_type: crate::Report,",
            "stringify!((Molecule,crate::Report))",
        ),
        (
            "report_type: crate::Report, report_result_type: crate::CanonicalResult, inplace: true,",
            "stringify!(crate::CanonicalResult)",
        ),
        ("inplace_result_type: usize, inplace: true,", "\"Molecule\""),
    ] {
        let registry =
            syn::parse_str::<declaration::MoleculeRegistry>(&operation("declared", fields))
                .unwrap();
        let output: String = matrices::expand_molecule_matrices(&registry)
            .unwrap()
            .to_string()
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect();
        assert!(output.contains(expected), "{output}");
    }
}
