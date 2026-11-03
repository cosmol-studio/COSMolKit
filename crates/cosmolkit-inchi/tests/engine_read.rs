use cosmolkit_inchi::engine::{self, Atom, Input, ReadInput};

#[test]
fn reads_owned_official_records_without_a_toolkit() {
    // GetStructFromINCHI's methane output, also covered at the source boundary
    // by source_port__inchi_dll__getstructfrominchi__line_2856.
    let output = engine::read(ReadInput {
        inchi: c"InChI=1S/CH4/h1H4",
        options: None,
    })
    .unwrap();
    assert_eq!(output.return_code, 0);
    let mut methane = Atom::default();
    methane.elname[0] = b'C' as i8;
    methane.num_iso_H = [4, 0, 0, 0];
    assert_eq!(output.atoms, std::slice::from_ref(&methane));
    assert!(output.stereo0d.is_empty());
    assert_eq!(output.warning_flags, [[0; 2]; 2]);
    assert_eq!(output.message.as_deref(), Some(b"".as_slice()));
    assert_eq!(
        output.log.as_deref(),
        Some(
            b"Input format: InChI (plain identifier)\nOutput format: Plain text\nAux. info suppressed\nNo timeout\nUp to 1024 atoms per structure"
                .as_slice()
        )
    );
    // A second invocation cannot invalidate the first result's owned buffers.
    let _ = engine::read(ReadInput {
        inchi: c"InChI=1S/H2O/h1H2",
        options: Some(c""),
    })
    .unwrap();
    assert_eq!(output.atoms, std::slice::from_ref(&methane));
    assert_eq!(
        engine::generate(Input {
            atoms: &output.atoms,
            stereo0d: &output.stereo0d,
            options: None,
        })
        .unwrap()
        .inchi
        .as_deref(),
        Some(b"InChI=1S/CH4/h1H4".as_slice())
    );
}

#[test]
fn preserves_official_tetrahedral_neighbor_order_and_explicit_hydrogen() {
    // Exact raw upstream record from the GetStructFromINCHI source regression;
    // unlike a toolkit reader, this entry does not remove the seventh atom.
    let output = engine::read(ReadInput {
        inchi: c"InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1",
        options: None,
    })
    .unwrap();
    assert_eq!(output.return_code, 0);
    assert_eq!(output.atoms.len(), 7);
    assert_eq!(output.atoms[6].elname[0], b'H' as i8);
    assert_eq!(output.stereo0d.len(), 1);
    let stereo = &output.stereo0d[0];
    assert_eq!(stereo.neighbor, [6, 0, 2, 3]);
    assert_eq!(stereo.central_atom, 1);
    assert_eq!(stereo.type_, 2);
    assert_eq!(stereo.parity, 1);
}

#[test]
fn preserves_allene_stereo_for_direct_engine_regeneration() {
    let inchi = c"InChI=1S/C5H8/c1-3-5-4-2/h3-4H,1-2H3/t5-/m0/s1";
    let output = engine::read(ReadInput {
        inchi,
        options: None,
    })
    .unwrap();
    assert_eq!(output.return_code, 0);
    assert_eq!(output.stereo0d.len(), 1);
    assert_eq!(output.stereo0d[0].type_, 3); // INCHI_StereoType_Allene
    let regenerated = engine::generate(Input {
        atoms: &output.atoms,
        stereo0d: &output.stereo0d,
        options: None,
    })
    .unwrap();
    assert_eq!(regenerated.return_code, 0);
    assert_eq!(regenerated.inchi.as_deref(), Some(inchi.to_bytes()));
}

#[test]
fn invalid_inchi_retains_native_status_without_a_fallback_structure() {
    for inchi in [c"", c"bad", c"InChI=2S/CH4"] {
        let output = engine::read(ReadInput {
            inchi,
            options: None,
        })
        .unwrap();
        assert_eq!(output.return_code, -1); // inchi_Ret_EOF
        assert!(output.atoms.is_empty());
        assert!(output.stereo0d.is_empty());
        assert_eq!(output.message, None);
        assert_eq!(output.log, None);
        assert_eq!(output.warning_flags, [[0; 2]; 2]);
    }
}

#[test]
fn passes_read_options_to_the_official_engine() {
    let output = engine::read(ReadInput {
        inchi: c"InChI=1S/CH4/h1H4",
        options: Some(c"-DoNotAddH"),
    })
    .unwrap();
    assert_eq!(output.return_code, 0);
    assert!(
        output
            .log
            .as_deref()
            .unwrap()
            .windows(b"Do not add H".len())
            .any(|part| part == b"Do not add H")
    );
}
