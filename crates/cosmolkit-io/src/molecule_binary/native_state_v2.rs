//! Explicit COS NativeStateV2, authorized by ROOT CK-a88d982c9d33462896acbb759460e4e8.
//! Raw4 and canonical3 share this single complete codec. Two wire payloads and
//! their full equality comparison cost O(state) twice; no RDKit wire claim.
use super::*;
use cosmolkit_model::{PropertyText, SourceAtomValenceFacts};

fn bad(message: &str) -> PickleError {
    PickleError::InvalidArchive(message.into())
}
fn enumeration(value: u8, type_name: &'static str) -> PickleError {
    PickleError::InvalidEnumValue { value, type_name }
}
fn bytes(w: &mut PickleWriter, value: &[u8]) {
    if value.len() > 10_000_000 {
        w.error
            .get_or_insert(PickleError::StringTooLong(value.len()));
        return;
    }
    w.write_count(value.len());
    w.buf.extend_from_slice(value);
}
fn text(w: &mut PickleWriter, value: &PropertyText) {
    bytes(w, value.as_bytes());
}
fn read_text(r: &mut PickleReader<'_>) -> Result<PropertyText, PickleError> {
    let len = r.read_u32()? as usize;
    if len > 10_000_000 {
        return Err(PickleError::StringTooLong(len));
    }
    Ok(r.read_exact_slice(len)?.into())
}
fn count(w: &mut PickleWriter, len: usize) {
    if len > 1_000_000 {
        w.error.get_or_insert(bad("native count exceeds 1000000"));
    }
    w.write_count(len);
}
fn id(w: &mut PickleWriter, index: usize) {
    match u64::try_from(index) {
        Ok(value) => w.write_u64(value),
        Err(error) => {
            w.error.get_or_insert(invalid(error));
        }
    }
}
fn opt<T>(w: &mut PickleWriter, value: Option<T>, write: impl FnOnce(&mut PickleWriter, T)) {
    w.write_bool(value.is_some());
    if let Some(value) = value {
        write(w, value);
    }
}
fn read_opt<T>(
    r: &mut PickleReader<'_>,
    read: impl FnOnce(&mut PickleReader<'_>) -> Result<T, PickleError>,
) -> Result<Option<T>, PickleError> {
    if r.read_bool()? {
        Ok(Some(read(r)?))
    } else {
        Ok(None)
    }
}
fn vec<T>(w: &mut PickleWriter, values: &[T], mut write: impl FnMut(&mut PickleWriter, &T)) {
    count(w, values.len());
    for value in values {
        write(w, value);
    }
}
fn read_vec<T>(
    r: &mut PickleReader<'_>,
    min: usize,
    mut read: impl FnMut(&mut PickleReader<'_>) -> Result<T, PickleError>,
) -> Result<Vec<T>, PickleError> {
    let count = r.read_count(min)?;
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(read(r)?);
    }
    Ok(values)
}
fn value(w: &mut PickleWriter, value: &PropertyValue) {
    // Native numeric codes are the reviewed COS wire definition, not casts of
    // a Rust enum or source-name inference. f64 is written as exact bits.
    match value {
        PropertyValue::String(s) => {
            w.write_u8(0);
            text(w, s);
        }
        PropertyValue::Int(n) => {
            w.write_u8(1);
            w.write_i32(*n);
        }
        PropertyValue::UInt(n) => {
            w.write_u8(2);
            w.write_u32(*n);
        }
        PropertyValue::IntVector(v) => {
            w.write_u8(3);
            vec(w, v, |w, n| w.write_i32(*n));
        }
        PropertyValue::Double(n) => {
            w.write_u8(4);
            w.write_f64(*n);
        }
        PropertyValue::Bool(b) => {
            w.write_u8(5);
            w.write_bool(*b);
        }
        PropertyValue::StringVector(v) => {
            w.write_u8(6);
            vec(w, v, text);
        }
    }
}
fn read_native_value(r: &mut PickleReader<'_>) -> Result<PropertyValue, PickleError> {
    Ok(match r.read_u8()? {
        0 => PropertyValue::String(read_text(r)?),
        1 => PropertyValue::Int(r.read_i32()?),
        2 => PropertyValue::UInt(r.read_u32()?),
        3 => PropertyValue::IntVector(read_vec(r, 4, |r| r.read_i32())?),
        4 => PropertyValue::Double(r.read_f64()?),
        5 => PropertyValue::Bool(r.read_bool()?),
        6 => PropertyValue::StringVector(read_vec(r, 4, read_text)?),
        value => return Err(enumeration(value, "PropertyValue")),
    })
}
fn props<'a>(
    w: &mut PickleWriter,
    len: usize,
    rows: impl Iterator<Item = (&'a PropertyText, &'a PropertyValue)>,
) {
    count(w, len);
    let mut actual = 0;
    for (key, val) in rows {
        if key.is_empty() {
            w.error.get_or_insert(bad("empty native property key"));
        }
        text(w, key);
        value(w, val);
        actual += 1;
    }
    if actual != len {
        w.error
            .get_or_insert(bad("native property record count disagreement"));
    }
}
fn read_props(r: &mut PickleReader<'_>) -> Result<Vec<(PropertyText, PropertyValue)>, PickleError> {
    let rows = read_vec(r, 5, |r| Ok((read_text(r)?, read_native_value(r)?)))?;
    let mut seen = BTreeSet::new();
    for (key, _) in &rows {
        if key.is_empty() || !seen.insert(key) {
            return Err(bad("empty or duplicate native property key"));
        }
    }
    Ok(rows)
}
fn map(w: &mut PickleWriter, rows: &BTreeMap<PropertyText, PropertyText>) {
    count(w, rows.len());
    for (key, val) in rows {
        if key.is_empty() {
            w.error.get_or_insert(bad("empty native property key"));
        }
        text(w, key);
        text(w, val);
    }
}
fn read_map(r: &mut PickleReader<'_>) -> Result<BTreeMap<PropertyText, PropertyText>, PickleError> {
    let rows = read_vec(r, 8, |r| Ok((read_text(r)?, read_text(r)?)))?;
    let mut result = BTreeMap::new();
    let mut previous: Option<PropertyText> = None;
    for (key, value) in rows {
        if previous.as_ref().is_some_and(|p| p >= &key) {
            return Err(bad("unordered or duplicate native scalar map key"));
        }
        previous = Some(key.clone());
        result.insert(key, value);
    }
    Ok(result)
}
fn dimension(w: &mut PickleWriter, d: CoordinateDimension) {
    w.write_u8(match d {
        CoordinateDimension::TwoD => 0,
        CoordinateDimension::ThreeD => 1,
    });
}
fn read_dimension(r: &mut PickleReader<'_>) -> Result<CoordinateDimension, PickleError> {
    match r.read_u8()? {
        0 => Ok(CoordinateDimension::TwoD),
        1 => Ok(CoordinateDimension::ThreeD),
        value => Err(enumeration(value, "CoordinateDimension")),
    }
}
fn atom(w: &mut PickleWriter, a: &Atom) {
    id(w, a.id().index());
    w.write_u8(a.atomic_number());
    w.write_i8(a.formal_charge());
    opt(w, a.isotope(), |w, n| {
        w.buf.extend_from_slice(&n.to_le_bytes())
    });
    write_chiral_tag(w, a.chiral_tag());
    opt(w, a.chiral_permutation(), |w, n| w.write_u32(n));
    w.write_bool(a.unknown_stereo());
    opt(w, a.mol_parity(), |w, n| w.write_i32(n));
    opt(w, a.mol_inversion_flag(), |w, n| w.write_i32(n));
    w.write_u8(a.radical_electrons());
    w.write_bool(a.is_aromatic());
    write_hybridization(w, a.hybridization());
    opt(w, a.atom_map(), |w, n| w.write_u32(n));
    w.write_bool(a.no_implicit());
    w.write_bool(a.implicit_hydrogen());
    w.write_u8(a.explicit_hydrogens());
    vec(w, a.tracked_isotopic_hydrogens(), |w, n| {
        w.buf.extend_from_slice(&n.to_le_bytes())
    });
    let facts = a.source_valence_facts();
    w.write_i8(facts.explicit_valence);
    w.write_i8(facts.implicit_valence);
    w.write_u64(a.temporary_flags());
    props(w, a.props().len(), ordered_atom_properties(a));
    write_pdb_info(w, a.pdb_residue_info());
    opt(w, a.template_attachment_order(), |w, order| {
        vec(w, order.entries(), |w, e| {
            id(w, e.target().index());
            w.write_string(e.label());
        });
    });
}
fn read_atom(r: &mut PickleReader<'_>) -> Result<Atom, PickleError> {
    let aid = read_native_id(r)?;
    let number = r.read_u8()?;
    let element = Element::from_atomic_number(number).ok_or(enumeration(number, "Element"))?;
    let mut a = Atom::from_spec(AtomId::new(aid), AtomSpec::new(element));
    a.set_formal_charge(r.read_i8()?);
    a.set_isotope(read_opt(r, read_u16_le)?);
    a.set_chiral_tag(read_chiral_tag(r)?);
    a.set_chiral_permutation(read_opt(r, |r| r.read_u32())?);
    a.set_unknown_stereo(r.read_bool()?);
    a.set_mol_parity(read_opt(r, |r| r.read_i32())?);
    a.set_mol_inversion_flag(read_opt(r, |r| r.read_i32())?);
    a.set_radical_electrons(r.read_u8()?);
    a.set_aromatic(r.read_bool()?);
    a.set_hybridization(read_hybridization(r)?);
    a.set_atom_map(read_opt(r, |r| r.read_u32())?);
    a.set_no_implicit(r.read_bool()?);
    a.set_implicit_hydrogen(r.read_bool()?);
    a.set_explicit_hydrogens(r.read_u8()?);
    a.set_tracked_isotopic_hydrogens(read_vec(r, 2, read_u16_le)?);
    a.set_source_valence_facts(SourceAtomValenceFacts {
        explicit_valence: r.read_i8()?,
        implicit_valence: r.read_i8()?,
    });
    a.set_temporary_flags(r.read_u64()?);
    // Ordinary replay into a fresh detached value never reads __computedProps.
    for (key, value) in read_props(r)? {
        a.set_prop(key, value).map_err(invalid)?;
    }
    a.set_pdb_residue_info(read_pdb_info(r)?);
    let order = read_opt(r, |r| {
        let entries = read_vec(r, 12, |r| {
            Ok(TemplateAttachment::new(
                AtomId::new(read_native_id(r)?),
                r.read_string()?,
            ))
        })?;
        TemplateAttachmentOrder::new(entries).map_err(invalid)
    })?;
    replace_atom_template_attachment_order(&mut a, order);
    Ok(a)
}
fn bond(w: &mut PickleWriter, b: &Bond) {
    id(w, b.id().index());
    id(w, b.begin().index());
    id(w, b.end().index());
    write_bond_order(w, b.order());
    write_bond_stereo(w, b.stereo());
    write_bond_direction(w, b.direction());
    w.write_bool(b.is_aromatic());
    w.write_bool(b.is_conjugated());
    opt(w, b.stereo_atoms(), |w, ids| {
        id(w, ids[0].index());
        id(w, ids[1].index());
    });
    w.write_bool(b.unknown_stereo());
    w.write_u64(b.temporary_flags());
    props(w, b.props().len(), ordered_bond_properties(b));
}
fn read_native_bond(r: &mut PickleReader<'_>) -> Result<Bond, PickleError> {
    let bid = BondId::new(read_native_id(r)?);
    let begin = AtomId::new(read_native_id(r)?);
    let end = AtomId::new(read_native_id(r)?);
    let spec = BondSpec::new(begin, end, read_bond_order(r)?)
        .with_stereo(read_bond_stereo(r)?)
        .with_direction(read_bond_direction(r)?)
        .with_aromatic(r.read_bool()?)
        .with_conjugated(r.read_bool()?);
    let mut b = Bond::from_spec(bid, spec);
    if let Some(ids) = read_opt(r, |r| {
        Ok([
            AtomId::new(read_native_id(r)?),
            AtomId::new(read_native_id(r)?),
        ])
    })? {
        b.set_stereo_atoms(Some(ids));
    }
    b.set_unknown_stereo(r.read_bool()?);
    b.set_temporary_flags(r.read_u64()?);
    for (key, value) in read_props(r)? {
        b.set_prop(key, value).map_err(invalid)?;
    }
    Ok(b)
}
fn coordinates(w: &mut PickleWriter, c: &CoordinateBlock) {
    if let Some(order) = &c.source_conformer_order {
        let two = order
            .iter()
            .filter(|&&d| d == CoordinateDimension::TwoD)
            .count();
        if two != c.conformers_2d.len() || order.len() - two != c.conformers_3d.len() {
            w.error
                .get_or_insert(bad("native conformer occurrence order disagreement"));
        }
    }

    vec(w, &c.conformers_2d, |w, c| {
        id(w, c.id());
        vec(w, c.coordinates(), |w, row| {
            for x in row {
                w.write_f64(*x)
            }
        });
        map(w, c.props());
    });
    vec(w, &c.conformers_3d, |w, c| {
        id(w, c.id());
        vec(w, c.coordinates(), |w, row| {
            for x in row {
                w.write_f64(*x)
            }
        });
        w.write_bool(c.is_3d());
        map(w, c.props());
    });
    opt(w, c.source_coordinate_dim, dimension);
    opt(w, c.source_conformer_order.as_ref(), |w, v| {
        vec(w, v, |w, d| dimension(w, *d))
    });
}
fn read_coordinates(r: &mut PickleReader<'_>) -> Result<CoordinateBlock, PickleError> {
    let conformers_2d = read_vec(r, 16, |r| {
        let id = read_native_id(r)?;
        let rows = read_vec(r, 16, |r| Ok([r.read_f64()?, r.read_f64()?]))?;
        let mut c = Conformer2D::new(id, rows);
        for (k, v) in read_map(r)? {
            c = c.with_prop(k, v);
        }
        Ok(c)
    })?;
    let conformers_3d = read_vec(r, 17, |r| {
        let id = read_native_id(r)?;
        let rows = read_vec(r, 24, |r| Ok([r.read_f64()?, r.read_f64()?, r.read_f64()?]))?;
        let mut c = Conformer3D::new(id, rows, r.read_bool()?);
        for (k, v) in read_map(r)? {
            c = c.with_prop(k, v);
        }
        Ok(c)
    })?;
    let source_coordinate_dim = read_opt(r, read_dimension)?;
    let source_conformer_order = read_opt(r, |r| read_vec(r, 1, read_dimension))?;
    let c = CoordinateBlock {
        conformers_2d,
        conformers_3d,
        source_coordinate_dim,
        source_conformer_order,
    };
    if let Some(order) = &c.source_conformer_order {
        let two = order
            .iter()
            .filter(|&&d| d == CoordinateDimension::TwoD)
            .count();
        if two != c.conformers_2d.len() || order.len() - two != c.conformers_3d.len() {
            return Err(bad("native conformer occurrence order disagreement"));
        }
    }
    Ok(c)
}
fn atom_ids(w: &mut PickleWriter, ids: &[AtomId]) {
    vec(w, ids, |w, idv| id(w, idv.index()));
}
fn bond_ids(w: &mut PickleWriter, ids: &[BondId]) {
    vec(w, ids, |w, idv| id(w, idv.index()));
}
fn read_atom_ids(r: &mut PickleReader<'_>) -> Result<Vec<AtomId>, PickleError> {
    read_vec(r, 8, |r| Ok(AtomId::new(read_native_id(r)?)))
}
fn read_bond_ids(r: &mut PickleReader<'_>) -> Result<Vec<BondId>, PickleError> {
    read_vec(r, 8, |r| Ok(BondId::new(read_native_id(r)?)))
}
fn kind(w: &mut PickleWriter, k: &SubstanceGroupKind) {
    let tag = match k {
        SubstanceGroupKind::Data => 0,
        SubstanceGroupKind::Superatom => 1,
        SubstanceGroupKind::MultipleGroup => 2,
        SubstanceGroupKind::StructuralRepeatUnit => 3,
        SubstanceGroupKind::Monomer => 4,
        SubstanceGroupKind::Copolymer => 5,
        SubstanceGroupKind::Crosslink => 6,
        SubstanceGroupKind::Graft => 7,
        SubstanceGroupKind::Modification => 8,
        SubstanceGroupKind::Mer => 9,
        SubstanceGroupKind::AnyPolymer => 10,
        SubstanceGroupKind::MixtureComponent => 11,
        SubstanceGroupKind::Mixture => 12,
        SubstanceGroupKind::Formulation => 13,
        SubstanceGroupKind::Generic(_) => 14,
    };
    w.write_u8(tag);
    if let SubstanceGroupKind::Generic(v) = k {
        text(w, v);
    }
}
fn read_kind(r: &mut PickleReader<'_>) -> Result<SubstanceGroupKind, PickleError> {
    Ok(match r.read_u8()? {
        0 => SubstanceGroupKind::Data,
        1 => SubstanceGroupKind::Superatom,
        2 => SubstanceGroupKind::MultipleGroup,
        3 => SubstanceGroupKind::StructuralRepeatUnit,
        4 => SubstanceGroupKind::Monomer,
        5 => SubstanceGroupKind::Copolymer,
        6 => SubstanceGroupKind::Crosslink,
        7 => SubstanceGroupKind::Graft,
        8 => SubstanceGroupKind::Modification,
        9 => SubstanceGroupKind::Mer,
        10 => SubstanceGroupKind::AnyPolymer,
        11 => SubstanceGroupKind::MixtureComponent,
        12 => SubstanceGroupKind::Mixture,
        13 => SubstanceGroupKind::Formulation,
        14 => SubstanceGroupKind::Generic(read_text(r)?),
        value => return Err(enumeration(value, "SubstanceGroupKind")),
    })
}
fn connection(w: &mut PickleWriter, c: &SGroupConnection) {
    match c {
        SGroupConnection::HeadToHead => w.write_u8(0),
        SGroupConnection::HeadToTail => w.write_u8(1),
        SGroupConnection::Either => w.write_u8(2),
        SGroupConnection::Unknown(v) => {
            w.write_u8(3);
            text(w, v);
        }
    }
}
fn read_connection(r: &mut PickleReader<'_>) -> Result<SGroupConnection, PickleError> {
    Ok(match r.read_u8()? {
        0 => SGroupConnection::HeadToHead,
        1 => SGroupConnection::HeadToTail,
        2 => SGroupConnection::Either,
        3 => SGroupConnection::Unknown(read_text(r)?),
        value => return Err(enumeration(value, "SGroupConnection")),
    })
}
fn bracket_style(w: &mut PickleWriter, c: &SGroupBracketStyle) {
    match c {
        SGroupBracketStyle::Bracket => w.write_u8(0),
        SGroupBracketStyle::Parenthesis => w.write_u8(1),
        SGroupBracketStyle::None => w.write_u8(2),
        SGroupBracketStyle::Unknown(v) => {
            w.write_u8(3);
            text(w, v);
        }
    }
}
fn read_bracket_style(r: &mut PickleReader<'_>) -> Result<SGroupBracketStyle, PickleError> {
    Ok(match r.read_u8()? {
        0 => SGroupBracketStyle::Bracket,
        1 => SGroupBracketStyle::Parenthesis,
        2 => SGroupBracketStyle::None,
        3 => SGroupBracketStyle::Unknown(read_text(r)?),
        value => return Err(enumeration(value, "SGroupBracketStyle")),
    })
}
fn group(w: &mut PickleWriter, s: &SubstanceGroup) {
    id(w, s.id().index());
    opt(w, s.rdkit_sequence_id(), |w, v| w.write_u32(v));
    opt(w, s.external_id(), |w, v| w.write_u32(v));
    kind(w, s.kind());
    atom_ids(w, s.atoms());
    bond_ids(w, s.bonds());
    // The reviewed constructors keep sparse roles within members. Sorting the
    // distinct existing member IDs exposes the actual BTreeMap order, without
    // introducing restoration rights outside members. O(B log B) extra cost.
    let roles: Vec<_> = s
        .bonds()
        .iter()
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .filter_map(|b| s.explicit_bond_role(b).map(|role| (b, role)))
        .collect();
    vec(w, &roles, |w, (b, role)| {
        id(w, b.index());
        w.write_u8(match role {
            SGroupBondRole::Crossing => 0,
            SGroupBondRole::Contained => 1,
        });
    });
    bond_ids(w, s.head_crossing_bonds());
    bond_ids(w, s.crossing_bond_correspondence());
    atom_ids(w, s.parent_atoms());
    opt(w, s.parent(), |w, p| id(w, p.index()));
    opt(w, s.label(), text);
    opt(w, s.connection(), connection);
    opt(w, s.subtype(), text);
    opt(w, s.bracket_style(), bracket_style);
    opt(w, s.expansion_state(), text);
    opt(w, s.class(), text);
    opt(w, s.component_number(), |w, n| w.write_u32(n));
    opt(w, s.display(), |w, d| {
        vec(w, &d.brackets, |w, b| {
            for row in b.points {
                for v in row {
                    w.write_f64(v)
                }
            }
        });
        opt(w, d.field_position, |w, p| {
            w.write_f64(p[0]);
            w.write_f64(p[1]);
        });
        opt(w, d.display_tag.as_ref(), text);
    });
    opt(w, s.data(), |w, d| {
        for field in [
            &d.field_name,
            &d.field_type,
            &d.field_info,
            &d.field_display,
            &d.units,
            &d.query_type,
            &d.query_op,
        ] {
            opt(w, field.as_ref(), text);
        }
        vec(w, &d.values, text);
    });
    vec(w, s.attach_points(), |w, p| {
        id(w, p.atom.index());
        opt(w, p.leaving_atom, |w, a| id(w, a.index()));
        opt(w, p.label.as_ref(), text);
        opt(w, p.order, |w, n| w.write_u32(n));
    });
    vec(w, s.cstates(), |w, c| {
        id(w, c.bond.index());
        for v in c.vector {
            w.write_f64(v);
        }
    });
    props(w, s.props().len(), s.property_records());
    vec(w, s.data_fields(), text);
}
fn read_group(r: &mut PickleReader<'_>) -> Result<SubstanceGroup, PickleError> {
    let gid = SubstanceGroupId::new(read_native_id(r)?);
    let sequence = read_opt(r, |r| r.read_u32())?;
    let external = read_opt(r, |r| r.read_u32())?;
    let mut s = SubstanceGroup::new(gid, read_kind(r)?)
        .with_atoms(read_atom_ids(r)?)
        .with_bonds(read_bond_ids(r)?);
    if let Some(v) = sequence {
        s = s.with_rdkit_sequence_id(v);
    }
    if let Some(v) = external {
        s = s.with_external_id(v);
    }
    let roles = read_vec(r, 9, |r| {
        let b = BondId::new(read_native_id(r)?);
        let role = match r.read_u8()? {
            0 => SGroupBondRole::Crossing,
            1 => SGroupBondRole::Contained,
            value => return Err(enumeration(value, "SGroupBondRole")),
        };
        Ok((b, role))
    })?;
    let mut previous = None;
    for (b, role) in roles {
        if !s.bonds().contains(&b) || previous.is_some_and(|p| p >= b) {
            return Err(bad("invalid native sparse role record"));
        }
        previous = Some(b);
        s = s.with_bond_role(b, role);
    }
    s = s
        .with_head_crossing_bonds(read_bond_ids(r)?)
        .with_crossing_bond_correspondence(read_bond_ids(r)?)
        .with_parent_atoms(read_atom_ids(r)?);
    if let Some(v) = read_opt(r, read_native_id)? {
        s = s.with_parent(SubstanceGroupId::new(v));
    }
    if let Some(v) = read_opt(r, read_text)? {
        s = s.with_label(v);
    }
    if let Some(v) = read_opt(r, read_connection)? {
        s = s.with_connection(v);
    }
    if let Some(v) = read_opt(r, read_text)? {
        s = s.with_subtype(v);
    }
    if let Some(v) = read_opt(r, read_bracket_style)? {
        s = s.with_bracket_style(v);
    }
    if let Some(v) = read_opt(r, read_text)? {
        s = s.with_expansion_state(v);
    }
    if let Some(v) = read_opt(r, read_text)? {
        s = s.with_class(v);
    }
    if let Some(v) = read_opt(r, |r| r.read_u32())? {
        s = s.with_component_number(v);
    }
    if let Some(v) = read_opt(r, |r| {
        let brackets = read_vec(r, 72, |r| {
            let mut points = [[0.0; 3]; 3];
            for row in &mut points {
                for v in row {
                    *v = r.read_f64()?;
                }
            }
            Ok(SGroupBracket::new(points))
        })?;
        let field_position = read_opt(r, |r| Ok([r.read_f64()?, r.read_f64()?]))?;
        let display_tag = read_opt(r, read_text)?;
        Ok(SGroupDisplay {
            brackets,
            field_position,
            display_tag,
        })
    })? {
        s = s.with_display(v);
    }
    if let Some(v) = read_opt(r, |r| {
        Ok(SGroupData {
            field_name: read_opt(r, read_text)?,
            field_type: read_opt(r, read_text)?,
            field_info: read_opt(r, read_text)?,
            field_display: read_opt(r, read_text)?,
            units: read_opt(r, read_text)?,
            query_type: read_opt(r, read_text)?,
            query_op: read_opt(r, read_text)?,
            values: read_vec(r, 4, read_text)?,
        })
    })? {
        s = s.with_data(v);
    }
    s = s.with_attach_points(read_vec(r, 11, |r| {
        Ok(SGroupAttachPoint {
            atom: AtomId::new(read_native_id(r)?),
            leaving_atom: read_opt(r, |r| Ok(AtomId::new(read_native_id(r)?)))?,
            label: read_opt(r, read_text)?,
            order: read_opt(r, |r| r.read_u32())?,
        })
    })?);
    s = s.with_cstates(read_vec(r, 32, |r| {
        Ok(SGroupCState::new(
            BondId::new(read_native_id(r)?),
            [r.read_f64()?, r.read_f64()?, r.read_f64()?],
        ))
    })?);
    for (k, v) in read_props(r)? {
        s.set_prop(k, v).map_err(invalid)?;
    }
    for v in read_vec(r, 4, read_text)? {
        s.push_data_field(v);
    }
    Ok(s)
}
fn stereo(w: &mut PickleWriter, s: &StereoGroup) {
    opt(w, s.id(), |w, v| w.write_u32(v));
    w.write_u32(s.write_id());
    write_stereo_group_kind(w, s.kind());
    atom_ids(w, s.atoms());
    bond_ids(w, s.bonds());
}
fn read_stereo(r: &mut PickleReader<'_>) -> Result<StereoGroup, PickleError> {
    let id = read_opt(r, |r| r.read_u32())?;
    let write_id = r.read_u32()?;
    let kind = read_stereo_group_kind(r)?;
    let mut s =
        StereoGroup::new(kind, read_atom_ids(r)?, read_bond_ids(r)?)?.with_write_id(write_id);
    if let Some(id) = id {
        s = s.with_id(id);
    }
    Ok(s)
}
fn molecule_properties(w: &mut PickleWriter, p: &MoleculeProperties) {
    opt(w, p.name(), text);
    props(w, p.props().len(), p.ordered_props());
    vec(w, p.sdf_data_fields(), |w, (k, v)| {
        text(w, k);
        text(w, v);
    });
    vec(w, p.sdf_property_lists(), |w, p| {
        w.write_u8(match p.target() {
            SdfPropertyListTarget::Atom => 0,
            SdfPropertyListTarget::Bond => 1,
        });
        text(w, p.name());
        vec(w, p.values(), |w, v| opt(w, v.as_ref(), value));
    });
}
fn read_molecule_properties(r: &mut PickleReader<'_>) -> Result<MoleculeProperties, PickleError> {
    let mut p = MoleculeProperties::default();
    if let Some(n) = read_opt(r, read_text)? {
        p = p.with_name(n);
    }
    for (k, v) in read_props(r)? {
        p = p.with_prop(k, v).map_err(invalid)?;
    }
    for (k, v) in read_vec(r, 8, |r| Ok((read_text(r)?, read_text(r)?)))? {
        p = p.with_sdf_data_field(k, v);
    }
    for list in read_vec(r, 9, |r| {
        let target = match r.read_u8()? {
            0 => SdfPropertyListTarget::Atom,
            1 => SdfPropertyListTarget::Bond,
            value => return Err(enumeration(value, "SdfPropertyListTarget")),
        };
        let name = read_text(r)?;
        let values = read_vec(r, 1, |r| read_opt(r, read_native_value))?;
        Ok(SdfPropertyList::new(target, name, values))
    })? {
        p = p.with_sdf_property_list(list);
    }
    Ok(p)
}
fn derived(
    w: &mut PickleWriter,
    d: &BinaryDerivedView<'_>,
    bits: Option<u16>,
) -> Result<(), PickleError> {
    opt(w, bits, |w, b| w.buf.extend_from_slice(&b.to_le_bytes()));
    for rings in [d.rings, d.ring_families] {
        w.write_bool(rings.is_some());
        if let Some(r) = rings {
            count(w, r.atom_row_count());
            count(w, r.bond_row_count());
            write_ring_info(w, r)?;
        }
    }
    opt(w, d.valence, |w, v| {
        vec(w, &v.explicit_valence, |w, n| w.write_i32(*n));
        vec(w, &v.implicit_hydrogens, |w, n| w.write_i32(*n));
    });
    w.write_bool(d.aromaticity_valid);
    w.write_bool(d.stereo_valid);
    Ok(())
}
fn read_native_derived(
    r: &mut PickleReader<'_>,
    atoms: usize,
    bonds: usize,
) -> Result<BinaryDerivedState, PickleError> {
    let valid_bits = read_opt(r, read_u16_le)?;
    let mut slots = Vec::with_capacity(2);
    for _ in 0..2 {
        slots.push(if r.read_bool()? {
            let ac = r.read_u32()? as usize;
            let bc = r.read_u32()? as usize;
            if ac > atoms || bc > bonds {
                return Err(bad("native ring extent exceeds topology"));
            }
            Some(read_ring_info(r, ac, bc)?)
        } else {
            None
        });
    }
    let valence = read_opt(r, |r| {
        let explicit_valence = read_vec(r, 4, |r| r.read_i32())?;
        let implicit_hydrogens = read_vec(r, 4, |r| r.read_i32())?;
        if explicit_valence.len() != atoms || implicit_hydrogens.len() != atoms {
            return Err(bad("native valence row count disagreement"));
        }
        Ok(ValenceAssignment {
            explicit_valence,
            implicit_hydrogens,
        })
    })?;
    let mut slots = slots.into_iter();
    Ok(BinaryDerivedState {
        valid_bits,
        rings: slots.next().unwrap(),
        ring_families: slots.next().unwrap(),
        valence,
        aromaticity_valid: r.read_bool()?,
        stereo_valid: r.read_bool()?,
    })
}
fn encode_parts(
    topology: &TopologyBlock,
    coordinates_value: &CoordinateBlock,
    properties_value: &MoleculeProperties,
    derived_value: &BinaryDerivedView<'_>,
    bits: Option<u16>,
) -> Result<Vec<u8>, PickleError> {
    // The authorized complete codec covers concrete state only. Preserve the
    // existing independent query boundary instead of dropping an attached tree.
    if topology.bonds.iter().any(|bond| bond.query().is_some()) {
        return Err(PickleError::InvalidMolecule(
            "concrete native state contains a query bond".into(),
        ));
    }
    let mut w = PickleWriter::new();
    vec(&mut w, &topology.atoms, atom);
    vec(&mut w, &topology.bonds, bond);
    count(&mut w, topology.atoms.len());
    for i in 0..topology.atoms.len() {
        vec(&mut w, topology.adjacency.neighbors_of(i), |w, n| {
            id(w, n.atom_index);
            id(w, n.bond.index());
        });
    }
    coordinates(&mut w, coordinates_value);
    vec(&mut w, &topology.substance_groups, group);
    vec(&mut w, &topology.stereo_groups, stereo);
    molecule_properties(&mut w, properties_value);
    derived(&mut w, derived_value, bits)?;
    w.into_inner()
}
pub(super) fn encode_input(input: &BinaryInput<'_>) -> Result<Vec<u8>, PickleError> {
    encode_parts(
        input.topology,
        input.coordinates,
        input.properties,
        &input.derived,
        Some(input.derived.valid_bits),
    )
}
pub(super) fn encode_record(record: &BinaryRecord) -> Result<Vec<u8>, PickleError> {
    let d = &record.derived;
    let view = BinaryDerivedView {
        rings: d.rings.as_ref(),
        ring_families: d.ring_families.as_ref(),
        valence: d.valence.as_ref(),
        aromaticity_valid: d.aromaticity_valid,
        stereo_valid: d.stereo_valid,
        valid_bits: d.valid_bits.unwrap_or(0),
    };
    encode_parts(
        &record.topology,
        &record.coordinates,
        &record.properties,
        &view,
        d.valid_bits,
    )
}
pub(super) fn decode(data: &[u8]) -> Result<BinaryRecord, PickleError> {
    let mut r = PickleReader::new(data);
    let atoms = read_vec(&mut r, 1, read_atom)?;
    let bonds = read_vec(&mut r, 1, read_native_bond)?;
    let adjacency = read_vec(&mut r, 4, |r| {
        read_vec(r, 16, |r| {
            Ok((read_native_id(r)?, BondId::new(read_native_id(r)?)))
        })
    })?;
    let coordinates = read_coordinates(&mut r)?;
    let groups = read_vec(&mut r, 1, read_group)?;
    let stereo = read_vec(&mut r, 1, read_stereo)?;
    let properties = read_molecule_properties(&mut r)?;
    let derived = read_native_derived(&mut r, atoms.len(), bonds.len())?;
    if r.remaining() != 0 {
        return Err(bad("trailing NativeStateV2 bytes"));
    }
    let topology = TopologyBlock::try_from_parts(atoms, bonds, groups, stereo).map_err(invalid)?;
    if adjacency.len() != topology.atoms.len() {
        return Err(bad("native adjacency row count disagreement"));
    }
    for (i, row) in adjacency.iter().enumerate() {
        let actual = topology.adjacency.neighbors_of(i);
        if row.len() != actual.len()
            || row
                .iter()
                .zip(actual)
                .any(|(&(atom, bond), actual)| atom != actual.atom_index || bond != actual.bond)
        {
            return Err(bad("native adjacency disagrees with validated topology"));
        }
    }
    let result = BinaryRecord {
        topology,
        coordinates,
        properties,
        derived,
    };
    result.validate()?;
    Ok(result)
}

#[cfg(test)]
mod recovery_chem29_native_codec {
    use super::*;
    use cosmolkit_model::StereoGroupError;
    fn malformed(atoms: &[AtomId], bonds: &[BondId]) -> Vec<u8> {
        let mut w = PickleWriter::new();
        w.write_bool(true);
        w.write_u32(7);
        w.write_u32(9);
        write_stereo_group_kind(&mut w, StereoGroupKind::Or);
        atom_ids(&mut w, atoms);
        bond_ids(&mut w, bonds);
        w.into_inner().unwrap()
    }
    #[test]
    fn wire_duplicate_atoms_and_bonds_retain_typed_constructor_cause() {
        for (atoms, bonds, expected) in [
            (
                vec![AtomId::new(0), AtomId::new(0)],
                vec![BondId::new(0), BondId::new(0)],
                StereoGroupError::DuplicateAtom,
            ),
            (
                vec![AtomId::new(0)],
                vec![BondId::new(0), BondId::new(0)],
                StereoGroupError::DuplicateBond,
            ),
        ] {
            let bytes = malformed(&atoms, &bonds);
            let mut reader = PickleReader::new(&bytes);
            let error = read_stereo(&mut reader).unwrap_err();
            assert_eq!(error, PickleError::StereoGroup(expected));
            assert_eq!(
                std::error::Error::source(&error)
                    .unwrap()
                    .downcast_ref::<StereoGroupError>(),
                Some(&expected)
            );
        }
    }
    #[test]
    fn truncated_member_array_precedes_duplicate_constructor_and_valid_ids_roundtrip() {
        let mut bytes = malformed(&[AtomId::new(0), AtomId::new(0)], &[BondId::new(0)]);
        bytes.pop();
        let mut reader = PickleReader::new(&bytes);
        // The unchanged bounded-vector reader rejects the truncated bond
        // array before the duplicate atom constructor is reached.
        assert_eq!(
            read_stereo(&mut reader),
            Err(PickleError::InvalidArchive(
                "invalid bounded count: 1".into()
            ))
        );
        let original = StereoGroup::new(
            StereoGroupKind::And,
            vec![AtomId::new(1), AtomId::new(0)],
            vec![BondId::new(0)],
        )
        .unwrap()
        .with_id(7)
        .with_write_id(9);
        let mut w = PickleWriter::new();
        stereo(&mut w, &original);
        let bytes = w.into_inner().unwrap();
        assert_eq!(
            read_stereo(&mut PickleReader::new(&bytes)).unwrap(),
            original
        );
    }
}
