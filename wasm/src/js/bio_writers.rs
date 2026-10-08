//! Canonical BIO writer calls and complete mutable source-defined options.
use crate::{
    bio_readers::BioStructure,
    host_values::{bool_value, integer},
};
use cosmolkit_wasm::rust as ck;
use std::path::Path;
use wasm_bindgen::prelude::*;
fn u16_value(v: &JsValue, name: &str) -> Result<u16, JsValue> {
    Ok(integer(v, name, 0., 65535.)? as u16)
}
#[wasm_bindgen]
pub struct BioPdbWriteParams {
    pub(crate) inner: ck::BioPdbWriteParams,
}
#[wasm_bindgen]
impl BioPdbWriteParams {
    #[wasm_bindgen(constructor)]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ter_records: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] numbered_ter: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ter_ignores_type: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] preserve_serial: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] end_record: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::BioPdbWriteParams::default();
        if !ter_records.is_undefined() {
            inner.ter_records = bool_value(&ter_records, "terRecords")?;
        }
        if !numbered_ter.is_undefined() {
            inner.numbered_ter = bool_value(&numbered_ter, "numberedTer")?;
        }
        if !ter_ignores_type.is_undefined() {
            inner.ter_ignores_type = bool_value(&ter_ignores_type, "terIgnoresType")?;
        }
        if !preserve_serial.is_undefined() {
            inner.preserve_serial = bool_value(&preserve_serial, "preserveSerial")?;
        }
        if !end_record.is_undefined() {
            inner.end_record = bool_value(&end_record, "endRecord")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=terRecords)]
    pub fn ter_records(&self) -> bool {
        self.inner.ter_records
    }
    #[wasm_bindgen(setter,js_name=terRecords)]
    pub fn set_ter_records(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.ter_records = value;
        self.inner.ter_records = bool_value(&value, "terRecords")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=numberedTer)]
    pub fn numbered_ter(&self) -> bool {
        self.inner.numbered_ter
    }
    #[wasm_bindgen(setter,js_name=numberedTer)]
    pub fn set_numbered_ter(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.numbered_ter = value;
        self.inner.numbered_ter = bool_value(&value, "numberedTer")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=terIgnoresType)]
    pub fn ter_ignores_type(&self) -> bool {
        self.inner.ter_ignores_type
    }
    #[wasm_bindgen(setter,js_name=terIgnoresType)]
    pub fn set_ter_ignores_type(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.ter_ignores_type = value;
        self.inner.ter_ignores_type = bool_value(&value, "terIgnoresType")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=preserveSerial)]
    pub fn preserve_serial(&self) -> bool {
        self.inner.preserve_serial
    }
    #[wasm_bindgen(setter,js_name=preserveSerial)]
    pub fn set_preserve_serial(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.preserve_serial = value;
        self.inner.preserve_serial = bool_value(&value, "preserveSerial")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=endRecord)]
    pub fn end_record(&self) -> bool {
        self.inner.end_record
    }
    #[wasm_bindgen(setter,js_name=endRecord)]
    pub fn set_end_record(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.end_record = value;
        self.inner.end_record = bool_value(&value, "endRecord")?;
        Ok(())
    }
}
#[wasm_bindgen]
pub struct BioMmcifWriteParams {
    pub(crate) inner: ck::BioMmcifWriteParams,
}
#[wasm_bindgen]
impl BioMmcifWriteParams {
    #[wasm_bindgen(constructor)]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] all_groups: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] block_name: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] entry: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] database_status: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] author: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] cell: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] symmetry: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] entity: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] entity_poly: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] struct_ref: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] chem_comp: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] exptl: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] diffrn: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] reflns: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] refine: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] title_keywords: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] ncs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] struct_asym: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] origx: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] struct_conf: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] struct_sheet: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] struct_biol: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] assembly: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] conn: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] cis: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] modres: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] scale: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] atom_type: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] entity_poly_seq: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] tls: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] software: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] group_pdb: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] auth_all: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] prefer_pairs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] compact: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] misuse_hash: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] align_pairs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] align_loops: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::BioMmcifWriteParams::default();
        // COSMolKit❗✔️: atoms: atoms.unwrap_or(all_groups),
        // Fixed scalar options: same Python constructor defaults and no chemistry logic.
        let all_groups = if all_groups.is_undefined() {
            true
        } else {
            bool_value(&all_groups, "allGroups")?
        };
        inner.atoms = if atoms.is_null() || atoms.is_undefined() {
            all_groups
        } else {
            bool_value(&atoms, "atoms")?
        };
        inner.block_name = if block_name.is_null() || block_name.is_undefined() {
            all_groups
        } else {
            bool_value(&block_name, "blockName")?
        };
        inner.entry = if entry.is_null() || entry.is_undefined() {
            all_groups
        } else {
            bool_value(&entry, "entry")?
        };
        inner.database_status = if database_status.is_null() || database_status.is_undefined() {
            all_groups
        } else {
            bool_value(&database_status, "databaseStatus")?
        };
        inner.author = if author.is_null() || author.is_undefined() {
            all_groups
        } else {
            bool_value(&author, "author")?
        };
        inner.cell = if cell.is_null() || cell.is_undefined() {
            all_groups
        } else {
            bool_value(&cell, "cell")?
        };
        inner.symmetry = if symmetry.is_null() || symmetry.is_undefined() {
            all_groups
        } else {
            bool_value(&symmetry, "symmetry")?
        };
        inner.entity = if entity.is_null() || entity.is_undefined() {
            all_groups
        } else {
            bool_value(&entity, "entity")?
        };
        inner.entity_poly = if entity_poly.is_null() || entity_poly.is_undefined() {
            all_groups
        } else {
            bool_value(&entity_poly, "entityPoly")?
        };
        inner.struct_ref = if struct_ref.is_null() || struct_ref.is_undefined() {
            all_groups
        } else {
            bool_value(&struct_ref, "structRef")?
        };
        inner.chem_comp = if chem_comp.is_null() || chem_comp.is_undefined() {
            all_groups
        } else {
            bool_value(&chem_comp, "chemComp")?
        };
        inner.exptl = if exptl.is_null() || exptl.is_undefined() {
            all_groups
        } else {
            bool_value(&exptl, "exptl")?
        };
        inner.diffrn = if diffrn.is_null() || diffrn.is_undefined() {
            all_groups
        } else {
            bool_value(&diffrn, "diffrn")?
        };
        inner.reflns = if reflns.is_null() || reflns.is_undefined() {
            all_groups
        } else {
            bool_value(&reflns, "reflns")?
        };
        inner.refine = if refine.is_null() || refine.is_undefined() {
            all_groups
        } else {
            bool_value(&refine, "refine")?
        };
        inner.title_keywords = if title_keywords.is_null() || title_keywords.is_undefined() {
            all_groups
        } else {
            bool_value(&title_keywords, "titleKeywords")?
        };
        inner.ncs = if ncs.is_null() || ncs.is_undefined() {
            all_groups
        } else {
            bool_value(&ncs, "ncs")?
        };
        inner.struct_asym = if struct_asym.is_null() || struct_asym.is_undefined() {
            all_groups
        } else {
            bool_value(&struct_asym, "structAsym")?
        };
        inner.origx = if origx.is_null() || origx.is_undefined() {
            all_groups
        } else {
            bool_value(&origx, "origx")?
        };
        inner.struct_conf = if struct_conf.is_null() || struct_conf.is_undefined() {
            all_groups
        } else {
            bool_value(&struct_conf, "structConf")?
        };
        inner.struct_sheet = if struct_sheet.is_null() || struct_sheet.is_undefined() {
            all_groups
        } else {
            bool_value(&struct_sheet, "structSheet")?
        };
        inner.struct_biol = if struct_biol.is_null() || struct_biol.is_undefined() {
            all_groups
        } else {
            bool_value(&struct_biol, "structBiol")?
        };
        inner.assembly = if assembly.is_null() || assembly.is_undefined() {
            all_groups
        } else {
            bool_value(&assembly, "assembly")?
        };
        inner.conn = if conn.is_null() || conn.is_undefined() {
            all_groups
        } else {
            bool_value(&conn, "conn")?
        };
        inner.cis = if cis.is_null() || cis.is_undefined() {
            all_groups
        } else {
            bool_value(&cis, "cis")?
        };
        inner.modres = if modres.is_null() || modres.is_undefined() {
            all_groups
        } else {
            bool_value(&modres, "modres")?
        };
        inner.scale = if scale.is_null() || scale.is_undefined() {
            all_groups
        } else {
            bool_value(&scale, "scale")?
        };
        inner.atom_type = if atom_type.is_null() || atom_type.is_undefined() {
            all_groups
        } else {
            bool_value(&atom_type, "atomType")?
        };
        inner.entity_poly_seq = if entity_poly_seq.is_null() || entity_poly_seq.is_undefined() {
            all_groups
        } else {
            bool_value(&entity_poly_seq, "entityPolySeq")?
        };
        inner.tls = if tls.is_null() || tls.is_undefined() {
            all_groups
        } else {
            bool_value(&tls, "tls")?
        };
        inner.software = if software.is_null() || software.is_undefined() {
            all_groups
        } else {
            bool_value(&software, "software")?
        };
        inner.group_pdb = if group_pdb.is_null() || group_pdb.is_undefined() {
            all_groups
        } else {
            bool_value(&group_pdb, "groupPdb")?
        };
        inner.auth_all = if auth_all.is_null() || auth_all.is_undefined() {
            false
        } else {
            bool_value(&auth_all, "authAll")?
        };
        if !prefer_pairs.is_undefined() {
            inner.prefer_pairs = bool_value(&prefer_pairs, "preferPairs")?;
        }
        if !compact.is_undefined() {
            inner.compact = bool_value(&compact, "compact")?;
        }
        if !misuse_hash.is_undefined() {
            inner.misuse_hash = bool_value(&misuse_hash, "misuseHash")?;
        }
        if !align_pairs.is_undefined() {
            inner.align_pairs = u16_value(&align_pairs, "alignPairs")?;
        }
        if !align_loops.is_undefined() {
            inner.align_loops = u16_value(&align_loops, "alignLoops")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=atoms)]
    pub fn atoms(&self) -> bool {
        self.inner.atoms
    }
    #[wasm_bindgen(setter,js_name=atoms)]
    pub fn set_atoms(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.atoms = value;
        self.inner.atoms = bool_value(&value, "atoms")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=blockName)]
    pub fn block_name(&self) -> bool {
        self.inner.block_name
    }
    #[wasm_bindgen(setter,js_name=blockName)]
    pub fn set_block_name(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.block_name = value;
        self.inner.block_name = bool_value(&value, "blockName")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=entry)]
    pub fn entry(&self) -> bool {
        self.inner.entry
    }
    #[wasm_bindgen(setter,js_name=entry)]
    pub fn set_entry(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.entry = value;
        self.inner.entry = bool_value(&value, "entry")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=databaseStatus)]
    pub fn database_status(&self) -> bool {
        self.inner.database_status
    }
    #[wasm_bindgen(setter,js_name=databaseStatus)]
    pub fn set_database_status(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.database_status = value;
        self.inner.database_status = bool_value(&value, "databaseStatus")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=author)]
    pub fn author(&self) -> bool {
        self.inner.author
    }
    #[wasm_bindgen(setter,js_name=author)]
    pub fn set_author(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.author = value;
        self.inner.author = bool_value(&value, "author")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=cell)]
    pub fn cell(&self) -> bool {
        self.inner.cell
    }
    #[wasm_bindgen(setter,js_name=cell)]
    pub fn set_cell(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.cell = value;
        self.inner.cell = bool_value(&value, "cell")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=symmetry)]
    pub fn symmetry(&self) -> bool {
        self.inner.symmetry
    }
    #[wasm_bindgen(setter,js_name=symmetry)]
    pub fn set_symmetry(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.symmetry = value;
        self.inner.symmetry = bool_value(&value, "symmetry")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=entity)]
    pub fn entity(&self) -> bool {
        self.inner.entity
    }
    #[wasm_bindgen(setter,js_name=entity)]
    pub fn set_entity(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.entity = value;
        self.inner.entity = bool_value(&value, "entity")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=entityPoly)]
    pub fn entity_poly(&self) -> bool {
        self.inner.entity_poly
    }
    #[wasm_bindgen(setter,js_name=entityPoly)]
    pub fn set_entity_poly(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.entity_poly = value;
        self.inner.entity_poly = bool_value(&value, "entityPoly")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=structRef)]
    pub fn struct_ref(&self) -> bool {
        self.inner.struct_ref
    }
    #[wasm_bindgen(setter,js_name=structRef)]
    pub fn set_struct_ref(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.struct_ref = value;
        self.inner.struct_ref = bool_value(&value, "structRef")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=chemComp)]
    pub fn chem_comp(&self) -> bool {
        self.inner.chem_comp
    }
    #[wasm_bindgen(setter,js_name=chemComp)]
    pub fn set_chem_comp(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.chem_comp = value;
        self.inner.chem_comp = bool_value(&value, "chemComp")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=exptl)]
    pub fn exptl(&self) -> bool {
        self.inner.exptl
    }
    #[wasm_bindgen(setter,js_name=exptl)]
    pub fn set_exptl(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.exptl = value;
        self.inner.exptl = bool_value(&value, "exptl")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=diffrn)]
    pub fn diffrn(&self) -> bool {
        self.inner.diffrn
    }
    #[wasm_bindgen(setter,js_name=diffrn)]
    pub fn set_diffrn(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.diffrn = value;
        self.inner.diffrn = bool_value(&value, "diffrn")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=reflns)]
    pub fn reflns(&self) -> bool {
        self.inner.reflns
    }
    #[wasm_bindgen(setter,js_name=reflns)]
    pub fn set_reflns(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.reflns = value;
        self.inner.reflns = bool_value(&value, "reflns")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=refine)]
    pub fn refine(&self) -> bool {
        self.inner.refine
    }
    #[wasm_bindgen(setter,js_name=refine)]
    pub fn set_refine(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.refine = value;
        self.inner.refine = bool_value(&value, "refine")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=titleKeywords)]
    pub fn title_keywords(&self) -> bool {
        self.inner.title_keywords
    }
    #[wasm_bindgen(setter,js_name=titleKeywords)]
    pub fn set_title_keywords(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.title_keywords = value;
        self.inner.title_keywords = bool_value(&value, "titleKeywords")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=ncs)]
    pub fn ncs(&self) -> bool {
        self.inner.ncs
    }
    #[wasm_bindgen(setter,js_name=ncs)]
    pub fn set_ncs(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.ncs = value;
        self.inner.ncs = bool_value(&value, "ncs")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=structAsym)]
    pub fn struct_asym(&self) -> bool {
        self.inner.struct_asym
    }
    #[wasm_bindgen(setter,js_name=structAsym)]
    pub fn set_struct_asym(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.struct_asym = value;
        self.inner.struct_asym = bool_value(&value, "structAsym")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=origx)]
    pub fn origx(&self) -> bool {
        self.inner.origx
    }
    #[wasm_bindgen(setter,js_name=origx)]
    pub fn set_origx(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.origx = value;
        self.inner.origx = bool_value(&value, "origx")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=structConf)]
    pub fn struct_conf(&self) -> bool {
        self.inner.struct_conf
    }
    #[wasm_bindgen(setter,js_name=structConf)]
    pub fn set_struct_conf(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.struct_conf = value;
        self.inner.struct_conf = bool_value(&value, "structConf")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=structSheet)]
    pub fn struct_sheet(&self) -> bool {
        self.inner.struct_sheet
    }
    #[wasm_bindgen(setter,js_name=structSheet)]
    pub fn set_struct_sheet(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.struct_sheet = value;
        self.inner.struct_sheet = bool_value(&value, "structSheet")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=structBiol)]
    pub fn struct_biol(&self) -> bool {
        self.inner.struct_biol
    }
    #[wasm_bindgen(setter,js_name=structBiol)]
    pub fn set_struct_biol(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.struct_biol = value;
        self.inner.struct_biol = bool_value(&value, "structBiol")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=assembly)]
    pub fn assembly(&self) -> bool {
        self.inner.assembly
    }
    #[wasm_bindgen(setter,js_name=assembly)]
    pub fn set_assembly(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.assembly = value;
        self.inner.assembly = bool_value(&value, "assembly")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=conn)]
    pub fn conn(&self) -> bool {
        self.inner.conn
    }
    #[wasm_bindgen(setter,js_name=conn)]
    pub fn set_conn(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.conn = value;
        self.inner.conn = bool_value(&value, "conn")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=cis)]
    pub fn cis(&self) -> bool {
        self.inner.cis
    }
    #[wasm_bindgen(setter,js_name=cis)]
    pub fn set_cis(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.cis = value;
        self.inner.cis = bool_value(&value, "cis")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=modres)]
    pub fn modres(&self) -> bool {
        self.inner.modres
    }
    #[wasm_bindgen(setter,js_name=modres)]
    pub fn set_modres(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.modres = value;
        self.inner.modres = bool_value(&value, "modres")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=scale)]
    pub fn scale(&self) -> bool {
        self.inner.scale
    }
    #[wasm_bindgen(setter,js_name=scale)]
    pub fn set_scale(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.scale = value;
        self.inner.scale = bool_value(&value, "scale")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=atomType)]
    pub fn atom_type(&self) -> bool {
        self.inner.atom_type
    }
    #[wasm_bindgen(setter,js_name=atomType)]
    pub fn set_atom_type(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.atom_type = value;
        self.inner.atom_type = bool_value(&value, "atomType")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=entityPolySeq)]
    pub fn entity_poly_seq(&self) -> bool {
        self.inner.entity_poly_seq
    }
    #[wasm_bindgen(setter,js_name=entityPolySeq)]
    pub fn set_entity_poly_seq(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.entity_poly_seq = value;
        self.inner.entity_poly_seq = bool_value(&value, "entityPolySeq")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=tls)]
    pub fn tls(&self) -> bool {
        self.inner.tls
    }
    #[wasm_bindgen(setter,js_name=tls)]
    pub fn set_tls(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.tls = value;
        self.inner.tls = bool_value(&value, "tls")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=software)]
    pub fn software(&self) -> bool {
        self.inner.software
    }
    #[wasm_bindgen(setter,js_name=software)]
    pub fn set_software(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.software = value;
        self.inner.software = bool_value(&value, "software")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=groupPdb)]
    pub fn group_pdb(&self) -> bool {
        self.inner.group_pdb
    }
    #[wasm_bindgen(setter,js_name=groupPdb)]
    pub fn set_group_pdb(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.group_pdb = value;
        self.inner.group_pdb = bool_value(&value, "groupPdb")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=authAll)]
    pub fn auth_all(&self) -> bool {
        self.inner.auth_all
    }
    #[wasm_bindgen(setter,js_name=authAll)]
    pub fn set_auth_all(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.auth_all = value;
        self.inner.auth_all = bool_value(&value, "authAll")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=preferPairs)]
    pub fn prefer_pairs(&self) -> bool {
        self.inner.prefer_pairs
    }
    #[wasm_bindgen(setter,js_name=preferPairs)]
    pub fn set_prefer_pairs(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.prefer_pairs = value;
        self.inner.prefer_pairs = bool_value(&value, "preferPairs")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=compact)]
    pub fn compact(&self) -> bool {
        self.inner.compact
    }
    #[wasm_bindgen(setter,js_name=compact)]
    pub fn set_compact(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.compact = value;
        self.inner.compact = bool_value(&value, "compact")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=misuseHash)]
    pub fn misuse_hash(&self) -> bool {
        self.inner.misuse_hash
    }
    #[wasm_bindgen(setter,js_name=misuseHash)]
    pub fn set_misuse_hash(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.misuse_hash = value;
        self.inner.misuse_hash = bool_value(&value, "misuseHash")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=alignPairs)]
    pub fn align_pairs(&self) -> u16 {
        self.inner.align_pairs
    }
    #[wasm_bindgen(setter,js_name=alignPairs)]
    pub fn set_align_pairs(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.align_pairs = value;
        self.inner.align_pairs = u16_value(&value, "alignPairs")?;
        Ok(())
    }
    #[wasm_bindgen(getter,js_name=alignLoops)]
    pub fn align_loops(&self) -> u16 {
        self.inner.align_loops
    }
    #[wasm_bindgen(setter,js_name=alignLoops)]
    pub fn set_align_loops(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.align_loops = value;
        self.inner.align_loops = u16_value(&value, "alignLoops")?;
        Ok(())
    }
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=toPdb)]
    pub fn to_pdb(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_pdb()
        self.inner
            .to_pdb()
            .map_err(|e| crate::bio_write_errors::pdb_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toPdbWithParams)]
    pub fn to_pdb_with_params(&self, params: &BioPdbWriteParams) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_pdb_with_params(&params.inner)
        self.inner
            .to_pdb_with_params(&params.inner)
            .map_err(|e| crate::bio_write_errors::pdb_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writePdb)]
    pub fn write_pdb(&self, path: &str) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_pdb(Path::new(path))
        self.inner
            .write_pdb(Path::new(path))
            .map_err(|e| crate::bio_write_errors::pdb_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writePdbWithParams)]
    pub fn write_pdb_with_params(
        &self,
        path: &str,
        params: &BioPdbWriteParams,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_pdb_with_params(Path::new(path), &params.inner)
        self.inner
            .write_pdb_with_params(Path::new(path), &params.inner)
            .map_err(|e| crate::bio_write_errors::pdb_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toMmcif)]
    pub fn to_mmcif(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_mmcif()
        self.inner
            .to_mmcif()
            .map_err(|e| crate::bio_write_errors::mmcif_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toMmcifWithParams)]
    pub fn to_mmcif_with_params(&self, params: &BioMmcifWriteParams) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_mmcif_with_params(&params.inner)
        self.inner
            .to_mmcif_with_params(&params.inner)
            .map_err(|e| crate::bio_write_errors::mmcif_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeMmcif)]
    pub fn write_mmcif(&self, path: &str) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_mmcif(Path::new(path))
        self.inner
            .write_mmcif(Path::new(path))
            .map_err(|e| crate::bio_write_errors::mmcif_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeMmcifWithParams)]
    pub fn write_mmcif_with_params(
        &self,
        path: &str,
        params: &BioMmcifWriteParams,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_mmcif_with_params(Path::new(path), &params.inner)
        self.inner
            .write_mmcif_with_params(Path::new(path), &params.inner)
            .map_err(|e| crate::bio_write_errors::mmcif_error(&e).unwrap_or_else(|e| e))
    }
}
