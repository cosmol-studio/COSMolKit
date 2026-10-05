//! QED: migrated existing modern implementation; pinned RDKit QED.py and CPython3.13 sum.
use crate::{DescriptorError, DescriptorInput, DescriptorResult};

#[derive(Debug, Clone, Copy)]
struct QedProperties {
    mw: f64,
    alogp: f64,
    hba: f64,
    hbd: f64,
    psa: f64,
    rotb: f64,
    arom: f64,
    alerts: f64,
}

impl QedProperties {
    fn values_in_rdkit_order(self) -> [f64; 8] {
        // RDKit✔️✔️: QEDproperties = namedtuple('QEDproperties', 'MW,ALOGP,HBA,HBD,PSA,ROTB,AROM,ALERTS')
        [
            self.mw,
            self.alogp,
            self.hba,
            self.hbd,
            self.psa,
            self.rotb,
            self.arom,
            self.alerts,
        ]
    }
}

#[derive(Debug, Clone, Copy)]
struct QedAdsParameter {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    dmax: f64,
}

const RDKIT_QED_WEIGHT_MEAN: QedProperties = QedProperties {
    // RDKit✔️✔️: QEDproperties = namedtuple('QEDproperties', 'MW,ALOGP,HBA,HBD,PSA,ROTB,AROM,ALERTS')
    // RDKit✔️✔️: WEIGHT_MEAN = QEDproperties(0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95)
    mw: 0.66,
    alogp: 0.46,
    hba: 0.05,
    hbd: 0.61,
    psa: 0.06,
    rotb: 0.65,
    arom: 0.48,
    alerts: 0.95,
};

const RDKIT_QED_ALIPHATIC_RINGS_SMARTS: &str = "[$([A;R][!a])]";
const RDKIT_QED_ACCEPTOR_SMARTS: &[&str] = &[
    // RDKit✔️✔️: AcceptorSmarts = [
    // RDKit✔️✔️:   '[oH0;X2]', '[OH1;X2;v2]', '[OH0;X2;v2]', '[OH0;X1;v2]', '[O-;X1]', '[SH0;X2;v2]', '[SH0;X1;v2]',
    // RDKit✔️✔️:   '[S-;X1]', '[nH0;X2]', '[NH0;X1;v3]', '[$([N;+0;X3;v3]);!$(N[C,S]=O)]'
    // RDKit✔️✔️: ]
    "[oH0;X2]",
    "[OH1;X2;v2]",
    "[OH0;X2;v2]",
    "[OH0;X1;v2]",
    "[O-;X1]",
    "[SH0;X2;v2]",
    "[SH0;X1;v2]",
    "[S-;X1]",
    "[nH0;X2]",
    "[NH0;X1;v3]",
    "[$([N;+0;X3;v3]);!$(N[C,S]=O)]",
];

const RDKIT_QED_STRUCTURAL_ALERT_SMARTS: &[&str] = &[
    // RDKit✔️✔️: StructuralAlertSmarts = [
    "*1[O,S,N]*1",
    "[S,C](=[O,S])[F,Br,Cl,I]",
    "[CX4][Cl,Br,I]",
    "[#6]S(=O)(=O)O[#6]",
    "[$([CH]),$(CC)]#CC(=O)[#6]",
    "[$([CH]),$(CC)]#CC(=O)O[#6]",
    "n[OH]",
    "[$([CH]),$(CC)]#CS(=O)(=O)[#6]",
    "C=C(C=O)C=O",
    "n1c([F,Cl,Br,I])cccc1",
    "[CH1](=O)",
    "[#8][#8]",
    "[C;!R]=[N;!R]",
    "[N!R]=[N!R]",
    "[#6](=O)[#6](=O)",
    "[#16][#16]",
    "[#7][NH2]",
    "C(=O)N[NH2]",
    "[#6]=S",
    "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
    "C1(=[O,N])C=CC(=[O,N])C=C1",
    "C1(=[O,N])C(=[O,N])C=CC=C1",
    "a21aa3a(aa1aaaa2)aaaa3",
    "a31a(a2a(aa1)aaaa2)aaaa3",
    "a1aa2a3a(a1)A=AA=A3=AA=A2",
    "c1cc([NH2])ccc1",
    // RDKit✔️✔️:   '[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,' +
    // RDKit✔️✔️:   'Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]', 'I', 'OS(=O)(=O)[O-]', '[N+](=O)[O-]', 'C(=O)N[OH]',
    "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
    "I",
    "OS(=O)(=O)[O-]",
    "[N+](=O)[O-]",
    "C(=O)N[OH]",
    "C1NC(=O)NC(=O)1",
    "[SH]",
    "[S-]",
    "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
    "[CR1]1[CR1][CR1]cc[CR1][CR1]1",
    "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
    "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "C#C",
    "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]",
    "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
    "[#6]=N[OH]",
    "[#6]=NOC=O",
    "[#6](=O)[CX4,CR0X3,O][#6](=O)",
    "c1ccc2c(c1)ccc(=O)o2",
    "[O+,o+,S+,s+]",
    "N=C=O",
    "[NX3,NX4][F,Cl,Br,I]",
    "c1ccccc1OC(=O)[#6]",
    "[CR0]=[CR0][CR0]=[CR0]",
    "[C+,c+,C-,c-]",
    "N=[N+]=[N-]",
    "C12C(NC(N1)=O)CSC2",
    "c1c([OH])c([OH,NH2,NH])ccc1",
    "P",
    "[N,O,S]C#N",
    "C=C=O",
    "[Si][F,Cl,Br,I]",
    "[SX2]O",
    "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
    "O1CCCCC1OC2CCC3CCCCC3C2",
    "N=[CR0][N,n,O,S]",
    "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
    "C=[C!r]C#N",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])",
    "[OH]c1ccc([OH,NH2,NH])cc1",
    "c1ccccc1OC(=O)O",
    "[SX2H0][N]",
    "c12ccccc1(SC(S)=N2)",
    "c12ccccc1(SC(=S)N2)",
    "c1nnnn1C=O",
    "s1c(S)nnc1NC=O",
    "S1C=CSC1=S",
    "C(=O)Onnn",
    "OS(=O)(=O)C(F)(F)F",
    "N#CC[OH]",
    "N#CC(=O)",
    "S(=O)(=O)C#N",
    "N[CH2]C#N",
    "C1(=O)NCC1",
    "S(=O)(=O)[O-,OH]",
    "NC[F,Cl,Br,I]",
    "C=[C!r]O",
    "[NX2+0]=[O+0]",
    "[OR0,NR0][OR0,NR0]",
    "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]",
    "[CX2R0][NX3R0]",
    "c1ccccc1[C;!R]=[C;!R]c2ccccc2",
    "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
    "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
    "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]",
    "[*]=[N+]=[*]",
    "[SX3](=O)[O-,OH]",
    "N#N",
    "F.F.F.F",
    "[R0;D2][R0;D2][R0;D2][R0;D2]",
    "[cR,CR]~C(=O)NC(=O)~[cR,CR]",
    "C=!@CC=[O,S]",
    "[#6,#8,#16][#6](=O)O[#6]",
    "c[C;R0](=[O,S])[#6]",
    "c[SX2][C;!R]",
    "C=C=C",
    "c1nc([F,Cl,Br,I,S])ncc1",
    "c1ncnc([F,Cl,Br,I,S])c1",
    "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
    "[#6]S(=O)(=O)c1ccc(cc1)F",
    "[15N]",
    "[13C]",
    "[18O]",
    "[34S]",
];

const RDKIT_QED_ADS_PARAMETERS: [QedAdsParameter; 8] = [
    // RDKit✔️✔️: adsParameters = {
    // RDKit✔️✔️:   'MW':
    QedAdsParameter {
        a: 2.817065973,
        b: 392.5754953,
        c: 290.7489764,
        d: 2.419764353,
        e: 49.22325677,
        f: 65.37051707,
        dmax: 104.9805561,
    },
    // RDKit✔️✔️:   'ALOGP':
    QedAdsParameter {
        a: 3.172690585,
        b: 137.8624751,
        c: 2.534937431,
        d: 4.581497897,
        e: 0.822739154,
        f: 0.576295591,
        dmax: 131.3186604,
    },
    // RDKit✔️✔️:   'HBA':
    QedAdsParameter {
        a: 2.948620388,
        b: 160.4605972,
        c: 3.615294657,
        d: 4.435986202,
        e: 0.290141953,
        f: 1.300669958,
        dmax: 148.7763046,
    },
    // RDKit✔️✔️:   'HBD':
    QedAdsParameter {
        a: 1.618662227,
        b: 1010.051101,
        c: 0.985094388,
        d: 0.000000001,
        e: 0.713820843,
        f: 0.920922555,
        dmax: 258.1632616,
    },
    // RDKit✔️✔️:   'PSA':
    QedAdsParameter {
        a: 1.876861559,
        b: 125.2232657,
        c: 62.90773554,
        d: 87.83366614,
        e: 12.01999824,
        f: 28.51324732,
        dmax: 104.5686167,
    },
    // RDKit✔️✔️:   'ROTB':
    QedAdsParameter {
        a: 0.010000000,
        b: 272.4121427,
        c: 2.558379970,
        d: 1.565547684,
        e: 1.271567166,
        f: 2.758063707,
        dmax: 105.4420403,
    },
    // RDKit✔️✔️:   'AROM':
    QedAdsParameter {
        a: 3.217788970,
        b: 957.7374108,
        c: 2.274627939,
        d: 0.000000001,
        e: 1.317690384,
        f: 0.375760881,
        dmax: 312.3372610,
    },
    // RDKit✔️✔️:   'ALERTS':
    QedAdsParameter {
        a: 0.010000000,
        b: 1199.094025,
        c: -0.09002883,
        d: 0.000000001,
        e: 0.185904477,
        f: 0.875193782,
        dmax: 417.7253140,
    },
];

fn rdkit_qed_ads(x: f64, p: QedAdsParameter) -> f64 {
    // RDKit✔️✔️: def ads(x, adsParameter):
    // RDKit✔️✔️:   """ ADS function """
    // RDKit✔️✔️:   p = adsParameter
    // RDKit✔️✔️:   exp1 = 1 + math.exp(-1 * (x - p.C + p.D / 2) / p.E)
    let exp1 = 1.0 + (-1.0 * (x - p.c + p.d / 2.0) / p.e).exp();
    // RDKit✔️✔️:   exp2 = 1 + math.exp(-1 * (x - p.C - p.D / 2) / p.F)
    let exp2 = 1.0 + (-1.0 * (x - p.c - p.d / 2.0) / p.f).exp();
    // RDKit✔️✔️:   dx = p.A + p.B / exp1 * (1 - 1 / exp2)
    let dx = p.a + p.b / exp1 * (1.0 - 1.0 / exp2);
    // RDKit✔️✔️:   return dx / p.DMAX
    dx / p.dmax
}

fn rdkit_qed_python313_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    // RDKit 2026.03.1 QED parity is pinned to CPython 3.13.12 because QED.py
    // delegates both reductions to builtins.sum(), whose float algorithm
    // changed in CPython 3.12.
    // CPython✔️✔️: if (PyFloat_CheckExact(result)) {
    // CPython✔️✔️:     double f_result = PyFloat_AS_DOUBLE(result);
    // CPython✔️✔️:     double c = 0.0;
    // CPython✔️✔️:     Py_SETREF(result, NULL);
    // CPython✔️✔️:     while(result == NULL) {
    // CPython✔️✔️:         item = PyIter_Next(iter);
    // CPython✔️✔️:         if (item == NULL) {
    // CPython✔️✔️:             Py_DECREF(iter);
    // CPython✔️✔️:             if (PyErr_Occurred())
    // CPython✔️✔️:                 return NULL;
    // CPython✔️✔️:             /* Avoid losing the sign on a negative result,
    // CPython✔️✔️:                and don't let adding the compensation convert
    // CPython✔️✔️:                an infinite or overflowed sum to a NaN. */
    // CPython✔️✔️:             if (c && Py_IS_FINITE(c)) {
    // CPython✔️✔️:                 f_result += c;
    // CPython✔️✔️:             }
    // CPython✔️✔️:             return PyFloat_FromDouble(f_result);
    // CPython✔️✔️:         }
    // CPython✔️✔️:         if (PyFloat_CheckExact(item)) {
    // CPython✔️✔️:             // Improved Kahan–Babuška algorithm by Arnold Neumaier
    // CPython✔️✔️:             double x = PyFloat_AS_DOUBLE(item);
    // CPython✔️✔️:             double t = f_result + x;
    // CPython✔️✔️:             if (fabs(f_result) >= fabs(x)) {
    // CPython✔️✔️:                 c += (f_result - t) + x;
    // CPython✔️✔️:             } else {
    // CPython✔️✔️:                 c += (x - t) + f_result;
    // CPython✔️✔️:             }
    // CPython✔️✔️:             f_result = t;
    // CPython✔️✔️:             _Py_DECREF_SPECIALIZED(item, _PyFloat_ExactDealloc);
    // CPython✔️✔️:             continue;
    // CPython✔️✔️:         }
    //
    // RDKit QED calls sum() with its default integer-zero start and an all-f64
    // generator. The first float is therefore produced by 0 + item before the
    // CPython float fast path processes the remaining values.
    let mut values = values.into_iter();
    let Some(first) = values.next() else {
        return 0.0;
    };
    let mut result = 0.0 + first;
    let mut compensation = 0.0;
    for value in values {
        let next = result + value;
        if result.abs() >= value.abs() {
            compensation += (result - next) + value;
        } else {
            compensation += (value - next) + result;
        }
        result = next;
    }
    if compensation != 0.0 && compensation.is_finite() {
        result += compensation;
    }
    result
}

fn qed_properties(original: &DescriptorInput<'_>) -> DescriptorResult<QedProperties> {
    // RDKit❗✔️: def properties(mol):
    // RDKit❗✔️:   mol = Chem.RemoveHs(mol)
    // Existing source RemoveHs performs all structural/state transitions on
    // detached owned blocks. No runtime capability or Molecule crosses here.
    let removed = cosmolkit_core::remove_hydrogens_with_params(
        original.topology().clone(),
        original.coordinates().clone(),
        original.properties().clone(),
        &cosmolkit_core::RemoveHsParams::default(),
    )
    .map_err(|source| DescriptorError::Hydrogens {
        function: "qed",
        source,
    })?;
    let valence = removed
        .final_valence
        .as_ref()
        .ok_or(DescriptorError::MissingFinalHydrogenState { field: "valence" })?;
    let rings = removed
        .final_rings
        .as_ref()
        .ok_or(DescriptorError::MissingFinalHydrogenState { field: "rings" })?;
    let input = DescriptorInput::new(
        &removed.topology,
        &removed.coordinates,
        &removed.properties,
        valence,
        rings,
    );
    let mut memo = crate::DescriptorComputedState::default();
    let context = crate::patterns::prepared_context(&input, "qed")?;
    // RDKit❗✔️:   qedProperties = QEDproperties(
    // RDKit❗✔️:     MW=rdmd._CalcMolWt(mol),
    let mw = crate::molecular_weight_with_valence(input.topology(), false, Some(input.valence()))?;
    // RDKit❗✔️:     ALOGP=Crippen.MolLogP(mol),
    let alogp = crate::crippen_totals(&input, true, false, &mut memo)?.logp;
    // RDKit❗✔️:     HBA=sum(
    // RDKit❗✔️:       len(mol.GetSubstructMatches(pattern)) for pattern in Acceptors
    // RDKit❗✔️:       if mol.HasSubstructMatch(pattern)),
    let mut hba = 0u32;
    for &pattern in RDKIT_QED_ACCEPTOR_SMARTS {
        if qed_has_match(&input, pattern, &context)? {
            let matches = crate::patterns::count_pattern_matches_with_context(
                &input, "qed", pattern, &context,
            )?;
            hba = hba
                .checked_add(matches)
                .ok_or(DescriptorError::CountOverflow {
                    function: "qed",
                    field: "HBA",
                })?;
        }
    }
    // RDKit❗✔️:     HBD=rdmd.CalcNumHBD(mol),
    let hbd = crate::num_hbd_prepared(&input)?;
    // RDKit❗✔️:     PSA=MolSurf.TPSA(mol),
    let psa = crate::tpsa(&input, false, false, &mut memo)?;
    // RDKit❗✔️:     ROTB=rdmd.CalcNumRotatableBonds(mol, rdmd.NumRotatableBondsOptions.Strict),
    let rotb = crate::num_rotatable_bonds_prepared(&input, crate::RotatableBondsOptions::Strict)?;
    // RDKit❗✔️:     AROM=len(Chem.GetSSSR(Chem.DeleteSubstructs(Chem.Mol(mol), AliphaticRings))),
    let arom = qed_arom(&input, &context)?;
    // RDKit❗✔️:     ALERTS=sum(1 for alert in StructuralAlerts if mol.HasSubstructMatch(alert)),
    let mut alerts = 0u32;
    for &pattern in RDKIT_QED_STRUCTURAL_ALERT_SMARTS {
        if qed_has_match(&input, pattern, &context)? {
            alerts += 1;
        }
    }
    // RDKit❗✔️:   )
    // RDKit❗✔️:   return qedProperties
    Ok(QedProperties {
        mw,
        alogp,
        hba: f64::from(hba),
        hbd: f64::from(hbd),
        psa,
        rotb: f64::from(rotb),
        arom: f64::from(arom),
        alerts: f64::from(alerts),
    })
}

fn qed_has_match(
    input: &DescriptorInput<'_>,
    pattern: &'static str,
    context: &cosmolkit_search::QueryMatchContext<'_>,
) -> DescriptorResult<bool> {
    // RDKit source (Wrap/substructmethods.h HasSubstructMatch):
    // RDKit❗✔️: ps.maxMatches = 1;
    // Use the one retained query and one matcher, including recursive matching;
    // maximum one match gives the source short-circuit behavior for alerts.
    let query = crate::patterns::retained_pattern("qed", pattern)?;
    let target = cosmolkit_search::SearchTarget::new(
        input.topology(),
        input.coordinates(),
        &input.topology().stereo_groups,
        Some(input.ring_info()),
        Some(input.valence()),
    );
    let params = cosmolkit_search::SubstructMatchParams {
        max_matches: 1,
        ..Default::default()
    };
    let matches = cosmolkit_search::try_get_substruct_matches_with_params_and_context(
        &target, &query, &params, context,
    )
    .map_err(|source| DescriptorError::Search {
        function: "qed",
        source: crate::DescriptorSearchCause::Match(source),
    })?;
    Ok(!matches.is_empty())
}

fn qed_arom(
    input: &DescriptorInput<'_>,
    context: &cosmolkit_search::QueryMatchContext<'_>,
) -> DescriptorResult<u32> {
    // Source QED calls DeleteSubstructs with onlyFrags/useChirality both false.
    // It only reads GetSSSR of the returned graph, so coordinate/property
    // transport is not an input here. Reuse the detached model batch editor.
    // RDKit❗❌: SubstructMatch(*res, query, fgpMatches, uniquify, recursionPossible,
    // RDKit❗❌:                useChirality);
    let matches = crate::patterns::pattern_matches_with_context(
        input,
        "qed",
        RDKIT_QED_ALIPHATIC_RINGS_SMARTS,
        context,
    )?;
    let topology = if matches.is_empty() {
        std::borrow::Cow::Borrowed(input.topology())
    } else {
        // RDKit❗❌: for (const auto &mxi : matches) {
        // RDKit❗❌:   INT_VECT tmp;
        // RDKit❗❌:   Union(mxi, delList, tmp);
        // RDKit❗❌:   delList = tmp;
        // RDKit❗❌: }
        let mut atoms = std::collections::BTreeSet::new();
        for matched in matches {
            atoms.extend(matched.atom_mapping);
        }
        // RDKit❗❌: res->beginBatchEdit();
        // RDKit❗❌: for (auto idx : delList) {
        // RDKit❗❌:   res->removeAtom(idx);
        // RDKit❗❌: }
        // RDKit❗❌: res->commitBatchEdit();
        // Batch edit has existing two-topology snapshot allocation debt; no
        // duplicate graph compaction algorithm is introduced in descriptors.
        let mut edit = input.topology().begin_batch_edit().map_err(|source| {
            DescriptorError::TopologyEdit {
                function: "qed",
                source,
            }
        })?;
        for atom in atoms {
            edit.remove_atom(cosmolkit_model::AtomId::new(atom))
                .map_err(|source| DescriptorError::TopologyEdit {
                    function: "qed",
                    source,
                })?;
        }
        let (topology, _) = edit
            .finish()
            .map_err(|source| DescriptorError::TopologyEdit {
                function: "qed",
                source,
            })?;
        std::borrow::Cow::Owned(topology)
    };
    // RDKit❗✔️: AROM=len(Chem.GetSSSR(Chem.DeleteSubstructs(Chem.Mol(mol), AliphaticRings))),
    // This MUST be ordinary SSSR, not symmetrized SSSR. The existing owner
    // borrows only graph rows, independent of computed valence properties.
    let rings = cosmolkit_core::find_sssr_from_parts(
        topology.atoms.len(),
        &topology.bonds,
        &topology.adjacency,
    )
    .map_err(|source| DescriptorError::Ring {
        function: "qed",
        source,
    })?;
    u32::try_from(rings.atom_rings().len()).map_err(|_| DescriptorError::CountOverflow {
        function: "qed",
        field: "AROM",
    })
}

pub fn qed(input: &DescriptorInput<'_>) -> DescriptorResult<f64> {
    // RDKit✔️✔️: def qed(mol, w=WEIGHT_MEAN, qedProperties=None):
    // RDKit✔️✔️:   """ Calculate the weighted sum of ADS mapped properties
    // RDKit✔️✔️:
    // RDKit✔️✔️:   some examples from the QED paper, reference values from Peter G's original implementation
    // RDKit✔️✔️:   >>> m = Chem.MolFromSmiles('N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O')
    // RDKit✔️✔️:   >>> qed(m)
    // RDKit✔️✔️:   0.253...
    // RDKit✔️✔️:   >>> m = Chem.MolFromSmiles('CNC(=NCCSCc1nc[nH]c1C)NC#N')
    // RDKit✔️✔️:   >>> qed(m)
    // RDKit✔️✔️:   0.234...
    // RDKit✔️✔️:   >>> m = Chem.MolFromSmiles('CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12')
    // RDKit✔️✔️:   >>> qed(m)
    // RDKit✔️✔️:   0.234...
    // RDKit✔️✔️:   """
    // RDKit✔️✔️:   if qedProperties is None:
    // RDKit✔️✔️:     qedProperties = properties(mol)
    let qed_properties = qed_properties(input)?;
    // RDKit✔️✔️:   d = [ads(pi, adsParameters[name]) for name, pi in qedProperties._asdict().items()]
    let values = qed_properties.values_in_rdkit_order();
    let d = values
        .into_iter()
        .zip(RDKIT_QED_ADS_PARAMETERS)
        .map(|(pi, parameter)| rdkit_qed_ads(pi, parameter))
        .collect::<Vec<_>>();
    // RDKit✔️✔️:   t = sum(wi * math.log(di) for wi, di in zip(w, d))
    let weights = RDKIT_QED_WEIGHT_MEAN.values_in_rdkit_order();
    let t = rdkit_qed_python313_sum(
        weights
            .iter()
            .copied()
            .zip(d.iter().copied())
            .map(|(wi, di)| wi * di.ln()),
    );
    // RDKit✔️✔️:   return math.exp(t / sum(w))
    Ok((t / rdkit_qed_python313_sum(weights)).exp())
}
