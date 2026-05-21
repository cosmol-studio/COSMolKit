pub(crate) const GEMMI_OP_DEN: i32 = 24;

#[derive(Clone, Copy, Debug)]
pub(crate) struct GemmiSymOp {
    pub rot: [[i32; 3]; 3],
    pub tran: [i32; 3],
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct GemmiSpaceGroupEntry {
    pub hm: &'static str,
    pub ext: u8,
    pub ops: &'static [GemmiSymOp],
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct GemmiAltNameEntry {
    pub hm: &'static str,
    pub ext: u8,
    pub pos: usize,
}

const OPS_0: &[GemmiSymOp] = &[GemmiSymOp {
    rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
    tran: [0, 0, 0],
}];

const OPS_1: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_2: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_3: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_4: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_5: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_6: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_7: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_8: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_9: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_10: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_11: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_12: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_13: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_14: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_15: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_16: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_17: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_18: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_19: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_20: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_21: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_22: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_23: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_24: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_25: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_26: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_27: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_28: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_29: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_30: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_31: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_32: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_33: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_34: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_35: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_36: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_37: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_38: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_39: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_40: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_41: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_42: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_43: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_44: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_45: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_46: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_47: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_48: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_49: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_50: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_51: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_52: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_53: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_54: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_55: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_56: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_57: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_58: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_59: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_60: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_61: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_62: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_63: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_64: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_65: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_66: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_67: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_68: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_69: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_70: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_71: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_72: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_73: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_74: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_75: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_76: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_77: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_78: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_79: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_80: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_81: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_82: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_83: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_84: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_85: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_86: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_87: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_88: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_89: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_90: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_91: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_92: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_93: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_94: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_95: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_96: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_97: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_98: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_99: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_100: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_101: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_102: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_103: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_104: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_105: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_106: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_107: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_108: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_109: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_110: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_111: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_112: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_113: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_114: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_115: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_116: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_117: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_118: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_119: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_120: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_121: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_122: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_123: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_124: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_125: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_126: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_127: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_128: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_129: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_130: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_131: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_132: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_133: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_134: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_135: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_136: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_137: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_138: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_139: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_140: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_141: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_142: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_143: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_144: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_145: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_146: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_147: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_148: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_149: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_150: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_151: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_152: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_153: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_154: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_155: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_156: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_157: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_158: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_159: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_160: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_161: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_162: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_163: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_164: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_165: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_166: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_167: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_168: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_169: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_170: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_171: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_172: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_173: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_174: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_175: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_176: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_177: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_178: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_179: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_180: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_181: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_182: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_183: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_184: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_185: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_186: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_187: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_188: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_189: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_190: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_191: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_192: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_193: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_194: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_195: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_196: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_197: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_198: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_199: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_200: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_201: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_202: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_203: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_204: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_205: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_206: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_207: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_208: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_209: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_210: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_211: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
];

const OPS_212: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
];

const OPS_213: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_214: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_215: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_216: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_217: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_218: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_219: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_220: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_221: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_222: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_223: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_224: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_225: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_226: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_227: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_228: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_229: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_230: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_231: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_232: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_233: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_234: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_235: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_236: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_237: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_238: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_239: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_240: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_241: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_242: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_243: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_244: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_245: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_246: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_247: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_248: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_249: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_250: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_251: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_252: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_253: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_254: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_255: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_256: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_257: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_258: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_259: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_260: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_261: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_262: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_263: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_264: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_265: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_266: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_267: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_268: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_269: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_270: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_271: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_272: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_273: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_274: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_275: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_276: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_277: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_278: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_279: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_280: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_281: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_282: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_283: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_284: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_285: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_286: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_287: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_288: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_289: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_290: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_291: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_292: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_293: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_294: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_295: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_296: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_297: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_298: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_299: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_300: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_301: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_302: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_303: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_304: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_305: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_306: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_307: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_308: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_309: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_310: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_311: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_312: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_313: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_314: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_315: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_316: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_317: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_318: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_319: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_320: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_321: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_322: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_323: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_324: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_325: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_326: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_327: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_328: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_329: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_330: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_331: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_332: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_333: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_334: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
];

const OPS_335: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
];

const OPS_336: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_337: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_338: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_339: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_340: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_341: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_342: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_343: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_344: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_345: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_346: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_347: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_348: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_349: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
];

const OPS_350: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_351: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
];

const OPS_352: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_353: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
];

const OPS_354: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_355: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_356: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_357: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_358: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_359: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_360: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_361: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_362: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_363: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_364: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
];

const OPS_365: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_366: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_367: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 18],
    },
];

const OPS_368: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_369: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_370: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_371: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 6],
    },
];

const OPS_372: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_373: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_374: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_375: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_376: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_377: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_378: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_379: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_380: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_381: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_382: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_383: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_384: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_385: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
];

const OPS_386: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
];

const OPS_387: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_388: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_389: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_390: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_391: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_392: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_393: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_394: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_395: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_396: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_397: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_398: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
];

const OPS_399: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_400: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_401: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_402: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_403: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_404: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_405: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_406: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_407: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_408: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_409: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_410: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_411: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_412: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_413: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_414: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_415: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_416: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_417: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_418: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_419: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_420: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_421: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_422: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_423: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_424: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_425: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
];

const OPS_426: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
];

const OPS_427: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
];

const OPS_428: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
];

const OPS_429: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_430: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
];

const OPS_431: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
];

const OPS_432: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
];

const OPS_433: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_434: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_435: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
];

const OPS_436: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_437: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_438: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_439: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_440: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
];

const OPS_441: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_442: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
];

const OPS_443: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
];

const OPS_444: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_445: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_446: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_447: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_448: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_449: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
];

const OPS_450: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_451: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
];

const OPS_452: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_453: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_454: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_455: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_456: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_457: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
];

const OPS_458: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_459: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 4],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 4],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 4],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
];

const OPS_460: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_461: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_462: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
];

const OPS_463: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
];

const OPS_464: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
];

const OPS_465: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
];

const OPS_466: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_467: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_468: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_469: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_470: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_471: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 20],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 4],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_472: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 4],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 20],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_473: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_474: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_475: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_476: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_477: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_478: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_479: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_480: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_481: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_482: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_483: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_484: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_485: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_486: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_487: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_488: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_489: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_490: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_491: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_492: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
];

const OPS_493: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_494: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_495: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
];

const OPS_496: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_497: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
];

const OPS_498: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 6, 18],
    },
];

const OPS_499: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_500: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_501: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 12],
    },
];

const OPS_502: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_503: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_504: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_505: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_506: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_507: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_508: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_509: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_510: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_511: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_512: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_513: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_514: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
];

const OPS_515: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_516: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_517: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_518: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_519: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_520: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_521: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_522: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_523: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
];

const OPS_524: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
];

const OPS_525: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
];

const OPS_526: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
];

const OPS_527: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 0, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 12, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 0, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 0, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 12, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
];

const OPS_528: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_529: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_530: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_531: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_532: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_533: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_534: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_535: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_536: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_537: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_538: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
];

const OPS_539: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_540: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_541: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_542: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_543: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_544: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_545: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_546: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_547: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_548: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_549: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_550: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_551: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_552: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_553: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_554: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_555: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
];

const OPS_556: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_557: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_558: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_559: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_560: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_561: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_562: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_563: &[GemmiSymOp] = &[
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    GemmiSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

pub(crate) const GEMMI_SPACEGROUPS: &[GemmiSpaceGroupEntry] = &[
    GemmiSpaceGroupEntry {
        hm: "P 1",
        ext: 0,
        ops: OPS_0,
    },
    GemmiSpaceGroupEntry {
        hm: "P -1",
        ext: 0,
        ops: OPS_1,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 2 1",
        ext: 0,
        ops: OPS_2,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 2",
        ext: 0,
        ops: OPS_3,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 1 1",
        ext: 0,
        ops: OPS_4,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 21 1",
        ext: 0,
        ops: OPS_5,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 21",
        ext: 0,
        ops: OPS_6,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 1 1",
        ext: 0,
        ops: OPS_7,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 2 1",
        ext: 0,
        ops: OPS_8,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 2 1",
        ext: 0,
        ops: OPS_9,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 2 1",
        ext: 0,
        ops: OPS_10,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 1 2",
        ext: 0,
        ops: OPS_11,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 1 2",
        ext: 0,
        ops: OPS_12,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 1 2",
        ext: 0,
        ops: OPS_13,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2 1 1",
        ext: 0,
        ops: OPS_14,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 1 1",
        ext: 0,
        ops: OPS_15,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 1 1",
        ext: 0,
        ops: OPS_16,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 m 1",
        ext: 0,
        ops: OPS_17,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 m",
        ext: 0,
        ops: OPS_18,
    },
    GemmiSpaceGroupEntry {
        hm: "P m 1 1",
        ext: 0,
        ops: OPS_19,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 c 1",
        ext: 0,
        ops: OPS_20,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 n 1",
        ext: 0,
        ops: OPS_21,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 a 1",
        ext: 0,
        ops: OPS_22,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 a",
        ext: 0,
        ops: OPS_23,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 n",
        ext: 0,
        ops: OPS_24,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 b",
        ext: 0,
        ops: OPS_25,
    },
    GemmiSpaceGroupEntry {
        hm: "P b 1 1",
        ext: 0,
        ops: OPS_26,
    },
    GemmiSpaceGroupEntry {
        hm: "P n 1 1",
        ext: 0,
        ops: OPS_27,
    },
    GemmiSpaceGroupEntry {
        hm: "P c 1 1",
        ext: 0,
        ops: OPS_28,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 m 1",
        ext: 0,
        ops: OPS_29,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 m 1",
        ext: 0,
        ops: OPS_30,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 m 1",
        ext: 0,
        ops: OPS_31,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 1 m",
        ext: 0,
        ops: OPS_32,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 1 m",
        ext: 0,
        ops: OPS_33,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 1 m",
        ext: 0,
        ops: OPS_34,
    },
    GemmiSpaceGroupEntry {
        hm: "B m 1 1",
        ext: 0,
        ops: OPS_35,
    },
    GemmiSpaceGroupEntry {
        hm: "C m 1 1",
        ext: 0,
        ops: OPS_36,
    },
    GemmiSpaceGroupEntry {
        hm: "I m 1 1",
        ext: 0,
        ops: OPS_37,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 c 1",
        ext: 0,
        ops: OPS_38,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 n 1",
        ext: 0,
        ops: OPS_39,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 a 1",
        ext: 0,
        ops: OPS_40,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 a 1",
        ext: 0,
        ops: OPS_41,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 n 1",
        ext: 0,
        ops: OPS_42,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 c 1",
        ext: 0,
        ops: OPS_43,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 1 a",
        ext: 0,
        ops: OPS_44,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 1 n",
        ext: 0,
        ops: OPS_45,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 1 b",
        ext: 0,
        ops: OPS_46,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 1 b",
        ext: 0,
        ops: OPS_47,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 1 n",
        ext: 0,
        ops: OPS_48,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 1 a",
        ext: 0,
        ops: OPS_49,
    },
    GemmiSpaceGroupEntry {
        hm: "B b 1 1",
        ext: 0,
        ops: OPS_50,
    },
    GemmiSpaceGroupEntry {
        hm: "C n 1 1",
        ext: 0,
        ops: OPS_51,
    },
    GemmiSpaceGroupEntry {
        hm: "I c 1 1",
        ext: 0,
        ops: OPS_52,
    },
    GemmiSpaceGroupEntry {
        hm: "C c 1 1",
        ext: 0,
        ops: OPS_53,
    },
    GemmiSpaceGroupEntry {
        hm: "B n 1 1",
        ext: 0,
        ops: OPS_54,
    },
    GemmiSpaceGroupEntry {
        hm: "I b 1 1",
        ext: 0,
        ops: OPS_55,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 2/m 1",
        ext: 0,
        ops: OPS_56,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 2/m",
        ext: 0,
        ops: OPS_57,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2/m 1 1",
        ext: 0,
        ops: OPS_58,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 21/m 1",
        ext: 0,
        ops: OPS_59,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 21/m",
        ext: 0,
        ops: OPS_60,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21/m 1 1",
        ext: 0,
        ops: OPS_61,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 2/m 1",
        ext: 0,
        ops: OPS_62,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 2/m 1",
        ext: 0,
        ops: OPS_63,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 2/m 1",
        ext: 0,
        ops: OPS_64,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 1 2/m",
        ext: 0,
        ops: OPS_65,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 1 2/m",
        ext: 0,
        ops: OPS_66,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 1 2/m",
        ext: 0,
        ops: OPS_67,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2/m 1 1",
        ext: 0,
        ops: OPS_68,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2/m 1 1",
        ext: 0,
        ops: OPS_69,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2/m 1 1",
        ext: 0,
        ops: OPS_70,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 2/c 1",
        ext: 0,
        ops: OPS_71,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 2/n 1",
        ext: 0,
        ops: OPS_72,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 2/a 1",
        ext: 0,
        ops: OPS_73,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 2/a",
        ext: 0,
        ops: OPS_74,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 2/n",
        ext: 0,
        ops: OPS_75,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 2/b",
        ext: 0,
        ops: OPS_76,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2/b 1 1",
        ext: 0,
        ops: OPS_77,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2/n 1 1",
        ext: 0,
        ops: OPS_78,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2/c 1 1",
        ext: 0,
        ops: OPS_79,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 21/c 1",
        ext: 0,
        ops: OPS_80,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 21/n 1",
        ext: 0,
        ops: OPS_81,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 21/a 1",
        ext: 0,
        ops: OPS_82,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 21/a",
        ext: 0,
        ops: OPS_83,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 21/n",
        ext: 0,
        ops: OPS_84,
    },
    GemmiSpaceGroupEntry {
        hm: "P 1 1 21/b",
        ext: 0,
        ops: OPS_85,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21/b 1 1",
        ext: 0,
        ops: OPS_86,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21/n 1 1",
        ext: 0,
        ops: OPS_87,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21/c 1 1",
        ext: 0,
        ops: OPS_88,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 2/c 1",
        ext: 0,
        ops: OPS_89,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 2/n 1",
        ext: 0,
        ops: OPS_90,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 2/a 1",
        ext: 0,
        ops: OPS_91,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 2/a 1",
        ext: 0,
        ops: OPS_92,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 2/n 1",
        ext: 0,
        ops: OPS_93,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 2/c 1",
        ext: 0,
        ops: OPS_94,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 1 2/a",
        ext: 0,
        ops: OPS_95,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 1 2/n",
        ext: 0,
        ops: OPS_96,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 1 2/b",
        ext: 0,
        ops: OPS_97,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 1 2/b",
        ext: 0,
        ops: OPS_98,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1 1 2/n",
        ext: 0,
        ops: OPS_99,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 1 2/a",
        ext: 0,
        ops: OPS_100,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2/b 1 1",
        ext: 0,
        ops: OPS_101,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2/n 1 1",
        ext: 0,
        ops: OPS_102,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2/c 1 1",
        ext: 0,
        ops: OPS_103,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2/c 1 1",
        ext: 0,
        ops: OPS_104,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2/n 1 1",
        ext: 0,
        ops: OPS_105,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2/b 1 1",
        ext: 0,
        ops: OPS_106,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 2 2",
        ext: 0,
        ops: OPS_107,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 2 21",
        ext: 0,
        ops: OPS_108,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 2 2",
        ext: 0,
        ops: OPS_109,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 21 2",
        ext: 0,
        ops: OPS_110,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 21 2",
        ext: 0,
        ops: OPS_111,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 21 21",
        ext: 0,
        ops: OPS_112,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 2 21",
        ext: 0,
        ops: OPS_113,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 21 21",
        ext: 0,
        ops: OPS_114,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 2 21",
        ext: 0,
        ops: OPS_115,
    },
    GemmiSpaceGroupEntry {
        hm: "A 21 2 2",
        ext: 0,
        ops: OPS_116,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2 21 2",
        ext: 0,
        ops: OPS_117,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 2 2",
        ext: 0,
        ops: OPS_118,
    },
    GemmiSpaceGroupEntry {
        hm: "A 2 2 2",
        ext: 0,
        ops: OPS_119,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2 2 2",
        ext: 0,
        ops: OPS_120,
    },
    GemmiSpaceGroupEntry {
        hm: "F 2 2 2",
        ext: 0,
        ops: OPS_121,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 2 2",
        ext: 0,
        ops: OPS_122,
    },
    GemmiSpaceGroupEntry {
        hm: "I 21 21 21",
        ext: 0,
        ops: OPS_123,
    },
    GemmiSpaceGroupEntry {
        hm: "P m m 2",
        ext: 0,
        ops: OPS_124,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 m m",
        ext: 0,
        ops: OPS_125,
    },
    GemmiSpaceGroupEntry {
        hm: "P m 2 m",
        ext: 0,
        ops: OPS_126,
    },
    GemmiSpaceGroupEntry {
        hm: "P m c 21",
        ext: 0,
        ops: OPS_127,
    },
    GemmiSpaceGroupEntry {
        hm: "P c m 21",
        ext: 0,
        ops: OPS_128,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 m a",
        ext: 0,
        ops: OPS_129,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 a m",
        ext: 0,
        ops: OPS_130,
    },
    GemmiSpaceGroupEntry {
        hm: "P b 21 m",
        ext: 0,
        ops: OPS_131,
    },
    GemmiSpaceGroupEntry {
        hm: "P m 21 b",
        ext: 0,
        ops: OPS_132,
    },
    GemmiSpaceGroupEntry {
        hm: "P c c 2",
        ext: 0,
        ops: OPS_133,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 a a",
        ext: 0,
        ops: OPS_134,
    },
    GemmiSpaceGroupEntry {
        hm: "P b 2 b",
        ext: 0,
        ops: OPS_135,
    },
    GemmiSpaceGroupEntry {
        hm: "P m a 2",
        ext: 0,
        ops: OPS_136,
    },
    GemmiSpaceGroupEntry {
        hm: "P b m 2",
        ext: 0,
        ops: OPS_137,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 m b",
        ext: 0,
        ops: OPS_138,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 c m",
        ext: 0,
        ops: OPS_139,
    },
    GemmiSpaceGroupEntry {
        hm: "P c 2 m",
        ext: 0,
        ops: OPS_140,
    },
    GemmiSpaceGroupEntry {
        hm: "P m 2 a",
        ext: 0,
        ops: OPS_141,
    },
    GemmiSpaceGroupEntry {
        hm: "P c a 21",
        ext: 0,
        ops: OPS_142,
    },
    GemmiSpaceGroupEntry {
        hm: "P b c 21",
        ext: 0,
        ops: OPS_143,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 a b",
        ext: 0,
        ops: OPS_144,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 c a",
        ext: 0,
        ops: OPS_145,
    },
    GemmiSpaceGroupEntry {
        hm: "P c 21 b",
        ext: 0,
        ops: OPS_146,
    },
    GemmiSpaceGroupEntry {
        hm: "P b 21 a",
        ext: 0,
        ops: OPS_147,
    },
    GemmiSpaceGroupEntry {
        hm: "P n c 2",
        ext: 0,
        ops: OPS_148,
    },
    GemmiSpaceGroupEntry {
        hm: "P c n 2",
        ext: 0,
        ops: OPS_149,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 n a",
        ext: 0,
        ops: OPS_150,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 a n",
        ext: 0,
        ops: OPS_151,
    },
    GemmiSpaceGroupEntry {
        hm: "P b 2 n",
        ext: 0,
        ops: OPS_152,
    },
    GemmiSpaceGroupEntry {
        hm: "P n 2 b",
        ext: 0,
        ops: OPS_153,
    },
    GemmiSpaceGroupEntry {
        hm: "P m n 21",
        ext: 0,
        ops: OPS_154,
    },
    GemmiSpaceGroupEntry {
        hm: "P n m 21",
        ext: 0,
        ops: OPS_155,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 m n",
        ext: 0,
        ops: OPS_156,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 n m",
        ext: 0,
        ops: OPS_157,
    },
    GemmiSpaceGroupEntry {
        hm: "P n 21 m",
        ext: 0,
        ops: OPS_158,
    },
    GemmiSpaceGroupEntry {
        hm: "P m 21 n",
        ext: 0,
        ops: OPS_159,
    },
    GemmiSpaceGroupEntry {
        hm: "P b a 2",
        ext: 0,
        ops: OPS_160,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 c b",
        ext: 0,
        ops: OPS_161,
    },
    GemmiSpaceGroupEntry {
        hm: "P c 2 a",
        ext: 0,
        ops: OPS_162,
    },
    GemmiSpaceGroupEntry {
        hm: "P n a 21",
        ext: 0,
        ops: OPS_163,
    },
    GemmiSpaceGroupEntry {
        hm: "P b n 21",
        ext: 0,
        ops: OPS_164,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 n b",
        ext: 0,
        ops: OPS_165,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 c n",
        ext: 0,
        ops: OPS_166,
    },
    GemmiSpaceGroupEntry {
        hm: "P c 21 n",
        ext: 0,
        ops: OPS_167,
    },
    GemmiSpaceGroupEntry {
        hm: "P n 21 a",
        ext: 0,
        ops: OPS_168,
    },
    GemmiSpaceGroupEntry {
        hm: "P n n 2",
        ext: 0,
        ops: OPS_169,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 n n",
        ext: 0,
        ops: OPS_170,
    },
    GemmiSpaceGroupEntry {
        hm: "P n 2 n",
        ext: 0,
        ops: OPS_171,
    },
    GemmiSpaceGroupEntry {
        hm: "C m m 2",
        ext: 0,
        ops: OPS_172,
    },
    GemmiSpaceGroupEntry {
        hm: "A 2 m m",
        ext: 0,
        ops: OPS_173,
    },
    GemmiSpaceGroupEntry {
        hm: "B m 2 m",
        ext: 0,
        ops: OPS_174,
    },
    GemmiSpaceGroupEntry {
        hm: "C m c 21",
        ext: 0,
        ops: OPS_175,
    },
    GemmiSpaceGroupEntry {
        hm: "C c m 21",
        ext: 0,
        ops: OPS_176,
    },
    GemmiSpaceGroupEntry {
        hm: "A 21 m a",
        ext: 0,
        ops: OPS_177,
    },
    GemmiSpaceGroupEntry {
        hm: "A 21 a m",
        ext: 0,
        ops: OPS_178,
    },
    GemmiSpaceGroupEntry {
        hm: "B b 21 m",
        ext: 0,
        ops: OPS_179,
    },
    GemmiSpaceGroupEntry {
        hm: "B m 21 b",
        ext: 0,
        ops: OPS_180,
    },
    GemmiSpaceGroupEntry {
        hm: "C c c 2",
        ext: 0,
        ops: OPS_181,
    },
    GemmiSpaceGroupEntry {
        hm: "A 2 a a",
        ext: 0,
        ops: OPS_182,
    },
    GemmiSpaceGroupEntry {
        hm: "B b 2 b",
        ext: 0,
        ops: OPS_183,
    },
    GemmiSpaceGroupEntry {
        hm: "A m m 2",
        ext: 0,
        ops: OPS_184,
    },
    GemmiSpaceGroupEntry {
        hm: "B m m 2",
        ext: 0,
        ops: OPS_185,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2 m m",
        ext: 0,
        ops: OPS_186,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 m m",
        ext: 0,
        ops: OPS_187,
    },
    GemmiSpaceGroupEntry {
        hm: "C m 2 m",
        ext: 0,
        ops: OPS_188,
    },
    GemmiSpaceGroupEntry {
        hm: "A m 2 m",
        ext: 0,
        ops: OPS_189,
    },
    GemmiSpaceGroupEntry {
        hm: "A b m 2",
        ext: 0,
        ops: OPS_190,
    },
    GemmiSpaceGroupEntry {
        hm: "B m a 2",
        ext: 0,
        ops: OPS_191,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2 c m",
        ext: 0,
        ops: OPS_192,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 m b",
        ext: 0,
        ops: OPS_193,
    },
    GemmiSpaceGroupEntry {
        hm: "C m 2 a",
        ext: 0,
        ops: OPS_194,
    },
    GemmiSpaceGroupEntry {
        hm: "A c 2 m",
        ext: 0,
        ops: OPS_195,
    },
    GemmiSpaceGroupEntry {
        hm: "A m a 2",
        ext: 0,
        ops: OPS_196,
    },
    GemmiSpaceGroupEntry {
        hm: "B b m 2",
        ext: 0,
        ops: OPS_197,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2 m b",
        ext: 0,
        ops: OPS_198,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 c m",
        ext: 0,
        ops: OPS_199,
    },
    GemmiSpaceGroupEntry {
        hm: "C c 2 m",
        ext: 0,
        ops: OPS_200,
    },
    GemmiSpaceGroupEntry {
        hm: "A m 2 a",
        ext: 0,
        ops: OPS_201,
    },
    GemmiSpaceGroupEntry {
        hm: "A b a 2",
        ext: 0,
        ops: OPS_202,
    },
    GemmiSpaceGroupEntry {
        hm: "B b a 2",
        ext: 0,
        ops: OPS_203,
    },
    GemmiSpaceGroupEntry {
        hm: "B 2 c b",
        ext: 0,
        ops: OPS_204,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 c b",
        ext: 0,
        ops: OPS_205,
    },
    GemmiSpaceGroupEntry {
        hm: "C c 2 a",
        ext: 0,
        ops: OPS_206,
    },
    GemmiSpaceGroupEntry {
        hm: "A c 2 a",
        ext: 0,
        ops: OPS_207,
    },
    GemmiSpaceGroupEntry {
        hm: "F m m 2",
        ext: 0,
        ops: OPS_208,
    },
    GemmiSpaceGroupEntry {
        hm: "F 2 m m",
        ext: 0,
        ops: OPS_209,
    },
    GemmiSpaceGroupEntry {
        hm: "F m 2 m",
        ext: 0,
        ops: OPS_210,
    },
    GemmiSpaceGroupEntry {
        hm: "F d d 2",
        ext: 0,
        ops: OPS_211,
    },
    GemmiSpaceGroupEntry {
        hm: "F 2 d d",
        ext: 0,
        ops: OPS_212,
    },
    GemmiSpaceGroupEntry {
        hm: "F d 2 d",
        ext: 0,
        ops: OPS_213,
    },
    GemmiSpaceGroupEntry {
        hm: "I m m 2",
        ext: 0,
        ops: OPS_214,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 m m",
        ext: 0,
        ops: OPS_215,
    },
    GemmiSpaceGroupEntry {
        hm: "I m 2 m",
        ext: 0,
        ops: OPS_216,
    },
    GemmiSpaceGroupEntry {
        hm: "I b a 2",
        ext: 0,
        ops: OPS_217,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 c b",
        ext: 0,
        ops: OPS_218,
    },
    GemmiSpaceGroupEntry {
        hm: "I c 2 a",
        ext: 0,
        ops: OPS_219,
    },
    GemmiSpaceGroupEntry {
        hm: "I m a 2",
        ext: 0,
        ops: OPS_220,
    },
    GemmiSpaceGroupEntry {
        hm: "I b m 2",
        ext: 0,
        ops: OPS_221,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 m b",
        ext: 0,
        ops: OPS_222,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 c m",
        ext: 0,
        ops: OPS_223,
    },
    GemmiSpaceGroupEntry {
        hm: "I c 2 m",
        ext: 0,
        ops: OPS_224,
    },
    GemmiSpaceGroupEntry {
        hm: "I m 2 a",
        ext: 0,
        ops: OPS_225,
    },
    GemmiSpaceGroupEntry {
        hm: "P m m m",
        ext: 0,
        ops: OPS_226,
    },
    GemmiSpaceGroupEntry {
        hm: "P n n n",
        ext: 49,
        ops: OPS_227,
    },
    GemmiSpaceGroupEntry {
        hm: "P n n n",
        ext: 50,
        ops: OPS_228,
    },
    GemmiSpaceGroupEntry {
        hm: "P c c m",
        ext: 0,
        ops: OPS_229,
    },
    GemmiSpaceGroupEntry {
        hm: "P m a a",
        ext: 0,
        ops: OPS_230,
    },
    GemmiSpaceGroupEntry {
        hm: "P b m b",
        ext: 0,
        ops: OPS_231,
    },
    GemmiSpaceGroupEntry {
        hm: "P b a n",
        ext: 49,
        ops: OPS_232,
    },
    GemmiSpaceGroupEntry {
        hm: "P b a n",
        ext: 50,
        ops: OPS_233,
    },
    GemmiSpaceGroupEntry {
        hm: "P n c b",
        ext: 49,
        ops: OPS_234,
    },
    GemmiSpaceGroupEntry {
        hm: "P n c b",
        ext: 50,
        ops: OPS_235,
    },
    GemmiSpaceGroupEntry {
        hm: "P c n a",
        ext: 49,
        ops: OPS_236,
    },
    GemmiSpaceGroupEntry {
        hm: "P c n a",
        ext: 50,
        ops: OPS_237,
    },
    GemmiSpaceGroupEntry {
        hm: "P m m a",
        ext: 0,
        ops: OPS_238,
    },
    GemmiSpaceGroupEntry {
        hm: "P m m b",
        ext: 0,
        ops: OPS_239,
    },
    GemmiSpaceGroupEntry {
        hm: "P b m m",
        ext: 0,
        ops: OPS_240,
    },
    GemmiSpaceGroupEntry {
        hm: "P c m m",
        ext: 0,
        ops: OPS_241,
    },
    GemmiSpaceGroupEntry {
        hm: "P m c m",
        ext: 0,
        ops: OPS_242,
    },
    GemmiSpaceGroupEntry {
        hm: "P m a m",
        ext: 0,
        ops: OPS_243,
    },
    GemmiSpaceGroupEntry {
        hm: "P n n a",
        ext: 0,
        ops: OPS_244,
    },
    GemmiSpaceGroupEntry {
        hm: "P n n b",
        ext: 0,
        ops: OPS_245,
    },
    GemmiSpaceGroupEntry {
        hm: "P b n n",
        ext: 0,
        ops: OPS_246,
    },
    GemmiSpaceGroupEntry {
        hm: "P c n n",
        ext: 0,
        ops: OPS_247,
    },
    GemmiSpaceGroupEntry {
        hm: "P n c n",
        ext: 0,
        ops: OPS_248,
    },
    GemmiSpaceGroupEntry {
        hm: "P n a n",
        ext: 0,
        ops: OPS_249,
    },
    GemmiSpaceGroupEntry {
        hm: "P m n a",
        ext: 0,
        ops: OPS_250,
    },
    GemmiSpaceGroupEntry {
        hm: "P n m b",
        ext: 0,
        ops: OPS_251,
    },
    GemmiSpaceGroupEntry {
        hm: "P b m n",
        ext: 0,
        ops: OPS_252,
    },
    GemmiSpaceGroupEntry {
        hm: "P c n m",
        ext: 0,
        ops: OPS_253,
    },
    GemmiSpaceGroupEntry {
        hm: "P n c m",
        ext: 0,
        ops: OPS_254,
    },
    GemmiSpaceGroupEntry {
        hm: "P m a n",
        ext: 0,
        ops: OPS_255,
    },
    GemmiSpaceGroupEntry {
        hm: "P c c a",
        ext: 0,
        ops: OPS_256,
    },
    GemmiSpaceGroupEntry {
        hm: "P c c b",
        ext: 0,
        ops: OPS_257,
    },
    GemmiSpaceGroupEntry {
        hm: "P b a a",
        ext: 0,
        ops: OPS_258,
    },
    GemmiSpaceGroupEntry {
        hm: "P c a a",
        ext: 0,
        ops: OPS_259,
    },
    GemmiSpaceGroupEntry {
        hm: "P b c b",
        ext: 0,
        ops: OPS_260,
    },
    GemmiSpaceGroupEntry {
        hm: "P b a b",
        ext: 0,
        ops: OPS_261,
    },
    GemmiSpaceGroupEntry {
        hm: "P b a m",
        ext: 0,
        ops: OPS_262,
    },
    GemmiSpaceGroupEntry {
        hm: "P m c b",
        ext: 0,
        ops: OPS_263,
    },
    GemmiSpaceGroupEntry {
        hm: "P c m a",
        ext: 0,
        ops: OPS_264,
    },
    GemmiSpaceGroupEntry {
        hm: "P c c n",
        ext: 0,
        ops: OPS_265,
    },
    GemmiSpaceGroupEntry {
        hm: "P n a a",
        ext: 0,
        ops: OPS_266,
    },
    GemmiSpaceGroupEntry {
        hm: "P b n b",
        ext: 0,
        ops: OPS_267,
    },
    GemmiSpaceGroupEntry {
        hm: "P b c m",
        ext: 0,
        ops: OPS_268,
    },
    GemmiSpaceGroupEntry {
        hm: "P c a m",
        ext: 0,
        ops: OPS_269,
    },
    GemmiSpaceGroupEntry {
        hm: "P m c a",
        ext: 0,
        ops: OPS_270,
    },
    GemmiSpaceGroupEntry {
        hm: "P m a b",
        ext: 0,
        ops: OPS_271,
    },
    GemmiSpaceGroupEntry {
        hm: "P b m a",
        ext: 0,
        ops: OPS_272,
    },
    GemmiSpaceGroupEntry {
        hm: "P c m b",
        ext: 0,
        ops: OPS_273,
    },
    GemmiSpaceGroupEntry {
        hm: "P n n m",
        ext: 0,
        ops: OPS_274,
    },
    GemmiSpaceGroupEntry {
        hm: "P m n n",
        ext: 0,
        ops: OPS_275,
    },
    GemmiSpaceGroupEntry {
        hm: "P n m n",
        ext: 0,
        ops: OPS_276,
    },
    GemmiSpaceGroupEntry {
        hm: "P m m n",
        ext: 49,
        ops: OPS_277,
    },
    GemmiSpaceGroupEntry {
        hm: "P m m n",
        ext: 50,
        ops: OPS_278,
    },
    GemmiSpaceGroupEntry {
        hm: "P n m m",
        ext: 49,
        ops: OPS_279,
    },
    GemmiSpaceGroupEntry {
        hm: "P n m m",
        ext: 50,
        ops: OPS_280,
    },
    GemmiSpaceGroupEntry {
        hm: "P m n m",
        ext: 49,
        ops: OPS_281,
    },
    GemmiSpaceGroupEntry {
        hm: "P m n m",
        ext: 50,
        ops: OPS_282,
    },
    GemmiSpaceGroupEntry {
        hm: "P b c n",
        ext: 0,
        ops: OPS_283,
    },
    GemmiSpaceGroupEntry {
        hm: "P c a n",
        ext: 0,
        ops: OPS_284,
    },
    GemmiSpaceGroupEntry {
        hm: "P n c a",
        ext: 0,
        ops: OPS_285,
    },
    GemmiSpaceGroupEntry {
        hm: "P n a b",
        ext: 0,
        ops: OPS_286,
    },
    GemmiSpaceGroupEntry {
        hm: "P b n a",
        ext: 0,
        ops: OPS_287,
    },
    GemmiSpaceGroupEntry {
        hm: "P c n b",
        ext: 0,
        ops: OPS_288,
    },
    GemmiSpaceGroupEntry {
        hm: "P b c a",
        ext: 0,
        ops: OPS_289,
    },
    GemmiSpaceGroupEntry {
        hm: "P c a b",
        ext: 0,
        ops: OPS_290,
    },
    GemmiSpaceGroupEntry {
        hm: "P n m a",
        ext: 0,
        ops: OPS_291,
    },
    GemmiSpaceGroupEntry {
        hm: "P m n b",
        ext: 0,
        ops: OPS_292,
    },
    GemmiSpaceGroupEntry {
        hm: "P b n m",
        ext: 0,
        ops: OPS_293,
    },
    GemmiSpaceGroupEntry {
        hm: "P c m n",
        ext: 0,
        ops: OPS_294,
    },
    GemmiSpaceGroupEntry {
        hm: "P m c n",
        ext: 0,
        ops: OPS_295,
    },
    GemmiSpaceGroupEntry {
        hm: "P n a m",
        ext: 0,
        ops: OPS_296,
    },
    GemmiSpaceGroupEntry {
        hm: "C m c m",
        ext: 0,
        ops: OPS_297,
    },
    GemmiSpaceGroupEntry {
        hm: "C c m m",
        ext: 0,
        ops: OPS_298,
    },
    GemmiSpaceGroupEntry {
        hm: "A m m a",
        ext: 0,
        ops: OPS_299,
    },
    GemmiSpaceGroupEntry {
        hm: "A m a m",
        ext: 0,
        ops: OPS_300,
    },
    GemmiSpaceGroupEntry {
        hm: "B b m m",
        ext: 0,
        ops: OPS_301,
    },
    GemmiSpaceGroupEntry {
        hm: "B m m b",
        ext: 0,
        ops: OPS_302,
    },
    GemmiSpaceGroupEntry {
        hm: "C m c a",
        ext: 0,
        ops: OPS_303,
    },
    GemmiSpaceGroupEntry {
        hm: "C c m b",
        ext: 0,
        ops: OPS_304,
    },
    GemmiSpaceGroupEntry {
        hm: "A b m a",
        ext: 0,
        ops: OPS_305,
    },
    GemmiSpaceGroupEntry {
        hm: "A c a m",
        ext: 0,
        ops: OPS_306,
    },
    GemmiSpaceGroupEntry {
        hm: "B b c m",
        ext: 0,
        ops: OPS_307,
    },
    GemmiSpaceGroupEntry {
        hm: "B m a b",
        ext: 0,
        ops: OPS_308,
    },
    GemmiSpaceGroupEntry {
        hm: "C m m m",
        ext: 0,
        ops: OPS_309,
    },
    GemmiSpaceGroupEntry {
        hm: "A m m m",
        ext: 0,
        ops: OPS_310,
    },
    GemmiSpaceGroupEntry {
        hm: "B m m m",
        ext: 0,
        ops: OPS_311,
    },
    GemmiSpaceGroupEntry {
        hm: "C c c m",
        ext: 0,
        ops: OPS_312,
    },
    GemmiSpaceGroupEntry {
        hm: "A m a a",
        ext: 0,
        ops: OPS_313,
    },
    GemmiSpaceGroupEntry {
        hm: "B b m b",
        ext: 0,
        ops: OPS_314,
    },
    GemmiSpaceGroupEntry {
        hm: "C m m a",
        ext: 0,
        ops: OPS_315,
    },
    GemmiSpaceGroupEntry {
        hm: "C m m b",
        ext: 0,
        ops: OPS_316,
    },
    GemmiSpaceGroupEntry {
        hm: "A b m m",
        ext: 0,
        ops: OPS_317,
    },
    GemmiSpaceGroupEntry {
        hm: "A c m m",
        ext: 0,
        ops: OPS_318,
    },
    GemmiSpaceGroupEntry {
        hm: "B m c m",
        ext: 0,
        ops: OPS_319,
    },
    GemmiSpaceGroupEntry {
        hm: "B m a m",
        ext: 0,
        ops: OPS_320,
    },
    GemmiSpaceGroupEntry {
        hm: "C c c a",
        ext: 49,
        ops: OPS_321,
    },
    GemmiSpaceGroupEntry {
        hm: "C c c a",
        ext: 50,
        ops: OPS_322,
    },
    GemmiSpaceGroupEntry {
        hm: "C c c b",
        ext: 49,
        ops: OPS_323,
    },
    GemmiSpaceGroupEntry {
        hm: "C c c b",
        ext: 50,
        ops: OPS_324,
    },
    GemmiSpaceGroupEntry {
        hm: "A b a a",
        ext: 49,
        ops: OPS_325,
    },
    GemmiSpaceGroupEntry {
        hm: "A b a a",
        ext: 50,
        ops: OPS_326,
    },
    GemmiSpaceGroupEntry {
        hm: "A c a a",
        ext: 49,
        ops: OPS_327,
    },
    GemmiSpaceGroupEntry {
        hm: "A c a a",
        ext: 50,
        ops: OPS_328,
    },
    GemmiSpaceGroupEntry {
        hm: "B b c b",
        ext: 49,
        ops: OPS_329,
    },
    GemmiSpaceGroupEntry {
        hm: "B b c b",
        ext: 50,
        ops: OPS_330,
    },
    GemmiSpaceGroupEntry {
        hm: "B b a b",
        ext: 49,
        ops: OPS_331,
    },
    GemmiSpaceGroupEntry {
        hm: "B b a b",
        ext: 50,
        ops: OPS_332,
    },
    GemmiSpaceGroupEntry {
        hm: "F m m m",
        ext: 0,
        ops: OPS_333,
    },
    GemmiSpaceGroupEntry {
        hm: "F d d d",
        ext: 49,
        ops: OPS_334,
    },
    GemmiSpaceGroupEntry {
        hm: "F d d d",
        ext: 50,
        ops: OPS_335,
    },
    GemmiSpaceGroupEntry {
        hm: "I m m m",
        ext: 0,
        ops: OPS_336,
    },
    GemmiSpaceGroupEntry {
        hm: "I b a m",
        ext: 0,
        ops: OPS_337,
    },
    GemmiSpaceGroupEntry {
        hm: "I m c b",
        ext: 0,
        ops: OPS_338,
    },
    GemmiSpaceGroupEntry {
        hm: "I c m a",
        ext: 0,
        ops: OPS_339,
    },
    GemmiSpaceGroupEntry {
        hm: "I b c a",
        ext: 0,
        ops: OPS_340,
    },
    GemmiSpaceGroupEntry {
        hm: "I c a b",
        ext: 0,
        ops: OPS_341,
    },
    GemmiSpaceGroupEntry {
        hm: "I m m a",
        ext: 0,
        ops: OPS_342,
    },
    GemmiSpaceGroupEntry {
        hm: "I m m b",
        ext: 0,
        ops: OPS_343,
    },
    GemmiSpaceGroupEntry {
        hm: "I b m m",
        ext: 0,
        ops: OPS_344,
    },
    GemmiSpaceGroupEntry {
        hm: "I c m m",
        ext: 0,
        ops: OPS_345,
    },
    GemmiSpaceGroupEntry {
        hm: "I m c m",
        ext: 0,
        ops: OPS_346,
    },
    GemmiSpaceGroupEntry {
        hm: "I m a m",
        ext: 0,
        ops: OPS_347,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4",
        ext: 0,
        ops: OPS_348,
    },
    GemmiSpaceGroupEntry {
        hm: "P 41",
        ext: 0,
        ops: OPS_349,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42",
        ext: 0,
        ops: OPS_350,
    },
    GemmiSpaceGroupEntry {
        hm: "P 43",
        ext: 0,
        ops: OPS_351,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4",
        ext: 0,
        ops: OPS_352,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41",
        ext: 0,
        ops: OPS_353,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4",
        ext: 0,
        ops: OPS_354,
    },
    GemmiSpaceGroupEntry {
        hm: "I -4",
        ext: 0,
        ops: OPS_355,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/m",
        ext: 0,
        ops: OPS_356,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/m",
        ext: 0,
        ops: OPS_357,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n",
        ext: 49,
        ops: OPS_358,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n",
        ext: 50,
        ops: OPS_359,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n",
        ext: 49,
        ops: OPS_360,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n",
        ext: 50,
        ops: OPS_361,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4/m",
        ext: 0,
        ops: OPS_362,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41/a",
        ext: 49,
        ops: OPS_363,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41/a",
        ext: 50,
        ops: OPS_364,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4 2 2",
        ext: 0,
        ops: OPS_365,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4 21 2",
        ext: 0,
        ops: OPS_366,
    },
    GemmiSpaceGroupEntry {
        hm: "P 41 2 2",
        ext: 0,
        ops: OPS_367,
    },
    GemmiSpaceGroupEntry {
        hm: "P 41 21 2",
        ext: 0,
        ops: OPS_368,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 2 2",
        ext: 0,
        ops: OPS_369,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 21 2",
        ext: 0,
        ops: OPS_370,
    },
    GemmiSpaceGroupEntry {
        hm: "P 43 2 2",
        ext: 0,
        ops: OPS_371,
    },
    GemmiSpaceGroupEntry {
        hm: "P 43 21 2",
        ext: 0,
        ops: OPS_372,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4 2 2",
        ext: 0,
        ops: OPS_373,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41 2 2",
        ext: 0,
        ops: OPS_374,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4 m m",
        ext: 0,
        ops: OPS_375,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4 b m",
        ext: 0,
        ops: OPS_376,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 c m",
        ext: 0,
        ops: OPS_377,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 n m",
        ext: 0,
        ops: OPS_378,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4 c c",
        ext: 0,
        ops: OPS_379,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4 n c",
        ext: 0,
        ops: OPS_380,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 m c",
        ext: 0,
        ops: OPS_381,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 b c",
        ext: 0,
        ops: OPS_382,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4 m m",
        ext: 0,
        ops: OPS_383,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4 c m",
        ext: 0,
        ops: OPS_384,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41 m d",
        ext: 0,
        ops: OPS_385,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41 c d",
        ext: 0,
        ops: OPS_386,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 2 m",
        ext: 0,
        ops: OPS_387,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 2 c",
        ext: 0,
        ops: OPS_388,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 21 m",
        ext: 0,
        ops: OPS_389,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 21 c",
        ext: 0,
        ops: OPS_390,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 m 2",
        ext: 0,
        ops: OPS_391,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 c 2",
        ext: 0,
        ops: OPS_392,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 b 2",
        ext: 0,
        ops: OPS_393,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 n 2",
        ext: 0,
        ops: OPS_394,
    },
    GemmiSpaceGroupEntry {
        hm: "I -4 m 2",
        ext: 0,
        ops: OPS_395,
    },
    GemmiSpaceGroupEntry {
        hm: "I -4 c 2",
        ext: 0,
        ops: OPS_396,
    },
    GemmiSpaceGroupEntry {
        hm: "I -4 2 m",
        ext: 0,
        ops: OPS_397,
    },
    GemmiSpaceGroupEntry {
        hm: "I -4 2 d",
        ext: 0,
        ops: OPS_398,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/m m m",
        ext: 0,
        ops: OPS_399,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/m c c",
        ext: 0,
        ops: OPS_400,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n b m",
        ext: 49,
        ops: OPS_401,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n b m",
        ext: 50,
        ops: OPS_402,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n n c",
        ext: 49,
        ops: OPS_403,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n n c",
        ext: 50,
        ops: OPS_404,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/m b m",
        ext: 0,
        ops: OPS_405,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/m n c",
        ext: 0,
        ops: OPS_406,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n m m",
        ext: 49,
        ops: OPS_407,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n m m",
        ext: 50,
        ops: OPS_408,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n c c",
        ext: 49,
        ops: OPS_409,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4/n c c",
        ext: 50,
        ops: OPS_410,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/m m c",
        ext: 0,
        ops: OPS_411,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/m c m",
        ext: 0,
        ops: OPS_412,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n b c",
        ext: 49,
        ops: OPS_413,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n b c",
        ext: 50,
        ops: OPS_414,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n n m",
        ext: 49,
        ops: OPS_415,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n n m",
        ext: 50,
        ops: OPS_416,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/m b c",
        ext: 0,
        ops: OPS_417,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/m n m",
        ext: 0,
        ops: OPS_418,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n m c",
        ext: 49,
        ops: OPS_419,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n m c",
        ext: 50,
        ops: OPS_420,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n c m",
        ext: 49,
        ops: OPS_421,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42/n c m",
        ext: 50,
        ops: OPS_422,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4/m m m",
        ext: 0,
        ops: OPS_423,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4/m c m",
        ext: 0,
        ops: OPS_424,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41/a m d",
        ext: 49,
        ops: OPS_425,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41/a m d",
        ext: 50,
        ops: OPS_426,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41/a c d",
        ext: 49,
        ops: OPS_427,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41/a c d",
        ext: 50,
        ops: OPS_428,
    },
    GemmiSpaceGroupEntry {
        hm: "P 3",
        ext: 0,
        ops: OPS_429,
    },
    GemmiSpaceGroupEntry {
        hm: "P 31",
        ext: 0,
        ops: OPS_430,
    },
    GemmiSpaceGroupEntry {
        hm: "P 32",
        ext: 0,
        ops: OPS_431,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3",
        ext: 72,
        ops: OPS_432,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3",
        ext: 82,
        ops: OPS_433,
    },
    GemmiSpaceGroupEntry {
        hm: "P -3",
        ext: 0,
        ops: OPS_434,
    },
    GemmiSpaceGroupEntry {
        hm: "R -3",
        ext: 72,
        ops: OPS_435,
    },
    GemmiSpaceGroupEntry {
        hm: "R -3",
        ext: 82,
        ops: OPS_436,
    },
    GemmiSpaceGroupEntry {
        hm: "P 3 1 2",
        ext: 0,
        ops: OPS_437,
    },
    GemmiSpaceGroupEntry {
        hm: "P 3 2 1",
        ext: 0,
        ops: OPS_438,
    },
    GemmiSpaceGroupEntry {
        hm: "P 31 1 2",
        ext: 0,
        ops: OPS_439,
    },
    GemmiSpaceGroupEntry {
        hm: "P 31 2 1",
        ext: 0,
        ops: OPS_440,
    },
    GemmiSpaceGroupEntry {
        hm: "P 32 1 2",
        ext: 0,
        ops: OPS_441,
    },
    GemmiSpaceGroupEntry {
        hm: "P 32 2 1",
        ext: 0,
        ops: OPS_442,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3 2",
        ext: 72,
        ops: OPS_443,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3 2",
        ext: 82,
        ops: OPS_444,
    },
    GemmiSpaceGroupEntry {
        hm: "P 3 m 1",
        ext: 0,
        ops: OPS_445,
    },
    GemmiSpaceGroupEntry {
        hm: "P 3 1 m",
        ext: 0,
        ops: OPS_446,
    },
    GemmiSpaceGroupEntry {
        hm: "P 3 c 1",
        ext: 0,
        ops: OPS_447,
    },
    GemmiSpaceGroupEntry {
        hm: "P 3 1 c",
        ext: 0,
        ops: OPS_448,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3 m",
        ext: 72,
        ops: OPS_449,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3 m",
        ext: 82,
        ops: OPS_450,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3 c",
        ext: 72,
        ops: OPS_451,
    },
    GemmiSpaceGroupEntry {
        hm: "R 3 c",
        ext: 82,
        ops: OPS_452,
    },
    GemmiSpaceGroupEntry {
        hm: "P -3 1 m",
        ext: 0,
        ops: OPS_453,
    },
    GemmiSpaceGroupEntry {
        hm: "P -3 1 c",
        ext: 0,
        ops: OPS_454,
    },
    GemmiSpaceGroupEntry {
        hm: "P -3 m 1",
        ext: 0,
        ops: OPS_455,
    },
    GemmiSpaceGroupEntry {
        hm: "P -3 c 1",
        ext: 0,
        ops: OPS_456,
    },
    GemmiSpaceGroupEntry {
        hm: "R -3 m",
        ext: 72,
        ops: OPS_457,
    },
    GemmiSpaceGroupEntry {
        hm: "R -3 m",
        ext: 82,
        ops: OPS_458,
    },
    GemmiSpaceGroupEntry {
        hm: "R -3 c",
        ext: 72,
        ops: OPS_459,
    },
    GemmiSpaceGroupEntry {
        hm: "R -3 c",
        ext: 82,
        ops: OPS_460,
    },
    GemmiSpaceGroupEntry {
        hm: "P 6",
        ext: 0,
        ops: OPS_461,
    },
    GemmiSpaceGroupEntry {
        hm: "P 61",
        ext: 0,
        ops: OPS_462,
    },
    GemmiSpaceGroupEntry {
        hm: "P 65",
        ext: 0,
        ops: OPS_463,
    },
    GemmiSpaceGroupEntry {
        hm: "P 62",
        ext: 0,
        ops: OPS_464,
    },
    GemmiSpaceGroupEntry {
        hm: "P 64",
        ext: 0,
        ops: OPS_465,
    },
    GemmiSpaceGroupEntry {
        hm: "P 63",
        ext: 0,
        ops: OPS_466,
    },
    GemmiSpaceGroupEntry {
        hm: "P -6",
        ext: 0,
        ops: OPS_467,
    },
    GemmiSpaceGroupEntry {
        hm: "P 6/m",
        ext: 0,
        ops: OPS_468,
    },
    GemmiSpaceGroupEntry {
        hm: "P 63/m",
        ext: 0,
        ops: OPS_469,
    },
    GemmiSpaceGroupEntry {
        hm: "P 6 2 2",
        ext: 0,
        ops: OPS_470,
    },
    GemmiSpaceGroupEntry {
        hm: "P 61 2 2",
        ext: 0,
        ops: OPS_471,
    },
    GemmiSpaceGroupEntry {
        hm: "P 65 2 2",
        ext: 0,
        ops: OPS_472,
    },
    GemmiSpaceGroupEntry {
        hm: "P 62 2 2",
        ext: 0,
        ops: OPS_473,
    },
    GemmiSpaceGroupEntry {
        hm: "P 64 2 2",
        ext: 0,
        ops: OPS_474,
    },
    GemmiSpaceGroupEntry {
        hm: "P 63 2 2",
        ext: 0,
        ops: OPS_475,
    },
    GemmiSpaceGroupEntry {
        hm: "P 6 m m",
        ext: 0,
        ops: OPS_476,
    },
    GemmiSpaceGroupEntry {
        hm: "P 6 c c",
        ext: 0,
        ops: OPS_477,
    },
    GemmiSpaceGroupEntry {
        hm: "P 63 c m",
        ext: 0,
        ops: OPS_478,
    },
    GemmiSpaceGroupEntry {
        hm: "P 63 m c",
        ext: 0,
        ops: OPS_479,
    },
    GemmiSpaceGroupEntry {
        hm: "P -6 m 2",
        ext: 0,
        ops: OPS_480,
    },
    GemmiSpaceGroupEntry {
        hm: "P -6 c 2",
        ext: 0,
        ops: OPS_481,
    },
    GemmiSpaceGroupEntry {
        hm: "P -6 2 m",
        ext: 0,
        ops: OPS_482,
    },
    GemmiSpaceGroupEntry {
        hm: "P -6 2 c",
        ext: 0,
        ops: OPS_483,
    },
    GemmiSpaceGroupEntry {
        hm: "P 6/m m m",
        ext: 0,
        ops: OPS_484,
    },
    GemmiSpaceGroupEntry {
        hm: "P 6/m c c",
        ext: 0,
        ops: OPS_485,
    },
    GemmiSpaceGroupEntry {
        hm: "P 63/m c m",
        ext: 0,
        ops: OPS_486,
    },
    GemmiSpaceGroupEntry {
        hm: "P 63/m m c",
        ext: 0,
        ops: OPS_487,
    },
    GemmiSpaceGroupEntry {
        hm: "P 2 3",
        ext: 0,
        ops: OPS_488,
    },
    GemmiSpaceGroupEntry {
        hm: "F 2 3",
        ext: 0,
        ops: OPS_489,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 3",
        ext: 0,
        ops: OPS_490,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21 3",
        ext: 0,
        ops: OPS_491,
    },
    GemmiSpaceGroupEntry {
        hm: "I 21 3",
        ext: 0,
        ops: OPS_492,
    },
    GemmiSpaceGroupEntry {
        hm: "P m -3",
        ext: 0,
        ops: OPS_493,
    },
    GemmiSpaceGroupEntry {
        hm: "P n -3",
        ext: 49,
        ops: OPS_494,
    },
    GemmiSpaceGroupEntry {
        hm: "P n -3",
        ext: 50,
        ops: OPS_495,
    },
    GemmiSpaceGroupEntry {
        hm: "F m -3",
        ext: 0,
        ops: OPS_496,
    },
    GemmiSpaceGroupEntry {
        hm: "F d -3",
        ext: 49,
        ops: OPS_497,
    },
    GemmiSpaceGroupEntry {
        hm: "F d -3",
        ext: 50,
        ops: OPS_498,
    },
    GemmiSpaceGroupEntry {
        hm: "I m -3",
        ext: 0,
        ops: OPS_499,
    },
    GemmiSpaceGroupEntry {
        hm: "P a -3",
        ext: 0,
        ops: OPS_500,
    },
    GemmiSpaceGroupEntry {
        hm: "I a -3",
        ext: 0,
        ops: OPS_501,
    },
    GemmiSpaceGroupEntry {
        hm: "P 4 3 2",
        ext: 0,
        ops: OPS_502,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 3 2",
        ext: 0,
        ops: OPS_503,
    },
    GemmiSpaceGroupEntry {
        hm: "F 4 3 2",
        ext: 0,
        ops: OPS_504,
    },
    GemmiSpaceGroupEntry {
        hm: "F 41 3 2",
        ext: 0,
        ops: OPS_505,
    },
    GemmiSpaceGroupEntry {
        hm: "I 4 3 2",
        ext: 0,
        ops: OPS_506,
    },
    GemmiSpaceGroupEntry {
        hm: "P 43 3 2",
        ext: 0,
        ops: OPS_507,
    },
    GemmiSpaceGroupEntry {
        hm: "P 41 3 2",
        ext: 0,
        ops: OPS_508,
    },
    GemmiSpaceGroupEntry {
        hm: "I 41 3 2",
        ext: 0,
        ops: OPS_509,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 3 m",
        ext: 0,
        ops: OPS_510,
    },
    GemmiSpaceGroupEntry {
        hm: "F -4 3 m",
        ext: 0,
        ops: OPS_511,
    },
    GemmiSpaceGroupEntry {
        hm: "I -4 3 m",
        ext: 0,
        ops: OPS_512,
    },
    GemmiSpaceGroupEntry {
        hm: "P -4 3 n",
        ext: 0,
        ops: OPS_513,
    },
    GemmiSpaceGroupEntry {
        hm: "F -4 3 c",
        ext: 0,
        ops: OPS_514,
    },
    GemmiSpaceGroupEntry {
        hm: "I -4 3 d",
        ext: 0,
        ops: OPS_515,
    },
    GemmiSpaceGroupEntry {
        hm: "P m -3 m",
        ext: 0,
        ops: OPS_516,
    },
    GemmiSpaceGroupEntry {
        hm: "P n -3 n",
        ext: 49,
        ops: OPS_517,
    },
    GemmiSpaceGroupEntry {
        hm: "P n -3 n",
        ext: 50,
        ops: OPS_518,
    },
    GemmiSpaceGroupEntry {
        hm: "P m -3 n",
        ext: 0,
        ops: OPS_519,
    },
    GemmiSpaceGroupEntry {
        hm: "P n -3 m",
        ext: 49,
        ops: OPS_520,
    },
    GemmiSpaceGroupEntry {
        hm: "P n -3 m",
        ext: 50,
        ops: OPS_521,
    },
    GemmiSpaceGroupEntry {
        hm: "F m -3 m",
        ext: 0,
        ops: OPS_522,
    },
    GemmiSpaceGroupEntry {
        hm: "F m -3 c",
        ext: 0,
        ops: OPS_523,
    },
    GemmiSpaceGroupEntry {
        hm: "F d -3 m",
        ext: 49,
        ops: OPS_524,
    },
    GemmiSpaceGroupEntry {
        hm: "F d -3 m",
        ext: 50,
        ops: OPS_525,
    },
    GemmiSpaceGroupEntry {
        hm: "F d -3 c",
        ext: 49,
        ops: OPS_526,
    },
    GemmiSpaceGroupEntry {
        hm: "F d -3 c",
        ext: 50,
        ops: OPS_527,
    },
    GemmiSpaceGroupEntry {
        hm: "I m -3 m",
        ext: 0,
        ops: OPS_528,
    },
    GemmiSpaceGroupEntry {
        hm: "I a -3 d",
        ext: 0,
        ops: OPS_529,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1 21 1",
        ext: 0,
        ops: OPS_530,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 21 1",
        ext: 0,
        ops: OPS_531,
    },
    GemmiSpaceGroupEntry {
        hm: "P 21212(a)",
        ext: 0,
        ops: OPS_532,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 2 21a)",
        ext: 0,
        ops: OPS_533,
    },
    GemmiSpaceGroupEntry {
        hm: "C 2 2 2a",
        ext: 0,
        ops: OPS_534,
    },
    GemmiSpaceGroupEntry {
        hm: "F 2 2 2a",
        ext: 0,
        ops: OPS_535,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 2 2a",
        ext: 0,
        ops: OPS_536,
    },
    GemmiSpaceGroupEntry {
        hm: "P 42 21 2a",
        ext: 0,
        ops: OPS_537,
    },
    GemmiSpaceGroupEntry {
        hm: "I 2 3a",
        ext: 0,
        ops: OPS_538,
    },
    GemmiSpaceGroupEntry {
        hm: "A 1",
        ext: 0,
        ops: OPS_539,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1",
        ext: 0,
        ops: OPS_540,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1",
        ext: 0,
        ops: OPS_541,
    },
    GemmiSpaceGroupEntry {
        hm: "F 1",
        ext: 0,
        ops: OPS_542,
    },
    GemmiSpaceGroupEntry {
        hm: "I 1",
        ext: 0,
        ops: OPS_543,
    },
    GemmiSpaceGroupEntry {
        hm: "A -1",
        ext: 0,
        ops: OPS_544,
    },
    GemmiSpaceGroupEntry {
        hm: "B -1",
        ext: 0,
        ops: OPS_545,
    },
    GemmiSpaceGroupEntry {
        hm: "C -1",
        ext: 0,
        ops: OPS_546,
    },
    GemmiSpaceGroupEntry {
        hm: "F -1",
        ext: 0,
        ops: OPS_547,
    },
    GemmiSpaceGroupEntry {
        hm: "I -1",
        ext: 0,
        ops: OPS_548,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 2 1",
        ext: 0,
        ops: OPS_549,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 1 2",
        ext: 0,
        ops: OPS_550,
    },
    GemmiSpaceGroupEntry {
        hm: "B 1 21 1",
        ext: 0,
        ops: OPS_551,
    },
    GemmiSpaceGroupEntry {
        hm: "C 1 1 21",
        ext: 0,
        ops: OPS_552,
    },
    GemmiSpaceGroupEntry {
        hm: "F 1 2 1",
        ext: 0,
        ops: OPS_553,
    },
    GemmiSpaceGroupEntry {
        hm: "F 1 m 1",
        ext: 0,
        ops: OPS_554,
    },
    GemmiSpaceGroupEntry {
        hm: "F 1 d 1",
        ext: 0,
        ops: OPS_555,
    },
    GemmiSpaceGroupEntry {
        hm: "F 1 2/m 1",
        ext: 0,
        ops: OPS_556,
    },
    GemmiSpaceGroupEntry {
        hm: "A b a m",
        ext: 0,
        ops: OPS_557,
    },
    GemmiSpaceGroupEntry {
        hm: "C 4 2 2",
        ext: 0,
        ops: OPS_558,
    },
    GemmiSpaceGroupEntry {
        hm: "C 4 2 21",
        ext: 0,
        ops: OPS_559,
    },
    GemmiSpaceGroupEntry {
        hm: "F 4 2 2",
        ext: 0,
        ops: OPS_560,
    },
    GemmiSpaceGroupEntry {
        hm: "C -4 2 m",
        ext: 0,
        ops: OPS_561,
    },
    GemmiSpaceGroupEntry {
        hm: "C -4 2 b",
        ext: 0,
        ops: OPS_562,
    },
    GemmiSpaceGroupEntry {
        hm: "F 4/m m m",
        ext: 0,
        ops: OPS_563,
    },
];

pub(crate) const GEMMI_ALT_NAMES: &[GemmiAltNameEntry] = &[
    GemmiAltNameEntry {
        hm: "A e m 2",
        ext: 0,
        pos: 190,
    },
    GemmiAltNameEntry {
        hm: "B m e 2",
        ext: 0,
        pos: 191,
    },
    GemmiAltNameEntry {
        hm: "B 2 e m",
        ext: 0,
        pos: 192,
    },
    GemmiAltNameEntry {
        hm: "C 2 m e",
        ext: 0,
        pos: 193,
    },
    GemmiAltNameEntry {
        hm: "C m 2 e",
        ext: 0,
        pos: 194,
    },
    GemmiAltNameEntry {
        hm: "A e 2 m",
        ext: 0,
        pos: 195,
    },
    GemmiAltNameEntry {
        hm: "A e a 2",
        ext: 0,
        pos: 202,
    },
    GemmiAltNameEntry {
        hm: "B b e 2",
        ext: 0,
        pos: 203,
    },
    GemmiAltNameEntry {
        hm: "B 2 e b",
        ext: 0,
        pos: 204,
    },
    GemmiAltNameEntry {
        hm: "C 2 c e",
        ext: 0,
        pos: 205,
    },
    GemmiAltNameEntry {
        hm: "C c 2 e",
        ext: 0,
        pos: 206,
    },
    GemmiAltNameEntry {
        hm: "A e 2 a",
        ext: 0,
        pos: 207,
    },
    GemmiAltNameEntry {
        hm: "C m c e",
        ext: 0,
        pos: 303,
    },
    GemmiAltNameEntry {
        hm: "C c m e",
        ext: 0,
        pos: 304,
    },
    GemmiAltNameEntry {
        hm: "A e m a",
        ext: 0,
        pos: 305,
    },
    GemmiAltNameEntry {
        hm: "A e a m",
        ext: 0,
        pos: 306,
    },
    GemmiAltNameEntry {
        hm: "B b e m",
        ext: 0,
        pos: 307,
    },
    GemmiAltNameEntry {
        hm: "B m e b",
        ext: 0,
        pos: 308,
    },
    GemmiAltNameEntry {
        hm: "C m m e",
        ext: 0,
        pos: 315,
    },
    GemmiAltNameEntry {
        hm: "A e m m",
        ext: 0,
        pos: 317,
    },
    GemmiAltNameEntry {
        hm: "B m e m",
        ext: 0,
        pos: 319,
    },
    GemmiAltNameEntry {
        hm: "C c c e",
        ext: 49,
        pos: 321,
    },
    GemmiAltNameEntry {
        hm: "C c c e",
        ext: 50,
        pos: 322,
    },
    GemmiAltNameEntry {
        hm: "A e a a",
        ext: 49,
        pos: 325,
    },
    GemmiAltNameEntry {
        hm: "A e a a",
        ext: 50,
        pos: 326,
    },
    GemmiAltNameEntry {
        hm: "B b e b",
        ext: 49,
        pos: 329,
    },
    GemmiAltNameEntry {
        hm: "B b e b",
        ext: 50,
        pos: 330,
    },
    GemmiAltNameEntry {
        hm: "P 21 21 2a",
        ext: 0,
        pos: 532,
    },
];
