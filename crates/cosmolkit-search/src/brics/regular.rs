//! Regular r-BRICS rules (not the reBRICS post-pass).
//  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//  Copyright (c) 2022, IBM LLC, Leili Zhang, Vasu Rao, Wendy Cornell
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//     * Neither the name of International Business Machine
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use super::RuleSet;

// rBRICS_public.py, environs and reactionDefs at 03938c2.
pub(super) const REGULAR: RuleSet = RuleSet {
    environments: &[
        ("L1", "[C;D3]([#0,#6,#7,#8])(=O)"),
        ("L3", "[O;D2]-;!@[#0,#6,#1]"),
        ("L4", "[C;!D1;!$(C=*)]-;!@[#6]"),
        (
            "L5",
            "[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]",
        ),
        ("L51", "[N;!R;!D1;$(N(!@[N,O]))]"),
        ("L6", "[C;D3;!R](=O)-;!@[#0,#6,#7,#8]"),
        ("L7a", "[C;D2,D3]-[#6]"),
        ("L7b", "[C;D2,D3]-[#6]"),
        ("#L8", "[C;!R;!D1]-;!@[#6]"),
        ("##L8", "[C;!R;!D1;!$(C!-*)]"),
        ("L8", "[C;!R;!D1;!$(C!-*);!$(C([H])([H])([H]))]"),
        ("L81", "[C;!R;!D1;$(C(-[C,N,O,S])(=[N,S]))]"),
        (
            "L9",
            "[RN,n;+0;$([RN,n](@[RC,RN,RO,RS,c,n,o,s])@[RC,RN,RO,RS,c,n,o,s])]",
        ),
        ("#L9", "[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]"),
        ("L10", "[N;R;$(N(@C(=O))@[C,N,O,S])]"),
        ("L11", "[S;D2](-;!@[#0,#6])"),
        ("L12", "[S;D4]([#6,#0])(=O)(=O)"),
        ("L12b", "[S;D4;!R](!@O)(!@O)"),
        ("L13", "[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]"),
        ("L14", "[c;$(c(:[c,n,o,s]):[n,o,s])]"),
        ("L14b", "[RC;$([RC](@[RC,RN,RO,RS])@[RN,RO,RS])]"),
        ("L15", "[C;$(C(-;@C)-;@C)]"),
        ("L16", "[c;$(c(:c):c)]"),
        ("L16b", "[RC;$([RC](@[RC])@[RC])]"),
        ("L17", "[C](-C)(-C)(-C)"),
        ("L18", "[R!#1;x3]"),
        ("L19", "[R!#1;x2]"),
        ("L182", "[R!#1;x3]"),
        ("L192", "[R!#1;x2]"),
        ("L20", "[CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]"),
        ("L21", "[CH2][CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]"),
        (
            "L22",
            "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]",
        ),
        (
            "L23",
            "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3,RC,c,$(C(~[!#6]))]",
        ),
        ("L30", "[C;D2]([#0,#6,#7,#8,#16])(#[N,C])"),
    ],
    groups: &[
        &[("L1", "L3", "-"), ("L1", "L5", "-"), ("L1", "L10", "-")],
        &[
            ("L30", "L30", "-"),
            ("L30", "L4", "-"),
            ("L30", "L5", "-"),
            ("L30", "L51", "-"),
            ("L30", "L6", "-"),
            ("L30", "L81", "-"),
            ("L30", "L9", "-"),
            ("L30", "L10", "-"),
            ("L30", "L11", "-"),
            ("L30", "L12", "-"),
            ("L30", "L12b", "-"),
            ("L30", "L13", "-"),
            ("L30", "L14", "-"),
            ("L30", "L14b", "-"),
            ("L30", "L15", "-"),
            ("L30", "L16", "-"),
            ("L30", "L16b", "-"),
        ],
        &[
            ("L3", "L4", "-"),
            ("L3", "L13", "-"),
            ("L3", "L14", "-"),
            ("L3", "L15", "-"),
            ("L3", "L16", "-"),
        ],
        &[("L4", "L5", "-"), ("L4", "L11", "-")],
        &[
            ("L5", "L12", "-"),
            ("L5", "L14", "-"),
            ("L5", "L16", "-"),
            ("L5", "L13", "-"),
            ("L5", "L15", "-"),
        ],
        &[
            ("L51", "L1", "-"),
            ("L51", "L4", "-"),
            ("L51", "L12", "-"),
            ("L51", "L12b", "-"),
            ("L51", "L14", "-"),
            ("L51", "L16", "-"),
            ("L51", "L13", "-"),
            ("L51", "L15", "-"),
        ],
        &[
            ("L6", "L13", "-"),
            ("L6", "L14", "-"),
            ("L6", "L15", "-"),
            ("L6", "L16", "-"),
        ],
        &[("L7a", "L7b", "=")],
        &[
            ("L8", "L9", "-"),
            ("L8", "L10", "-"),
            ("L8", "L13", "-"),
            ("L8", "L14", "-"),
            ("L8", "L15", "-"),
            ("L8", "L16", "-"),
        ],
        &[
            ("L81", "L8", "-"),
            ("L81", "L9", "-"),
            ("L81", "L10", "-"),
            ("L81", "L13", "-"),
            ("L81", "L14", "-"),
            ("L81", "L15", "-"),
            ("L81", "L16", "-"),
        ],
        &[
            ("L9", "L13", "-"),
            ("L9", "L14", "-"),
            ("L9", "L15", "-"),
            ("L9", "L16", "-"),
        ],
        &[
            ("L10", "L13", "-"),
            ("L10", "L14", "-"),
            ("L10", "L15", "-"),
            ("L10", "L16", "-"),
        ],
        &[
            ("L11", "L13", "-"),
            ("L11", "L14", "-"),
            ("L11", "L15", "-"),
            ("L11", "L16", "-"),
        ],
        &[
            ("L12b", "L12b", "-"),
            ("L12b", "L5", "-"),
            ("L12b", "L4", "-"),
            ("L12b", "L13", "-"),
            ("L12b", "L14", "-"),
            ("L12b", "L15", "-"),
            ("L12b", "L16", "-"),
        ],
        &[
            ("L13", "L14", "-;@,!@"),
            ("L13", "L15", "-;!@"),
            ("L13", "L16", "-;@,!@"),
        ],
        &[
            ("L14", "L14", "-;@,!@"),
            ("L14", "L15", "-;@,!@"),
            ("L14", "L16", "-;@,!@"),
        ],
        &[
            ("L3", "L14b", "-"),
            ("L5", "L14b", "-"),
            ("L51", "L14b", "-"),
            ("L6", "L14b", "-"),
            ("L8", "L14b", "-"),
            ("L81", "L14b", "-"),
            ("L9", "L14b", "-"),
            ("L10", "L14b", "-"),
            ("L11", "L14b", "-"),
            ("L12b", "L14b", "-"),
            ("L13", "L14b", "-"),
            ("L14", "L14b", "-"),
            ("L14b", "L14b", "-"),
            ("L14b", "L16", "-"),
            ("L14b", "L16b", "-"),
            ("L14b", "L15", "-"),
            ("L14b", "L17", "-"),
        ],
        &[("L15", "L16", "-;@,!@")],
        &[("L16", "L16", "-;@,!@")],
        &[
            ("L3", "L16b", "-"),
            ("L5", "L16b", "-"),
            ("L51", "L16b", "-"),
            ("L6", "L16b", "-"),
            ("L8", "L16b", "-"),
            ("L81", "L16b", "-"),
            ("L9", "L16b", "-"),
            ("L10", "L16b", "-"),
            ("L11", "L16b", "-"),
            ("L12b", "L16b", "-"),
            ("L13", "L16b", "-"),
            ("L15", "L16b", "-"),
            ("L16", "L16b", "-"),
            ("L16b", "L16b", "-"),
            ("L17", "L16b", "-"),
        ],
        &[
            ("L17", "L17", "-"),
            ("L17", "L16", "-"),
            ("L17", "L15", "-"),
            ("L17", "L14", "-"),
            ("L17", "L13", "-"),
            ("L17", "L12b", "-"),
            ("L17", "L12", "-"),
            ("L17", "L11", "-"),
            ("L17", "L10", "-"),
            ("L17", "L9", "-"),
            ("L17", "L8", "-"),
            ("L17", "L81", "-"),
            ("L17", "L51", "-"),
            ("L17", "L5", "-"),
        ],
        &[("L18", "L19", "-;@"), ("L182", "L192", "=;@")],
        &[("L20", "L21", "-")],
        &[("L22", "L23", "-")],
    ],
};
