pub(crate) fn base26_triplet_1(a: &[u8]) -> [u8; 3] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1173 base26_triplet_1
    // INCHI✔️✔️: const char* base26_triplet_1( const unsigned char *a )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     UINT32 b0, b1, h;
    // INCHI✔️✔️:     b0 = (UINT32) a[0];             /* 1111 1111  */
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     b1 = (UINT32) ( a[1] & 0x3f );  /* 0011 1111  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 | b1 << 8 );
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     b1 = (UINT32) ( a[1] & 0xfc );  /* 1111 1100  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 << 8 | b1 ) >> 2;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     return t26[h];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: base26_triplet_1
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_1
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: typedef unsigned int UINT32;
    // INCHI✔️✔️:      As the 2^14 (16384) is very close to 26^3 (17576), a triplet of uppercase
    // INCHI✔️✔️:      letters A..Z encodes 14 bits with good efficiency.
    // INCHI✔️✔️:      For speed, we just tabulate triplets below.
    // INCHI✔️✔️:
    // INCHI✔️✔️:      We should throw away 17576-16384= 1192 triplets.
    // INCHI✔️✔️:      These are 676 triplets starting from 'E', the most frequent letter in English
    // INCHI✔️✔️:      texts (the other 516 are those started at 'T' , "TAA" to "TTV").
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_1

    let b0 = u32::from(a[0]);
    let b1 = u32::from(a[1] & 0x3f);
    let table_index = b0 | (b1 << 8);
    let alphabet_index = if table_index < 2_704 {
        table_index
    } else if table_index < 12_168 {
        table_index + 676
    } else {
        table_index + 1_192
    };
    [
        b'A' + (alphabet_index / 676) as u8,
        b'A' + ((alphabet_index / 26) % 26) as u8,
        b'A' + (alphabet_index % 26) as u8,
    ]
}

pub(crate) fn base26_triplet_2(a: &[u8]) -> [u8; 3] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1191 base26_triplet_2
    // INCHI✔️✔️: const char* base26_triplet_2( const unsigned char *a )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     UINT32 b0, b1, b2, h;
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     b0 = (UINT32) ( a[1] & 0xc0 );   /* 1100 0000  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[2] );         /* 1111 1111  */
    // INCHI✔️✔️:     b2 = (UINT32) ( a[3] & 0x0f );  /* 0000 1111  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 | b1 << 8 | b2 << 16 ) >> 6;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     b0 = (UINT32) ( a[1] & 0x03 );   /* 0000 0011 */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[2] );         /* 1111 1111 */
    // INCHI✔️✔️:     b2 = (UINT32) ( a[3] & 0xf0 );  /* 1111 0000 */
    // INCHI✔️✔️:     h = (UINT32) ( b0 << 16 | b1 << 8 | b2 ) >> 4;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     return t26[h];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: base26_triplet_2
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_2
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: typedef unsigned int UINT32;
    // INCHI✔️✔️:      As the 2^14 (16384) is very close to 26^3 (17576), a triplet of uppercase
    // INCHI✔️✔️:      letters A..Z encodes 14 bits with good efficiency.
    // INCHI✔️✔️:      For speed, we just tabulate triplets below.
    // INCHI✔️✔️:
    // INCHI✔️✔️:      We should throw away 17576-16384= 1192 triplets.
    // INCHI✔️✔️:      These are 676 triplets starting from 'E', the most frequent letter in English
    // INCHI✔️✔️:      texts (the other 516 are those started at 'T' , "TAA" to "TTV").
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_2

    let b0 = u32::from(a[1] & 0xc0);
    let b1 = u32::from(a[2]);
    let b2 = u32::from(a[3] & 0x0f);
    let table_index = (b0 | (b1 << 8) | (b2 << 16)) >> 6;
    let alphabet_index = if table_index < 2_704 {
        table_index
    } else if table_index < 12_168 {
        table_index + 676
    } else {
        table_index + 1_192
    };
    [
        b'A' + (alphabet_index / 676) as u8,
        b'A' + ((alphabet_index / 26) % 26) as u8,
        b'A' + (alphabet_index % 26) as u8,
    ]
}

pub(crate) fn base26_triplet_3(a: &[u8]) -> [u8; 3] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1212 base26_triplet_3
    // INCHI✔️✔️: const char* base26_triplet_3( const unsigned char *a )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     UINT32 b0, b1, b2, h;
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     b0 = (UINT32) ( a[3] & 0xf0 );   /* 1111 0000  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[4] );         /* 1111 1111  */
    // INCHI✔️✔️:     b2 = (UINT32) ( a[5] & 0x03 );  /* 0000 0011  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 | b1 << 8 | b2 << 16 ) >> 4;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     b0 = (UINT32) ( a[3] & 0x0f );   /* 0000 1111 */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[4] );         /* 1111 1111 */
    // INCHI✔️✔️:     b2 = (UINT32) ( a[5] & 0xc0 );  /* 1100 0000 */
    // INCHI✔️✔️:     h = (UINT32) ( b0 << 16 | b1 << 8 | b2 ) >> 6;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     return t26[h];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: base26_triplet_3
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_3
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: typedef unsigned int UINT32;
    // INCHI✔️✔️:      As the 2^14 (16384) is very close to 26^3 (17576), a triplet of uppercase
    // INCHI✔️✔️:      letters A..Z encodes 14 bits with good efficiency.
    // INCHI✔️✔️:      For speed, we just tabulate triplets below.
    // INCHI✔️✔️:
    // INCHI✔️✔️:      We should throw away 17576-16384= 1192 triplets.
    // INCHI✔️✔️:      These are 676 triplets starting from 'E', the most frequent letter in English
    // INCHI✔️✔️:      texts (the other 516 are those started at 'T' , "TAA" to "TTV").
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_3

    let b0 = u32::from(a[3] & 0xf0);
    let b1 = u32::from(a[4]);
    let b2 = u32::from(a[5] & 0x03);
    let table_index = (b0 | (b1 << 8) | (b2 << 16)) >> 4;
    let alphabet_index = if table_index < 2_704 {
        table_index
    } else if table_index < 12_168 {
        table_index + 676
    } else {
        table_index + 1_192
    };
    [
        b'A' + (alphabet_index / 676) as u8,
        b'A' + ((alphabet_index / 26) % 26) as u8,
        b'A' + (alphabet_index % 26) as u8,
    ]
}

pub(crate) fn base26_triplet_4(a: &[u8]) -> [u8; 3] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1233 base26_triplet_4
    // INCHI✔️✔️: const char* base26_triplet_4( const unsigned char *a )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     UINT32 b0, b1, h;
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     b0 = (UINT32) ( a[5] & 0xfc );   /* 1111 1100  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[6] );         /* 1111 1111  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 | b1 << 8 ) >> 2;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     b0 = (UINT32) ( a[5] & 0x3f );   /* 0011 1111  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[6] );         /* 1111 1111  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 << 8 | b1 );
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     return t26[h];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: base26_triplet_4
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_4
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: typedef unsigned int UINT32;
    // INCHI✔️✔️:      As the 2^14 (16384) is very close to 26^3 (17576), a triplet of uppercase
    // INCHI✔️✔️:      letters A..Z encodes 14 bits with good efficiency.
    // INCHI✔️✔️:      For speed, we just tabulate triplets below.
    // INCHI✔️✔️:
    // INCHI✔️✔️:      We should throw away 17576-16384= 1192 triplets.
    // INCHI✔️✔️:      These are 676 triplets starting from 'E', the most frequent letter in English
    // INCHI✔️✔️:      texts (the other 516 are those started at 'T' , "TAA" to "TTV").
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_triplet_4

    let b0 = u32::from(a[5] & 0xfc);
    let b1 = u32::from(a[6]);
    let table_index = (b0 | (b1 << 8)) >> 2;
    let alphabet_index = if table_index < 2_704 {
        table_index
    } else if table_index < 12_168 {
        table_index + 676
    } else {
        table_index + 1_192
    };
    [
        b'A' + (alphabet_index / 676) as u8,
        b'A' + ((alphabet_index / 26) % 26) as u8,
        b'A' + (alphabet_index % 26) as u8,
    ]
}

pub(crate) fn base26_dublet_for_bits_28_to_36(a: &[u8]) -> [u8; 2] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1262 base26_dublet_for_bits_28_to_36
    // INCHI✔️✔️: const char* base26_dublet_for_bits_28_to_36( unsigned char *a )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     UINT32 b0, b1, h;
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     b0 = (UINT32) ( a[3] & 0xf0 );    /* 1111 0000  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[4] & 0x1f );    /* 0001 1111  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 | b1 << 8 ) >> 4;
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     b0 = (UINT32) ( a[3] & 0x0f );    /* 0000 1111  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[4] & 0xf8 );    /* 1111 1000  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 << 8 | b1 ) >> 3;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     return d26[h];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: base26_dublet_for_bits_28_to_36
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_dublet_for_bits_28_to_36
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: typedef unsigned int UINT32;
    // INCHI✔️✔️: /*
    // INCHI✔️✔️:     Doublets
    // INCHI✔️✔️: */
    // INCHI✔️✔️: static const char d26[][3] =
    // INCHI✔️✔️: {
    // INCHI✔️✔️: "AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP",
    // INCHI✔️✔️: "AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD","BE","BF",
    // INCHI✔️✔️: "BG","BH","BI","BJ","BK","BL","BM","BN","BO","BP","BQ","BR","BS","BT","BU","BV",
    // INCHI✔️✔️: "BW","BX","BY","BZ","CA","CB","CC","CD","CE","CF","CG","CH","CI","CJ","CK","CL",
    // INCHI✔️✔️: "CM","CN","CO","CP","CQ","CR","CS","CT","CU","CV","CW","CX","CY","CZ","DA","DB",
    // INCHI✔️✔️: "DC","DD","DE","DF","DG","DH","DI","DJ","DK","DL","DM","DN","DO","DP","DQ","DR",
    // INCHI✔️✔️: "DS","DT","DU","DV","DW","DX","DY","DZ","EA","EB","EC","ED","EE","EF","EG","EH",
    // INCHI✔️✔️: "EI","EJ","EK","EL","EM","EN","EO","EP","EQ","ER","ES","ET","EU","EV","EW","EX",
    // INCHI✔️✔️: "EY","EZ","FA","FB","FC","FD","FE","FF","FG","FH","FI","FJ","FK","FL","FM","FN",
    // INCHI✔️✔️: "FO","FP","FQ","FR","FS","FT","FU","FV","FW","FX","FY","FZ","GA","GB","GC","GD",
    // INCHI✔️✔️: "GE","GF","GG","GH","GI","GJ","GK","GL","GM","GN","GO","GP","GQ","GR","GS","GT",
    // INCHI✔️✔️: "GU","GV","GW","GX","GY","GZ","HA","HB","HC","HD","HE","HF","HG","HH","HI","HJ",
    // INCHI✔️✔️: "HK","HL","HM","HN","HO","HP","HQ","HR","HS","HT","HU","HV","HW","HX","HY","HZ",
    // INCHI✔️✔️: "IA","IB","IC","ID","IE","IF","IG","IH","II","IJ","IK","IL","IM","IN","IO","IP",
    // INCHI✔️✔️: "IQ","IR","IS","IT","IU","IV","IW","IX","IY","IZ","JA","JB","JC","JD","JE","JF",
    // INCHI✔️✔️: "JG","JH","JI","JJ","JK","JL","JM","JN","JO","JP","JQ","JR","JS","JT","JU","JV",
    // INCHI✔️✔️: "JW","JX","JY","JZ","KA","KB","KC","KD","KE","KF","KG","KH","KI","KJ","KK","KL",
    // INCHI✔️✔️: "KM","KN","KO","KP","KQ","KR","KS","KT","KU","KV","KW","KX","KY","KZ","LA","LB",
    // INCHI✔️✔️: "LC","LD","LE","LF","LG","LH","LI","LJ","LK","LL","LM","LN","LO","LP","LQ","LR",
    // INCHI✔️✔️: "LS","LT","LU","LV","LW","LX","LY","LZ","MA","MB","MC","MD","ME","MF","MG","MH",
    // INCHI✔️✔️: "MI","MJ","MK","ML","MM","MN","MO","MP","MQ","MR","MS","MT","MU","MV","MW","MX",
    // INCHI✔️✔️: "MY","MZ","NA","NB","NC","ND","NE","NF","NG","NH","NI","NJ","NK","NL","NM","NN",
    // INCHI✔️✔️: "NO","NP","NQ","NR","NS","NT","NU","NV","NW","NX","NY","NZ","OA","OB","OC","OD",
    // INCHI✔️✔️: "OE","OF","OG","OH","OI","OJ","OK","OL","OM","ON","OO","OP","OQ","OR","OS","OT",
    // INCHI✔️✔️: "OU","OV","OW","OX","OY","OZ","PA","PB","PC","PD","PE","PF","PG","PH","PI","PJ",
    // INCHI✔️✔️: "PK","PL","PM","PN","PO","PP","PQ","PR","PS","PT","PU","PV","PW","PX","PY","PZ",
    // INCHI✔️✔️: "QA","QB","QC","QD","QE","QF","QG","QH","QI","QJ","QK","QL","QM","QN","QO","QP",
    // INCHI✔️✔️: "QQ","QR","QS","QT","QU","QV","QW","QX","QY","QZ","RA","RB","RC","RD","RE","RF",
    // INCHI✔️✔️: "RG","RH","RI","RJ","RK","RL","RM","RN","RO","RP","RQ","RR","RS","RT","RU","RV",
    // INCHI✔️✔️: "RW","RX","RY","RZ","SA","SB","SC","SD","SE","SF","SG","SH","SI","SJ","SK","SL",
    // INCHI✔️✔️: "SM","SN","SO","SP","SQ","SR","SS","ST","SU","SV","SW","SX","SY","SZ","TA","TB",
    // INCHI✔️✔️: "TC","TD","TE","TF","TG","TH","TI","TJ","TK","TL","TM","TN","TO","TP","TQ","TR",
    // INCHI✔️✔️: "TS","TT","TU","TV","TW","TX","TY","TZ","UA","UB","UC","UD","UE","UF","UG","UH",
    // INCHI✔️✔️: "UI","UJ","UK","UL","UM","UN","UO","UP","UQ","UR","US","UT","UU","UV","UW","UX",
    // INCHI✔️✔️: "UY","UZ","VA","VB","VC","VD","VE","VF","VG","VH","VI","VJ","VK","VL","VM","VN",
    // INCHI✔️✔️: "VO","VP","VQ","VR","VS","VT","VU","VV","VW","VX","VY","VZ","WA","WB","WC","WD",
    // INCHI✔️✔️: "WE","WF","WG","WH","WI","WJ","WK","WL","WM","WN","WO","WP","WQ","WR","WS","WT",
    // INCHI✔️✔️: "WU","WV","WW","WX","WY","WZ","XA","XB","XC","XD","XE","XF","XG","XH","XI","XJ",
    // INCHI✔️✔️: "XK","XL","XM","XN","XO","XP","XQ","XR","XS","XT","XU","XV","XW","XX","XY","XZ",
    // INCHI✔️✔️: "YA","YB","YC","YD","YE","YF","YG","YH","YI","YJ","YK","YL","YM","YN","YO","YP",
    // INCHI✔️✔️: "YQ","YR","YS","YT","YU","YV","YW","YX","YY","YZ","ZA","ZB","ZC","ZD","ZE","ZF",
    // INCHI✔️✔️: "ZG","ZH","ZI","ZJ","ZK","ZL","ZM","ZN","ZO","ZP","ZQ","ZR","ZS","ZT","ZU","ZV",
    // INCHI✔️✔️: "ZW","ZX","ZY","ZZ"
    // INCHI✔️✔️: };
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_dublet_for_bits_28_to_36

    let b0 = u32::from(a[3] & 0xf0);
    let b1 = u32::from(a[4] & 0x1f);
    let table_index = (b0 | (b1 << 8)) >> 4;
    [
        b'A' + (table_index / 26) as u8,
        b'A' + (table_index % 26) as u8,
    ]
}

pub(crate) fn base26_dublet_for_bits_56_to_64(a: &[u8]) -> [u8; 2] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1287 base26_dublet_for_bits_56_to_64
    // INCHI✔️✔️: const char* base26_dublet_for_bits_56_to_64( unsigned char *a )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     UINT32 b0, b1, h;
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     b0 = (UINT32) ( a[7] );            /* 1111 1111  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[8] & 0x01 );    /* 0000 0001  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 | b1 << 8 );
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     b0 = (UINT32) ( a[7] );            /* 1111 1111  */
    // INCHI✔️✔️:     b1 = (UINT32) ( a[8] & 0x80 );    /* 1000 0000  */
    // INCHI✔️✔️:     h = (UINT32) ( b0 << 8 | b1 ) >> 7;
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     return d26[h];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: base26_dublet_for_bits_56_to_64
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_dublet_for_bits_56_to_64
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: typedef unsigned int UINT32;
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // The complete verbatim `d26` definition shared by both active doublet
    // functions is reproduced above from ikey_base26.c:1111-1160.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: base26_dublet_for_bits_56_to_64

    let b0 = u32::from(a[7]);
    let b1 = u32::from(a[8] & 0x01);
    let table_index = b0 | (b1 << 8);
    [
        b'A' + (table_index / 26) as u8,
        b'A' + (table_index % 26) as u8,
    ]
}

pub(crate) fn get_xtra_hash_major_hex(a: &[u8], sz_xtra: &mut [u8]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1307 get_xtra_hash_major_hex
    // INCHI✔️✔️: void get_xtra_hash_major_hex( const unsigned char *a, char* szXtra )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     unsigned char c;
    // INCHI✔️✔️:     int i, j, start_byte = 8;
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     c = a[start_byte] & 0xfe;  /* 1111 1110  */
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     c = a[start_byte] & 0x7f;  /* 0111 1111  */
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     j = sprintf(szXtra, "%02x", c);
    // INCHI✔️✔️:     for (i = start_byte + 1; i < 32; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         j += sprintf(szXtra + j, "%02x", a[i]);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: get_xtra_hash_major_hex
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_xtra_hash_major_hex
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: #include <stdio.h>
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_xtra_hash_major_hex

    const HEX_LOWER: &[u8; 16] = b"0123456789abcdef";
    let mut output_index = 0_usize;
    for input_index in 8_usize..32 {
        let byte = if input_index == 8 {
            a[input_index] & 0xfe
        } else {
            a[input_index]
        };
        sz_xtra[output_index] = HEX_LOWER[usize::from(byte >> 4)];
        sz_xtra[output_index + 1] = HEX_LOWER[usize::from(byte & 0x0f)];
        output_index += 2;
    }
    sz_xtra[output_index] = 0;
}

pub(crate) fn get_xtra_hash_minor_hex(a: &[u8], sz_xtra: &mut [u8]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1328 get_xtra_hash_minor_hex
    // INCHI✔️✔️: void get_xtra_hash_minor_hex( const unsigned char *a, char* szXtra )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     unsigned char c;
    // INCHI✔️✔️:     int i, j, start_byte = 4;
    // INCHI✔️✔️: #ifndef FIX_BASE26_ENC_BUG
    // INCHI✔️✔️:     c = a[start_byte] & 0xe0;  /* 1110 0000  */
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     c = a[start_byte] & 0x07;  /* 0000 0111  */
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:     j = sprintf(szXtra, "%02x", c);
    // INCHI✔️✔️:     for (i = start_byte + 1; i < 32; i++)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         j += sprintf(szXtra + j, "%02x", a[i]);
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: get_xtra_hash_minor_hex
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_xtra_hash_minor_hex
    // INCHI✔️✔️: /*#define FIX_BASE26_ENC_BUG 1*/
    // INCHI✔️✔️: #include <stdio.h>
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; FIX_BASE26_ENC_BUG is undefined.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: get_xtra_hash_minor_hex

    const HEX_LOWER: &[u8; 16] = b"0123456789abcdef";
    let mut output_index = 0_usize;
    for input_index in 4_usize..32 {
        let byte = if input_index == 4 {
            a[input_index] & 0xe0
        } else {
            a[input_index]
        };
        sz_xtra[output_index] = HEX_LOWER[usize::from(byte >> 4)];
        sz_xtra[output_index + 1] = HEX_LOWER[usize::from(byte & 0x0f)];
        output_index += 2;
    }
    sz_xtra[output_index] = 0;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ikey_base26__base26_triplet_1__line_1173() {
        let input_for_index = |index: u32| [index as u8, ((index >> 8) & 0x3f) as u8];
        let cases = [
            (0, *b"AAA"),
            (255, *b"AJV"),
            (256, *b"AJW"),
            (2_703, *b"DZZ"),
            (2_704, *b"FAA"),
            (2_705, *b"FAB"),
            (12_167, *b"SZZ"),
            (12_168, *b"TTW"),
            (12_169, *b"TTX"),
            (16_383, *b"ZZZ"),
        ];
        for (index, expected) in cases {
            assert_eq!(base26_triplet_1(&input_for_index(index)), expected);
        }

        for low_byte in u8::MIN..=u8::MAX {
            for low_six_bits in 0_u8..=0x3f {
                let expected = base26_triplet_1(&[low_byte, low_six_bits]);
                assert_eq!(base26_triplet_1(&[low_byte, low_six_bits | 0x40]), expected);
                assert_eq!(base26_triplet_1(&[low_byte, low_six_bits | 0x80]), expected);
                assert_eq!(base26_triplet_1(&[low_byte, low_six_bits | 0xc0]), expected);
            }
        }
    }

    #[test]
    fn official_c_oracle__base26_triplet_1__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--base26-triplet-1-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "base26_triplet_1");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 2);
            let expected = parse_bytes(&official["output"]["bytes"], 3);
            let rust = base26_triplet_1(&input);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["nul"], 0, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 16_388);
    }

    #[test]
    fn source_port__ikey_base26__base26_triplet_2__line_1191() {
        let input_for_index = |index: u32| {
            [
                0,
                ((index & 0x03) << 6) as u8,
                ((index >> 2) & 0xff) as u8,
                ((index >> 10) & 0x0f) as u8,
            ]
        };
        let cases = [
            (0, *b"AAA"),
            (255, *b"AJV"),
            (256, *b"AJW"),
            (2_703, *b"DZZ"),
            (2_704, *b"FAA"),
            (2_705, *b"FAB"),
            (12_167, *b"SZZ"),
            (12_168, *b"TTW"),
            (12_169, *b"TTX"),
            (16_383, *b"ZZZ"),
        ];
        for (index, expected) in cases {
            assert_eq!(base26_triplet_2(&input_for_index(index)), expected);
        }

        for index in 0_u32..16_384 {
            let input = input_for_index(index);
            let expected = base26_triplet_2(&input);
            assert_eq!(
                base26_triplet_2(&[0xff, input[1] | 0x3f, input[2], input[3] | 0xf0]),
                expected
            );
        }
    }

    #[test]
    fn official_c_oracle__base26_triplet_2__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--base26-triplet-2-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "base26_triplet_2");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 4);
            let expected = parse_bytes(&official["output"]["bytes"], 3);
            let rust = base26_triplet_2(&input);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["nul"], 0, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 16_388);
    }

    #[test]
    fn source_port__ikey_base26__base26_triplet_3__line_1212() {
        let input_for_index = |index: u32| {
            [
                0,
                0,
                0,
                ((index & 0x0f) << 4) as u8,
                ((index >> 4) & 0xff) as u8,
                ((index >> 12) & 0x03) as u8,
            ]
        };
        let cases = [
            (0, *b"AAA"),
            (255, *b"AJV"),
            (256, *b"AJW"),
            (2_703, *b"DZZ"),
            (2_704, *b"FAA"),
            (2_705, *b"FAB"),
            (12_167, *b"SZZ"),
            (12_168, *b"TTW"),
            (12_169, *b"TTX"),
            (16_383, *b"ZZZ"),
        ];
        for (index, expected) in cases {
            assert_eq!(base26_triplet_3(&input_for_index(index)), expected);
        }

        for index in 0_u32..16_384 {
            let input = input_for_index(index);
            let expected = base26_triplet_3(&input);
            assert_eq!(
                base26_triplet_3(&[0xff, 0xa5, 0x5a, input[3] | 0x0f, input[4], input[5] | 0xfc,]),
                expected
            );
        }
    }

    #[test]
    fn official_c_oracle__base26_triplet_3__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--base26-triplet-3-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "base26_triplet_3");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 6);
            let expected = parse_bytes(&official["output"]["bytes"], 3);
            let rust = base26_triplet_3(&input);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["nul"], 0, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 16_388);
    }

    #[test]
    fn source_port__ikey_base26__base26_triplet_4__line_1233() {
        let input_for_index = |index: u32| {
            [
                0,
                0,
                0,
                0,
                0,
                ((index & 0x3f) << 2) as u8,
                ((index >> 6) & 0xff) as u8,
            ]
        };
        let cases = [
            (0, *b"AAA"),
            (255, *b"AJV"),
            (256, *b"AJW"),
            (2_703, *b"DZZ"),
            (2_704, *b"FAA"),
            (2_705, *b"FAB"),
            (12_167, *b"SZZ"),
            (12_168, *b"TTW"),
            (12_169, *b"TTX"),
            (16_383, *b"ZZZ"),
        ];
        for (index, expected) in cases {
            assert_eq!(base26_triplet_4(&input_for_index(index)), expected);
        }

        for index in 0_u32..16_384 {
            let input = input_for_index(index);
            let expected = base26_triplet_4(&input);
            assert_eq!(
                base26_triplet_4(&[0xff, 0xa5, 0x5a, 0x3c, 0xc3, input[5] | 0x03, input[6],]),
                expected
            );
        }
    }

    #[test]
    fn official_c_oracle__base26_triplet_4__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--base26-triplet-4-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "base26_triplet_4");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 7);
            let expected = parse_bytes(&official["output"]["bytes"], 3);
            let rust = base26_triplet_4(&input);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["nul"], 0, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 16_388);
    }

    #[test]
    fn source_port__ikey_base26__base26_dublet_for_bits_28_to_36__line_1262() {
        let input_for_index = |index: u32| {
            [
                0,
                0,
                0,
                ((index & 0x0f) << 4) as u8,
                ((index >> 4) & 0x1f) as u8,
            ]
        };
        let cases = [
            (0, *b"AA"),
            (1, *b"AB"),
            (25, *b"AZ"),
            (26, *b"BA"),
            (255, *b"JV"),
            (256, *b"JW"),
            (510, *b"TQ"),
            (511, *b"TR"),
        ];
        for (index, expected) in cases {
            assert_eq!(
                base26_dublet_for_bits_28_to_36(&input_for_index(index)),
                expected
            );
        }

        for index in 0_u32..512 {
            let input = input_for_index(index);
            let expected = [b'A' + (index / 26) as u8, b'A' + (index % 26) as u8];
            assert_eq!(base26_dublet_for_bits_28_to_36(&input), expected);
            assert_eq!(
                base26_dublet_for_bits_28_to_36(&[
                    0xff,
                    0xa5,
                    0x5a,
                    input[3] | 0x0f,
                    input[4] | 0xe0,
                ]),
                expected
            );
        }
    }

    #[test]
    fn official_c_oracle__base26_dublet_for_bits_28_to_36__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--base26-dublet-for-bits-28-to-36-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "base26_dublet_for_bits_28_to_36");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 5);
            let expected = parse_bytes(&official["output"]["bytes"], 2);
            let rust = base26_dublet_for_bits_28_to_36(&input);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["nul"], 0, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 516);
    }

    #[test]
    fn source_port__ikey_base26__base26_dublet_for_bits_56_to_64__line_1287() {
        let input_for_index = |index: u32| {
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                (index & 0xff) as u8,
                ((index >> 8) & 0x01) as u8,
            ]
        };
        let cases = [
            (0, *b"AA"),
            (1, *b"AB"),
            (25, *b"AZ"),
            (26, *b"BA"),
            (255, *b"JV"),
            (256, *b"JW"),
            (510, *b"TQ"),
            (511, *b"TR"),
        ];
        for (index, expected) in cases {
            assert_eq!(
                base26_dublet_for_bits_56_to_64(&input_for_index(index)),
                expected
            );
        }

        for index in 0_u32..512 {
            let input = input_for_index(index);
            let expected = [b'A' + (index / 26) as u8, b'A' + (index % 26) as u8];
            assert_eq!(base26_dublet_for_bits_56_to_64(&input), expected);
            assert_eq!(
                base26_dublet_for_bits_56_to_64(&[
                    0xff,
                    0xa5,
                    0x5a,
                    0x3c,
                    0xc3,
                    0x69,
                    0x96,
                    input[7],
                    input[8] | 0xfe,
                ]),
                expected
            );
        }
    }

    #[test]
    fn official_c_oracle__base26_dublet_for_bits_56_to_64__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--base26-dublet-for-bits-56-to-64-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "base26_dublet_for_bits_56_to_64");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 9);
            let expected = parse_bytes(&official["output"]["bytes"], 2);
            let rust = base26_dublet_for_bits_56_to_64(&input);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["nul"], 0, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 516);
    }

    #[test]
    fn source_port__ikey_base26__get_xtra_hash_major_hex__line_1307() {
        const HEX_LOWER: &[u8; 16] = b"0123456789abcdef";

        for start_byte in u8::MIN..=u8::MAX {
            let mut input = [0_u8; 32];
            for (index, byte) in input.iter_mut().enumerate() {
                *byte = start_byte.wrapping_add((index as u8).wrapping_mul(17));
            }
            input[8] = start_byte;
            let input_before = input;
            let mut output = [0xa5_u8; 64];

            get_xtra_hash_major_hex(&input, &mut output);

            for input_index in 8_usize..32 {
                let byte = if input_index == 8 {
                    input[input_index] & 0xfe
                } else {
                    input[input_index]
                };
                let output_index = (input_index - 8) * 2;
                assert_eq!(
                    output[output_index],
                    HEX_LOWER[usize::from(byte >> 4)],
                    "start_byte={start_byte}, input_index={input_index}"
                );
                assert_eq!(
                    output[output_index + 1],
                    HEX_LOWER[usize::from(byte & 0x0f)],
                    "start_byte={start_byte}, input_index={input_index}"
                );
            }
            assert_eq!(output[48], 0, "start_byte={start_byte}");
            assert_eq!(&output[49..], &[0xa5; 15], "start_byte={start_byte}");
            assert_eq!(input, input_before, "start_byte={start_byte}");
        }
    }

    #[test]
    fn official_c_oracle__get_xtra_hash_major_hex__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--get-xtra-hash-major-hex-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "get_xtra_hash_major_hex");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 32);
            let expected = parse_bytes(&official["output"]["bytes"], 64);
            let mut rust = [0xa5_u8; 64];
            get_xtra_hash_major_hex(&input, &mut rust);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 256);
    }

    #[test]
    fn source_port__ikey_base26__get_xtra_hash_minor_hex__line_1328() {
        const HEX_LOWER: &[u8; 16] = b"0123456789abcdef";

        for start_byte in u8::MIN..=u8::MAX {
            let mut input = [0_u8; 32];
            for (index, byte) in input.iter_mut().enumerate() {
                *byte = start_byte.wrapping_add((index as u8).wrapping_mul(17));
            }
            input[4] = start_byte;
            let input_before = input;
            let mut output = [0xa5_u8; 64];

            get_xtra_hash_minor_hex(&input, &mut output);

            for input_index in 4_usize..32 {
                let byte = if input_index == 4 {
                    input[input_index] & 0xe0
                } else {
                    input[input_index]
                };
                let output_index = (input_index - 4) * 2;
                assert_eq!(
                    output[output_index],
                    HEX_LOWER[usize::from(byte >> 4)],
                    "start_byte={start_byte}, input_index={input_index}"
                );
                assert_eq!(
                    output[output_index + 1],
                    HEX_LOWER[usize::from(byte & 0x0f)],
                    "start_byte={start_byte}, input_index={input_index}"
                );
            }
            assert_eq!(output[56], 0, "start_byte={start_byte}");
            assert_eq!(&output[57..], &[0xa5; 7], "start_byte={start_byte}");
            assert_eq!(input, input_before, "start_byte={start_byte}");
        }
    }

    #[test]
    fn official_c_oracle__get_xtra_hash_minor_hex__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tests/tools/inchi_official_c_oracle/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--get-xtra-hash-minor-hex-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "get_xtra_hash_minor_hex");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_bytes = |path: &Value, expected_length: usize| {
                let values = path
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_bytes(&official["input"]["bytes"], 32);
            let expected = parse_bytes(&official["output"]["bytes"], 64);
            let mut rust = [0xa5_u8; 64];
            get_xtra_hash_minor_hex(&input, &mut rust);
            assert_eq!(rust.as_slice(), expected.as_slice(), "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 256);
    }
}
