#[allow(non_snake_case)]
pub(crate) fn EquString(equivalence_value: i32) -> &'static str {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:564 EquString
    // INCHI✔️✔️: const char *EquString( int EquVal )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int bFrom = EquVal & ( iiSTEREO | iiSTEREO_INV | iiNUMB | iiEQU );
    // INCHI✔️✔️:     int bType = EquVal & ( iitISO | iitNONTAUT );
    // INCHI✔️✔️:     int bEq2 = EquVal & ( iiEq2NONTAUT | iiEq2ISO | iiEq2INV );
    // INCHI✔️✔️:     const char *r = "";
    // INCHI✔️✔️:
    // INCHI✔️✔️: #if ( FIX_EMPTY_LAYER_BUG == 1 )
    // INCHI✔️✔️:     int bEmpty = EquVal & iiEmpty;
    // INCHI✔️✔️:     if (bEmpty)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         r = "e";
    // INCHI✔️✔️:         return r;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:
    // INCHI✔️✔️:     switch (bFrom)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:
    // INCHI✔️✔️:         case iiSTEREO:  /* ------------ Stereo --------------------*/
    // INCHI✔️✔️:             switch (bType)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 case iitISO:  /* iso main stereo =... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";            /* iso main stereo = main stereo */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";           /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case iitNONTAUT: /* non-taut stereo =... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";            /* non-taut stereo = main stereo */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";           /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case ( iitNONTAUT | iitISO ): /* iso non-taut stereo = ... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";            /* iso non-taut stereo = main stereo */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2ISO:
    // INCHI✔️✔️:                             r = "M";            /* iso non-taut stereo = main iso stereo */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2NONTAUT:
    // INCHI✔️✔️:                             r = "n";            /* iso non-taut stereo = non-taut stereo */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";           /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 default:
    // INCHI✔️✔️:                     r = "??";           /* should not happen */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         case iiSTEREO_INV: /*---------- Inverted Aux Stereo ------*/
    // INCHI✔️✔️:             if (bEq2 & iiEq2INV)
    // INCHI✔️✔️:             { /* stereo = Inverted(another stereo) */
    // INCHI✔️✔️:                 bEq2 &= ~iiEq2INV;
    // INCHI✔️✔️:                 switch (bType)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     case 0: /* main = ...*/
    // INCHI✔️✔️:                         switch (bEq2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 r = "im";       /* main         = Inv(main) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2ISO:
    // INCHI✔️✔️:                                 r = "iM";       /* main         = Inv(main iso) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2NONTAUT:
    // INCHI✔️✔️:                                 r = "in";       /* maim         = Inv(non-taut) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case ( iiEq2NONTAUT | iiEq2ISO ):
    // INCHI✔️✔️:                                 r = "iN";       /* maim         = Inv(non-taut iso ) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             default:
    // INCHI✔️✔️:                                 r = "??";           /* should not happen */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     case iitISO: /* main iso = ...*/
    // INCHI✔️✔️:                         switch (bEq2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 r = "im";       /* main iso     = Inv(main) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2ISO:
    // INCHI✔️✔️:                                 r = "iM";       /* main iso     = Inv(main iso) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2NONTAUT:
    // INCHI✔️✔️:                                 r = "in";       /* maim iso     = Inv(non-taut) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case ( iiEq2NONTAUT | iiEq2ISO ):
    // INCHI✔️✔️:                                 r = "iN";       /* maim         = Inv(non-taut iso ) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             default:
    // INCHI✔️✔️:                                 r = "??";           /* should not happen */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     case iitNONTAUT: /* non-taut = ... */
    // INCHI✔️✔️:                         switch (bEq2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 r = "im";       /* non-taut     = Inv(main) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2ISO:
    // INCHI✔️✔️:                                 r = "iM";       /* non-taut     = Inv(main iso) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2NONTAUT:
    // INCHI✔️✔️:                                 r = "in";       /* non-taut     = Inv(non-taut) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case ( iiEq2NONTAUT | iiEq2ISO ):
    // INCHI✔️✔️:                                 r = "iN";       /* non-taut     = Inv(non-taut iso ) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             default:
    // INCHI✔️✔️:                                 r = "??";           /* should not happen */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     case ( iitNONTAUT | iitISO ):
    // INCHI✔️✔️:                         switch (bEq2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 r = "im";       /* non-taut iso = Inv(main) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2ISO:
    // INCHI✔️✔️:                                 r = "iM";       /* non-taut iso = Inv(main iso) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2NONTAUT:
    // INCHI✔️✔️:                                 r = "in";       /* non-taut iso = Inv(non-taut) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case ( iiEq2NONTAUT | iiEq2ISO ):
    // INCHI✔️✔️:                                 r = "iN";       /* non-taut iso = Inv(non-taut iso ) */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             default:
    // INCHI✔️✔️:                                 r = "??";           /* should not happen */
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     default:
    // INCHI✔️✔️:                         r = "??";           /* should not happen */
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             else
    // INCHI✔️✔️:             {  /* Inv stereo = another (non-inverted) stereo */
    // INCHI✔️✔️:
    // INCHI✔️✔️:                 switch (bType)
    // INCHI✔️✔️:                 {
    // INCHI✔️✔️:                     case iitISO: /* main iso = ...*/
    // INCHI✔️✔️:                         switch (bEq2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 r = "m";       /* main         = (inverted aux) main */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             default:
    // INCHI✔️✔️:                                 r = "??";           /* should not happen */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     case iitNONTAUT: /* non-taut = ... */
    // INCHI✔️✔️:                         switch (bEq2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 r = "m";       /* non-taut     = (inverted aux) main */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             default:
    // INCHI✔️✔️:                                 r = "??";           /* should not happen */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     case ( iitNONTAUT | iitISO ): /* non-taut iso = ...*/
    // INCHI✔️✔️:                         switch (bEq2)
    // INCHI✔️✔️:                         {
    // INCHI✔️✔️:                             case 0:
    // INCHI✔️✔️:                                 r = "m";        /* non-taut iso  = (inverted aux) main */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2ISO:
    // INCHI✔️✔️:                                 r = "M";       /* non-taut iso  = (inverted aux) main iso */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             case iiEq2NONTAUT:
    // INCHI✔️✔️:                                 r = "n";       /* non-taut iso  = (inverted aux) non-taut */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                             default:
    // INCHI✔️✔️:                                 r = "??";           /* should not happen */
    // INCHI✔️✔️:                                 break;
    // INCHI✔️✔️:                         }
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                     default:
    // INCHI✔️✔️:                         r = "??";           /* should not happen */
    // INCHI✔️✔️:                         break;
    // INCHI✔️✔️:                 }
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         case ( iiNUMB | iiSTEREO_INV ): /*------------- Inv Stereo Numbering ------------*/
    // INCHI✔️✔️:             switch (bType)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 case 0: /* inv stereo numb main = ...*/
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* inv stereo numb main     = main numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case iitISO: /* inv stereo iso numb main = ...*/
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* inv stereo iso numb main = main numb  */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2INV:
    // INCHI✔️✔️:                             r = "im";      /* inv stereo iso numb main = InvStereo(main) numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2ISO:
    // INCHI✔️✔️:                             r = "M";      /* inv stereo iso numb main = isotopic main numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case iitNONTAUT: /* inv stereo numb non-taut = ... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* inv stereo numb non-taut = main numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2NONTAUT:
    // INCHI✔️✔️:                             r = "n";       /* inv stereo numb non-taut = non-taut numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2INV:
    // INCHI✔️✔️:                             r = "im";      /* inv stereo numb non-taut =  InvStereo(main) numb  */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case ( iitNONTAUT | iitISO ): /* inv stereo numb non-taut iso = ... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* inv stereo numb non-taut iso = main numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2ISO:
    // INCHI✔️✔️:                             r = "M";       /* inv stereo numb non-taut iso = main numb iso */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case ( iiEq2ISO | iiEq2INV ):
    // INCHI✔️✔️:                             r = "iM";       /* inv stereo numb non-taut iso = InvStereo(main iso) numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2NONTAUT:
    // INCHI✔️✔️:                             r = "n";       /* inv stereo numb non-taut iso = non-taut numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case ( iiEq2NONTAUT | iiEq2ISO ):
    // INCHI✔️✔️:                             r = "N";       /* inv stereo numb non-taut iso = non-taut iso numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2INV:
    // INCHI✔️✔️:                             r = "im";      /* inv stereo numb non-taut iso = InvStereo(main) numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case ( iiEq2NONTAUT | iiEq2INV ):
    // INCHI✔️✔️:                             r = "in";      /* inv stereo numb non-taut iso = InvStereo(non-taut) numb ) */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";           /* should not happen  */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 default:
    // INCHI✔️✔️:                     r = "??";           /* should not happen */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         case iiNUMB:           /*------------- Canonical Numbering ------------*/
    // INCHI✔️✔️:             switch (bType)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 case 0:         /* numb main = ...*/
    // INCHI✔️✔️:                     r = "??";      /* should not happen */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case iitISO:     /* iso numb main = ...*/
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* iso numb main = main numb  */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case iitNONTAUT: /* numb non-taut = ... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* numb non-taut = main numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case ( iitNONTAUT | iitISO ): /* numb non-taut iso = ... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* numb non-taut iso = main numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2ISO:
    // INCHI✔️✔️:                             r = "M";       /* numb non-taut iso = main numb iso */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2NONTAUT:
    // INCHI✔️✔️:                             r = "n";       /* numb non-taut iso = non-taut numb */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";           /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 default:
    // INCHI✔️✔️:                     r = "??";           /* should not happen */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         case iiEQU:         /*------------- Atom Equivalence ------------*/
    // INCHI✔️✔️:             switch (bType)
    // INCHI✔️✔️:             {
    // INCHI✔️✔️:                 case 0:         /* equivalence main = ...*/
    // INCHI✔️✔️:                     r = "??";      /* should not happen */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case iitISO:     /* equivalence main iso = ...*/
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* equivalence main = main equ  */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case iitNONTAUT: /* equivalence non-taut = ... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* equivalence non-taut = main equ */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 case ( iitNONTAUT | iitISO ): /*  equivalence non-taut iso = ... */
    // INCHI✔️✔️:                     switch (bEq2)
    // INCHI✔️✔️:                     {
    // INCHI✔️✔️:                         case 0:
    // INCHI✔️✔️:                             r = "m";       /* equivalence non-taut iso = main equ */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2ISO:
    // INCHI✔️✔️:                             r = "M";       /* equivalence non-taut iso = main iso equ */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         case iiEq2NONTAUT:
    // INCHI✔️✔️:                             r = "n";       /* equivalence non-taut iso = non-taut equ */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                         default:
    // INCHI✔️✔️:                             r = "??";      /* should not happen */
    // INCHI✔️✔️:                             break;
    // INCHI✔️✔️:                     }
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:                 default:
    // INCHI✔️✔️:                     r = "??";          /* should not happen */
    // INCHI✔️✔️:                     break;
    // INCHI✔️✔️:             }
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:         default:
    // INCHI✔️✔️:             r = "??";      /* should not happen */
    // INCHI✔️✔️:             break;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     return r;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: EquString

    const STEREO: i32 = crate::source_types::iiSTEREO as i32;
    const STEREO_INVERTED: i32 = crate::source_types::iiSTEREO_INV as i32;
    const NUMBERING: i32 = crate::source_types::iiNUMB as i32;
    const EQUIVALENCE: i32 = crate::source_types::iiEQU as i32;
    const ISOTOPIC: i32 = crate::source_types::iitISO as i32;
    const NONTAUTOMERIC: i32 = crate::source_types::iitNONTAUT as i32;
    const EQUAL_NONTAUTOMERIC: i32 = crate::source_types::iiEq2NONTAUT as i32;
    const EQUAL_ISOTOPIC: i32 = crate::source_types::iiEq2ISO as i32;
    const EQUAL_INVERTED: i32 = crate::source_types::iiEq2INV as i32;

    let from = equivalence_value & (STEREO | STEREO_INVERTED | NUMBERING | EQUIVALENCE);
    let layer_type = equivalence_value & (ISOTOPIC | NONTAUTOMERIC);
    let equal_to = equivalence_value & (EQUAL_NONTAUTOMERIC | EQUAL_ISOTOPIC | EQUAL_INVERTED);

    match from {
        STEREO => match layer_type {
            ISOTOPIC | NONTAUTOMERIC => {
                if equal_to == 0 {
                    "m"
                } else {
                    "??"
                }
            }
            value if value == ISOTOPIC | NONTAUTOMERIC => match equal_to {
                0 => "m",
                EQUAL_ISOTOPIC => "M",
                EQUAL_NONTAUTOMERIC => "n",
                _ => "??",
            },
            _ => "??",
        },
        STEREO_INVERTED => {
            if equal_to & EQUAL_INVERTED != 0 {
                match equal_to & !EQUAL_INVERTED {
                    0 => "im",
                    EQUAL_ISOTOPIC => "iM",
                    EQUAL_NONTAUTOMERIC => "in",
                    value if value == EQUAL_NONTAUTOMERIC | EQUAL_ISOTOPIC => "iN",
                    _ => "??",
                }
            } else {
                match layer_type {
                    ISOTOPIC | NONTAUTOMERIC => {
                        if equal_to == 0 {
                            "m"
                        } else {
                            "??"
                        }
                    }
                    value if value == ISOTOPIC | NONTAUTOMERIC => match equal_to {
                        0 => "m",
                        EQUAL_ISOTOPIC => "M",
                        EQUAL_NONTAUTOMERIC => "n",
                        _ => "??",
                    },
                    _ => "??",
                }
            }
        }
        value if value == NUMBERING | STEREO_INVERTED => match layer_type {
            0 => {
                if equal_to == 0 {
                    "m"
                } else {
                    "??"
                }
            }
            ISOTOPIC => match equal_to {
                0 => "m",
                EQUAL_INVERTED => "im",
                EQUAL_ISOTOPIC => "M",
                _ => "??",
            },
            NONTAUTOMERIC => match equal_to {
                0 => "m",
                EQUAL_NONTAUTOMERIC => "n",
                EQUAL_INVERTED => "im",
                _ => "??",
            },
            value if value == ISOTOPIC | NONTAUTOMERIC => match equal_to {
                0 => "m",
                EQUAL_ISOTOPIC => "M",
                value if value == EQUAL_ISOTOPIC | EQUAL_INVERTED => "iM",
                EQUAL_NONTAUTOMERIC => "n",
                value if value == EQUAL_NONTAUTOMERIC | EQUAL_ISOTOPIC => "N",
                EQUAL_INVERTED => "im",
                value if value == EQUAL_NONTAUTOMERIC | EQUAL_INVERTED => "in",
                _ => "??",
            },
            _ => "??",
        },
        NUMBERING | EQUIVALENCE => match layer_type {
            ISOTOPIC | NONTAUTOMERIC => {
                if equal_to == 0 {
                    "m"
                } else {
                    "??"
                }
            }
            value if value == ISOTOPIC | NONTAUTOMERIC => match equal_to {
                0 => "m",
                EQUAL_ISOTOPIC => "M",
                EQUAL_NONTAUTOMERIC => "n",
                _ => "??",
            },
            _ => "??",
        },
        _ => "??",
    }
}

#[allow(non_snake_case)]
pub(crate) fn szGetTag(
    heap: &mut SourceHeap,
    tags: &[INCHI_TAG],
    n_tag: i32,
    tag_bits: i32,
    output: SourceMutPointer<i8>,
    always: &mut i32,
    tag_flag: i16,
) -> Result<SourceMutPointer<i8>, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2087 szGetTag
    // INCHI✔️❌: char *szGetTag( const INCHI_TAG *Tag,
    // INCHI✔️❌:                 int             nTag,
    // INCHI✔️❌:                 int             bTag,
    // INCHI✔️❌:                 char            *szTag,
    // INCHI✔️❌:                 int             *bAlways,
    // INCHI✔️❌:                 short           tag_flag)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, bit, num, len;
    // INCHI✔️❌:     const int MAX_TAG_NUM = tag_flag ? (int)IL_MAX_ORD : (int)AL_MAX_ORD; /* djb-rwth: fixing GHI #160 */
    // INCHI✔️❌:     if (0 < nTag && nTag < 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* no plain text comments: pick up the last tag */
    // INCHI✔️❌:         for (i = 0, j = -1, bit = 1; i < MAX_TAG_NUM; i++, bit <<= 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (bTag & bit)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 j = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (j >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌: #if USE_BCF
    // INCHI✔️❌:             int stl1, stl2, dstsz;
    // INCHI✔️❌:             stl1 = strlen(Tag[j].szXmlLabel) + 1;
    // INCHI✔️❌:             stl2 = strlen(Tag[j].szPlainLabel) + 1;
    // INCHI✔️❌:             dstsz = max_3(stl1, stl2, 5);
    // INCHI✔️❌:             strcpy_s( szTag, dstsz, nTag == 1 ? Tag[j].szXmlLabel : nTag == 2 ? Tag[j].szPlainLabel : "???" ); /* djb-rwth: function replaced with its safe C11 variant */
    // INCHI✔️❌: #else
    // INCHI✔️❌:             strcpy(szTag, nTag == 1 ? Tag[j].szXmlLabel : nTag == 2 ? Tag[j].szPlainLabel : "???"); /* djb-rwth: addressing coverity ID #499488 -- when nTag == 2, the "???" is avoided, which is correct */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             if (nTag != 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *bAlways = Tag[j].bAlwaysOutput;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return szTag;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:         if (nTag == 3)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* plain text with comments */
    // INCHI✔️❌:             szTag[0] = '{';
    // INCHI✔️❌:             szTag[1] = '\0';
    // INCHI✔️❌:             for (i = 0, j = -1, bit = 1, num = 0; i < MAX_TAG_NUM; i++, bit <<= 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bTag & bit)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j = i;
    // INCHI✔️❌:                     if (num++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         strcat(szTag, ":");
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     strcat(szTag, Tag[i].szPlainComment);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(szTag, "}");
    // INCHI✔️❌:                 num = (int) strlen( Tag[j].szPlainLabel );
    // INCHI✔️❌:                 len = (int) strlen( szTag );
    // INCHI✔️❌:                 if (len)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     memmove(szTag + num, szTag, (long long)len + 1); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                     memcpy(szTag, Tag[j].szPlainLabel, num);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     strcpy(szTag, Tag[j].szPlainLabel);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 *bAlways = Tag[j].bAlwaysOutput;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcpy(szTag, "???");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return szTag;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     strcpy(szTag, "???");
    // INCHI✔️❌:     return szTag;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: szGetTag

    fn c_bytes(
        heap: &SourceHeap,
        pointer: SourceConstPointer<i8>,
    ) -> Result<Vec<i8>, SourceHeapError> {
        let bytes = heap.slice(pointer)?;
        let length = bytes
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        Ok(bytes[..length].to_vec())
    }

    fn write_c_bytes(
        heap: &mut SourceHeap,
        output: SourceMutPointer<i8>,
        bytes: &[i8],
    ) -> Result<(), SourceHeapError> {
        let target = heap.slice_mut(output)?;
        let required = bytes
            .len()
            .checked_add(1)
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        if target.len() < required {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        target[..bytes.len()].copy_from_slice(bytes);
        target[bytes.len()] = 0;
        Ok(())
    }

    let max_tag_num = if tag_flag != 0 {
        crate::source_types::local_ichiprt1::tagIdentLblOrd_IL_MAX_ORD as usize
    } else {
        crate::source_types::local_ichiprt1::tagAuxLblOrd_AL_MAX_ORD as usize
    };
    if tags.len() < max_tag_num {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    if 0 < n_tag && n_tag < 3 {
        let mut selected = None;
        let mut bit = 1_i32;
        for index in 0..max_tag_num {
            if tag_bits & bit != 0 {
                selected = Some(index);
            }
            bit = bit.wrapping_shl(1);
        }
        if let Some(index) = selected {
            let source = if n_tag == 1 {
                tags[index].szXmlLabel
            } else {
                tags[index].szPlainLabel
            };
            let bytes = c_bytes(heap, source)?;
            write_c_bytes(heap, output, &bytes)?;
            if n_tag != 2 {
                *always = tags[index].bAlwaysOutput;
            }
            return Ok(output);
        }
    } else if n_tag == 3 {
        let mut bytes = vec![b'{' as i8];
        let mut selected = None;
        let mut bit = 1_i32;
        for (index, tag) in tags.iter().enumerate().take(max_tag_num) {
            if tag_bits & bit != 0 {
                if selected.is_some() {
                    bytes.push(b':' as i8);
                }
                bytes.extend(c_bytes(heap, tag.szPlainComment)?);
                selected = Some(index);
            }
            bit = bit.wrapping_shl(1);
        }
        if let Some(index) = selected {
            bytes.push(b'}' as i8);
            let label = c_bytes(heap, tags[index].szPlainLabel)?;
            let mut complete = label;
            complete.extend(bytes);
            write_c_bytes(heap, output, &complete)?;
            *always = tags[index].bAlwaysOutput;
        } else {
            write_c_bytes(heap, output, &[b'?' as i8, b'?' as i8, b'?' as i8])?;
        }
        return Ok(output);
    }

    write_c_bytes(heap, output, &[b'?' as i8, b'?' as i8, b'?' as i8])?;
    Ok(output)
}

#[allow(non_snake_case)]
pub(crate) fn str_LineEnd(
    heap: &mut SourceHeap,
    tag: SourceConstPointer<i8>,
    overflow: &mut i32,
    buffer: &mut INCHI_IOS_STRING,
    ind: i32,
    plain_text_tags: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2183 str_LineEnd
    // INCHI✔️❌: int str_LineEnd( const char       *tag,
    // INCHI✔️❌:                  int              *bOverflow,
    // INCHI✔️❌:                  INCHI_IOS_STRING *buf,
    // INCHI✔️❌:                  int               ind,
    // INCHI✔️❌:                  int               bPlainTextTags )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int tag_len;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* check buffer overflow */
    // INCHI✔️❌:     if (*bOverflow)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ind < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Plain text */
    // INCHI✔️❌:         /* insert plain text tag if:
    // INCHI✔️❌:            (a) pStr has non-zero length, or
    // INCHI✔️❌:            (b) ind < -1
    // INCHI✔️❌:         */
    // INCHI✔️❌:         if (buf->pStr[0] || ind < -1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             tag_len = bPlainTextTags ? (int) strlen( tag ) : 0;
    // INCHI✔️❌:             if (tag_len > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int n_added = tag_len + 2 + 2;
    // INCHI✔️❌:                 inchi_strbuf_update( buf, n_added );
    // INCHI✔️❌:
    // INCHI✔️❌:                 memmove(buf->pStr + tag_len, buf->pStr, (long long)buf->nUsedLength + 1); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                 /* NB: trailing 0 is also memmoved */
    // INCHI✔️❌:                 memcpy(buf->pStr, tag, tag_len);
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* to be sure...  */
    // INCHI✔️❌:                 buf->nUsedLength = strlen( buf->pStr );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: str_LineEnd

    if *overflow != 0 {
        return Ok(1);
    }
    if ind < 0 {
        let first = *heap
            .slice(buffer.pStr.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if first != 0 || ind < -1 {
            let tag_bytes = heap.slice(tag)?;
            let tag_length = if plain_text_tags != 0 {
                tag_bytes
                    .iter()
                    .position(|byte| *byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?
            } else {
                0
            };
            let tag_bytes = tag_bytes[..tag_length].to_vec();
            if tag_length > 0 {
                let addition = i32::try_from(tag_length)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                    .wrapping_add(4);
                if inchi_strbuf_update(heap, Some(buffer), addition)? < 0 {
                    return Err(SourceHeapError::AllocationFailed);
                }
                let used = usize::try_from(buffer.nUsedLength)
                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                let source = heap
                    .slice(buffer.pStr.as_const())?
                    .get(..=used)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .to_vec();
                let target = heap.slice_mut(buffer.pStr)?;
                let end = tag_length
                    .checked_add(source.len())
                    .ok_or(SourceHeapError::AllocationSizeOverflow)?;
                target
                    .get_mut(tag_length..end)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(&source);
                target[..tag_length].copy_from_slice(&tag_bytes);
                buffer.nUsedLength = i32::try_from(
                    target
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or(SourceHeapError::MissingNulTerminator)?,
                )
                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn MergeZzInStrHillFormulaComponent(
    heap: &mut SourceHeap,
    string: SourceMutPointer<i8>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5495 MergeZzInStrHillFormulaComponent
    // INCHI✔️❌: void MergeZzInStrHillFormulaComponent(char *s)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char *pz = strstr(s, "Zz");
    // INCHI✔️❌:     if (pz)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int n = 1, offset;
    // INCHI✔️❌:         char *pz2 = pz + 2, *pd = pz + 2;
    // INCHI✔️❌:         if (*pd && (isdigit(*pd)))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             n = strtol(pd, &pz2, 10);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pz2 = strstr(pz2, "Zz");
    // INCHI✔️❌:         if (pz2)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int n2 = 1;
    // INCHI✔️❌:             char *pd2 = pz2 + 2;
    // INCHI✔️❌:             if (*pd2 && (isdigit(*pd2)))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 n2 = strtol(pd2, &pz2, 10);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             n += n2;
    // INCHI✔️❌:             offset = (int)(pd - s);
    // INCHI✔️❌:             sprintf(s + offset, "%d", n);
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MergeZzInStrHillFormulaComponent

    let bytes = heap.slice(string.as_const())?;
    let nul = bytes
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let Some(first_z) = bytes[..nul]
        .windows(2)
        .position(|window| window == [b'Z' as i8, b'z' as i8])
    else {
        return Ok(());
    };
    let digit_position = first_z + 2;
    let mut next_search = digit_position;
    let mut count = 1_i32;
    if digit_position < nul && (bytes[digit_position] as u8).is_ascii_digit() {
        let start = string.offset(
            i64::try_from(digit_position).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
        )?;
        let mut end = start.as_const();
        count = inchi_strtol(heap, start.as_const(), Some(&mut end), 10)? as i32;
        next_search = usize::try_from(end.as_mut().difference(string)?)
            .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?;
    }
    let bytes = heap.slice(string.as_const())?;
    let Some(relative_second) = bytes[next_search..nul]
        .windows(2)
        .position(|window| window == [b'Z' as i8, b'z' as i8])
    else {
        return Ok(());
    };
    let second_z = next_search + relative_second;
    let second_digit = second_z + 2;
    let mut second_count = 1_i32;
    if second_digit < nul && (bytes[second_digit] as u8).is_ascii_digit() {
        let start = string.offset(
            i64::try_from(second_digit).map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
        )?;
        second_count = inchi_strtol(heap, start.as_const(), None, 10)? as i32;
    }
    count = count.wrapping_add(second_count);
    let rendered = count.to_string();
    let target = heap.slice_mut(string)?;
    let end = digit_position
        .checked_add(rendered.len())
        .and_then(|value| value.checked_add(1))
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let destination = target
        .get_mut(digit_position..end)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    for (target, source) in destination.iter_mut().zip(rendered.bytes()) {
        *target = source as i8;
    }
    destination[rendered.len()] = 0;
    Ok(())
}

#[allow(non_snake_case)]
pub(crate) fn MergeZzInHillFormula(
    heap: &mut SourceHeap,
    string_buffer: &mut INCHI_IOS_STRING,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5435 MergeZzInHillFormula
    // INCHI✔️❌: int MergeZzInHillFormula(INCHI_IOS_STRING *strbuf)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char *p, *scopy = NULL, *stmp=NULL, *pend=NULL, *p0 = NULL; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     size_t sublen; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!strbuf->pStr || strbuf->nUsedLength < 1)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     scopy = (char *)inchi_calloc((long long)strbuf->nAllocatedLength+1, sizeof(char)); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     if (!scopy)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(scopy); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         return -1; /* failed */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     memcpy(scopy, strbuf->pStr, strbuf->nAllocatedLength);
    // INCHI✔️❌:     stmp = (char *)inchi_calloc((long long)strbuf->nAllocatedLength + 1, sizeof(char)); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     if (!stmp)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(scopy); /* djb-rwth: avoiding memory leak */
    // INCHI✔️❌:         return -1; /* failed */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_reset(strbuf);
    // INCHI✔️❌:     p0 = scopy;
    // INCHI✔️❌:     p = p0;
    // INCHI✔️❌:     do
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         pend = strchr(p, '.');
    // INCHI✔️❌:         if (!pend)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             pend = strchr(p, '\0');
    // INCHI✔️❌:         }
    // INCHI✔️❌:         sublen = pend - p;
    // INCHI✔️❌:         memcpy(stmp, p, sublen);
    // INCHI✔️❌:         stmp[sublen] = '\0';
    // INCHI✔️❌:         MergeZzInStrHillFormulaComponent(stmp);
    // INCHI✔️❌:         if (stmp)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_strbuf_printf(strbuf, "%-s%-c", stmp, *pend);
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     } while ( *pend && (p=pend+1));
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (scopy)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(scopy);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (stmp)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_free(stmp);
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MergeZzInHillFormula

    if string_buffer.pStr.is_null() || string_buffer.nUsedLength < 1 {
        return Ok(0);
    }
    let allocated = usize::try_from(string_buffer.nAllocatedLength)
        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let temporary_length = u64::try_from(allocated)
        .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
        .checked_add(1)
        .ok_or(SourceHeapError::AllocationSizeOverflow)?;
    let copy = match inchi_calloc::<i8>(heap, temporary_length, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => return Ok(-1),
        Err(error) => return Err(error),
    };
    let original = heap
        .slice(string_buffer.pStr.as_const())?
        .get(..allocated)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    heap.slice_mut(copy)?[..allocated].copy_from_slice(&original);
    let temporary = match inchi_calloc::<i8>(heap, temporary_length, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed) => {
            inchi_free(heap, copy)?;
            return Ok(-1);
        }
        Err(error) => {
            inchi_free(heap, copy)?;
            return Err(error);
        }
    };

    inchi_strbuf_reset(heap, Some(string_buffer))?;
    let format =
        heap.allocate_model_storage(b"%-s%-c\0".iter().map(|byte| *byte as i8).collect())?;
    let mut position = 0_usize;
    loop {
        let copy_bytes = heap.slice(copy.as_const())?;
        let end = copy_bytes[position..]
            .iter()
            .position(|byte| matches!(*byte as u8, b'.' | 0))
            .map(|relative| position + relative)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let delimiter = copy_bytes[end];
        let component = copy_bytes[position..end].to_vec();
        let temporary_bytes = heap.slice_mut(temporary)?;
        temporary_bytes[..component.len()].copy_from_slice(&component);
        temporary_bytes[component.len()] = 0;
        MergeZzInStrHillFormulaComponent(heap, temporary)?;
        inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(temporary.as_const()),
                    SourceFormatArgument::Byte(delimiter),
                ],
                ..SourceVaList::default()
            },
        )?;
        if delimiter == 0 {
            break;
        }
        position = end + 1;
    }
    inchi_free(heap, copy)?;
    inchi_free(heap, temporary)?;
    Ok(0)
}

fn formula_ident_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
    // BEGIN INCHI C GLOBAL ENTRIES: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:297 IdentLbl[IL_FML__ORD]
    // INCHI✔️❌:                                                                        /* IL_FML__ORD, */    { "/",   "formula",        "formula",        1 }, /* basic part formula */
    // END INCHI C GLOBAL ENTRIES: IdentLbl[IL_FML__ORD]
    // BEGIN INCHI C GLOBAL ENTRIES: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:315 IdentLbl[IL_REC__ORD]
    // INCHI✔️❌:                                                                                                                                                                                                                                                                                                                                                        /* IL_REC__ORD, */    { "/r",  "reconnected bond(s) to metal(s) formula",  "formula",  0 }
    // END INCHI C GLOBAL ENTRIES: IdentLbl[IL_REC__ORD]

    let allocate = |heap: &mut SourceHeap, text: &[u8]| {
        heap.allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())
            .map(SourceMutPointer::as_const)
    };
    let mut tags = vec![INCHI_TAG::default(); 19];
    tags[4] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/\0")?,
        szPlainComment: allocate(heap, b"formula\0")?,
        szXmlLabel: allocate(heap, b"formula\0")?,
        bAlwaysOutput: 1,
    };
    tags[18] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/r\0")?,
        szPlainComment: allocate(heap, b"reconnected bond(s) to metal(s) formula\0")?,
        szXmlLabel: allocate(heap, b"formula\0")?,
        bAlwaysOutput: 0,
    };
    Ok(tags)
}

fn connections_ident_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
    // BEGIN INCHI C GLOBAL ENTRY: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:298 IdentLbl[IL_CONN_ORD]
    // INCHI✔️❌:                                                                        /* IL_CONN_ORD, */    { "/c",  "connections",    "connections",    1 },
    // END INCHI C GLOBAL ENTRY: IdentLbl[IL_CONN_ORD]
    let allocate = |heap: &mut SourceHeap, text: &[u8]| {
        heap.allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())
            .map(SourceMutPointer::as_const)
    };
    let mut tags = vec![INCHI_TAG::default(); 19];
    tags[5] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/c\0")?,
        szPlainComment: allocate(heap, b"connections\0")?,
        szXmlLabel: allocate(heap, b"connections\0")?,
        bAlwaysOutput: 1,
    };
    Ok(tags)
}

fn hydrogens_ident_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
    // BEGIN INCHI C GLOBAL ENTRY: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:299 IdentLbl[IL_ALLH_ORD]
    // INCHI✔️❌:                                                                        /* IL_ALLH_ORD, */    { "/h",  "H_atoms",        "H",              1 },
    // END INCHI C GLOBAL ENTRY: IdentLbl[IL_ALLH_ORD]
    let allocate = |heap: &mut SourceHeap, text: &[u8]| {
        heap.allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())
            .map(SourceMutPointer::as_const)
    };
    let mut tags = vec![INCHI_TAG::default(); 19];
    tags[6] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/h\0")?,
        szPlainComment: allocate(heap, b"H_atoms\0")?,
        szXmlLabel: allocate(heap, b"H\0")?,
        bAlwaysOutput: 1,
    };
    Ok(tags)
}

fn charge_proton_ident_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
    // BEGIN INCHI C GLOBAL ENTRIES: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:290 IdentLbl
    // INCHI✔️❌:     /* IL_FIXH_ORD, */    { "/",   "fixed_H",        "fixed-H",        0 }, /* fixed H */
    // INCHI✔️❌:                                                                        /* IL_CHRG_ORD, */    { "/q",  "charge",         "charge",         1 },
    // INCHI✔️❌:                                                                        /* IL_PROT_ORD, */    { "/p",  "protons",        "protons",        0 },
    // END INCHI C GLOBAL ENTRIES: IdentLbl
    let allocate = |heap: &mut SourceHeap, text: &[u8]| {
        heap.allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())
            .map(SourceMutPointer::as_const)
    };
    let mut tags = vec![INCHI_TAG::default(); 19];
    tags[0] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/\0")?,
        szPlainComment: allocate(heap, b"fixed_H\0")?,
        szXmlLabel: allocate(heap, b"fixed-H\0")?,
        bAlwaysOutput: 0,
    };
    tags[7] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/q\0")?,
        szPlainComment: allocate(heap, b"charge\0")?,
        szXmlLabel: allocate(heap, b"charge\0")?,
        bAlwaysOutput: 1,
    };
    tags[8] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/p\0")?,
        szPlainComment: allocate(heap, b"protons\0")?,
        szXmlLabel: allocate(heap, b"protons\0")?,
        bAlwaysOutput: 0,
    };
    Ok(tags)
}

fn stereo_ident_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
    // BEGIN INCHI C GLOBAL ENTRIES: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:289 IdentLbl
    // INCHI✔️❌:     /* IL_FIXH_ORD, */    { "/",   "fixed_H",        "fixed-H",        0 }, /* fixed H */
    // INCHI✔️❌:     /* IL_STER_ORD, */    { "/",   "stereo",         "stereo",         0 }, /* stereo */
    // INCHI✔️❌:     /* IL_DBND_ORD, */    { "/b",  "dbond",          "dbond",          0 },
    // INCHI✔️❌:     /* IL_SP3S_ORD, */    { "/t",  "sp3",            "sp3",            0 },
    // INCHI✔️❌:     /* IL_INVS_ORD, */    { "/m",  "sp3:inverted",   "abs.inverted",   0 }, /* mirrored */
    // INCHI✔️❌:     /* IL_TYPS_ORD, */    { "/s",  "type (1=abs, 2=rel, 3=rac)", "type",           0 }, /* stereo type */
    // END INCHI C GLOBAL ENTRIES: IdentLbl
    let allocate = |heap: &mut SourceHeap, text: &[u8]| {
        heap.allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())
            .map(SourceMutPointer::as_const)
    };
    let mut tags = vec![INCHI_TAG::default(); 19];
    tags[0] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/\0")?,
        szPlainComment: allocate(heap, b"fixed_H\0")?,
        szXmlLabel: allocate(heap, b"fixed-H\0")?,
        bAlwaysOutput: 0,
    };
    tags[2] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/\0")?,
        szPlainComment: allocate(heap, b"stereo\0")?,
        szXmlLabel: allocate(heap, b"stereo\0")?,
        bAlwaysOutput: 0,
    };
    tags[9] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/b\0")?,
        szPlainComment: allocate(heap, b"dbond\0")?,
        szXmlLabel: allocate(heap, b"dbond\0")?,
        bAlwaysOutput: 0,
    };
    tags[10] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/t\0")?,
        szPlainComment: allocate(heap, b"sp3\0")?,
        szXmlLabel: allocate(heap, b"sp3\0")?,
        bAlwaysOutput: 0,
    };
    tags[11] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/m\0")?,
        szPlainComment: allocate(heap, b"sp3:inverted\0")?,
        szXmlLabel: allocate(heap, b"abs.inverted\0")?,
        bAlwaysOutput: 0,
    };
    tags[12] = INCHI_TAG {
        szPlainLabel: allocate(heap, b"/s\0")?,
        szPlainComment: allocate(heap, b"type (1=abs, 2=rel, 3=rac)\0")?,
        szXmlLabel: allocate(heap, b"type\0")?,
        bAlwaysOutput: 0,
    };
    Ok(tags)
}

fn isotopic_ident_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
    // BEGIN INCHI C GLOBAL ENTRIES: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:289 IdentLbl
    // INCHI✔️❌:     /* IL_FIXH_ORD, */    { "/",   "fixed_H",        "fixed-H",        0 }, /* fixed H */
    // INCHI✔️❌:     /* IL_ISOT_ORD, */    { "/",   "isotopic",       "isotopic",       0 }, /* isotopic */
    // INCHI✔️❌:     /* IL_STER_ORD, */    { "/",   "stereo",         "stereo",         0 }, /* stereo */
    // INCHI✔️❌:     /* IL_DBND_ORD, */    { "/b",  "dbond",          "dbond",          0 },
    // INCHI✔️❌:     /* IL_SP3S_ORD, */    { "/t",  "sp3",            "sp3",            0 },
    // INCHI✔️❌:     /* IL_INVS_ORD, */    { "/m",  "sp3:inverted",   "abs.inverted",   0 }, /* mirrored */
    // INCHI✔️❌:     /* IL_TYPS_ORD, */    { "/s",  "type (1=abs, 2=rel, 3=rac)", "type",           0 }, /* stereo type */
    // INCHI✔️❌:     /* IL_ATMS_ORD, */    { "/i",  "atoms",          "atoms",          1 },
    // INCHI✔️❌:     /* IL_XCGA_ORD, */    { "/h",  "exchangeable_H", "H-isotopic",     1 },
    // INCHI✔️❌:     /* IL_FMLF_ORD, */    { "/f",  "formula",        "formula",        1 }, /* fixed H formula */
    // INCHI✔️❌:     /* IL_HFIX_ORD, */    { "/h",  "H_fixed" ,       "H-fixed" ,       1 }, /* fixed-H */
    // INCHI✔️❌:     /* IL_TRNS_ORD, */    { "/o",  "transposition",  "transposition",  0 }, /* order */
    // END INCHI C GLOBAL ENTRIES: IdentLbl
    let allocate = |heap: &mut SourceHeap, text: &[u8]| {
        heap.allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())
            .map(SourceMutPointer::as_const)
    };
    let mut tags = vec![INCHI_TAG::default(); 19];
    for (index, plain, comment, xml, always) in [
        (
            0,
            b"/\0".as_slice(),
            b"fixed_H\0".as_slice(),
            b"fixed-H\0".as_slice(),
            0,
        ),
        (
            1,
            b"/\0".as_slice(),
            b"isotopic\0".as_slice(),
            b"isotopic\0".as_slice(),
            0,
        ),
        (
            2,
            b"/\0".as_slice(),
            b"stereo\0".as_slice(),
            b"stereo\0".as_slice(),
            0,
        ),
        (
            9,
            b"/b\0".as_slice(),
            b"dbond\0".as_slice(),
            b"dbond\0".as_slice(),
            0,
        ),
        (
            10,
            b"/t\0".as_slice(),
            b"sp3\0".as_slice(),
            b"sp3\0".as_slice(),
            0,
        ),
        (
            11,
            b"/m\0".as_slice(),
            b"sp3:inverted\0".as_slice(),
            b"abs.inverted\0".as_slice(),
            0,
        ),
        (
            12,
            b"/s\0".as_slice(),
            b"type (1=abs, 2=rel, 3=rac)\0".as_slice(),
            b"type\0".as_slice(),
            0,
        ),
        (
            13,
            b"/i\0".as_slice(),
            b"atoms\0".as_slice(),
            b"atoms\0".as_slice(),
            1,
        ),
        (
            14,
            b"/h\0".as_slice(),
            b"exchangeable_H\0".as_slice(),
            b"H-isotopic\0".as_slice(),
            1,
        ),
        (
            15,
            b"/f\0".as_slice(),
            b"formula\0".as_slice(),
            b"formula\0".as_slice(),
            1,
        ),
        (
            16,
            b"/h\0".as_slice(),
            b"H_fixed\0".as_slice(),
            b"H-fixed\0".as_slice(),
            1,
        ),
        (
            17,
            b"/o\0".as_slice(),
            b"transposition\0".as_slice(),
            b"transposition\0".as_slice(),
            0,
        ),
    ] {
        tags[index] = INCHI_TAG {
            szPlainLabel: allocate(heap, plain)?,
            szPlainComment: allocate(heap, comment)?,
            szXmlLabel: allocate(heap, xml)?,
            bAlwaysOutput: always,
        };
    }
    Ok(tags)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINChI2(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut string_buffer: Option<&mut INCHI_IOS_STRING>,
    sort_rows: SourceMutPointer<SourceMutPointer<crate::source_types::INCHI_SORT>>,
    basic_or_reconnected: i32,
    original_atom_data: SourceConstPointer<ORIG_ATOM_DATA>,
    original_structure: SourceConstPointer<ORIG_STRUCT>,
    input_parameters: &crate::source_types::INPUT_PARMS,
    disconnected_coordinates: i32,
    output_type: i32,
    output_options: i32,
    component_counts: &[i32; 2],
    non_tautomeric_counts: &[i32; 2],
    tautomeric_counts: &[i32; 2],
    output: &mut INCHI_IOSTREAM,
    log: &mut INCHI_IOSTREAM,
    input_structure_number: i32,
    sort_print_flags: SourceMutPointer<i32>,
    save_option_bits: u8,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:966 OutputINChI2
    // INCHI✔️❌: int OutputINChI2( CANON_GLOBALS     *pCG,
    // INCHI✔️❌:                   INCHI_IOS_STRING  *strbuf,
    // INCHI✔️❌:                   INCHI_SORT        *pINChISortTautAndNonTaut2[][TAUT_NUM],
    // INCHI✔️❌:                   int               INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                   ORIG_ATOM_DATA    *orig_inp_data,
    // INCHI✔️❌:                   ORIG_STRUCT       *pOrigStruct,
    // INCHI✔️❌:                   INPUT_PARMS       *ip,
    // INCHI✔️❌:                   int               bDisconnectedCoord,
    // INCHI✔️❌:                   int               bOutputType,
    // INCHI✔️❌:                   int               bINChIOutputOptions,
    // INCHI✔️❌:                   int               num_components2[],
    // INCHI✔️❌:                   int               num_non_taut2[],
    // INCHI✔️❌:                   int               num_taut2[],
    // INCHI✔️❌:                   INCHI_IOSTREAM    *output_file,
    // INCHI✔️❌:                   INCHI_IOSTREAM    *log_file,
    // INCHI✔️❌:                   int               num_input_struct,
    // INCHI✔️❌:                   int               *pSortPrintINChIFlags,
    // INCHI✔️❌:                   unsigned char     save_opt_bits )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int bINChIOutputOptions0 = bINChIOutputOptions & ~( INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS );
    // INCHI✔️❌:     int bINChIOutputOptionsCur;
    // INCHI✔️❌:     int bCurOption, ret, i;
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < 3; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         switch (i)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             case 1:
    // INCHI✔️❌:                 bCurOption = INCHI_OUT_PLAIN_TEXT;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             case 2:
    // INCHI✔️❌:                 bCurOption = INCHI_OUT_PLAIN_TEXT_COMMENTS;
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             default:
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (bINChIOutputOptions & bCurOption)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             bINChIOutputOptionsCur = bINChIOutputOptions0 | bCurOption;
    // INCHI✔️❌:             if (i != 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bINChIOutputOptionsCur &= ~INCHI_OUT_TABBED_OUTPUT;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             ret |= OutputINChI1( pCG,
    // INCHI✔️❌:                                  strbuf,
    // INCHI✔️❌:                                  pINChISortTautAndNonTaut2,
    // INCHI✔️❌:                                  INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                  orig_inp_data,
    // INCHI✔️❌:                                  pOrigStruct,
    // INCHI✔️❌:                                  ip,
    // INCHI✔️❌:                                  bDisconnectedCoord,
    // INCHI✔️❌:                                  bOutputType,
    // INCHI✔️❌:                                  bINChIOutputOptionsCur,
    // INCHI✔️❌:                                  num_components2,
    // INCHI✔️❌:                                  num_non_taut2,
    // INCHI✔️❌:                                  num_taut2,
    // INCHI✔️❌:                                  output_file,
    // INCHI✔️❌:                                  log_file,
    // INCHI✔️❌:                                  num_input_struct,
    // INCHI✔️❌:                                  pSortPrintINChIFlags,
    // INCHI✔️❌:                                  save_opt_bits );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINChI2

    let base_options = output_options
        & !((crate::source_types::INCHI_OUT_PLAIN_TEXT
            | crate::source_types::INCHI_OUT_PLAIN_TEXT_COMMENTS) as i32);
    let mut result = 0_i32;
    for (index, current_option) in [
        crate::source_types::INCHI_OUT_PLAIN_TEXT as i32,
        crate::source_types::INCHI_OUT_PLAIN_TEXT_COMMENTS as i32,
    ]
    .into_iter()
    .enumerate()
    {
        if output_options & current_option == 0 {
            continue;
        }
        let mut current_options = base_options | current_option;
        if index != 0 {
            current_options &= !(crate::source_types::INCHI_OUT_TABBED_OUTPUT as i32);
        }
        result |= OutputINChI1(
            heap,
            canonical_globals,
            string_buffer.as_deref_mut(),
            sort_rows,
            basic_or_reconnected,
            original_atom_data,
            original_structure,
            input_parameters,
            disconnected_coordinates,
            output_type,
            current_options,
            component_counts,
            non_tautomeric_counts,
            tautomeric_counts,
            output,
            log,
            input_structure_number,
            sort_print_flags,
            save_option_bits,
            stdout,
        )?;
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINChI1(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    string_buffer: Option<&mut INCHI_IOS_STRING>,
    sort_rows: SourceMutPointer<SourceMutPointer<crate::source_types::INCHI_SORT>>,
    mut basic_or_reconnected: i32,
    original_atom_data: SourceConstPointer<ORIG_ATOM_DATA>,
    original_structure: SourceConstPointer<ORIG_STRUCT>,
    input_parameters: &crate::source_types::INPUT_PARMS,
    disconnected_coordinates: i32,
    output_type: i32,
    output_options: i32,
    component_counts: &[i32; 2],
    non_tautomeric_counts: &[i32; 2],
    tautomeric_counts: &[i32; 2],
    output: &mut INCHI_IOSTREAM,
    log: &mut INCHI_IOSTREAM,
    input_structure_number: i32,
    sort_print_flags: SourceMutPointer<i32>,
    save_option_bits: u8,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:1043 OutputINChI1
    // INCHI✔️❌: complete configured source frame follows verbatim; the second axis remains
    // INCHI✔️❌: open because SourceHeap lookup and temporary clones add material overhead.
    // INCHI✔️❌: int OutputINChI1( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                   INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                   INCHI_SORT *pINChISortTautAndNonTaut2[][TAUT_NUM],
    // INCHI✔️❌:                   int INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                   ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                   ORIG_STRUCT *pOrigStruct,
    // INCHI✔️❌:                   INPUT_PARMS *ip,
    // INCHI✔️❌:                   int bDisconnectedCoord,
    // INCHI✔️❌:                   int bOutputType,
    // INCHI✔️❌:                   int bINChIOutputOptions,
    // INCHI✔️❌:                   int num_components2[],
    // INCHI✔️❌:                   int num_non_taut2[],
    // INCHI✔️❌:                   int num_taut2[],
    // INCHI✔️❌:                   INCHI_IOSTREAM *out_file,
    // INCHI✔️❌:                   INCHI_IOSTREAM *log_file,
    // INCHI✔️❌:                   int num_input_struct,
    // INCHI✔️❌:                   int *pSortPrintINChIFlags,
    // INCHI✔️❌:                   unsigned char save_opt_bits )
    // INCHI✔️❌: {
    // INCHI✔️❌:     INCHI_OUT_CTL io;
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         bINChIOutputOptions bits:
    // INCHI✔️❌:           INCHI_OUT_NO_AUX_INFO           0x0001    do not output Aux Info
    // INCHI✔️❌:           INCHI_OUT_SHORT_AUX_INFO        0x0002    output short version of Aux Info
    // INCHI✔️❌:           INCHI_OUT_ONLY_AUX_INFO         0x0004    output only Aux Info
    // INCHI✔️❌:           INCHI_OUT_EMBED_REC             0x0008    embed reconnected INChI into disconnected INChI
    // INCHI✔️❌:
    // INCHI✔️❌:         bOutputType =
    // INCHI✔️❌:          TAUT_YES  => tautomeric only (if no tautomeric components then no output;
    // INCHI✔️❌:          TAUT_NON  => only non-tautomeric output (if no non-taut present then no output;
    // INCHI✔️❌:          TAUT_BOTH => tautomeric and non-tautomeric
    // INCHI✔️❌:     */
    // INCHI✔️❌:     int  i, j, ii, jj, /*ii2, jj2,*/ bEmbeddedOutputCalled = 0;
    // INCHI✔️❌:     int  bTautIsoHNum, bTautIsoAt, bHasIsotopicAtoms[TAUT_NUM];
    // INCHI✔️❌:     int  bStereoSp2[TAUT_NUM], bStereoSp3[TAUT_NUM];
    // INCHI✔️❌:     int  bIsotopicStereoSp2[TAUT_NUM], bIsotopicStereoSp3[TAUT_NUM];
    // INCHI✔️❌:     int  bStereoAbsInverted[TAUT_NUM], bIsotopicStereoAbsInverted[TAUT_NUM];
    // INCHI✔️❌:     int  bStereoAbs[TAUT_NUM], bIsotopicStereoAbs[TAUT_NUM];
    // INCHI✔️❌:     int  bTautomericAcid, bHardAddRemProton;
    // INCHI✔️❌:     int  bRequestedRacemicStereo = 0, bRequestedRelativeStereo = 0;
    // INCHI✔️❌:     int  npass = 0; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:
    // INCHI✔️❌:     INCHI_SORT   *is, *is2;
    // INCHI✔️❌:     INChI        *pINChI /*, *pINChI2*/;
    // INCHI✔️❌:     INChI_Aux    *pINChI_Aux = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     int  ret = 0;        /*  0 failed, 1 success */
    // INCHI✔️❌:     int  intermediate_result = 0;
    // INCHI✔️❌:     int then_goto_repeat = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int  bHasIsoH;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:     int  bTautAndNonTaut, bTautIsNonTaut;
    // INCHI✔️❌:     int nAtomsAllComp1, nAtomsAllComp2;    /* v. 1.05 Total atoms in all components */
    // INCHI✔️❌:
    // INCHI✔️❌:     int  bPlainText = 0 != ( bINChIOutputOptions & ( INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS ) );
    // INCHI✔️❌:     int  bPlainTextCommnts = 0 != ( bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT_COMMENTS );
    // INCHI✔️❌:
    // INCHI✔️❌:     char *pLF, *pTAB;
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:     int silent = 1;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     int bFixTranspChargeBug = 0;
    // INCHI✔️❌: #if ( FIX_TRANSPOSITION_CHARGE_BUG == 1 ) /* 2008-01-02 */
    // INCHI✔️❌:     if (INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG & bINChIOutputOptions)
    // INCHI✔️❌:         bFixTranspChargeBug = 1;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     io.bAbcNumbers = ip->bAbcNumbers;
    // INCHI✔️❌:
    // INCHI✔️❌:     io.ATOM_MODE = ( ( io.bAbcNumbers ? CT_MODE_ABC_NUMBERS : 0 )
    // INCHI✔️❌:                     | CT_MODE_ATOM_COUNTS
    // INCHI✔️❌:                     | CT_MODE_NO_ORPHANS
    // INCHI✔️❌: #if ( EQL_H_NUM_TOGETHER == 1 )
    // INCHI✔️❌:                     | CT_MODE_EQL_H_TOGETHER
    // INCHI✔️❌: #endif
    // INCHI✔️❌: #if ( ABC_CT_NUM_CLOSURES == 1 )
    // INCHI✔️❌:                     | ( io.bAbcNumbers && ip->bCtPredecessors ? CT_MODE_ABC_NUM_CLOSURES : 0 )
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                     | ( ip->bCtPredecessors ? CT_MODE_PREDECESSORS : 0 ) );
    // INCHI✔️❌:
    // INCHI✔️❌:     io.TAUT_MODE = ( io.bAbcNumbers ? CT_MODE_ABC_NUMBERS : 0 );
    // INCHI✔️❌:     io.pSortPrintINChIFlags = pSortPrintINChIFlags;
    // INCHI✔️❌:     io.num_components = num_components2[INCHI_basic_or_INCHI_reconnected];
    // INCHI✔️❌:     io.pINChISortTautAndNonTaut = pINChISortTautAndNonTaut2[INCHI_basic_or_INCHI_reconnected];
    // INCHI✔️❌:     io.pINChISort = io.pINChISortTautAndNonTaut[TAUT_YES];
    // INCHI✔️❌:     io.pINChISort2 = io.pINChISortTautAndNonTaut[TAUT_YES];
    // INCHI✔️❌:     io.bAlways = 0;
    // INCHI✔️❌:     io.bUseMulipliers = 1;
    // INCHI✔️❌:     io.bOmitRepetitions = 1;
    // INCHI✔️❌:     io.bPlainTextTags = 2;  /* 0 => no plain tags, 1=> plain text tags, 2=>plaintext tags without consecutive // */
    // INCHI✔️❌:     io.bOutputType = bOutputType; /* remains constant */
    // INCHI✔️❌:     io.bOutType = bOutputType;    /* will change! */
    // INCHI✔️❌:     io.bOverflow = 0;
    // INCHI✔️❌:     io.bSecondNonTautPass = 0;
    // INCHI✔️❌:     io.bNonTautIsoIdentifierNotEmpty = 0;
    // INCHI✔️❌:     io.bNonTautNonIsoIdentifierNotEmpty = 0;
    // INCHI✔️❌:     io.bNonTautIsIdenticalToTaut = 1;
    // INCHI✔️❌:     io.bFhTag = 0;
    // INCHI✔️❌:     io.nTag = bPlainTextCommnts ? 3 : bPlainText ? 2 : 0; /* tag type */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     if (NULL==orig_inp_data)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*intermediate_result = 1;
    // INCHI✔️❌:         goto exit_function;*/
    // INCHI✔️❌:         io.n_zy     = 0;
    // INCHI✔️❌:         io.n_pzz    = 0;
    // INCHI✔️❌:         io.n_pzz    = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         io.n_zy = orig_inp_data->n_zy;
    // INCHI✔️❌:         io.n_pzz = 0;
    // INCHI✔️❌:         if (orig_inp_data->polymer)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             io.n_pzz = orig_inp_data->polymer->n_pzz;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     io.bPolymers = ip->bPolymers;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:     /* @@@ From now on we will go with silent output; it ends on @@@ below */
    // INCHI✔️❌:     silent = 1;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Analyze layers, make adjustments and fixes, etc. */
    // INCHI✔️❌:
    // INCHI✔️❌:     set_line_separators( bINChIOutputOptions, &pLF, &pTAB );
    // INCHI✔️❌:     memset( io.sDifSegs, DIFV_BOTH_EMPTY, sizeof( io.sDifSegs ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     if (!strbuf || !( strbuf->pStr ) || strbuf->nAllocatedLength <= 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "Cannot allocate output buffer. No output for structure #%d.%s%s%s%s\n",
    // INCHI✔️❌:                          num_input_struct, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:         return ret;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    -- commented out to allow empty InChI --
    // INCHI✔️❌:     if (!io.num_components ) return 0;
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < TAUT_NUM; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         bHasIsotopicAtoms[i] =
    // INCHI✔️❌:             io.num_comp[i] =
    // INCHI✔️❌:             bStereoSp2[i] =
    // INCHI✔️❌:             bStereoSp3[i] =
    // INCHI✔️❌:             bIsotopicStereoSp2[i] =
    // INCHI✔️❌:             bIsotopicStereoSp3[i] =
    // INCHI✔️❌:             io.bIsotopicOrigNumb[i] =
    // INCHI✔️❌:             bStereoAbs[i] =
    // INCHI✔️❌:             bIsotopicStereoAbs[i] =
    // INCHI✔️❌:             bStereoAbsInverted[i] =
    // INCHI✔️❌:             bIsotopicStereoAbsInverted[i] =
    // INCHI✔️❌:             io.bRacemicStereo[i] =
    // INCHI✔️❌:             io.bRelativeStereo[i] =
    // INCHI✔️❌:             io.bIsotopicRacemicStereo[i] =
    // INCHI✔️❌:             io.bIsotopicRelativeStereo[i] =
    // INCHI✔️❌:             io.bAtomEqu[i] =
    // INCHI✔️❌:             io.bTautEqu[i] =
    // INCHI✔️❌:             io.bIsotopicAtomEqu[i] =
    // INCHI✔️❌:             io.bIsotopicTautEqu[i] =
    // INCHI✔️❌:             io.bInvStereo[i] =
    // INCHI✔️❌:             io.bInvIsotopicStereo[i] =
    // INCHI✔️❌:             io.bInvStereoOrigNumb[i] =
    // INCHI✔️❌:             io.bInvIsotopicStereoOrigNumb[i] =
    // INCHI✔️❌:             io.bIgn_UU_Sp3[i] =
    // INCHI✔️❌:             io.bIgn_UU_Sp2[i] =
    // INCHI✔️❌:             io.bIgn_UU_Sp3_Iso[i] =
    // INCHI✔️❌:             io.bIgn_UU_Sp2_Iso[i] =
    // INCHI✔️❌:             io.bChargesRadVal[i] =
    // INCHI✔️❌:             io.bOrigCoord[i] = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*    Find if it is isotopic */
    // INCHI✔️❌:         io.bIsotopic =
    // INCHI✔️❌:         io.bTautomeric =
    // INCHI✔️❌:         io.bNonTautomeric =
    // INCHI✔️❌:         bTautomericAcid =
    // INCHI✔️❌:         bHardAddRemProton =
    // INCHI✔️❌:         bTautIsoHNum =
    // INCHI✔️❌:         bTautIsoAt =
    // INCHI✔️❌:         bTautAndNonTaut =
    // INCHI✔️❌:         bTautIsNonTaut = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         /*
    // INCHI✔️❌:              x = bStereo, bStereoSp2, bStereoSp3, bStereoAbsInverted,
    // INCHI✔️❌:                  bIsotopicStereo, bIsotopicStereoSp2, bIsotopicStereoSp3, bIsotopicStereoAbsInverted
    // INCHI✔️❌:
    // INCHI✔️❌:              OUT_N1: x[TAUT_NON] refers to non-tautomeric only
    // INCHI✔️❌:              OUT_T1: x[TAUT_YES] refers to tautomeric if exists otherwise non-tautomeric
    // INCHI✔️❌:              OUT_NT: x[TAUT_NON] refers to non-taut representations of tautomeric
    // INCHI✔️❌:              OUT_TN: x[TAUT_YES] refers to tautomeric if exists otherwise non-tautomeric
    // INCHI✔️❌:                      x[TAUT_NON] refers to non-taut representations of tautomeric
    // INCHI✔️❌:          */
    // INCHI✔️❌:
    // INCHI✔️❌:     memset( io.num_iso_H, 0, sizeof( io.num_iso_H ) ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     io.nNumRemovedProtons = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     bHasIsoH = 0;
    // INCHI✔️❌:     io.bTautomericOutputAllowed = ( io.bOutType == OUT_T1 || io.bOutType == OUT_TN );
    // INCHI✔️❌:     io.pINChISort = io.pINChISortTautAndNonTaut[io.bTautomericOutputAllowed ? TAUT_YES : TAUT_NON];
    // INCHI✔️❌:     is = io.pINChISort;
    // INCHI✔️❌:     /* djb-rwth: removing redundant variables/code */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0, is2 = io.pINChISortTautAndNonTaut[TAUT_NON]; i < io.num_components; i++, is++, is2 ? is2++ : NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         CompINChILayers( is, is2, io.sDifSegs, bFixTranspChargeBug );
    // INCHI✔️❌:
    // INCHI✔️❌:         io.bNonTautIsIdenticalToTaut = io.bNonTautIsIdenticalToTaut && !CompINChITautVsNonTaut( is, is2, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (is && ( pINChI_Aux = is->pINChI_Aux[TAUT_YES] ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = 0; j < NUM_H_ISOTOPES; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bHasIsoH += abs( pINChI_Aux->nNumRemovedIsotopicH[j] );
    // INCHI✔️❌:                 io.num_iso_H[j] += pINChI_Aux->nNumRemovedIsotopicH[j];
    // INCHI✔️❌:             }
    // INCHI✔️❌:             io.nNumRemovedProtons += pINChI_Aux->nNumRemovedProtons;
    // INCHI✔️❌:             /* djb-rwth: removing redundant code */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (io.bTautomericOutputAllowed && is) /* djb-rwth: fixing a NULL pointer dereference */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Check for removed isotopic H */
    // INCHI✔️❌:             for (j = TAUT_YES; j < TAUT_NUM; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 switch (io.bOutType)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     case OUT_N1: /* x[TAUT_NON]: non-tautomeric only -- never happens */
    // INCHI✔️❌:                         jj = GET_II( io.bOutType, is );
    // INCHI✔️❌:                         if (jj != j)
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case OUT_T1: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric */
    // INCHI✔️❌:                         jj = GET_II( io.bOutType, is );
    // INCHI✔️❌:                         if (jj != j)
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case OUT_NT: /* x[TAUT_NON]: only non-taut representations of tautomeric -- never happens */
    // INCHI✔️❌:                         jj = GET_II( io.bOutType, is );
    // INCHI✔️❌:                         if (jj != j)
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     /* main path of control flow */
    // INCHI✔️❌:                     case OUT_TN: /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric;
    // INCHI✔️❌:                                   * x[TAUT_NON]: non-taut only if tautomeric is present */
    // INCHI✔️❌:                         jj = ( j == TAUT_YES ) ? GET_II( OUT_T1, is ) : ( j == TAUT_NON ) ? GET_II( OUT_NT, is ) : -1;
    // INCHI✔️❌:                         if (jj == TAUT_YES)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* Fix12 */
    // INCHI✔️❌:                             if (is->pINChI[jj]->lenTautomer > 0)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 bTautAndNonTaut += ( !is->pINChI[jj]->bDeleted && HAS_N( is ) );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 bTautIsNonTaut++;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (jj < 0)
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     default:
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (jj != j)
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:                 if (( pINChI = is->pINChI[jj] ) && pINChI->nNumberOfAtoms > 0 && ( pINChI_Aux = is->pINChI_Aux[jj] ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bTautIsoHNum += ( pINChI_Aux->nNumRemovedIsotopicH[0] +
    // INCHI✔️❌:                                      pINChI_Aux->nNumRemovedIsotopicH[1] +
    // INCHI✔️❌:                                      pINChI_Aux->nNumRemovedIsotopicH[2] );
    // INCHI✔️❌:                     bTautIsoAt += ( pINChI->nNumberOfIsotopicAtoms > 0 || pINChI->nNumberOfIsotopicTGroups > 0 );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     io.sDifSegs[DIFL_M][DIFS_p_PROTONS] = io.nNumRemovedProtons ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:     io.sDifSegs[DIFL_MI][DIFS_h_H_ATOMS] = bHasIsoH ? DIFV_NEQ2PRECED : DIFV_BOTH_EMPTY;
    // INCHI✔️❌:
    // INCHI✔️❌:     MarkUnusedAndEmptyLayers( io.sDifSegs );
    // INCHI✔️❌:
    // INCHI✔️❌:     io.bNonTautIsIdenticalToTaut = io.bNonTautIsIdenticalToTaut && !bTautIsoHNum;
    // INCHI✔️❌:     nAtomsAllComp1 = nAtomsAllComp2 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0, is = io.pINChISort; i < io.num_components; i++, is++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int bCurIso, bCurHasIsoStereo /* Fix14 */, bCurTaut /*, bCurTaut2*/; /* djb-rwth: removing redundant variables */
    // INCHI✔️❌:         int bCompExists, bCurIsoHPos, bCurIsoHStereo;
    // INCHI✔️❌:         int bCurStereoSp2, bCurIsoStereoSp2, bCurStereoSp3, bCurIsoStereoSp3, bCurIsoStereoSp3Inv;
    // INCHI✔️❌:         int bCurRacemic, bCurRelative, bCurIsoRacemic, bCurIsoRelative;
    // INCHI✔️❌:         bCompExists = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:         for (j = TAUT_NON; j < TAUT_NUM; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             switch (io.bOutType)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 case OUT_N1:
    // INCHI✔️❌:                     /* x[TAUT_NON]: non-tautomeric only */
    // INCHI✔️❌:                     jj = GET_II( io.bOutType, is );
    // INCHI✔️❌:                     if (jj != j)
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     ii = TAUT_NON;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case OUT_T1:
    // INCHI✔️❌:                     /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric */
    // INCHI✔️❌:                     jj = GET_II( io.bOutType, is );
    // INCHI✔️❌:                     if (jj != j)
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     ii = TAUT_YES;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 case OUT_NT:
    // INCHI✔️❌:                     /* x[TAUT_NON]: only non-taut representations of tautomeric */
    // INCHI✔️❌:                     jj = GET_II( io.bOutType, is );
    // INCHI✔️❌:                     if (jj != j)
    // INCHI✔️❌:                         continue;
    // INCHI✔️❌:                     ii = TAUT_NON;
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:                 /* main control flow comes here: requested both mobile and fixed H results */
    // INCHI✔️❌:                 case OUT_TN:
    // INCHI✔️❌:                     /* x[TAUT_YES]: tautomeric if present otherwise non-tautomeric; */
    // INCHI✔️❌:                     /* x[TAUT_NON]: non-taut only if tautomeric is present          */
    // INCHI✔️❌:                     jj = ( j == TAUT_YES ) ? GET_II( OUT_T1, is ) : ( j == TAUT_NON ) ? GET_II( OUT_NT, is ) : -1;
    // INCHI✔️❌:                     if (jj < 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* Fix12 */
    // INCHI✔️❌:                         if (bTautAndNonTaut && bTautIsNonTaut &&
    // INCHI✔️❌:                              j == TAUT_NON && 0 <= ( jj = GET_II( OUT_T1, is ) ) &&
    // INCHI✔️❌:                              !is->pINChI[jj]->bDeleted && !is->pINChI[jj]->lenTautomer)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             ; /* the requested non-tautomeric component is in tautomeric position   */
    // INCHI✔️❌:                               /*   (is->pINChI[TAUT_YES]);                                          */
    // INCHI✔️❌:                               /*   process it also as non-tautomeric if Fixed-H layer was requested */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             continue;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     ii = j; /* ii is what we wanted; jj is what we found (0 = TAUT_NON: fixed_H, 1 = TAUT_YES: mobile_H) */
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* -- not used 2004-09-16 ---
    // INCHI✔️❌:                     if ( is2 ) {
    // INCHI✔️❌:                         jj2 = ( j == TAUT_YES )? GET_II(OUT_T1,is2) : ( j == TAUT_NON )? GET_II(OUT_NT,is2) : -1;
    // INCHI✔️❌:                         if ( jj2 >= 0 ) {
    // INCHI✔️❌:                             ii2 = j;
    // INCHI✔️❌:                         } else {
    // INCHI✔️❌:                             ii2 = -1;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     } else {
    // INCHI✔️❌:                         jj2 = ii2 = -1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     -----------------------------*/
    // INCHI✔️❌:                     break;
    // INCHI✔️❌:
    // INCHI✔️❌:                 default:
    // INCHI✔️❌:                     continue;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (( pINChI = is->pINChI[jj] ) && pINChI->nNumberOfAtoms > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /*pINChI_Aux = is->pINChI_Aux[jj];*/
    // INCHI✔️❌:                 bCompExists++;
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (j == TAUT_NON)
    // INCHI✔️❌:                     nAtomsAllComp1 += pINChI->nNumberOfAtoms;
    // INCHI✔️❌:                 else if (j == TAUT_YES)
    // INCHI✔️❌:                     nAtomsAllComp2 += pINChI->nNumberOfAtoms;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                 bCurTaut = ( pINChI->lenTautomer > 0 );
    // INCHI✔️❌:                 bCurIso = ( pINChI->nNumberOfIsotopicAtoms > 0 || pINChI->nNumberOfIsotopicTGroups > 0 );
    // INCHI✔️❌:                 bCurIsoHPos = ( (pINChI->nPossibleLocationsOfIsotopicH && pINChI->nPossibleLocationsOfIsotopicH[0] > 1) || pINChI->lenTautomer > 1 ); /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 /* present isotopic H + their possible positions AND/OR isotopic atoms */
    // INCHI✔️❌:                 bCurIsoHStereo = (bCurIsoHPos && ( bTautIsoHNum || bTautIsoAt )) || bCurIso; /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                 if (jj == j && pINChI->bDeleted)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     io.num_comp[j] --;
    // INCHI✔️❌:                     if (bCurTaut)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         io.bTautomeric |= 1; /* tautomeric representation is present */
    // INCHI✔️❌:                         io.bNonTautomeric |= HAS_N( is );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     io.bIsotopic |= bCurIso;
    // INCHI✔️❌:                     continue; /* deleted H(+) in tautomeric representation */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 bCurStereoSp2 = pINChI->Stereo && ( pINChI->Stereo->nNumberOfStereoBonds > 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:                 bCurHasIsoStereo =
    // INCHI✔️❌:                     bCurStereoSp3 = pINChI->Stereo && ( pINChI->Stereo->nNumberOfStereoCenters > 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:                 bCurIsoStereoSp2 = bCurIsoHStereo && pINChI->StereoIsotopic && ( pINChI->StereoIsotopic->nNumberOfStereoBonds > 0 );
    // INCHI✔️❌:                 bCurIsoStereoSp3 = bCurIsoHStereo && pINChI->StereoIsotopic && ( pINChI->StereoIsotopic->nNumberOfStereoCenters > 0 );
    // INCHI✔️❌:                 bCurIsoStereoSp3Inv = bCurIsoStereoSp3 && pINChI->StereoIsotopic->nCompInv2Abs; /* inversion changes sp3 stereo */
    // INCHI✔️❌:                 bRequestedRacemicStereo |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_RAC_STEREO ) );
    // INCHI✔️❌:                 bRequestedRelativeStereo |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_REL_STEREO ) );
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Check whether isotopic stereo is same as non-isotopic; if same than do not output isotopic stereo */
    // INCHI✔️❌:                 if (bCurStereoSp2 && bCurIsoStereoSp2)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bCurIsoStereoSp2 = !Eql_INChI_Stereo( pINChI->Stereo, EQL_SP2, pINChI->StereoIsotopic, EQL_SP2, 0 );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (bCurStereoSp3 && bCurIsoStereoSp3)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* bCurIsoStereoSp3=0 means (iso stereo sp3) = (non-iso stereo sp3) or (iso stereo sp3) = Inv(non-iso stereo sp3) */
    // INCHI✔️❌:                     bCurIsoStereoSp3 = !Eql_INChI_Stereo( pINChI->Stereo, EQL_SP3, pINChI->StereoIsotopic, EQL_SP3,
    // INCHI✔️❌:                         ( pINChI->nFlags & INCHI_FLAG_RAC_STEREO ) || ( pINChI->nFlags & INCHI_FLAG_REL_STEREO ) );
    // INCHI✔️❌:                     if (!bCurIsoStereoSp3)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* Inversion changes iso sp3 differently from non-iso sp3 Fix11 */
    // INCHI✔️❌:                         bCurIsoStereoSp3Inv &= ( pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 bCurRelative = bRequestedRelativeStereo && bCurStereoSp3;
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 bCurRelative = bCurRelative &&
    // INCHI✔️❌:                     ( pINChI->Stereo->nNumberOfStereoCenters > 1 ) &&
    // INCHI✔️❌:                     ( pINChI->Stereo->nCompInv2Abs != 0 ) &&
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                     bCurIsoRelative = bRequestedRelativeStereo && ( bCurIsoStereoSp3 || bCurIsoStereoSp3Inv );
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 bCurIsoRelative = bCurIsoRelative &&
    // INCHI✔️❌:                     ( pINChI->StereoIsotopic->nNumberOfStereoCenters > 1 ) &&
    // INCHI✔️❌:                     ( pINChI->StereoIsotopic->nCompInv2Abs != 0 ) &&
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                     bCurRacemic = bRequestedRacemicStereo && bCurStereoSp3;
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 bCurRacemic = bCurRacemic &&
    // INCHI✔️❌:                     ( pINChI->Stereo->nCompInv2Abs != 0 ) &&
    // INCHI✔️❌:                     ( pINChI->Stereo->nNumberOfStereoCenters > 0 ) ?
    // INCHI✔️❌:                     pINChI->Stereo->nNumberOfStereoCenters : 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:                 bCurIsoRacemic = bRequestedRacemicStereo && ( bCurIsoStereoSp3 || bCurIsoStereoSp3Inv );
    // INCHI✔️❌: #if ( REL_RAC_STEREO_IGN_1_SC == 1 )
    // INCHI✔️❌:                 bCurIsoRacemic = bCurIsoRacemic &
    // INCHI✔️❌:                     ( pINChI->StereoIsotopic->nCompInv2Abs != 0 ) &&
    // INCHI✔️❌:                     ( pINChI->StereoIsotopic->nNumberOfStereoCenters > 0 ) ?
    // INCHI✔️❌:                     pINChI->StereoIsotopic->nNumberOfStereoCenters : 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                 if (bRequestedRelativeStereo)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     bCurStereoSp3 = bCurRelative || (bCurStereoSp3 && ( pINChI->Stereo->nNumberOfStereoCenters > 1 )); /* Fix11 */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     bCurIsoStereoSp3 = bCurIsoRelative ? bCurIsoStereoSp3 : 0;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (bRequestedRacemicStereo)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         bCurStereoSp3 = bCurRacemic > 1 || (bCurStereoSp3 && ( pINChI->Stereo->nNumberOfStereoCenters > 1 )); /* Fix11 */ /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                         bCurIsoStereoSp3 = bCurIsoRacemic > 1 ? bCurIsoStereoSp3 : 0;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:                 io.bIsotopic |= bCurIso;
    // INCHI✔️❌:                 bHasIsotopicAtoms[ii] |= bCurIso;
    // INCHI✔️❌:                 bStereoSp2[ii] |= bCurStereoSp2;
    // INCHI✔️❌:                 bStereoSp3[ii] |= bCurStereoSp3;
    // INCHI✔️❌:                 io.bIgn_UU_Sp3[ii] |= !bCurStereoSp3 && ( pINChI->nFlags & INCHI_FLAG_SC_IGN_ALL_UU );
    // INCHI✔️❌:                 io.bIgn_UU_Sp2[ii] |= !bCurStereoSp2 && ( pINChI->nFlags & INCHI_FLAG_SB_IGN_ALL_UU );
    // INCHI✔️❌:                 bIsotopicStereoSp2[ii] |= bCurIsoStereoSp2;
    // INCHI✔️❌:                 bIsotopicStereoSp3[ii] |= bCurIsoStereoSp3;
    // INCHI✔️❌:                 io.bIgn_UU_Sp3_Iso[ii] |= !bCurIsoStereoSp3 && ( pINChI->nFlags & INCHI_FLAG_SC_IGN_ALL_ISO_UU );
    // INCHI✔️❌:                 io.bIgn_UU_Sp2_Iso[ii] |= !bCurIsoStereoSp2 && ( pINChI->nFlags & INCHI_FLAG_SB_IGN_ALL_ISO_UU );
    // INCHI✔️❌:                 bStereoAbs[ii] |= bCurStereoSp3 && ( pINChI->Stereo->nCompInv2Abs != 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:                 bStereoAbsInverted[ii] |= bCurStereoSp3 && ( pINChI->Stereo->nCompInv2Abs < 0 );
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Fix08: missing isotopic inverted flag if isotopic = inverted non-isotopic */
    // INCHI✔️❌:                 bIsotopicStereoAbsInverted[ii] |= (bCurIsoStereoSp3 && ( pINChI->StereoIsotopic->nCompInv2Abs < 0 )) ||
    // INCHI✔️❌:                     (!bCurIsoStereoSp3  && pINChI->StereoIsotopic  && pINChI->Stereo &&
    // INCHI✔️❌:                     pINChI->StereoIsotopic->nCompInv2Abs &&
    // INCHI✔️❌:                     pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs); /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Fix 11: missing /s1 if only isotopic stereo is inverted */
    // INCHI✔️❌:                 bIsotopicStereoAbs[ii] |= (bCurIsoStereoSp3 && ( pINChI->StereoIsotopic->nCompInv2Abs != 0 )) ||
    // INCHI✔️❌:                     (!bCurIsoStereoSp3  && pINChI->StereoIsotopic  && pINChI->Stereo &&
    // INCHI✔️❌:                     pINChI->StereoIsotopic->nCompInv2Abs &&
    // INCHI✔️❌:                     pINChI->StereoIsotopic->nCompInv2Abs != pINChI->Stereo->nCompInv2Abs); /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:
    // INCHI✔️❌:                 io.bRelativeStereo[ii] |= bCurRelative;
    // INCHI✔️❌:                 io.bIsotopicRelativeStereo[ii] |= bCurIsoRelative;
    // INCHI✔️❌:                 io.bRacemicStereo[ii] |= bCurRacemic;
    // INCHI✔️❌:                 io.bIsotopicRacemicStereo[ii] |= bCurIsoRacemic;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:                 bTautomericAcid |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_ACID_TAUT ) );
    // INCHI✔️❌:                 bHardAddRemProton |= ( 0 != ( pINChI->nFlags & INCHI_FLAG_HARD_ADD_REM_PROTON ) );
    // INCHI✔️❌:                 if (bCurTaut)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     io.bTautomeric |= 1; /* tautomeric representation is present */
    // INCHI✔️❌:                     /* does tautomeric structure have also a non-tautomeric repesentation? */
    // INCHI✔️❌:                     io.bNonTautomeric |= HAS_N( is );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Auxiliary info */
    // INCHI✔️❌:                 if (!( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) && ( pINChI_Aux = is->pINChI_Aux[jj] ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* detect presence of constitutional equivalence onfo */
    // INCHI✔️❌:                     int bCurEqu, bCurTautEqu = 0, bCurIsoEqu = 0, bCurIsoTautEqu = 0; /* Fix15-disabled */
    // INCHI✔️❌:                     io.bAtomEqu[ii] |= ( bCurEqu = bHasEquString( pINChI_Aux->nConstitEquNumbers,
    // INCHI✔️❌:                         pINChI_Aux->nNumberOfAtoms ) ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                     if (bCurTaut)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         io.bTautEqu[ii] |= ( bCurTautEqu = bHasEquString( pINChI_Aux->nConstitEquTGroupNumbers,
    // INCHI✔️❌:                             pINChI_Aux->nNumberOfTGroups ) ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (bCurIso)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         io.bIsotopicAtomEqu[ii] |= ( bCurIsoEqu = bHasEquString( pINChI_Aux->nConstitEquIsotopicNumbers,
    // INCHI✔️❌:                             pINChI_Aux->nNumberOfAtoms ) ) /*|| bCurEqu*/; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                         if (bCurTaut)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             io.bIsotopicTautEqu[ii] |= ( bCurIsoTautEqu = bHasEquString( pINChI_Aux->nConstitEquIsotopicTGroupNumbers,
    // INCHI✔️❌:                                 pINChI_Aux->nNumberOfTGroups ) ) /*|| bCurTautEqu*/; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         /* non-zero if isotopic numbering for inverted isotopic stereo is different */
    // INCHI✔️❌:                         io.bIsotopicOrigNumb[ii] |= bCurHasIsoStereo && /* Fix14 */
    // INCHI✔️❌:                             pINChI_Aux->nOrigAtNosInCanonOrdInv &&
    // INCHI✔️❌:                             pINChI_Aux->nIsotopicOrigAtNosInCanonOrd &&
    // INCHI✔️❌:                             ( 0 != memcmp( pINChI_Aux->nOrigAtNosInCanonOrdInv,
    // INCHI✔️❌:                                 pINChI_Aux->nIsotopicOrigAtNosInCanonOrd,
    // INCHI✔️❌:                                 sizeof( pINChI_Aux->nOrigAtNosInCanonOrdInv[0] ) * pINChI_Aux->nNumberOfAtoms ) );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* Inverted stereo */
    // INCHI✔️❌:                     if (bCurStereoSp3 && pINChI->Stereo->nCompInv2Abs)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         io.bInvStereo[ii] |= 1;
    // INCHI✔️❌:                         io.bInvStereoOrigNumb[ii] |= pINChI_Aux->nOrigAtNosInCanonOrd &&
    // INCHI✔️❌:                             pINChI_Aux->nOrigAtNosInCanonOrdInv &&
    // INCHI✔️❌:                             ( 0 != memcmp( pINChI_Aux->nOrigAtNosInCanonOrd,
    // INCHI✔️❌:                                 pINChI_Aux->nOrigAtNosInCanonOrdInv,
    // INCHI✔️❌:                                 sizeof( pINChI_Aux->nOrigAtNosInCanonOrd[0] ) * pINChI_Aux->nNumberOfAtoms ) );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* Inverted isotopic stereo */
    // INCHI✔️❌:                     if (bCurIsoStereoSp3 && pINChI->StereoIsotopic->nCompInv2Abs)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         io.bInvIsotopicStereo[ii] |= 1;
    // INCHI✔️❌:
    // INCHI✔️❌:                         io.bInvIsotopicStereoOrigNumb[ii]
    // INCHI✔️❌:                             |= pINChI_Aux->nIsotopicOrigAtNosInCanonOrd &&
    // INCHI✔️❌:                             pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv &&
    // INCHI✔️❌:                             ( 0 != memcmp( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd,
    // INCHI✔️❌:                                 pINChI_Aux->nIsotopicOrigAtNosInCanonOrdInv,
    // INCHI✔️❌:                                 sizeof( pINChI_Aux->nIsotopicOrigAtNosInCanonOrd[0] ) * pINChI_Aux->nNumberOfAtoms ) );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (pINChI_Aux->OrigInfo && bHasOrigInfo( pINChI_Aux->OrigInfo, pINChI_Aux->nNumberOfAtoms ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         io.bChargesRadVal[ii] |= 1;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (bCompExists)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (j = TAUT_NON; j < TAUT_NUM; j++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 io.num_comp[j] ++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (io.bTautomeric /*&& bTautomericAcid*/) /* "&& bTautomericAcid" commented out 2004-06-02 */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         io.bTautomeric += bTautomericAcid; /* long-range tautomerism */
    // INCHI✔️❌:         io.bTautomeric += ( bHardAddRemProton ? 4 : 0 );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (bRequestedRacemicStereo || bRequestedRelativeStereo)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* do not output inverted stereo info */
    // INCHI✔️❌:         for (i = 0; i < TAUT_NUM; i++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Fix11 */
    // INCHI✔️❌:             bStereoAbsInverted[i] =
    // INCHI✔️❌:                 bStereoAbs[i] =
    // INCHI✔️❌:                 io.bInvStereo[i] =
    // INCHI✔️❌:                 io.bInvStereoOrigNumb[i] = 0;
    // INCHI✔️❌:                 /* io.bIsotopicRelativeStereo[i]=0 may happen because iso stereo is same or inverted non-iso stereo */
    // INCHI✔️❌:             bIsotopicStereoAbsInverted[i] =
    // INCHI✔️❌:                 bIsotopicStereoAbs[i] =
    // INCHI✔️❌:                 io.bInvIsotopicStereo[i] =
    // INCHI✔️❌:                 io.bInvIsotopicStereoOrigNumb[i] = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     io.iCurTautMode = io.bOutType == OUT_N1 ? TAUT_NON :  /*  only non-taut */
    // INCHI✔️❌:
    // INCHI✔️❌:         io.bOutType == OUT_T1 ? TAUT_YES :      /*  tautomeric if present, otherwise non-tautomeric     */
    // INCHI✔️❌:         io.bOutType == OUT_NT ? TAUT_NON :      /*  only non-taut representations of tautomeric         */
    // INCHI✔️❌:         io.bOutType == OUT_TN ? TAUT_YES :       /*  tautomeric if present otherwise non-tautomeric;     */
    // INCHI✔️❌:         -1; /*  separately output non-taut representations of tautomeric if present */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (io.iCurTautMode < 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0;  /* error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Now print out */
    // INCHI✔️❌:
    // INCHI✔️❌:     io.bOverflow = 0;
    // INCHI✔️❌:     io.num_components = io.num_comp[io.iCurTautMode];
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (bINChIOutputOptions & INCHI_OUT_ONLY_AUX_INFO)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto output_aux_info;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     io.nCurINChISegment = DIFL_M;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: version and kind */
    // INCHI✔️❌:     if (INCHI_basic_or_INCHI_reconnected == INCHI_BAS || !( bINChIOutputOptions & INCHI_OUT_EMBED_REC ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int is_beta = 0;
    // INCHI✔️❌:         int nAtomsAllComp = inchi_max( nAtomsAllComp1, nAtomsAllComp2 );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (nAtomsAllComp > NORMALLY_ALLOWED_INP_MAX_ATOMS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* v. 1.05 for LargeMolecules */
    // INCHI✔️❌:             is_beta = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (pOrigStruct && pOrigStruct->polymer)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* v. 1.05 for Polymers */
    // INCHI✔️❌:             is_beta = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* specifically put 'B' for empty structure InChI    */
    // INCHI✔️❌:         /* if "Polymers" or "LargeMolecules" requested        */
    // INCHI✔️❌:         else if (!pOrigStruct && ( ip->bLargeMolecules || ip->bPolymers ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             is_beta = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (pOrigStruct && pOrigStruct->n_zy)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             is_beta = 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         OutputINCHI_VersionAndKind( out_file, strbuf, bINChIOutputOptions, is_beta, pLF, pTAB );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: atoms */
    // INCHI✔️❌:     intermediate_result = OutputINCHI_MainLayerFormula( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                         num_components2,
    // INCHI✔️❌:                                                         &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                         &io, pLF, pTAB );
    // INCHI✔️❌:     if (intermediate_result != 0)
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: connection table */
    // INCHI✔️❌:     intermediate_result = OutputINCHI_MainLayerConnections( pCG, out_file, strbuf, num_components2,
    // INCHI✔️❌:                                                                 &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                                 &io, pLF, pTAB );
    // INCHI✔️❌:     if (intermediate_result != 0)
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: hydrogens (with tautomeric info) */
    // INCHI✔️❌:     intermediate_result = OutputINCHI_MainLayerHydrogens( pCG, out_file, strbuf, num_components2,
    // INCHI✔️❌:                                                               &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                               &io, pLF, pTAB );
    // INCHI✔️❌:     if (intermediate_result != 0)
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:     io.bFhTag = 0;
    // INCHI✔️❌:     npass = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: repeat_INChI_output:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: charge and  removed protons */
    // INCHI✔️❌:     intermediate_result = OutputINCHI_ChargeAndRemovedAddedProtonsLayers( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                                               &io, pLF, pTAB );
    // INCHI✔️❌:     if (intermediate_result != 0)
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: polymer layer */
    // INCHI✔️❌:     if (npass == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         intermediate_result = OutputINCHI_PolymerLayer( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                         &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                         orig_inp_data, pOrigStruct,
    // INCHI✔️❌:                                                         &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: stereo (non-isotopic) */
    // INCHI✔️❌:     intermediate_result = OutputINCHI_StereoLayer( pCG, out_file, strbuf, &io, pLF, pTAB );
    // INCHI✔️❌:     if (intermediate_result != 0)
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Switch from M to MI or from F to FI */
    // INCHI✔️❌:     io.nCurINChISegment++;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: isotopic */
    // INCHI✔️❌:     intermediate_result = OutputINCHI_IsotopicLayer( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                          &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                          &io, pLF, pTAB );
    // INCHI✔️❌:
    // INCHI✔️❌:     if (intermediate_result != 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         At this point the INChI part of the output has been done.
    // INCHI✔️❌:         If this INChI is tautomeric and non-tautomeric results exist,
    // INCHI✔️❌:         then we need to output non-tautomeric data:
    // INCHI✔️❌:             fixed H
    // INCHI✔️❌:             stereo
    // INCHI✔️❌:             isotopic
    // INCHI✔️❌:             isotopic stereo
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: FixedH and sublayers */
    // INCHI✔️❌:     intermediate_result = OutputINCHI_FixedHLayerWithSublayers( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                                 &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                                 &io, pLF, pTAB,
    // INCHI✔️❌:                                                                 &then_goto_repeat );
    // INCHI✔️❌:     if (intermediate_result != 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (then_goto_repeat)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         npass++;
    // INCHI✔️❌:         goto repeat_INChI_output;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         InChI output:  reconnected structure
    // INCHI✔️❌:     */
    // INCHI✔️❌:
    // INCHI✔️❌:     bEmbeddedOutputCalled = 0;
    // INCHI✔️❌:     if (bDisconnectedCoord && INCHI_basic_or_INCHI_reconnected == INCHI_BAS &&
    // INCHI✔️❌:         ( bINChIOutputOptions & INCHI_OUT_EMBED_REC ) && num_components2[INCHI_REC])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int nRet;
    // INCHI✔️❌:         bEmbeddedOutputCalled = 1;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* output blank line before /R: in case of bPlainTextCommnts=1 */
    // INCHI✔️❌:         inchi_ios_print_nodisplay( out_file, "%s", pLF );
    // INCHI✔️❌:         /* end of disconnected INChI output */
    // INCHI✔️❌:
    // INCHI✔️❌:         nRet = OutputINChI1( pCG,
    // INCHI✔️❌:                              strbuf,
    // INCHI✔️❌:                              pINChISortTautAndNonTaut2,
    // INCHI✔️❌:                              INCHI_REC,
    // INCHI✔️❌:                              orig_inp_data,
    // INCHI✔️❌:                              pOrigStruct,
    // INCHI✔️❌:                              ip,
    // INCHI✔️❌:                              0 /*bDisconnectedCoord*/,
    // INCHI✔️❌:                              bOutputType,
    // INCHI✔️❌:                              bINChIOutputOptions | INCHI_OUT_NO_AUX_INFO,
    // INCHI✔️❌:                              num_components2,
    // INCHI✔️❌:                              num_non_taut2,
    // INCHI✔️❌:                              num_taut2,
    // INCHI✔️❌:                              out_file,
    // INCHI✔️❌:                              log_file,
    // INCHI✔️❌:                              num_input_struct,
    // INCHI✔️❌:                              pSortPrintINChIFlags,
    // INCHI✔️❌:                              save_opt_bits );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!nRet)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function; /* error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* InChI output: save InChI creation options if requested */
    // INCHI✔️❌:     if (!bEmbeddedOutputCalled &&
    // INCHI✔️❌:         ( bINChIOutputOptions & INCHI_OUT_SAVEOPT ) &&
    // INCHI✔️❌:         ( 0 == ( bINChIOutputOptions & INCHI_OUT_STDINCHI ) )    /* not std-InChI output */
    // INCHI✔️❌:         )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         char let1, let2;
    // INCHI✔️❌:         GetSaveOptLetters( save_opt_bits, &let1, &let2 );
    // INCHI✔️❌:         inchi_ios_print_nodisplay( out_file, "\\%c%c", let1, let2 );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (!bEmbeddedOutputCalled && !bPlainTextCommnts)
    // INCHI✔️❌:     { /* plain text comment earlier ended with LF */
    // INCHI✔️❌:         inchi_ios_print_nodisplay( out_file, "%s%s",
    // INCHI✔️❌:             ( !num_components2[0] && !num_components2[1] ) ? "//" : "", /* empty InChI=// */
    // INCHI✔️❌:             ( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) ? "\n" : pTAB );
    // INCHI✔️❌: /* end of INChI= output */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:     /* @@@ Here we end up with silent output: display previously hidden output */
    // INCHI✔️❌:     if (inchi_ios_flush_not_displayed( out_file ) != -1)
    // INCHI✔️❌:         silent = 0;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: output_aux_info:
    // INCHI✔️❌:
    // INCHI✔️❌:     /*  Output Aux Info */
    // INCHI✔️❌:
    // INCHI✔️❌:     io.bFhTag = 0;
    // INCHI✔️❌:     if (( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ) == 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         io.num_components = io.num_comp[io.iCurTautMode];
    // INCHI✔️❌:
    // INCHI✔️❌:         /* AuxInfo: header and normalization type */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_HeaderAndNormalization_type( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                                          bINChIOutputOptions,
    // INCHI✔️❌:                                                                          &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                                          num_components2,
    // INCHI✔️❌:                                                                          &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:     repeat_INChI_Aux_output:
    // INCHI✔️❌:
    // INCHI✔️❌:         /* AuxInfo: original atom numbers and symmetry numbers (constit. equivalence /E: )    */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_OriginalNumbersAndEquivalenceClasses( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                                                   num_components2,
    // INCHI✔️❌:                                                                                   &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* AuxInfo: tautomeric groups equivalence */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_TautomericGroupsEquivalence( pCG, out_file, strbuf, &io );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* AuxInfo: stereo data */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_Stereo( pCG, out_file, strbuf, &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:
    // INCHI✔️❌:     repeat_INChI_Aux_Iso_output:
    // INCHI✔️❌:             /* AuxInfo: isotopic info */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_IsotopicInfo( pCG, out_file, strbuf,
    // INCHI✔️❌:                                                           &INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                           &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /*
    // INCHI✔️❌:           At this point the INChI_Aux part of the output has been completed.
    // INCHI✔️❌:           If this INChI is tautomeric and non-tautomeric results exist,
    // INCHI✔️❌:           then we need to output non-tautomeric auxilialy data
    // INCHI✔️❌:           (same as above excluding tautomeric information).
    // INCHI✔️❌:           Currently, this is enabled for xml output only
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         if (io.bOutType == OUT_TN && io.bTautomeric && io.bNonTautomeric &&
    // INCHI✔️❌:             /* Check whether the Fixed-H layer is empty */
    // INCHI✔️❌:             ( *pSortPrintINChIFlags & ( ( INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
    // INCHI✔️❌:                 FLAG_SORT_PRINT_NO_NFIX_H_REC ) ) &&
    // INCHI✔️❌:                 ( *pSortPrintINChIFlags & ( ( INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
    // INCHI✔️❌:                     FLAG_SORT_PRINT_NO_IFIX_H_REC ) )
    // INCHI✔️❌:               )
    // INCHI✔️❌:         {
    // INCHI✔️❌:             io.bNonTautomeric = 0; /* bNonTautIdentifierNotEmpty == 0 => no fixed H info 02-10-2995 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (io.bOutType == OUT_TN && io.bTautomeric && io.bNonTautomeric)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* add the second (non-tautomeric) output */
    // INCHI✔️❌:             io.bOutType = OUT_NONTAUT;
    // INCHI✔️❌:             io.iCurTautMode = TAUT_NON;
    // INCHI✔️❌:             io.pINChISort = io.pINChISortTautAndNonTaut[TAUT_NON];
    // INCHI✔️❌:             io.bSecondNonTautPass = 1;
    // INCHI✔️❌:             io.num_components = io.num_comp[io.iCurTautMode];
    // INCHI✔️❌:             io.bFhTag = AL_FIXH;
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf ); /*pStr[io.tot_len=0] = '\0';*/
    // INCHI✔️❌:
    // INCHI✔️❌:             /* if InChI Fixed-H isotopic is empty then do not output corresponding AuxInfo */
    // INCHI✔️❌:             if (!( *pSortPrintINChIFlags &
    // INCHI✔️❌:                 ( ( INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
    // INCHI✔️❌:                     FLAG_SORT_PRINT_NO_NFIX_H_REC ) )
    // INCHI✔️❌:                )
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 npass++;
    // INCHI✔️❌:                 goto repeat_INChI_Aux_output;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 npass++;
    // INCHI✔️❌:                 goto repeat_INChI_Aux_Iso_output;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io.bOutType == OUT_NONTAUT && io.bOutputType == OUT_TN && io.bTautomeric && io.bNonTautomeric)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* the second (non-taut) output has been done; restore variables */
    // INCHI✔️❌:                 io.bOutType = OUT_TN;
    // INCHI✔️❌:                 io.iCurTautMode = TAUT_YES;
    // INCHI✔️❌:                 io.pINChISort = io.pINChISortTautAndNonTaut[TAUT_YES];
    // INCHI✔️❌:                 io.bSecondNonTautPass = 0;
    // INCHI✔️❌:                 /* set correct num components for the reversibility info 02-10-2005 */
    // INCHI✔️❌:                 io.num_components = io.num_comp[io.iCurTautMode];
    // INCHI✔️❌:                 io.bFhTag = 0;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*    Charges, radicals, unusual valences */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_ChargesRadicalsAndUnusualValences( pCG, out_file, strbuf, &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Output the original input structure -- quick fix */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_ReversibilityInfo( pCG, out_file, strbuf, pOrigStruct, &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Output polymeric Aux Info */
    // INCHI✔️❌:         intermediate_result = OutputAUXINFO_PolymerInfo( pCG, out_file, strbuf, pOrigStruct, &io, pLF, pTAB );
    // INCHI✔️❌:         if (intermediate_result != 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*
    // INCHI✔️❌:             Output INChI_Aux of the reconnected structure
    // INCHI✔️❌:         */
    // INCHI✔️❌:
    // INCHI✔️❌:         bEmbeddedOutputCalled = 0;
    // INCHI✔️❌:         if (bDisconnectedCoord && INCHI_basic_or_INCHI_reconnected == INCHI_BAS && ( bINChIOutputOptions & INCHI_OUT_EMBED_REC ) &&
    // INCHI✔️❌:              num_components2[INCHI_REC] && !( bINChIOutputOptions & INCHI_OUT_NO_AUX_INFO ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int nRet;
    // INCHI✔️❌:             bEmbeddedOutputCalled = 1;
    // INCHI✔️❌:             inchi_ios_print( out_file, "%s", pLF );
    // INCHI✔️❌:
    // INCHI✔️❌:             nRet = OutputINChI1( pCG,
    // INCHI✔️❌:                                  strbuf,
    // INCHI✔️❌:                                  pINChISortTautAndNonTaut2,
    // INCHI✔️❌:                                  INCHI_REC,
    // INCHI✔️❌:                                  NULL,
    // INCHI✔️❌:                                  NULL,
    // INCHI✔️❌:                                  ip,
    // INCHI✔️❌:                                  0 /*bDisconnectedCoord*/,
    // INCHI✔️❌:                                  bOutputType,
    // INCHI✔️❌:                                  INCHI_OUT_ONLY_AUX_INFO | bINChIOutputOptions,
    // INCHI✔️❌:                                  num_components2,
    // INCHI✔️❌:                                  num_non_taut2,
    // INCHI✔️❌:                                  num_taut2,
    // INCHI✔️❌:                                  out_file,
    // INCHI✔️❌:                                  log_file,
    // INCHI✔️❌:                                  num_input_struct,
    // INCHI✔️❌:                                  pSortPrintINChIFlags,
    // INCHI✔️❌:                                  save_opt_bits );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!nRet)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 goto exit_function; /* error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Close INChI_Aux */
    // INCHI✔️❌:         if (!bEmbeddedOutputCalled && !bPlainTextCommnts)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print( out_file, "%s\n", ( !num_components2[0] && !num_components2[1] ) ? "//" : "" );
    // INCHI✔️❌:             /* plain text comment earlier ended with LF */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* in wINChI window, separate AuxInfo: from InChIKey: with blank line */
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s",
    // INCHI✔️❌:             ( bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) ? "\n" : "" );
    // INCHI✔️❌:     } /* end of output AuxInfo */
    // INCHI✔️❌:
    // INCHI✔️❌:     ret = 1;
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: #ifdef TARGET_LIB_FOR_WINCHI
    // INCHI✔️❌:     /* @@@ If for any error we get here silent, display previously hidden output */
    // INCHI✔️❌:     if (silent)
    // INCHI✔️❌:     {
    // INCHI✔️❌:      /*
    // INCHI✔️❌:         if ( !inchi_ios_flush_not_displayed( out_file ) != -1  )
    // INCHI✔️❌:             silent = 0;
    // INCHI✔️❌:     */
    // INCHI✔️❌:         silent = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     if (io.bOverflow)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_ios_print( out_file, "\nFATAL ERROR: Output buffer overflow\n" );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (intermediate_result)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ret = 0;
    // INCHI✔️❌:         inchi_ios_eprint( log_file, "InChI serialization error for structure #%d.%s%s%s%s\n",
    // INCHI✔️❌:                                     num_input_struct, SDF_LBL_VAL( ip->pSdfLabel, ip->pSdfValue ) );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: } /* OutputINChI1 */
    // INCHI✔️❌:
    // INCHI✔️❌:
    // INCHI✔️❌: /****************************************************************************/
    // INCHI✔️❌: char *szGetTag( const INCHI_TAG *Tag,
    // INCHI✔️❌:                 int             nTag,
    // INCHI✔️❌:                 int             bTag,
    // INCHI✔️❌:                 char            *szTag,
    // INCHI✔️❌:                 int             *bAlways,
    // INCHI✔️❌:                 short           tag_flag)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, j, bit, num, len;
    // INCHI✔️❌:     const int MAX_TAG_NUM = tag_flag ? (int)IL_MAX_ORD : (int)AL_MAX_ORD; /* djb-rwth: fixing GHI #160 */
    // INCHI✔️❌:     if (0 < nTag && nTag < 3)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* no plain text comments: pick up the last tag */
    // INCHI✔️❌:         for (i = 0, j = -1, bit = 1; i < MAX_TAG_NUM; i++, bit <<= 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (bTag & bit)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 j = i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (j >= 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌: #if USE_BCF
    // INCHI✔️❌:             int stl1, stl2, dstsz;
    // INCHI✔️❌:             stl1 = strlen(Tag[j].szXmlLabel) + 1;
    // INCHI✔️❌:             stl2 = strlen(Tag[j].szPlainLabel) + 1;
    // INCHI✔️❌:             dstsz = max_3(stl1, stl2, 5);
    // INCHI✔️❌:             strcpy_s( szTag, dstsz, nTag == 1 ? Tag[j].szXmlLabel : nTag == 2 ? Tag[j].szPlainLabel : "???" ); /* djb-rwth: function replaced with its safe C11 variant */
    // INCHI✔️❌: #else
    // INCHI✔️❌:             strcpy(szTag, nTag == 1 ? Tag[j].szXmlLabel : nTag == 2 ? Tag[j].szPlainLabel : "???"); /* djb-rwth: addressing coverity ID #499488 -- when nTag == 2, the "???" is avoided, which is correct */
    // INCHI✔️❌: #endif
    // INCHI✔️❌:             if (nTag != 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 *bAlways = Tag[j].bAlwaysOutput;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return szTag;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:         if (nTag == 3)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* plain text with comments */
    // INCHI✔️❌:             szTag[0] = '{';
    // INCHI✔️❌:             szTag[1] = '\0';
    // INCHI✔️❌:             for (i = 0, j = -1, bit = 1, num = 0; i < MAX_TAG_NUM; i++, bit <<= 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (bTag & bit)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     j = i;
    // INCHI✔️❌:                     if (num++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         strcat(szTag, ":");
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     strcat(szTag, Tag[i].szPlainComment);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (num)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcat(szTag, "}");
    // INCHI✔️❌:                 num = (int) strlen( Tag[j].szPlainLabel );
    // INCHI✔️❌:                 len = (int) strlen( szTag );
    // INCHI✔️❌:                 if (len)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     memmove(szTag + num, szTag, (long long)len + 1); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                     memcpy(szTag, Tag[j].szPlainLabel, num);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     strcpy(szTag, Tag[j].szPlainLabel);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 *bAlways = Tag[j].bAlwaysOutput;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 strcpy(szTag, "???");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return szTag;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:     strcpy(szTag, "???");
    // INCHI✔️❌:     return szTag;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINChI1

    fn inchi(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI>,
    ) -> Result<Option<crate::source_types::INChI>, SourceHeapError> {
        if pointer.is_null() {
            return Ok(None);
        }
        Ok(Some(
            heap.slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        ))
    }
    fn aux(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI_Aux>,
    ) -> Result<Option<crate::source_types::INChI_Aux>, SourceHeapError> {
        if pointer.is_null() {
            return Ok(None);
        }
        Ok(Some(
            heap.slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        ))
    }
    fn stereo(
        heap: &SourceHeap,
        pointer: SourceMutPointer<crate::source_types::INChI_Stereo>,
    ) -> Result<Option<crate::source_types::INChI_Stereo>, SourceHeapError> {
        if pointer.is_null() {
            return Ok(None);
        }
        Ok(Some(
            heap.slice(pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone(),
        ))
    }
    fn selected_representation(
        heap: &SourceHeap,
        sort: &crate::source_types::INCHI_SORT,
        mode: i32,
    ) -> Result<i32, SourceHeapError> {
        let mobile = inchi(heap, sort.pINChI[crate::source_types::TAUT_YES as usize])?;
        let fixed = inchi(heap, sort.pINChI[crate::source_types::TAUT_NON as usize])?;
        let mobile_exists = mobile
            .as_ref()
            .is_some_and(|value| value.nNumberOfAtoms != 0);
        let fixed_exists = fixed
            .as_ref()
            .is_some_and(|value| value.nNumberOfAtoms != 0);
        let mobile_non_taut =
            mobile_exists && mobile.as_ref().is_some_and(|value| value.lenTautomer == 0);
        let fixed_non_taut =
            fixed_exists && fixed.as_ref().is_some_and(|value| value.lenTautomer == 0);
        Ok(match mode {
            value if value == crate::source_types::OUT_N1 as i32 => {
                if mobile_non_taut {
                    crate::source_types::TAUT_YES as i32
                } else if fixed_non_taut {
                    crate::source_types::TAUT_NON as i32
                } else {
                    -1
                }
            }
            value
                if value == crate::source_types::OUT_T1 as i32
                    || value == crate::source_types::OUT_TN as i32 =>
            {
                if mobile_exists {
                    crate::source_types::TAUT_YES as i32
                } else if fixed_exists {
                    crate::source_types::TAUT_NON as i32
                } else {
                    -1
                }
            }
            value if value == crate::source_types::OUT_NN as i32 => {
                if fixed_non_taut {
                    crate::source_types::TAUT_NON as i32
                } else if mobile_non_taut {
                    crate::source_types::TAUT_YES as i32
                } else {
                    -1
                }
            }
            value if value == crate::source_types::OUT_NT as i32 => {
                if mobile
                    .as_ref()
                    .is_some_and(|value| value.nNumberOfAtoms != 0 && value.lenTautomer > 0)
                    && fixed_non_taut
                {
                    crate::source_types::TAUT_NON as i32
                } else {
                    -1
                }
            }
            _ => -1,
        })
    }
    fn has_fixed(
        heap: &SourceHeap,
        sort: &crate::source_types::INCHI_SORT,
    ) -> Result<i32, SourceHeapError> {
        Ok(i32::from(
            inchi(heap, sort.pINChI[crate::source_types::TAUT_NON as usize])?
                .is_some_and(|value| value.nNumberOfAtoms != 0),
        ))
    }
    fn pointers_differ(
        heap: &SourceHeap,
        first: SourceMutPointer<u16>,
        second: SourceMutPointer<u16>,
        count: i32,
    ) -> Result<i32, SourceHeapError> {
        if first.is_null() || second.is_null() {
            return Ok(0);
        }
        let count = usize::try_from(count).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        let first = heap
            .slice(first.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let second = heap
            .slice(second.as_const())?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        Ok(i32::from(first != second))
    }
    fn print_values(
        heap: &mut SourceHeap,
        output: &mut INCHI_IOSTREAM,
        stdout: SourceMutPointer<FILE>,
        display: bool,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect::<Vec<_>>())?;
        let arguments = SourceVaList {
            arguments,
            ..SourceVaList::default()
        };
        if display {
            inchi_ios_print(heap, Some(output), stdout, format.as_const(), &arguments)?;
        } else {
            inchi_ios_print_nodisplay(heap, Some(output), stdout, format.as_const(), &arguments)?;
        }
        Ok(())
    }
    fn c_nonempty(
        heap: &SourceHeap,
        pointer: SourceMutPointer<i8>,
    ) -> Result<bool, SourceHeapError> {
        if pointer.is_null() {
            return Ok(false);
        }
        Ok(*heap
            .slice(pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            != 0)
    }
    fn report_serialization_error(
        heap: &mut SourceHeap,
        log: &mut INCHI_IOSTREAM,
        input_parameters: &crate::source_types::INPUT_PARMS,
        input_structure_number: i32,
        allocation_message: bool,
    ) -> Result<(), SourceHeapError> {
        let empty = heap.allocate_model_storage(vec![0_i8])?;
        let space = heap.allocate_model_storage(vec![b' ' as i8, 0])?;
        let equal = heap.allocate_model_storage(vec![b'=' as i8, 0])?;
        let missing = heap.allocate_model_storage(
            b"Missing SDfile data value\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        )?;
        let has_label = c_nonempty(heap, input_parameters.pSdfLabel)?;
        let has_value = c_nonempty(heap, input_parameters.pSdfValue)?;
        let first = if has_label { space } else { empty };
        let label = if has_label {
            input_parameters.pSdfLabel
        } else {
            empty
        };
        let separator = if has_label {
            if has_value { equal } else { space }
        } else {
            empty
        };
        let value = if has_value {
            input_parameters.pSdfValue
        } else if has_label {
            missing
        } else {
            empty
        };
        let format = if allocation_message {
            b"Cannot allocate output buffer. No output for structure #%d.%s%s%s%s\n\0".as_slice()
        } else {
            b"InChI serialization error for structure #%d.%s%s%s%s\n\0".as_slice()
        };
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect::<Vec<_>>())?;
        crate::source::base::ichi_io::inchi_ios_eprint(
            heap,
            Some(log),
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Signed(i64::from(input_structure_number)),
                    SourceFormatArgument::Bytes(first.as_const()),
                    SourceFormatArgument::Bytes(label.as_const()),
                    SourceFormatArgument::Bytes(separator.as_const()),
                    SourceFormatArgument::Bytes(value.as_const()),
                ],
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    let Some(string_buffer) = string_buffer else {
        report_serialization_error(heap, log, input_parameters, input_structure_number, true)?;
        return Ok(0);
    };
    if string_buffer.pStr.is_null() || string_buffer.nAllocatedLength <= 0 {
        report_serialization_error(heap, log, input_parameters, input_structure_number, true)?;
        return Ok(0);
    }

    let basic_index =
        usize::try_from(basic_or_reconnected).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let number_of_components = *component_counts
        .get(basic_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let row_offset = i64::from(basic_or_reconnected)
        .checked_mul(i64::from(crate::source_types::TAUT_NUM))
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    let selected_row = sort_rows.offset(row_offset)?;
    let selected_sorts = heap.slice(selected_row.as_const())?;
    let mobile_sorts = *selected_sorts
        .get(crate::source_types::TAUT_YES as usize)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let fixed_sorts = *selected_sorts
        .get(crate::source_types::TAUT_NON as usize)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let plain_text = i32::from(
        output_options
            & (crate::source_types::INCHI_OUT_PLAIN_TEXT
                | crate::source_types::INCHI_OUT_PLAIN_TEXT_COMMENTS) as i32
            != 0,
    );
    let plain_comments =
        i32::from(output_options & crate::source_types::INCHI_OUT_PLAIN_TEXT_COMMENTS as i32 != 0);
    let mut io = INCHI_OUT_CTL {
        bAbcNumbers: input_parameters.bAbcNumbers,
        ATOM_MODE: (if input_parameters.bAbcNumbers != 0 {
            2
        } else {
            0
        }) | 4
            | 1
            | 16
            | if input_parameters.bAbcNumbers != 0 && input_parameters.bCtPredecessors != 0 {
                32
            } else {
                0
            }
            | if input_parameters.bCtPredecessors != 0 {
                8
            } else {
                0
            },
        TAUT_MODE: if input_parameters.bAbcNumbers != 0 {
            2
        } else {
            0
        },
        pSortPrintINChIFlags: sort_print_flags,
        num_components: number_of_components,
        pINChISortTautAndNonTaut: selected_row,
        pINChISort: mobile_sorts,
        pINChISort2: mobile_sorts,
        bUseMulipliers: 1,
        bOmitRepetitions: 1,
        bPlainTextTags: 2,
        bOutputType: output_type,
        bOutType: output_type,
        bNonTautIsIdenticalToTaut: 1,
        nTag: if plain_comments != 0 {
            3
        } else if plain_text != 0 {
            2
        } else {
            0
        },
        bPolymers: input_parameters.bPolymers,
        ..INCHI_OUT_CTL::default()
    };
    if original_atom_data.is_null() {
        io.n_zy = 0;
        io.n_pzz = 0;
    } else {
        let data = heap
            .slice(original_atom_data)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        io.n_zy = data.n_zy;
        if !data.polymer.is_null() {
            io.n_pzz = heap
                .slice(data.polymer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .n_pzz;
        }
    }

    let mut line_feed: &'static [i8] = &[0];
    let mut tab: &'static [i8] = &[b'\n' as i8, 0];
    set_line_separators(output_options, &mut line_feed, &mut tab);
    let line_feed = heap.allocate_model_storage(line_feed.to_vec())?;
    let tab = heap.allocate_model_storage(tab.to_vec())?;
    io.sDifSegs = [[crate::source_types::tagMarkDiff_DIFV_BOTH_EMPTY as i8; 11]; 4];

    let mut has_isotopic_atoms = [0_i32; 2];
    let mut stereo_sp2 = [0_i32; 2];
    let mut stereo_sp3 = [0_i32; 2];
    let mut isotopic_stereo_sp2 = [0_i32; 2];
    let mut isotopic_stereo_sp3 = [0_i32; 2];
    let mut stereo_abs_inverted = [0_i32; 2];
    let mut isotopic_stereo_abs_inverted = [0_i32; 2];
    let mut stereo_abs = [0_i32; 2];
    let mut isotopic_stereo_abs = [0_i32; 2];
    let mut tautomeric_acid = 0_i32;
    let mut hard_add_remove_proton = 0_i32;
    let mut requested_racemic = 0_i32;
    let mut requested_relative = 0_i32;
    let mut taut_iso_h_number = 0_i32;
    let mut taut_iso_atom = 0_i32;
    let mut has_iso_h = 0_i32;
    let mut taut_and_non_taut = 0_i32;
    let mut taut_is_non_taut = 0_i32;
    let fix_transposed_charge = i32::from(
        output_options & crate::source_types::INCHI_OUT_FIX_TRANSPOSITION_CHARGE_BUG as i32 != 0,
    );
    io.bTautomericOutputAllowed = i32::from(
        io.bOutType == crate::source_types::OUT_T1 as i32
            || io.bOutType == crate::source_types::OUT_TN as i32,
    );
    io.pINChISort = if io.bTautomericOutputAllowed != 0 {
        mobile_sorts
    } else {
        fixed_sorts
    };

    for component in 0..number_of_components {
        let mobile_sort = heap
            .slice(io.pINChISort.offset(i64::from(component))?.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let fixed_sort = heap
            .slice(fixed_sorts.offset(i64::from(component))?.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        crate::source::base::ichimake::CompINChILayers(
            heap,
            &mobile_sort,
            &fixed_sort,
            &mut io.sDifSegs,
            fix_transposed_charge,
        )?;
        io.bNonTautIsIdenticalToTaut &= i32::from(
            crate::source::base::ichimake::CompINChITautVsNonTaut(
                heap,
                &mobile_sort,
                &fixed_sort,
                1,
            )? == 0,
        );
        if let Some(value) = aux(
            heap,
            mobile_sort.pINChI_Aux[crate::source_types::TAUT_YES as usize],
        )? {
            for isotope in 0..3 {
                has_iso_h = has_iso_h.wrapping_add(i32::from(
                    value.nNumRemovedIsotopicH[isotope]
                        .checked_abs()
                        .ok_or(SourceHeapError::SourceIntegerOverflow)?,
                ));
                io.num_iso_H[isotope] = io.num_iso_H[isotope]
                    .wrapping_add(i32::from(value.nNumRemovedIsotopicH[isotope]));
            }
            io.nNumRemovedProtons = io
                .nNumRemovedProtons
                .wrapping_add(i32::from(value.nNumRemovedProtons));
        }
        if io.bTautomericOutputAllowed != 0 {
            let j = crate::source_types::TAUT_YES as i32;
            let jj = if io.bOutType == crate::source_types::OUT_TN as i32 {
                selected_representation(heap, &mobile_sort, crate::source_types::OUT_T1 as i32)?
            } else {
                selected_representation(heap, &mobile_sort, io.bOutType)?
            };
            if io.bOutType == crate::source_types::OUT_TN as i32
                && jj == crate::source_types::TAUT_YES as i32
            {
                let value = inchi(heap, mobile_sort.pINChI[jj as usize])?
                    .ok_or(SourceHeapError::NullPointer)?;
                if value.lenTautomer > 0 {
                    taut_and_non_taut = taut_and_non_taut.wrapping_add(i32::from(
                        value.bDeleted == 0 && has_fixed(heap, &mobile_sort)? != 0,
                    ));
                } else {
                    taut_is_non_taut = taut_is_non_taut.wrapping_add(1);
                }
            }
            if jj == j {
                if let Some(value) = inchi(heap, mobile_sort.pINChI[jj as usize])? {
                    if value.nNumberOfAtoms > 0 {
                        if let Some(value_aux) = aux(heap, mobile_sort.pINChI_Aux[jj as usize])? {
                            taut_iso_h_number = taut_iso_h_number.wrapping_add(
                                i32::from(value_aux.nNumRemovedIsotopicH[0])
                                    .wrapping_add(i32::from(value_aux.nNumRemovedIsotopicH[1]))
                                    .wrapping_add(i32::from(value_aux.nNumRemovedIsotopicH[2])),
                            );
                            taut_iso_atom = taut_iso_atom.wrapping_add(i32::from(
                                value.nNumberOfIsotopicAtoms > 0
                                    || value.nNumberOfIsotopicTGroups > 0,
                            ));
                        }
                    }
                }
            }
        }
    }
    io.sDifSegs[crate::source_types::tagDiffINChILayers_DIFL_M as usize]
        [crate::source_types::tagDiffINChISegments_DIFS_p_PROTONS as usize] =
        if io.nNumRemovedProtons != 0 {
            crate::source_types::tagMarkDiff_DIFV_NEQ2PRECED as i8
        } else {
            crate::source_types::tagMarkDiff_DIFV_BOTH_EMPTY as i8
        };
    io.sDifSegs[crate::source_types::tagDiffINChILayers_DIFL_MI as usize]
        [crate::source_types::tagDiffINChISegments_DIFS_h_H_ATOMS as usize] = if has_iso_h != 0 {
        crate::source_types::tagMarkDiff_DIFV_NEQ2PRECED as i8
    } else {
        crate::source_types::tagMarkDiff_DIFV_BOTH_EMPTY as i8
    };
    crate::source::base::ichimake::MarkUnusedAndEmptyLayers(&mut io.sDifSegs);
    io.bNonTautIsIdenticalToTaut &= i32::from(taut_iso_h_number == 0);

    let mut all_fixed_atoms = 0_i32;
    let mut all_mobile_atoms = 0_i32;
    for component in 0..number_of_components {
        let sort = heap
            .slice(io.pINChISort.offset(i64::from(component))?.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        let mut component_exists = 0_i32;
        for requested in 0..2_i32 {
            let mut found = if io.bOutType == crate::source_types::OUT_TN as i32 {
                if requested == crate::source_types::TAUT_YES as i32 {
                    selected_representation(heap, &sort, crate::source_types::OUT_T1 as i32)?
                } else {
                    selected_representation(heap, &sort, crate::source_types::OUT_NT as i32)?
                }
            } else {
                selected_representation(heap, &sort, io.bOutType)?
            };
            let target = match io.bOutType {
                value if value == crate::source_types::OUT_N1 as i32 => {
                    crate::source_types::TAUT_NON as i32
                }
                value if value == crate::source_types::OUT_T1 as i32 => {
                    crate::source_types::TAUT_YES as i32
                }
                value if value == crate::source_types::OUT_NT as i32 => {
                    crate::source_types::TAUT_NON as i32
                }
                value if value == crate::source_types::OUT_TN as i32 => requested,
                _ => continue,
            };
            if found != requested && io.bOutType != crate::source_types::OUT_TN as i32 {
                continue;
            }
            if io.bOutType == crate::source_types::OUT_TN as i32 && found < 0 {
                if taut_and_non_taut != 0
                    && taut_is_non_taut != 0
                    && requested == crate::source_types::TAUT_NON as i32
                {
                    found =
                        selected_representation(heap, &sort, crate::source_types::OUT_T1 as i32)?;
                    if found < 0 {
                        continue;
                    }
                    let fallback = inchi(heap, sort.pINChI[found as usize])?
                        .ok_or(SourceHeapError::NullPointer)?;
                    if fallback.bDeleted != 0 || fallback.lenTautomer != 0 {
                        continue;
                    }
                } else {
                    continue;
                }
            }
            if found < 0 || found > 1 {
                continue;
            }
            let Some(value) = inchi(heap, sort.pINChI[found as usize])? else {
                continue;
            };
            if value.nNumberOfAtoms <= 0 {
                continue;
            }
            component_exists = component_exists.wrapping_add(1);
            if requested == crate::source_types::TAUT_NON as i32 {
                all_fixed_atoms = all_fixed_atoms.wrapping_add(value.nNumberOfAtoms);
            } else {
                all_mobile_atoms = all_mobile_atoms.wrapping_add(value.nNumberOfAtoms);
            }
            let current_taut = i32::from(value.lenTautomer > 0);
            let current_iso =
                i32::from(value.nNumberOfIsotopicAtoms > 0 || value.nNumberOfIsotopicTGroups > 0);
            let possible_iso_h = i32::from(
                (!value.nPossibleLocationsOfIsotopicH.is_null()
                    && heap.slice(value.nPossibleLocationsOfIsotopicH.as_const())?[0] > 1)
                    || value.lenTautomer > 1,
            );
            let iso_h_stereo = i32::from(
                (possible_iso_h != 0 && (taut_iso_h_number != 0 || taut_iso_atom != 0))
                    || current_iso != 0,
            );
            if found == requested && value.bDeleted != 0 {
                io.num_comp[requested as usize] = io.num_comp[requested as usize].wrapping_sub(1);
                if current_taut != 0 {
                    io.bTautomeric |= 1;
                    io.bNonTautomeric |= has_fixed(heap, &sort)?;
                }
                io.bIsotopic |= current_iso;
                continue;
            }
            let normal_stereo = stereo(heap, value.Stereo)?;
            let iso_stereo = stereo(heap, value.StereoIsotopic)?;
            let current_sp2 = i32::from(
                normal_stereo
                    .as_ref()
                    .is_some_and(|stereo| stereo.nNumberOfStereoBonds > 0),
            );
            let current_has_iso_stereo = i32::from(
                normal_stereo
                    .as_ref()
                    .is_some_and(|stereo| stereo.nNumberOfStereoCenters > 0),
            );
            let mut current_sp3 = current_has_iso_stereo;
            let mut current_iso_sp2 = i32::from(
                iso_h_stereo != 0
                    && iso_stereo
                        .as_ref()
                        .is_some_and(|stereo| stereo.nNumberOfStereoBonds > 0),
            );
            let mut current_iso_sp3 = i32::from(
                iso_h_stereo != 0
                    && iso_stereo
                        .as_ref()
                        .is_some_and(|stereo| stereo.nNumberOfStereoCenters > 0),
            );
            let mut current_iso_sp3_inverted = i32::from(
                current_iso_sp3 != 0
                    && iso_stereo
                        .as_ref()
                        .is_some_and(|stereo| stereo.nCompInv2Abs != 0),
            );
            requested_racemic |= i32::from(
                value.nFlags & u64::from(crate::source_types::INCHI_FLAG_RAC_STEREO) != 0,
            );
            requested_relative |= i32::from(
                value.nFlags & u64::from(crate::source_types::INCHI_FLAG_REL_STEREO) != 0,
            );
            if current_sp2 != 0 && current_iso_sp2 != 0 {
                current_iso_sp2 = i32::from(
                    crate::source::base::ichiprt2::Eql_INChI_Stereo(
                        heap,
                        normal_stereo.as_ref(),
                        crate::source_types::EQL_SP2 as i32,
                        iso_stereo.as_ref(),
                        crate::source_types::EQL_SP2 as i32,
                        0,
                    )? == 0,
                );
            }
            if current_sp3 != 0 && current_iso_sp3 != 0 {
                current_iso_sp3 = i32::from(
                    crate::source::base::ichiprt2::Eql_INChI_Stereo(
                        heap,
                        normal_stereo.as_ref(),
                        crate::source_types::EQL_SP3 as i32,
                        iso_stereo.as_ref(),
                        crate::source_types::EQL_SP3 as i32,
                        i32::from(
                            value.nFlags
                                & u64::from(
                                    crate::source_types::INCHI_FLAG_RAC_STEREO
                                        | crate::source_types::INCHI_FLAG_REL_STEREO,
                                )
                                != 0,
                        ),
                    )? == 0,
                );
                if current_iso_sp3 == 0 {
                    current_iso_sp3_inverted &= i32::from(
                        iso_stereo.as_ref().map(|value| value.nCompInv2Abs)
                            != normal_stereo.as_ref().map(|value| value.nCompInv2Abs),
                    );
                }
            }
            let current_relative = i32::from(requested_relative != 0 && current_sp3 != 0);
            let current_iso_relative = i32::from(
                requested_relative != 0 && (current_iso_sp3 != 0 || current_iso_sp3_inverted != 0),
            );
            let current_racemic = i32::from(requested_racemic != 0 && current_sp3 != 0);
            let current_iso_racemic = i32::from(
                requested_racemic != 0 && (current_iso_sp3 != 0 || current_iso_sp3_inverted != 0),
            );
            if requested_relative != 0 {
                current_sp3 = i32::from(
                    current_relative != 0
                        || (current_sp3 != 0
                            && normal_stereo
                                .as_ref()
                                .is_some_and(|value| value.nNumberOfStereoCenters > 1)),
                );
                current_iso_sp3 = if current_iso_relative != 0 {
                    current_iso_sp3
                } else {
                    0
                };
            } else if requested_racemic != 0 {
                current_sp3 = i32::from(
                    current_racemic > 1
                        || (current_sp3 != 0
                            && normal_stereo
                                .as_ref()
                                .is_some_and(|value| value.nNumberOfStereoCenters > 1)),
                );
                current_iso_sp3 = if current_iso_racemic > 1 {
                    current_iso_sp3
                } else {
                    0
                };
            }
            let target =
                usize::try_from(target).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            io.bIsotopic |= current_iso;
            has_isotopic_atoms[target] |= current_iso;
            stereo_sp2[target] |= current_sp2;
            stereo_sp3[target] |= current_sp3;
            io.bIgn_UU_Sp3[target] |= i32::from(
                current_sp3 == 0
                    && value.nFlags & u64::from(crate::source_types::INCHI_FLAG_SC_IGN_ALL_UU) != 0,
            );
            io.bIgn_UU_Sp2[target] |= i32::from(
                current_sp2 == 0
                    && value.nFlags & u64::from(crate::source_types::INCHI_FLAG_SB_IGN_ALL_UU) != 0,
            );
            isotopic_stereo_sp2[target] |= current_iso_sp2;
            isotopic_stereo_sp3[target] |= current_iso_sp3;
            io.bIgn_UU_Sp3_Iso[target] |= i32::from(
                current_iso_sp3 == 0
                    && value.nFlags & u64::from(crate::source_types::INCHI_FLAG_SC_IGN_ALL_ISO_UU)
                        != 0,
            );
            io.bIgn_UU_Sp2_Iso[target] |= i32::from(
                current_iso_sp2 == 0
                    && value.nFlags & u64::from(crate::source_types::INCHI_FLAG_SB_IGN_ALL_ISO_UU)
                        != 0,
            );
            stereo_abs[target] |= i32::from(
                current_sp3 != 0
                    && normal_stereo
                        .as_ref()
                        .is_some_and(|value| value.nCompInv2Abs != 0),
            );
            stereo_abs_inverted[target] |= i32::from(
                current_sp3 != 0
                    && normal_stereo
                        .as_ref()
                        .is_some_and(|value| value.nCompInv2Abs < 0),
            );
            let inversion_differs =
                iso_stereo
                    .as_ref()
                    .zip(normal_stereo.as_ref())
                    .is_some_and(|(iso, normal)| {
                        iso.nCompInv2Abs != 0 && iso.nCompInv2Abs != normal.nCompInv2Abs
                    });
            isotopic_stereo_abs_inverted[target] |= i32::from(
                (current_iso_sp3 != 0
                    && iso_stereo
                        .as_ref()
                        .is_some_and(|value| value.nCompInv2Abs < 0))
                    || (current_iso_sp3 == 0 && inversion_differs),
            );
            isotopic_stereo_abs[target] |= i32::from(
                (current_iso_sp3 != 0
                    && iso_stereo
                        .as_ref()
                        .is_some_and(|value| value.nCompInv2Abs != 0))
                    || (current_iso_sp3 == 0 && inversion_differs),
            );
            io.bRelativeStereo[target] |= current_relative;
            io.bIsotopicRelativeStereo[target] |= current_iso_relative;
            io.bRacemicStereo[target] |= current_racemic;
            io.bIsotopicRacemicStereo[target] |= current_iso_racemic;
            tautomeric_acid |=
                i32::from(value.nFlags & u64::from(crate::source_types::INCHI_FLAG_ACID_TAUT) != 0);
            hard_add_remove_proton |= i32::from(
                value.nFlags & u64::from(crate::source_types::INCHI_FLAG_HARD_ADD_REM_PROTON) != 0,
            );
            if current_taut != 0 {
                io.bTautomeric |= 1;
                io.bNonTautomeric |= has_fixed(heap, &sort)?;
            }
            if output_options & crate::source_types::INCHI_OUT_NO_AUX_INFO as i32 == 0 {
                if let Some(value_aux) = aux(heap, sort.pINChI_Aux[found as usize])? {
                    io.bAtomEqu[target] |= crate::source::base::ichiprt2::bHasEquString(
                        heap,
                        value_aux.nConstitEquNumbers.as_const(),
                        value_aux.nNumberOfAtoms,
                    )?;
                    if current_taut != 0 {
                        io.bTautEqu[target] |= crate::source::base::ichiprt2::bHasEquString(
                            heap,
                            value_aux.nConstitEquTGroupNumbers.as_const(),
                            value_aux.nNumberOfTGroups,
                        )?;
                    }
                    if current_iso != 0 {
                        io.bIsotopicAtomEqu[target] |=
                            crate::source::base::ichiprt2::bHasEquString(
                                heap,
                                value_aux.nConstitEquIsotopicNumbers.as_const(),
                                value_aux.nNumberOfAtoms,
                            )?;
                        if current_taut != 0 {
                            io.bIsotopicTautEqu[target] |=
                                crate::source::base::ichiprt2::bHasEquString(
                                    heap,
                                    value_aux.nConstitEquIsotopicTGroupNumbers.as_const(),
                                    value_aux.nNumberOfTGroups,
                                )?;
                        }
                        io.bIsotopicOrigNumb[target] |= current_has_iso_stereo
                            & pointers_differ(
                                heap,
                                value_aux.nOrigAtNosInCanonOrdInv,
                                value_aux.nIsotopicOrigAtNosInCanonOrd,
                                value_aux.nNumberOfAtoms,
                            )?;
                    }
                    if current_sp3 != 0
                        && normal_stereo
                            .as_ref()
                            .is_some_and(|value| value.nCompInv2Abs != 0)
                    {
                        io.bInvStereo[target] |= 1;
                        io.bInvStereoOrigNumb[target] |= pointers_differ(
                            heap,
                            value_aux.nOrigAtNosInCanonOrd,
                            value_aux.nOrigAtNosInCanonOrdInv,
                            value_aux.nNumberOfAtoms,
                        )?;
                    }
                    if current_iso_sp3 != 0
                        && iso_stereo
                            .as_ref()
                            .is_some_and(|value| value.nCompInv2Abs != 0)
                    {
                        io.bInvIsotopicStereo[target] |= 1;
                        io.bInvIsotopicStereoOrigNumb[target] |= pointers_differ(
                            heap,
                            value_aux.nIsotopicOrigAtNosInCanonOrd,
                            value_aux.nIsotopicOrigAtNosInCanonOrdInv,
                            value_aux.nNumberOfAtoms,
                        )?;
                    }
                    if crate::source::base::ichiprt2::bHasOrigInfo(
                        heap,
                        value_aux.OrigInfo.as_const(),
                        value_aux.nNumberOfAtoms,
                    )? != 0
                    {
                        io.bChargesRadVal[target] |= 1;
                    }
                }
            }
        }
        if component_exists != 0 {
            io.num_comp[0] = io.num_comp[0].wrapping_add(1);
            io.num_comp[1] = io.num_comp[1].wrapping_add(1);
        }
    }
    if io.bTautomeric != 0 {
        io.bTautomeric = io.bTautomeric.wrapping_add(tautomeric_acid);
        io.bTautomeric =
            io.bTautomeric
                .wrapping_add(if hard_add_remove_proton != 0 { 4 } else { 0 });
    }
    if requested_racemic != 0 || requested_relative != 0 {
        for index in 0..2 {
            stereo_abs_inverted[index] = 0;
            stereo_abs[index] = 0;
            io.bInvStereo[index] = 0;
            io.bInvStereoOrigNumb[index] = 0;
            isotopic_stereo_abs_inverted[index] = 0;
            isotopic_stereo_abs[index] = 0;
            io.bInvIsotopicStereo[index] = 0;
            io.bInvIsotopicStereoOrigNumb[index] = 0;
        }
    }
    io.iCurTautMode = match io.bOutType {
        value if value == crate::source_types::OUT_N1 as i32 => 0,
        value if value == crate::source_types::OUT_T1 as i32 => 1,
        value if value == crate::source_types::OUT_NT as i32 => 0,
        value if value == crate::source_types::OUT_TN as i32 => 1,
        _ => -1,
    };
    if io.iCurTautMode < 0 {
        return Ok(0);
    }
    io.num_components = io.num_comp[io.iCurTautMode as usize];
    let component_counts_pointer = heap.allocate_model_storage(component_counts.to_vec())?;

    let mut intermediate_result = 0_i32;
    let mut embedded_output_called = false;
    if output_options & crate::source_types::INCHI_OUT_ONLY_AUX_INFO as i32 == 0 {
        io.nCurINChISegment = crate::source_types::tagDiffINChILayers_DIFL_M as i32;
        if basic_or_reconnected == crate::source_types::INCHI_BAS as i32
            || output_options & crate::source_types::INCHI_OUT_EMBED_REC as i32 == 0
        {
            let structure = if original_structure.is_null() {
                None
            } else {
                Some(
                    heap.slice(original_structure)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                )
            };
            let is_beta = i32::from(
                all_fixed_atoms.max(all_mobile_atoms) > 1024
                    || structure
                        .as_ref()
                        .is_some_and(|value| !value.polymer.is_null() || value.n_zy != 0)
                    || (structure.is_none()
                        && (input_parameters.bLargeMolecules != 0
                            || input_parameters.bPolymers != 0)),
            );
            OutputINCHI_VersionAndKind(
                heap,
                output,
                string_buffer,
                output_options,
                is_beta,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            )?;
        }
        intermediate_result = OutputINCHI_MainLayerFormula(
            heap,
            canonical_globals,
            output,
            string_buffer,
            component_counts,
            &basic_or_reconnected,
            &mut io,
            line_feed.as_const(),
            tab.as_const(),
            stdout,
        )?;
        if intermediate_result == 0 {
            intermediate_result = OutputINCHI_MainLayerConnections(
                heap,
                canonical_globals,
                output,
                string_buffer,
                component_counts,
                &basic_or_reconnected,
                &mut io,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            )?;
        }
        if intermediate_result == 0 {
            intermediate_result = OutputINCHI_MainLayerHydrogens(
                heap,
                canonical_globals,
                output,
                string_buffer,
                component_counts,
                &basic_or_reconnected,
                &mut io,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            )?;
        }
        io.bFhTag = 0;
        let mut pass = 0_i32;
        while intermediate_result == 0 {
            intermediate_result = OutputINCHI_ChargeAndRemovedAddedProtonsLayers(
                heap,
                canonical_globals,
                output,
                string_buffer,
                &mut io,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            )?;
            if intermediate_result != 0 {
                break;
            }
            if pass == 0 {
                intermediate_result = OutputINCHI_PolymerLayer(
                    heap,
                    canonical_globals,
                    Some(&mut *output),
                    string_buffer,
                    &mut basic_or_reconnected,
                    original_atom_data,
                    original_structure,
                    &mut io,
                    line_feed.as_const(),
                    tab.as_const(),
                    stdout,
                )?;
                if intermediate_result != 0 {
                    break;
                }
            }
            intermediate_result = OutputINCHI_StereoLayer(
                heap,
                canonical_globals,
                output,
                string_buffer,
                &mut io,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            )?;
            if intermediate_result != 0 {
                break;
            }
            io.nCurINChISegment = io.nCurINChISegment.wrapping_add(1);
            intermediate_result = OutputINCHI_IsotopicLayer(
                heap,
                canonical_globals,
                output,
                string_buffer,
                &basic_or_reconnected,
                &mut io,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            )?;
            if intermediate_result != 0 {
                break;
            }
            let mut repeat = 0_i32;
            intermediate_result = OutputINCHI_FixedHLayerWithSublayers(
                heap,
                canonical_globals,
                output,
                string_buffer,
                &basic_or_reconnected,
                &mut io,
                line_feed.as_const(),
                tab.as_const(),
                &mut repeat,
                stdout,
            )?;
            if intermediate_result != 0 || repeat == 0 {
                break;
            }
            pass = pass.wrapping_add(1);
        }
        if intermediate_result == 0
            && disconnected_coordinates != 0
            && basic_or_reconnected == crate::source_types::INCHI_BAS as i32
            && output_options & crate::source_types::INCHI_OUT_EMBED_REC as i32 != 0
            && component_counts[crate::source_types::INCHI_REC as usize] != 0
        {
            embedded_output_called = true;
            print_values(
                heap,
                output,
                stdout,
                false,
                b"%s\0",
                vec![SourceFormatArgument::Bytes(line_feed.as_const())],
            )?;
            let nested = OutputINChI1(
                heap,
                canonical_globals,
                Some(&mut *string_buffer),
                sort_rows,
                crate::source_types::INCHI_REC as i32,
                original_atom_data,
                original_structure,
                input_parameters,
                0,
                output_type,
                output_options | crate::source_types::INCHI_OUT_NO_AUX_INFO as i32,
                component_counts,
                non_tautomeric_counts,
                tautomeric_counts,
                output,
                log,
                input_structure_number,
                sort_print_flags,
                save_option_bits,
                stdout,
            )?;
            if nested == 0 {
                intermediate_result = 1;
            }
        }
        if intermediate_result == 0
            && !embedded_output_called
            && output_options & crate::source_types::INCHI_OUT_SAVEOPT as i32 != 0
            && output_options & crate::source_types::INCHI_OUT_STDINCHI as i32 == 0
        {
            let mut first = 0_i8;
            let mut second = 0_i8;
            GetSaveOptLetters(save_option_bits, &mut first, &mut second);
            print_values(
                heap,
                output,
                stdout,
                false,
                b"\\%c%c\0",
                vec![
                    SourceFormatArgument::Signed(i64::from(first)),
                    SourceFormatArgument::Signed(i64::from(second)),
                ],
            )?;
        }
        if intermediate_result == 0 && !embedded_output_called && plain_comments == 0 {
            let empty = heap.allocate_model_storage(
                if component_counts[0] == 0 && component_counts[1] == 0 {
                    vec![b'/' as i8, b'/' as i8, 0]
                } else {
                    vec![0_i8]
                },
            )?;
            let ending = if output_options & crate::source_types::INCHI_OUT_NO_AUX_INFO as i32 != 0
            {
                heap.allocate_model_storage(vec![b'\n' as i8, 0])?
            } else {
                tab
            };
            print_values(
                heap,
                output,
                stdout,
                false,
                b"%s%s\0",
                vec![
                    SourceFormatArgument::Bytes(empty.as_const()),
                    SourceFormatArgument::Bytes(ending.as_const()),
                ],
            )?;
        }
        inchi_strbuf_reset(heap, Some(&mut *string_buffer))?;
    }

    if intermediate_result == 0 {
        io.bFhTag = 0;
        if output_options & crate::source_types::INCHI_OUT_NO_AUX_INFO as i32 == 0 {
            io.num_components = io.num_comp[io.iCurTautMode as usize];
            intermediate_result = OutputAUXINFO_HeaderAndNormalization_type(
                heap,
                canonical_globals,
                Some(&mut *output),
                string_buffer,
                output_options,
                &mut basic_or_reconnected,
                component_counts_pointer.as_const(),
                &mut io,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            )?;
            let mut aux_pass = 0_i32;
            let mut iso_only = false;
            while intermediate_result == 0 {
                if !iso_only {
                    intermediate_result = OutputAUXINFO_OriginalNumbersAndEquivalenceClasses(
                        heap,
                        canonical_globals,
                        Some(&mut *output),
                        string_buffer,
                        component_counts_pointer.as_const(),
                        &mut io,
                        line_feed.as_const(),
                        tab.as_const(),
                        stdout,
                    )?;
                    if intermediate_result != 0 {
                        break;
                    }
                    intermediate_result = OutputAUXINFO_TautomericGroupsEquivalence(
                        heap,
                        canonical_globals,
                        Some(&mut *output),
                        string_buffer,
                        &mut io,
                        stdout,
                    )?;
                    if intermediate_result != 0 {
                        break;
                    }
                    intermediate_result = OutputAUXINFO_Stereo(
                        heap,
                        canonical_globals,
                        Some(&mut *output),
                        string_buffer,
                        &mut io,
                        line_feed.as_const(),
                        tab.as_const(),
                        stdout,
                    )?;
                    if intermediate_result != 0 {
                        break;
                    }
                }
                intermediate_result = OutputAUXINFO_IsotopicInfo(
                    heap,
                    canonical_globals,
                    Some(&mut *output),
                    string_buffer,
                    &mut basic_or_reconnected,
                    &mut io,
                    line_feed.as_const(),
                    tab.as_const(),
                    stdout,
                )?;
                if intermediate_result != 0 {
                    break;
                }
                let flags = *heap
                    .slice(sort_print_flags.as_const())?
                    .first()
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let no_non_iso = if basic_or_reconnected == crate::source_types::INCHI_BAS as i32 {
                    crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_BAS as i32
                } else {
                    crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_REC as i32
                };
                let no_iso = if basic_or_reconnected == crate::source_types::INCHI_BAS as i32 {
                    crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_BAS as i32
                } else {
                    crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_REC as i32
                };
                if io.bOutType == crate::source_types::OUT_TN as i32
                    && io.bTautomeric != 0
                    && io.bNonTautomeric != 0
                    && flags & no_non_iso != 0
                    && flags & no_iso != 0
                {
                    io.bNonTautomeric = 0;
                }
                if io.bOutType == crate::source_types::OUT_TN as i32
                    && io.bTautomeric != 0
                    && io.bNonTautomeric != 0
                {
                    io.bOutType = crate::source_types::OUT_NN as i32;
                    io.iCurTautMode = crate::source_types::TAUT_NON as i32;
                    io.pINChISort = fixed_sorts;
                    io.bSecondNonTautPass = 1;
                    io.num_components = io.num_comp[crate::source_types::TAUT_NON as usize];
                    io.bFhTag = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_FIXH as i32;
                    inchi_strbuf_reset(heap, Some(&mut *string_buffer))?;
                    aux_pass = aux_pass.wrapping_add(1);
                    iso_only = flags & no_non_iso != 0;
                    continue;
                }
                if io.bOutType == crate::source_types::OUT_NN as i32
                    && io.bOutputType == crate::source_types::OUT_TN as i32
                    && io.bTautomeric != 0
                    && io.bNonTautomeric != 0
                {
                    io.bOutType = crate::source_types::OUT_TN as i32;
                    io.iCurTautMode = crate::source_types::TAUT_YES as i32;
                    io.pINChISort = mobile_sorts;
                    io.bSecondNonTautPass = 0;
                    io.num_components = io.num_comp[crate::source_types::TAUT_YES as usize];
                    io.bFhTag = 0;
                }
                let _ = aux_pass;
                break;
            }
            if intermediate_result == 0 {
                intermediate_result = OutputAUXINFO_ChargesRadicalsAndUnusualValences(
                    heap,
                    canonical_globals,
                    Some(&mut *output),
                    string_buffer,
                    &mut io,
                    line_feed.as_const(),
                    tab.as_const(),
                    stdout,
                )?;
            }
            let structure = if original_structure.is_null() {
                None
            } else {
                Some(
                    heap.slice(original_structure)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        .clone(),
                )
            };
            if intermediate_result == 0 {
                intermediate_result = OutputAUXINFO_ReversibilityInfo(
                    heap,
                    canonical_globals,
                    Some(&mut *output),
                    string_buffer,
                    structure.as_ref(),
                    &mut io,
                    line_feed.as_const(),
                    tab.as_const(),
                    stdout,
                )?;
            }
            if intermediate_result == 0 {
                intermediate_result = OutputAUXINFO_PolymerInfo(
                    heap,
                    canonical_globals,
                    Some(&mut *output),
                    string_buffer,
                    structure.as_ref(),
                    &mut io,
                    line_feed.as_const(),
                    tab.as_const(),
                    stdout,
                )?;
            }
            if intermediate_result == 0
                && disconnected_coordinates != 0
                && basic_or_reconnected == crate::source_types::INCHI_BAS as i32
                && output_options & crate::source_types::INCHI_OUT_EMBED_REC as i32 != 0
                && component_counts[crate::source_types::INCHI_REC as usize] != 0
            {
                embedded_output_called = true;
                print_values(
                    heap,
                    output,
                    stdout,
                    true,
                    b"%s\0",
                    vec![SourceFormatArgument::Bytes(line_feed.as_const())],
                )?;
                let nested = OutputINChI1(
                    heap,
                    canonical_globals,
                    Some(&mut *string_buffer),
                    sort_rows,
                    crate::source_types::INCHI_REC as i32,
                    SourceConstPointer::null(),
                    SourceConstPointer::null(),
                    input_parameters,
                    0,
                    output_type,
                    crate::source_types::INCHI_OUT_ONLY_AUX_INFO as i32 | output_options,
                    component_counts,
                    non_tautomeric_counts,
                    tautomeric_counts,
                    output,
                    log,
                    input_structure_number,
                    sort_print_flags,
                    save_option_bits,
                    stdout,
                )?;
                if nested == 0 {
                    intermediate_result = 1;
                }
            }
            if intermediate_result == 0 && !embedded_output_called && plain_comments == 0 {
                let empty = heap.allocate_model_storage(
                    if component_counts[0] == 0 && component_counts[1] == 0 {
                        vec![b'/' as i8, b'/' as i8, 0]
                    } else {
                        vec![0_i8]
                    },
                )?;
                print_values(
                    heap,
                    output,
                    stdout,
                    true,
                    b"%s\n\0",
                    vec![SourceFormatArgument::Bytes(empty.as_const())],
                )?;
            }
            if intermediate_result == 0 {
                let ending = heap.allocate_model_storage(
                    if output_options & crate::source_types::INCHI_OUT_WINCHI_WINDOW as i32 != 0 {
                        vec![b'\n' as i8, 0]
                    } else {
                        vec![0_i8]
                    },
                )?;
                print_values(
                    heap,
                    output,
                    stdout,
                    true,
                    b"%s\0",
                    vec![SourceFormatArgument::Bytes(ending.as_const())],
                )?;
            }
        }
    }

    let mut result = i32::from(intermediate_result == 0);
    if io.bOverflow != 0 {
        print_values(
            heap,
            output,
            stdout,
            true,
            b"\nFATAL ERROR: Output buffer overflow\n\0",
            Vec::new(),
        )?;
    }
    if intermediate_result != 0 {
        result = 0;
        report_serialization_error(heap, log, input_parameters, input_structure_number, false)?;
    }
    Ok(result)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_MainLayerFormula(
    heap: &mut SourceHeap,
    _canon_globals: SourceMutPointer<CANON_GLOBALS>,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    number_of_components: &[i32; 2],
    basic_or_reconnected: &i32,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3170 OutputINCHI_MainLayerFormula
    // INCHI✔️❌: int OutputINCHI_MainLayerFormula( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                   INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                   INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                   int              num_components2[],
    // INCHI✔️❌:                                   int              *INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                   INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                   char             *pLF,
    // INCHI✔️❌:                                   char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     /* constitution ( dot-disconnected Hill formulas: <formula> ) */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (num_components2[0] || num_components2[1])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szGetTag( IdentLbl, io->nTag, io->bTag1 = *INCHI_basic_or_INCHI_reconnected == INCHI_REC ? IL_REC_ : IL_FML_, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:         io->tot_len = str_HillFormula( io->pINChISort, strbuf, &io->bOverflow, io->bOutType,
    // INCHI✔️❌:                                    io->num_components, io->bUseMulipliers );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, 1 ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (io->n_pzz > 0 && io->n_zy > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int retm = MergeZzInHillFormula(strbuf);
    // INCHI✔️❌:             if (0 != retm)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return -1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     LOG_NO_ARGS("\n#################### (L3318:ichiprt1.c) ##########################\n");
    // INCHI✔️❌:     LOG_MULT_ARGS("This is the Chemical formula : %s\n", strbuf->pStr);
    // INCHI✔️❌:     LOG_NO_ARGS("####################################################################\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_MainLayerFormula

    if number_of_components[0] != 0 || number_of_components[1] != 0 {
        io.bTag1 = if *basic_or_reconnected == INCHI_REC as i32 {
            1 << 18
        } else {
            1 << 4
        };
        let tags = formula_ident_labels(heap)?;
        let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
        szGetTag(heap, &tags, io.nTag, io.bTag1, tag, &mut io.bAlways, 1)?;
        io.szTag1.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = str_HillFormula(
            heap,
            io.pINChISort,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.num_components,
            io.bUseMulipliers,
        )?;
        if str_LineEnd(
            heap,
            tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            1,
        )? != 0
        {
            return Ok(1);
        }
        if io.n_pzz > 0 && io.n_zy > 0 && MergeZzInHillFormula(heap, string_buffer)? != 0 {
            return Ok(-1);
        }
        let format =
            heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
                ..SourceVaList::default()
            },
        )?;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_MainLayerConnections(
    heap: &mut SourceHeap,
    canon_globals: SourceMutPointer<CANON_GLOBALS>,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    _number_of_components: &[i32; 2],
    _basic_or_reconnected: &i32,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3213 OutputINCHI_MainLayerConnections
    // INCHI✔️❌: int OutputINCHI_MainLayerConnections( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                       INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                       INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                       int              num_components2[],
    // INCHI✔️❌:                                       int              *INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                       INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                       char             *pLF,
    // INCHI✔️❌:                                       char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* connections ( semicolon/dot-disconnected connection tables ) */
    // INCHI✔️❌:
    // INCHI✔️❌:     szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_CONN, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:     inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:     io->tot_len = 0;
    // INCHI✔️❌:     io->tot_len2 = str_Connections( pCG, io->pINChISort, strbuf, &io->bOverflow, io->bOutType,
    // INCHI✔️❌:                                     io->ATOM_MODE, io->num_components, io->bUseMulipliers );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* current version does not output empty (";;;;") connectivity */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (io->tot_len != io->tot_len2)
    // INCHI✔️❌:     { /* 2004-06-30: never output empty connection table */
    // INCHI✔️❌:         io->tot_len = io->tot_len2;
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -2, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1; /* pStr overfow */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     LOG_NO_ARGS("\n##################### (L3357:ichiprt1.c) #########################\n");
    // INCHI✔️❌:     LOG_MULT_ARGS("This is the Connection Layer : %s\n", strbuf->pStr);
    // INCHI✔️❌:     LOG_NO_ARGS("####################################################################\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_MainLayerConnections

    io.bTag1 = 1 << 5;
    let tags = connections_ident_labels(heap)?;
    let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
    szGetTag(heap, &tags, io.nTag, io.bTag1, tag, &mut io.bAlways, 1)?;
    io.szTag1.copy_from_slice(
        heap.slice(tag.as_const())?
            .get(..64)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    inchi_strbuf_reset(heap, Some(string_buffer))?;
    io.tot_len = 0;
    io.tot_len2 = crate::source::base::ichiprt3::str_Connections(
        heap,
        io.pINChISort,
        string_buffer,
        &mut io.bOverflow,
        io.bOutType,
        io.ATOM_MODE,
        io.num_components,
        io.bUseMulipliers,
    )?;
    if io.tot_len != io.tot_len2 {
        io.tot_len = io.tot_len2;
        if str_LineEnd(
            heap,
            tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -2,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        let format =
            heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
                ..SourceVaList::default()
            },
        )?;
    }
    let _ = canon_globals;
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_MainLayerHydrogens(
    heap: &mut SourceHeap,
    _canon_globals: SourceMutPointer<CANON_GLOBALS>,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    _number_of_components: &[i32; 2],
    _basic_or_reconnected: &i32,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3251 OutputINCHI_MainLayerHydrogens
    // INCHI✔️❌: int OutputINCHI_MainLayerHydrogens( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                     INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                     INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                     int              num_components2[],
    // INCHI✔️❌:                                     int              *INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                     INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                     char             *pLF,
    // INCHI✔️❌:                                     char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     /* hydrogen atoms (do not output empty) */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (INCHI_SEGM_FILL == INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_h_H_ATOMS] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_ALLH, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:         io->tot_len = 0;
    // INCHI✔️❌:         io->tot_len2 = str_H_atoms( io->pINChISort, strbuf, &io->bOverflow, io->bOutType,
    // INCHI✔️❌:                                 io->ATOM_MODE, io->TAUT_MODE,
    // INCHI✔️❌:                                 io->num_components, io->bUseMulipliers );
    // INCHI✔️❌:         if (io->tot_len != io->tot_len2)
    // INCHI✔️❌:         { /* 2004-06-21: never output empty */
    // INCHI✔️❌:             io->tot_len = io->tot_len2;
    // INCHI✔️❌:             if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -2, 1 ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     LOG_NO_ARGS("\n###################### (L3396:ichiprt1.c) ########################\n");
    // INCHI✔️❌:     LOG_MULT_ARGS("This is the Hydrogen Layer : %s\n", strbuf->pStr);
    // INCHI✔️❌:     LOG_NO_ARGS("####################################################################\n");
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_MainLayerHydrogens

    let layer =
        usize::try_from(io.nCurINChISegment).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let segment = crate::source_types::tagDiffINChISegments_DIFS_h_H_ATOMS as usize;
    let difference = *io
        .sDifSegs
        .get(layer)
        .and_then(|segments| segments.get(segment))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if crate::source::base::ichimake::INChI_SegmentAction(difference)
        == crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32
    {
        io.bTag1 = crate::source_types::local_ichiprt1::tagIdentLblBit_IL_ALLH as i32;
        let tags = hydrogens_ident_labels(heap)?;
        let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
        szGetTag(heap, &tags, io.nTag, io.bTag1, tag, &mut io.bAlways, 1)?;
        io.szTag1.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len2 = crate::source::base::ichiprt3::str_H_atoms(
            heap,
            io.pINChISort,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.ATOM_MODE,
            io.TAUT_MODE,
            io.num_components,
            io.bUseMulipliers,
        )?;
        if io.tot_len != io.tot_len2 {
            io.tot_len = io.tot_len2;
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                -2,
                1,
            )? != 0
            {
                return Ok(1);
            }
            let format =
                heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
            inchi_ios_print_nodisplay(
                heap,
                Some(output),
                stdout,
                format.as_const(),
                &SourceVaList {
                    arguments: vec![
                        SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                        SourceFormatArgument::Bytes(line_feed),
                    ],
                    ..SourceVaList::default()
                },
            )?;
        }
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_ChargeAndRemovedAddedProtonsLayers(
    heap: &mut SourceHeap,
    _canon_globals: SourceMutPointer<CANON_GLOBALS>,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3293 OutputINCHI_ChargeAndRemovedAddedProtonsLayers
    // INCHI✔️❌: int OutputINCHI_ChargeAndRemovedAddedProtonsLayers( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                                     INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                                     INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                                     INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                                     char             *pLF,
    // INCHI✔️❌:                                                     char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     /* charge  */
    // INCHI✔️❌:
    // INCHI✔️❌:     io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_q_CHARGE] );
    // INCHI✔️❌:     if (io->nSegmAction)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_CHRG | io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:         io->tot_len = 0;
    // INCHI✔️❌:         if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             io->tot_len = str_Charge2( io->pINChISort, io->pINChISort2,
    // INCHI✔️❌:                                    strbuf, &io->bOverflow, io->bOutType, io->num_components,
    // INCHI✔️❌:                                    io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:             io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* removed protons */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (io->iCurTautMode == TAUT_YES && !io->bSecondNonTautPass)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_p_PROTONS] );
    // INCHI✔️❌:         if (io->nSegmAction)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_PROT | io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:             io->tot_len = 0;
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%+d", io->nNumRemovedProtons );
    // INCHI✔️❌:             if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "/" );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_ChargeAndRemovedAddedProtonsLayers

    fn selected_difference(io: &INCHI_OUT_CTL, segment: usize) -> Result<i8, SourceHeapError> {
        let layer = usize::try_from(io.nCurINChISegment)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        io.sDifSegs
            .get(layer)
            .and_then(|segments| segments.get(segment))
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    fn set_tag(
        heap: &mut SourceHeap,
        io: &mut INCHI_OUT_CTL,
        tags: &[INCHI_TAG],
        tag_bit: i32,
    ) -> Result<SourceMutPointer<i8>, SourceHeapError> {
        io.bTag1 = tag_bit | io.bFhTag;
        let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
        szGetTag(heap, tags, io.nTag, io.bTag1, tag, &mut io.bAlways, 1)?;
        io.szTag1.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        Ok(tag)
    }

    fn print_buffer(
        heap: &mut SourceHeap,
        output: &mut INCHI_IOSTREAM,
        string_buffer: &INCHI_IOS_STRING,
        line_feed: SourceConstPointer<i8>,
        stdout: SourceMutPointer<FILE>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    let tags = charge_proton_ident_labels(heap)?;
    io.nSegmAction = crate::source::base::ichimake::INChI_SegmentAction(selected_difference(
        io,
        crate::source_types::tagDiffINChISegments_DIFS_q_CHARGE as usize,
    )?);
    if io.nSegmAction != 0 {
        let tag = set_tag(
            heap,
            io,
            &tags,
            crate::source_types::local_ichiprt1::tagIdentLblBit_IL_CHRG as i32,
        )?;
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        if io.nSegmAction == crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32 {
            io.tot_len = crate::source::base::ichiprt3::str_Charge2(
                heap,
                io.pINChISort,
                io.pINChISort2,
                string_buffer,
                &mut io.bOverflow,
                io.bOutType,
                io.num_components,
                io.bSecondNonTautPass,
                io.bOmitRepetitions,
                io.bUseMulipliers,
            )?;
            io.bNonTautNonIsoIdentifierNotEmpty = io
                .bNonTautNonIsoIdentifierNotEmpty
                .wrapping_add(io.bSecondNonTautPass);
        }
        if str_LineEnd(
            heap,
            tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            io.nSegmAction.wrapping_neg(),
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_buffer(heap, output, string_buffer, line_feed, stdout)?;
    }

    if io.iCurTautMode == crate::source_types::TAUT_YES as i32 && io.bSecondNonTautPass == 0 {
        io.nSegmAction = crate::source::base::ichimake::INChI_SegmentAction(selected_difference(
            io,
            crate::source_types::tagDiffINChISegments_DIFS_p_PROTONS as usize,
        )?);
        if io.nSegmAction != 0 {
            let tag = set_tag(
                heap,
                io,
                &tags,
                crate::source_types::local_ichiprt1::tagIdentLblBit_IL_PROT as i32,
            )?;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            let format =
                heap.allocate_model_storage(b"%+d\0".iter().map(|byte| *byte as i8).collect())?;
            match inchi_strbuf_printf(
                heap,
                Some(string_buffer),
                format.as_const(),
                &SourceVaList {
                    arguments: vec![SourceFormatArgument::Signed(i64::from(
                        io.nNumRemovedProtons,
                    ))],
                    ..SourceVaList::default()
                },
            ) {
                Ok(_) | Err(SourceHeapError::AllocationFailed) => {}
                Err(error) => return Err(error),
            }
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                io.nSegmAction.wrapping_neg(),
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(1);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
        } else if io.bPlainTextTags == 1 {
            let slash = heap.allocate_model_storage(vec![b'/' as i8, 0])?;
            inchi_ios_print_nodisplay(
                heap,
                Some(output),
                stdout,
                slash.as_const(),
                &SourceVaList::default(),
            )?;
        }
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_StereoLayer(
    heap: &mut SourceHeap,
    _canon_globals: SourceMutPointer<CANON_GLOBALS>,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3354 OutputINCHI_StereoLayer
    // INCHI✔️❌: int OutputINCHI_StereoLayer( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                              INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                              INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                              INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                              char             *pLF,
    // INCHI✔️❌:                              char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int i;
    // INCHI✔️❌:         i = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ) ||
    // INCHI✔️❌:          INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ) ||
    // INCHI✔️❌:          INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ) ||
    // INCHI✔️❌:          INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  stereo */
    // INCHI✔️❌:
    // INCHI✔️❌:         szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_STER | io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  sp2 */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*if ( bStereoSp2[io->iCurTautMode]  )*/
    // INCHI✔️❌:         if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_DBND, io->szTag2, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:             io->tot_len = 0;
    // INCHI✔️❌:             if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 io->tot_len = str_Sp2( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
    // INCHI✔️❌:                                         io->bOutType, io->TAUT_MODE, io->num_components,
    // INCHI✔️❌:                                         io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:
    // INCHI✔️❌:                 io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io->bPlainTextTags == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "/" ); /* sp2 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  sp3 */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*if ( bStereoSp3[io->iCurTautMode]  )*/
    // INCHI✔️❌:         if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             io->bRelRac = io->bRelativeStereo[io->iCurTautMode] || io->bRacemicStereo[io->iCurTautMode];
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_SP3S, io->szTag2, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:             io->tot_len = 0;
    // INCHI✔️❌:             if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 io->tot_len = str_Sp3( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
    // INCHI✔️❌:                                        io->bOutType, io->TAUT_MODE, io->num_components, io->bRelRac,
    // INCHI✔️❌:                                    io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:
    // INCHI✔️❌:                 io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "/" ); /* sp3 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* bStereoAbsInverted[io->iCurTautMode]  */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* if ( bStereoAbs[io->iCurTautMode]  ) */
    // INCHI✔️❌:         if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_INVS, io->szTag2, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:             if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 io->tot_len = str_StereoAbsInv( io->pINChISort, strbuf,
    // INCHI✔️❌:                                             &io->bOverflow, io->bOutType, io->num_components );
    // INCHI✔️❌:                 io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 3;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io->bPlainTextTags == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "/" ); /* stereo-abs-inv */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* stereo type */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*if ( io->bRacemicStereo[io->iCurTautMode] || io->bRelativeStereo[io->iCurTautMode] || bStereoAbs[io->iCurTautMode] )*/
    // INCHI✔️❌:         if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             const char *p_stereo = io->bRelativeStereo[io->iCurTautMode] ? x_rel :
    // INCHI✔️❌:                 io->bRacemicStereo[io->iCurTautMode] ? x_rac : x_abs;
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_TYPS, io->szTag2, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:             if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ( io->tot_len ) += MakeDelim( p_stereo, strbuf, &io->bOverflow );
    // INCHI✔️❌:                 io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (io->bPlainTextTags == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "/" );  /* no abs, inv or racemic stereo */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "////" ); /* sp3, sp2, abs-inv, stereo.type */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_StereoLayer
    fn selected_difference(io: &INCHI_OUT_CTL, segment: usize) -> Result<i8, SourceHeapError> {
        let layer = usize::try_from(io.nCurINChISegment)
            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        io.sDifSegs
            .get(layer)
            .and_then(|segments| segments.get(segment))
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    fn make_tag(
        heap: &mut SourceHeap,
        tags: &[INCHI_TAG],
        n_tag: i32,
        bits: i32,
        source: &[i8; 64],
        always: &mut i32,
    ) -> Result<(SourceMutPointer<i8>, [i8; 64]), SourceHeapError> {
        let tag = heap.allocate_model_storage(source.to_vec())?;
        szGetTag(heap, tags, n_tag, bits, tag, always, 1)?;
        let mut value = [0_i8; 64];
        value.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        Ok((tag, value))
    }

    fn print_buffer(
        heap: &mut SourceHeap,
        output: &mut INCHI_IOSTREAM,
        string_buffer: &INCHI_IOS_STRING,
        line_feed: SourceConstPointer<i8>,
        stdout: SourceMutPointer<FILE>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    fn print_literal(
        heap: &mut SourceHeap,
        output: &mut INCHI_IOSTREAM,
        literal: &[u8],
        stdout: SourceMutPointer<FILE>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(literal.iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList::default(),
        )?;
        Ok(())
    }

    let bond_segment = crate::source_types::tagDiffINChISegments_DIFS_b_SBONDS as usize;
    let center_segment = crate::source_types::tagDiffINChISegments_DIFS_t_SATOMS as usize;
    let inversion_segment = crate::source_types::tagDiffINChISegments_DIFS_m_SP3INV as usize;
    let type_segment = crate::source_types::tagDiffINChISegments_DIFS_s_STYPE as usize;
    let fill = crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32;
    let segment_action = |value| crate::source::base::ichimake::INChI_SegmentAction(value);

    let _ = segment_action(selected_difference(io, center_segment)?);
    let bond_action = segment_action(selected_difference(io, bond_segment)?);
    let center_action = segment_action(selected_difference(io, center_segment)?);
    let inversion_action = segment_action(selected_difference(io, inversion_segment)?);
    let type_action = segment_action(selected_difference(io, type_segment)?);

    if bond_action != 0 || center_action != 0 || inversion_action != 0 || type_action != 0 {
        let tags = stereo_ident_labels(heap)?;
        io.bTag1 = crate::source_types::local_ichiprt1::tagIdentLblBit_IL_STER as i32 | io.bFhTag;
        let (_, tag1) = make_tag(heap, &tags, io.nTag, io.bTag1, &io.szTag1, &mut io.bAlways)?;
        io.szTag1 = tag1;

        io.nSegmAction = bond_action;
        if io.nSegmAction != 0 {
            io.bTag2 =
                io.bTag1 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_DBND as i32;
            let (tag, value) =
                make_tag(heap, &tags, io.nTag, io.bTag2, &io.szTag2, &mut io.bAlways)?;
            io.szTag2 = value;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            if io.nSegmAction == fill {
                io.tot_len = str_Sp2(
                    heap,
                    io.pINChISort,
                    io.pINChISort2,
                    string_buffer,
                    &mut io.bOverflow,
                    io.bOutType,
                    io.TAUT_MODE,
                    io.num_components,
                    io.bSecondNonTautPass,
                    io.bOmitRepetitions,
                    io.bUseMulipliers,
                )?;
                io.bNonTautNonIsoIdentifierNotEmpty = io
                    .bNonTautNonIsoIdentifierNotEmpty
                    .wrapping_add(io.bSecondNonTautPass);
            }
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                io.nSegmAction.wrapping_neg(),
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(1);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
        } else if io.bPlainTextTags == 1 {
            print_literal(heap, output, b"/\0", stdout)?;
        }

        io.nSegmAction = center_action;
        if io.nSegmAction != 0 {
            let taut_index = usize::try_from(io.iCurTautMode)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            io.bRelRac = i32::from(
                *io.bRelativeStereo
                    .get(taut_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
                    || *io
                        .bRacemicStereo
                        .get(taut_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0,
            );
            io.bTag2 =
                io.bTag1 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_SP3S as i32;
            let (tag, value) =
                make_tag(heap, &tags, io.nTag, io.bTag2, &io.szTag2, &mut io.bAlways)?;
            io.szTag2 = value;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            if io.nSegmAction == fill {
                io.tot_len = str_Sp3(
                    heap,
                    io.pINChISort,
                    io.pINChISort2,
                    string_buffer,
                    &mut io.bOverflow,
                    io.bOutType,
                    io.TAUT_MODE,
                    io.num_components,
                    io.bRelRac,
                    io.bSecondNonTautPass,
                    io.bOmitRepetitions,
                    io.bUseMulipliers,
                )?;
                io.bNonTautNonIsoIdentifierNotEmpty = io
                    .bNonTautNonIsoIdentifierNotEmpty
                    .wrapping_add(io.bSecondNonTautPass);
            }
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                io.nSegmAction.wrapping_neg(),
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(2);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
        } else if io.bPlainTextTags == 1 {
            print_literal(heap, output, b"/\0", stdout)?;
        }

        io.nSegmAction = inversion_action;
        if io.nSegmAction != 0 {
            io.bTag2 =
                io.bTag1 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_INVS as i32;
            let (tag, value) =
                make_tag(heap, &tags, io.nTag, io.bTag2, &io.szTag2, &mut io.bAlways)?;
            io.szTag2 = value;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            if io.nSegmAction == fill {
                io.tot_len = str_StereoAbsInv(
                    heap,
                    io.pINChISort,
                    string_buffer,
                    &mut io.bOverflow,
                    io.bOutType,
                    io.num_components,
                )?;
                io.bNonTautNonIsoIdentifierNotEmpty = io
                    .bNonTautNonIsoIdentifierNotEmpty
                    .wrapping_add(io.bSecondNonTautPass);
            }
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                io.nSegmAction.wrapping_neg(),
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(3);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
        } else if io.bPlainTextTags == 1 {
            print_literal(heap, output, b"/\0", stdout)?;
        }

        io.nSegmAction = type_action;
        if io.nSegmAction != 0 {
            let taut_index = usize::try_from(io.iCurTautMode)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let relative = *io
                .bRelativeStereo
                .get(taut_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                != 0;
            let racemic = *io
                .bRacemicStereo
                .get(taut_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                != 0;
            io.bTag2 =
                io.bTag1 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_TYPS as i32;
            let (tag, value) =
                make_tag(heap, &tags, io.nTag, io.bTag2, &io.szTag2, &mut io.bAlways)?;
            io.szTag2 = value;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            if io.nSegmAction == fill {
                let stereo_text = if relative {
                    b"2\0".as_slice()
                } else if racemic {
                    b"3\0".as_slice()
                } else {
                    b"1\0".as_slice()
                };
                let stereo = heap
                    .allocate_model_storage(stereo_text.iter().map(|byte| *byte as i8).collect())?;
                io.tot_len = io.tot_len.wrapping_add(MakeDelim(
                    heap,
                    stereo.as_const(),
                    string_buffer,
                    &mut io.bOverflow,
                )?);
                io.bNonTautNonIsoIdentifierNotEmpty = io
                    .bNonTautNonIsoIdentifierNotEmpty
                    .wrapping_add(io.bSecondNonTautPass);
            }
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                io.nSegmAction.wrapping_neg(),
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(1);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
        }
        if io.bPlainTextTags == 1 {
            print_literal(heap, output, b"/\0", stdout)?;
        }
    } else if io.bPlainTextTags == 1 {
        print_literal(heap, output, b"////\0", stdout)?;
    }

    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_IsotopicLayer(
    heap: &mut SourceHeap,
    canon_globals: SourceMutPointer<CANON_GLOBALS>,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    basic_or_reconnected: &i32,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3502 OutputINCHI_IsotopicLayer
    // INCHI✔️❌: int OutputINCHI_IsotopicLayer( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                int              *INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                char             *pLF,
    // INCHI✔️❌:                                char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     if (INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_i_IATOMS] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  isotopic #1:  composition -- atoms -- do not output in xml if empty */
    // INCHI✔️❌:         szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_ISOT | io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:         /* isotopic atoms without mobile H.
    // INCHI✔️❌:          * Fixed 2004-06-15: always output if not bXml. Note:
    // INCHI✔️❌:          * Previous condition if( bHasIsotopicAtoms[io->iCurTautMode] || bIsotopic && !bXml)
    // INCHI✔️❌:          * did not optput /i in case of only mobile isotopic H
    // INCHI✔️❌:          */
    // INCHI✔️❌:         if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_i_IATOMS] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_ATMS, io->szTag2, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:             io->tot_len = 0;
    // INCHI✔️❌:             /*if ( bHasIsotopicAtoms[io->iCurTautMode] )*/
    // INCHI✔️❌:             if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 io->tot_len2 = str_IsoAtoms( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
    // INCHI✔️❌:                                              io->bOutType, io->TAUT_MODE, io->num_components, io->bAbcNumbers,
    // INCHI✔️❌:                                              io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:                 io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 io->tot_len2 = io->tot_len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             io->tot_len = io->tot_len2;
    // INCHI✔️❌:             if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  isotopic #1a:  composition -- exchangeable isotopic H (mobile H only) */
    // INCHI✔️❌:         /*if ( !io->bSecondNonTautPass && bHasIsoH )*/
    // INCHI✔️❌:         if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_h_H_ATOMS] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_XCGA, io->szTag2, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:             io->tot_len = 0;
    // INCHI✔️❌:             ( io->tot_len ) += MakeIsoHString( io->num_iso_H, strbuf, io->TAUT_MODE, &io->bOverflow );
    // INCHI✔️❌:             io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:             if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /***************************************************
    // INCHI✔️❌:          *
    // INCHI✔️❌:          *       Isotopic stereo
    // INCHI✔️❌:          *
    // INCHI✔️❌:          ***************************************************/
    // INCHI✔️❌:
    // INCHI✔️❌:         /*if ( bIsotopicStereo[io->iCurTautMode] )*/
    // INCHI✔️❌:         if (INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ) ||
    // INCHI✔️❌:              INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ) ||
    // INCHI✔️❌:              INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ) ||
    // INCHI✔️❌:              INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /*  stereo */
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag2 = io->bTag1 | IL_STER, io->szTag2, &io->bAlways, 1 );
    // INCHI✔️❌:
    // INCHI✔️❌:             /************************
    // INCHI✔️❌:               isotopic #2:  sp2
    // INCHI✔️❌:              ************************/
    // INCHI✔️❌:             /*if ( bIsotopicStereoSp2[io->iCurTautMode]  )*/
    // INCHI✔️❌:             if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_b_SBONDS] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_DBND, io->szTag3, &io->bAlways, 1 );
    // INCHI✔️❌:                 inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:                 io->tot_len = 0;
    // INCHI✔️❌:                 if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     io->tot_len = str_IsoSp2( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
    // INCHI✔️❌:                                               io->bOutType, io->TAUT_MODE, io->num_components,
    // INCHI✔️❌:                                           io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:                     io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 3;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (io->bPlainTextTags == 1) inchi_ios_print_nodisplay( out_file, "/" ); /* iso sp2 */
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /************************
    // INCHI✔️❌:               isotopic #3:  sp3
    // INCHI✔️❌:              ************************/
    // INCHI✔️❌:             /*if ( bIsotopicStereoSp3[io->iCurTautMode]  )*/
    // INCHI✔️❌:             if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_t_SATOMS] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 io->bRelRac = io->bIsotopicRelativeStereo[io->iCurTautMode] || io->bIsotopicRacemicStereo[io->iCurTautMode];
    // INCHI✔️❌:
    // INCHI✔️❌:                 szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_SP3S, io->szTag3, &io->bAlways, 1 );
    // INCHI✔️❌:                 inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:                 io->tot_len = 0;
    // INCHI✔️❌:                 if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     io->tot_len = str_IsoSp3( io->pINChISort, io->pINChISort2, strbuf, &io->bOverflow,
    // INCHI✔️❌:                                               io->bOutType, io->TAUT_MODE, io->num_components, io->bRelRac,
    // INCHI✔️❌:                                               io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:                     io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 5;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (io->bPlainTextTags == 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print_nodisplay( out_file, "/" ); /* iso-sp3 */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* isotopic #4: abs inverted */
    // INCHI✔️❌:             if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_m_SP3INV] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_INVS, io->szTag3, &io->bAlways, 1 );
    // INCHI✔️❌:                 inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:                 io->tot_len = 0;
    // INCHI✔️❌:                 if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     io->tot_len = str_IsoStereoAbsInv( io->pINChISort, strbuf,
    // INCHI✔️❌:                                                    &io->bOverflow, io->bOutType, io->num_components );
    // INCHI✔️❌:                     io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 5;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (io->bPlainTextTags == 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     inchi_ios_print_nodisplay( out_file, "/" );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             /* isotopic #5: stereo type. Do not output if it has already been output in non-iso */
    // INCHI✔️❌:             if ((io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_s_STYPE] ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 const char *p_stereo = io->bIsotopicRelativeStereo[io->iCurTautMode] ? x_rel :
    // INCHI✔️❌:                     io->bIsotopicRacemicStereo[io->iCurTautMode] ? x_rac : x_abs;
    // INCHI✔️❌:                 szGetTag( IdentLbl, io->nTag, io->bTag3 = io->bTag2 | IL_TYPS, io->szTag3, &io->bAlways, 1 );
    // INCHI✔️❌:                 inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:                 io->tot_len = 0;
    // INCHI✔️❌:                 if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     io->tot_len += MakeDelim( p_stereo, strbuf, &io->bOverflow );
    // INCHI✔️❌:                     io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     return 6;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (io->bPlainTextTags == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "/" );  /* no abs, inv or racemic stereo */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* no isotopic stereo */
    // INCHI✔️❌:             if (io->bPlainTextTags == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "////" ); /* sp3, sp2, abs-inv, stereo.type */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (io->bPlainTextTags == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "///" ); /* isotopic composition, sp2, sp3 */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (io->bPlainTextTags == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "//" );   /* inv or racemic stereo */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: #if ( CANON_FIXH_TRANS == 1 )
    // INCHI✔️❌:     if (io->bOutType == OUT_NONTAUT && io->bOutputType == OUT_TN && io->bSecondNonTautPass &&
    // INCHI✔️❌:          INCHI_SEGM_FILL == INChI_SegmentAction( io->sDifSegs[DIFL_F][DIFS_o_TRANSP] ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* find and print non-tautomeric components transposition, if non-trivial */
    // INCHI✔️❌:         AT_NUMB *nTrans_n, *nTrans_s;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (0 < bin_AuxTautTrans( io->pINChISort, io->pINChISort2, &nTrans_n, &nTrans_s, io->bOutType, io->num_components ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* a non-trivial transposition does exist; output start tag */
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_TRNS | io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:             io->tot_len = 0;
    // INCHI✔️❌:             /* print the transposition, cycle after cycle */
    // INCHI✔️❌:             io->tot_len = str_AuxTautTrans( pCG, nTrans_n, nTrans_s, strbuf,
    // INCHI✔️❌:                                             &io->bOverflow, io->TAUT_MODE, io->num_components );
    // INCHI✔️❌:             io->bNonTautIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:             if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 7;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:              /* detected transposition */
    // INCHI✔️❌:             ( *io->pSortPrintINChIFlags ) |=
    // INCHI✔️❌:                 ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_TRANSPOS_BAS : FLAG_SORT_PRINT_TRANSPOS_REC;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io->bPlainTextTags == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print_nodisplay( out_file, "/" );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_IsotopicLayer

    fn selected_difference(
        io: &INCHI_OUT_CTL,
        layer: usize,
        segment: usize,
    ) -> Result<i8, SourceHeapError> {
        io.sDifSegs
            .get(layer)
            .and_then(|segments| segments.get(segment))
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }
    fn make_tag(
        heap: &mut SourceHeap,
        tags: &[INCHI_TAG],
        n_tag: i32,
        bits: i32,
        source: &[i8; 64],
        always: &mut i32,
    ) -> Result<(SourceMutPointer<i8>, [i8; 64]), SourceHeapError> {
        let tag = heap.allocate_model_storage(source.to_vec())?;
        szGetTag(heap, tags, n_tag, bits, tag, always, 1)?;
        let mut value = [0_i8; 64];
        value.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        Ok((tag, value))
    }
    fn print_buffer(
        heap: &mut SourceHeap,
        output: &mut INCHI_IOSTREAM,
        buffer: &INCHI_IOS_STRING,
        line_feed: SourceConstPointer<i8>,
        stdout: SourceMutPointer<FILE>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }
    fn print_literal(
        heap: &mut SourceHeap,
        output: &mut INCHI_IOSTREAM,
        literal: &[u8],
        stdout: SourceMutPointer<FILE>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(literal.iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList::default(),
        )?;
        Ok(())
    }

    let current_layer =
        usize::try_from(io.nCurINChISegment).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let action = |value| crate::source::base::ichimake::INChI_SegmentAction(value);
    let fill = crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32;
    let atom_segment = crate::source_types::tagDiffINChISegments_DIFS_i_IATOMS as usize;
    let hydrogen_segment = crate::source_types::tagDiffINChISegments_DIFS_h_H_ATOMS as usize;
    let bond_segment = crate::source_types::tagDiffINChISegments_DIFS_b_SBONDS as usize;
    let center_segment = crate::source_types::tagDiffINChISegments_DIFS_t_SATOMS as usize;
    let inversion_segment = crate::source_types::tagDiffINChISegments_DIFS_m_SP3INV as usize;
    let type_segment = crate::source_types::tagDiffINChISegments_DIFS_s_STYPE as usize;
    let atom_action = action(selected_difference(io, current_layer, atom_segment)?);
    let tags = isotopic_ident_labels(heap)?;

    if atom_action != 0 {
        io.bTag1 = crate::source_types::local_ichiprt1::tagIdentLblBit_IL_ISOT as i32 | io.bFhTag;
        let (_, value) = make_tag(heap, &tags, io.nTag, io.bTag1, &io.szTag1, &mut io.bAlways)?;
        io.szTag1 = value;

        io.nSegmAction = atom_action;
        io.bTag2 = io.bTag1 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_ATMS as i32;
        let (tag, value) = make_tag(heap, &tags, io.nTag, io.bTag2, &io.szTag2, &mut io.bAlways)?;
        io.szTag2 = value;
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        if io.nSegmAction == fill {
            io.tot_len2 = str_IsoAtoms(
                heap,
                io.pINChISort,
                io.pINChISort2,
                string_buffer,
                &mut io.bOverflow,
                io.bOutType,
                io.TAUT_MODE,
                io.num_components,
                io.bAbcNumbers,
                io.bSecondNonTautPass,
                io.bOmitRepetitions,
                io.bUseMulipliers,
            )?;
            io.bNonTautIsoIdentifierNotEmpty = io
                .bNonTautIsoIdentifierNotEmpty
                .wrapping_add(io.bSecondNonTautPass);
        } else {
            io.tot_len2 = io.tot_len;
        }
        io.tot_len = io.tot_len2;
        if str_LineEnd(
            heap,
            tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            io.nSegmAction.wrapping_neg(),
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_buffer(heap, output, string_buffer, line_feed, stdout)?;

        io.nSegmAction = action(selected_difference(io, current_layer, hydrogen_segment)?);
        if io.nSegmAction != 0 {
            io.bTag2 =
                io.bTag1 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_XCGA as i32;
            let (tag, value) =
                make_tag(heap, &tags, io.nTag, io.bTag2, &io.szTag2, &mut io.bAlways)?;
            io.szTag2 = value;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            io.tot_len = io.tot_len.wrapping_add(MakeIsoHString(
                heap,
                &io.num_iso_H,
                string_buffer,
                io.TAUT_MODE,
                &mut io.bOverflow,
            )?);
            io.bNonTautIsoIdentifierNotEmpty = io
                .bNonTautIsoIdentifierNotEmpty
                .wrapping_add(io.bSecondNonTautPass);
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                io.nSegmAction.wrapping_neg(),
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(2);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
        }

        let bond_action = action(selected_difference(io, current_layer, bond_segment)?);
        let center_action = action(selected_difference(io, current_layer, center_segment)?);
        let inversion_action = action(selected_difference(io, current_layer, inversion_segment)?);
        let type_action = action(selected_difference(io, current_layer, type_segment)?);
        if bond_action != 0 || center_action != 0 || inversion_action != 0 || type_action != 0 {
            io.bTag2 =
                io.bTag1 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_STER as i32;
            let (_, value) = make_tag(heap, &tags, io.nTag, io.bTag2, &io.szTag2, &mut io.bAlways)?;
            io.szTag2 = value;

            io.nSegmAction = bond_action;
            if io.nSegmAction != 0 {
                io.bTag3 =
                    io.bTag2 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_DBND as i32;
                let (tag, value) =
                    make_tag(heap, &tags, io.nTag, io.bTag3, &io.szTag3, &mut io.bAlways)?;
                io.szTag3 = value;
                inchi_strbuf_reset(heap, Some(string_buffer))?;
                io.tot_len = 0;
                if io.nSegmAction == fill {
                    io.tot_len = str_IsoSp2(
                        heap,
                        io.pINChISort,
                        io.pINChISort2,
                        string_buffer,
                        &mut io.bOverflow,
                        io.bOutType,
                        io.TAUT_MODE,
                        io.num_components,
                        io.bSecondNonTautPass,
                        io.bOmitRepetitions,
                        io.bUseMulipliers,
                    )?;
                    io.bNonTautIsoIdentifierNotEmpty = io
                        .bNonTautIsoIdentifierNotEmpty
                        .wrapping_add(io.bSecondNonTautPass);
                }
                if str_LineEnd(
                    heap,
                    tag.as_const(),
                    &mut io.bOverflow,
                    string_buffer,
                    io.nSegmAction.wrapping_neg(),
                    io.bPlainTextTags,
                )? != 0
                {
                    return Ok(3);
                }
                print_buffer(heap, output, string_buffer, line_feed, stdout)?;
            } else if io.bPlainTextTags == 1 {
                print_literal(heap, output, b"/\0", stdout)?;
            }

            io.nSegmAction = center_action;
            if io.nSegmAction != 0 {
                let taut = usize::try_from(io.iCurTautMode)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                io.bRelRac = i32::from(
                    *io.bIsotopicRelativeStereo
                        .get(taut)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        != 0
                        || *io
                            .bIsotopicRacemicStereo
                            .get(taut)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            != 0,
                );
                io.bTag3 =
                    io.bTag2 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_SP3S as i32;
                let (tag, value) =
                    make_tag(heap, &tags, io.nTag, io.bTag3, &io.szTag3, &mut io.bAlways)?;
                io.szTag3 = value;
                inchi_strbuf_reset(heap, Some(string_buffer))?;
                io.tot_len = 0;
                if io.nSegmAction == fill {
                    io.tot_len = str_IsoSp3(
                        heap,
                        io.pINChISort,
                        io.pINChISort2,
                        string_buffer,
                        &mut io.bOverflow,
                        io.bOutType,
                        io.TAUT_MODE,
                        io.num_components,
                        io.bRelRac,
                        io.bSecondNonTautPass,
                        io.bOmitRepetitions,
                        io.bUseMulipliers,
                    )?;
                    io.bNonTautIsoIdentifierNotEmpty = io
                        .bNonTautIsoIdentifierNotEmpty
                        .wrapping_add(io.bSecondNonTautPass);
                }
                if str_LineEnd(
                    heap,
                    tag.as_const(),
                    &mut io.bOverflow,
                    string_buffer,
                    io.nSegmAction.wrapping_neg(),
                    io.bPlainTextTags,
                )? != 0
                {
                    return Ok(5);
                }
                print_buffer(heap, output, string_buffer, line_feed, stdout)?;
            } else if io.bPlainTextTags == 1 {
                print_literal(heap, output, b"/\0", stdout)?;
            }

            io.nSegmAction = inversion_action;
            if io.nSegmAction != 0 {
                io.bTag3 =
                    io.bTag2 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_INVS as i32;
                let (tag, value) =
                    make_tag(heap, &tags, io.nTag, io.bTag3, &io.szTag3, &mut io.bAlways)?;
                io.szTag3 = value;
                inchi_strbuf_reset(heap, Some(string_buffer))?;
                io.tot_len = 0;
                if io.nSegmAction == fill {
                    io.tot_len = str_IsoStereoAbsInv(
                        heap,
                        io.pINChISort,
                        string_buffer,
                        &mut io.bOverflow,
                        io.bOutType,
                        io.num_components,
                    )?;
                    io.bNonTautIsoIdentifierNotEmpty = io
                        .bNonTautIsoIdentifierNotEmpty
                        .wrapping_add(io.bSecondNonTautPass);
                }
                if str_LineEnd(
                    heap,
                    tag.as_const(),
                    &mut io.bOverflow,
                    string_buffer,
                    io.nSegmAction.wrapping_neg(),
                    io.bPlainTextTags,
                )? != 0
                {
                    return Ok(5);
                }
                print_buffer(heap, output, string_buffer, line_feed, stdout)?;
            } else if io.bPlainTextTags == 1 {
                print_literal(heap, output, b"/\0", stdout)?;
            }

            io.nSegmAction = type_action;
            if io.nSegmAction != 0 {
                let taut = usize::try_from(io.iCurTautMode)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let stereo_text = if *io
                    .bIsotopicRelativeStereo
                    .get(taut)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
                {
                    b"2\0".as_slice()
                } else if *io
                    .bIsotopicRacemicStereo
                    .get(taut)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    != 0
                {
                    b"3\0".as_slice()
                } else {
                    b"1\0".as_slice()
                };
                io.bTag3 =
                    io.bTag2 | crate::source_types::local_ichiprt1::tagIdentLblBit_IL_TYPS as i32;
                let (tag, value) =
                    make_tag(heap, &tags, io.nTag, io.bTag3, &io.szTag3, &mut io.bAlways)?;
                io.szTag3 = value;
                inchi_strbuf_reset(heap, Some(string_buffer))?;
                io.tot_len = 0;
                if io.nSegmAction == fill {
                    let stereo = heap.allocate_model_storage(
                        stereo_text.iter().map(|byte| *byte as i8).collect(),
                    )?;
                    io.tot_len = io.tot_len.wrapping_add(MakeDelim(
                        heap,
                        stereo.as_const(),
                        string_buffer,
                        &mut io.bOverflow,
                    )?);
                    io.bNonTautIsoIdentifierNotEmpty = io
                        .bNonTautIsoIdentifierNotEmpty
                        .wrapping_add(io.bSecondNonTautPass);
                }
                if str_LineEnd(
                    heap,
                    tag.as_const(),
                    &mut io.bOverflow,
                    string_buffer,
                    io.nSegmAction.wrapping_neg(),
                    io.bPlainTextTags,
                )? != 0
                {
                    return Ok(6);
                }
                print_buffer(heap, output, string_buffer, line_feed, stdout)?;
            }
            if io.bPlainTextTags == 1 {
                print_literal(heap, output, b"/\0", stdout)?;
            }
        } else if io.bPlainTextTags == 1 {
            print_literal(heap, output, b"////\0", stdout)?;
        }
    } else if io.bPlainTextTags == 1 {
        print_literal(heap, output, b"///\0", stdout)?;
        print_literal(heap, output, b"//\0", stdout)?;
    }

    let transposition_layer = crate::source_types::tagDiffINChILayers_DIFL_F as usize;
    let transposition_segment = crate::source_types::tagDiffINChISegments_DIFS_o_TRANSP as usize;
    if io.bOutType == crate::source_types::local_ichiprt1::OUT_NONTAUT as i32
        && io.bOutputType == crate::source_types::OUT_TN as i32
        && io.bSecondNonTautPass != 0
        && action(selected_difference(
            io,
            transposition_layer,
            transposition_segment,
        )?) == fill
    {
        let mut non_taut = SourceMutPointer::null();
        let mut taut = SourceMutPointer::null();
        if bin_AuxTautTrans(
            heap,
            io.pINChISort,
            io.pINChISort2,
            &mut non_taut,
            &mut taut,
            io.bOutType,
            io.num_components,
        )? > 0
        {
            io.bTag1 =
                crate::source_types::local_ichiprt1::tagIdentLblBit_IL_TRNS as i32 | io.bFhTag;
            let (tag, value) =
                make_tag(heap, &tags, io.nTag, io.bTag1, &io.szTag1, &mut io.bAlways)?;
            io.szTag1 = value;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            io.tot_len = str_AuxTautTrans(
                heap,
                non_taut,
                taut,
                string_buffer,
                &mut io.bOverflow,
                io.TAUT_MODE,
                io.num_components,
            )?;
            io.bNonTautIsoIdentifierNotEmpty = io
                .bNonTautIsoIdentifierNotEmpty
                .wrapping_add(io.bSecondNonTautPass);
            if str_LineEnd(
                heap,
                tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                -1,
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(7);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
            let flags = heap
                .slice_mut(io.pSortPrintINChIFlags)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *flags |= if *basic_or_reconnected == crate::source_types::INCHI_BAS as i32 {
                crate::source_types::FLAG_SORT_PRINT_TRANSPOS_BAS as i32
            } else {
                crate::source_types::FLAG_SORT_PRINT_TRANSPOS_REC as i32
            };
        } else if io.bPlainTextTags == 1 {
            print_literal(heap, output, b"/\0", stdout)?;
        }
    }
    let _ = canon_globals;
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_FixedHLayerWithSublayers(
    heap: &mut SourceHeap,
    _canon_globals: SourceMutPointer<CANON_GLOBALS>,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    basic_or_reconnected: &i32,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    then_goto_repeat: &mut i32,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3750 OutputINCHI_FixedHLayerWithSublayers
    // INCHI✔️❌: int OutputINCHI_FixedHLayerWithSublayers( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                           INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                           INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                           int              *INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                           INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                           char             *pLF,
    // INCHI✔️❌:                                           char             *pTAB,
    // INCHI✔️❌:                                           int              *then_goto_repeat )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     *then_goto_repeat = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (io->bOutType == OUT_TN &&
    // INCHI✔️❌:          !( io->bSecondNonTautPass ) &&
    // INCHI✔️❌:          io->bNonTautIsIdenticalToTaut &&
    // INCHI✔️❌:          io->bTautomeric &&
    // INCHI✔️❌:          io->bNonTautomeric)
    // INCHI✔️❌:     {
    // INCHI✔️❌:             /* Fixed-H layer is empty in the Identifier */
    // INCHI✔️❌:         ( *io->pSortPrintINChIFlags ) |=
    // INCHI✔️❌:             ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
    // INCHI✔️❌:             FLAG_SORT_PRINT_NO_NFIX_H_REC;
    // INCHI✔️❌:         ( *io->pSortPrintINChIFlags ) |=
    // INCHI✔️❌:             ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
    // INCHI✔️❌:             FLAG_SORT_PRINT_NO_IFIX_H_REC;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (io->bOutType == OUT_TN &&
    // INCHI✔️❌:          !io->bNonTautIsIdenticalToTaut && /* added 2004-10-04 Fix16 */
    // INCHI✔️❌: #ifdef OLD_ITEM_DISCOVERY
    // INCHI✔️❌:          io->bTautomeric &&
    // INCHI✔️❌:          io->bNonTautomeric &&
    // INCHI✔️❌: #endif
    // INCHI✔️❌:          INChI_SegmentAction( io->sDifSegs[DIFL_F][DIFS_f_FORMULA] )
    // INCHI✔️❌:                                     /* special case: removed isolated H(+): */
    // INCHI✔️❌:                                     /* || io->iCurTautMode == TAUT_YES && num_comp[TAUT_YES] < num_comp[TAUT_NON] &&
    // INCHI✔️❌:                                         0 < num_comp[TAUT_NON]*/
    // INCHI✔️❌:        )
    // INCHI✔️❌:
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* add the second (non-tautomeric) output */
    // INCHI✔️❌:         io->bOutType = OUT_NONTAUT;    /* pick up only non-tautomeric representation of tautomeric */
    // INCHI✔️❌:         io->iCurTautMode = TAUT_NON;
    // INCHI✔️❌:         io->pINChISort = io->pINChISortTautAndNonTaut[TAUT_NON];
    // INCHI✔️❌:         io->bSecondNonTautPass = 1;
    // INCHI✔️❌:         io->nCurINChISegment = DIFL_F;
    // INCHI✔️❌:         io->num_components = io->num_comp[io->iCurTautMode]; /* number of components could change due to removal of isolated H(+) from tautomeric */
    // INCHI✔️❌:         io->bFhTag = IL_FIXH;
    // INCHI✔️❌:         szGetTag( IdentLbl, io->nTag, io->bTag1 = io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:         /***** constitution non-taut: dot-disconnected Hill formulas: <formula> -- only if different */
    // INCHI✔️❌:         szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_FMLF | io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:         io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_f_FORMULA] );
    // INCHI✔️❌:         if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             io->tot_len2 = str_HillFormula2( io->pINChISort, io->pINChISort2,
    // INCHI✔️❌:                                              strbuf, &io->bOverflow, io->bOutType,
    // INCHI✔️❌:                                              io->num_components, io->bUseMulipliers );
    // INCHI✔️❌:             if (io->n_pzz > 0 && io->n_zy > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 MergeZzInHillFormula(strbuf);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             io->tot_len2 = io->tot_len;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         io->tot_len = io->tot_len2;
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:
    // INCHI✔️❌:         io->nSegmAction = INChI_SegmentAction( io->sDifSegs[io->nCurINChISegment][DIFS_h_H_ATOMS] );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (INCHI_SEGM_FILL == io->nSegmAction)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( IdentLbl, io->nTag, io->bTag1 = IL_HFIX | io->bFhTag, io->szTag1, &io->bAlways, 1 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf ); io->tot_len = 0; /* open H-fixed */
    // INCHI✔️❌:             /* output the second non-tautomeric item: fixed H -- do not output in xml if empty */
    // INCHI✔️❌:             io->tot_len2 = str_FixedH_atoms( io->pINChISort, strbuf,
    // INCHI✔️❌:                                              &io->bOverflow, io->bOutType, io->ATOM_MODE,
    // INCHI✔️❌:                                              io->num_components, io->bUseMulipliers );
    // INCHI✔️❌:             io->tot_len = io->tot_len2;
    // INCHI✔️❌:             if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -io->nSegmAction, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 2;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:             io->bNonTautNonIsoIdentifierNotEmpty += io->bSecondNonTautPass;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         *then_goto_repeat = 1;
    // INCHI✔️❌:         return 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (io->bOutType == OUT_NONTAUT && io->bOutputType == OUT_TN && io->bSecondNonTautPass /* && io->bTautomeric && io->bNonTautomeric*/)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* the second (non-taut) output has been done; restore variables */
    // INCHI✔️❌:             io->bOutType = OUT_TN;
    // INCHI✔️❌:             io->iCurTautMode = TAUT_YES;
    // INCHI✔️❌:             io->pINChISort = io->pINChISortTautAndNonTaut[TAUT_YES];
    // INCHI✔️❌:             io->bSecondNonTautPass = 0;
    // INCHI✔️❌:             io->num_components = io->num_comp[io->iCurTautMode];
    // INCHI✔️❌:             if (!io->bNonTautNonIsoIdentifierNotEmpty)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Fixed-H layer is empty in the Identifier */
    // INCHI✔️❌:                 ( *io->pSortPrintINChIFlags ) |= ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_NFIX_H_BAS :
    // INCHI✔️❌:                     FLAG_SORT_PRINT_NO_NFIX_H_REC;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (!io->bNonTautIsoIdentifierNotEmpty)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Fixed-H layer is empty in the Identifier */
    // INCHI✔️❌:                 ( *io->pSortPrintINChIFlags ) |= ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
    // INCHI✔️❌:                     FLAG_SORT_PRINT_NO_IFIX_H_REC;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             io->bFhTag = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_FixedHLayerWithSublayers

    fn selected_sort(
        heap: &SourceHeap,
        io: &INCHI_OUT_CTL,
        taut_mode: usize,
    ) -> Result<SourceMutPointer<crate::source_types::INCHI_SORT>, SourceHeapError> {
        heap.slice(io.pINChISortTautAndNonTaut.as_const())?
            .get(taut_mode)
            .copied()
            .ok_or(SourceHeapError::PointerOutOfBounds)
    }

    fn update_tag(
        heap: &mut SourceHeap,
        tags: &[INCHI_TAG],
        io: &mut INCHI_OUT_CTL,
        tag_bits: i32,
    ) -> Result<SourceMutPointer<i8>, SourceHeapError> {
        io.bTag1 = tag_bits;
        let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
        szGetTag(heap, tags, io.nTag, io.bTag1, tag, &mut io.bAlways, 1)?;
        let tag_length = io.szTag1.len();
        io.szTag1.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..tag_length)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        Ok(tag)
    }

    fn print_buffer(
        heap: &mut SourceHeap,
        output: &mut INCHI_IOSTREAM,
        string_buffer: &INCHI_IOS_STRING,
        line_feed: SourceConstPointer<i8>,
        stdout: SourceMutPointer<FILE>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print_nodisplay(
            heap,
            Some(output),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    *then_goto_repeat = 0;
    let out_tn = crate::source_types::OUT_TN as i32;
    let out_non_taut = crate::source_types::local_ichiprt1::OUT_NONTAUT as i32;
    let taut_non = crate::source_types::TAUT_NON as usize;
    let taut_yes = crate::source_types::TAUT_YES as usize;
    let basic = *basic_or_reconnected == crate::source_types::INCHI_BAS as i32;

    if io.bOutType == out_tn
        && io.bSecondNonTautPass == 0
        && io.bNonTautIsIdenticalToTaut != 0
        && io.bTautomeric != 0
        && io.bNonTautomeric != 0
    {
        let flags = heap
            .slice_mut(io.pSortPrintINChIFlags)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        *flags |= if basic {
            crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_BAS as i32
        } else {
            crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_REC as i32
        };
        *flags |= if basic {
            crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_BAS as i32
        } else {
            crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_REC as i32
        };
    }

    let fixed_layer = crate::source_types::tagDiffINChILayers_DIFL_F as usize;
    let formula_segment = crate::source_types::tagDiffINChISegments_DIFS_f_FORMULA as usize;
    let fixed_h_segment = crate::source_types::tagDiffINChISegments_DIFS_h_H_ATOMS as usize;
    let formula_difference = io.sDifSegs[fixed_layer][formula_segment];
    if io.bOutType == out_tn
        && io.bNonTautIsIdenticalToTaut == 0
        && crate::source::base::ichimake::INChI_SegmentAction(formula_difference) != 0
    {
        io.bOutType = out_non_taut;
        io.iCurTautMode = crate::source_types::TAUT_NON as i32;
        io.pINChISort = selected_sort(heap, io, taut_non)?;
        io.bSecondNonTautPass = 1;
        io.nCurINChISegment = crate::source_types::tagDiffINChILayers_DIFL_F as i32;
        io.num_components = io.num_comp[taut_non];
        io.bFhTag = crate::source_types::local_ichiprt1::tagIdentLblBit_IL_FIXH as i32;

        let tags = isotopic_ident_labels(heap)?;
        update_tag(heap, &tags, io, io.bFhTag)?;
        let formula_tag = update_tag(
            heap,
            &tags,
            io,
            crate::source_types::local_ichiprt1::tagIdentLblBit_IL_FMLF as i32 | io.bFhTag,
        )?;
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.nSegmAction = crate::source::base::ichimake::INChI_SegmentAction(
            io.sDifSegs[fixed_layer][formula_segment],
        );
        if io.nSegmAction == crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32 {
            io.tot_len2 = str_HillFormula2(
                heap,
                io.pINChISort,
                io.pINChISort2,
                string_buffer,
                &mut io.bOverflow,
                io.bOutType,
                io.num_components,
                io.bUseMulipliers,
            )?;
            if io.n_pzz > 0 && io.n_zy > 0 {
                let _ = MergeZzInHillFormula(heap, string_buffer)?;
            }
            io.bNonTautNonIsoIdentifierNotEmpty = io
                .bNonTautNonIsoIdentifierNotEmpty
                .wrapping_add(io.bSecondNonTautPass);
        } else {
            io.tot_len2 = io.tot_len;
        }
        io.tot_len = io.tot_len2;
        if str_LineEnd(
            heap,
            formula_tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            io.nSegmAction.wrapping_neg(),
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_buffer(heap, output, string_buffer, line_feed, stdout)?;

        io.nSegmAction = crate::source::base::ichimake::INChI_SegmentAction(
            io.sDifSegs[fixed_layer][fixed_h_segment],
        );
        if io.nSegmAction == crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32 {
            let fixed_h_tag = update_tag(
                heap,
                &tags,
                io,
                crate::source_types::local_ichiprt1::tagIdentLblBit_IL_HFIX as i32 | io.bFhTag,
            )?;
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            io.tot_len2 = str_FixedH_atoms(
                heap,
                io.pINChISort,
                string_buffer,
                &mut io.bOverflow,
                io.bOutType,
                io.ATOM_MODE,
                io.num_components,
                io.bUseMulipliers,
            )?;
            io.tot_len = io.tot_len2;
            if str_LineEnd(
                heap,
                fixed_h_tag.as_const(),
                &mut io.bOverflow,
                string_buffer,
                io.nSegmAction.wrapping_neg(),
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(2);
            }
            print_buffer(heap, output, string_buffer, line_feed, stdout)?;
            io.bNonTautNonIsoIdentifierNotEmpty = io
                .bNonTautNonIsoIdentifierNotEmpty
                .wrapping_add(io.bSecondNonTautPass);
        }
        *then_goto_repeat = 1;
        return Ok(0);
    } else if io.bOutType == out_non_taut && io.bOutputType == out_tn && io.bSecondNonTautPass != 0
    {
        io.bOutType = out_tn;
        io.iCurTautMode = crate::source_types::TAUT_YES as i32;
        io.pINChISort = selected_sort(heap, io, taut_yes)?;
        io.bSecondNonTautPass = 0;
        io.num_components = io.num_comp[taut_yes];
        let flags = heap
            .slice_mut(io.pSortPrintINChIFlags)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if io.bNonTautNonIsoIdentifierNotEmpty == 0 {
            *flags |= if basic {
                crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_BAS as i32
            } else {
                crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_REC as i32
            };
        }
        if io.bNonTautIsoIdentifierNotEmpty == 0 {
            *flags |= if basic {
                crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_BAS as i32
            } else {
                crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_REC as i32
            };
        }
        io.bFhTag = 0;
    }

    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_PolymerLayer(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    _inchi_basic_or_reconnected: &mut i32,
    original_atom_data: SourceConstPointer<ORIG_ATOM_DATA>,
    original_structure: SourceConstPointer<ORIG_STRUCT>,
    output_control: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3880 OutputINCHI_PolymerLayer
    // INCHI✔️❌: complete verbatim source frame follows; the second axis remains open because
    // INCHI✔️❌: SourceHeap allocation lookup and temporary Rust clones add material overhead.
    /*
    static int OutputINCHI_PolymerLayer( CANON_GLOBALS *pCG,
                                         INCHI_IOSTREAM *out_file,
                                         INCHI_IOS_STRING *strbuf,
                                         int *INCHI_basic_or_INCHI_reconnected,
                                         ORIG_ATOM_DATA *orig_inp_data,
                                         ORIG_STRUCT *pOrigStruct,
                                         INCHI_OUT_CTL *io,
                                         char *pLF,
                                         char *pTAB )
    {
        int i, err = 0;
        int nunits2 = 0;
        int n_used_stars = 0;
        int *cano_nums = NULL, *compnt_nums = NULL, *unum = NULL, *old_stars = NULL;
        OAD_PolymerUnit *u = NULL;
        OAD_PolymerUnit **units2 = NULL;
        OAD_Polymer *p = NULL;
        OAD_AtProps *aprops = NULL;
        int nat,num_inp_bonds;
        inp_ATOM    *at = NULL;
        int is_inchi2inchi = 0;

        if (!orig_inp_data)
        {
            goto exit_function;
        }

        at = orig_inp_data->at;
        nat = orig_inp_data->num_inp_atoms;
        num_inp_bonds = orig_inp_data->num_inp_bonds;

        if (pOrigStruct && !pOrigStruct->polymer)
        {
            return 0;
        }

        if (pOrigStruct)
        {
            p = pOrigStruct->polymer;
            is_inchi2inchi = !pOrigStruct->szAtoms && !pOrigStruct->szBonds && !pOrigStruct->szCoord;

            if (is_inchi2inchi)
            {
                err = NOT_YET_I2I_FOR_POLYMERS;
                goto exit_function;
            }

            /*OAD_Polymer_DebugTrace( p );*/

            /* Get canonical numbers and numbers-of-components for each original atom */
            cano_nums = (int*)inchi_calloc((long long)pOrigStruct->num_atoms + 1, sizeof(int)); /* djb-rwth: cast operator added */
            if (!cano_nums)
            {
                err = 1;
                goto exit_function;
            }
            compnt_nums = (int*)inchi_calloc((long long)pOrigStruct->num_atoms + 1, sizeof(int)); /* djb-rwth: cast operator added */
            if (!compnt_nums)
            {
                err = 2;
                goto exit_function;
            }
            err = InternallyGetCanoNumsAndComponentNums(pCG,
                strbuf,
                io,
                pOrigStruct->num_atoms,
                cano_nums,
                compnt_nums);
            if (err != 0)
            {
                err = 3;
                goto exit_function;
            }

            /* Set atom properties for sorting */
            aprops = (OAD_AtProps*)inchi_calloc((long long)nat + 1, sizeof(OAD_AtProps)); /* djb-rwth: cast operator added */
            /* nat + 1: add extra element for possibe 1-based indexing */
            if (!aprops)
            {
                /* djb-rwth: avoiding memory leak */
                if (cano_nums)
                {
                    inchi_free(cano_nums);
                }
                if (compnt_nums)
                {
                    inchi_free(compnt_nums);
                }
                if (aprops)
                {
                    inchi_free(aprops);
                }
                return 0;
            }

            /* Note that aprops[] is in orig_atoms domain (0-based) and      */
            /* u (from units) are in cano_nums domain (1-based)              */
            /* Supply non-NULL cano_nums to adjust the domains (base will be adjusted at place) */
            OAD_Polymer_SetAtProps(p, at, nat, &num_inp_bonds, aprops, cano_nums);

            /* Make a working copy of polymer units data: units2 is a copy        */
            /* of original polymer units (p->units) with atomic numbers changed    */
            /* to curr canonical ones; atoms in alists sorted; atoms in blists    */
            /* and blists themselves sorted                                     */
            units2 = (OAD_PolymerUnit**)inchi_calloc(p->n, sizeof(OAD_PolymerUnit*));

            if (NULL == units2)
            {
                err = 3;
                goto exit_function;
            }
            memset(units2, 0, sizeof(*units2)); /* djb-rwth: memset_s C11/Annex K variant? */

            old_stars = (int*)inchi_calloc(pOrigStruct->polymer->n_pzz, sizeof(int));
            if (NULL == old_stars)
            {
                err = 3;
                goto exit_function;
            }
            for (i = 0; i < pOrigStruct->polymer->n_pzz; i++)
            {
                old_stars[i] = pOrigStruct->polymer->pzz[i];
            }

            for (i = 0; i < p->n; i++)
            {
                units2[i] = OAD_PolymerUnit_CreateCopy(p->units[i]);
                if (NULL == units2[i]) /* djb-rwth: unresolved issue -- revision required? -- units2 properly allocated, and loop index well defined */
                {
                    err = 4;
                    goto exit_function;
                }
                nunits2 = i + 1;
            }

            /* unum contains numbers of units (0..p->n) as they go  */
            /* when sorted by alist's in lexicographic order        */
            unum = (int*)inchi_calloc(p->n, sizeof(int));
            if (NULL == unum)
            {
                err = 4;
                goto exit_function;
            }

            err = OAD_Polymer_PrepareWorkingSet(p, cano_nums, compnt_nums, units2, unum);

            if (err != 0)
            {
                err = 5;
                goto exit_function;
            }

            /* Prepare polymer substring */

            /* Mark layer beginning */
            inchi_strbuf_printf(strbuf, "%s", "/z");

            /* Print polymer units data */
            n_used_stars = 0;
            for (i = 0; i < p->n; i++)
            {
                /* For each unit u ... */
                u = units2[unum[i]];
                /* djb-rwth: addressing coverity ID #499574 -- all NULL checks already done above */
                err = OutputINCHI_PolymerLayer_SingleUnit(u,
                    io->bPolymers,
                    pOrigStruct->polymer->n_pzz,
                    &n_used_stars, aprops,
                    cano_nums,
                    orig_inp_data,
                    pOrigStruct, strbuf);
                if (err)
                {
                    goto exit_function;
                }
                if (i < p->n - 1)
                {
                    inchi_strbuf_printf(strbuf, ";");
                }
            }
            inchi_ios_print_nodisplay(out_file, "%s%s", strbuf->pStr, pLF);

            LOG_NO_ARGS("\n******************* (L4184:ichiprt1.c) ********************\n");
            LOG_MULT_ARGS("Polymer Layer start: %s\n", strbuf->pStr);
            LOG_NO_ARGS("\n***********************************************************\n");

        exit_function:
            if (cano_nums)
            {
                inchi_free(cano_nums);
            }
            if (compnt_nums)
            {
                inchi_free(compnt_nums);
            }
            if (aprops)
            {
                inchi_free(aprops);
            }
            if (unum)
            {
                inchi_free(unum);
            }
            if (units2)
            {
                for (i = 0; i < nunits2; i++)
                {
                    OAD_PolymerUnit_Free(units2[i]);
                }
                inchi_free(units2);
            }
            if (old_stars)
            {
                for (i = 0; i < pOrigStruct->polymer->n_pzz; i++)
                    pOrigStruct->polymer->pzz[i] = old_stars[i];
                inchi_free(old_stars);
            }

        }
        return err;
    }
        */
    // END INCHI C FUNCTION: OutputINCHI_PolymerLayer

    if original_atom_data.is_null() {
        return Ok(0);
    }
    let original_data = heap
        .slice(original_atom_data)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let mut number_of_input_bonds = original_data.num_inp_bonds;

    if original_structure.is_null() {
        return Ok(0);
    }
    let structure = heap
        .slice(original_structure)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if structure.polymer.is_null() {
        return Ok(0);
    }
    if structure.szAtoms.is_null() && structure.szBonds.is_null() && structure.szCoord.is_null() {
        return Ok(crate::source_types::local_ichiprt1::NOT_YET_I2I_FOR_POLYMERS as i32);
    }
    let mut polymer = heap
        .slice(structure.polymer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();

    let mut error = 0_i32;
    let mut copied_units = 0_i32;
    let mut canonical_numbers = SourceMutPointer::<i32>::null();
    let mut component_numbers = SourceMutPointer::<i32>::null();
    let mut unit_numbers = SourceMutPointer::<i32>::null();
    let mut old_stars = SourceMutPointer::<i32>::null();
    let mut units2 =
        SourceMutPointer::<SourceMutPointer<crate::source_types::OAD_PolymerUnit>>::null();
    let mut atom_properties = SourceMutPointer::<crate::source_types::OAD_AtProps>::null();

    let body_result = (|| -> Result<(), SourceHeapError> {
        let atom_count = u64::try_from(i64::from(structure.num_atoms).wrapping_add(1))
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        canonical_numbers = match inchi_calloc::<i32>(heap, atom_count, 4) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                error = 1;
                return Ok(());
            }
            Err(other) => return Err(other),
        };
        component_numbers = match inchi_calloc::<i32>(heap, atom_count, 4) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                error = 2;
                return Ok(());
            }
            Err(other) => return Err(other),
        };
        if InternallyGetCanoNumsAndComponentNums(
            heap,
            canonical_globals,
            string_buffer,
            output_control,
            structure.num_atoms,
            canonical_numbers,
            component_numbers,
        )? != 0
        {
            error = 3;
            return Ok(());
        }

        let property_count = u64::try_from(i64::from(original_data.num_inp_atoms).wrapping_add(1))
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        atom_properties = match inchi_calloc::<crate::source_types::OAD_AtProps>(
            heap,
            property_count,
            std::mem::size_of::<crate::source_types::OAD_AtProps>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                inchi_free(heap, canonical_numbers)?;
                canonical_numbers = SourceMutPointer::null();
                inchi_free(heap, component_numbers)?;
                component_numbers = SourceMutPointer::null();
                error = 0;
                return Ok(());
            }
            Err(other) => return Err(other),
        };
        OAD_Polymer_SetAtProps(
            heap,
            structure.polymer.as_const(),
            original_data.at,
            original_data.num_inp_atoms,
            &mut number_of_input_bonds,
            atom_properties,
            canonical_numbers.as_const(),
        )?;

        let unit_count = u64::try_from(polymer.n)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        units2 = match inchi_calloc(
            heap,
            unit_count,
            std::mem::size_of::<SourceMutPointer<crate::source_types::OAD_PolymerUnit>>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                error = 3;
                return Ok(());
            }
            Err(other) => return Err(other),
        };
        let star_count = u64::try_from(polymer.n_pzz)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?;
        old_stars = match inchi_calloc::<i32>(heap, star_count, 4) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                error = 3;
                return Ok(());
            }
            Err(other) => return Err(other),
        };
        for index in 0..polymer.n_pzz {
            let star = *heap
                .slice(polymer.pzz.as_const().offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(old_stars.offset(i64::from(index))?)?[0] = star;
        }

        for index in 0..polymer.n {
            let source_unit = *heap
                .slice(polymer.units.as_const().offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let copy = OAD_PolymerUnit_CreateCopy(heap, source_unit)?;
            heap.slice_mut(units2.offset(i64::from(index))?)?[0] = copy;
            if copy.is_null() {
                error = 4;
                return Ok(());
            }
            copied_units = index.wrapping_add(1);
        }
        unit_numbers = match inchi_calloc::<i32>(heap, unit_count, 4) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => {
                error = 4;
                return Ok(());
            }
            Err(other) => return Err(other),
        };

        if OAD_Polymer_PrepareWorkingSet(
            heap,
            &mut polymer,
            canonical_numbers.as_const(),
            component_numbers.as_const(),
            units2,
            unit_numbers,
        )? != 0
        {
            error = 5;
            return Ok(());
        }

        let layer =
            heap.allocate_model_storage(b"%s\0".iter().map(|byte| *byte as i8).collect())?;
        let layer_value =
            heap.allocate_model_storage(b"/z\0".iter().map(|byte| *byte as i8).collect())?;
        inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            layer.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(layer_value.as_const())],
                ..SourceVaList::default()
            },
        )?;
        let semicolon = heap.allocate_model_storage(vec![b';' as i8, 0])?;
        let mut used_stars = 0_i32;
        let properties = heap.slice(atom_properties.as_const())?.to_vec();
        let canonical = heap.slice(canonical_numbers.as_const())?.to_vec();
        for index in 0..polymer.n {
            let order = *heap
                .slice(unit_numbers.as_const().offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let unit_pointer = *heap
                .slice(units2.as_const().offset(i64::from(order))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let mut unit = heap
                .slice(unit_pointer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            error = OutputINCHI_PolymerLayer_SingleUnit(
                heap,
                &mut unit,
                output_control.bPolymers,
                polymer.n_pzz,
                &mut used_stars,
                &properties,
                &canonical,
                &original_data,
                &structure,
                string_buffer,
            )?;
            if error != 0 {
                return Ok(());
            }
            if index < polymer.n.wrapping_sub(1) {
                inchi_strbuf_printf(
                    heap,
                    Some(string_buffer),
                    semicolon.as_const(),
                    &SourceVaList::default(),
                )?;
            }
        }
        let format =
            heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
        let _ = inchi_ios_print_nodisplay(
            heap,
            output.as_deref_mut(),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    })();

    if !canonical_numbers.is_null() {
        inchi_free(heap, canonical_numbers)?;
    }
    if !component_numbers.is_null() {
        inchi_free(heap, component_numbers)?;
    }
    if !atom_properties.is_null() {
        inchi_free(heap, atom_properties)?;
    }
    if !unit_numbers.is_null() {
        inchi_free(heap, unit_numbers)?;
    }
    if !units2.is_null() {
        for index in 0..copied_units {
            let unit = *heap
                .slice(units2.as_const().offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            OAD_PolymerUnit_Free(heap, unit)?;
        }
        inchi_free(heap, units2)?;
    }
    if !old_stars.is_null() {
        for index in 0..polymer.n_pzz {
            let old = *heap
                .slice(old_stars.as_const().offset(i64::from(index))?)?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(polymer.pzz.offset(i64::from(index))?)?[0] = old;
        }
        inchi_free(heap, old_stars)?;
    }
    body_result?;
    Ok(error)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_HeaderAndNormalization_type(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    output_options: i32,
    basic_or_reconnected: &mut i32,
    component_counts: SourceConstPointer<i32>,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4339 OutputAUXINFO_HeaderAndNormalization_type
    // INCHI✔️❌: int OutputAUXINFO_HeaderAndNormalization_type( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                                INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                                INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                                int              bINChIOutputOptions,
    // INCHI✔️❌:                                                int              *INCHI_basic_or_INCHI_reconnected,
    // INCHI✔️❌:                                                int              num_components2[],
    // INCHI✔️❌:                                                INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                                char             *pLF,
    // INCHI✔️❌:                                                char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* AuxInfo header  */
    // INCHI✔️❌:     if (*INCHI_basic_or_INCHI_reconnected == INCHI_BAS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "AuxInfo=" ); /* in wINChI window, separate INChI: from AuxInfo: with blank line */
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s%s%s",
    // INCHI✔️❌:             ( bINChIOutputOptions & INCHI_OUT_WINCHI_WINDOW ) ? "\n" : "",
    // INCHI✔️❌:             strbuf->pStr, pLF );
    // INCHI✔️❌:         szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_VERS, io->szTag1, &io->bAlways, 0 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "%s", x_curr_ver );
    // INCHI✔️❌:         /* avoid leading slash in plain output */
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (*INCHI_basic_or_INCHI_reconnected == INCHI_REC)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_REC_, io->szTag1, &io->bAlways, 0 );
    // INCHI✔️❌:             inchi_ios_print( out_file, "%s%s", io->szTag1, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* AuxInfo normalization type */
    // INCHI✔️❌:     if (num_components2[0] || num_components2[1])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_NORM, io->szTag1, &io->bAlways, 0 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "%d", ( io->bTautomeric && io->bTautomericOutputAllowed ) ? io->bTautomeric : 0 );
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputAUXINFO_HeaderAndNormalization_type

    fn aux_allocate(
        heap: &mut SourceHeap,
        text: &[u8],
    ) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        Ok(heap
            .allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())?
            .as_const())
    }
    fn aux_append(
        heap: &mut SourceHeap,
        buffer: &mut INCHI_IOS_STRING,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect())?;
        inchi_strbuf_printf(
            heap,
            Some(buffer),
            format.as_const(),
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }
    fn aux_print(
        heap: &mut SourceHeap,
        output: &mut Option<&mut INCHI_IOSTREAM>,
        stdout: SourceMutPointer<FILE>,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect())?;
        let _ = inchi_ios_print_nodisplay(
            heap,
            output.as_deref_mut(),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    let mut aux_labels = vec![INCHI_TAG::default(); 18];
    for (index, plain, xml, always) in [
        (0, b"fixed_H\0".as_slice(), b"fixed-H\0".as_slice(), 0),
        (1, b"isotopic\0".as_slice(), b"isotopic\0".as_slice(), 0),
        (
            2,
            b"abs_stereo_inverted\0".as_slice(),
            b"stereo.abs.inverted\0".as_slice(),
            0,
        ),
        (
            3,
            b"reversibility\0".as_slice(),
            b"reversibility\0".as_slice(),
            0,
        ),
        (4, b"version\0".as_slice(), b"version\0".as_slice(), 1),
        (
            5,
            b"normalization_type\0".as_slice(),
            b"norm-type\0".as_slice(),
            1,
        ),
        (
            17,
            b"reconnected bond(s) to metal(s) part\0".as_slice(),
            b"\0".as_slice(),
            1,
        ),
    ] {
        let plain_label = match index {
            4 => b"\0".as_slice(),
            5 => b"/\0".as_slice(),
            17 => b"/R:\0".as_slice(),
            _ => b"/\0".as_slice(),
        };
        aux_labels[index] = INCHI_TAG {
            szPlainLabel: aux_allocate(heap, plain_label)?,
            szPlainComment: aux_allocate(heap, plain)?,
            szXmlLabel: aux_allocate(heap, xml)?,
            bAlwaysOutput: always,
        };
    }
    let tag_pointer = heap.allocate_model_storage(io.szTag1.to_vec())?;
    let line_feed_empty = aux_allocate(heap, b"\0")?;
    if *basic_or_reconnected == crate::source_types::INCHI_BAS as i32 {
        aux_append(heap, string_buffer, b"AuxInfo=\0", Vec::new())?;
        let blank = if output_options & crate::source_types::INCHI_OUT_WINCHI_WINDOW as i32 != 0 {
            aux_allocate(heap, b"\n\0")?
        } else {
            line_feed_empty
        };
        aux_print(
            heap,
            &mut output,
            stdout,
            b"%s%s%s\0",
            vec![
                SourceFormatArgument::Bytes(blank),
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
        io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_VERS as i32;
        szGetTag(
            heap,
            &aux_labels,
            io.nTag,
            io.bTag1,
            tag_pointer,
            &mut io.bAlways,
            0,
        )?;
        io.szTag1.copy_from_slice(
            heap.slice(tag_pointer.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        let version = heap.allocate_model_storage(
            crate::source_types::local_ichiprt1::x_curr_ver
                .iter()
                .map(|byte| *byte as i8)
                .collect(),
        )?;
        aux_append(
            heap,
            string_buffer,
            b"%s\0",
            vec![SourceFormatArgument::Bytes(version.as_const())],
        )?;
        if str_LineEnd(
            heap,
            tag_pointer.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        aux_print(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    } else if *basic_or_reconnected == crate::source_types::INCHI_REC as i32 {
        io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_REC_ as i32;
        szGetTag(
            heap,
            &aux_labels,
            io.nTag,
            io.bTag1,
            tag_pointer,
            &mut io.bAlways,
            0,
        )?;
        io.szTag1.copy_from_slice(
            heap.slice(tag_pointer.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        aux_print(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(tag_pointer.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    }
    let first_count = *heap
        .slice(component_counts)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let second_count = *heap
        .slice(component_counts.offset(1)?)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if first_count != 0 || second_count != 0 {
        io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_NORM as i32;
        szGetTag(
            heap,
            &aux_labels,
            io.nTag,
            io.bTag1,
            tag_pointer,
            &mut io.bAlways,
            0,
        )?;
        io.szTag1.copy_from_slice(
            heap.slice(tag_pointer.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        let normalization = if io.bTautomeric != 0 && io.bTautomericOutputAllowed != 0 {
            io.bTautomeric
        } else {
            0
        };
        aux_append(
            heap,
            string_buffer,
            b"%d\0",
            vec![SourceFormatArgument::Signed(i64::from(normalization))],
        )?;
        if str_LineEnd(
            heap,
            tag_pointer.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        aux_print(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_OriginalNumbersAndEquivalenceClasses(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    component_counts: SourceConstPointer<i32>,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4396 OutputAUXINFO_OriginalNumbersAndEquivalenceClasses
    // INCHI✔️❌: int OutputAUXINFO_OriginalNumbersAndEquivalenceClasses( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                                         INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                                         INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                                         int              num_components2[],
    // INCHI✔️❌:                                                         INCHI_OUT_CTL   *io,
    // INCHI✔️❌:                                                         char            *pLF,
    // INCHI✔️❌:                                                         char            *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* Original atom numbers in order of canonical numbers */
    // INCHI✔️❌:     if (num_components2[0] || num_components2[1])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szGetTag( AuxLbl, io->nTag,
    // INCHI✔️❌:                  io->bTag1 = ( io->bSecondNonTautPass ? AL_FIXN : AL_ANBR ) | io->bFhTag, io->szTag1, &io->bAlways, 0 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:         io->tot_len = 0;
    // INCHI✔️❌:         /* Original numbering output */
    // INCHI✔️❌:         io->tot_len = str_AuxNumb( pCG, io->pINChISort, io->pINChISort2,
    // INCHI✔️❌:                                    strbuf, &io->bOverflow,
    // INCHI✔️❌:                                    io->bOutType, io->TAUT_MODE, io->num_components,
    // INCHI✔️❌:                                    io->bSecondNonTautPass, io->bOmitRepetitions );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /*
    // INCHI✔️❌:         Symmetry numbers (constit. equivalence)    /E:
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (io->bAtomEqu[io->iCurTautMode])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*  aux equ atoms */
    // INCHI✔️❌:         /* 1. Compare to tautomeric equivalence (in case of second, non-taut, pass only) */
    // INCHI✔️❌:         /* 2. Compare to the previous component if (1) failed to find equivalence */
    // INCHI✔️❌:         szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_AEQU | io->bFhTag, io->szTag1, &io->bAlways, 0 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:         io->tot_len = 0;
    // INCHI✔️❌:         io->tot_len = str_AuxEqu( io->pINChISort, io->pINChISort2,
    // INCHI✔️❌:                               strbuf, &io->bOverflow, io->bOutType, io->TAUT_MODE,
    // INCHI✔️❌:                               io->num_components, io->bSecondNonTautPass,
    // INCHI✔️❌:                               io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (io->bPlainTextTags == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print( out_file, "/" );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputAUXINFO_OriginalNumbersAndEquivalenceClasses

    fn allocate_text(
        heap: &mut SourceHeap,
        text: &[u8],
    ) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        Ok(heap
            .allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())?
            .as_const())
    }
    fn aux_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
        let rows: [(&[u8], &[u8], &[u8], i32); 18] = [
            (b"/\0", b"fixed_H\0", b"fixed-H\0", 0),
            (b"/\0", b"isotopic\0", b"isotopic\0", 0),
            (
                b"/\0",
                b"abs_stereo_inverted\0",
                b"stereo.abs.inverted\0",
                0,
            ),
            (b"/\0", b"reversibility\0", b"reversibility\0", 0),
            (b"\0", b"version\0", b"version\0", 1),
            (b"/\0", b"normalization_type\0", b"norm-type\0", 1),
            (b"/N:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/E:\0", b"atom_equivalence\0", b"atom.equivalence\0", 0),
            (b"/gE:\0", b"group_equivalence\0", b"group.equivalence\0", 0),
            (b"/it:\0", b"sp3\0", b"sp3\0", 0),
            (b"/iN:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 0),
            (
                b"/CRV:\0",
                b"charge_radical_valence\0",
                b"charges-rad-val\0",
                0,
            ),
            (b"/rA:\0", b"atoms\0", b"atoms\0", 0),
            (b"/rB:\0", b"bonds\0", b"bonds\0", 0),
            (b"/rC:\0", b"xyz\0", b"xyz\0", 0),
            (b"/F:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/I:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (
                b"/R:\0",
                b"reconnected bond(s) to metal(s) part\0",
                b"\0",
                1,
            ),
        ];
        rows.into_iter()
            .map(|(plain, comment, xml, always)| {
                Ok(INCHI_TAG {
                    szPlainLabel: allocate_text(heap, plain)?,
                    szPlainComment: allocate_text(heap, comment)?,
                    szXmlLabel: allocate_text(heap, xml)?,
                    bAlwaysOutput: always,
                })
            })
            .collect()
    }
    fn print_values(
        heap: &mut SourceHeap,
        output: &mut Option<&mut INCHI_IOSTREAM>,
        stdout: SourceMutPointer<FILE>,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print(
            heap,
            output.as_deref_mut(),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    let labels = aux_labels(heap)?;
    let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
    let count0 = heap
        .slice(component_counts)?
        .first()
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let count1 = heap
        .slice(component_counts.offset(1)?)?
        .first()
        .copied()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if count0 != 0 || count1 != 0 {
        io.bTag1 = if io.bSecondNonTautPass != 0 {
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_FIXN as i32
        } else {
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_ANBR as i32
        } | io.bFhTag;
        szGetTag(heap, &labels, io.nTag, io.bTag1, tag, &mut io.bAlways, 0)?;
        io.szTag1.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxNumb(
            heap,
            canonical_globals,
            io.pINChISort,
            io.pINChISort2,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bSecondNonTautPass,
            io.bOmitRepetitions,
        )?;
        if str_LineEnd(
            heap,
            tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    }
    let taut_index =
        usize::try_from(io.iCurTautMode).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atom_equivalence = *io
        .bAtomEqu
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom_equivalence != 0 {
        io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_AEQU as i32 | io.bFhTag;
        szGetTag(heap, &labels, io.nTag, io.bTag1, tag, &mut io.bAlways, 0)?;
        io.szTag1.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxEqu(
            heap,
            io.pINChISort,
            io.pINChISort2,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bSecondNonTautPass,
            io.bOmitRepetitions,
            io.bUseMulipliers,
        )?;
        if str_LineEnd(
            heap,
            tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    } else if io.bPlainTextTags == 1 {
        print_values(heap, &mut output, stdout, b"/\0", Vec::new())?;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_TautomericGroupsEquivalence(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    io: &mut INCHI_OUT_CTL,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4461 OutputAUXINFO_TautomericGroupsEquivalence
    // INCHI✔️❌: int OutputAUXINFO_TautomericGroupsEquivalence( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                                INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                                INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                                INCHI_OUT_CTL    *io )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (io->bTautomericOutputAllowed && io->bTautomeric && io->bTautEqu[io->iCurTautMode] && !io->bSecondNonTautPass)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /*-- Tautomeric groups constitutional equivalence */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*-- aux tgroup equ */
    // INCHI✔️❌:         szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_GEQU | io->bFhTag, io->szTag1, &io->bAlways, 0 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:         io->tot_len = str_AuxTgroupEqu( io->pINChISort,
    // INCHI✔️❌:                                     strbuf, &io->bOverflow, io->bOutType, io->TAUT_MODE,
    // INCHI✔️❌:                                     io->num_components, io->bUseMulipliers );
    // INCHI✔️❌:         if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s", strbuf->pStr );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (io->bTautomericOutputAllowed && io->bTautomeric)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io->bPlainTextTags == 1)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_ios_print( out_file, "/" );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputAUXINFO_TautomericGroupsEquivalence

    fn allocate_text(
        heap: &mut SourceHeap,
        text: &[u8],
    ) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        Ok(heap
            .allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())?
            .as_const())
    }

    fn aux_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
        let rows: [(&[u8], &[u8], &[u8], i32); 18] = [
            (b"/\0", b"fixed_H\0", b"fixed-H\0", 0),
            (b"/\0", b"isotopic\0", b"isotopic\0", 0),
            (
                b"/\0",
                b"abs_stereo_inverted\0",
                b"stereo.abs.inverted\0",
                0,
            ),
            (b"/\0", b"reversibility\0", b"reversibility\0", 0),
            (b"\0", b"version\0", b"version\0", 1),
            (b"/\0", b"normalization_type\0", b"norm-type\0", 1),
            (b"/N:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/E:\0", b"atom_equivalence\0", b"atom.equivalence\0", 0),
            (b"/gE:\0", b"group_equivalence\0", b"group.equivalence\0", 0),
            (b"/it:\0", b"sp3\0", b"sp3\0", 0),
            (b"/iN:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 0),
            (
                b"/CRV:\0",
                b"charge_radical_valence\0",
                b"charges-rad-val\0",
                0,
            ),
            (b"/rA:\0", b"atoms\0", b"atoms\0", 0),
            (b"/rB:\0", b"bonds\0", b"bonds\0", 0),
            (b"/rC:\0", b"xyz\0", b"xyz\0", 0),
            (b"/F:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/I:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (
                b"/R:\0",
                b"reconnected bond(s) to metal(s) part\0",
                b"\0",
                1,
            ),
        ];
        rows.into_iter()
            .map(|(plain, comment, xml, always)| {
                Ok(INCHI_TAG {
                    szPlainLabel: allocate_text(heap, plain)?,
                    szPlainComment: allocate_text(heap, comment)?,
                    szXmlLabel: allocate_text(heap, xml)?,
                    bAlwaysOutput: always,
                })
            })
            .collect()
    }

    fn print_values(
        heap: &mut SourceHeap,
        output: &mut Option<&mut INCHI_IOSTREAM>,
        stdout: SourceMutPointer<FILE>,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print(
            heap,
            output.as_deref_mut(),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    let taut_index =
        usize::try_from(io.iCurTautMode).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let taut_equivalence = *io
        .bTautEqu
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if io.bTautomericOutputAllowed != 0
        && io.bTautomeric != 0
        && taut_equivalence != 0
        && io.bSecondNonTautPass == 0
    {
        let labels = aux_labels(heap)?;
        let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
        io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_GEQU as i32 | io.bFhTag;
        szGetTag(heap, &labels, io.nTag, io.bTag1, tag, &mut io.bAlways, 0)?;
        io.szTag1.copy_from_slice(
            heap.slice(tag.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxTgroupEqu(
            heap,
            io.pINChISort,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bUseMulipliers,
        )?;
        if str_LineEnd(
            heap,
            tag.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s\0",
            vec![SourceFormatArgument::Bytes(string_buffer.pStr.as_const())],
        )?;
    } else if io.bTautomericOutputAllowed != 0 && io.bTautomeric != 0 && io.bPlainTextTags == 1 {
        print_values(heap, &mut output, stdout, b"/\0", Vec::new())?;
    }

    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_Stereo(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4500 OutputAUXINFO_Stereo
    // INCHI✔️❌: int OutputAUXINFO_Stereo( CANON_GLOBALS     *pCG,
    // INCHI✔️❌:                            INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                            INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                            INCHI_OUT_CTL   *io,
    // INCHI✔️❌:                            char             *pLF,
    // INCHI✔️❌:                            char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /*--    Inverted stereo -- sp3 only + canonical numbering
    // INCHI✔️❌:     */
    // INCHI✔️❌:     if (io->bInvStereo[io->iCurTautMode])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_STER | io->bFhTag, io->szTag1, &io->bAlways, 0 );
    // INCHI✔️❌:         /*-- inverted sp3 start tag */
    // INCHI✔️❌:         szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_SP3I, io->szTag2, &io->bAlways, 0 );
    // INCHI✔️❌:         inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:         io->tot_len = str_AuxInvSp3( io->pINChISort, io->pINChISort2, strbuf,
    // INCHI✔️❌:                                  &io->bOverflow, io->bOutType, io->TAUT_MODE, io->num_components,
    // INCHI✔️❌:                                  io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
    // INCHI✔️❌:         if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:
    // INCHI✔️❌:         /*-- inverted sp3  canonical numbering */
    // INCHI✔️❌:         if (io->bInvStereoOrigNumb[io->iCurTautMode])
    // INCHI✔️❌:         {
    // INCHI✔️❌:             szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_SP3N, io->szTag2, &io->bAlways, 0 );
    // INCHI✔️❌:             inchi_strbuf_reset( strbuf ); io->tot_len = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:             io->tot_len = str_AuxInvSp3Numb( pCG, io->pINChISort, io->pINChISort2,
    // INCHI✔️❌:                                          strbuf, &io->bOverflow, io->bOutType,
    // INCHI✔️❌:                                          io->TAUT_MODE, io->num_components,
    // INCHI✔️❌:                                          io->bSecondNonTautPass, io->bOmitRepetitions );
    // INCHI✔️❌:
    // INCHI✔️❌:             if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 1;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (io->bPlainTextTags == 1) inchi_ios_print( out_file, "/" );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (io->bPlainTextTags == 1)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_ios_print( out_file, "//" );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* Inverted stereo -- sp3 only + canonical numbering */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputAUXINFO_Stereo

    fn allocate_text(
        heap: &mut SourceHeap,
        text: &[u8],
    ) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        Ok(heap
            .allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())?
            .as_const())
    }
    fn aux_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
        let rows: [(&[u8], &[u8], &[u8], i32); 18] = [
            (b"/\0", b"fixed_H\0", b"fixed-H\0", 0),
            (b"/\0", b"isotopic\0", b"isotopic\0", 0),
            (
                b"/\0",
                b"abs_stereo_inverted\0",
                b"stereo.abs.inverted\0",
                0,
            ),
            (b"/\0", b"reversibility\0", b"reversibility\0", 0),
            (b"\0", b"version\0", b"version\0", 1),
            (b"/\0", b"normalization_type\0", b"norm-type\0", 1),
            (b"/N:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/E:\0", b"atom_equivalence\0", b"atom.equivalence\0", 0),
            (b"/gE:\0", b"group_equivalence\0", b"group.equivalence\0", 0),
            (b"/it:\0", b"sp3\0", b"sp3\0", 0),
            (b"/iN:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 0),
            (
                b"/CRV:\0",
                b"charge_radical_valence\0",
                b"charges-rad-val\0",
                0,
            ),
            (b"/rA:\0", b"atoms\0", b"atoms\0", 0),
            (b"/rB:\0", b"bonds\0", b"bonds\0", 0),
            (b"/rC:\0", b"xyz\0", b"xyz\0", 0),
            (b"/F:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/I:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (
                b"/R:\0",
                b"reconnected bond(s) to metal(s) part\0",
                b"\0",
                1,
            ),
        ];
        rows.into_iter()
            .map(|(plain, comment, xml, always)| {
                Ok(INCHI_TAG {
                    szPlainLabel: allocate_text(heap, plain)?,
                    szPlainComment: allocate_text(heap, comment)?,
                    szXmlLabel: allocate_text(heap, xml)?,
                    bAlwaysOutput: always,
                })
            })
            .collect()
    }
    fn print_values(
        heap: &mut SourceHeap,
        output: &mut Option<&mut INCHI_IOSTREAM>,
        stdout: SourceMutPointer<FILE>,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print(
            heap,
            output.as_deref_mut(),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    let taut_index =
        usize::try_from(io.iCurTautMode).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let inverted = *io
        .bInvStereo
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if inverted != 0 {
        let labels = aux_labels(heap)?;
        let tag1 = heap.allocate_model_storage(io.szTag1.to_vec())?;
        let tag2 = heap.allocate_model_storage(io.szTag2.to_vec())?;
        io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_STER as i32 | io.bFhTag;
        szGetTag(heap, &labels, io.nTag, io.bTag1, tag1, &mut io.bAlways, 0)?;
        io.szTag1.copy_from_slice(
            heap.slice(tag1.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        io.bTag2 = io.bTag1 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_SP3I as i32;
        szGetTag(heap, &labels, io.nTag, io.bTag2, tag2, &mut io.bAlways, 0)?;
        io.szTag2.copy_from_slice(
            heap.slice(tag2.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxInvSp3(
            heap,
            io.pINChISort,
            io.pINChISort2,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bSecondNonTautPass,
            io.bOmitRepetitions,
            io.bUseMulipliers,
        )?;
        if str_LineEnd(
            heap,
            tag2.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;

        let inverted_numbers = *io
            .bInvStereoOrigNumb
            .get(taut_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if inverted_numbers != 0 {
            io.bTag2 = io.bTag1 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_SP3N as i32;
            szGetTag(heap, &labels, io.nTag, io.bTag2, tag2, &mut io.bAlways, 0)?;
            io.szTag2.copy_from_slice(
                heap.slice(tag2.as_const())?
                    .get(..64)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            io.tot_len = crate::source::base::ichiprt3::str_AuxInvSp3Numb(
                heap,
                canonical_globals,
                io.pINChISort,
                io.pINChISort2,
                string_buffer,
                &mut io.bOverflow,
                io.bOutType,
                io.TAUT_MODE,
                io.num_components,
                io.bSecondNonTautPass,
                io.bOmitRepetitions,
            )?;
            if str_LineEnd(
                heap,
                tag2.as_const(),
                &mut io.bOverflow,
                string_buffer,
                -1,
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(1);
            }
            print_values(
                heap,
                &mut output,
                stdout,
                b"%s%s\0",
                vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
            )?;
        } else if io.bPlainTextTags == 1 {
            print_values(heap, &mut output, stdout, b"/\0", Vec::new())?;
        }
    } else if io.bPlainTextTags == 1 {
        print_values(heap, &mut output, stdout, b"//\0", Vec::new())?;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_IsotopicInfo(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    basic_or_reconnected: &mut i32,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4560 OutputAUXINFO_IsotopicInfo
    // INCHI✔️❌: complete configured source frame follows verbatim.
    /*
    int OutputAUXINFO_IsotopicInfo( CANON_GLOBALS    *pCG,
                                    INCHI_IOSTREAM   *out_file,
                                    INCHI_IOS_STRING *strbuf,
                                    int              *INCHI_basic_or_INCHI_reconnected,
                                    INCHI_OUT_CTL    *io,
                                    char             *pLF,
                                    char             *pTAB )
    {
        int i;

        /* if InChI Fixed-H isotopic is empty, then do not output corresponding AuxInfo */

        i = io->bSecondNonTautPass &&
            ( *io->pSortPrintINChIFlags & ( ( *INCHI_basic_or_INCHI_reconnected == INCHI_BAS ) ? FLAG_SORT_PRINT_NO_IFIX_H_BAS :
                FLAG_SORT_PRINT_NO_IFIX_H_REC ) );

        if (io->bIsotopic && !i &&
            ( io->bIsotopicOrigNumb[io->iCurTautMode] ||
                io->bIsotopicAtomEqu[io->iCurTautMode] ||
                (io->bTautomericOutputAllowed && io->bTautomeric && io->bIsotopicTautEqu[io->iCurTautMode]) ||
                (io->bInvIsotopicStereo[io->iCurTautMode]
                && ( io->bIgn_UU_Sp3_Iso[io->iCurTautMode])) || io->bIgn_UU_Sp2_Iso[io->iCurTautMode] ) ) /* djb-rwth: addressing LLVM warnings */
        {
            /*-- isotopic aux info header */
            szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_ISOT | io->bFhTag, io->szTag1, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf ); /* pStr[io->tot_len = 0] = '\0'; */
            /*-- Original atom numbers in order of isotopic canonical numbers */
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_ISON, io->szTag2, &io->bAlways, 0 );
            if (io->bIsotopicOrigNumb[io->iCurTautMode])
            {
                inchi_strbuf_reset( strbuf );
                io->tot_len = 0;
                io->tot_len = str_AuxIsoNumb( pCG, io->pINChISort, io->pINChISort2,
                                          strbuf, &io->bOverflow, io->bOutType,
                                          io->TAUT_MODE, io->num_components,
                                          io->bSecondNonTautPass, io->bOmitRepetitions );
                if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
                {
                    return 1;
                }
                inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
            }
            else
            {
                /*if ( io->bPlainTextTags == 1 ) inchi_ios_print( out_file, "/" );*/
                inchi_ios_print( out_file, "%s%s", io->szTag2, pLF ); /* mark isotopic output */
            }

            /*-- Isotopic symmetry */
            if (io->bIsotopicAtomEqu[io->iCurTautMode])
            {
                /*-- atoms */
                szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_AEQU, io->szTag2, &io->bAlways, 0 );
                inchi_strbuf_reset( strbuf ); io->tot_len = 0;
                io->tot_len = str_AuxIsoEqu( io->pINChISort, io->pINChISort2,
                                         strbuf,
                                         &io->bOverflow, io->bOutType, io->TAUT_MODE, io->num_components,
                                         io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
                if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -2/*was -1: Fix15*/, io->bPlainTextTags ))
                {
                    return 1;
                }
                inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
            }
            else
            {
                if (io->bPlainTextTags == 1)
                {
                    inchi_ios_print( out_file, "/" );
                }
            }

            /*-- Tautomeric groups, isotopic */
            if (io->bTautomericOutputAllowed && io->bTautomeric && io->bIsotopicTautEqu[io->iCurTautMode])
            {
                /*-- Isotopic tautomeric groups equivalence */
                szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_GEQU, io->szTag2, &io->bAlways, 0 );
                inchi_strbuf_reset( strbuf ); io->tot_len = 0;
                io->tot_len = str_AuxIsoTgroupEqu( io->pINChISort,
                                               strbuf, &io->bOverflow,
                                               io->bOutType, io->TAUT_MODE, io->num_components,
                                               io->bOmitRepetitions, io->bUseMulipliers );
                if (str_LineEnd( io->szTag2, &io->bOverflow, strbuf, -2/*was -1: Fix15*/, io->bPlainTextTags ))
                {
                    return 1;
                }
                inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
            }
            else
            {
                if (io->bTautomericOutputAllowed && io->bTautomeric)
                {
                    if (io->bPlainTextTags == 1)
                    {
                        inchi_ios_print( out_file, "/" );
                    }
                }
            }
            /*-- Isotopic inverted stereo */
            if (io->bInvIsotopicStereo[io->iCurTautMode])
            {
                szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_STER, io->szTag2, &io->bAlways, 0 );
                /*-- inverted isotopic sp3 start tag */
                szGetTag( AuxLbl, io->nTag, io->bTag3 = io->bTag2 | AL_SP3I, io->szTag3, &io->bAlways, 0 );
                inchi_strbuf_reset( strbuf ); io->tot_len = 0;
                io->tot_len = str_AuxInvIsoSp3( io->pINChISort, io->pINChISort2,
                                            strbuf, &io->bOverflow,
                                            io->bOutType, io->TAUT_MODE, io->num_components,
                                            io->bSecondNonTautPass, io->bOmitRepetitions, io->bUseMulipliers );
                if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
                {
                    return 1;
                }
                inchi_ios_print( out_file, "%s", strbuf->pStr );
                /*-- inverted isotopic sp3  canonical numbering */
                if (io->bInvIsotopicStereoOrigNumb[io->iCurTautMode])
                {
                    szGetTag( AuxLbl, io->nTag, io->bTag3 = io->bTag2 | AL_SP3N, io->szTag3, &io->bAlways, 0 );
                    inchi_strbuf_reset( strbuf ); io->tot_len = 0;
                    io->tot_len = str_AuxInvIsoSp3Numb( pCG, io->pINChISort, io->pINChISort2,
                                                    strbuf, &io->bOverflow,
                                                    io->bOutType, io->TAUT_MODE,
                                                    io->num_components,
                                                    io->bSecondNonTautPass,
                                                    io->bOmitRepetitions );

                    if (str_LineEnd( io->szTag3, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
                    {
                        return 1;
                    }
                    inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
                }
                else
                {
                    if (io->bPlainTextTags == 1)
                    {
                        inchi_ios_print( out_file, "/" );
                    }
                }
            }
            else
            {
                if (io->bPlainTextTags == 1)
                {
                    inchi_ios_print( out_file, "//" );
                }
            }
            /*-- totally omitted undefined/unknown isotopic stereo */
        } /* Aux info isotopic */

        return 0;
    }
    */
    // END INCHI C FUNCTION: OutputAUXINFO_IsotopicInfo

    fn allocate_text(
        heap: &mut SourceHeap,
        text: &[u8],
    ) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        Ok(heap
            .allocate_model_storage(text.iter().map(|byte| *byte as i8).collect())?
            .as_const())
    }
    fn aux_labels(heap: &mut SourceHeap) -> Result<Vec<INCHI_TAG>, SourceHeapError> {
        let rows: [(&[u8], &[u8], &[u8], i32); 18] = [
            (b"/\0", b"fixed_H\0", b"fixed-H\0", 0),
            (b"/\0", b"isotopic\0", b"isotopic\0", 0),
            (
                b"/\0",
                b"abs_stereo_inverted\0",
                b"stereo.abs.inverted\0",
                0,
            ),
            (b"/\0", b"reversibility\0", b"reversibility\0", 0),
            (b"\0", b"version\0", b"version\0", 1),
            (b"/\0", b"normalization_type\0", b"norm-type\0", 1),
            (b"/N:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/E:\0", b"atom_equivalence\0", b"atom.equivalence\0", 0),
            (b"/gE:\0", b"group_equivalence\0", b"group.equivalence\0", 0),
            (b"/it:\0", b"sp3\0", b"sp3\0", 0),
            (b"/iN:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 0),
            (
                b"/CRV:\0",
                b"charge_radical_valence\0",
                b"charges-rad-val\0",
                0,
            ),
            (b"/rA:\0", b"atoms\0", b"atoms\0", 0),
            (b"/rB:\0", b"bonds\0", b"bonds\0", 0),
            (b"/rC:\0", b"xyz\0", b"xyz\0", 0),
            (b"/F:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (b"/I:\0", b"original_atom_numbers\0", b"atom.orig-nbr\0", 1),
            (
                b"/R:\0",
                b"reconnected bond(s) to metal(s) part\0",
                b"\0",
                1,
            ),
        ];
        rows.into_iter()
            .map(|(plain, comment, xml, always)| {
                Ok(INCHI_TAG {
                    szPlainLabel: allocate_text(heap, plain)?,
                    szPlainComment: allocate_text(heap, comment)?,
                    szXmlLabel: allocate_text(heap, xml)?,
                    bAlwaysOutput: always,
                })
            })
            .collect()
    }
    fn print_values(
        heap: &mut SourceHeap,
        output: &mut Option<&mut INCHI_IOSTREAM>,
        stdout: SourceMutPointer<FILE>,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().map(|byte| *byte as i8).collect())?;
        inchi_ios_print(
            heap,
            output.as_deref_mut(),
            stdout,
            format.as_const(),
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    let fixed_h_isotopic_empty = if io.bSecondNonTautPass != 0 {
        let flags = *heap
            .slice(io.pSortPrintINChIFlags.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mask = if *basic_or_reconnected == crate::source_types::INCHI_BAS as i32 {
            crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_BAS as i32
        } else {
            crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_REC as i32
        };
        flags & mask
    } else {
        0
    };
    if io.bIsotopic == 0 || fixed_h_isotopic_empty != 0 {
        return Ok(0);
    }
    let taut_index =
        usize::try_from(io.iCurTautMode).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let isotopic_numbers = *io
        .bIsotopicOrigNumb
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let isotopic_equivalence = *io
        .bIsotopicAtomEqu
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let isotopic_taut_equivalence = *io
        .bIsotopicTautEqu
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let inverted_isotopic_stereo = *io
        .bInvIsotopicStereo
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ignore_undefined_sp3 = *io
        .bIgn_UU_Sp3_Iso
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let ignore_undefined_sp2 = *io
        .bIgn_UU_Sp2_Iso
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let taut_output = io.bTautomericOutputAllowed != 0 && io.bTautomeric != 0;
    if isotopic_numbers == 0
        && isotopic_equivalence == 0
        && !(taut_output && isotopic_taut_equivalence != 0)
        && !(inverted_isotopic_stereo != 0 && ignore_undefined_sp3 != 0)
        && ignore_undefined_sp2 == 0
    {
        return Ok(0);
    }

    let labels = aux_labels(heap)?;
    let tag1 = heap.allocate_model_storage(io.szTag1.to_vec())?;
    let tag2 = heap.allocate_model_storage(io.szTag2.to_vec())?;
    let tag3 = heap.allocate_model_storage(io.szTag3.to_vec())?;
    io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_ISOT as i32 | io.bFhTag;
    szGetTag(heap, &labels, io.nTag, io.bTag1, tag1, &mut io.bAlways, 0)?;
    io.szTag1.copy_from_slice(
        heap.slice(tag1.as_const())?
            .get(..64)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    inchi_strbuf_reset(heap, Some(string_buffer))?;
    io.bTag2 = io.bTag1 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_ISON as i32;
    szGetTag(heap, &labels, io.nTag, io.bTag2, tag2, &mut io.bAlways, 0)?;
    io.szTag2.copy_from_slice(
        heap.slice(tag2.as_const())?
            .get(..64)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    if isotopic_numbers != 0 {
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxIsoNumb(
            heap,
            canonical_globals,
            io.pINChISort,
            io.pINChISort2,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bSecondNonTautPass,
            io.bOmitRepetitions,
        )?;
        if str_LineEnd(
            heap,
            tag2.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    } else {
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(tag2.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    }

    if isotopic_equivalence != 0 {
        io.bTag2 = io.bTag1 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_AEQU as i32;
        szGetTag(heap, &labels, io.nTag, io.bTag2, tag2, &mut io.bAlways, 0)?;
        io.szTag2.copy_from_slice(
            heap.slice(tag2.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxIsoEqu(
            heap,
            io.pINChISort,
            io.pINChISort2,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bSecondNonTautPass,
            io.bOmitRepetitions,
            io.bUseMulipliers,
        )?;
        if str_LineEnd(
            heap,
            tag2.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -2,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    } else if io.bPlainTextTags == 1 {
        print_values(heap, &mut output, stdout, b"/\0", Vec::new())?;
    }

    if taut_output && isotopic_taut_equivalence != 0 {
        io.bTag2 = io.bTag1 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_GEQU as i32;
        szGetTag(heap, &labels, io.nTag, io.bTag2, tag2, &mut io.bAlways, 0)?;
        io.szTag2.copy_from_slice(
            heap.slice(tag2.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxIsoTgroupEqu(
            heap,
            io.pINChISort,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bOmitRepetitions,
            io.bUseMulipliers,
        )?;
        if str_LineEnd(
            heap,
            tag2.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -2,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s%s\0",
            vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
        )?;
    } else if taut_output && io.bPlainTextTags == 1 {
        print_values(heap, &mut output, stdout, b"/\0", Vec::new())?;
    }

    if inverted_isotopic_stereo != 0 {
        io.bTag2 = io.bTag1 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_STER as i32;
        szGetTag(heap, &labels, io.nTag, io.bTag2, tag2, &mut io.bAlways, 0)?;
        io.szTag2.copy_from_slice(
            heap.slice(tag2.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        io.bTag3 = io.bTag2 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_SP3I as i32;
        szGetTag(heap, &labels, io.nTag, io.bTag3, tag3, &mut io.bAlways, 0)?;
        io.szTag3.copy_from_slice(
            heap.slice(tag3.as_const())?
                .get(..64)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        io.tot_len = 0;
        io.tot_len = crate::source::base::ichiprt3::str_AuxInvIsoSp3(
            heap,
            io.pINChISort,
            io.pINChISort2,
            string_buffer,
            &mut io.bOverflow,
            io.bOutType,
            io.TAUT_MODE,
            io.num_components,
            io.bSecondNonTautPass,
            io.bOmitRepetitions,
            io.bUseMulipliers,
        )?;
        if str_LineEnd(
            heap,
            tag3.as_const(),
            &mut io.bOverflow,
            string_buffer,
            -1,
            io.bPlainTextTags,
        )? != 0
        {
            return Ok(1);
        }
        print_values(
            heap,
            &mut output,
            stdout,
            b"%s\0",
            vec![SourceFormatArgument::Bytes(string_buffer.pStr.as_const())],
        )?;
        let inverted_numbers = *io
            .bInvIsotopicStereoOrigNumb
            .get(taut_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if inverted_numbers != 0 {
            io.bTag3 = io.bTag2 | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_SP3N as i32;
            szGetTag(heap, &labels, io.nTag, io.bTag3, tag3, &mut io.bAlways, 0)?;
            io.szTag3.copy_from_slice(
                heap.slice(tag3.as_const())?
                    .get(..64)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            inchi_strbuf_reset(heap, Some(string_buffer))?;
            io.tot_len = 0;
            io.tot_len = crate::source::base::ichiprt3::str_AuxInvIsoSp3Numb(
                heap,
                canonical_globals,
                io.pINChISort,
                io.pINChISort2,
                string_buffer,
                &mut io.bOverflow,
                io.bOutType,
                io.TAUT_MODE,
                io.num_components,
                io.bSecondNonTautPass,
                io.bOmitRepetitions,
            )?;
            if str_LineEnd(
                heap,
                tag3.as_const(),
                &mut io.bOverflow,
                string_buffer,
                -1,
                io.bPlainTextTags,
            )? != 0
            {
                return Ok(1);
            }
            print_values(
                heap,
                &mut output,
                stdout,
                b"%s%s\0",
                vec![
                    SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                    SourceFormatArgument::Bytes(line_feed),
                ],
            )?;
        } else if io.bPlainTextTags == 1 {
            print_values(heap, &mut output, stdout, b"/\0", Vec::new())?;
        }
    } else if io.bPlainTextTags == 1 {
        print_values(heap, &mut output, stdout, b"//\0", Vec::new())?;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_ChargesRadicalsAndUnusualValences(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4717 OutputAUXINFO_ChargesRadicalsAndUnusualValences
    // INCHI✔️❌: complete configured source frame follows verbatim.
    /*
    int OutputAUXINFO_ChargesRadicalsAndUnusualValences( CANON_GLOBALS    *pCG,
                                                         INCHI_IOSTREAM   *out_file,
                                                         INCHI_IOS_STRING *strbuf,
                                                         INCHI_OUT_CTL    *io,
                                                         char             *pLF,
                                                         char             *pTAB )
    {
        if (!io->bSecondNonTautPass && io->bChargesRadVal[io->iCurTautMode])
        {
            /*  aux equ atoms */
            /* 1. Compare to tautomeric equivalence (in case of second, non-taut, pass only) */
            /* 2. Compare to the previous component if (1) failed to find equivalence */
            szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_CRV_ | io->bFhTag, io->szTag1, &io->bAlways, 0 );

            inchi_strbuf_reset( strbuf );
            io->tot_len = 0;

            io->tot_len = str_AuxChargeRadVal( io->pINChISort, strbuf,
                                               &io->bOverflow, io->bOutType, io->TAUT_MODE,
                                               io->num_components, io->bUseMulipliers );

            if (str_LineEnd( io->szTag1, &io->bOverflow, strbuf, -1, io->bPlainTextTags ))
            {
                return 1;
            }

            inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: OutputAUXINFO_ChargesRadicalsAndUnusualValences

    if io.bSecondNonTautPass != 0 {
        return Ok(0);
    }
    let taut_index =
        usize::try_from(io.iCurTautMode).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if *io
        .bChargesRadVal
        .get(taut_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        == 0
    {
        return Ok(0);
    }
    let allocate_text = |heap: &mut SourceHeap, bytes: &[u8]| {
        Ok::<_, SourceHeapError>(
            heap.allocate_model_storage(bytes.iter().map(|byte| *byte as i8).collect())?
                .as_const(),
        )
    };
    let mut labels = vec![
        INCHI_TAG::default();
        crate::source_types::local_ichiprt1::tagAuxLblOrd_AL_MAX_ORD as usize
    ];
    labels[0] = INCHI_TAG {
        szPlainLabel: allocate_text(heap, b"/\0")?,
        szPlainComment: allocate_text(heap, b"fixed_H\0")?,
        szXmlLabel: allocate_text(heap, b"fixed-H\0")?,
        bAlwaysOutput: 0,
    };
    labels[11] = INCHI_TAG {
        szPlainLabel: allocate_text(heap, b"/CRV:\0")?,
        szPlainComment: allocate_text(heap, b"charge_radical_valence\0")?,
        szXmlLabel: allocate_text(heap, b"charges-rad-val\0")?,
        bAlwaysOutput: 0,
    };
    let tag = heap.allocate_model_storage(io.szTag1.to_vec())?;
    io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_CRV_ as i32 | io.bFhTag;
    szGetTag(heap, &labels, io.nTag, io.bTag1, tag, &mut io.bAlways, 0)?;
    io.szTag1
        .copy_from_slice(&heap.slice(tag.as_const())?[..64]);
    inchi_strbuf_reset(heap, Some(string_buffer))?;
    io.tot_len = 0;
    io.tot_len = crate::source::base::ichiprt3::str_AuxChargeRadVal(
        heap,
        io.pINChISort,
        string_buffer,
        &mut io.bOverflow,
        io.bOutType,
        io.TAUT_MODE,
        io.num_components,
        io.bUseMulipliers,
    )?;
    if str_LineEnd(
        heap,
        tag.as_const(),
        &mut io.bOverflow,
        string_buffer,
        -1,
        io.bPlainTextTags,
    )? != 0
    {
        return Ok(1);
    }
    let format =
        heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, b'%' as i8, b's' as i8, 0])?;
    inchi_ios_print(
        heap,
        output.as_deref_mut(),
        stdout,
        format.as_const(),
        &SourceVaList {
            arguments: vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
            ..SourceVaList::default()
        },
    )?;
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_ReversibilityInfo(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    original_structure: Option<&ORIG_STRUCT>,
    io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4753 OutputAUXINFO_ReversibilityInfo
    // INCHI✔️❌: complete configured source frame follows verbatim.
    /*
    int OutputAUXINFO_ReversibilityInfo( CANON_GLOBALS    *pCG,
                                         INCHI_IOSTREAM   *out_file,
                                         INCHI_IOS_STRING *strbuf,
                                         ORIG_STRUCT      *pOrigStruct,
                                         INCHI_OUT_CTL    *io,
                                         char             *pLF,
                                         char             *pTAB )
    {
        if (!io->bSecondNonTautPass &&
             pOrigStruct && pOrigStruct->num_atoms &&
             pOrigStruct->szAtoms
             && pOrigStruct->szBonds
             && pOrigStruct->szCoord)
        {
            int length, cur_pos, line_len, last_pos, nMaxLineLen;
            char *p;
            nMaxLineLen = inchi_min( 80, strbuf->nAllocatedLength ); /* restrict line length to 80 characters */

            szGetTag( AuxLbl, io->nTag, io->bTag1 = AL_REVR | io->bFhTag, io->szTag1, &io->bAlways, 0 );

            /* Atoms /A: */
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_ATMR, io->szTag2, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf );
            inchi_ios_print( out_file, "%s%s", io->szTag2, strbuf->pStr );
            p = pOrigStruct->szAtoms;
            length = (int) strlen( p );
            io->tot_len = strbuf->nUsedLength;
            line_len = nMaxLineLen - io->tot_len;
            for (cur_pos = 0; cur_pos < length; cur_pos = last_pos)
            {
                if (length - cur_pos >= line_len)
                {
                    last_pos = cur_pos + line_len;
                    /* search backward for the nearest first atom letter (always uppercase) */
                    while (cur_pos < last_pos && !isupper( UCINT p[last_pos] ))
                    {
                        last_pos--;
                    }
                }
                else
                {
                    last_pos = length;
                }
                if (last_pos > cur_pos)
                {
                    memcpy(strbuf->pStr + strbuf->nUsedLength, p + cur_pos, (long long)last_pos - (long long)cur_pos); /* djb-rwth: cast operators added */
                    strbuf->pStr[strbuf->nUsedLength + last_pos - cur_pos] = '\0';
                    /*strbuf->nUsedLength = strbuf->nUsedLength + last_pos - cur_pos;*/

                    if (1) /* always show "Zy" as "Zz" */
                    {
                        char *pzy, *pstart=strbuf->pStr + strbuf->nUsedLength;
                        while ((pzy = strstr( pstart, "Zy" ))) /* djb-rwth: addressing LLVM warning */
                        {
                            *(++pzy) = 'z';
                            pstart = pzy;
                        }
                    }

                    inchi_ios_print( out_file, "%s%s", strbuf->pStr, io->bPlainTextTags ? "" : "\n" );
                }
                else
                {
                    break;
                }
            }
            if (pLF[0])
            {
                inchi_ios_print( out_file, "%s", pLF );
            }

            inchi_strbuf_reset( strbuf );

            /* Bonds /B: */
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_BNDR, io->szTag2, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf );
            inchi_ios_print( out_file, "%s%s", io->szTag2, strbuf->pStr );

            p = pOrigStruct->szBonds;
            length = (int) strlen( p );
            line_len = nMaxLineLen - io->tot_len;
            for (cur_pos = 0; cur_pos < length; cur_pos = last_pos)
            {
                if (length - cur_pos >= line_len)
                {
                    last_pos = cur_pos + line_len - 1;
                    /* search backward for the nearest first bond delimiter ";" */
                    while (cur_pos < last_pos && p[last_pos] != ';')
                    {
                        last_pos--;
                    }
                    if (cur_pos < last_pos)
                    {
                        last_pos++; /* include ';' at the end of the line */
                    }
                }
                else
                {
                    last_pos = length;
                }
                if (last_pos > cur_pos)
                {
                    memcpy(strbuf->pStr, p + cur_pos, (long long)last_pos - (long long)cur_pos); /* djb-rwth: cast operators added */
                    strbuf->pStr[last_pos - cur_pos] = '\0';
                    strbuf->nUsedLength = last_pos - cur_pos;
                    inchi_ios_print( out_file, "%s%s", strbuf->pStr, io->bPlainTextTags ? "" : "\n" );
                    inchi_strbuf_reset( strbuf );
                }
                else
                {
                    break;
                }
            }
            if (pLF[0])
            {
                inchi_ios_print( out_file, "%s", pLF );
            }

            /* Coordinates /C:    */
            szGetTag( AuxLbl, io->nTag, io->bTag2 = io->bTag1 | AL_XYZR, io->szTag2, &io->bAlways, 0 );
            inchi_strbuf_reset( strbuf );
            inchi_ios_print( out_file, "%s%s", io->szTag2, strbuf->pStr );

            p = pOrigStruct->szCoord;
            length = (int) strlen( p );
            line_len = nMaxLineLen - io->tot_len;
            for (cur_pos = 0; cur_pos < length; cur_pos = last_pos)
            {
                if (length - cur_pos >= line_len)
                {
                    last_pos = cur_pos + line_len - 1;
                    /* search backward for the nearest first coord. delimiter ";" */
                    while (cur_pos < last_pos && p[last_pos] != ';')
                    {
                        last_pos--;
                    }
                    if (cur_pos < last_pos)
                    {
                        last_pos++; /* include ';' at the end of the line */
                    }
                }
                else
                {
                    last_pos = length;
                }
                if (last_pos > cur_pos)
                {
                    memcpy(strbuf->pStr, p + cur_pos, (long long)last_pos - (long long)cur_pos); /* djb-rwth: cast operator added */
                    strbuf->pStr[last_pos - cur_pos] = '\0';
                    strbuf->nUsedLength = last_pos - cur_pos;
                    inchi_ios_print( out_file, "%s%s", strbuf->pStr, io->bPlainTextTags ? "" : "\n" );
                    inchi_strbuf_reset( strbuf );
                }
                else
                {
                    break;
                }
            }

            if (pLF[0])
            {
                inchi_ios_print( out_file, "%s", pLF );
            }
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: OutputAUXINFO_ReversibilityInfo

    let Some(structure) = original_structure else {
        return Ok(0);
    };
    if io.bSecondNonTautPass != 0
        || structure.num_atoms == 0
        || structure.szAtoms.is_null()
        || structure.szBonds.is_null()
        || structure.szCoord.is_null()
    {
        return Ok(0);
    }
    let max_line = 80_i32.min(string_buffer.nAllocatedLength);
    let text =
        |heap: &SourceHeap, ptr: SourceConstPointer<i8>| -> Result<Vec<i8>, SourceHeapError> {
            let bytes = heap.slice(ptr)?;
            let n = bytes
                .iter()
                .position(|b| *b == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            Ok(bytes[..n].to_vec())
        };
    let alloc = |heap: &mut SourceHeap, bytes: &[u8]| {
        Ok::<_, SourceHeapError>(
            heap.allocate_model_storage(bytes.iter().map(|b| *b as i8).collect())?
                .as_const(),
        )
    };
    let mut labels = vec![
        INCHI_TAG::default();
        crate::source_types::local_ichiprt1::tagAuxLblOrd_AL_MAX_ORD as usize
    ];
    for (i, plain, comment, xml) in [
        (
            0,
            b"/\0".as_slice(),
            b"fixed_H\0".as_slice(),
            b"fixed-H\0".as_slice(),
        ),
        (3, b"/\0", b"reversibility\0", b"reversibility\0"),
        (12, b"/rA:\0", b"atoms\0", b"atoms\0"),
        (13, b"/rB:\0", b"bonds\0", b"bonds\0"),
        (14, b"/rC:\0", b"xyz\0", b"xyz\0"),
    ] {
        labels[i] = INCHI_TAG {
            szPlainLabel: alloc(heap, plain)?,
            szPlainComment: alloc(heap, comment)?,
            szXmlLabel: alloc(heap, xml)?,
            bAlwaysOutput: 0,
        };
    }
    let tags = heap.allocate_model_storage(vec![0_i8; 128])?;
    io.bTag1 = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_REVR as i32 | io.bFhTag;
    szGetTag(heap, &labels, io.nTag, io.bTag1, tags, &mut io.bAlways, 0)?;
    io.szTag1
        .copy_from_slice(&heap.slice(tags.as_const())?[..64]);
    let lf = text(heap, line_feed)?;
    fn emit(
        heap: &mut SourceHeap,
        out: &mut Option<&mut INCHI_IOSTREAM>,
        stdout: SourceMutPointer<FILE>,
        bytes: &[i8],
    ) -> Result<(), SourceHeapError> {
        let s =
            heap.allocate_model_storage(bytes.iter().copied().chain(std::iter::once(0)).collect())?;
        let f = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
        inchi_ios_print(
            heap,
            out.as_deref_mut(),
            stdout,
            f.as_const(),
            &SourceVaList {
                arguments: vec![SourceFormatArgument::Bytes(s.as_const())],
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }
    for (field, bit, atom_mode) in [
        (
            structure.szAtoms.as_const(),
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_ATMR as i32,
            true,
        ),
        (
            structure.szBonds.as_const(),
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_BNDR as i32,
            false,
        ),
        (
            structure.szCoord.as_const(),
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_XYZR as i32,
            false,
        ),
    ] {
        io.bTag2 = io.bTag1 | bit;
        let tag2 = tags.offset(64)?;
        szGetTag(heap, &labels, io.nTag, io.bTag2, tag2, &mut io.bAlways, 0)?;
        io.szTag2
            .copy_from_slice(&heap.slice(tag2.as_const())?[..64]);
        inchi_strbuf_reset(heap, Some(string_buffer))?;
        emit(heap, &mut output, stdout, &text(heap, tag2.as_const())?)?;
        let source = text(heap, field)?;
        io.tot_len = string_buffer.nUsedLength;
        let line_len = max_line.wrapping_sub(io.tot_len);
        let mut cur = 0usize;
        while cur < source.len() {
            let mut last = source.len();
            if (source.len() - cur) as i32 >= line_len {
                last = cur
                    + usize::try_from(line_len)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                if !atom_mode {
                    last = last.saturating_sub(1);
                }
                while cur < last
                    && if atom_mode {
                        !source
                            .get(last)
                            .is_some_and(|byte| (*byte as u8).is_ascii_uppercase())
                    } else {
                        source[last] != b';' as i8
                    }
                {
                    last -= 1;
                }
                if !atom_mode && cur < last {
                    last += 1;
                }
            }
            if last <= cur {
                break;
            }
            let mut chunk = source[cur..last].to_vec();
            if atom_mode {
                for i in 0..chunk.len().saturating_sub(1) {
                    if chunk[i] == b'Z' as i8 && chunk[i + 1] == b'y' as i8 {
                        chunk[i + 1] = b'z' as i8;
                    }
                }
            }
            let n = chunk.len();
            let dst = heap.slice_mut(string_buffer.pStr)?;
            dst[..n].copy_from_slice(&chunk);
            dst[n] = 0;
            if !atom_mode {
                string_buffer.nUsedLength =
                    i32::try_from(n).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
            }
            emit(heap, &mut output, stdout, &chunk)?;
            if io.bPlainTextTags == 0 {
                emit(heap, &mut output, stdout, &[b'\n' as i8])?;
            }
            if !atom_mode {
                inchi_strbuf_reset(heap, Some(string_buffer))?;
            }
            cur = last;
        }
        if !lf.is_empty() {
            emit(heap, &mut output, stdout, &lf)?;
        }
        inchi_strbuf_reset(heap, Some(string_buffer))?;
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn OutputAUXINFO_PolymerInfo(
    heap: &mut SourceHeap,
    _canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    mut output: Option<&mut INCHI_IOSTREAM>,
    string_buffer: &mut INCHI_IOS_STRING,
    original_structure: Option<&ORIG_STRUCT>,
    _io: &mut INCHI_OUT_CTL,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4923 OutputAUXINFO_PolymerInfo
    // INCHI✔️❌: complete configured source frame follows verbatim.
    /*
    int OutputAUXINFO_PolymerInfo( CANON_GLOBALS    *pCG,
                                   INCHI_IOSTREAM   *out_file,
                                   INCHI_IOS_STRING *strbuf,
                                   ORIG_STRUCT      *pOrigStruct,
                                   INCHI_OUT_CTL    *io,
                                   char             *pLF,
                                   char             *pTAB )
    {
        int k, i, q;
        OAD_Polymer *p;
        OAD_PolymerUnit *u;


        if (!pOrigStruct)
        {
            return 0;
        }
        p = pOrigStruct->polymer;
        if (!p)
        {
            return 0;
        }

        inchi_strbuf_reset( strbuf );

        inchi_ios_print( out_file, "/Z:" );


        /* Print polymer units data */
        for (i = 0; i < p->n; i++)
        {
            /* For each unit u ... */
            u = p->units[i];

            /* print kinds of unit */
            inchi_strbuf_printf( strbuf, "%-d%-d%-d-", u->type, u->subtype, u->conn );
            inchi_strbuf_printf( strbuf, "%-s-", u->smt[0] ? u->smt : "n" );

            /* Print unit atoms */
            print_sequence_of_nums_compressing_ranges( u->na, u->alist, strbuf );

            /* Print bonds from unit to otside */
            if (u->nb > 0)
            {
                inchi_strbuf_printf( strbuf, "(" );
                for (k = 0; k < 2 * u->nb - 1; k++)
                {
                    inchi_strbuf_printf( strbuf, "%-d,", u->blist[k] );
                }
                inchi_strbuf_printf( strbuf, "%-d)", u->blist[2 * u->nb - 1] );
            }

            if (fabs( -fabs( u->xbr1[0] ) + 777777.777 ) > 1.e-7)
            {
                inchi_strbuf_printf( strbuf, "[" );
                for (q = 0; q < 3; q++)
                {
                    inchi_strbuf_printf( strbuf, "%-f,", u->xbr1[q] );
                }
                inchi_strbuf_printf( strbuf, "%-f]", u->xbr1[3] );
            }
            if (fabs( -fabs( u->xbr2[0] ) + 777777.777 ) > 1.e-7)
            {
                inchi_strbuf_printf( strbuf, "[" );
                for (q = 0; q < 3; q++)
                {
                    inchi_strbuf_printf( strbuf, "%-f,", u->xbr2[q] );
                }
                inchi_strbuf_printf( strbuf, "%-f]", u->xbr2[3] );
            }

            if (i < p->n - 1)
            {
                inchi_strbuf_printf( strbuf, ";" );
            }
        }

        inchi_ios_print( out_file, "%s%s", strbuf->pStr, pLF );

        return 0;
    }
    */
    // END INCHI C FUNCTION: OutputAUXINFO_PolymerInfo

    let Some(original_structure) = original_structure else {
        return Ok(0);
    };
    if original_structure.polymer.is_null() {
        return Ok(0);
    }
    let polymer = heap
        .slice(original_structure.polymer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();

    inchi_strbuf_reset(heap, Some(string_buffer))?;

    fn format_pointer(
        heap: &mut SourceHeap,
        format: &[u8],
    ) -> Result<SourceConstPointer<i8>, SourceHeapError> {
        Ok(heap
            .allocate_model_storage(format.iter().map(|byte| *byte as i8).collect())?
            .as_const())
    }
    fn append(
        heap: &mut SourceHeap,
        buffer: &mut INCHI_IOS_STRING,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format = format_pointer(heap, format)?;
        inchi_strbuf_printf(
            heap,
            Some(buffer),
            format,
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }
    fn print(
        heap: &mut SourceHeap,
        output: &mut Option<&mut INCHI_IOSTREAM>,
        stdout: SourceMutPointer<FILE>,
        format: &[u8],
        arguments: Vec<SourceFormatArgument>,
    ) -> Result<(), SourceHeapError> {
        let format = format_pointer(heap, format)?;
        inchi_ios_print(
            heap,
            output.as_deref_mut(),
            stdout,
            format,
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    }

    print(heap, &mut output, stdout, b"/Z:\0", Vec::new())?;

    for index in 0..polymer.n {
        let index = usize::try_from(index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let unit_pointer = *heap
            .slice(polymer.units.as_const())?
            .get(index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let unit = heap
            .slice(unit_pointer.as_const())?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();

        append(
            heap,
            string_buffer,
            b"%-d%-d%-d-\0",
            vec![
                SourceFormatArgument::Signed(i64::from(unit.type_)),
                SourceFormatArgument::Signed(i64::from(unit.subtype)),
                SourceFormatArgument::Signed(i64::from(unit.conn)),
            ],
        )?;
        let smt = if unit.smt[0] != 0 {
            let end = unit
                .smt
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            heap.allocate_model_storage(unit.smt[..=end].to_vec())?
        } else {
            heap.allocate_model_storage(vec![b'n' as i8, 0])?
        };
        append(
            heap,
            string_buffer,
            b"%-s-\0",
            vec![SourceFormatArgument::Bytes(smt.as_const())],
        )?;

        let atom_count =
            usize::try_from(unit.na).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if atom_count == 0 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let atoms = heap
            .slice(unit.alist.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec();
        crate::source::base::ichiprt2::print_sequence_of_nums_compressing_ranges(
            heap,
            &atoms,
            string_buffer,
        )?;

        if unit.nb > 0 {
            append(heap, string_buffer, b"(\0", Vec::new())?;
            let bond_value_count = unit
                .nb
                .checked_mul(2)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            let last = bond_value_count
                .checked_sub(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            for bond_index in 0..last {
                let value = *heap
                    .slice(unit.blist.as_const())?
                    .get(
                        usize::try_from(bond_index)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                    )
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                append(
                    heap,
                    string_buffer,
                    b"%-d,\0",
                    vec![SourceFormatArgument::Signed(i64::from(value))],
                )?;
            }
            let value = *heap
                .slice(unit.blist.as_const())?
                .get(usize::try_from(last).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            append(
                heap,
                string_buffer,
                b"%-d)\0",
                vec![SourceFormatArgument::Signed(i64::from(value))],
            )?;
        }

        for coordinates in [unit.xbr1, unit.xbr2] {
            if (-coordinates[0].abs() + 777_777.777).abs() > 1.0e-7 {
                append(heap, string_buffer, b"[\0", Vec::new())?;
                for value in coordinates.iter().take(3) {
                    append(
                        heap,
                        string_buffer,
                        b"%-f,\0",
                        vec![SourceFormatArgument::Float(*value)],
                    )?;
                }
                append(
                    heap,
                    string_buffer,
                    b"%-f]\0",
                    vec![SourceFormatArgument::Float(coordinates[3])],
                )?;
            }
        }

        if index + 1
            < usize::try_from(polymer.n).map_err(|_| SourceHeapError::PointerOutOfBounds)?
        {
            append(heap, string_buffer, b";\0", Vec::new())?;
        }
    }

    print(
        heap,
        &mut output,
        stdout,
        b"%s%s\0",
        vec![
            SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
            SourceFormatArgument::Bytes(line_feed),
        ],
    )?;
    Ok(0)
}

#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn OutputINCHI_PolymerLayer_SingleUnit(
    heap: &mut SourceHeap,
    unit: &mut crate::source_types::OAD_PolymerUnit,
    polymers_mode: i32,
    total_star_atoms: i32,
    used_stars: &mut i32,
    atom_properties: &[crate::source_types::OAD_AtProps],
    _canonical_numbers: &[i32],
    original_atom_data: &crate::source_types::ORIG_ATOM_DATA,
    original_structure: &crate::source_types::ORIG_STRUCT,
    string_buffer: &mut INCHI_IOS_STRING,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4111 OutputINCHI_PolymerLayer_SingleUnit
    // INCHI✔️❌: static int OutputINCHI_PolymerLayer_SingleUnit( OAD_PolymerUnit *u,
    // INCHI✔️❌:                                                 int bPolymers,
    // INCHI✔️❌:                                                 int total_star_atoms,
    // INCHI✔️❌:                                                 int *n_used_stars,
    // INCHI✔️❌:                                                 OAD_AtProps *aprops,
    // INCHI✔️❌:                                                 int *cano_nums,
    // INCHI✔️❌:                                                 ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                                                 ORIG_STRUCT *pOrigStruct,
    // INCHI✔️❌:                                                 INCHI_IOS_STRING *strbuf )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, k, tmp, a1 = 0, a2 = 0, a3 = 0, a4 = 0, b, curr_star_num;
    // INCHI✔️❌:     int err = 0;
    // INCHI✔️❌:     OAD_Polymer *p = orig_inp_data->polymer;
    // INCHI✔️❌:     inp_ATOM    *at = orig_inp_data->at;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* print unit type and subtype */
    // INCHI✔️❌:     inchi_strbuf_printf( strbuf, "%-d%-d%-d-", u->type, u->subtype, u->conn );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* print unit atoms */
    // INCHI✔️❌:     print_sequence_of_nums_compressing_ranges( u->na, u->alist, strbuf );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Print the crossing bonds or frame-shiftable pattern */
    // INCHI✔️❌:     if (u->nb > 2)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* not supported yet, too many bonds in SBL */
    // INCHI✔️❌:         err = 12;
    // INCHI✔️❌:         goto exit_function;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Print crossing bonds "(cap1-partner1,cap2-partner2)"    */
    // INCHI✔️❌:     if (u->nb == 2 && ( !u->cyclizable || !u->cyclized ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         int swap = 0;
    // INCHI✔️❌:         a1 = u->blist[0];
    // INCHI✔️❌:         a2 = u->blist[1];
    // INCHI✔️❌:         a3 = u->blist[2];
    // INCHI✔️❌:         a4 = u->blist[3];
    // INCHI✔️❌:
    // INCHI✔️❌:         if (is_in_the_ilist( u->alist, a1, u->na ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             tmp = a2;
    // INCHI✔️❌:             a2 = a1;
    // INCHI✔️❌:             a1 = tmp;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (is_in_the_ilist( u->alist, a3, u->na ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             tmp = a4;
    // INCHI✔️❌:             a4 = a3;
    // INCHI✔️❌:             a3 = tmp;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Always print first the crossing bond pointing to more senior CRU end ("head")    */
    // INCHI✔️❌:         if (bPolymers==POLYMERS_LEGACY)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* old, v. 1.05 */
    // INCHI✔️❌:             /* The first printed is the crossing bond with higher canonical number of the cap */
    // INCHI✔️❌:             swap = a3 < a1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* new in v. 1.06 */
    // INCHI✔️❌:             /* The first printed is the crossing bond pointing to more senior CRU end ("head")    */
    // INCHI✔️❌:             swap = (OAD_Polymer_IsFirstAtomRankLower(a2, a4, aprops) == 1);
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (swap)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "(%-d-%-d,%-d-%-d)", a3, a4, a1, a2 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "(%-d-%-d,%-d-%-d)", a1, a2, a3, a4 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     else if (u->nb <= 2 && ( u->cyclizable || u->nbkbonds > 0 ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Print frame-shiftable pattern "cap1,cap2-(b1a1,b1a2, b2a1,b2a2, ... )"     */
    // INCHI✔️❌:         /* where b1, b2, ... are CRU bonds potentially invilved in frame shift          */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*  Get actual star atoms numbers from all-units pool                           */
    // INCHI✔️❌:         /*  NB: ordered according to already established in units2/unum                 */
    // INCHI✔️❌:         if (u->cap1 > 0 || u->cap2 > 0)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             int n_expl_H = 0, pos = 0;
    // INCHI✔️❌:             char *sza = pOrigStruct->szAtoms;
    // INCHI✔️❌:             while (sza[pos])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 if (sza[pos] == 'H')
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (isupper( UCINT sza[pos + 1] ) || !sza[pos + 1])        /* if ( next_c is Uppercase or NUL ) */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         n_expl_H++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     if (!sza[pos + 1])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 pos++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (u->cap1 > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 curr_star_num = pOrigStruct->num_atoms - n_expl_H - total_star_atoms + *n_used_stars + 1;
    // INCHI✔️❌:                 if (curr_star_num > pOrigStruct->num_atoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     err = 11; goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 a1 = curr_star_num;
    // INCHI✔️❌:                 (*n_used_stars)++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (u->cap2 > 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 curr_star_num = pOrigStruct->num_atoms - n_expl_H - total_star_atoms + *n_used_stars + 1;
    // INCHI✔️❌:                 if (curr_star_num > pOrigStruct->num_atoms)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     err = 11; goto exit_function;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 a2 = curr_star_num;
    // INCHI✔️❌:                 (*n_used_stars)++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         /* a1 and a2 are number of star atoms associated (but actually  */
    // INCHI✔️❌:         /* disconnected at this moment ) with SRU head and tail atoms   */
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "(%-d,%-d-", a1, a2 );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (u->cyclizable == CLOSING_SRU_DIRADICAL)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%-d)", u->end_atom1 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (u->cyclizable == CLOSING_SRU_HIGHER_ORDER_BOND)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             a3 = u->end_atom1;
    // INCHI✔️❌:             a4 = u->end_atom2;
    // INCHI✔️❌:             inchi_sort_int_pair_ascending( &a3, &a4 );
    // INCHI✔️❌:             /* if ( a3 > a4 )   { tmp = a4; a4 = a3;  a3 = tmp;                  }*/
    // INCHI✔️❌:             inchi_strbuf_printf( strbuf, "%-d.%-d)", a3, a4 );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (u->cyclizable == CLOSING_SRU_RING)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (u->nbkbonds == 0)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* last resort */
    // INCHI✔️❌:                 a3 = u->end_atom1;
    // INCHI✔️❌:                 a4 = u->end_atom2;
    // INCHI✔️❌:                 inchi_sort_int_pair_ascending( &a3, &a4 );
    // INCHI✔️❌:                 /* if ( a3 > a4 ) { tmp = a4; a4 = a3; a3 = tmp; } */
    // INCHI✔️❌:                 inchi_strbuf_printf( strbuf, "%-d,%-d)", a3, a4 );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Sort all backbone bonds in min-at-number order */
    // INCHI✔️❌:                 for (b = 1; b < u->nbkbonds; b++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     int *tmp_psbond = u->bkbonds[b];
    // INCHI✔️❌:                     j = b - 1;
    // INCHI✔️❌:                     while (j >= 0 && IsBondAtomNumsLesser( u->bkbonds[j], tmp_psbond ) > 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         u->bkbonds[j + 1] = u->bkbonds[j];
    // INCHI✔️❌:                         j--;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     u->bkbonds[j + 1] = tmp_psbond;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 if (p->treat==POLYMERS_MODERN || p->treat==POLYMERS_LEGACY_PLUS)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* was #if DO_POLYMER_FRAME_SHIFT_AT_STRUCT_TO_INCHI_CONVERSION==1 */
    // INCHI✔️❌:                     /* was: OAD_PolymerUnit_ReorderPolymerFrameShiftLinks( u, orig_inp_data, aprops, cano_nums );
    // INCHI✔️❌:                     OAD_Polymer_DebugTrace(p);*/
    // INCHI✔️❌:
    // INCHI✔️❌:                     /*\tFind senior link and move it to the beginning of list
    // INCHI✔️❌:                         If necessary, swap atoms in link so that the first
    // INCHI✔️❌:                         points to SRU head and the next to SRU tail
    // INCHI✔️❌:
    // INCHI✔️❌:                         The net result of ordering is as follows:
    // INCHI✔️❌:                         at1,at2,  at3,at4,  at5,at6, ...
    // INCHI✔️❌:                         here
    // INCHI✔️❌:                         at1, at2 is the most senior bond, and at1 is more senior than at2
    // INCHI✔️❌:                         all other pairs at3,at4,  at5,at6, ... are sorted just in increasing
    // INCHI✔️❌:                         order of first number in pair, then second one, e.g.: at3<at4; at3<=at5
    // INCHI✔️❌:                         (and if at3==at5 then at4<at6)
    // INCHI✔️❌:                     */
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (p->frame_shift_scheme != FSS_NONE && u->nbkbonds >= 1 && u->cap1 >= 1 && u->cap2 >= 1)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (OAD_PolymerUnit_SetReopeningDetails(u, at))
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             /* Find senior backbone bond and move it to the beginning of list */
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 int senior_bond, bond0at1, bond0at2;
    // INCHI✔️❌:                                 OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(u, at, aprops, &senior_bond);
    // INCHI✔️❌:                                 if (senior_bond)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     bond0at1 = u->bkbonds[0][0];
    // INCHI✔️❌:                                     bond0at2 = u->bkbonds[0][1];
    // INCHI✔️❌:                                     u->bkbonds[0][0] = u->bkbonds[senior_bond][0];
    // INCHI✔️❌:                                     u->bkbonds[0][1] = u->bkbonds[senior_bond][1];
    // INCHI✔️❌:                                     u->end_atom1 = u->bkbonds[0][0];
    // INCHI✔️❌:                                     u->end_atom2 = u->bkbonds[0][1];
    // INCHI✔️❌:                                     u->bkbonds[senior_bond][0] = bond0at1;
    // INCHI✔️❌:                                     u->bkbonds[senior_bond][1] = bond0at2;
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* p->really_do_frame_shift = 0; */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 for (k = 0; k < u->nbkbonds; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     a3 = u->bkbonds[k][0]; a4 = u->bkbonds[k][1];
    // INCHI✔️❌:                     /*if ( a3 > a4 )    { tmp = a4; a4 = a3; a3 = tmp; }*/
    // INCHI✔️❌:                     inchi_strbuf_printf( strbuf, "%-d,%-d%-c", a3, a4, k == u->nbkbonds - 1 ? ')' : ',' );
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:
    // INCHI✔️❌:     return err;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_PolymerLayer_SingleUnit

    let write = |heap: &mut SourceHeap,
                 buffer: &mut INCHI_IOS_STRING,
                 format: &[u8],
                 arguments: Vec<SourceFormatArgument>|
     -> Result<(), SourceHeapError> {
        let format =
            heap.allocate_model_storage(format.iter().copied().map(|byte| byte as i8).collect())?;
        inchi_strbuf_printf(
            heap,
            Some(buffer),
            format.as_const(),
            &SourceVaList {
                arguments,
                ..SourceVaList::default()
            },
        )?;
        Ok(())
    };
    let signed = |value: i32| SourceFormatArgument::Signed(i64::from(value));

    let polymer = heap
        .slice(original_atom_data.polymer.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    write(
        heap,
        string_buffer,
        b"%-d%-d%-d-\0",
        vec![signed(unit.type_), signed(unit.subtype), signed(unit.conn)],
    )?;
    let atom_count =
        usize::try_from(unit.na).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    let unit_atoms = heap
        .slice(unit.alist.as_const())?
        .get(..atom_count)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .to_vec();
    crate::source::base::ichiprt2::print_sequence_of_nums_compressing_ranges(
        heap,
        &unit_atoms,
        string_buffer,
    )?;

    if unit.nb > 2 {
        return Ok(12);
    }

    let mut atom1 = 0_i32;
    let mut atom2 = 0_i32;
    if unit.nb == 2 && (unit.cyclizable == 0 || unit.cyclized == 0) {
        let crossing = heap.slice(unit.blist.as_const())?;
        let mut atom3 = *crossing
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut atom4 = *crossing.get(1).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut atom5 = *crossing.get(2).ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut atom6 = *crossing.get(3).ok_or(SourceHeapError::PointerOutOfBounds)?;
        if crate::source::base::util::is_in_the_ilist(Some(&unit_atoms), atom3, unit.na)?.is_some()
        {
            std::mem::swap(&mut atom3, &mut atom4);
        }
        if crate::source::base::util::is_in_the_ilist(Some(&unit_atoms), atom5, unit.na)?.is_some()
        {
            std::mem::swap(&mut atom5, &mut atom6);
        }
        let swap = if polymers_mode == crate::source_types::POLYMERS_LEGACY as i32 {
            atom5 < atom3
        } else {
            crate::source::base::runichi3::OAD_Polymer_IsFirstAtomRankLower(
                atom4,
                atom6,
                atom_properties,
            )? == 1
        };
        let values = if swap {
            [atom5, atom6, atom3, atom4]
        } else {
            [atom3, atom4, atom5, atom6]
        };
        write(
            heap,
            string_buffer,
            b"(%-d-%-d,%-d-%-d)\0",
            values.into_iter().map(signed).collect(),
        )?;
    } else if unit.nb <= 2 && (unit.cyclizable != 0 || unit.nbkbonds > 0) {
        if unit.cap1 > 0 || unit.cap2 > 0 {
            let atom_text = heap.slice(original_structure.szAtoms.as_const())?;
            let terminator = atom_text
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            let mut explicit_hydrogens = 0_i32;
            let mut position = 0_usize;
            while position < terminator {
                if atom_text[position] == b'H' as i8 {
                    let next = *atom_text
                        .get(position + 1)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?
                        as u8;
                    if next.is_ascii_uppercase() || next == 0 {
                        explicit_hydrogens = explicit_hydrogens
                            .checked_add(1)
                            .ok_or(SourceHeapError::SourceIntegerOverflow)?;
                    }
                    if next == 0 {
                        break;
                    }
                }
                position += 1;
            }
            let next_star = |used: i32| -> Result<i32, SourceHeapError> {
                original_structure
                    .num_atoms
                    .checked_sub(explicit_hydrogens)
                    .and_then(|value| value.checked_sub(total_star_atoms))
                    .and_then(|value| value.checked_add(used))
                    .and_then(|value| value.checked_add(1))
                    .ok_or(SourceHeapError::SourceIntegerOverflow)
            };
            if unit.cap1 > 0 {
                atom1 = next_star(*used_stars)?;
                if atom1 > original_structure.num_atoms {
                    return Ok(11);
                }
                *used_stars = used_stars
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
            if unit.cap2 > 0 {
                atom2 = next_star(*used_stars)?;
                if atom2 > original_structure.num_atoms {
                    return Ok(11);
                }
                *used_stars = used_stars
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }
        write(
            heap,
            string_buffer,
            b"(%-d,%-d-\0",
            vec![signed(atom1), signed(atom2)],
        )?;

        if unit.cyclizable == crate::source_types::CLOSING_SRU_DIRADICAL as i32 {
            write(heap, string_buffer, b"%-d)\0", vec![signed(unit.end_atom1)])?;
        } else if unit.cyclizable == crate::source_types::CLOSING_SRU_HIGHER_ORDER_BOND as i32 {
            let (mut atom3, mut atom4) = (unit.end_atom1, unit.end_atom2);
            inchi_sort_int_pair_ascending(&mut atom3, &mut atom4);
            write(
                heap,
                string_buffer,
                b"%-d.%-d)\0",
                vec![signed(atom3), signed(atom4)],
            )?;
        } else if unit.cyclizable == crate::source_types::CLOSING_SRU_RING as i32 {
            if unit.nbkbonds == 0 {
                let (mut atom3, mut atom4) = (unit.end_atom1, unit.end_atom2);
                inchi_sort_int_pair_ascending(&mut atom3, &mut atom4);
                write(
                    heap,
                    string_buffer,
                    b"%-d,%-d)\0",
                    vec![signed(atom3), signed(atom4)],
                )?;
            } else {
                for bond in 1..unit.nbkbonds {
                    let bond_index = usize::try_from(bond)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    let temporary = *heap
                        .slice(unit.bkbonds.as_const())?
                        .get(bond_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let mut previous = bond - 1;
                    while previous >= 0 {
                        let previous_index = usize::try_from(previous)
                            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                        let current = *heap
                            .slice(unit.bkbonds.as_const())?
                            .get(previous_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        let current_row = heap.slice(current.as_const())?;
                        let temporary_row = heap.slice(temporary.as_const())?;
                        let current_pair = [current_row[0], current_row[1]];
                        let temporary_pair = [temporary_row[0], temporary_row[1]];
                        if IsBondAtomNumsLesser(&current_pair, &temporary_pair) <= 0 {
                            break;
                        }
                        heap.slice_mut(unit.bkbonds)?[previous_index + 1] = current;
                        previous -= 1;
                    }
                    let destination = usize::try_from(previous + 1)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    heap.slice_mut(unit.bkbonds)?[destination] = temporary;
                }

                if polymer.treat == crate::source_types::POLYMERS_MODERN as i32
                    || polymer.treat == crate::source_types::POLYMERS_LEGACY_PLUS as i32
                {
                    if polymer.frame_shift_scheme
                        != crate::source_types::tagFrameShifScheme_FSS_NONE as i32
                        && unit.nbkbonds >= 1
                        && unit.cap1 >= 1
                        && unit.cap2 >= 1
                    {
                        let atoms = heap.slice(original_atom_data.at.as_const())?;
                        if crate::source::base::runichi3::OAD_PolymerUnit_SetReopeningDetails(
                            heap, unit, atoms,
                        )? != 0
                        {
                            let mut senior_bond = 0_i32;
                            crate::source::base::runichi3::OAD_PolymerUnit_SortBackboneBondsAndSetSeniors(
                                heap,
                                unit,
                                original_atom_data.at,
                                atom_properties,
                                &mut senior_bond,
                            )?;
                            if senior_bond != 0 {
                                let rows = heap.slice(unit.bkbonds.as_const())?;
                                let first_pointer = rows[0];
                                let senior_pointer = rows[usize::try_from(senior_bond)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?];
                                let first = {
                                    let row = heap.slice(first_pointer.as_const())?;
                                    [row[0], row[1]]
                                };
                                let senior = {
                                    let row = heap.slice(senior_pointer.as_const())?;
                                    [row[0], row[1]]
                                };
                                heap.slice_mut(first_pointer)?[..2].copy_from_slice(&senior);
                                unit.end_atom1 = senior[0];
                                unit.end_atom2 = senior[1];
                                heap.slice_mut(senior_pointer)?[..2].copy_from_slice(&first);
                            }
                        }
                    }
                }

                for index in 0..unit.nbkbonds {
                    let row_pointer = *heap
                        .slice(unit.bkbonds.as_const())?
                        .get(
                            usize::try_from(index)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let row = heap.slice(row_pointer.as_const())?;
                    let delimiter = if index == unit.nbkbonds - 1 {
                        b')'
                    } else {
                        b','
                    };
                    write(
                        heap,
                        string_buffer,
                        b"%-d,%-d%-c\0",
                        vec![
                            signed(row[0]),
                            signed(row[1]),
                            SourceFormatArgument::Byte(delimiter as i8),
                        ],
                    )?;
                }
            }
        }
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn IsBondAtomNumsLesser(bond1: &[i32; 2], bond2: &[i32; 2]) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5007 IsBondAtomNumsLesser
    // INCHI✔️✔️: int IsBondAtomNumsLesser( int *bond1, int* bond2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int min1 = inchi_min( bond1[0], bond1[1] );
    // INCHI✔️✔️:     int min2 = inchi_min( bond2[0], bond2[1] );
    // INCHI✔️✔️:     int max1 = inchi_max( bond1[0], bond1[1] );
    // INCHI✔️✔️:     int max2 = inchi_max( bond2[0], bond2[1] );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     if (min1 < min2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return -1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (min1 > min2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         return 1;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:     if (min1 == min2)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         if (max1 < max2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return -1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:         if (max1 > max2)
    // INCHI✔️✔️:         {
    // INCHI✔️✔️:             return 1;
    // INCHI✔️✔️:         }
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return 0;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: IsBondAtomNumsLesser

    let minimum1 = bond1[0].min(bond1[1]);
    let minimum2 = bond2[0].min(bond2[1]);
    let maximum1 = bond1[0].max(bond1[1]);
    let maximum2 = bond2[0].max(bond2[1]);
    if minimum1 < minimum2 {
        return -1;
    }
    if minimum1 > minimum2 {
        return 1;
    }
    if maximum1 < maximum2 {
        return -1;
    }
    if maximum1 > maximum2 {
        return 1;
    }
    0
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub(crate) fn EditINCHI_HidePolymerZz(
    heap: &mut SourceHeap,
    output: &mut INCHI_IOSTREAM,
    number_of_polymer_zz: i32,
    number_of_placeholder_zy: i32,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5039 EditINCHI_HidePolymerZz
    // INCHI✔️❌: complete source frame follows verbatim.
    /*
void EditINCHI_HidePolymerZz(INCHI_IOSTREAM *out, int n_pzz, int n_zy)
{
    char *s = out->s.pStr, *s0, *buf = NULL;
    char prev_layer_symbol = '0';
    int i, j, ii, nzz,
        nzz1 = 0, nslash = 0, ncopied = 0,
        start = 0, skip = 0, is_in_z_layer = 0,
        eol_was_consumed = 0, pre_eol = 0,
        nonprt_sym = 0, nonprt_prev = 0;

    if (n_zy > 0) 
    {
        /* We have some placeholder pseudo atoms which should not be removed below (if anyway they are allowed) */
        if (n_pzz == 0)
        {
            /* Have nothing to remove, just exit */
            return;
        }
        /* Have both polymer-related and placeholder pseudo atoms */
        if (n_pzz < 2)
        {
            /* Something strange, should not arrive here, cowardly exit */
            return;
        }
    }

    /* Ensure that polymeric layer is present */
    if (!strstr(s, "/z"))
    {
        return;
    }
    s0 = strstr(s, "InChI=1B/");
    if (!s0)
    {
        return;
    }

#if 0
    nzz1 = CountPseudoElementInFormula("Zz", s0 + strlen("InChI=1B/"));
    if (nzz1 == 0)
    {
        return;
    }
    if (nzz1 != (n_pzz + n_zy))
    {
        /* Something strange, should not arrive here, cowardly exit */
        return;
    }
#endif
    nzz1 = n_pzz;

    /* OK, we must hide n_pzz Zz's*/
    buf = (char *) inchi_calloc( (long long)out->s.nUsedLength + 1, sizeof( char ) ); /* djb-rwth: cast operator added */
    if (!buf)
    {
        return;
    }

    /* Consume '\n' temporarily */
    if (s[out->s.nUsedLength - 1] == '\n')
    {
        s[out->s.nUsedLength - 1] = '\0';
        out->s.nUsedLength--;
        eol_was_consumed = 1;
    }

    start = s0 - s;
    nzz = nzz1;
    is_in_z_layer = skip = 0;
    for (i = start; i < out->s.nUsedLength; i++)
    {
        pre_eol = (i == out->s.nUsedLength - 1);
        nonprt_sym = s[i] == '\n' || s[i] == '\r' || s[i] == '\t';

        if (!skip)
        {
            buf[ncopied] = s[i];
            ncopied++;
            if (nonprt_sym && nonprt_prev)
            {
                continue;
            }
        }
        nonprt_prev = nonprt_sym;

        if (is_in_z_layer && !skip)
        {
            if (s[i] == '(')
            {
                /* Software version 1.07 : skip pattern "(cap,cap-bkbonds)" but not "(cap-end, cap-end)" */
                const char *q;
                const char *p = out->s.pStr + i + 1;
                AT_NUMB ia = (AT_NUMB) inchi_strtol(p, &q, 10); /* make compiler happy: */ /* djb-rwth: removing redundant code; ignoring LLVM warning: variable used to store function return value */
                if (*q != '-')
                {
                    skip = 1; 
                }
            }
        }
        else if (is_in_z_layer && skip)
        {
            if (s[i] == '-')
            {
                skip = 0;
            }
        }

        if (s[i] == '/' || pre_eol || nonprt_sym )
        {
            if (is_in_z_layer)
            {
                is_in_z_layer = 0;
            }

            if (s[i] == '/')
            {
                nslash++;
            }

            if (nslash == 2 ||
                ( nslash == 1 && pre_eol ) ||
                ( prev_layer_symbol == 'f' )
                )
            {
                if (nzz)
                {
                    /* eat Zz's */
                    ii = i;
                    if (pre_eol)
                    {
                        ii = i + 1;
                    }
                    if (s[ii - 1] == 'z' && s[ii - 2] == 'Z')
                    {
                        ncopied -= 2;
                        for (j = ii - 3; j >= 0; j--)
                        {
                            if (s[j] == '.')
                            {
                                break;
                            }
                            ncopied--;
                        }
                        ncopied--;
                        if (!pre_eol)
                        {
                            buf[ncopied - 1] = '/';
                        }
                        else
                        {
                            buf[ncopied - 1] = '\0';
                        }
                    }
                }
            }
            else if (nslash > 2 || pre_eol || s[i] == '\n')
            {
                if (nzz)
                {
                    if (prev_layer_symbol != 'p' &&
                         prev_layer_symbol != 's' &&
                         prev_layer_symbol != 'f' &&
                         prev_layer_symbol != 'z'
                        )
                    {
                        /* eat nzz last ; if any */
                        int n_eaten = 0;
                        char eatable = ';';
                        if (prev_layer_symbol == 'm')
                        {
                            eatable = '.';
                        }
                        ii = i;
                        if (pre_eol)
                        {
                            ii = i + 1;
                        }
                        for (j = ii - 1; j >= 0; j--)
                        {
                            if (s[j] == eatable  && n_eaten < nzz)
                            {
                                ncopied--;
                                n_eaten++;
                            }
                            else
                            {
                                break;
                            }
                        }
                        if (!pre_eol)
                        {
                            if (s[i] == '\n')
                            {
                                buf[ncopied - 1] = '\n';
                            }
                            else
                            {
                                buf[ncopied - 1] = '/';
                            }
                        }
                        else
                        {
                            buf[ncopied] = '\0';
                            break;
                        }
                    }
                }
            }

            prev_layer_symbol = s[i + 1];
            if (s[i + 1] == 'r')
            {
                if (i + 3 < out->s.nUsedLength && s[i + 3] == ':')
                {
                    /* ra: rB; rC: in AuxInfo */
                    ;
                }
                else
                {
                    /* reconnected layer */
                    nslash = 1;
                    /* nzz = nzz2;    paranoidal */
                }
            }
            else if (s[i + 1] == 'z')
            {
                is_in_z_layer = 1;
            }
        }
    }

    out->s.nUsedLength = 0;
    inchi_ios_print_nodisplay( out, "%s%s", buf, eol_was_consumed ? "\n" : "" );
    inchi_free( buf );

    return;
}
    */
    // END INCHI C FUNCTION: EditINCHI_HidePolymerZz
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: EditINCHI_HidePolymerZz
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; the literal #if 0 block is inactive.
    // INCHI✔️❌: Rust uses checked model pointers for C undefined pointer/index domains and retains the
    // INCHI✔️❌: source allocation, in-place newline mutation, output reset, print, and cleanup ordering.
    // END INCHI ACTIVE MACRO CONFIGURATION: EditINCHI_HidePolymerZz

    if number_of_placeholder_zy > 0 {
        if number_of_polymer_zz == 0 {
            return Ok(());
        }
        if number_of_polymer_zz < 2 {
            return Ok(());
        }
    }

    if output.s.pStr.is_null() {
        return Err(SourceHeapError::NullPointer);
    }
    let source = heap.slice(output.s.pStr.as_const())?.to_vec();
    let nul = source
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let visible = &source[..nul];
    let find = |needle: &[u8]| {
        visible
            .windows(needle.len())
            .position(|window| window.iter().map(|byte| *byte as u8).eq(needle.iter().copied()))
    };
    if find(b"/z").is_none() {
        return Ok(());
    }
    let Some(start) = find(b"InChI=1B/") else {
        return Ok(());
    };

    let allocation_count = i64::from(output.s.nUsedLength) + 1;
    let buffer = match inchi_calloc::<i8>(heap, allocation_count as u64, 1) {
        Ok(pointer) => pointer,
        Err(SourceHeapError::AllocationFailed)
        | Err(SourceHeapError::AllocationSizeOverflow)
        | Err(SourceHeapError::AllocationElementCountOutOfRange) => return Ok(()),
        Err(error) => return Err(error),
    };

    let operation = (|| -> Result<(), SourceHeapError> {
        let mut input = source;
        let mut used = usize::try_from(output.s.nUsedLength)
            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
        if used == 0 || used > input.len() {
            return Err(SourceHeapError::UnsupportedSourceBehavior);
        }
        let mut consumed_eol = false;
        if input[used - 1] as u8 == b'\n' {
            input[used - 1] = 0;
            heap.slice_mut(output.s.pStr)?[used - 1] = 0;
            output.s.nUsedLength = output.s.nUsedLength.wrapping_sub(1);
            used -= 1;
            consumed_eol = true;
        }

        let mut copied = vec![0_i8; usize::try_from(allocation_count)
            .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?];
        let mut copied_count = 0_i32;
        let mut slash_count = 0_i32;
        let mut skip = false;
        let mut in_z_layer = false;
        let mut previous_nonprinting = false;
        let mut previous_layer_symbol = b'0';
        let polymer_zz = number_of_polymer_zz;
        let mut index = start;
        while index < used {
            let pre_eol = index == used - 1;
            let character = *input
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
            let nonprinting = matches!(character, b'\n' | b'\r' | b'\t');

            if !skip {
                let target = usize::try_from(copied_count)
                    .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                *copied
                    .get_mut(target)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? = character as i8;
                copied_count = copied_count.wrapping_add(1);
                if nonprinting && previous_nonprinting {
                    index += 1;
                    continue;
                }
            }
            previous_nonprinting = nonprinting;

            if in_z_layer && !skip {
                if character == b'(' {
                    let parse_pointer = output
                        .s
                        .pStr
                        .as_const()
                        .offset(i64::try_from(index + 1).map_err(|_| SourceHeapError::PointerOffsetOverflow)?)?;
                    let mut end = SourceConstPointer::null();
                    let _atom_number = inchi_strtol(heap, parse_pointer, Some(&mut end), 10)? as u16;
                    let end_character = *heap
                        .slice(end)?
                        .first()
                        .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
                    if end_character != b'-' {
                        skip = true;
                    }
                }
            } else if in_z_layer && skip && character == b'-' {
                skip = false;
            }

            if character == b'/' || pre_eol || nonprinting {
                if in_z_layer {
                    in_z_layer = false;
                }
                if character == b'/' {
                    slash_count = slash_count.wrapping_add(1);
                }

                if slash_count == 2
                    || (slash_count == 1 && pre_eol)
                    || previous_layer_symbol == b'f'
                {
                    if polymer_zz != 0 {
                        let boundary = if pre_eol { index + 1 } else { index };
                        if boundary < 2 {
                            return Err(SourceHeapError::UnsupportedSourceBehavior);
                        }
                        if input[boundary - 1] as u8 == b'z'
                            && input[boundary - 2] as u8 == b'Z'
                        {
                            copied_count = copied_count.wrapping_sub(2);
                            let mut scan = i64::try_from(boundary)
                                .map_err(|_| SourceHeapError::PointerOffsetOverflow)? - 3;
                            while scan >= 0 {
                                if input[scan as usize] as u8 == b'.' {
                                    break;
                                }
                                copied_count = copied_count.wrapping_sub(1);
                                scan -= 1;
                            }
                            copied_count = copied_count.wrapping_sub(1);
                            let target = if pre_eol {
                                copied_count.wrapping_sub(1)
                            } else {
                                copied_count.wrapping_sub(1)
                            };
                            let target = usize::try_from(target)
                                .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                            copied[target] = if pre_eol { 0 } else { b'/' as i8 };
                        }
                    }
                } else if slash_count > 2 || pre_eol || character == b'\n' {
                    if polymer_zz != 0
                        && !matches!(previous_layer_symbol, b'p' | b's' | b'f' | b'z')
                    {
                        let eatable = if previous_layer_symbol == b'm' { b'.' } else { b';' };
                        let boundary = if pre_eol { index + 1 } else { index };
                        let mut scan = i64::try_from(boundary)
                            .map_err(|_| SourceHeapError::PointerOffsetOverflow)? - 1;
                        let mut eaten = 0_i32;
                        while scan >= 0
                            && input[scan as usize] as u8 == eatable
                            && eaten < polymer_zz
                        {
                            copied_count = copied_count.wrapping_sub(1);
                            eaten = eaten.wrapping_add(1);
                            scan -= 1;
                        }
                        if !pre_eol {
                            let target = usize::try_from(copied_count.wrapping_sub(1))
                                .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                            copied[target] = if character == b'\n' { b'\n' as i8 } else { b'/' as i8 };
                        } else {
                            let target = usize::try_from(copied_count)
                                .map_err(|_| SourceHeapError::UnsupportedSourceBehavior)?;
                            copied[target] = 0;
                            break;
                        }
                    }
                }

                let next = *input
                    .get(index + 1)
                    .ok_or(SourceHeapError::PointerOutOfBounds)? as u8;
                previous_layer_symbol = next;
                if next == b'r' {
                    if !(index + 3 < used && input[index + 3] as u8 == b':') {
                        slash_count = 1;
                    }
                } else if next == b'z' {
                    in_z_layer = true;
                }
            }
            index += 1;
        }

        heap.slice_mut(buffer)?.copy_from_slice(&copied);
        output.s.nUsedLength = 0;
        let format = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, b'%' as i8, b's' as i8, 0])?;
        let ending = heap.allocate_model_storage(if consumed_eol {
            vec![b'\n' as i8, 0]
        } else {
            vec![0]
        })?;
        let print_result = inchi_ios_print_nodisplay(
            heap,
            Some(output),
            SourceMutPointer::null(),
            format.as_const(),
            &SourceVaList {
                arguments: vec![
                    SourceFormatArgument::Bytes(buffer.as_const()),
                    SourceFormatArgument::Bytes(ending.as_const()),
                ],
                ..SourceVaList::default()
            },
        );
        let ending_cleanup = heap.free(ending);
        let format_cleanup = heap.free(format);
        match (print_result, ending_cleanup, format_cleanup) {
            (Err(error), _, _)
            | (Ok(_), Err(error), _)
            | (Ok(_), Ok(()), Err(error)) => Err(error),
            (Ok(_), Ok(()), Ok(())) => Ok(()),
        }
    })();
    let cleanup = inchi_free(heap, buffer);
    match (operation, cleanup) {
        (Err(error), _) | (Ok(()), Err(error)) => Err(error),
        (Ok(()), Ok(())) => Ok(()),
    }
}

pub(crate) fn inchi_sort_int_pair_ascending(a: &mut i32, b: &mut i32) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5525 inchi_sort_int_pair_ascending
    // INCHI✔️✔️: static void inchi_sort_int_pair_ascending(int *a, int *b)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int tmp;
    // INCHI✔️✔️:     if (*a > *b)
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         tmp = *b;
    // INCHI✔️✔️:         *b = *a;
    // INCHI✔️✔️:         *a = tmp;
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: inchi_sort_int_pair_ascending

    if *a > *b {
        let temporary = *b;
        *b = *a;
        *a = temporary;
    }
}

#[allow(non_snake_case)]
pub(crate) fn GetSaveOptLetters(save_option_bits: u8, letter1: &mut i8, letter2: &mut i8) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3100 GetSaveOptLetters
    // INCHI✔️✔️: void GetSaveOptLetters( unsigned char save_opt_bits, char* let1, char* let2 )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     const char a2p[] = "ABCDEFGHIJKLMNOP";
    // INCHI✔️✔️:     /* SaveOptBits layout: {unused|unused|Ket|15T|RecMet|FixedH|SUU|SLUUD} */
    // INCHI✔️✔️:     *let1 = a2p[(size_t) ( save_opt_bits & 0x0f )];
    // INCHI✔️✔️:     *let2 = a2p[(size_t) ( ( save_opt_bits & 0x30 ) >> 4 )];
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: GetSaveOptLetters

    const LETTERS: &[u8; 16] = b"ABCDEFGHIJKLMNOP";
    *letter1 = LETTERS[usize::from(save_option_bits & 0x0f)] as i8;
    *letter2 = LETTERS[usize::from((save_option_bits & 0x30) >> 4)] as i8;
}

pub(crate) fn set_line_separators(
    output_options: i32,
    line_feed: &mut &'static [i8],
    tab: &mut &'static [i8],
) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3112 set_line_separators
    // INCHI✔️✔️: void set_line_separators( int bINChIOutputOptions, char **pLF, char **pTAB )
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     int  bPlainTextCommnts = 0 != ( bINChIOutputOptions & INCHI_OUT_PLAIN_TEXT_COMMENTS );
    // INCHI✔️✔️:
    // INCHI✔️✔️:     *pLF = (char *)(bPlainTextCommnts ? "\n" : "\0");
    // INCHI✔️✔️:
    // INCHI✔️✔️: #if ( !defined(TARGET_API_LIB) && !defined(TARGET_LIB_FOR_WINCHI) )
    // INCHI✔️✔️:     {
    // INCHI✔️✔️:         int  bPlainText = 0 != ( bINChIOutputOptions & ( INCHI_OUT_PLAIN_TEXT | INCHI_OUT_PLAIN_TEXT_COMMENTS ) );
    // INCHI✔️✔️:         int  bPlainTabbedOutput = 0 != ( bINChIOutputOptions & INCHI_OUT_TABBED_OUTPUT ) &&
    // INCHI✔️✔️:             bPlainText && !bPlainTextCommnts;
    // INCHI✔️✔️:
    // INCHI✔️✔️:         *pTAB = bPlainTabbedOutput ? "\t" : "\n";
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #else
    // INCHI✔️✔️:     *pTAB = "\n";
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️:
    // INCHI✔️✔️:     return;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: set_line_separators

    const EMPTY: &[i8] = &[0];
    const NEWLINE: &[i8] = &[b'\n' as i8, 0];
    *line_feed = if output_options & INCHI_OUT_PLAIN_TEXT_COMMENTS as i32 != 0 {
        NEWLINE
    } else {
        EMPTY
    };
    *tab = NEWLINE;
}

#[allow(non_snake_case)]
pub(crate) fn OutputINCHI_VersionAndKind(
    heap: &mut SourceHeap,
    output: &mut INCHI_IOSTREAM,
    string_buffer: &mut INCHI_IOS_STRING,
    output_options: i32,
    is_beta: i32,
    line_feed: SourceConstPointer<i8>,
    _tab: SourceConstPointer<i8>,
    stdout: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3137 OutputINCHI_VersionAndKind
    // INCHI✔️❌: int OutputINCHI_VersionAndKind( INCHI_IOSTREAM   *out_file,
    // INCHI✔️❌:                                 INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                 int              bINChIOutputOptions,
    // INCHI✔️❌:                                 int              is_beta,
    // INCHI✔️❌:                                 char             *pLF,
    // INCHI✔️❌:                                 char             *pTAB )
    // INCHI✔️❌: {
    // INCHI✔️❌:     inchi_ios_print_nodisplay( out_file, "%s%s=%s", pLF, INCHI_NAME, pLF );
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:     inchi_strbuf_printf( strbuf, "%s", x_curr_ver );
    // INCHI✔️❌:
    // INCHI✔️❌:     /* - add 'Beta' flag if applicable */
    // INCHI✔️❌:     if (is_beta)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "B" );
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* - add 'Standard' flag if applicable */
    // INCHI✔️❌:     else if (bINChIOutputOptions & INCHI_OUT_STDINCHI)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         inchi_strbuf_printf( strbuf, "S" );
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_ios_print_nodisplay( out_file, "%s%s", strbuf->pStr, pLF );
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OutputINCHI_VersionAndKind

    let format_header =
        heap.allocate_model_storage(b"%s%s=%s\0".iter().map(|byte| *byte as i8).collect())?;
    let format_pair =
        heap.allocate_model_storage(b"%s%s\0".iter().map(|byte| *byte as i8).collect())?;
    let format_string = heap.allocate_model_storage(vec![b'%' as i8, b's' as i8, 0])?;
    let beta = heap.allocate_model_storage(vec![b'B' as i8, 0])?;
    let standard = heap.allocate_model_storage(vec![b'S' as i8, 0])?;
    let name = heap.allocate_model_storage(
        crate::source_types::INCHI_NAME
            .iter()
            .map(|byte| *byte as i8)
            .collect(),
    )?;
    let version = heap.allocate_model_storage(
        crate::source_types::local_ichiprt1::x_curr_ver
            .iter()
            .map(|byte| *byte as i8)
            .collect(),
    )?;

    inchi_ios_print_nodisplay(
        heap,
        Some(output),
        stdout,
        format_header.as_const(),
        &SourceVaList {
            arguments: vec![
                SourceFormatArgument::Bytes(line_feed),
                SourceFormatArgument::Bytes(name.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
            ..SourceVaList::default()
        },
    )?;
    inchi_strbuf_reset(heap, Some(string_buffer))?;
    inchi_strbuf_printf(
        heap,
        Some(string_buffer),
        format_string.as_const(),
        &SourceVaList {
            arguments: vec![SourceFormatArgument::Bytes(version.as_const())],
            ..SourceVaList::default()
        },
    )?;
    if is_beta != 0 {
        inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            beta.as_const(),
            &SourceVaList::default(),
        )?;
    } else if output_options & INCHI_OUT_STDINCHI as i32 != 0 {
        inchi_strbuf_printf(
            heap,
            Some(string_buffer),
            standard.as_const(),
            &SourceVaList::default(),
        )?;
    }
    inchi_ios_print_nodisplay(
        heap,
        Some(output),
        stdout,
        format_pair.as_const(),
        &SourceVaList {
            arguments: vec![
                SourceFormatArgument::Bytes(string_buffer.pStr.as_const()),
                SourceFormatArgument::Bytes(line_feed),
            ],
            ..SourceVaList::default()
        },
    )?;
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn InternallyGetCanoNumsAndComponentNums(
    heap: &mut SourceHeap,
    canonical_globals: SourceMutPointer<CANON_GLOBALS>,
    string_buffer: &mut INCHI_IOS_STRING,
    output_control: &mut INCHI_OUT_CTL,
    number_of_atoms: i32,
    canonical_numbers: SourceMutPointer<i32>,
    component_numbers: SourceMutPointer<i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5362 InternallyGetCanoNumsAndComponentNums
    // INCHI✔️❌: int InternallyGetCanoNumsAndComponentNums( CANON_GLOBALS    *pCG,
    // INCHI✔️❌:                                            INCHI_IOS_STRING *strbuf,
    // INCHI✔️❌:                                            INCHI_OUT_CTL    *io,
    // INCHI✔️❌:                                            int              nat,
    // INCHI✔️❌:                                            int              *cano_nums,
    // INCHI✔️❌:                                            int              *compnt_nums )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int orig_num, cano_num, icompnt, i, k, ndigit, err;
    // INCHI✔️❌:     char c, cnum[8];
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!cano_nums || !compnt_nums || !strbuf->pStr)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:     io->tot_len = str_AuxNumb( pCG, io->pINChISort, io->pINChISort2,
    // INCHI✔️❌:                                strbuf, &io->bOverflow, io->bOutType,
    // INCHI✔️❌:                                io->TAUT_MODE, io->num_components,
    // INCHI✔️❌:                                io->bSecondNonTautPass, io->bOmitRepetitions );
    // INCHI✔️❌:     for (i = 0; i < nat; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         compnt_nums[i] = -1;
    // INCHI✔️❌:         cano_nums[i + 1] = -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     ndigit = 0;
    // INCHI✔️❌:     err = 0;
    // INCHI✔️❌:     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:     icompnt = 1;
    // INCHI✔️❌:     cano_num = 0;
    // INCHI✔️❌:     for (k = 0; k <= strbuf->nUsedLength; k++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         c = strbuf->pStr[k];
    // INCHI✔️❌:         if (c == ',' || c == ';' || c == '\0')
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cnum[ndigit] = '\0';
    // INCHI✔️❌:             orig_num = atoi( cnum );
    // INCHI✔️❌:             cano_nums[orig_num] = cano_num;
    // INCHI✔️❌:             compnt_nums[cano_num] = icompnt;
    // INCHI✔️❌:             cnum[0] = '\0';
    // INCHI✔️❌:             ndigit = 0;
    // INCHI✔️❌:             cano_num++;
    // INCHI✔️❌:             if (c == ';')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 icompnt++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (c == '\0')
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else if (isdigit( c ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cnum[ndigit] = c;
    // INCHI✔️❌:             ndigit++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             err = 2;
    // INCHI✔️❌:             goto exit_function;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌: exit_function:
    // INCHI✔️❌:     inchi_strbuf_reset( strbuf );
    // INCHI✔️❌:
    // INCHI✔️❌:     return err;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: InternallyGetCanoNumsAndComponentNums

    if canonical_numbers.is_null() || component_numbers.is_null() || string_buffer.pStr.is_null() {
        return Ok(1);
    }

    inchi_strbuf_reset(heap, Some(string_buffer))?;
    let sorted_inchi = output_control.pINChISort;
    let sorted_inchi2 = output_control.pINChISort2;
    let output_type = output_control.bOutType;
    let taut_mode = output_control.TAUT_MODE;
    let number_of_components = output_control.num_components;
    let second_non_taut_pass = output_control.bSecondNonTautPass;
    let omit_repetitions = output_control.bOmitRepetitions;
    output_control.tot_len = crate::source::base::ichiprt3::str_AuxNumb(
        heap,
        canonical_globals,
        sorted_inchi,
        sorted_inchi2,
        string_buffer,
        &mut output_control.bOverflow,
        output_type,
        taut_mode,
        number_of_components,
        second_non_taut_pass,
        omit_repetitions,
    )?;

    for index in 0..number_of_atoms {
        let component = component_numbers.offset(i64::from(index))?;
        *heap
            .slice_mut(component)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = -1;
        let canonical = canonical_numbers.offset(i64::from(index.wrapping_add(1)))?;
        *heap
            .slice_mut(canonical)?
            .first_mut()
            .ok_or(SourceHeapError::PointerOutOfBounds)? = -1;
    }

    let mut digits = [0_i8; 8];
    let mut digit_count = 0_usize;
    let mut error = 0_i32;
    let mut component_number = 1_i32;
    let mut canonical_number = 0_i32;
    let mut index = 0_i32;
    while index <= string_buffer.nUsedLength {
        let character_pointer = string_buffer.pStr.as_const().offset(i64::from(index))?;
        let character = *heap
            .slice(character_pointer)?
            .first()
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if character == b',' as i8 || character == b';' as i8 || character == 0 {
            *digits
                .get_mut(digit_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
            let original_number = digits[..digit_count].iter().fold(0_i32, |value, digit| {
                value
                    .wrapping_mul(10)
                    .wrapping_add(i32::from(*digit - b'0' as i8))
            });
            let canonical = canonical_numbers.offset(i64::from(original_number))?;
            *heap
                .slice_mut(canonical)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = canonical_number;
            let component = component_numbers.offset(i64::from(canonical_number))?;
            *heap
                .slice_mut(component)?
                .first_mut()
                .ok_or(SourceHeapError::PointerOutOfBounds)? = component_number;
            digits[0] = 0;
            digit_count = 0;
            canonical_number = canonical_number.wrapping_add(1);
            if character == b';' as i8 {
                component_number = component_number.wrapping_add(1);
            }
            if character == 0 {
                break;
            }
        } else if (b'0' as i8..=b'9' as i8).contains(&character) {
            *digits
                .get_mut(digit_count)
                .ok_or(SourceHeapError::PointerOutOfBounds)? = character;
            digit_count += 1;
        } else {
            error = 2;
            break;
        }
        index = index.wrapping_add(1);
    }

    inchi_strbuf_reset(heap, Some(string_buffer))?;
    Ok(error)
}

#[allow(non_snake_case)]
pub(crate) fn CleanOrigCoord(
    heap: &mut SourceHeap,
    coordinates: &mut MOL_COORD,
    delimiter: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2228 CleanOrigCoord
    // INCHI✔️❌: int CleanOrigCoord( MOL_COORD szCoord, int delim )
    // INCHI✔️❌: {
    // INCHI✔️❌: #define MIN_BOND_LENGTH   (1.0e-6)
    // INCHI✔️❌:     char szVal[LEN_COORD + 1];
    // INCHI✔️❌:     MOL_COORD szBuf;
    // INCHI✔️❌:     char *q;
    // INCHI✔️❌:     int len, last, fst, dec_pnt, num_zer = 0, len_buf = 0, e;
    // INCHI✔️❌:     int k, i;
    // INCHI✔️❌:     double coord;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (k = 0; k < NUM_COORD*LEN_COORD; k += LEN_COORD)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memcpy(szVal, szCoord + k, LEN_COORD);
    // INCHI✔️❌:         szVal[LEN_COORD] = '\0';
    // INCHI✔️❌:         lrtrim( szVal, &len );
    // INCHI✔️❌:         coord = strtod( szVal, &q );
    // INCHI✔️❌:         if (MIN_BOND_LENGTH > fabs( coord ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             strcpy(szVal, "0");
    // INCHI✔️❌:             len = 1;
    // INCHI✔️❌:             num_zer++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len = (int) ( q - szVal );
    // INCHI✔️❌:             /* last = (last mantissa digit position + 1)  */
    // INCHI✔️❌:             if (( q = strchr( szVal, 'e' ) ) || ( q = strchr( szVal, 'E' ) ) ||
    // INCHI✔️❌:                 ( q = strchr( szVal, 'd' ) ) || ( q = strchr( szVal, 'D' ) ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* floating point */
    // INCHI✔️❌:                 last = q - szVal;
    // INCHI✔️❌:                 /* remove (+) and leading zeroes from the exponent */
    // INCHI✔️❌:                 e = (int) strtol( szVal + last + 1, &q, 10 ); /* exponent */
    // INCHI✔️❌:                 if (e)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* new exp; update the length */
    // INCHI✔️❌:                     len = last + 1 + sprintf(szVal + last + 1, "%d", e); /* print exp without leading zeroes and '+' */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* exponent is zero */
    // INCHI✔️❌:                     len = last;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 last = len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* fst = (first mantissa digit); fst=1 if the sign is present, otherwise 0 */
    // INCHI✔️❌:             fst = ( szVal[0] != '.' && !isdigit( UCINT szVal[0] ) );
    // INCHI✔️❌:             /* dec_pnt = (decimal point position) or last */
    // INCHI✔️❌:             if ((q = strchr( szVal, '.' ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 dec_pnt = (int) ( q - szVal );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 dec_pnt = last;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             last -= 1; /* last mantissa digit position */
    // INCHI✔️❌:             /* remove trailing zeroes in the range dec_pnt+1..last-1 */
    // INCHI✔️❌:             for (i = last; dec_pnt < i && '0' == szVal[i]; i--)
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             if (i == dec_pnt)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 i--; /* remove decimal point, too */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (i < last)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memmove(szVal + i + 1, szVal + last + 1, (long long)len - (long long)last); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                 len -= last - i;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* remove leading zeroes */
    // INCHI✔️❌:             for (i = fst; i < len && '0' == szVal[i]; i++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if ((i > fst) && (len - fst <= LEN_COORD + 1 - i) && (len - fst <= LEN_COORD + 1 - fst)) /* djb-rwth: fixing GHI #138 */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memmove(szVal + fst, szVal + i, (long long)len - (long long)fst); /* djb-rwth: cast operator added */
    // INCHI✔️❌:                 len -= i - fst;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (len_buf && (len_buf < (int)sizeof(MOL_COORD)))
    // INCHI✔️❌:         {
    // INCHI✔️❌: #pragma warning (push)
    // INCHI✔️❌: #pragma warning (disable: 6386)
    // INCHI✔️❌:             szBuf[len_buf++] = delim;
    // INCHI✔️❌: #pragma warning (pop)
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (len_buf >= (int)sizeof(MOL_COORD)) /* djb-rwth: fixing coverity ID #499520 */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len_buf = (int)sizeof(MOL_COORD) - 1;
    // INCHI✔️❌:             len = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         memcpy(szBuf + len_buf, szVal, len); /* does not copy zero termination*/
    // INCHI✔️❌:         len_buf += len;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     /* zero termination */
    // INCHI✔️❌:     if (len_buf < ( int )sizeof( MOL_COORD ))
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memset( szBuf + len_buf, 0, sizeof( MOL_COORD ) - len_buf ); /* djb-rwth: memset_s C11/Annex K variant? */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     memcpy(szCoord, szBuf, sizeof(MOL_COORD));
    // INCHI✔️❌:
    // INCHI✔️❌:     return num_zer;
    // INCHI✔️❌: #undef MIN_BOND_LENGTH
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: CleanOrigCoord
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: CleanOrigCoord
    // INCHI✔️❌: #define LEN_COORD 10
    // INCHI✔️❌: #define NUM_COORD 3
    // INCHI✔️❌: typedef char MOL_COORD[LEN_COORD*NUM_COORD + NUM_COORD - 1]
    // INCHI✔️❌: #define UCINT (unsigned char)
    // INCHI✔️❌: C locale; COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: CleanOrigCoord

    const MIN_BOND_LENGTH: f64 = 1.0e-6;
    let mut output = [0_i8; 32];
    let mut output_length = 0_usize;
    let mut zero_count = 0_i32;

    for offset in (0..(NUM_COORD * LEN_COORD) as usize).step_by(LEN_COORD as usize) {
        let mut field = coordinates[offset..offset + LEN_COORD as usize].to_vec();
        field.push(0);
        let field_pointer = heap.allocate_model_storage(field)?;
        let parsed =
            (|| -> Result<(Vec<i8>, usize, f64, usize, Option<(usize, i32)>), SourceHeapError> {
                let mut trimmed_length = 0_i32;
                let _ = lrtrim(heap, field_pointer, Some(&mut trimmed_length))?;
                let (coordinate, end) = source_strtod_with_end(heap, field_pointer.as_const())?;
                let value = heap.slice(field_pointer.as_const())?.to_vec();
                let exponent = [b'e', b'E', b'd', b'D']
                    .into_iter()
                    .find_map(|target| value.iter().position(|byte| *byte as u8 == target));
                let exponent = match exponent {
                    Some(position) => Some((
                        position,
                        inchi_strtol(
                            heap,
                            field_pointer.as_const().offset(
                                i64::try_from(position + 1)
                                    .map_err(|_| SourceHeapError::PointerOffsetOverflow)?,
                            )?,
                            None,
                            10,
                        )? as i32,
                    )),
                    None => None,
                };
                Ok((
                    value,
                    usize::try_from(trimmed_length)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                    coordinate,
                    end,
                    exponent,
                ))
            })();
        let cleanup = heap.free(field_pointer);
        let (mut value, mut length, coordinate, parsed_end, exponent) = match (parsed, cleanup) {
            (Ok(value), Ok(())) => value,
            (Err(error), Ok(())) | (Ok(_), Err(error)) => return Err(error),
            (Err(error), Err(_)) => return Err(error),
        };

        if coordinate.abs() < MIN_BOND_LENGTH {
            value[0] = b'0' as i8;
            value[1] = 0;
            length = 1;
            zero_count = zero_count.wrapping_add(1);
        } else {
            length = parsed_end;
            let mut mantissa_end = if let Some((exponent, exponent_value)) = exponent {
                if exponent_value != 0 {
                    let digits = exponent_value.to_string();
                    for (index, byte) in digits.bytes().enumerate() {
                        value[exponent + 1 + index] = byte as i8;
                    }
                    length = exponent + 1 + digits.len();
                    value[length] = 0;
                } else {
                    length = exponent;
                }
                exponent
            } else {
                length
            };

            let first = usize::from(value[0] as u8 != b'.' && !(value[0] as u8).is_ascii_digit());
            let value_end = value
                .iter()
                .position(|byte| *byte == 0)
                .unwrap_or(value.len());
            let decimal_point = value[..value_end]
                .iter()
                .position(|byte| *byte as u8 == b'.')
                .unwrap_or(mantissa_end);
            if mantissa_end != 0 {
                mantissa_end -= 1;
                let mut last_nonzero = mantissa_end;
                while decimal_point < last_nonzero && value[last_nonzero] as u8 == b'0' {
                    last_nonzero -= 1;
                }
                if last_nonzero == decimal_point {
                    last_nonzero = last_nonzero.saturating_sub(1);
                }
                if last_nonzero < mantissa_end {
                    value.copy_within(mantissa_end + 1..=length, last_nonzero + 1);
                    length -= mantissa_end - last_nonzero;
                }
            }

            let mut leading_end = first;
            while leading_end < length && value[leading_end] as u8 == b'0' {
                leading_end += 1;
            }
            if leading_end > first
                && length - first <= LEN_COORD as usize + 1 - leading_end
                && length - first <= LEN_COORD as usize + 1 - first
            {
                let count = length - first;
                value.copy_within(leading_end..leading_end + count, first);
                length -= leading_end - first;
            }
        }

        if output_length != 0 && output_length < output.len() {
            output[output_length] = delimiter as i8;
            output_length += 1;
        }
        if output_length >= output.len() {
            output_length = output.len() - 1;
            length = 0;
        }
        output[output_length..output_length + length].copy_from_slice(&value[..length]);
        output_length += length;
    }
    if output_length < output.len() {
        output[output_length..].fill(0);
    }
    *coordinates = output;
    Ok(zero_count)
}

#[allow(non_snake_case)]
pub(crate) fn WriteOrigCoord(
    heap: &mut SourceHeap,
    num_inp_atoms: i32,
    molecular_coordinates: SourceConstPointer<MOL_COORD>,
    index: &mut i32,
    buffer: SourceMutPointer<i8>,
    buffer_length: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2340 WriteOrigCoord
    // INCHI✔️❌: int WriteOrigCoord( int       num_inp_atoms,
    // INCHI✔️❌:                     MOL_COORD *szMolCoord,
    // INCHI✔️❌:                     int       *i,
    // INCHI✔️❌:                     char      *szBuf,
    // INCHI✔️❌:                     int       buf_len )
    // INCHI✔️❌: {
    // INCHI✔️❌:
    // INCHI✔️❌:     int j, num_zer, len, cur_len;
    // INCHI✔️❌:     char *p;
    // INCHI✔️❌:     MOL_COORD szCurCoord;
    // INCHI✔️❌:     cur_len = 0;
    // INCHI✔️❌:     for (j = *i; j < num_inp_atoms; )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memcpy(szCurCoord, szMolCoord[j], sizeof(szCurCoord));
    // INCHI✔️❌:         num_zer = CleanOrigCoord( szCurCoord, ',' );
    // INCHI✔️❌:         if (NUM_COORD == num_zer)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((p = (char *) memchr( szCurCoord, '\0', sizeof( szCurCoord ) ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len = (int) ( p - szCurCoord );
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len = sizeof( szCurCoord );
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (len + cur_len + 1 < buf_len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (len)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 memcpy(szBuf + cur_len, szCurCoord, len * sizeof(szBuf[0]));
    // INCHI✔️❌:             }
    // INCHI✔️❌:             szBuf[cur_len += len] = ';';
    // INCHI✔️❌:             cur_len++;
    // INCHI✔️❌:             j++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     szBuf[cur_len] = '\0';
    // INCHI✔️❌:     *i = j; /* next item */
    // INCHI✔️❌:
    // INCHI✔️❌:     return cur_len;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: WriteOrigCoord
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: WriteOrigCoord
    // INCHI✔️❌: #define NUM_COORD 3
    // INCHI✔️❌: typedef char MOL_COORD[LEN_COORD*NUM_COORD + NUM_COORD - 1]
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: WriteOrigCoord

    let mut current_index = *index;
    if current_index < 0 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let coordinates = if current_index < num_inp_atoms {
        let count =
            usize::try_from(num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        heap.slice(molecular_coordinates)?
            .get(..count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    } else {
        Vec::new()
    };
    let mut current_length = 0_i32;
    while current_index < num_inp_atoms {
        let mut coordinate = *coordinates
            .get(usize::try_from(current_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let zero_count = CleanOrigCoord(heap, &mut coordinate, b',' as i32)?;
        let length = if zero_count == NUM_COORD as i32 {
            0_i32
        } else {
            i32::try_from(
                coordinate
                    .iter()
                    .position(|byte| *byte == 0)
                    .unwrap_or(coordinate.len()),
            )
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
        };
        if length.wrapping_add(current_length).wrapping_add(1) < buffer_length {
            let length_usize =
                usize::try_from(length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let current_usize =
                usize::try_from(current_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            if length_usize != 0 {
                heap.slice_mut(buffer)?
                    .get_mut(current_usize..current_usize + length_usize)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .copy_from_slice(&coordinate[..length_usize]);
            }
            current_length = current_length.wrapping_add(length);
            *heap
                .slice_mut(buffer)?
                .get_mut(
                    usize::try_from(current_length)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)? = b';' as i8;
            current_length = current_length.wrapping_add(1);
            current_index = current_index.wrapping_add(1);
        } else {
            break;
        }
    }
    *heap
        .slice_mut(buffer)?
        .get_mut(usize::try_from(current_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    *index = current_index;
    Ok(current_length)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn WriteOrigAtoms(
    heap: &mut SourceHeap,
    _canon_globals: &mut CANON_GLOBALS,
    num_inp_atoms: i32,
    atoms: SourceConstPointer<inp_ATOM>,
    index: &mut i32,
    buffer: SourceMutPointer<i8>,
    buffer_length: i32,
    structure_data: &STRUCT_DATA,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2410 WriteOrigAtoms
    // INCHI✔️❌: int WriteOrigAtoms( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                     int           num_inp_atoms,
    // INCHI✔️❌:                     inp_ATOM      *at,
    // INCHI✔️❌:                     int           *i,
    // INCHI✔️❌:                     char          *szBuf,
    // INCHI✔️❌:                     int           buf_len,
    // INCHI✔️❌:                     STRUCT_DATA   *sd )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, k, n, len, len0, cur_len, val, bonds_val, mw, parity, num_trans, is_ok, b_self;
    // INCHI✔️❌:     static const char szIsoH[] = "hdt";
    // INCHI✔️❌:     char szCurAtom[32];
    // INCHI✔️❌:     AT_NUMB nNeighOrder[MAXVAL], neigh;
    // INCHI✔️❌:
    // INCHI✔️❌:     cur_len = 0;
    // INCHI✔️❌:     if (0 == *i)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         cur_len = sprintf(szBuf, "%d%s", num_inp_atoms,
    // INCHI✔️❌:             (sd->bChiralFlag & FLAG_INP_AT_CHIRAL) ? "c" :
    // INCHI✔️❌:             (sd->bChiralFlag & FLAG_INP_AT_NONCHIRAL) ? "n" : "");
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (j = *i; j < num_inp_atoms; )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* tetrahedral parity treatment */
    // INCHI✔️❌:         parity = 0;
    // INCHI✔️❌:         num_trans = 0;
    // INCHI✔️❌:         if (at[j].p_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* verify neighbors */
    // INCHI✔️❌:             is_ok = 1;
    // INCHI✔️❌:             b_self = 0;
    // INCHI✔️❌:             for (n = 0, k = 0; n < MAX_NUM_STEREO_ATOM_NEIGH; n++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 neigh = at[j].p_orig_at_num[n] - 1;
    // INCHI✔️❌:                 if (is_in_the_list( at[j].neighbor, neigh, at[j].valence ) &&
    // INCHI✔️❌:                      at[neigh].orig_at_number == at[j].p_orig_at_num[n])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* real neighbor */
    // INCHI✔️❌:                     nNeighOrder[k++] = at[j].p_orig_at_num[n];
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if ((int) neigh == j && at[neigh].orig_at_number == at[j].p_orig_at_num[n])
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         /* central atom is a neighbor */
    // INCHI✔️❌:                         num_trans = n; /* move this neighbor to 0 position permutation parity */
    // INCHI✔️❌:                         b_self++;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         is_ok = 0;
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             if (is_ok && b_self <= 1 && b_self + k == MAX_NUM_STEREO_ATOM_NEIGH)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 num_trans += insertions_sort( pCG, nNeighOrder, k, sizeof( nNeighOrder[0] ), comp_AT_RANK );
    // INCHI✔️❌:                 if (ATOM_PARITY_WELL_DEF( at[j].p_parity ))
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     parity = 2 - ( num_trans + at[j].p_parity ) % 2;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (ATOM_PARITY_ILL_DEF( at[j].p_parity ))
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         parity = at[j].p_parity;
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         ; /* invalid atom parity */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             else
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 ;/* add error message here */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         len = len0 = (int) strlen( at[j].elname );
    // INCHI✔️❌:
    // INCHI✔️❌:         memcpy(szCurAtom, at[j].elname, len);
    // INCHI✔️❌:         bonds_val = nBondsValenceInpAt( at + j, NULL, NULL );
    // INCHI✔️❌:
    // INCHI✔️❌:         if (( val = needed_unusual_el_valence( at[j].el_number, at[j].charge, at[j].radical,
    // INCHI✔️❌:             at[j].chem_bonds_valence, bonds_val, at[j].num_H, at[j].valence ) ) ||
    // INCHI✔️❌:              at[j].charge || at[j].radical || at[j].iso_atw_diff || NUM_ISO_H( at, j ) || parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* valence */
    // INCHI✔️❌:             if (val)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len += sprintf(szCurAtom + len, "%d", val > 0 ? val : 0);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* charge */
    // INCHI✔️❌:             if ((val = at[j].charge)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 szCurAtom[len++] = val > 0 ? '+' : '-';
    // INCHI✔️❌:                 if (( val = abs( val ) ) > 1)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     len += sprintf(szCurAtom + len, "%d", val);
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* radical */
    // INCHI✔️❌:             if ((val = at[j].radical)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len += sprintf(szCurAtom + len, ".%d", val);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* isotopic shift */
    // INCHI✔️❌:             if ((val = at[j].iso_atw_diff)) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 mw = get_atomic_mass_from_elnum( at[j].el_number );
    // INCHI✔️❌:                 if (val == 1)
    // INCHI✔️❌:                     val = mw;
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                     if (val > 0)
    // INCHI✔️❌:                         val = mw + val - 1;
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                         val = mw + val;
    // INCHI✔️❌:
    // INCHI✔️❌:                 len += sprintf(szCurAtom + len, "%si%d", len == len0 ? "." : "", val);
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* parity */
    // INCHI✔️❌:             if (parity)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 len += sprintf(szCurAtom + len, "%s%s", len == len0 ? "." : "",
    // INCHI✔️❌:                     parity == AB_PARITY_ODD ? "o" :
    // INCHI✔️❌:                     parity == AB_PARITY_EVEN ? "e" :
    // INCHI✔️❌:                     parity == AB_PARITY_UNKN ? "u" :
    // INCHI✔️❌:                     parity == AB_PARITY_UNDF ? "?" : "");
    // INCHI✔️❌:             }
    // INCHI✔️❌:             /* implicit isotopic H */
    // INCHI✔️❌:             if (NUM_ISO_H( at, j ))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 for (k = 0; k < NUM_H_ISOTOPES; k++)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if ((val = at[j].num_iso_H[k])) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         len += sprintf(szCurAtom + len, "%s%c", len == len0 ? "." : "", szIsoH[k]);
    // INCHI✔️❌:                         if (val > 1)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             len += sprintf(szCurAtom + len, "%d", val);
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (len + cur_len < buf_len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memcpy(szBuf + cur_len, szCurAtom, len);
    // INCHI✔️❌:             cur_len += len;
    // INCHI✔️❌:             j++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         szBuf[cur_len] = '\0';
    // INCHI✔️❌:         *i = j;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return cur_len;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: WriteOrigAtoms
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: WriteOrigAtoms
    // INCHI✔️❌: #define MAXVAL 20
    // INCHI✔️❌: #define MAX_NUM_STEREO_ATOM_NEIGH 4
    // INCHI✔️❌: #define NUM_H_ISOTOPES 3
    // INCHI✔️❌: #define NUM_ISO_H(AT, N) (AT[N].num_iso_H[0] + AT[N].num_iso_H[1] + AT[N].num_iso_H[2])
    // INCHI✔️❌: #define ATOM_PARITY_WELL_DEF(X) (1 <= (X) && (X) <= 2)
    // INCHI✔️❌: #define ATOM_PARITY_ILL_DEF(X) (3 <= (X) && (X) <= 4)
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: WriteOrigAtoms

    if *index < 0 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let atom_count =
        usize::try_from(num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = if atom_count == 0 {
        Vec::new()
    } else {
        heap.slice(atoms)?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let mut current_length = 0_i32;
    if *index == 0 {
        let suffix = if structure_data.bChiralFlag & FLAG_INP_AT_CHIRAL as i32 != 0 {
            "c"
        } else if structure_data.bChiralFlag & FLAG_INP_AT_NONCHIRAL as i32 != 0 {
            "n"
        } else {
            ""
        };
        let header = format!("{num_inp_atoms}{suffix}");
        let destination = heap
            .slice_mut(buffer)?
            .get_mut(..=header.len())
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        for (target, source) in destination.iter_mut().zip(header.bytes().chain([0])) {
            *target = source as i8;
        }
        current_length =
            i32::try_from(header.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    }

    let mut current_index = *index;
    while current_index < num_inp_atoms {
        let atom = atoms
            .get(usize::try_from(current_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut parity = 0_i32;
        let mut num_trans = 0_i32;
        if atom.p_parity != 0 {
            let mut is_ok = true;
            let mut self_count = 0_i32;
            let mut neighbor_order = [0_u16; MAXVAL as usize];
            let mut neighbor_count = 0_usize;
            for n in 0..MAX_NUM_STEREO_ATOM_NEIGH as usize {
                let neighbor = atom.p_orig_at_num[n].wrapping_sub(1);
                let neighbor_atom = atoms
                    .get(usize::from(neighbor))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if is_in_the_list(Some(&atom.neighbor), neighbor, i32::from(atom.valence))?
                    .is_some()
                    && neighbor_atom.orig_at_number == atom.p_orig_at_num[n]
                {
                    neighbor_order[neighbor_count] = atom.p_orig_at_num[n];
                    neighbor_count += 1;
                } else if i32::from(neighbor) == current_index
                    && neighbor_atom.orig_at_number == atom.p_orig_at_num[n]
                {
                    num_trans = n as i32;
                    self_count = self_count.wrapping_add(1);
                } else {
                    is_ok = false;
                    break;
                }
            }
            if is_ok
                && self_count <= 1
                && self_count + neighbor_count as i32 == MAX_NUM_STEREO_ATOM_NEIGH as i32
            {
                num_trans = num_trans.wrapping_add(insertions_sort_AT_RANK(
                    &mut neighbor_order,
                    neighbor_count as i32,
                )?);
                let source_parity = i32::from(atom.p_parity);
                if (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
                    .contains(&source_parity)
                {
                    parity = 2 - num_trans.wrapping_add(source_parity) % 2;
                } else if (AB_MIN_ILL_DEFINED_PARITY as i32..=AB_MAX_ILL_DEFINED_PARITY as i32)
                    .contains(&source_parity)
                {
                    parity = source_parity;
                }
            }
        }

        let element_length = atom
            .elname
            .iter()
            .position(|byte| *byte == 0)
            .ok_or(SourceHeapError::MissingNulTerminator)?;
        let mut current_atom = atom.elname[..element_length]
            .iter()
            .map(|byte| *byte as u8)
            .collect::<Vec<_>>();
        let bonds_valence = nBondsValenceInpAt(atom, None, None);
        let unusual_valence = needed_unusual_el_valence(
            i32::from(atom.el_number),
            i32::from(atom.charge),
            i32::from(atom.radical),
            i32::from(atom.chem_bonds_valence),
            bonds_valence,
            i32::from(atom.num_H),
            i32::from(atom.valence),
        )?;
        let isotope_hydrogens = atom
            .num_iso_H
            .iter()
            .map(|value| i32::from(*value))
            .sum::<i32>();
        if unusual_valence != 0
            || atom.charge != 0
            || atom.radical != 0
            || atom.iso_atw_diff != 0
            || isotope_hydrogens != 0
            || parity != 0
        {
            if unusual_valence != 0 {
                current_atom.extend(
                    (if unusual_valence > 0 {
                        unusual_valence
                    } else {
                        0
                    })
                    .to_string()
                    .bytes(),
                );
            }
            let charge = i32::from(atom.charge);
            if charge != 0 {
                current_atom.push(if charge > 0 { b'+' } else { b'-' });
                if charge.abs() > 1 {
                    current_atom.extend(charge.abs().to_string().bytes());
                }
            }
            let radical = i32::from(atom.radical);
            if radical != 0 {
                current_atom.push(b'.');
                current_atom.extend(radical.to_string().bytes());
            }
            let isotope_difference = i32::from(atom.iso_atw_diff);
            if isotope_difference != 0 {
                let mass = get_atomic_mass_from_elnum(i32::from(atom.el_number))?;
                let isotope_mass = if isotope_difference == 1 {
                    mass
                } else if isotope_difference > 0 {
                    mass.wrapping_add(isotope_difference).wrapping_sub(1)
                } else {
                    mass.wrapping_add(isotope_difference)
                };
                if current_atom.len() == element_length {
                    current_atom.push(b'.');
                }
                current_atom.push(b'i');
                current_atom.extend(isotope_mass.to_string().bytes());
            }
            if parity != 0 {
                if current_atom.len() == element_length {
                    current_atom.push(b'.');
                }
                if let Some(symbol) = match parity {
                    value if value == AB_PARITY_ODD as i32 => Some(b'o'),
                    value if value == AB_PARITY_EVEN as i32 => Some(b'e'),
                    value if value == AB_PARITY_UNKN as i32 => Some(b'u'),
                    value if value == AB_PARITY_UNDF as i32 => Some(b'?'),
                    _ => None,
                } {
                    current_atom.push(symbol);
                }
            }
            for (isotope, symbol) in atom.num_iso_H.iter().zip(*b"hdt") {
                let count = i32::from(*isotope);
                if count != 0 {
                    if current_atom.len() == element_length {
                        current_atom.push(b'.');
                    }
                    current_atom.push(symbol);
                    if count > 1 {
                        current_atom.extend(count.to_string().bytes());
                    }
                }
            }
        }
        if current_atom.len() > 32 {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let atom_length = i32::try_from(current_atom.len())
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if atom_length.wrapping_add(current_length) < buffer_length {
            let start =
                usize::try_from(current_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            heap.slice_mut(buffer)?
                .get_mut(start..start + current_atom.len())
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(
                    &current_atom
                        .iter()
                        .map(|byte| *byte as i8)
                        .collect::<Vec<_>>(),
                );
            current_length = current_length.wrapping_add(atom_length);
            current_index = current_index.wrapping_add(1);
        } else {
            break;
        }
        *heap
            .slice_mut(buffer)?
            .get_mut(
                usize::try_from(current_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?,
            )
            .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
        *index = current_index;
    }
    Ok(current_length)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
// Preserve official source-frame trailing whitespace that rustfmt otherwise removes.
#[rustfmt::skip]
pub(crate) fn WriteOrigBonds(
    heap: &mut SourceHeap,
    canon_globals: &mut CANON_GLOBALS,
    num_inp_atoms: i32,
    atoms: SourceMutPointer<inp_ATOM>,
    index: &mut i32,
    buffer: SourceMutPointer<i8>,
    buffer_length: i32,
    structure_data: &mut STRUCT_DATA,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2610 WriteOrigBonds
    // INCHI✔️❌: int WriteOrigBonds( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                     int           num_inp_atoms,
    // INCHI✔️❌:                     inp_ATOM      *at,
    // INCHI✔️❌:                     int           *i,
    // INCHI✔️❌:                     char          *szBuf,
    // INCHI✔️❌:                     int           buf_len,
    // INCHI✔️❌:                     STRUCT_DATA   *sd )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int j, k, k2, kk, len, cur_len, j2 = 0, bond_stereo, bond_char, bond_parity, bond_parityNM, num_trans; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:     char szCurBonds[7 * MAXVAL + 2]; /* num_neigh*(1 byte bond type + 2 bytes for bond parity up to 4 digits per neighbor number) + at the end one ';' */
    // INCHI✔️❌:     AT_RANK nNeighOrder[MAXVAL];
    // INCHI✔️❌:     int  chain_len, pnxt_atom, pinxt2cur, pinxt_sb_parity_ord;
    // INCHI✔️❌:     int  chain_len2, pnxt_atom2, pinxt2cur2, pinxt_sb_parity_ord2, m1, m2;
    // INCHI✔️❌:     int  pcur_atom, picur2nxt, picur_sb_parity_ord;
    // INCHI✔️❌:
    // INCHI✔️❌:     cur_len = 0;
    // INCHI✔️❌:     for (j = *i; j < num_inp_atoms; )
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len = 0;
    // INCHI✔️❌:         if (at[j].valence >= 1) /* djb-rwth: changing condition to avoid garbage values */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             for (k = 0; k < at[j].valence; k++)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 nNeighOrder[k] = k;
    // INCHI✔️❌:             }
    // INCHI✔️❌:             pCG->m_pn_RankForSort = at[j].neighbor;
    // INCHI✔️❌:             num_trans = insertions_sort( pCG, nNeighOrder, at[j].valence, sizeof( nNeighOrder[0] ), CompRank ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             num_trans = 0; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
    // INCHI✔️❌:             nNeighOrder[0] = 0;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         for (kk = 0; kk < at[j].valence; kk++) 
    // INCHI✔️❌:         {
    // INCHI✔️❌:             k = nNeighOrder[kk];
    // INCHI✔️❌:             j2 = at[j].neighbor[k];
    // INCHI✔️❌:             bond_parity = 0;
    // INCHI✔️❌:             bond_parityNM = 0;
    // INCHI✔️❌:             if (j2 < j)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 bond_stereo = at[j].bond_stereo[k];
    // INCHI✔️❌:                 switch (at[j].bond_type[k])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     case BOND_TYPE_SINGLE:
    // INCHI✔️❌:                         switch (bond_stereo)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case  STEREO_SNGL_UP:
    // INCHI✔️❌:                                 bond_char = 'p';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case -STEREO_SNGL_UP:
    // INCHI✔️❌:                                 bond_char = 'P';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case  STEREO_SNGL_DOWN:
    // INCHI✔️❌:                                 bond_char = 'n';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case -STEREO_SNGL_DOWN:
    // INCHI✔️❌:                                 bond_char = 'N';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌: #if ( FIX_EITHER_STEREO_IN_AUX_INFO == 1 )
    // INCHI✔️❌:                             case  STEREO_SNGL_EITHER:
    // INCHI✔️❌:                                 bond_char = 'v';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             case -STEREO_SNGL_EITHER:
    // INCHI✔️❌:                                 bond_char = 'V';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌: #else
    // INCHI✔️❌:                             case  STEREO_SNGL_EITHER:
    // INCHI✔️❌:                             case -STEREO_SNGL_EITHER:
    // INCHI✔️❌:                                 bond_char = 'v';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌: #endif
    // INCHI✔️❌:                             default:
    // INCHI✔️❌:                                 bond_char = 's';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case BOND_TYPE_DOUBLE:
    // INCHI✔️❌:                         switch (bond_stereo)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             case  STEREO_DBLE_EITHER:
    // INCHI✔️❌:                             case -STEREO_DBLE_EITHER:
    // INCHI✔️❌:                                 bond_char = 'w';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                             default:
    // INCHI✔️❌:                                 bond_char = 'd';
    // INCHI✔️❌:                                 break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case BOND_TYPE_TRIPLE:
    // INCHI✔️❌:                         bond_char = 't';
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     case BOND_TYPE_ALTERN:
    // INCHI✔️❌:                         bond_char = 'a';
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                     default:
    // INCHI✔️❌:                         bond_char = 's';
    // INCHI✔️❌:                         break;
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /* check for allene/cumulene */
    // INCHI✔️❌:                 k2 = (int) ( is_in_the_list( at[j2].neighbor, (AT_NUMB) j, at[j2].valence ) - at[j2].neighbor );
    // INCHI✔️❌:                 chain_len = chain_len2 = 0;
    // INCHI✔️❌:                 if (at[j].sb_parity[0])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (m1 = 0; m1 < MAX_NUM_STEREO_BONDS && at[j].sb_parity[m1]; m1++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (k == at[j].sb_ord[m1])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             chain_len = get_opposite_sb_atom( at, j, k,
    // INCHI✔️❌:                                           &pnxt_atom, &pinxt2cur, &pinxt_sb_parity_ord );
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if (at[j2].sb_parity[0])
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     for (m2 = 0; m2 < MAX_NUM_STEREO_BONDS && at[j2].sb_parity[m2]; m2++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (k2 == at[j2].sb_ord[m2])
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             chain_len2 = get_opposite_sb_atom( at, j2, k2,
    // INCHI✔️❌:                                            &pnxt_atom2, &pinxt2cur2, &pinxt_sb_parity_ord2 );
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 if ((chain_len == 1 && chain_len2 == 1) ||  /* regular stereobond */
    // INCHI✔️❌:                      (chain_len > 1 && j > pnxt_atom)) /* djb-rwth: addressing LLVM warnings */
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* j  is a cumulene endpoint */
    // INCHI✔️❌:                     int m;
    // INCHI✔️❌:                     pcur_atom = j;  /* pcur_atom > pnxt_atom */
    // INCHI✔️❌:                     picur2nxt = k;
    // INCHI✔️❌:                     picur_sb_parity_ord = -1;
    // INCHI✔️❌:                     for (m = 0; m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m]; m++)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (at[pcur_atom].sb_ord[m] == k)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             picur_sb_parity_ord = m;
    // INCHI✔️❌:                             break;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 else
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     if (chain_len2 > 1 && j2 > pnxt_atom2)
    // INCHI✔️❌:                     { /* j2 is a cumulene endpoint */
    // INCHI✔️❌:                         int m;
    // INCHI✔️❌:                         pcur_atom = j2;
    // INCHI✔️❌:                         picur2nxt = k2;
    // INCHI✔️❌:                         pnxt_atom = pnxt_atom2;
    // INCHI✔️❌:                         pinxt2cur = pinxt2cur2;
    // INCHI✔️❌:                         pinxt_sb_parity_ord = pinxt_sb_parity_ord2;
    // INCHI✔️❌:                         picur_sb_parity_ord = -1;
    // INCHI✔️❌:                         for (m = 0; m < MAX_NUM_STEREO_BONDS && at[pcur_atom].sb_parity[m]; m++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (at[pcur_atom].sb_ord[m] == k2)
    // INCHI✔️❌:                                 picur_sb_parity_ord = m;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         chain_len = chain_len2;
    // INCHI✔️❌:                         /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         chain_len = 0; /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:                 /*len += sprintf( szCurBonds + len, "%c%d", bond_char, val+1);*/
    // INCHI✔️❌:                 if (chain_len)
    // INCHI✔️❌:                 {
    // INCHI✔️❌:                     /* both atoms belong to a stereo bond */
    // INCHI✔️❌:                     int kc;
    // INCHI✔️❌:                     int p1 = 0, p2, p1NM = 0, p2NM, neigh, neigh1, neigh2, bHasMetal, bWellDef; /* djb-rwth: initialising p1 and p1NM */
    // INCHI✔️❌:                     int     bNeighSwitched1, bNeighSwitched2;
    // INCHI✔️❌:
    // INCHI✔️❌:                     /* djb-rwth: avoiding buffer overrun as picur_sb_parity_ord == -1 is possible */
    // INCHI✔️❌:                     if (picur_sb_parity_ord >= 0)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         p1 = SB_PARITY_1( at[pcur_atom].sb_parity[picur_sb_parity_ord] );
    // INCHI✔️❌:                         p1NM = SB_PARITY_2( at[pcur_atom].sb_parity[picur_sb_parity_ord] );
    // INCHI✔️❌:                     }
    // INCHI✔️❌:
    // INCHI✔️❌:                     p2 = SB_PARITY_1( at[pnxt_atom].sb_parity[pinxt_sb_parity_ord] );
    // INCHI✔️❌:                     p2NM = SB_PARITY_2( at[pnxt_atom].sb_parity[pinxt_sb_parity_ord] );
    // INCHI✔️❌:
    // INCHI✔️❌:                     bWellDef = ATOM_PARITY_WELL_DEF( p1 ) && ATOM_PARITY_WELL_DEF( p2 );
    // INCHI✔️❌:                     bHasMetal = ATOM_PARITY_WELL_DEF( p1NM ) && ATOM_PARITY_WELL_DEF( p2NM );
    // INCHI✔️❌:
    // INCHI✔️❌:                     bNeighSwitched1 = bNeighSwitched2 = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:                     if (bWellDef || bHasMetal)
    // INCHI✔️❌:                     {
    // INCHI✔️❌:
    // INCHI✔️❌:                         neigh1 = num_inp_atoms;
    // INCHI✔️❌:                         for (kc = 0; kc < at[pcur_atom].valence; kc++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (kc == picur2nxt)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             neigh = at[pcur_atom].neighbor[kc];
    // INCHI✔️❌:                             if (bHasMetal && is_el_a_metal( at[neigh].el_number ))
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             if (neigh < neigh1)
    // INCHI✔️❌:                                 neigh1 = neigh;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (neigh1 < num_inp_atoms)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bNeighSwitched1 = ( neigh1 != at[pcur_atom].neighbor[(int) at[pcur_atom].sn_ord[picur_sb_parity_ord]] );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             AddErrorMessage( sd->pStrErrStruct, "Cannot find 0D stereobond neighbor" );
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             sd->nStructReadError =  99;
    // INCHI✔️❌:                             sd->nErrorType = _IS_ERROR;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         neigh2 = num_inp_atoms;
    // INCHI✔️❌:                         for (kc = 0; kc < at[pnxt_atom].valence; kc++)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (kc == pinxt2cur)
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             neigh = at[pnxt_atom].neighbor[kc];
    // INCHI✔️❌:                             if (bHasMetal && is_el_a_metal( at[neigh].el_number ))
    // INCHI✔️❌:                                 continue;
    // INCHI✔️❌:                             if (neigh < neigh2)
    // INCHI✔️❌:                                 neigh2 = neigh;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (neigh2 < num_inp_atoms)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bNeighSwitched2 = ( neigh2 != at[pnxt_atom].neighbor[(int) at[pnxt_atom].sn_ord[pinxt_sb_parity_ord]] );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         else
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             AddErrorMessage( sd->pStrErrStruct, "Cannot find 0D stereobond neighbor" );
    // INCHI✔️❌:                             /*
    // INCHI✔️❌:                             sd->nStructReadError =  99;
    // INCHI✔️❌:                             sd->nErrorType = _IS_ERROR;
    // INCHI✔️❌:                             */
    // INCHI✔️❌:                         }
    // INCHI✔️❌:
    // INCHI✔️❌:                         if (neigh1 < num_inp_atoms && neigh2 < num_inp_atoms)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             if (ATOM_PARITY_WELL_DEF( p1 ) && ATOM_PARITY_WELL_DEF( p2 ))
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 bond_parity = 2 - ( p1 + p2 + bNeighSwitched1 + bNeighSwitched2 ) % 2;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 bond_parity = inchi_min( p1, p2 );
    // INCHI✔️❌:                             }
    // INCHI✔️❌:
    // INCHI✔️❌:                             if (bHasMetal)
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 bond_parityNM = 2 - ( p1NM + p2NM + bNeighSwitched1 + bNeighSwitched2 ) % 2;
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                             else
    // INCHI✔️❌:                             {
    // INCHI✔️❌:                                 if (p1NM && p2NM)
    // INCHI✔️❌:                                 {
    // INCHI✔️❌:                                     bond_parityNM = inchi_min( p1NM, p2NM );
    // INCHI✔️❌:                                 }
    // INCHI✔️❌:                             }
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                     else
    // INCHI✔️❌:                     {
    // INCHI✔️❌:                         if (p1 && p2)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bond_parity = inchi_min( p1, p2 );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (p1NM && p2NM)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bond_parityNM = inchi_min( p1NM, p2NM );
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                         if (bond_parityNM && !bond_parity)
    // INCHI✔️❌:                         {
    // INCHI✔️❌:                             bond_parity = AB_PARITY_UNDF;
    // INCHI✔️❌:                         }
    // INCHI✔️❌:                     }
    // INCHI✔️❌:                 }
    // INCHI✔️❌:
    // INCHI✔️❌:                 len += sprintf(szCurBonds + len, "%c%s%s%d",
    // INCHI✔️❌:                     bond_char,
    // INCHI✔️❌:
    // INCHI✔️❌:                     (bond_parity == AB_PARITY_ODD) ? "-" :
    // INCHI✔️❌:                     (bond_parity == AB_PARITY_EVEN) ? "+" :
    // INCHI✔️❌:                     (bond_parity == AB_PARITY_UNKN) ? "u" :
    // INCHI✔️❌:                     (bond_parity == AB_PARITY_UNDF) ? "?" : "",
    // INCHI✔️❌:
    // INCHI✔️❌:                     (bond_parityNM == AB_PARITY_ODD) ? "-" :
    // INCHI✔️❌:                     (bond_parityNM == AB_PARITY_EVEN) ? "+" :
    // INCHI✔️❌:                     (bond_parityNM == AB_PARITY_UNKN) ? "u" :
    // INCHI✔️❌:                     (bond_parityNM == AB_PARITY_UNDF) ? "?" : "",
    // INCHI✔️❌:
    // INCHI✔️❌:                     j2 + 1);
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (len + cur_len + 2 < buf_len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             memcpy(szBuf + cur_len, szCurBonds, len);
    // INCHI✔️❌:             cur_len += len;
    // INCHI✔️❌:             szBuf[cur_len++] = ';';
    // INCHI✔️❌:             j++;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     szBuf[cur_len] = '\0';
    // INCHI✔️❌:     *i = num_inp_atoms > 0 ? j : 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     return cur_len;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: WriteOrigBonds
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: WriteOrigBonds
    // INCHI✔️❌: #define MAXVAL 20
    // INCHI✔️❌: #define MAX_NUM_STEREO_BONDS 3
    // INCHI✔️❌: #define FIX_EITHER_STEREO_IN_AUX_INFO 1
    // INCHI✔️❌: #define SB_PARITY_1(X) (X & SB_PARITY_MASK)
    // INCHI✔️❌: #define SB_PARITY_2(X) (((X) >> SB_PARITY_SHFT) & SB_PARITY_MASK)
    // INCHI✔️❌: #define SB_PARITY_SHFT 3
    // INCHI✔️❌: #define SB_PARITY_MASK 0x07
    // INCHI✔️❌: #define ATOM_PARITY_WELL_DEF(X) (1 <= (X) && (X) <= 2)
    // INCHI✔️❌: #define inchi_min(a,b) ((a)<(b)?(a):(b))
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: WriteOrigBonds

    if *index < 0 {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let atom_count =
        usize::try_from(num_inp_atoms).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let atoms = if atom_count == 0 {
        Vec::new()
    } else {
        heap.slice(atoms.as_const())?
            .get(..atom_count)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .to_vec()
    };
    let parity_part = |value: i8| i32::from(value) & SB_PARITY_MASK as i32;
    let parity_part_non_metal =
        |value: i8| (i32::from(value) >> SB_PARITY_SHFT) & SB_PARITY_MASK as i32;
    let is_well_defined = |value: i32| {
        (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32).contains(&value)
    };
    let parity_symbol = |value: i32| match value {
        value if value == AB_PARITY_ODD as i32 => Some(b'-'),
        value if value == AB_PARITY_EVEN as i32 => Some(b'+'),
        value if value == AB_PARITY_UNKN as i32 => Some(b'u'),
        value if value == AB_PARITY_UNDF as i32 => Some(b'?'),
        _ => None,
    };

    let mut current_length = 0_i32;
    let mut current_index = *index;
    while current_index < num_inp_atoms {
        let atom_index =
            usize::try_from(current_index).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let atom = atoms
            .get(atom_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let valence = i32::from(atom.valence);
        let valence_count =
            usize::try_from(valence.max(0)).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if valence_count > MAXVAL as usize {
            return Err(SourceHeapError::PointerOutOfBounds);
        }

        let mut neighbor_order = [0 as AT_RANK; MAXVAL as usize];
        if valence >= 1 {
            for (ordinal, target) in neighbor_order[..valence_count].iter_mut().enumerate() {
                *target = ordinal as AT_RANK;
            }
            let rank_pointer = heap.allocate_model_storage(atom.neighbor.to_vec())?;
            canon_globals.m_pn_RankForSort = rank_pointer.as_const();

            let ranks = atom.neighbor;
            let width = std::mem::size_of::<AT_RANK>();
            let mut order_bytes = Vec::with_capacity(valence_count * width);
            for ordinal in &neighbor_order[..valence_count] {
                order_bytes.extend_from_slice(&ordinal.to_ne_bytes());
            }
            let mut compare = |left: &[u8], right: &[u8]| {
                let left_rank = AT_RANK::from_ne_bytes(
                    left.try_into()
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                );
                let right_rank = AT_RANK::from_ne_bytes(
                    right
                        .try_into()
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                );
                let left_value = *ranks
                    .get(usize::from(left_rank))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let right_value = *ranks
                    .get(usize::from(right_rank))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                Ok(i32::from(left_value) - i32::from(right_value))
            };
            let _num_trans = insertions_sort(&mut order_bytes, valence_count, width, &mut compare)?;
            for (target, bytes) in neighbor_order[..valence_count]
                .iter_mut()
                .zip(order_bytes.chunks_exact(width))
            {
                *target = AT_RANK::from_ne_bytes(
                    bytes
                        .try_into()
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                );
            }
        } else {
            neighbor_order[0] = 0;
        }

        let mut current_bonds = Vec::<u8>::new();
        for ordinal in &neighbor_order[..valence_count] {
            let bond_ordinal = usize::from(*ordinal);
            let next_atom = i32::from(
                *atom
                    .neighbor
                    .get(bond_ordinal)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
            );
            let mut bond_parity = 0_i32;
            let mut bond_parity_non_metal = 0_i32;
            if next_atom < current_index {
                let bond_stereo = i32::from(
                    *atom
                        .bond_stereo
                        .get(bond_ordinal)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let bond_type = u32::from(
                    *atom
                        .bond_type
                        .get(bond_ordinal)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?,
                );
                let bond_character = match bond_type {
                    BOND_TYPE_SINGLE => match bond_stereo {
                        value if value == STEREO_SNGL_UP as i32 => b'p',
                        value if value == -(STEREO_SNGL_UP as i32) => b'P',
                        value if value == STEREO_SNGL_DOWN as i32 => b'n',
                        value if value == -(STEREO_SNGL_DOWN as i32) => b'N',
                        value if value == STEREO_SNGL_EITHER as i32 => b'v',
                        value
                            if value == -(STEREO_SNGL_EITHER as i32)
                                && FIX_EITHER_STEREO_IN_AUX_INFO == 1 =>
                        {
                            b'V'
                        }
                        value if value == -(STEREO_SNGL_EITHER as i32) => b'v',
                        _ => b's',
                    },
                    BOND_TYPE_DOUBLE => {
                        if bond_stereo == STEREO_DBLE_EITHER as i32
                            || bond_stereo == -(STEREO_DBLE_EITHER as i32)
                        {
                            b'w'
                        } else {
                            b'd'
                        }
                    }
                    BOND_TYPE_TRIPLE => b't',
                    BOND_TYPE_ALTERN => b'a',
                    _ => b's',
                };

                let next_index =
                    usize::try_from(next_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let next = atoms
                    .get(next_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let reverse_ordinal = is_in_the_list(
                    Some(&next.neighbor),
                    current_index as u16,
                    i32::from(next.valence),
                )?
                .ok_or(SourceHeapError::PointerOutOfBounds)?;

                let mut chain_length = 0_i32;
                let mut next_endpoint = 0_i32;
                let mut next_to_current = 0_i32;
                let mut next_parity_ordinal = 0_i32;
                if atom.sb_parity[0] != 0 {
                    for parity_ordinal in 0..MAX_NUM_STEREO_BONDS as usize {
                        if atom.sb_parity[parity_ordinal] == 0 {
                            break;
                        }
                        if i32::from(atom.sb_ord[parity_ordinal])
                            == i32::try_from(bond_ordinal)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                        {
                            chain_length = get_opposite_sb_atom_slice(
                                &atoms,
                                current_index,
                                i32::try_from(bond_ordinal)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                                Some(&mut next_endpoint),
                                Some(&mut next_to_current),
                                Some(&mut next_parity_ordinal),
                            )?;
                            break;
                        }
                    }
                }

                let mut chain_length2 = 0_i32;
                let mut next_endpoint2 = 0_i32;
                let mut next_to_current2 = 0_i32;
                let mut next_parity_ordinal2 = 0_i32;
                if next.sb_parity[0] != 0 {
                    for parity_ordinal in 0..MAX_NUM_STEREO_BONDS as usize {
                        if next.sb_parity[parity_ordinal] == 0 {
                            break;
                        }
                        if i32::from(next.sb_ord[parity_ordinal])
                            == i32::try_from(reverse_ordinal)
                                .map_err(|_| SourceHeapError::SourceIntegerOverflow)?
                        {
                            chain_length2 = get_opposite_sb_atom_slice(
                                &atoms,
                                next_atom,
                                i32::try_from(reverse_ordinal)
                                    .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                                Some(&mut next_endpoint2),
                                Some(&mut next_to_current2),
                                Some(&mut next_parity_ordinal2),
                            )?;
                            break;
                        }
                    }
                }

                let mut current_endpoint = 0_i32;
                let mut current_to_next = 0_i32;
                let mut current_parity_ordinal = -1_i32;
                if (chain_length == 1 && chain_length2 == 1)
                    || (chain_length > 1 && current_index > next_endpoint)
                {
                    current_endpoint = current_index;
                    current_to_next = i32::try_from(bond_ordinal)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    for parity_ordinal in 0..MAX_NUM_STEREO_BONDS as usize {
                        let endpoint = atoms
                            .get(
                                usize::try_from(current_endpoint)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        if endpoint.sb_parity[parity_ordinal] == 0 {
                            break;
                        }
                        if i32::from(endpoint.sb_ord[parity_ordinal]) == current_to_next {
                            current_parity_ordinal = parity_ordinal as i32;
                            break;
                        }
                    }
                } else if chain_length2 > 1 && next_atom > next_endpoint2 {
                    current_endpoint = next_atom;
                    current_to_next = i32::try_from(reverse_ordinal)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
                    next_endpoint = next_endpoint2;
                    next_to_current = next_to_current2;
                    next_parity_ordinal = next_parity_ordinal2;
                    current_parity_ordinal = -1;
                    let endpoint = atoms
                        .get(
                            usize::try_from(current_endpoint)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    for parity_ordinal in 0..MAX_NUM_STEREO_BONDS as usize {
                        if endpoint.sb_parity[parity_ordinal] == 0 {
                            break;
                        }
                        if i32::from(endpoint.sb_ord[parity_ordinal]) == current_to_next {
                            current_parity_ordinal = parity_ordinal as i32;
                        }
                    }
                    chain_length = chain_length2;
                } else {
                    chain_length = 0;
                }

                if chain_length != 0 {
                    let current_endpoint_index = usize::try_from(current_endpoint)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let next_endpoint_index = usize::try_from(next_endpoint)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let current_endpoint_atom = atoms
                        .get(current_endpoint_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let next_endpoint_atom = atoms
                        .get(next_endpoint_index)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;

                    let mut p1 = 0_i32;
                    let mut p1_non_metal = 0_i32;
                    if current_parity_ordinal >= 0 {
                        let value = *current_endpoint_atom
                            .sb_parity
                            .get(
                                usize::try_from(current_parity_ordinal)
                                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                            )
                            .ok_or(SourceHeapError::PointerOutOfBounds)?;
                        p1 = parity_part(value);
                        p1_non_metal = parity_part_non_metal(value);
                    }
                    let next_parity = *next_endpoint_atom
                        .sb_parity
                        .get(
                            usize::try_from(next_parity_ordinal)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    let p2 = parity_part(next_parity);
                    let p2_non_metal = parity_part_non_metal(next_parity);
                    let well_defined = is_well_defined(p1) && is_well_defined(p2);
                    let has_metal = is_well_defined(p1_non_metal) && is_well_defined(p2_non_metal);

                    if well_defined || has_metal {
                        let mut neighbor1 = num_inp_atoms;
                        for neighbor_ordinal in 0..i32::from(current_endpoint_atom.valence).max(0) {
                            if neighbor_ordinal == current_to_next {
                                continue;
                            }
                            let neighbor = i32::from(
                                *current_endpoint_atom
                                    .neighbor
                                    .get(
                                        usize::try_from(neighbor_ordinal)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            );
                            let neighbor_atom = atoms
                                .get(
                                    usize::try_from(neighbor)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if has_metal && is_el_a_metal(i32::from(neighbor_atom.el_number))? != 0
                            {
                                continue;
                            }
                            neighbor1 = neighbor1.min(neighbor);
                        }
                        let mut neighbor_switched1 = 0_i32;
                        if neighbor1 < num_inp_atoms {
                            let parity_ordinal = usize::try_from(current_parity_ordinal)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            let stereo_neighbor_ordinal = usize::try_from(
                                *current_endpoint_atom
                                    .sn_ord
                                    .get(parity_ordinal)
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            )
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            neighbor_switched1 = i32::from(
                                neighbor1
                                    != i32::from(
                                        *current_endpoint_atom
                                            .neighbor
                                            .get(stereo_neighbor_ordinal)
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                    ),
                            );
                        } else {
                            let message = b"Cannot find 0D stereobond neighbor\0"
                                .iter()
                                .map(|byte| *byte as i8)
                                .collect::<Vec<_>>();
                            AddErrorMessage(
                                Some(&mut structure_data.pStrErrStruct),
                                Some(&message),
                            )?;
                        }

                        let mut neighbor2 = num_inp_atoms;
                        for neighbor_ordinal in 0..i32::from(next_endpoint_atom.valence).max(0) {
                            if neighbor_ordinal == next_to_current {
                                continue;
                            }
                            let neighbor = i32::from(
                                *next_endpoint_atom
                                    .neighbor
                                    .get(
                                        usize::try_from(neighbor_ordinal)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            );
                            let neighbor_atom = atoms
                                .get(
                                    usize::try_from(neighbor)
                                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                )
                                .ok_or(SourceHeapError::PointerOutOfBounds)?;
                            if has_metal && is_el_a_metal(i32::from(neighbor_atom.el_number))? != 0
                            {
                                continue;
                            }
                            neighbor2 = neighbor2.min(neighbor);
                        }
                        let mut neighbor_switched2 = 0_i32;
                        if neighbor2 < num_inp_atoms {
                            let stereo_neighbor_ordinal = usize::try_from(
                                *next_endpoint_atom
                                    .sn_ord
                                    .get(
                                        usize::try_from(next_parity_ordinal)
                                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                                    )
                                    .ok_or(SourceHeapError::PointerOutOfBounds)?,
                            )
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                            neighbor_switched2 = i32::from(
                                neighbor2
                                    != i32::from(
                                        *next_endpoint_atom
                                            .neighbor
                                            .get(stereo_neighbor_ordinal)
                                            .ok_or(SourceHeapError::PointerOutOfBounds)?,
                                    ),
                            );
                        } else {
                            let message = b"Cannot find 0D stereobond neighbor\0"
                                .iter()
                                .map(|byte| *byte as i8)
                                .collect::<Vec<_>>();
                            AddErrorMessage(
                                Some(&mut structure_data.pStrErrStruct),
                                Some(&message),
                            )?;
                        }

                        if neighbor1 < num_inp_atoms && neighbor2 < num_inp_atoms {
                            if is_well_defined(p1) && is_well_defined(p2) {
                                bond_parity =
                                    2 - (p1 + p2 + neighbor_switched1 + neighbor_switched2) % 2;
                            } else {
                                bond_parity = p1.min(p2);
                            }
                            if has_metal {
                                bond_parity_non_metal = 2
                                    - (p1_non_metal
                                        + p2_non_metal
                                        + neighbor_switched1
                                        + neighbor_switched2)
                                        % 2;
                            } else if p1_non_metal != 0 && p2_non_metal != 0 {
                                bond_parity_non_metal = p1_non_metal.min(p2_non_metal);
                            }
                        }
                    } else {
                        if p1 != 0 && p2 != 0 {
                            bond_parity = p1.min(p2);
                        }
                        if p1_non_metal != 0 && p2_non_metal != 0 {
                            bond_parity_non_metal = p1_non_metal.min(p2_non_metal);
                        }
                        if bond_parity_non_metal != 0 && bond_parity == 0 {
                            bond_parity = AB_PARITY_UNDF as i32;
                        }
                    }
                }

                current_bonds.push(bond_character);
                if let Some(symbol) = parity_symbol(bond_parity) {
                    current_bonds.push(symbol);
                }
                if let Some(symbol) = parity_symbol(bond_parity_non_metal) {
                    current_bonds.push(symbol);
                }
                current_bonds.extend(next_atom.wrapping_add(1).to_string().bytes());
            }
        }

        if current_bonds.len() > (7 * MAXVAL + 2) as usize {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        let bond_length = i32::try_from(current_bonds.len())
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if bond_length.wrapping_add(current_length).wrapping_add(2) < buffer_length {
            let start =
                usize::try_from(current_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let output = heap.slice_mut(buffer)?;
            output
                .get_mut(start..start + current_bonds.len())
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .copy_from_slice(
                    &current_bonds
                        .iter()
                        .map(|byte| *byte as i8)
                        .collect::<Vec<_>>(),
                );
            current_length = current_length.wrapping_add(bond_length);
            *output
                .get_mut(
                    usize::try_from(current_length)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                )
                .ok_or(SourceHeapError::PointerOutOfBounds)? = b';' as i8;
            current_length = current_length.wrapping_add(1);
            current_index = current_index.wrapping_add(1);
        } else {
            break;
        }
    }
    *heap
        .slice_mut(buffer)?
        .get_mut(usize::try_from(current_length).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
        .ok_or(SourceHeapError::PointerOutOfBounds)? = 0;
    *index = if num_inp_atoms > 0 { current_index } else { 0 };
    Ok(current_length)
}

#[allow(non_snake_case)]
pub(crate) fn OrigStruct_FillOut(
    heap: &mut SourceHeap,
    canon_globals: &mut CANON_GLOBALS,
    original_atom_data: &mut ORIG_ATOM_DATA,
    original_structure: &mut ORIG_STRUCT,
    structure_data: &mut STRUCT_DATA,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2933 OrigStruct_FillOut
    // INCHI✔️❌: int OrigStruct_FillOut( CANON_GLOBALS *pCG,
    // INCHI✔️❌:                        ORIG_ATOM_DATA *orig_inp_data,
    // INCHI✔️❌:                        ORIG_STRUCT    *pOrigStruct,
    // INCHI✔️❌:                        STRUCT_DATA    *sd )
    // INCHI✔️❌: {
    // INCHI✔️❌:     char szBuf[ORIG_STR_BUFLEN];
    // INCHI✔️❌:     int  i, len, len_coord, len_atoms, len_bonds;
    // INCHI✔️❌:
    // INCHI✔️❌:     pOrigStruct->polymer = NULL;
    // INCHI✔️❌:     pOrigStruct->v3000 = NULL;
    // INCHI✔️❌:
    // INCHI✔️❌:     pOrigStruct->n_zy = orig_inp_data->n_zy;
    // INCHI✔️❌:     /* Coordinates */
    // INCHI✔️❌:     len_coord = i = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (orig_inp_data->szCoord)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         while ((len = WriteOrigCoord( orig_inp_data->num_inp_atoms,
    // INCHI✔️❌:             orig_inp_data->szCoord, &i, szBuf, sizeof( szBuf ) ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:         {
    // INCHI✔️❌:             len_coord += len;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         pOrigStruct->szCoord = (char*) inchi_malloc( ( (long long)len_coord + 1 ) * sizeof( pOrigStruct->szCoord[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:         i = 0;
    // INCHI✔️❌:         if (pOrigStruct->szCoord &&
    // INCHI✔️❌:              len_coord == WriteOrigCoord( orig_inp_data->num_inp_atoms,
    // INCHI✔️❌:                  orig_inp_data->szCoord, &i, pOrigStruct->szCoord, len_coord + 1 ) &&
    // INCHI✔️❌:              i == orig_inp_data->num_inp_atoms)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* success */
    // INCHI✔️❌:             if (orig_inp_data->szCoord)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 inchi_free( orig_inp_data->szCoord );
    // INCHI✔️❌:                 orig_inp_data->szCoord = NULL;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return -1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Atoms */
    // INCHI✔️❌:     len_atoms = i = 0;
    // INCHI✔️❌:     while ((len = WriteOrigAtoms( pCG, orig_inp_data->num_inp_atoms,
    // INCHI✔️❌:         orig_inp_data->at, &i, szBuf, sizeof( szBuf ), sd ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len_atoms += len;
    // INCHI✔️❌:         if (!orig_inp_data->num_inp_atoms)
    // INCHI✔️❌:             break;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     pOrigStruct->szAtoms = (char*) inchi_malloc( ( (long long)len_atoms + 1 ) * sizeof( pOrigStruct->szAtoms[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     i = 0;
    // INCHI✔️❌:     if (pOrigStruct->szAtoms &&
    // INCHI✔️❌:          len_atoms == WriteOrigAtoms( pCG, orig_inp_data->num_inp_atoms,
    // INCHI✔️❌:              orig_inp_data->at, &i, pOrigStruct->szAtoms, len_atoms + 1, sd ) &&
    // INCHI✔️❌:          i == orig_inp_data->num_inp_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ; /* success */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Bonds */
    // INCHI✔️❌:     len_bonds = 0;
    // INCHI✔️❌:     i = 1;
    // INCHI✔️❌:     while ((len = WriteOrigBonds( pCG, orig_inp_data->num_inp_atoms,
    // INCHI✔️❌: #if ( FIX_CURE53_ISSUE_OOB_ALREADY_HAVE_THIS_MESSAGE==1 )
    // INCHI✔️❌:         orig_inp_data->at, &i, szBuf, sizeof(szBuf), sd))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌: #else
    // INCHI✔️❌:         orig_inp_data->at, &i, szBuf, sizeof(szBuf), NULL)))
    // INCHI✔️❌: #endif
    // INCHI✔️❌:     {
    // INCHI✔️❌:         len_bonds += len;
    // INCHI✔️❌:         if (!orig_inp_data->num_inp_atoms)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             break;
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     pOrigStruct->szBonds = (char*) inchi_malloc( ( (long long)len_bonds + 2 ) * sizeof( pOrigStruct->szBonds[0] ) ); /* djb-rwth: cast operator added */
    // INCHI✔️❌:     i = 1;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (pOrigStruct->szBonds &&
    // INCHI✔️❌:          len_bonds == WriteOrigBonds( pCG, orig_inp_data->num_inp_atoms,
    // INCHI✔️❌:              orig_inp_data->at, &i, pOrigStruct->szBonds, len_bonds + 2, sd ) &&
    // INCHI✔️❌:          i == orig_inp_data->num_inp_atoms)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ; /* success */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     pOrigStruct->num_atoms = orig_inp_data->num_inp_atoms;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Extensions of v. 1.05 */
    // INCHI✔️❌:     if (orig_inp_data->polymer != NULL
    // INCHI✔️❌:          && orig_inp_data->polymer->n > 0
    // INCHI✔️❌:          && orig_inp_data->valid_polymer)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pOrigStruct->polymer = orig_inp_data->polymer;
    // INCHI✔️❌:                                 /* pointer copy, do not free after use! */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     if (orig_inp_data->v3000 != NULL)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         pOrigStruct->v3000 = orig_inp_data->v3000;
    // INCHI✔️❌:                                 /* pointer copy, do not free after use! */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigStruct_FillOut
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigStruct_FillOut
    // INCHI✔️❌: #define ORIG_STR_BUFLEN (7*MAXVAL+2)
    // INCHI✔️❌: #define MAXVAL 20
    // INCHI✔️❌: #define FIX_CURE53_ISSUE_OOB_ALREADY_HAVE_THIS_MESSAGE 1
    // INCHI✔️❌: #define inchi_malloc malloc
    // INCHI✔️❌: #define inchi_free(X) do{ if(X) free(X); }while(0)
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: OrigStruct_FillOut

    let scratch = heap.allocate_model_storage(vec![0_i8; ORIG_STR_BUFLEN as usize])?;
    let result = (|| {
        original_structure.polymer = SourceMutPointer::null();
        original_structure.v3000 = SourceMutPointer::null();
        original_structure.n_zy = original_atom_data.n_zy;

        let mut index = 0_i32;
        let mut coordinate_length = 0_i32;
        if !original_atom_data.szCoord.is_null() {
            loop {
                let length = WriteOrigCoord(
                    heap,
                    original_atom_data.num_inp_atoms,
                    original_atom_data.szCoord.as_const(),
                    &mut index,
                    scratch,
                    ORIG_STR_BUFLEN as i32,
                )?;
                if length == 0 {
                    break;
                }
                coordinate_length = coordinate_length
                    .checked_add(length)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
            let allocation_length = u64::try_from(coordinate_length)
                .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
                .checked_add(1)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?;
            original_structure.szCoord = match inchi_malloc(heap, allocation_length) {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => return Err(error),
            };
            index = 0;
            if original_structure.szCoord.is_null() {
                return Ok(-1);
            }
            let written = WriteOrigCoord(
                heap,
                original_atom_data.num_inp_atoms,
                original_atom_data.szCoord.as_const(),
                &mut index,
                original_structure.szCoord,
                coordinate_length
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?,
            )?;
            if written != coordinate_length || index != original_atom_data.num_inp_atoms {
                return Ok(-1);
            }
            inchi_free(heap, original_atom_data.szCoord)?;
            original_atom_data.szCoord = SourceMutPointer::null();
        }

        index = 0;
        let mut atom_length = 0_i32;
        loop {
            let length = WriteOrigAtoms(
                heap,
                canon_globals,
                original_atom_data.num_inp_atoms,
                original_atom_data.at.as_const(),
                &mut index,
                scratch,
                ORIG_STR_BUFLEN as i32,
                structure_data,
            )?;
            if length == 0 {
                break;
            }
            atom_length = atom_length
                .checked_add(length)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            if original_atom_data.num_inp_atoms == 0 {
                break;
            }
        }
        let allocation_length = u64::try_from(atom_length)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
            .checked_add(1)
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        original_structure.szAtoms = match inchi_malloc(heap, allocation_length) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        index = 0;
        if original_structure.szAtoms.is_null() {
            return Ok(-1);
        }
        let written = WriteOrigAtoms(
            heap,
            canon_globals,
            original_atom_data.num_inp_atoms,
            original_atom_data.at.as_const(),
            &mut index,
            original_structure.szAtoms,
            atom_length
                .checked_add(1)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?,
            structure_data,
        )?;
        if written != atom_length || index != original_atom_data.num_inp_atoms {
            return Ok(-1);
        }

        let mut bond_length = 0_i32;
        index = 1;
        loop {
            let length = WriteOrigBonds(
                heap,
                canon_globals,
                original_atom_data.num_inp_atoms,
                original_atom_data.at,
                &mut index,
                scratch,
                ORIG_STR_BUFLEN as i32,
                structure_data,
            )?;
            if length == 0 {
                break;
            }
            bond_length = bond_length
                .checked_add(length)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            if original_atom_data.num_inp_atoms == 0 {
                break;
            }
        }
        let allocation_length = u64::try_from(bond_length)
            .map_err(|_| SourceHeapError::AllocationElementCountOutOfRange)?
            .checked_add(2)
            .ok_or(SourceHeapError::AllocationSizeOverflow)?;
        original_structure.szBonds = match inchi_malloc(heap, allocation_length) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => return Err(error),
        };
        index = 1;
        if original_structure.szBonds.is_null() {
            return Ok(-1);
        }
        let written = WriteOrigBonds(
            heap,
            canon_globals,
            original_atom_data.num_inp_atoms,
            original_atom_data.at,
            &mut index,
            original_structure.szBonds,
            bond_length
                .checked_add(2)
                .ok_or(SourceHeapError::SourceIntegerOverflow)?,
            structure_data,
        )?;
        if written != bond_length || index != original_atom_data.num_inp_atoms {
            return Ok(-1);
        }

        original_structure.num_atoms = original_atom_data.num_inp_atoms;
        if !original_atom_data.polymer.is_null()
            && heap
                .slice(original_atom_data.polymer.as_const())?
                .first()
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .n
                > 0
            && original_atom_data.valid_polymer != 0
        {
            original_structure.polymer = original_atom_data.polymer;
        }
        if !original_atom_data.v3000.is_null() {
            original_structure.v3000 = original_atom_data.v3000;
        }
        Ok(0)
    })();
    let free_result = heap.free(scratch);
    match result {
        Err(error) => Err(error),
        Ok(value) => {
            free_result?;
            Ok(value)
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn OrigStruct_Free(
    heap: &mut SourceHeap,
    original_structure: Option<&mut ORIG_STRUCT>,
) -> Result<(), SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3051 OrigStruct_Free
    // INCHI✔️❌: void OrigStruct_Free( ORIG_STRUCT *pOrigStruct )
    // INCHI✔️❌: {
    // INCHI✔️❌:     if (pOrigStruct)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (pOrigStruct->szAtoms)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( pOrigStruct->szAtoms );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (pOrigStruct->szBonds)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( pOrigStruct->szBonds );
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (pOrigStruct->szCoord)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             inchi_free( pOrigStruct->szCoord );
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* For
    // INCHI✔️❌:
    // INCHI✔️❌:             OAD_Polymer *polymer;
    // INCHI✔️❌:             OAD_V3000    *v3000;
    // INCHI✔️❌:
    // INCHI✔️❌:             we used shallow (pointer) copy of analogs from orig_inp_data, so do not free these here */
    // INCHI✔️❌:
    // INCHI✔️❌:         /*memset( pOrigStruct, 0, sizeof(*pOrigStruct) );*/
    // INCHI✔️❌:         pOrigStruct->szAtoms = NULL;
    // INCHI✔️❌:         pOrigStruct->szBonds = NULL;
    // INCHI✔️❌:         pOrigStruct->szCoord = NULL;
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: OrigStruct_Free

    let Some(original_structure) = original_structure else {
        return Ok(());
    };
    if !original_structure.szAtoms.is_null() {
        inchi_free(heap, original_structure.szAtoms)?;
    }
    if !original_structure.szBonds.is_null() {
        inchi_free(heap, original_structure.szBonds)?;
    }
    if !original_structure.szCoord.is_null() {
        inchi_free(heap, original_structure.szCoord)?;
    }
    original_structure.szAtoms = SourceMutPointer::null();
    original_structure.szBonds = SourceMutPointer::null();
    original_structure.szCoord = SourceMutPointer::null();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn original_coordinates(fields: [&[u8]; 3]) -> MOL_COORD {
        let mut coordinates = [b' ' as i8; 32];
        for (index, field) in fields.into_iter().enumerate() {
            assert!(field.len() <= LEN_COORD as usize);
            for (target, source) in coordinates
                [index * LEN_COORD as usize..index * LEN_COORD as usize + field.len()]
                .iter_mut()
                .zip(field.iter().copied())
            {
                *target = source as i8;
            }
        }
        coordinates
    }

    fn coordinate_bytes(coordinates: &MOL_COORD) -> Vec<u8> {
        coordinates
            .iter()
            .take_while(|byte| **byte != 0)
            .map(|byte| *byte as u8)
            .collect()
    }

    #[test]
    fn source_port__ichiprt1__cleanorigcoord__line_2228() {
        let mut heap = SourceHeap::default();

        let mut decimal = original_coordinates([b" 0.000000", b"001.2300", b"-0002.500"]);
        assert_eq!(CleanOrigCoord(&mut heap, &mut decimal, b',' as i32), Ok(1));
        assert_eq!(coordinate_bytes(&decimal), b"0,1.23,-2.5");
        assert!(
            decimal[b"0,1.23,-2.5".len()..]
                .iter()
                .all(|byte| *byte == 0)
        );

        let mut threshold = original_coordinates([b"0.000001", b"0.0000009", b"-0.000000"]);
        assert_eq!(
            CleanOrigCoord(&mut heap, &mut threshold, b'|' as i32),
            Ok(2)
        );
        assert_eq!(coordinate_bytes(&threshold), b".000001|0|0");

        let mut exponent = original_coordinates([b"01.200e+03", b"02.500D-02", b"1.000E+000"]);
        assert_eq!(CleanOrigCoord(&mut heap, &mut exponent, b',' as i32), Ok(0));
        assert_eq!(coordinate_bytes(&exponent), b"1.2e3,2.5D-2,1");

        let mut special = original_coordinates([b"nan", b"+inf", b"1.0e+"]);
        assert_eq!(CleanOrigCoord(&mut heap, &mut special, b';' as i32), Ok(0));
        assert_eq!(coordinate_bytes(&special), b"nan;+inf;1");

        let mut invalid = original_coordinates([b"nonnumeric", b"\t\r\n", b"-0"]);
        assert_eq!(CleanOrigCoord(&mut heap, &mut invalid, b',' as i32), Ok(3));
        assert_eq!(coordinate_bytes(&invalid), b"0,0,0");

        let mut full = original_coordinates([b"1234567890", b"1234567890", b"1234567890"]);
        assert_eq!(CleanOrigCoord(&mut heap, &mut full, b',' as i32), Ok(0));
        assert_eq!(
            full.map(|byte| byte as u8),
            *b"1234567890,1234567890,1234567890"
        );
    }

    #[test]
    fn source_port__ichiprt1__writeorigcoord__line_2340() {
        let mut heap = SourceHeap::default();
        let coordinates = heap
            .allocate_model_storage(vec![
                original_coordinates([b"1", b"2", b"3"]),
                original_coordinates([b"0", b"-0", b"0.0000001"]),
                original_coordinates([b"1234567890", b"1234567890", b"1234567890"]),
            ])
            .unwrap();

        let too_small = heap.allocate_model_storage(vec![0x55_i8; 8]).unwrap();
        let mut index = 0_i32;
        assert_eq!(
            WriteOrigCoord(
                &mut heap,
                3,
                coordinates.as_const(),
                &mut index,
                too_small,
                6,
            ),
            Ok(0)
        );
        assert_eq!(index, 0);
        assert_eq!(&heap.slice(too_small.as_const()).unwrap()[..2], &[0, 0x55]);

        let first_chunk = heap.allocate_model_storage(vec![0x55_i8; 8]).unwrap();
        assert_eq!(
            WriteOrigCoord(
                &mut heap,
                3,
                coordinates.as_const(),
                &mut index,
                first_chunk,
                7,
            ),
            Ok(6)
        );
        assert_eq!(index, 1);
        assert_eq!(
            &heap.slice(first_chunk.as_const()).unwrap()[..8],
            &[
                b'1' as i8, b',' as i8, b'2' as i8, b',' as i8, b'3' as i8, b';' as i8, 0, 0x55
            ]
        );

        let zero_chunk = heap.allocate_model_storage(vec![0x55_i8; 4]).unwrap();
        assert_eq!(
            WriteOrigCoord(
                &mut heap,
                3,
                coordinates.as_const(),
                &mut index,
                zero_chunk,
                2,
            ),
            Ok(1)
        );
        assert_eq!(index, 2);
        assert_eq!(
            &heap.slice(zero_chunk.as_const()).unwrap()[..3],
            &[b';' as i8, 0, 0x55]
        );

        let full_chunk = heap.allocate_model_storage(vec![0x55_i8; 35]).unwrap();
        assert_eq!(
            WriteOrigCoord(
                &mut heap,
                3,
                coordinates.as_const(),
                &mut index,
                full_chunk,
                34,
            ),
            Ok(33)
        );
        assert_eq!(index, 3);
        assert_eq!(
            heap.slice(full_chunk.as_const()).unwrap()[..34]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            b"1234567890,1234567890,1234567890;\0"
        );
        assert_eq!(heap.slice(full_chunk.as_const()).unwrap()[34], 0x55);

        let empty = heap.allocate_model_storage(vec![0x55_i8; 2]).unwrap();
        index = 0;
        assert_eq!(
            WriteOrigCoord(
                &mut heap,
                0,
                SourceConstPointer::null(),
                &mut index,
                empty,
                0,
            ),
            Ok(0)
        );
        assert_eq!(index, 0);
        assert_eq!(&heap.slice(empty.as_const()).unwrap()[..2], &[0, 0x55]);
    }

    fn original_atom(
        symbol: &[u8],
        element_number: u8,
        bond_types: &[u8],
        chemical_bonds_valence: i8,
        hydrogens: i8,
    ) -> inp_ATOM {
        let mut atom = inp_ATOM {
            el_number: element_number,
            valence: bond_types.len() as i8,
            chem_bonds_valence: chemical_bonds_valence,
            num_H: hydrogens,
            ..inp_ATOM::default()
        };
        for (target, source) in atom
            .elname
            .iter_mut()
            .zip(symbol.iter().copied().chain(std::iter::once(0)))
        {
            *target = source as i8;
        }
        atom.bond_type[..bond_types.len()].copy_from_slice(bond_types);
        atom
    }

    fn write_single_original_atom(
        atom: inp_ATOM,
        chiral_flag: i32,
        buffer_length: i32,
    ) -> (Result<i32, SourceHeapError>, i32, Vec<u8>) {
        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
        let buffer = heap.allocate_model_storage(vec![0x55_i8; 64]).unwrap();
        let mut index = 0_i32;
        let result = WriteOrigAtoms(
            &mut heap,
            &mut CANON_GLOBALS::default(),
            1,
            atoms.as_const(),
            &mut index,
            buffer,
            buffer_length,
            &STRUCT_DATA {
                bChiralFlag: chiral_flag,
                ..STRUCT_DATA::default()
            },
        );
        let bytes = heap
            .slice(buffer.as_const())
            .unwrap()
            .iter()
            .take_while(|byte| **byte != 0)
            .map(|byte| *byte as u8)
            .collect();
        (result, index, bytes)
    }

    #[test]
    fn source_port__ichiprt1__writeorigatoms__line_2410() {
        let carbon = original_atom(b"C", 6, &[2, 2], 4, 0);
        assert_eq!(
            write_single_original_atom(
                carbon.clone(),
                (FLAG_INP_AT_CHIRAL | FLAG_INP_AT_NONCHIRAL) as i32,
                64,
            ),
            (Ok(3), 1, b"1cC".to_vec())
        );
        assert_eq!(
            write_single_original_atom(carbon.clone(), FLAG_INP_AT_NONCHIRAL as i32, 64,),
            (Ok(3), 1, b"1nC".to_vec())
        );
        assert_eq!(
            write_single_original_atom(carbon.clone(), 0, 2),
            (Ok(1), 0, b"1".to_vec())
        );

        let zero_valence = original_atom(b"C", 6, &[], 0, 0);
        assert_eq!(
            write_single_original_atom(zero_valence, 0, 64),
            (Ok(3), 1, b"1C0".to_vec())
        );

        let mut charged = original_atom(b"C", 6, &[2], 2, 1);
        charged.charge = 3;
        assert_eq!(
            write_single_original_atom(charged, 0, 64),
            (Ok(5), 1, b"1C3+3".to_vec())
        );

        let mut radical = original_atom(b"C", 6, &[3], 3, 0);
        radical.radical = RADICAL_DOUBLET as i8;
        assert_eq!(
            write_single_original_atom(radical, 0, 64),
            (Ok(4), 1, b"1C.2".to_vec())
        );

        let mut isotope = carbon.clone();
        isotope.iso_atw_diff = 1;
        assert_eq!(
            write_single_original_atom(isotope, 0, 64),
            (Ok(6), 1, b"1C.i12".to_vec())
        );
        let mut shifted_isotope = carbon.clone();
        shifted_isotope.iso_atw_diff = -2;
        assert_eq!(
            write_single_original_atom(shifted_isotope, 0, 64),
            (Ok(6), 1, b"1C.i10".to_vec())
        );

        let mut isotope_hydrogens = carbon.clone();
        isotope_hydrogens.num_iso_H = [1, 2, 3];
        assert_eq!(
            write_single_original_atom(isotope_hydrogens, 0, 64),
            (Ok(8), 1, b"1C.hd2t3".to_vec())
        );

        let mut parity_atoms = vec![
            original_atom(b"C", 6, &[1, 1, 1], 3, 1),
            original_atom(b"C", 6, &[1], 1, 3),
            original_atom(b"C", 6, &[1], 1, 3),
            original_atom(b"C", 6, &[1], 1, 3),
        ];
        for (atom_number, atom) in parity_atoms.iter_mut().enumerate() {
            atom.orig_at_number = atom_number as u16 + 1;
        }
        parity_atoms[0].neighbor[..3].copy_from_slice(&[1, 2, 3]);
        parity_atoms[0].p_orig_at_num = [4, 2, 1, 3];
        parity_atoms[0].p_parity = AB_PARITY_ODD as i8;

        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(parity_atoms.clone()).unwrap();
        let buffer = heap.allocate_model_storage(vec![0x55_i8; 64]).unwrap();
        let mut index = 0_i32;
        assert_eq!(
            WriteOrigAtoms(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                4,
                atoms.as_const(),
                &mut index,
                buffer,
                64,
                &STRUCT_DATA::default(),
            ),
            Ok(7)
        );
        assert_eq!(index, 4);
        assert_eq!(
            heap.slice(buffer.as_const()).unwrap()[..8]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            b"4C.oCCC\0"
        );

        parity_atoms[0].p_parity = AB_PARITY_UNKN as i8;
        let atoms = heap.allocate_model_storage(parity_atoms.clone()).unwrap();
        let buffer = heap.allocate_model_storage(vec![0x55_i8; 64]).unwrap();
        index = 0;
        assert_eq!(
            WriteOrigAtoms(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                4,
                atoms.as_const(),
                &mut index,
                buffer,
                64,
                &STRUCT_DATA::default(),
            ),
            Ok(7)
        );
        assert_eq!(
            heap.slice(buffer.as_const()).unwrap()[..8]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            b"4C.uCCC\0"
        );

        parity_atoms[0].p_parity = AB_PARITY_IISO as i8;
        let atoms = heap.allocate_model_storage(parity_atoms).unwrap();
        let buffer = heap.allocate_model_storage(vec![0x55_i8; 64]).unwrap();
        index = 0;
        assert_eq!(
            WriteOrigAtoms(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                4,
                atoms.as_const(),
                &mut index,
                buffer,
                64,
                &STRUCT_DATA::default(),
            ),
            Ok(5)
        );
        assert_eq!(
            heap.slice(buffer.as_const()).unwrap()[..6]
                .iter()
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>(),
            b"4CCCC\0"
        );
    }

    #[test]
    fn source_port__ichiprt1__writeorigbonds__line_2610() {
        fn run(
            atoms: Vec<inp_ATOM>,
            buffer_length: i32,
        ) -> (Result<i32, SourceHeapError>, i32, Vec<u8>, String, bool) {
            let atom_count = atoms.len() as i32;
            let mut heap = SourceHeap::default();
            let atom_pointer = if atoms.is_empty() {
                SourceMutPointer::null()
            } else {
                heap.allocate_model_storage(atoms).unwrap()
            };
            let buffer = heap.allocate_model_storage(vec![0x55_i8; 512]).unwrap();
            let mut index = 0_i32;
            let mut canon_globals = CANON_GLOBALS::default();
            let mut structure_data = STRUCT_DATA::default();
            let result = WriteOrigBonds(
                &mut heap,
                &mut canon_globals,
                atom_count,
                atom_pointer,
                &mut index,
                buffer,
                buffer_length,
                &mut structure_data,
            );
            let output = heap
                .slice(buffer.as_const())
                .unwrap()
                .iter()
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>();
            let error = structure_data
                .pStrErrStruct
                .iter()
                .take_while(|byte| **byte != 0)
                .map(|byte| *byte as u8)
                .collect::<Vec<_>>();
            (
                result,
                index,
                output,
                String::from_utf8(error).unwrap(),
                !canon_globals.m_pn_RankForSort.is_null(),
            )
        }

        assert_eq!(
            run(Vec::new(), 512),
            (Ok(0), 0, Vec::new(), String::new(), false)
        );

        let cases = [
            (BOND_TYPE_SINGLE as u8, STEREO_SNGL_UP as i8, b'p'),
            (BOND_TYPE_SINGLE as u8, -(STEREO_SNGL_UP as i8), b'P'),
            (BOND_TYPE_SINGLE as u8, STEREO_SNGL_DOWN as i8, b'n'),
            (BOND_TYPE_SINGLE as u8, -(STEREO_SNGL_DOWN as i8), b'N'),
            (BOND_TYPE_SINGLE as u8, STEREO_SNGL_EITHER as i8, b'v'),
            (BOND_TYPE_SINGLE as u8, -(STEREO_SNGL_EITHER as i8), b'V'),
            (BOND_TYPE_DOUBLE as u8, STEREO_DBLE_EITHER as i8, b'w'),
            (BOND_TYPE_DOUBLE as u8, 0, b'd'),
            (BOND_TYPE_TRIPLE as u8, 0, b't'),
            (BOND_TYPE_ALTERN as u8, 0, b'a'),
            (99, 0, b's'),
        ];
        let mut all_codes = vec![inp_ATOM::default(); cases.len() + 1];
        all_codes[0].valence = cases.len() as i8;
        for (ordinal, (bond_type, stereo, _)) in cases.iter().copied().enumerate() {
            let leaf = ordinal + 1;
            all_codes[0].neighbor[ordinal] = leaf as u16;
            all_codes[0].bond_type[ordinal] = bond_type;
            all_codes[0].bond_stereo[ordinal] = -stereo;
            all_codes[leaf].valence = 1;
            all_codes[leaf].neighbor[0] = 0;
            all_codes[leaf].bond_type[0] = bond_type;
            all_codes[leaf].bond_stereo[0] = stereo;
        }
        let mut expected = vec![b';'];
        for (_, _, code) in cases {
            expected.extend([code, b'1', b';']);
        }
        assert_eq!(
            run(all_codes, 512),
            (Ok(expected.len() as i32), 12, expected, String::new(), true,)
        );

        let mut sorted = vec![inp_ATOM::default(); 3];
        sorted[0].valence = 1;
        sorted[0].neighbor[0] = 2;
        sorted[0].bond_type[0] = BOND_TYPE_TRIPLE as u8;
        sorted[1].valence = 1;
        sorted[1].neighbor[0] = 2;
        sorted[1].bond_type[0] = BOND_TYPE_DOUBLE as u8;
        sorted[2].valence = 2;
        sorted[2].neighbor[..2].copy_from_slice(&[1, 0]);
        sorted[2].bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_TYPE_TRIPLE as u8]);
        assert_eq!(
            run(sorted, 512),
            (Ok(7), 3, b";;t1d2;".to_vec(), String::new(), true,)
        );

        assert_eq!(
            run(vec![inp_ATOM::default()], 3),
            (Ok(1), 1, b";".to_vec(), String::new(), false)
        );
        assert_eq!(
            run(vec![inp_ATOM::default()], 2),
            (Ok(0), 0, Vec::new(), String::new(), false)
        );

        fn stereo_fixture(
            first_parity: i8,
            second_parity: i8,
            substituent_element: u8,
        ) -> Vec<inp_ATOM> {
            let mut atoms = vec![inp_ATOM::default(); 4];
            atoms[0].valence = 2;
            atoms[0].neighbor[..2].copy_from_slice(&[1, 2]);
            atoms[0].bond_type[..2]
                .copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_TYPE_SINGLE as u8]);
            atoms[0].sb_ord[0] = 0;
            atoms[0].sn_ord[0] = 1;
            atoms[0].sb_parity[0] = first_parity;

            atoms[1].valence = 2;
            atoms[1].neighbor[..2].copy_from_slice(&[0, 3]);
            atoms[1].bond_type[..2]
                .copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_TYPE_SINGLE as u8]);
            atoms[1].sb_ord[0] = 0;
            atoms[1].sn_ord[0] = 1;
            atoms[1].sb_parity[0] = second_parity;

            atoms[2].valence = 1;
            atoms[2].neighbor[0] = 0;
            atoms[2].bond_type[0] = BOND_TYPE_SINGLE as u8;
            atoms[2].el_number = substituent_element;
            atoms[3].valence = 1;
            atoms[3].neighbor[0] = 1;
            atoms[3].bond_type[0] = BOND_TYPE_SINGLE as u8;
            atoms[3].el_number = substituent_element;
            atoms
        }

        assert_eq!(
            run(
                stereo_fixture(AB_PARITY_ODD as i8, AB_PARITY_ODD as i8, 6),
                512,
            ),
            (Ok(11), 4, b";d+1;s1;s2;".to_vec(), String::new(), true,)
        );
        assert_eq!(
            run(
                stereo_fixture(AB_PARITY_UNKN as i8, AB_PARITY_UNDF as i8, 6),
                512,
            ),
            (Ok(11), 4, b";du1;s1;s2;".to_vec(), String::new(), true,)
        );
        let connected_and_non_metal = (AB_PARITY_ODD | (AB_PARITY_EVEN << SB_PARITY_SHFT)) as i8;
        assert_eq!(
            run(
                stereo_fixture(connected_and_non_metal, connected_and_non_metal, 6,),
                512,
            ),
            (Ok(12), 4, b";d++1;s1;s2;".to_vec(), String::new(), true,)
        );
        let non_metal_only = (AB_PARITY_ODD << SB_PARITY_SHFT) as i8;
        assert_eq!(
            run(stereo_fixture(non_metal_only, non_metal_only, 26), 512,),
            (
                Ok(10),
                4,
                b";d1;s1;s2;".to_vec(),
                "Cannot find 0D stereobond neighbor".to_owned(),
                true,
            )
        );
    }

    #[test]
    fn source_port__ichiprt1__origstruct_fillout__line_2933() {
        fn text(heap: &SourceHeap, pointer: SourceMutPointer<i8>) -> Vec<u8> {
            let bytes = heap.slice(pointer.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..=end].iter().map(|byte| *byte as u8).collect()
        }

        fn input(heap: &mut SourceHeap, with_coordinates: bool) -> ORIG_ATOM_DATA {
            let atoms = heap
                .allocate_model_storage(vec![original_atom(b"C", 6, &[], 0, 0)])
                .unwrap();
            let coordinates = if with_coordinates {
                heap.allocate_model_storage(vec![original_coordinates([b"1", b"2", b"3"])])
                    .unwrap()
            } else {
                SourceMutPointer::null()
            };
            ORIG_ATOM_DATA {
                at: atoms,
                num_inp_atoms: 1,
                szCoord: coordinates,
                n_zy: 9,
                ..ORIG_ATOM_DATA::default()
            }
        }

        let mut heap = SourceHeap::default();
        let polymer = heap
            .allocate_model_storage(vec![crate::source_types::OAD_Polymer {
                n: 1,
                ..crate::source_types::OAD_Polymer::default()
            }])
            .unwrap();
        let v3000 = heap
            .allocate_model_storage(vec![crate::source_types::OAD_V3000::default()])
            .unwrap();
        let mut original_atom_data = input(&mut heap, true);
        original_atom_data.polymer = polymer;
        original_atom_data.valid_polymer = 1;
        original_atom_data.v3000 = v3000;
        let old_coordinates = original_atom_data.szCoord;
        let mut original_structure = ORIG_STRUCT::default();
        heap.trace_source_allocations();
        assert_eq!(
            OrigStruct_FillOut(
                &mut heap,
                &mut CANON_GLOBALS::default(),
                &mut original_atom_data,
                &mut original_structure,
                &mut STRUCT_DATA::default(),
            ),
            Ok(0)
        );
        assert_eq!(heap.source_allocation_calls(), 3);
        assert!(original_atom_data.szCoord.is_null());
        assert_eq!(
            heap.slice(old_coordinates.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(original_structure.num_atoms, 1);
        assert_eq!(original_structure.n_zy, 9);
        assert_eq!(original_structure.polymer, polymer);
        assert_eq!(original_structure.v3000, v3000);
        assert_eq!(text(&heap, original_structure.szCoord), b"1,2,3;\0");
        assert_eq!(text(&heap, original_structure.szAtoms), b"1C0\0");
        assert_eq!(text(&heap, original_structure.szBonds), b"\0");
        OrigStruct_Free(&mut heap, Some(&mut original_structure)).unwrap();
        assert_eq!(heap.slice(polymer.as_const()).unwrap()[0].n, 1);
        assert_eq!(heap.slice(v3000.as_const()).unwrap().len(), 1);

        let mut coordinate_failure_heap = SourceHeap::default();
        let mut coordinate_failure_input = input(&mut coordinate_failure_heap, true);
        let retained_coordinates = coordinate_failure_input.szCoord;
        let mut coordinate_failure_output = ORIG_STRUCT {
            num_atoms: 77,
            ..ORIG_STRUCT::default()
        };
        coordinate_failure_heap.fail_after_allocations(0);
        assert_eq!(
            OrigStruct_FillOut(
                &mut coordinate_failure_heap,
                &mut CANON_GLOBALS::default(),
                &mut coordinate_failure_input,
                &mut coordinate_failure_output,
                &mut STRUCT_DATA::default(),
            ),
            Ok(-1)
        );
        assert!(coordinate_failure_output.szCoord.is_null());
        assert_eq!(coordinate_failure_output.num_atoms, 77);
        assert_eq!(coordinate_failure_output.n_zy, 9);
        assert_eq!(coordinate_failure_input.szCoord, retained_coordinates);
        assert!(
            coordinate_failure_heap
                .slice(retained_coordinates.as_const())
                .is_ok()
        );

        let mut atom_failure_heap = SourceHeap::default();
        let mut atom_failure_input = input(&mut atom_failure_heap, true);
        let transferred_coordinates = atom_failure_input.szCoord;
        let mut atom_failure_output = ORIG_STRUCT {
            num_atoms: 77,
            ..ORIG_STRUCT::default()
        };
        atom_failure_heap.fail_after_allocations(1);
        assert_eq!(
            OrigStruct_FillOut(
                &mut atom_failure_heap,
                &mut CANON_GLOBALS::default(),
                &mut atom_failure_input,
                &mut atom_failure_output,
                &mut STRUCT_DATA::default(),
            ),
            Ok(-1)
        );
        assert!(atom_failure_input.szCoord.is_null());
        assert_eq!(
            atom_failure_heap.slice(transferred_coordinates.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            text(&atom_failure_heap, atom_failure_output.szCoord),
            b"1,2,3;\0"
        );
        assert!(atom_failure_output.szAtoms.is_null());
        assert_eq!(atom_failure_output.num_atoms, 77);
        OrigStruct_Free(&mut atom_failure_heap, Some(&mut atom_failure_output)).unwrap();

        let mut bond_failure_heap = SourceHeap::default();
        let mut bond_failure_input = input(&mut bond_failure_heap, false);
        let mut bond_failure_output = ORIG_STRUCT {
            num_atoms: 77,
            ..ORIG_STRUCT::default()
        };
        bond_failure_heap.fail_after_allocations(1);
        assert_eq!(
            OrigStruct_FillOut(
                &mut bond_failure_heap,
                &mut CANON_GLOBALS::default(),
                &mut bond_failure_input,
                &mut bond_failure_output,
                &mut STRUCT_DATA::default(),
            ),
            Ok(-1)
        );
        assert_eq!(
            text(&bond_failure_heap, bond_failure_output.szAtoms),
            b"1C0\0"
        );
        assert!(bond_failure_output.szBonds.is_null());
        assert_eq!(bond_failure_output.num_atoms, 77);
        OrigStruct_Free(&mut bond_failure_heap, Some(&mut bond_failure_output)).unwrap();

        let mut empty_heap = SourceHeap::default();
        let mut empty_input = ORIG_ATOM_DATA::default();
        let mut empty_output = ORIG_STRUCT::default();
        assert_eq!(
            OrigStruct_FillOut(
                &mut empty_heap,
                &mut CANON_GLOBALS::default(),
                &mut empty_input,
                &mut empty_output,
                &mut STRUCT_DATA::default(),
            ),
            Ok(0)
        );
        assert_eq!(text(&empty_heap, empty_output.szAtoms), b"0\0");
        assert_eq!(text(&empty_heap, empty_output.szBonds), b"\0");
        assert_eq!(empty_output.num_atoms, 0);
    }

    #[test]
    fn source_port__ichiprt1__origstruct_free__line_3051() {
        let mut heap = SourceHeap::default();
        assert_eq!(OrigStruct_Free(&mut heap, None), Ok(()));

        let mut empty = ORIG_STRUCT {
            num_atoms: 7,
            n_zy: 9,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(OrigStruct_Free(&mut heap, Some(&mut empty)), Ok(()));
        assert_eq!(empty.num_atoms, 7);
        assert_eq!(empty.n_zy, 9);

        let atoms = heap.allocate_model_storage(vec![1_i8, 2, 0]).unwrap();
        let bonds = heap.allocate_model_storage(vec![3_i8, 4, 0]).unwrap();
        let coordinates = heap.allocate_model_storage(vec![5_i8, 6, 0]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![crate::source_types::OAD_Polymer {
                n: 11,
                ..crate::source_types::OAD_Polymer::default()
            }])
            .unwrap();
        let v3000 = heap
            .allocate_model_storage(vec![crate::source_types::OAD_V3000::default()])
            .unwrap();
        let mut structure = ORIG_STRUCT {
            num_atoms: 13,
            szAtoms: atoms,
            szBonds: bonds,
            szCoord: coordinates,
            polymer,
            v3000,
            n_zy: 17,
        };
        assert_eq!(OrigStruct_Free(&mut heap, Some(&mut structure)), Ok(()));
        assert!(structure.szAtoms.is_null());
        assert!(structure.szBonds.is_null());
        assert!(structure.szCoord.is_null());
        assert_eq!(structure.num_atoms, 13);
        assert_eq!(structure.polymer, polymer);
        assert_eq!(structure.v3000, v3000);
        assert_eq!(structure.n_zy, 17);
        assert_eq!(
            heap.slice(atoms.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(bonds.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(
            heap.slice(coordinates.as_const()),
            Err(SourceHeapError::MissingAllocation)
        );
        assert_eq!(heap.slice(polymer.as_const()).unwrap()[0].n, 11);
        assert_eq!(heap.slice(v3000.as_const()).unwrap().len(), 1);

        assert_eq!(OrigStruct_Free(&mut heap, Some(&mut structure)), Ok(()));
        assert_eq!(heap.slice(polymer.as_const()).unwrap()[0].n, 11);
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_chargesradicalsandunusualvalences__line_4717() {
        fn buffer(heap: &mut SourceHeap, size: usize) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; size]).unwrap(),
                nAllocatedLength: size as i32,
                nPtr: 32,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn output_text(heap: &SourceHeap, output: &INCHI_IOSTREAM) -> String {
            let bytes = heap.slice(output.s.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|byte| *byte as u8 as char)
                .collect()
        }
        fn clear(heap: &mut SourceHeap, output: &mut INCHI_IOSTREAM) {
            output.s.nUsedLength = 0;
            heap.slice_mut(output.s.pStr).unwrap()[0] = 0;
        }

        let mut heap = SourceHeap::default();
        let line_feed = heap.allocate_model_storage(vec![b'\n' as i8, 0]).unwrap();
        let original = heap
            .allocate_model_storage(vec![crate::source_types::ORIG_INFO {
                cCharge: 1,
                cRadical: 0,
                cUnusualValence: 0,
            }])
            .unwrap();
        let aux = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Aux {
                nNumberOfAtoms: 1,
                OrigInfo: original,
                ..crate::source_types::INChI_Aux::default()
            }])
            .unwrap();
        let inchi = heap
            .allocate_model_storage(vec![crate::source_types::INChI {
                nNumberOfAtoms: 1,
                ..crate::source_types::INChI::default()
            }])
            .unwrap();
        let sort = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        let mut string_buffer = buffer(&mut heap, 256);
        let mut output = INCHI_IOSTREAM {
            s: buffer(&mut heap, 256),
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let mut inactive = INCHI_OUT_CTL::default();
        assert_eq!(
            OutputAUXINFO_ChargesRadicalsAndUnusualValences(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut inactive,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null()
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "");

        let mut control = INCHI_OUT_CTL {
            bChargesRadVal: [1, 0],
            nTag: 2,
            bPlainTextTags: 1,
            bOutType: crate::source_types::OUT_T1 as i32,
            num_components: 1,
            pINChISort: sort,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_ChargesRadicalsAndUnusualValences(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null()
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "/CRV:1+1\n");
        assert_eq!(control.tot_len, 3);
        assert_eq!(
            control.bTag1,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_CRV_ as i32
        );

        clear(&mut heap, &mut output);
        control.nTag = 3;
        control.bFhTag = crate::source_types::local_ichiprt1::tagAuxLblBit_AL_FIXH as i32;
        assert_eq!(
            OutputAUXINFO_ChargesRadicalsAndUnusualValences(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null()
            ),
            Ok(0)
        );
        assert_eq!(
            output_text(&heap, &output),
            "/CRV:{fixed_H:charge_radical_valence}1+1\n"
        );

        clear(&mut heap, &mut output);
        control.bSecondNonTautPass = 1;
        assert_eq!(
            OutputAUXINFO_ChargesRadicalsAndUnusualValences(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null()
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "");

        control.bSecondNonTautPass = 0;
        control.bOverflow = 1;
        assert_eq!(
            OutputAUXINFO_ChargesRadicalsAndUnusualValences(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null()
            ),
            Ok(1)
        );
        assert_eq!(output_text(&heap, &output), "");
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_reversibilityinfo__line_4753() {
        fn buffer(heap: &mut SourceHeap, size: usize) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; size]).unwrap(),
                nAllocatedLength: size as i32,
                nPtr: size as i32,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn c_text(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .map(|byte| byte as i8)
                    .chain(std::iter::once(0))
                    .collect(),
            )
            .unwrap()
        }
        fn output_text(heap: &SourceHeap, output: &INCHI_IOSTREAM) -> String {
            let bytes = heap.slice(output.s.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|byte| *byte as u8 as char)
                .collect()
        }
        fn run(
            heap: &mut SourceHeap,
            structure: Option<&ORIG_STRUCT>,
            size: usize,
            plain: i32,
            second: i32,
        ) -> (Result<i32, SourceHeapError>, String, INCHI_OUT_CTL) {
            let mut string_buffer = buffer(heap, size);
            let mut output = INCHI_IOSTREAM {
                s: buffer(heap, 512),
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            };
            let line_feed = c_text(heap, "\n");
            let mut io = INCHI_OUT_CTL {
                nTag: 2,
                bPlainTextTags: plain,
                bSecondNonTautPass: second,
                ..INCHI_OUT_CTL::default()
            };
            let result = OutputAUXINFO_ReversibilityInfo(
                heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                structure,
                &mut io,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            );
            (result, output_text(heap, &output), io)
        }

        let mut heap = SourceHeap::default();
        assert_eq!(run(&mut heap, None, 80, 1, 0).0, Ok(0));
        assert_eq!(run(&mut heap, None, 80, 1, 0).1, "");

        let atoms = c_text(&mut heap, "CZyO");
        let bonds = c_text(&mut heap, "1-2;");
        let coords = c_text(&mut heap, "0,0;");
        let structure = ORIG_STRUCT {
            num_atoms: 2,
            szAtoms: atoms,
            szBonds: bonds,
            szCoord: coords,
            ..ORIG_STRUCT::default()
        };
        let (result, output, io) = run(&mut heap, Some(&structure), 80, 1, 0);
        assert_eq!(result, Ok(0));
        assert_eq!(output, "/rA:CZzO\n/rB:1-2;\n/rC:0,0;\n");
        assert_eq!(
            io.bTag1,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_REVR as i32
        );
        assert_eq!(
            io.bTag2,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_REVR as i32
                | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_XYZR as i32
        );

        let split_atoms = c_text(&mut heap, "CaaaaBbbbb");
        let split_bonds = c_text(&mut heap, "a;b;c;");
        let split_coords = c_text(&mut heap, "x;y;");
        let split = ORIG_STRUCT {
            num_atoms: 2,
            szAtoms: split_atoms,
            szBonds: split_bonds,
            szCoord: split_coords,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(
            run(&mut heap, Some(&split), 8, 0, 0).1,
            "/rA:Caaaa\nBbbbb\n\n/rB:a;b;c;\n\n/rC:x;y;\n\n"
        );
        assert_eq!(run(&mut heap, Some(&structure), 80, 1, 1).1, "");

        let missing = ORIG_STRUCT {
            num_atoms: 1,
            szAtoms: atoms,
            szBonds: SourceMutPointer::null(),
            szCoord: coords,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(run(&mut heap, Some(&missing), 80, 1, 0).1, "");
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_polymerinfo__line_4923() {
        use crate::source_types::{OAD_Polymer, OAD_PolymerUnit};

        fn buffer(heap: &mut SourceHeap, size: usize) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; size]).unwrap(),
                nAllocatedLength: size as i32,
                nPtr: size as i32,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn output_text(heap: &SourceHeap, output: &INCHI_IOSTREAM) -> String {
            let bytes = heap.slice(output.s.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|byte| *byte as u8 as char)
                .collect()
        }
        fn smt(value: &str) -> [i8; 80] {
            let mut result = [0_i8; 80];
            for (destination, source) in result.iter_mut().zip(value.bytes()) {
                *destination = source as i8;
            }
            result
        }
        fn run(
            heap: &mut SourceHeap,
            structure: Option<&ORIG_STRUCT>,
        ) -> (Result<i32, SourceHeapError>, String) {
            let mut string_buffer = buffer(heap, 32);
            let mut output = INCHI_IOSTREAM {
                s: buffer(heap, 32),
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            };
            let line_feed = heap.allocate_model_storage(vec![b'\n' as i8, 0]).unwrap();
            let mut io = INCHI_OUT_CTL::default();
            let result = OutputAUXINFO_PolymerInfo(
                heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                structure,
                &mut io,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            );
            (result, output_text(heap, &output))
        }

        let mut heap = SourceHeap::default();
        assert_eq!(run(&mut heap, None), (Ok(0), String::new()));
        assert_eq!(
            run(&mut heap, Some(&ORIG_STRUCT::default())),
            (Ok(0), String::new())
        );

        let empty_polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let empty_structure = ORIG_STRUCT {
            polymer: empty_polymer,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(
            run(&mut heap, Some(&empty_structure)),
            (Ok(0), "/Z:\n".to_owned())
        );

        let first_atoms = heap.allocate_model_storage(vec![1, 2, 3, 5]).unwrap();
        let second_atoms = heap.allocate_model_storage(vec![7]).unwrap();
        let second_bonds = heap.allocate_model_storage(vec![1, 4, 2, 9]).unwrap();
        let first = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                type_: 1,
                subtype: 2,
                conn: 3,
                na: 4,
                alist: first_atoms,
                xbr1: [777_777.777, 9.0, 9.0, 9.0],
                xbr2: [-777_777.777, 8.0, 8.0, 8.0],
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let second = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                type_: 4,
                subtype: 5,
                conn: 6,
                smt: smt("eu"),
                na: 1,
                nb: 2,
                alist: second_atoms,
                blist: second_bonds,
                xbr1: [1.25, -0.0, 2.0, 3.5],
                xbr2: [f64::NAN, 1.0, 2.0, 3.0],
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let unit_pointers = heap.allocate_model_storage(vec![first, second]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units: unit_pointers,
                n: 2,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let structure = ORIG_STRUCT {
            polymer,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(
            run(&mut heap, Some(&structure)),
            (
                Ok(0),
                "/Z:123-n-1-3,5;456-eu-7(1,4,2,9)[1.250000,-0.000000,2.000000,3.500000]\n"
                    .to_owned()
            )
        );

        let missing_units = heap
            .allocate_model_storage(vec![OAD_Polymer {
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let missing_units_structure = ORIG_STRUCT {
            polymer: missing_units,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(
            run(&mut heap, Some(&missing_units_structure)),
            (Err(SourceHeapError::NullPointer), "/Z:".to_owned())
        );

        let no_atoms_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                type_: 1,
                subtype: 1,
                conn: 1,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let no_atoms_units = heap.allocate_model_storage(vec![no_atoms_unit]).unwrap();
        let no_atoms_polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units: no_atoms_units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let no_atoms_structure = ORIG_STRUCT {
            polymer: no_atoms_polymer,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(
            run(&mut heap, Some(&no_atoms_structure)),
            (Err(SourceHeapError::PointerOutOfBounds), "/Z:".to_owned())
        );

        let overflow_atom = heap.allocate_model_storage(vec![1]).unwrap();
        let overflow_unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 1,
                nb: i32::MAX,
                alist: overflow_atom,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let overflow_units = heap.allocate_model_storage(vec![overflow_unit]).unwrap();
        let overflow_polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units: overflow_units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        let overflow_structure = ORIG_STRUCT {
            polymer: overflow_polymer,
            ..ORIG_STRUCT::default()
        };
        assert_eq!(
            run(&mut heap, Some(&overflow_structure)).0,
            Err(SourceHeapError::SourceIntegerOverflow)
        );
    }

    #[test]
    fn source_port__ichiprt1__outputinchi2__line_966() {
        use crate::source_types::{
            INCHI_OUT_NO_AUX_INFO, INCHI_OUT_PLAIN_TEXT, INCHI_OUT_PLAIN_TEXT_COMMENTS,
            INCHI_OUT_TABBED_OUTPUT, INPUT_PARMS,
        };

        fn buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 4096]).unwrap(),
                nAllocatedLength: 4096,
                nPtr: 4096,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn stream(heap: &mut SourceHeap) -> INCHI_IOSTREAM {
            INCHI_IOSTREAM {
                s: buffer(heap),
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            }
        }
        fn text(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> String {
            let bytes = heap.slice(stream.s.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|byte| *byte as u8 as char)
                .collect()
        }
        fn run(options: i32, provide_buffer: bool) -> (i32, String, String) {
            let mut heap = SourceHeap::default();
            let rows = heap
                .allocate_model_storage(vec![SourceMutPointer::null(); 4])
                .unwrap();
            let flags = heap.allocate_model_storage(vec![0_i32]).unwrap();
            let mut string_buffer = buffer(&mut heap);
            let mut output = stream(&mut heap);
            let mut log = stream(&mut heap);
            let result = OutputINChI2(
                &mut heap,
                SourceMutPointer::null(),
                provide_buffer.then_some(&mut string_buffer),
                rows,
                crate::source_types::INCHI_BAS as i32,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &INPUT_PARMS::default(),
                0,
                crate::source_types::OUT_T1 as i32,
                options,
                &[0, 0],
                &[0, 0],
                &[0, 0],
                &mut output,
                &mut log,
                23,
                flags,
                0,
                SourceMutPointer::null(),
            )
            .unwrap();
            (result, text(&heap, &output), text(&heap, &log))
        }

        let no_aux = INCHI_OUT_NO_AUX_INFO as i32;
        assert_eq!(run(no_aux, false), (0, String::new(), String::new()));
        assert_eq!(
            run(no_aux | INCHI_OUT_PLAIN_TEXT as i32, true),
            (1, "InChI=1//\n".to_owned(), String::new())
        );
        assert_eq!(
            run(
                no_aux | INCHI_OUT_PLAIN_TEXT_COMMENTS as i32 | INCHI_OUT_TABBED_OUTPUT as i32,
                true,
            ),
            (1, "\nInChI=\n1\n".to_owned(), String::new())
        );
        assert_eq!(
            run(
                no_aux
                    | INCHI_OUT_PLAIN_TEXT as i32
                    | INCHI_OUT_PLAIN_TEXT_COMMENTS as i32
                    | INCHI_OUT_TABBED_OUTPUT as i32,
                true,
            ),
            (1, "InChI=1//\n\nInChI=\n1\n".to_owned(), String::new())
        );
        assert_eq!(
            run(
                no_aux | INCHI_OUT_PLAIN_TEXT as i32 | INCHI_OUT_PLAIN_TEXT_COMMENTS as i32,
                false,
            ),
            (
                0,
                String::new(),
                concat!(
                    "Cannot allocate output buffer. No output for structure #23.\n",
                    "Cannot allocate output buffer. No output for structure #23.\n"
                )
                .to_owned()
            )
        );
    }

    #[test]
    fn source_port__ichiprt1__outputinchi1__line_1043() {
        use crate::source_types::{
            INCHI_OUT_EMBED_REC, INCHI_OUT_NO_AUX_INFO, INCHI_OUT_ONLY_AUX_INFO,
            INCHI_OUT_PLAIN_TEXT, INCHI_OUT_PLAIN_TEXT_COMMENTS, INCHI_OUT_SAVEOPT,
            INCHI_OUT_STDINCHI, INCHI_SORT, INChI, INChI_Aux, INPUT_PARMS,
        };

        fn buffer(heap: &mut SourceHeap, size: usize) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; size]).unwrap(),
                nAllocatedLength: size as i32,
                nPtr: size as i32,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn stream(heap: &mut SourceHeap) -> INCHI_IOSTREAM {
            INCHI_IOSTREAM {
                s: buffer(heap, 4096),
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            }
        }
        fn text(heap: &SourceHeap, stream: &INCHI_IOSTREAM) -> String {
            let bytes = heap.slice(stream.s.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|byte| *byte as u8 as char)
                .collect()
        }
        fn c_text(heap: &mut SourceHeap, value: &str) -> SourceMutPointer<i8> {
            heap.allocate_model_storage(
                value
                    .bytes()
                    .map(|byte| byte as i8)
                    .chain(std::iter::once(0))
                    .collect(),
            )
            .unwrap()
        }
        fn empty_rows(heap: &mut SourceHeap) -> SourceMutPointer<SourceMutPointer<INCHI_SORT>> {
            heap.allocate_model_storage(vec![SourceMutPointer::null(); 4])
                .unwrap()
        }
        fn run_empty(
            output_options: i32,
            output_type: i32,
            save_bits: u8,
            input: INPUT_PARMS,
            provide_buffer: bool,
        ) -> (Result<i32, SourceHeapError>, String, String, i32) {
            let mut heap = SourceHeap::default();
            let rows = empty_rows(&mut heap);
            let flags = heap.allocate_model_storage(vec![0_i32]).unwrap();
            let mut output = stream(&mut heap);
            let mut log = stream(&mut heap);
            let mut string_buffer = buffer(&mut heap, 4096);
            let result = OutputINChI1(
                &mut heap,
                SourceMutPointer::null(),
                provide_buffer.then_some(&mut string_buffer),
                rows,
                crate::source_types::INCHI_BAS as i32,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &input,
                0,
                output_type,
                output_options,
                &[0, 0],
                &[0, 0],
                &[0, 0],
                &mut output,
                &mut log,
                17,
                flags,
                save_bits,
                SourceMutPointer::null(),
            );
            let flag_value = heap.slice(flags.as_const()).unwrap()[0];
            (result, text(&heap, &output), text(&heap, &log), flag_value)
        }

        let no_aux = (INCHI_OUT_NO_AUX_INFO | INCHI_OUT_PLAIN_TEXT) as i32;
        assert_eq!(
            run_empty(
                no_aux,
                crate::source_types::OUT_T1 as i32,
                0,
                INPUT_PARMS::default(),
                true,
            ),
            (Ok(1), "InChI=1//\n".to_owned(), String::new(), 0)
        );
        assert_eq!(
            run_empty(
                no_aux | INCHI_OUT_STDINCHI as i32,
                crate::source_types::OUT_T1 as i32,
                0,
                INPUT_PARMS::default(),
                true,
            ),
            (Ok(1), "InChI=1S//\n".to_owned(), String::new(), 0)
        );
        assert_eq!(
            run_empty(
                no_aux | INCHI_OUT_SAVEOPT as i32,
                crate::source_types::OUT_T1 as i32,
                0x2f,
                INPUT_PARMS::default(),
                true,
            ),
            (Ok(1), "InChI=1\\PC//\n".to_owned(), String::new(), 0)
        );
        assert_eq!(
            run_empty(
                no_aux | INCHI_OUT_PLAIN_TEXT_COMMENTS as i32,
                crate::source_types::OUT_T1 as i32,
                0,
                INPUT_PARMS::default(),
                true,
            ),
            (Ok(1), "\nInChI=\n1\n".to_owned(), String::new(), 0)
        );
        assert_eq!(
            run_empty(
                INCHI_OUT_ONLY_AUX_INFO as i32,
                crate::source_types::OUT_T1 as i32,
                0,
                INPUT_PARMS::default(),
                true,
            ),
            (Ok(1), "AuxInfo=???1//\n".to_owned(), String::new(), 0)
        );
        assert_eq!(
            run_empty(
                (INCHI_OUT_ONLY_AUX_INFO | INCHI_OUT_PLAIN_TEXT) as i32,
                crate::source_types::OUT_T1 as i32,
                0,
                INPUT_PARMS::default(),
                true,
            ),
            (Ok(1), "AuxInfo=1//\n".to_owned(), String::new(), 0)
        );
        assert_eq!(
            run_empty(no_aux, 99, 0, INPUT_PARMS::default(), true,),
            (Ok(0), String::new(), String::new(), 0)
        );

        let mut heap = SourceHeap::default();
        let label = c_text(&mut heap, "ID");
        let value = c_text(&mut heap, "42");
        let input = INPUT_PARMS {
            pSdfLabel: label,
            pSdfValue: value,
            ..INPUT_PARMS::default()
        };
        let rows = empty_rows(&mut heap);
        let flags = heap.allocate_model_storage(vec![0_i32]).unwrap();
        let mut output = stream(&mut heap);
        let mut log = stream(&mut heap);
        assert_eq!(
            OutputINChI1(
                &mut heap,
                SourceMutPointer::null(),
                None,
                rows,
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &input,
                0,
                crate::source_types::OUT_T1 as i32,
                no_aux,
                &[0, 0],
                &[0, 0],
                &[0, 0],
                &mut output,
                &mut log,
                17,
                flags,
                0,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(text(&heap, &output), "");
        assert_eq!(
            text(&heap, &log),
            "Cannot allocate output buffer. No output for structure #17. ID=42\n"
        );

        let mut heap = SourceHeap::default();
        let formula = c_text(&mut heap, "CH4");
        let atoms = heap.allocate_model_storage(vec![6_u8]).unwrap();
        let hydrogens = heap.allocate_model_storage(vec![4_i8]).unwrap();
        let mobile_inchi = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                szHillFormula: formula,
                nAtom: atoms,
                nNum_H: hydrogens,
                ..INChI::default()
            }])
            .unwrap();
        let atom_numbers = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let mobile_aux = heap
            .allocate_model_storage(vec![INChI_Aux {
                nNumberOfAtoms: 1,
                nOrigAtNosInCanonOrd: atom_numbers,
                ..INChI_Aux::default()
            }])
            .unwrap();
        let mobile_sort = heap
            .allocate_model_storage(vec![INCHI_SORT {
                pINChI: [SourceMutPointer::null(), mobile_inchi],
                pINChI_Aux: [SourceMutPointer::null(), mobile_aux],
                ..INCHI_SORT::default()
            }])
            .unwrap();
        let fixed_sort = heap
            .allocate_model_storage(vec![INCHI_SORT::default()])
            .unwrap();
        let rows = heap
            .allocate_model_storage(vec![
                fixed_sort,
                mobile_sort,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
            ])
            .unwrap();
        let flags = heap.allocate_model_storage(vec![0_i32]).unwrap();
        let mut output = stream(&mut heap);
        let mut log = stream(&mut heap);
        let mut string_buffer = buffer(&mut heap, 4096);
        assert_eq!(
            OutputINChI1(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut string_buffer),
                rows,
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &INPUT_PARMS::default(),
                0,
                crate::source_types::OUT_T1 as i32,
                no_aux,
                &[1, 0],
                &[0, 0],
                &[1, 0],
                &mut output,
                &mut log,
                1,
                flags,
                0,
                SourceMutPointer::null(),
            ),
            Ok(1)
        );
        assert_eq!(text(&heap, &output), "InChI=1/CH4/h1H4\n");
        assert_eq!(text(&heap, &log), "");

        let recursive_rows = heap
            .allocate_model_storage(vec![fixed_sort, mobile_sort, fixed_sort, mobile_sort])
            .unwrap();
        let recursive_flags = heap.allocate_model_storage(vec![0_i32]).unwrap();
        let mut recursive_output = stream(&mut heap);
        let mut recursive_log = stream(&mut heap);
        let mut recursive_buffer = buffer(&mut heap, 4096);
        assert_eq!(
            OutputINChI1(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut recursive_buffer),
                recursive_rows,
                0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &INPUT_PARMS::default(),
                1,
                crate::source_types::OUT_T1 as i32,
                no_aux | INCHI_OUT_EMBED_REC as i32,
                &[1, 1],
                &[0, 0],
                &[1, 1],
                &mut recursive_output,
                &mut recursive_log,
                1,
                recursive_flags,
                0,
                SourceMutPointer::null(),
            ),
            Ok(1)
        );
        assert_eq!(
            text(&heap, &recursive_output),
            "InChI=1/CH4/h1H4/rCH4/h1H4\n"
        );
        assert_eq!(text(&heap, &recursive_log), "");
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_headerandnormalization_type__line_4339() {
        fn buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 128]).unwrap(),
                nAllocatedLength: 128,
                ..INCHI_IOS_STRING::default()
            }
        }

        let mut heap = SourceHeap::default();
        let line_feed = heap.allocate_model_storage(vec![b'\n' as i8, 0]).unwrap();
        let counts = heap.allocate_model_storage(vec![0_i32, 0]).unwrap();
        let mut string_buffer = buffer(&mut heap);
        let mut basic = crate::source_types::INCHI_BAS as i32;
        let mut control = INCHI_OUT_CTL {
            nTag: 2,
            bPlainTextTags: 1,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_HeaderAndNormalization_type(
                &mut heap,
                SourceMutPointer::null(),
                None,
                &mut string_buffer,
                0,
                &mut basic,
                counts.as_const(),
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            &heap.slice(string_buffer.pStr.as_const()).unwrap()
                [..string_buffer.nUsedLength as usize],
            &[b'1' as i8]
        );
        assert_eq!(
            control.bTag1,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_VERS as i32
        );

        let mut normalized = buffer(&mut heap);
        let normalized_counts = heap.allocate_model_storage(vec![1_i32, 0]).unwrap();
        let mut reconnected = crate::source_types::INCHI_REC as i32;
        let mut normalized_control = INCHI_OUT_CTL {
            nTag: 2,
            bPlainTextTags: 1,
            bTautomeric: 1,
            bTautomericOutputAllowed: 1,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_HeaderAndNormalization_type(
                &mut heap,
                SourceMutPointer::null(),
                None,
                &mut normalized,
                0,
                &mut reconnected,
                normalized_counts.as_const(),
                &mut normalized_control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            &heap.slice(normalized.pStr.as_const()).unwrap()[..normalized.nUsedLength as usize],
            &[b'/' as i8, b'1' as i8]
        );
        assert_eq!(
            normalized_control.bTag1,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_NORM as i32
        );
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_originalnumbersandequivalenceclasses__line_4396() {
        let mut heap = SourceHeap::default();
        let line_feed = heap.allocate_model_storage(vec![b'\n' as i8, 0]).unwrap();
        let counts = heap.allocate_model_storage(vec![0_i32, 0]).unwrap();
        let mut string_buffer = INCHI_IOS_STRING {
            pStr: heap.allocate_model_storage(vec![0_i8; 128]).unwrap(),
            nAllocatedLength: 128,
            ..INCHI_IOS_STRING::default()
        };
        let mut output = INCHI_IOSTREAM {
            s: INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 128]).unwrap(),
                nAllocatedLength: 128,
                ..INCHI_IOS_STRING::default()
            },
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let mut control = INCHI_OUT_CTL {
            nTag: 2,
            bPlainTextTags: 1,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_OriginalNumbersAndEquivalenceClasses(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                counts.as_const(),
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(output.s.nUsedLength, 1);
        assert_eq!(
            &heap.slice(output.s.pStr.as_const()).unwrap()[..2],
            &[b'/' as i8, 0]
        );

        control.bAtomEqu[0] = 1;
        assert_eq!(
            OutputAUXINFO_OriginalNumbersAndEquivalenceClasses(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                counts.as_const(),
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            control.bTag1,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_AEQU as i32
        );
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_tautomericgroupsequivalence__line_4461() {
        fn buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 128]).unwrap(),
                nAllocatedLength: 128,
                ..INCHI_IOS_STRING::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut string_buffer = buffer(&mut heap);
        let mut output = INCHI_IOSTREAM {
            s: buffer(&mut heap),
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let mut control = INCHI_OUT_CTL {
            bTautomericOutputAllowed: 1,
            bTautomeric: 1,
            bPlainTextTags: 1,
            nTag: 2,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_TautomericGroupsEquivalence(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            &heap.slice(output.s.pStr.as_const()).unwrap()[..2],
            &[b'/' as i8, 0]
        );

        let numbers = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let aux = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Aux {
                nNumberOfAtoms: 2,
                nNumberOfTGroups: 2,
                nConstitEquTGroupNumbers: numbers,
                ..crate::source_types::INChI_Aux::default()
            }])
            .unwrap();
        let inchi = heap
            .allocate_model_storage(vec![crate::source_types::INChI {
                nNumberOfAtoms: 2,
                ..crate::source_types::INChI::default()
            }])
            .unwrap();
        control.pINChISort = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        control.bTautEqu[0] = 1;
        control.bOutType = crate::source_types::OUT_T1 as i32;
        control.num_components = 1;
        output.s.nUsedLength = 0;
        heap.slice_mut(output.s.pStr).unwrap()[0] = 0;
        assert_eq!(
            OutputAUXINFO_TautomericGroupsEquivalence(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        let output_bytes = heap.slice(output.s.pStr.as_const()).unwrap();
        let output_end = output_bytes.iter().position(|byte| *byte == 0).unwrap();
        assert_eq!(
            &output_bytes[..output_end],
            &b"/gE:(1,2)"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            control.bTag1,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_GEQU as i32
        );
        assert_eq!(control.tot_len, 5);

        control.bOverflow = 1;
        assert_eq!(
            OutputAUXINFO_TautomericGroupsEquivalence(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                SourceMutPointer::null(),
            ),
            Ok(1)
        );
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_stereo__line_4500() {
        fn buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 256]).unwrap(),
                nAllocatedLength: 256,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn output_text(heap: &SourceHeap, output: &INCHI_IOSTREAM) -> String {
            let bytes = heap.slice(output.s.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|byte| *byte as u8 as char)
                .collect()
        }

        let mut heap = SourceHeap::default();
        let mut string_buffer = buffer(&mut heap);
        let mut output = INCHI_IOSTREAM {
            s: buffer(&mut heap),
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let line_feed = heap.allocate_model_storage(vec![b'\n' as i8, 0]).unwrap();
        let mut control = INCHI_OUT_CTL {
            nTag: 2,
            bPlainTextTags: 1,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_Stereo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "//");

        let number = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let parity = heap.allocate_model_storage(vec![1_i8]).unwrap();
        let number_inv = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let parity_inv = heap.allocate_model_storage(vec![2_i8]).unwrap();
        let stereo = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: number,
                t_parity: parity,
                nNumberInv: number_inv,
                t_parityInv: parity_inv,
                nCompInv2Abs: 1,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let inchi = heap
            .allocate_model_storage(vec![crate::source_types::INChI {
                nNumberOfAtoms: 1,
                Stereo: stereo,
                ..crate::source_types::INChI::default()
            }])
            .unwrap();
        let aux = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Aux {
                nNumberOfAtoms: 1,
                nOrigAtNosInCanonOrd: number,
                nOrigAtNosInCanonOrdInv: number_inv,
                ..crate::source_types::INChI_Aux::default()
            }])
            .unwrap();
        control.pINChISort = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        control.bInvStereo[0] = 1;
        control.bInvStereoOrigNumb[0] = 1;
        control.bOutType = crate::source_types::OUT_T1 as i32;
        control.num_components = 1;
        output.s.nUsedLength = 0;
        heap.slice_mut(output.s.pStr).unwrap()[0] = 0;
        assert_eq!(
            OutputAUXINFO_Stereo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "/it:1+\n/iN:1\n");
        assert_eq!(
            control.bTag2,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_STER as i32
                | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_SP3N as i32
        );
        assert_eq!(control.tot_len, 1);

        control.bInvStereoOrigNumb[0] = 0;
        output.s.nUsedLength = 0;
        heap.slice_mut(output.s.pStr).unwrap()[0] = 0;
        assert_eq!(
            OutputAUXINFO_Stereo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "/it:1+\n/");
    }

    #[test]
    fn source_port__ichiprt1__outputauxinfo_isotopicinfo__line_4560() {
        fn buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 512]).unwrap(),
                nAllocatedLength: 512,
                ..INCHI_IOS_STRING::default()
            }
        }
        fn output_text(heap: &SourceHeap, output: &INCHI_IOSTREAM) -> String {
            let bytes = heap.slice(output.s.pStr.as_const()).unwrap();
            let end = bytes.iter().position(|byte| *byte == 0).unwrap();
            bytes[..end]
                .iter()
                .map(|byte| *byte as u8 as char)
                .collect()
        }
        fn clear(heap: &mut SourceHeap, output: &mut INCHI_IOSTREAM) {
            output.s.nUsedLength = 0;
            heap.slice_mut(output.s.pStr).unwrap()[0] = 0;
        }

        let mut heap = SourceHeap::default();
        let mut string_buffer = buffer(&mut heap);
        let mut output = INCHI_IOSTREAM {
            s: buffer(&mut heap),
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        let line_feed = heap.allocate_model_storage(vec![b'\n' as i8, 0]).unwrap();
        let mut basic = crate::source_types::INCHI_BAS as i32;

        let mut inactive = INCHI_OUT_CTL {
            nTag: 2,
            bPlainTextTags: 1,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_IsotopicInfo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut basic,
                &mut inactive,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "");

        let mut marker_only = INCHI_OUT_CTL {
            nTag: 2,
            bPlainTextTags: 1,
            bIsotopic: 1,
            bIgn_UU_Sp2_Iso: [1, 0],
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            OutputAUXINFO_IsotopicInfo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut basic,
                &mut marker_only,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "/I:\n///");

        let normal_numbers = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let isotopic_numbers = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let inverted_numbers = heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let isotopic_inverted_numbers = heap.allocate_model_storage(vec![2_u16, 1]).unwrap();
        let normal_equivalence = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let isotopic_equivalence = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let normal_tgroup_equivalence = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let isotopic_tgroup_equivalence = heap.allocate_model_storage(vec![1_u16, 1]).unwrap();
        let stereo_number = heap.allocate_model_storage(vec![2_u16]).unwrap();
        let stereo_parity = heap.allocate_model_storage(vec![1_i8]).unwrap();
        let stereo_number_inverted = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let stereo_parity_inverted = heap.allocate_model_storage(vec![2_i8]).unwrap();
        let stereo = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                nNumberOfStereoCenters: 1,
                nNumber: stereo_number,
                t_parity: stereo_parity,
                nNumberInv: stereo_number_inverted,
                t_parityInv: stereo_parity_inverted,
                nCompInv2Abs: 1,
                ..crate::source_types::INChI_Stereo::default()
            }])
            .unwrap();
        let inchi = heap
            .allocate_model_storage(vec![crate::source_types::INChI {
                nNumberOfAtoms: 2,
                nNumberOfIsotopicAtoms: 1,
                StereoIsotopic: stereo,
                ..crate::source_types::INChI::default()
            }])
            .unwrap();
        let aux = heap
            .allocate_model_storage(vec![crate::source_types::INChI_Aux {
                nNumberOfAtoms: 2,
                nNumberOfTGroups: 2,
                bIsIsotopic: 1,
                nOrigAtNosInCanonOrd: normal_numbers,
                nIsotopicOrigAtNosInCanonOrd: isotopic_numbers,
                nOrigAtNosInCanonOrdInv: inverted_numbers,
                nIsotopicOrigAtNosInCanonOrdInv: isotopic_inverted_numbers,
                nConstitEquNumbers: normal_equivalence,
                nConstitEquIsotopicNumbers: isotopic_equivalence,
                nConstitEquTGroupNumbers: normal_tgroup_equivalence,
                nConstitEquIsotopicTGroupNumbers: isotopic_tgroup_equivalence,
                ..crate::source_types::INChI_Aux::default()
            }])
            .unwrap();
        let sort = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [inchi, SourceMutPointer::null()],
                pINChI_Aux: [aux, SourceMutPointer::null()],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        let mut control = INCHI_OUT_CTL {
            nTag: 2,
            bPlainTextTags: 1,
            bIsotopic: 1,
            bIsotopicOrigNumb: [1, 0],
            bIsotopicAtomEqu: [1, 0],
            bIsotopicTautEqu: [1, 0],
            bInvIsotopicStereo: [1, 0],
            bInvIsotopicStereoOrigNumb: [1, 0],
            bTautomericOutputAllowed: 1,
            bTautomeric: 1,
            bOutType: crate::source_types::OUT_T1 as i32,
            num_components: 1,
            pINChISort: sort,
            ..INCHI_OUT_CTL::default()
        };
        clear(&mut heap, &mut output);
        assert_eq!(
            OutputAUXINFO_IsotopicInfo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut basic,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(
            output_text(&heap, &output),
            "/I:1,2\n/E:(1,2)\n/gE:(1,2)\n/it:1+/iN:2,1\n"
        );
        assert_eq!(
            control.bTag3,
            crate::source_types::local_ichiprt1::tagAuxLblBit_AL_ISOT as i32
                | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_STER as i32
                | crate::source_types::local_ichiprt1::tagAuxLblBit_AL_SP3N as i32
        );
        assert_eq!(control.tot_len, 3);

        let flags = heap
            .allocate_model_storage(vec![
                crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_BAS as i32,
            ])
            .unwrap();
        control.bSecondNonTautPass = 1;
        control.pSortPrintINChIFlags = flags;
        clear(&mut heap, &mut output);
        assert_eq!(
            OutputAUXINFO_IsotopicInfo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut basic,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        assert_eq!(output_text(&heap, &output), "");

        control.bSecondNonTautPass = 0;
        control.bOverflow = 1;
        clear(&mut heap, &mut output);
        assert_eq!(
            OutputAUXINFO_IsotopicInfo(
                &mut heap,
                SourceMutPointer::null(),
                Some(&mut output),
                &mut string_buffer,
                &mut basic,
                &mut control,
                line_feed.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(1)
        );
        assert_eq!(output_text(&heap, &output), "");
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_polymerlayer__line_3880() {
        use crate::source_types::{
            INCHI_SORT, INChI, INChI_Aux, OAD_Polymer, OAD_PolymerUnit, ORIG_ATOM_DATA,
            ORIG_STRUCT, OUT_T1, TAUT_YES, inp_ATOM,
        };

        fn output_buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![0_i8; 128]).unwrap(),
                nAllocatedLength: 128,
                ..INCHI_IOS_STRING::default()
            }
        }

        let mut heap = SourceHeap::default();
        let mut buffer = output_buffer(&mut heap);
        assert_eq!(
            OutputINCHI_PolymerLayer(
                &mut heap,
                SourceMutPointer::null(),
                None,
                &mut buffer,
                &mut 0,
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                &mut INCHI_OUT_CTL::default(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer::default()])
            .unwrap();
        let empty_structure = heap
            .allocate_model_storage(vec![ORIG_STRUCT::default()])
            .unwrap();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let original = heap
            .allocate_model_storage(vec![ORIG_ATOM_DATA {
                at: atoms,
                num_inp_atoms: 1,
                ..ORIG_ATOM_DATA::default()
            }])
            .unwrap();
        assert_eq!(
            OutputINCHI_PolymerLayer(
                &mut heap,
                SourceMutPointer::null(),
                None,
                &mut buffer,
                &mut 0,
                original.as_const(),
                empty_structure.as_const(),
                &mut INCHI_OUT_CTL::default(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );

        let inchi2inchi = heap
            .allocate_model_storage(vec![ORIG_STRUCT {
                polymer,
                ..ORIG_STRUCT::default()
            }])
            .unwrap();
        assert_eq!(
            OutputINCHI_PolymerLayer(
                &mut heap,
                SourceMutPointer::null(),
                None,
                &mut buffer,
                &mut 0,
                original.as_const(),
                inchi2inchi.as_const(),
                &mut INCHI_OUT_CTL::default(),
                SourceConstPointer::null(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(crate::source_types::local_ichiprt1::NOT_YET_I2I_FOR_POLYMERS as i32)
        );

        let numbers = heap.allocate_model_storage(vec![1_u16]).unwrap();
        let aux = heap
            .allocate_model_storage(vec![INChI_Aux {
                nNumberOfAtoms: 1,
                nOrigAtNosInCanonOrd: numbers,
                ..INChI_Aux::default()
            }])
            .unwrap();
        let inchi = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 1,
                lenTautomer: 1,
                ..INChI::default()
            }])
            .unwrap();
        let mut sort = INCHI_SORT::default();
        sort.pINChI[TAUT_YES as usize] = inchi;
        sort.pINChI_Aux[TAUT_YES as usize] = aux;
        let sort = heap.allocate_model_storage(vec![sort]).unwrap();
        let unit_atoms = heap.allocate_model_storage(vec![0_i32]).unwrap();
        let unit = heap
            .allocate_model_storage(vec![OAD_PolymerUnit {
                na: 1,
                alist: unit_atoms,
                ..OAD_PolymerUnit::default()
            }])
            .unwrap();
        let units = heap.allocate_model_storage(vec![unit]).unwrap();
        let polymer = heap
            .allocate_model_storage(vec![OAD_Polymer {
                units,
                n: 1,
                ..OAD_Polymer::default()
            }])
            .unwrap();
        heap.slice_mut(empty_structure).unwrap()[0].polymer = polymer;
        heap.slice_mut(original).unwrap()[0].polymer = polymer;
        heap.slice_mut(atoms).unwrap()[0].orig_at_number = 0;
        heap.slice_mut(atoms).unwrap()[0].el_number = 6;
        let atom_text = heap.allocate_model_storage(vec![b'C' as i8, 0]).unwrap();
        heap.slice_mut(empty_structure).unwrap()[0].szAtoms = atom_text;
        let newline = heap.allocate_model_storage(vec![b'\n' as i8, 0]).unwrap();
        let mut control = INCHI_OUT_CTL {
            pINChISort: sort,
            bOutType: OUT_T1 as i32,
            num_components: 0,
            ..INCHI_OUT_CTL::default()
        };
        let mut basic = 0;
        assert_eq!(
            OutputINCHI_PolymerLayer(
                &mut heap,
                SourceMutPointer::null(),
                None,
                &mut buffer,
                &mut basic,
                original.as_const(),
                empty_structure.as_const(),
                &mut control,
                newline.as_const(),
                SourceConstPointer::null(),
                SourceMutPointer::null(),
            ),
            Ok(0)
        );
        let text = heap.slice(buffer.pStr.as_const()).unwrap();
        assert_eq!(
            &text[..buffer.nUsedLength as usize],
            &[
                b'/' as i8, b'z' as i8, b'0' as i8, b'0' as i8, b'0' as i8, b'-' as i8, b'1' as i8
            ]
        );
        assert_eq!(heap.slice(polymer.as_const()).unwrap()[0].n, 1);
    }

    #[test]
    fn source_port__ichiprt1__internallygetcanonumsandcomponentnums__line_5362() {
        use crate::source_types::{INCHI_SORT, INChI, INChI_Aux, OUT_NT, OUT_T1, TAUT_YES};

        fn output_buffer(heap: &mut SourceHeap) -> INCHI_IOS_STRING {
            INCHI_IOS_STRING {
                pStr: heap.allocate_model_storage(vec![81_i8; 128]).unwrap(),
                nAllocatedLength: 128,
                nUsedLength: 7,
                nPtr: 17,
            }
        }

        fn sort_with_numbers(
            heap: &mut SourceHeap,
            numbers: &[u16],
        ) -> SourceMutPointer<INCHI_SORT> {
            let number_pointer = heap.allocate_model_storage(numbers.to_vec()).unwrap();
            let aux = heap
                .allocate_model_storage(vec![INChI_Aux {
                    nNumberOfAtoms: numbers.len() as i32,
                    nOrigAtNosInCanonOrd: number_pointer,
                    ..INChI_Aux::default()
                }])
                .unwrap();
            let inchi = heap
                .allocate_model_storage(vec![INChI {
                    nNumberOfAtoms: numbers.len() as i32,
                    lenTautomer: 1,
                    ..INChI::default()
                }])
                .unwrap();
            let mut sort = INCHI_SORT::default();
            sort.pINChI[TAUT_YES as usize] = inchi;
            sort.pINChI_Aux[TAUT_YES as usize] = aux;
            heap.allocate_model_storage(vec![sort]).unwrap()
        }

        let mut heap = SourceHeap::default();
        let sorted = sort_with_numbers(&mut heap, &[2, 1]);
        let second = sort_with_numbers(&mut heap, &[3]);
        let first_sort = heap.slice(sorted.as_const()).unwrap()[0].clone();
        let second_sort = heap.slice(second.as_const()).unwrap()[0].clone();
        let sorts = heap
            .allocate_model_storage(vec![first_sort, second_sort])
            .unwrap();
        let canonical = heap.allocate_model_storage(vec![77_i32; 4]).unwrap();
        let components = heap.allocate_model_storage(vec![77_i32; 4]).unwrap();
        let mut buffer = output_buffer(&mut heap);
        let mut control = INCHI_OUT_CTL {
            pINChISort: sorts,
            bOutType: OUT_T1 as i32,
            num_components: 2,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            InternallyGetCanoNumsAndComponentNums(
                &mut heap,
                SourceMutPointer::null(),
                &mut buffer,
                &mut control,
                3,
                canonical,
                components,
            ),
            Ok(0)
        );
        assert_eq!(control.tot_len, 5);
        assert_eq!(heap.slice(canonical.as_const()).unwrap(), &[77, 1, 0, 2]);
        assert_eq!(heap.slice(components.as_const()).unwrap(), &[1, 1, 2, 77]);
        assert_eq!(buffer.nUsedLength, 0);
        assert_eq!(buffer.nPtr, 0);
        assert_eq!(heap.slice(buffer.pStr.as_const()).unwrap()[0], 0);

        let mut invalid_buffer = output_buffer(&mut heap);
        let invalid_canonical = heap.allocate_model_storage(vec![55_i32; 2]).unwrap();
        let invalid_components = heap.allocate_model_storage(vec![66_i32; 1]).unwrap();
        let equal = sort_with_numbers(&mut heap, &[1, 2]);
        let equal_sort = heap.slice(equal.as_const()).unwrap()[0].clone();
        let non_numbers = heap.allocate_model_storage(vec![1_u16, 2]).unwrap();
        let non_aux = heap
            .allocate_model_storage(vec![INChI_Aux {
                nNumberOfAtoms: 2,
                nOrigAtNosInCanonOrd: non_numbers,
                ..INChI_Aux::default()
            }])
            .unwrap();
        let non_inchi = heap
            .allocate_model_storage(vec![INChI {
                nNumberOfAtoms: 2,
                lenTautomer: 0,
                ..INChI::default()
            }])
            .unwrap();
        let mut equal_sort = equal_sort;
        equal_sort.pINChI[0] = non_inchi;
        equal_sort.pINChI_Aux[0] = non_aux;
        let equal_sorts = heap.allocate_model_storage(vec![equal_sort]).unwrap();
        let mut invalid_control = INCHI_OUT_CTL {
            pINChISort: equal_sorts,
            bOutType: OUT_NT as i32,
            TAUT_MODE: 0,
            bSecondNonTautPass: 1,
            bOmitRepetitions: 1,
            num_components: 1,
            ..INCHI_OUT_CTL::default()
        };
        assert_eq!(
            InternallyGetCanoNumsAndComponentNums(
                &mut heap,
                SourceMutPointer::null(),
                &mut invalid_buffer,
                &mut invalid_control,
                1,
                invalid_canonical,
                invalid_components,
            ),
            Ok(2)
        );
        assert_eq!(invalid_buffer.nUsedLength, 0);
        assert_eq!(heap.slice(invalid_canonical.as_const()).unwrap(), &[55, -1]);
        assert_eq!(heap.slice(invalid_components.as_const()).unwrap(), &[-1]);

        let mut null_buffer = INCHI_IOS_STRING::default();
        let mut null_control = INCHI_OUT_CTL::default();
        assert_eq!(
            InternallyGetCanoNumsAndComponentNums(
                &mut heap,
                SourceMutPointer::null(),
                &mut null_buffer,
                &mut null_control,
                0,
                canonical,
                components,
            ),
            Ok(1)
        );
        assert_eq!(
            InternallyGetCanoNumsAndComponentNums(
                &mut heap,
                SourceMutPointer::null(),
                &mut buffer,
                &mut control,
                0,
                SourceMutPointer::null(),
                components,
            ),
            Ok(1)
        );
    }

    #[test]
    fn source_port__ichiprt1__equstring__line_564() {
        const S: i32 = 0x001;
        const I: i32 = 0x002;
        const N: i32 = 0x004;
        const E: i32 = 0x008;
        const ISO: i32 = 0x010;
        const NONTAUT: i32 = 0x020;
        const EQ_NONTAUT: i32 = 0x040;
        const EQ_ISO: i32 = 0x080;
        const EQ_INV: i32 = 0x100;
        const MASK: i32 = 0x1ff;

        let mut cases = vec![
            (S, ISO, 0, "m"),
            (S, NONTAUT, 0, "m"),
            (S, ISO | NONTAUT, 0, "m"),
            (S, ISO | NONTAUT, EQ_ISO, "M"),
            (S, ISO | NONTAUT, EQ_NONTAUT, "n"),
            (I, ISO, 0, "m"),
            (I, NONTAUT, 0, "m"),
            (I, ISO | NONTAUT, 0, "m"),
            (I, ISO | NONTAUT, EQ_ISO, "M"),
            (I, ISO | NONTAUT, EQ_NONTAUT, "n"),
            (N | I, 0, 0, "m"),
            (N | I, ISO, 0, "m"),
            (N | I, ISO, EQ_INV, "im"),
            (N | I, ISO, EQ_ISO, "M"),
            (N | I, NONTAUT, 0, "m"),
            (N | I, NONTAUT, EQ_NONTAUT, "n"),
            (N | I, NONTAUT, EQ_INV, "im"),
            (N | I, ISO | NONTAUT, 0, "m"),
            (N | I, ISO | NONTAUT, EQ_ISO, "M"),
            (N | I, ISO | NONTAUT, EQ_ISO | EQ_INV, "iM"),
            (N | I, ISO | NONTAUT, EQ_NONTAUT, "n"),
            (N | I, ISO | NONTAUT, EQ_NONTAUT | EQ_ISO, "N"),
            (N | I, ISO | NONTAUT, EQ_INV, "im"),
            (N | I, ISO | NONTAUT, EQ_NONTAUT | EQ_INV, "in"),
        ];
        for layer_type in [0, ISO, NONTAUT, ISO | NONTAUT] {
            for (equal_without_inversion, expected) in [
                (0, "im"),
                (EQ_ISO, "iM"),
                (EQ_NONTAUT, "in"),
                (EQ_NONTAUT | EQ_ISO, "iN"),
            ] {
                cases.push((I, layer_type, equal_without_inversion | EQ_INV, expected));
            }
        }
        for from in [N, E] {
            cases.extend([
                (from, ISO, 0, "m"),
                (from, NONTAUT, 0, "m"),
                (from, ISO | NONTAUT, 0, "m"),
                (from, ISO | NONTAUT, EQ_ISO, "M"),
                (from, ISO | NONTAUT, EQ_NONTAUT, "n"),
            ]);
        }

        for value in 0_i32..=0x3ff {
            let active = value & MASK;
            let from = active & 0x00f;
            let layer_type = active & (ISO | NONTAUT);
            let equal_to = active & (EQ_NONTAUT | EQ_ISO | EQ_INV);
            let expected = cases
                .iter()
                .find_map(|case| {
                    (case.0 == from && case.1 == layer_type && case.2 == equal_to).then_some(case.3)
                })
                .unwrap_or("??");
            assert_eq!(EquString(value), expected, "value={value:#05x}");
            assert_eq!(
                EquString(value | i32::MIN),
                expected,
                "signed value={:#010x}",
                value | i32::MIN
            );
        }
        assert_eq!(EquString(crate::source_types::iiEmpty as i32), "??");
    }

    fn tag_text(heap: &mut SourceHeap, text: &str) -> SourceConstPointer<i8> {
        heap.allocate_model_storage(
            text.bytes()
                .map(|byte| byte as i8)
                .chain(std::iter::once(0))
                .collect(),
        )
        .unwrap()
        .as_const()
    }

    fn test_tags(heap: &mut SourceHeap) -> Vec<INCHI_TAG> {
        (0..19)
            .map(|index| INCHI_TAG {
                szPlainLabel: tag_text(heap, &format!("P{index}")),
                szPlainComment: tag_text(heap, &format!("C{index}")),
                szXmlLabel: tag_text(heap, &format!("X{index}")),
                bAlwaysOutput: 100 + index,
            })
            .collect()
    }

    fn tag_result(heap: &SourceHeap, output: SourceMutPointer<i8>) -> Vec<i8> {
        let bytes = heap.slice(output.as_const()).unwrap();
        let nul = bytes.iter().position(|byte| *byte == 0).unwrap();
        bytes[..=nul].to_vec()
    }

    #[test]
    fn source_port__ichiprt1__szgettag__line_2087() {
        let mut heap = SourceHeap::default();
        let tags = test_tags(&mut heap);
        let output = heap.allocate_model_storage(vec![77_i8; 128]).unwrap();

        for invalid in [i32::MIN, -1, 0, 4, i32::MAX] {
            let mut always = 41;
            assert_eq!(
                szGetTag(&mut heap, &tags, invalid, -1, output, &mut always, 1),
                Ok(output)
            );
            assert_eq!(
                tag_result(&heap, output),
                b"???\0".iter().map(|byte| *byte as i8).collect::<Vec<_>>()
            );
            assert_eq!(always, 41);
        }

        for n_tag in [1, 2, 3] {
            let mut always = 42;
            assert_eq!(
                szGetTag(&mut heap, &tags, n_tag, 0, output, &mut always, 1),
                Ok(output)
            );
            assert_eq!(
                tag_result(&heap, output),
                b"???\0".iter().map(|byte| *byte as i8).collect::<Vec<_>>()
            );
            assert_eq!(always, 42);
        }

        let mut always = -7;
        assert_eq!(
            szGetTag(
                &mut heap,
                &tags,
                1,
                (1 << 2) | (1 << 17),
                output,
                &mut always,
                1
            ),
            Ok(output)
        );
        assert_eq!(
            tag_result(&heap, output),
            b"X17\0".iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(always, 117);

        always = -8;
        assert_eq!(
            szGetTag(&mut heap, &tags, 2, 1 << 18, output, &mut always, 1),
            Ok(output)
        );
        assert_eq!(
            tag_result(&heap, output),
            b"P18\0".iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(always, -8);

        always = -9;
        assert_eq!(
            szGetTag(
                &mut heap,
                &tags,
                3,
                (1 << 0) | (1 << 2) | (1 << 18),
                output,
                &mut always,
                1
            ),
            Ok(output)
        );
        assert_eq!(
            tag_result(&heap, output),
            b"P18{C0:C2:C18}\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>()
        );
        assert_eq!(always, 118);

        always = 55;
        assert_eq!(
            szGetTag(&mut heap, &tags, 1, 1 << 18, output, &mut always, 0),
            Ok(output)
        );
        assert_eq!(
            tag_result(&heap, output),
            b"???\0".iter().map(|byte| *byte as i8).collect::<Vec<_>>()
        );
        assert_eq!(always, 55);

        heap.slice_mut(output).unwrap().fill(77);
        assert_eq!(
            szGetTag(&mut heap, &tags, 3, 1 << 17, output, &mut always, 0),
            Ok(output)
        );
        assert_eq!(
            tag_result(&heap, output),
            b"P17{C17}\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>()
        );
        assert_eq!(always, 117);

        let tail = heap.slice(output.as_const()).unwrap();
        let nul = tail.iter().position(|byte| *byte == 0).unwrap();
        assert_eq!(tail[nul + 1], 77);
    }

    fn line_end_case(
        text: &[u8],
        allocated: usize,
        increment: i32,
        overflow: i32,
        ind: i32,
        plain_text_tags: i32,
    ) -> (i32, Vec<i8>, INCHI_IOS_STRING) {
        let mut heap = SourceHeap::default();
        let tag = tag_text(&mut heap, "/c");
        let mut storage = vec![88_i8; allocated];
        for (target, source) in storage.iter_mut().zip(text.iter().copied()) {
            *target = source as i8;
        }
        storage[text.len()] = 0;
        let pointer = heap.allocate_model_storage(storage).unwrap();
        let mut buffer = INCHI_IOS_STRING {
            pStr: pointer,
            nAllocatedLength: i32::try_from(allocated).unwrap(),
            nUsedLength: i32::try_from(text.len()).unwrap(),
            nPtr: increment,
        };
        let mut overflow_value = overflow;
        let result = str_LineEnd(
            &mut heap,
            tag,
            &mut overflow_value,
            &mut buffer,
            ind,
            plain_text_tags,
        )
        .unwrap();
        let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
        let nul = bytes.iter().position(|byte| *byte == 0).unwrap();
        (result, bytes[..=nul].to_vec(), buffer)
    }

    #[test]
    fn source_port__ichiprt1__str_lineend__line_2183() {
        let (result, bytes, buffer) = line_end_case(b"abc", 16, 8, 1, -1, 1);
        assert_eq!(result, 1);
        assert_eq!(bytes, vec![b'a' as i8, b'b' as i8, b'c' as i8, 0]);
        assert_eq!(buffer.nUsedLength, 3);

        for (text, ind, tags, expected) in [
            (&b""[..], -1, 1, &b""[..]),
            (&b""[..], -2, 1, &b"/c"[..]),
            (&b"abc"[..], -1, 0, &b"abc"[..]),
            (&b"abc"[..], 0, 1, &b"abc"[..]),
            (&b"abc"[..], -1, 1, &b"/cabc"[..]),
        ] {
            let (result, bytes, buffer) = line_end_case(text, 16, 8, 0, ind, tags);
            assert_eq!(result, 0);
            let mut expected_bytes: Vec<i8> = expected.iter().map(|byte| *byte as i8).collect();
            expected_bytes.push(0);
            assert_eq!(bytes, expected_bytes);
            assert_eq!(buffer.nUsedLength, i32::try_from(expected.len()).unwrap());
        }

        let (result, bytes, buffer) = line_end_case(b"abcdef", 7, 16, 0, -1, 1);
        assert_eq!(result, 0);
        assert_eq!(
            bytes,
            b"/cabcdef\0"
                .iter()
                .map(|byte| *byte as i8)
                .collect::<Vec<_>>()
        );
        assert_eq!(buffer.nUsedLength, 8);
        assert_eq!(buffer.nAllocatedLength, 23);
    }

    #[test]
    fn source_port__ichiprt1__mergezzinstrhillformulacomponent__line_5495() {
        let cases = [
            ("", ""),
            ("C2H6", "C2H6"),
            ("Zz", "Zz"),
            ("Zz2", "Zz2"),
            ("ZzZz", "Zz2"),
            ("C2ZzH3Zz", "C2Zz2"),
            ("C2Zz2H3Zz4tail", "C2Zz6"),
            ("Zz0Zz0", "Zz0"),
            ("Zz2147483647Zz1", "Zz-2147483648"),
            ("Zz4294967295Zz2", "Zz1"),
            ("Zz2Zz3Zz4", "Zz5"),
        ];
        for (input, expected) in cases {
            let mut heap = SourceHeap::default();
            let mut bytes: Vec<i8> = input.bytes().map(|byte| byte as i8).collect();
            bytes.resize(64, 91);
            bytes[input.len()] = 0;
            let original = bytes.clone();
            let string = heap.allocate_model_storage(bytes).unwrap();
            assert_eq!(MergeZzInStrHillFormulaComponent(&mut heap, string), Ok(()));
            let bytes = heap.slice(string.as_const()).unwrap();
            let nul = bytes.iter().position(|byte| *byte == 0).unwrap();
            assert_eq!(
                &bytes[..nul],
                &expected.bytes().map(|byte| byte as i8).collect::<Vec<_>>(),
                "input={input}"
            );
            if expected != input {
                assert_eq!(bytes[nul + 1], original[nul + 1], "input={input}");
            }
        }
    }

    fn merge_hill_buffer(heap: &mut SourceHeap, text: &str) -> INCHI_IOS_STRING {
        let allocated = 64;
        let mut bytes = vec![73_i8; allocated];
        for (target, source) in bytes.iter_mut().zip(text.bytes()) {
            *target = source as i8;
        }
        bytes[text.len()] = 0;
        INCHI_IOS_STRING {
            pStr: heap.allocate_model_storage(bytes).unwrap(),
            nAllocatedLength: allocated as i32,
            nUsedLength: text.len() as i32,
            nPtr: 16,
        }
    }

    #[test]
    fn source_port__ichiprt1__mergezzinhillformula__line_5435() {
        let mut heap = SourceHeap::default();
        let mut null = INCHI_IOS_STRING::default();
        assert_eq!(MergeZzInHillFormula(&mut heap, &mut null), Ok(0));

        let mut empty = merge_hill_buffer(&mut heap, "");
        let empty_before = empty.clone();
        assert_eq!(MergeZzInHillFormula(&mut heap, &mut empty), Ok(0));
        assert_eq!(empty, empty_before);

        for (input, expected) in [
            ("C2H6", "C2H6"),
            ("Zz2Zz3", "Zz5"),
            ("Zz2Zz3.C2H6.ZzZz", "Zz5.C2H6.Zz2"),
            (".ZzZz.", ".Zz2."),
        ] {
            let mut buffer = merge_hill_buffer(&mut heap, input);
            assert_eq!(MergeZzInHillFormula(&mut heap, &mut buffer), Ok(0));
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            let visible = bytes.iter().position(|byte| *byte == 0).unwrap();
            assert_eq!(
                &bytes[..visible],
                &expected.bytes().map(|byte| byte as i8).collect::<Vec<_>>(),
                "input={input}"
            );
            assert_eq!(buffer.nUsedLength, expected.len() as i32 + 1);
            assert_eq!(buffer.nPtr, 0);
        }

        for successful_allocations in [0, 1] {
            let mut heap = SourceHeap::default();
            let mut buffer = merge_hill_buffer(&mut heap, "ZzZz");
            let before = buffer.clone();
            heap.fail_after_allocations(successful_allocations);
            assert_eq!(MergeZzInHillFormula(&mut heap, &mut buffer), Ok(-1));
            assert_eq!(buffer, before);
            assert_eq!(
                &heap.slice(buffer.pStr.as_const()).unwrap()[..5],
                &[b'Z' as i8, b'z' as i8, b'Z' as i8, b'z' as i8, 0]
            );
        }
    }

    fn version_case(output_options: i32, is_beta: i32, newline: bool) -> (Vec<i8>, Vec<i8>) {
        let mut heap = SourceHeap::default();
        let output_storage = heap.allocate(vec![0_i8; 64]).unwrap();
        let mut output = INCHI_IOSTREAM {
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: output_storage,
                nAllocatedLength: 64,
                nUsedLength: 0,
                nPtr: 16,
            },
            ..INCHI_IOSTREAM::default()
        };
        let string_storage = heap.allocate(vec![99_i8; 8]).unwrap();
        let mut string_buffer = INCHI_IOS_STRING {
            pStr: string_storage,
            nAllocatedLength: 8,
            nUsedLength: 4,
            nPtr: 4,
        };
        let line_feed = heap
            .allocate(if newline {
                vec![b'\n' as i8, 0]
            } else {
                vec![0]
            })
            .unwrap();
        let tab = heap.allocate(vec![b'\n' as i8, 0]).unwrap();
        let stdout = heap
            .allocate(vec![crate::source_types::SourceFile::default()])
            .unwrap();
        assert_eq!(
            OutputINCHI_VersionAndKind(
                &mut heap,
                &mut output,
                &mut string_buffer,
                output_options,
                is_beta,
                line_feed.as_const(),
                tab.as_const(),
                stdout,
            ),
            Ok(0)
        );
        let output_bytes = heap.slice(output.s.pStr.as_const()).unwrap()
            [..usize::try_from(output.s.nUsedLength).unwrap()]
            .to_vec();
        let version_bytes = heap.slice(string_buffer.pStr.as_const()).unwrap()
            [..usize::try_from(string_buffer.nUsedLength).unwrap()]
            .to_vec();
        (output_bytes, version_bytes)
    }

    fn formula_wrapper_case(
        n_tag: i32,
        reconnected: i32,
        formula: &str,
        merge: bool,
        overflow: i32,
    ) -> (i32, String, String, INCHI_OUT_CTL) {
        let mut heap = SourceHeap::default();
        let formula_pointer = tag_text(&mut heap, formula).as_mut();
        let inchi = heap
            .allocate_model_storage(vec![crate::source_types::INChI {
                nNumberOfAtoms: 1,
                lenTautomer: 1,
                szHillFormula: formula_pointer,
                ..crate::source_types::INChI::default()
            }])
            .unwrap();
        let sorts = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [SourceMutPointer::null(), inchi],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        let string_storage = heap.allocate_model_storage(vec![0_i8; 128]).unwrap();
        let mut string_buffer = INCHI_IOS_STRING {
            pStr: string_storage,
            nAllocatedLength: 128,
            nUsedLength: 0,
            nPtr: 32,
        };
        let output_storage = heap.allocate_model_storage(vec![0_i8; 128]).unwrap();
        let mut output = INCHI_IOSTREAM {
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: output_storage,
                nAllocatedLength: 128,
                nUsedLength: 0,
                nPtr: 32,
            },
            ..INCHI_IOSTREAM::default()
        };
        let line_feed = tag_text(&mut heap, "\n");
        let tab = tag_text(&mut heap, "\n");
        let stdout = heap
            .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
            .unwrap();
        let mut io = INCHI_OUT_CTL {
            nTag: n_tag,
            bOutType: crate::source_types::OUT_T1 as i32,
            num_components: 1,
            bUseMulipliers: 1,
            bOverflow: overflow,
            n_pzz: i32::from(merge),
            n_zy: i32::from(merge),
            pINChISort: sorts,
            ..INCHI_OUT_CTL::default()
        };
        let result = OutputINCHI_MainLayerFormula(
            &mut heap,
            SourceMutPointer::null(),
            &mut output,
            &mut string_buffer,
            &[1, 0],
            &reconnected,
            &mut io,
            line_feed,
            tab,
            stdout,
        )
        .unwrap();
        let string = tag_result(&heap, string_buffer.pStr);
        let output_bytes = tag_result(&heap, output.s.pStr);
        (
            result,
            String::from_utf8(
                string[..string.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap(),
            String::from_utf8(
                output_bytes[..output_bytes.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap(),
            io,
        )
    }

    fn connections_wrapper_case(
        n_tag: i32,
        table: &[u16],
        plain_text_tags: i32,
        overflow: i32,
        number_of_components: i32,
        fail_after_allocations: Option<u64>,
    ) -> (i32, String, String, INCHI_OUT_CTL) {
        let mut heap = SourceHeap::default();
        let table_pointer = heap.allocate_model_storage(table.to_vec()).unwrap();
        let inchi = heap
            .allocate_model_storage(vec![crate::source_types::INChI {
                nNumberOfAtoms: 3,
                nConnTable: table_pointer,
                lenConnTable: i32::try_from(table.len()).unwrap(),
                ..crate::source_types::INChI::default()
            }])
            .unwrap();
        let sort = crate::source_types::INCHI_SORT {
            pINChI: [SourceMutPointer::null(), inchi],
            ..crate::source_types::INCHI_SORT::default()
        };
        let sorts = heap
            .allocate_model_storage(vec![sort; usize::try_from(number_of_components).unwrap()])
            .unwrap();
        let string_storage = heap.allocate_model_storage(vec![b'x' as i8; 128]).unwrap();
        let mut string_buffer = INCHI_IOS_STRING {
            pStr: string_storage,
            nAllocatedLength: 128,
            nUsedLength: 1,
            nPtr: 32,
        };
        let output_storage = heap.allocate_model_storage(vec![0_i8; 128]).unwrap();
        let mut output = INCHI_IOSTREAM {
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: output_storage,
                nAllocatedLength: 128,
                nUsedLength: 0,
                nPtr: 32,
            },
            ..INCHI_IOSTREAM::default()
        };
        let line_feed = tag_text(&mut heap, "\n");
        let stdout = heap
            .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
            .unwrap();
        let mut io = INCHI_OUT_CTL {
            nTag: n_tag,
            bOutType: crate::source_types::OUT_T1 as i32,
            ATOM_MODE: 0,
            num_components: number_of_components,
            bUseMulipliers: 1,
            bPlainTextTags: plain_text_tags,
            bOverflow: overflow,
            pINChISort: sorts,
            ..INCHI_OUT_CTL::default()
        };
        if let Some(successful_allocations) = fail_after_allocations {
            heap.fail_after_allocations(successful_allocations);
        }
        let result = OutputINCHI_MainLayerConnections(
            &mut heap,
            SourceMutPointer::null(),
            &mut output,
            &mut string_buffer,
            &[1, 0],
            &0,
            &mut io,
            line_feed,
            line_feed,
            stdout,
        )
        .unwrap();
        let string = tag_result(&heap, string_buffer.pStr);
        let output_bytes = tag_result(&heap, output.s.pStr);
        (
            result,
            String::from_utf8(
                string[..string.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap(),
            String::from_utf8(
                output_bytes[..output_bytes.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap(),
            io,
        )
    }

    fn hydrogens_wrapper_case(
        segment_difference: i8,
        n_tag: i32,
        hydrogens: &[i8],
    ) -> (i32, String, String, INCHI_OUT_CTL) {
        let mut heap = SourceHeap::default();
        let hydrogen_pointer = heap.allocate_model_storage(hydrogens.to_vec()).unwrap();
        let inchi = heap
            .allocate_model_storage(vec![crate::source_types::INChI {
                nNumberOfAtoms: i32::try_from(hydrogens.len()).unwrap(),
                nNum_H: hydrogen_pointer,
                ..crate::source_types::INChI::default()
            }])
            .unwrap();
        let sorts = heap
            .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                pINChI: [SourceMutPointer::null(), inchi],
                ..crate::source_types::INCHI_SORT::default()
            }])
            .unwrap();
        let string_storage = heap.allocate_model_storage(vec![b'x' as i8; 128]).unwrap();
        heap.slice_mut(string_storage).unwrap()[1] = 0;
        let mut string_buffer = INCHI_IOS_STRING {
            pStr: string_storage,
            nAllocatedLength: 128,
            nUsedLength: 1,
            nPtr: 32,
        };
        let output_storage = heap.allocate_model_storage(vec![0_i8; 128]).unwrap();
        let mut output = INCHI_IOSTREAM {
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: output_storage,
                nAllocatedLength: 128,
                nUsedLength: 0,
                nPtr: 32,
            },
            ..INCHI_IOSTREAM::default()
        };
        let line_feed = tag_text(&mut heap, "\n");
        let stdout = heap
            .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
            .unwrap();
        let mut io = INCHI_OUT_CTL {
            nTag: n_tag,
            bOutType: crate::source_types::OUT_T1 as i32,
            num_components: 1,
            bUseMulipliers: 1,
            pINChISort: sorts,
            bTag1: -71,
            bAlways: -72,
            tot_len: -73,
            tot_len2: -74,
            ..INCHI_OUT_CTL::default()
        };
        io.sDifSegs[0][crate::source_types::tagDiffINChISegments_DIFS_h_H_ATOMS as usize] =
            segment_difference;
        let result = OutputINCHI_MainLayerHydrogens(
            &mut heap,
            SourceMutPointer::null(),
            &mut output,
            &mut string_buffer,
            &[1, 0],
            &0,
            &mut io,
            line_feed,
            line_feed,
            stdout,
        )
        .unwrap();
        let string = tag_result(&heap, string_buffer.pStr);
        let output_bytes = tag_result(&heap, output.s.pStr);
        (
            result,
            String::from_utf8(
                string[..string.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap(),
            String::from_utf8(
                output_bytes[..output_bytes.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap(),
            io,
        )
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_mainlayerformula__line_3170() {
        for (n_tag, reconnect, expected, always) in [
            (2, 0, "/C2H6", 0),
            (2, INCHI_REC as i32, "/rC2H6", 0),
            (1, 0, "formulaC2H6", 1),
            (3, 0, "/{formula}C2H6", 1),
        ] {
            let (result, string, output, io) =
                formula_wrapper_case(n_tag, reconnect, "C2H6", false, 0);
            assert_eq!(result, 0);
            assert_eq!(string, expected);
            assert_eq!(output, format!("{expected}\n"));
            assert_eq!(io.tot_len, 4);
            assert_eq!(io.bAlways, always);
        }
        let (result, string, output, _) = formula_wrapper_case(2, 0, "ZzZz", true, 0);
        assert_eq!(result, 0);
        assert_eq!(string, "/Zz2");
        assert_eq!(output, "/Zz2\n");

        let (result, _, output, io) = formula_wrapper_case(2, 0, "C2H6", false, 1);
        assert_eq!(result, 1);
        assert_eq!(output, "");
        assert_eq!(io.tot_len, 0);

        let mut heap = SourceHeap::default();
        let storage = heap
            .allocate_model_storage(vec![b'x' as i8, 0, 0, 0])
            .unwrap();
        let mut string_buffer = INCHI_IOS_STRING {
            pStr: storage,
            nAllocatedLength: 4,
            nUsedLength: 1,
            nPtr: 2,
        };
        let output_storage = heap.allocate_model_storage(vec![0_i8; 8]).unwrap();
        let mut output = INCHI_IOSTREAM {
            type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
            s: INCHI_IOS_STRING {
                pStr: output_storage,
                nAllocatedLength: 8,
                nPtr: 4,
                ..INCHI_IOS_STRING::default()
            },
            ..INCHI_IOSTREAM::default()
        };
        let mut io = INCHI_OUT_CTL::default();
        let lf = tag_text(&mut heap, "\n");
        let stdout = heap
            .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
            .unwrap();
        assert_eq!(
            OutputINCHI_MainLayerFormula(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut string_buffer,
                &[0, 0],
                &0,
                &mut io,
                lf,
                lf,
                stdout,
            ),
            Ok(0)
        );
        assert_eq!(string_buffer.nUsedLength, 1);
        assert_eq!(
            heap.slice(string_buffer.pStr.as_const()).unwrap()[0],
            b'x' as i8
        );
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_mainlayerconnections__line_3213() {
        for (n_tag, expected, always) in [
            (2, "/c1-2-3", 0),
            (1, "connections1-2-3", 1),
            (3, "/c{connections}1-2-3", 1),
        ] {
            let (result, string, output, io) =
                connections_wrapper_case(n_tag, &[2, 1, 3, 2], 1, 0, 1, None);
            assert_eq!(result, 0);
            assert_eq!(string, expected);
            assert_eq!(output, format!("{expected}\n"));
            assert_eq!(io.tot_len, 5);
            assert_eq!(io.tot_len2, 5);
            assert_eq!(io.bTag1, 1 << 5);
            assert_eq!(io.bAlways, always);
        }

        let (result, string, output, io) =
            connections_wrapper_case(2, &[2, 1, 3, 2], 0, 0, 1, None);
        assert_eq!(
            (result, string.as_str(), output.as_str()),
            (0, "1-2-3", "1-2-3\n")
        );
        assert_eq!((io.tot_len, io.tot_len2), (5, 5));

        let (result, string, output, io) = connections_wrapper_case(2, &[], 1, 0, 1, None);
        assert_eq!((result, string.as_str(), output.as_str()), (0, "", ""));
        assert_eq!((io.tot_len, io.tot_len2), (0, 0));

        let (result, string, output, io) =
            connections_wrapper_case(2, &[2, 1, 3, 2], 1, 0, 2, Some(0));
        assert_eq!(
            (
                result,
                string.as_str(),
                output.as_str(),
                io.tot_len,
                io.tot_len2,
                io.bOverflow
            ),
            (1, "2*", "", 2, 2, 1)
        );
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_mainlayerhydrogens__line_3251() {
        for (n_tag, expected, always) in [(2, "/h1H", -72), (1, "H1H", 1), (3, "/h{H_atoms}1H", 1)]
        {
            let (result, string, output, io) = hydrogens_wrapper_case(11, n_tag, &[1, 0]);
            assert_eq!(result, 0);
            assert_eq!(string, expected);
            assert_eq!(output, format!("{expected}\n"));
            assert_eq!((io.tot_len, io.tot_len2), (2, 2));
            assert_eq!(io.bTag1, 1 << 6);
            assert_eq!(io.bAlways, always);
        }

        let (result, string, output, io) = hydrogens_wrapper_case(11, 2, &[0, 0]);
        assert_eq!((result, string.as_str(), output.as_str()), (0, "", ""));
        assert_eq!((io.tot_len, io.tot_len2, io.bTag1), (0, 0, 1 << 6));

        for segment_difference in [0, 4] {
            let (result, string, output, io) =
                hydrogens_wrapper_case(segment_difference, 2, &[1, 0]);
            assert_eq!((result, string.as_str(), output.as_str()), (0, "x", ""));
            assert_eq!((io.tot_len, io.tot_len2), (-73, -74));
            assert_eq!((io.bTag1, io.bAlways), (-71, -72));
        }
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_chargeandremovedaddedprotonslayers__line_3293() {
        fn run_case(
            charge_difference: i8,
            proton_difference: i8,
            n_tag: i32,
            plain_text_tags: i32,
            taut_mode: i32,
            second_pass: i32,
            overflow: i32,
            fixed_h_tag: i32,
        ) -> (i32, String, String, INCHI_OUT_CTL) {
            let mut heap = SourceHeap::default();
            let inchi = heap
                .allocate_model_storage(vec![crate::source_types::INChI {
                    nNumberOfAtoms: 1,
                    lenTautomer: 1,
                    nTotalCharge: 2,
                    ..crate::source_types::INChI::default()
                }])
                .unwrap();
            let sorts = heap
                .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), inchi],
                    ..crate::source_types::INCHI_SORT::default()
                }])
                .unwrap();
            let string_storage = heap.allocate_model_storage(vec![b'x' as i8; 128]).unwrap();
            heap.slice_mut(string_storage).unwrap()[1] = 0;
            let mut string_buffer = INCHI_IOS_STRING {
                pStr: string_storage,
                nAllocatedLength: 128,
                nUsedLength: 1,
                nPtr: 32,
            };
            let output_storage = heap.allocate_model_storage(vec![0_i8; 256]).unwrap();
            let mut output = INCHI_IOSTREAM {
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                s: INCHI_IOS_STRING {
                    pStr: output_storage,
                    nAllocatedLength: 256,
                    nUsedLength: 0,
                    nPtr: 32,
                },
                ..INCHI_IOSTREAM::default()
            };
            let line_feed = tag_text(&mut heap, "\n");
            let stdout = heap
                .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
                .unwrap();
            let mut io = INCHI_OUT_CTL {
                bOutType: crate::source_types::OUT_T1 as i32,
                bPlainTextTags: plain_text_tags,
                bOmitRepetitions: 0,
                bUseMulipliers: 1,
                bNonTautNonIsoIdentifierNotEmpty: 5,
                bSecondNonTautPass: second_pass,
                bFhTag: fixed_h_tag,
                iCurTautMode: taut_mode,
                num_components: 1,
                nNumRemovedProtons: -3,
                nTag: n_tag,
                bOverflow: overflow,
                pINChISort: sorts,
                pINChISort2: sorts,
                ..INCHI_OUT_CTL::default()
            };
            io.sDifSegs[0][crate::source_types::tagDiffINChISegments_DIFS_q_CHARGE as usize] =
                charge_difference;
            io.sDifSegs[0][crate::source_types::tagDiffINChISegments_DIFS_p_PROTONS as usize] =
                proton_difference;
            let result = OutputINCHI_ChargeAndRemovedAddedProtonsLayers(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut string_buffer,
                &mut io,
                line_feed,
                line_feed,
                stdout,
            )
            .unwrap();
            let string = tag_result(&heap, string_buffer.pStr);
            let output_bytes = tag_result(&heap, output.s.pStr);
            (
                result,
                String::from_utf8(
                    string[..string.len() - 1]
                        .iter()
                        .map(|byte| *byte as u8)
                        .collect(),
                )
                .unwrap(),
                String::from_utf8(
                    output_bytes[..output_bytes.len() - 1]
                        .iter()
                        .map(|byte| *byte as u8)
                        .collect(),
                )
                .unwrap(),
                io,
            )
        }

        let fill = 11;
        let empty = 4;
        let omit = 0;
        let taut = crate::source_types::TAUT_YES as i32;

        let (result, string, output, io) = run_case(fill, fill, 2, 1, taut, 0, 0, 0);
        assert_eq!((result, string.as_str()), (0, "/p-3"));
        assert_eq!(output, "/q+2\n/p-3\n");
        assert_eq!(
            (io.nSegmAction, io.bTag1, io.tot_len, io.bAlways),
            (1, 1 << 8, 0, 0)
        );
        assert_eq!(io.bNonTautNonIsoIdentifierNotEmpty, 5);

        let (result, string, output, io) = run_case(empty, omit, 2, 1, 0, 0, 0, 0);
        assert_eq!(
            (result, string.as_str(), output.as_str()),
            (0, "/q", "/q\n")
        );
        assert_eq!((io.nSegmAction, io.tot_len), (2, 0));

        let (result, string, output, io) = run_case(omit, omit, 2, 1, taut, 0, 0, 0);
        assert_eq!((result, string.as_str(), output.as_str()), (0, "x", "/"));
        assert_eq!(io.nSegmAction, 0);

        let (result, string, output, io) = run_case(omit, empty, 2, 1, taut, 0, 0, 0);
        assert_eq!(
            (result, string.as_str(), output.as_str()),
            (0, "/p-3", "/p-3\n")
        );
        assert_eq!((io.nSegmAction, io.tot_len), (2, 0));

        let fixed_h = crate::source_types::local_ichiprt1::tagIdentLblBit_IL_FIXH as i32;
        let (result, string, output, io) = run_case(fill, omit, 3, 1, 0, 0, 0, fixed_h);
        assert_eq!(
            (result, string.as_str(), output.as_str()),
            (0, "/q{fixed_H:charge}+2", "/q{fixed_H:charge}+2\n")
        );
        assert_eq!(io.bTag1, fixed_h | (1 << 7));

        let (result, string, output, io) = run_case(fill, fill, 2, 1, taut, 1, 0, 0);
        assert_eq!(
            (result, string.as_str(), output.as_str()),
            (0, "/q+2", "/q+2\n")
        );
        assert_eq!(io.bNonTautNonIsoIdentifierNotEmpty, 6);
        assert_eq!(io.nSegmAction, 1);

        let (result, string, output, io) = run_case(fill, fill, 2, 1, taut, 0, 1, 0);
        assert_eq!((result, string.as_str(), output.as_str()), (1, "+2", ""));
        assert_eq!(io.bOverflow, 1);
        assert_eq!(io.bNonTautNonIsoIdentifierNotEmpty, 5);
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_stereolayer__line_3354() {
        fn run_case(
            actions: [i8; 4],
            plain_text_tags: i32,
            relative: i32,
            racemic: i32,
            overflow: i32,
        ) -> (i32, String, String, INCHI_OUT_CTL) {
            let mut heap = SourceHeap::default();
            let bond_atom1 = heap.allocate_model_storage(vec![1_u16]).unwrap();
            let bond_atom2 = heap.allocate_model_storage(vec![2_u16]).unwrap();
            let bond_parity = heap.allocate_model_storage(vec![1_i8]).unwrap();
            let centers = heap.allocate_model_storage(vec![3_u16]).unwrap();
            let center_parity = heap.allocate_model_storage(vec![2_i8]).unwrap();
            let stereo = heap
                .allocate_model_storage(vec![crate::source_types::INChI_Stereo {
                    nNumberOfStereoBonds: 1,
                    nBondAtom1: bond_atom1,
                    nBondAtom2: bond_atom2,
                    b_parity: bond_parity,
                    nNumberOfStereoCenters: 1,
                    nNumber: centers,
                    t_parity: center_parity,
                    nCompInv2Abs: -1,
                    ..crate::source_types::INChI_Stereo::default()
                }])
                .unwrap();
            let inchi = heap
                .allocate_model_storage(vec![crate::source_types::INChI {
                    nNumberOfAtoms: 1,
                    lenTautomer: 1,
                    Stereo: stereo,
                    ..crate::source_types::INChI::default()
                }])
                .unwrap();
            let sorts = heap
                .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), inchi],
                    ..crate::source_types::INCHI_SORT::default()
                }])
                .unwrap();
            let string_storage = heap.allocate_model_storage(vec![0_i8; 256]).unwrap();
            let mut string_buffer = INCHI_IOS_STRING {
                pStr: string_storage,
                nAllocatedLength: 256,
                nUsedLength: 0,
                nPtr: 32,
            };
            let output_storage = heap.allocate_model_storage(vec![0_i8; 1024]).unwrap();
            let mut output = INCHI_IOSTREAM {
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                s: INCHI_IOS_STRING {
                    pStr: output_storage,
                    nAllocatedLength: 1024,
                    nUsedLength: 0,
                    nPtr: 32,
                },
                ..INCHI_IOSTREAM::default()
            };
            let line_feed = tag_text(&mut heap, "\n");
            let stdout = heap
                .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
                .unwrap();
            let mut io = INCHI_OUT_CTL {
                TAUT_MODE: 0,
                bOutType: crate::source_types::OUT_T1 as i32,
                bPlainTextTags: plain_text_tags,
                bOmitRepetitions: 0,
                bUseMulipliers: 1,
                bNonTautNonIsoIdentifierNotEmpty: 5,
                bSecondNonTautPass: 1,
                iCurTautMode: crate::source_types::TAUT_YES as i32,
                num_components: 1,
                nTag: 2,
                bOverflow: overflow,
                pINChISort: sorts,
                pINChISort2: sorts,
                ..INCHI_OUT_CTL::default()
            };
            io.bRelativeStereo[crate::source_types::TAUT_YES as usize] = relative;
            io.bRacemicStereo[crate::source_types::TAUT_YES as usize] = racemic;
            for (segment, action) in [
                crate::source_types::tagDiffINChISegments_DIFS_b_SBONDS as usize,
                crate::source_types::tagDiffINChISegments_DIFS_t_SATOMS as usize,
                crate::source_types::tagDiffINChISegments_DIFS_m_SP3INV as usize,
                crate::source_types::tagDiffINChISegments_DIFS_s_STYPE as usize,
            ]
            .into_iter()
            .zip(actions)
            {
                io.sDifSegs[0][segment] = action;
            }
            let result = OutputINCHI_StereoLayer(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut string_buffer,
                &mut io,
                line_feed,
                line_feed,
                stdout,
            )
            .unwrap();
            let decode = |bytes: Vec<i8>| {
                String::from_utf8(
                    bytes[..bytes.len() - 1]
                        .iter()
                        .map(|byte| *byte as u8)
                        .collect(),
                )
                .unwrap()
            };
            (
                result,
                decode(tag_result(&heap, string_buffer.pStr)),
                decode(tag_result(&heap, output.s.pStr)),
                io,
            )
        }

        let fill = 11;
        let empty = 4;
        let omit = 0;
        let (result, string, output, io) = run_case([fill; 4], 1, 1, 0, 0);
        assert_eq!(result, 0);
        assert_eq!(string, "/s2");
        assert_eq!(output, "/b1-2-\n/t3+\n/m1\n/s2\n/");
        assert_eq!(io.bNonTautNonIsoIdentifierNotEmpty, 9);
        assert_eq!(io.bRelRac, 1);
        assert_eq!(io.nSegmAction, 1);

        assert_eq!(
            run_case([fill; 4], 1, 0, 1, 0).2,
            "/b1-2-\n/t3+\n/m1\n/s3\n/"
        );
        assert_eq!(
            run_case([fill; 4], 1, 0, 0, 0).2,
            "/b1-2-\n/t3+\n/m1\n/s1\n/"
        );
        assert_eq!(run_case([omit; 4], 1, 0, 0, 0).2, "////");
        assert_eq!(
            run_case([fill, omit, omit, omit], 1, 0, 0, 0).2,
            "/b1-2-\n///"
        );
        let empty_case = run_case([empty; 4], 1, 0, 0, 0);
        assert_eq!(empty_case.2, "/b\n/t\n/m\n/s\n/");
        assert_eq!(empty_case.3.bNonTautNonIsoIdentifierNotEmpty, 5);

        for (actions, expected) in [
            ([fill, omit, omit, omit], 1),
            ([omit, fill, omit, omit], 2),
            ([omit, omit, fill, omit], 3),
            ([omit, omit, omit, fill], 1),
        ] {
            assert_eq!(run_case(actions, 1, 0, 0, 1).0, expected);
        }
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_isotopiclayer__line_3502() {
        fn run_case(
            actions: [i8; 6],
            plain_text_tags: i32,
            relative: i32,
            racemic: i32,
            overflow: i32,
            transposition: bool,
        ) -> (i32, String, INCHI_OUT_CTL, i32) {
            let mut heap = SourceHeap::default();
            let mut sorts = SourceMutPointer::null();
            let mut sorts2 = SourceMutPointer::null();
            let mut components = 0;
            if transposition {
                let non = heap
                    .allocate_model_storage(vec![crate::source_types::INChI {
                        nNumberOfAtoms: 1,
                        lenTautomer: 0,
                        ..crate::source_types::INChI::default()
                    }])
                    .unwrap();
                let taut = heap
                    .allocate_model_storage(vec![crate::source_types::INChI {
                        nNumberOfAtoms: 1,
                        lenTautomer: 1,
                        ..crate::source_types::INChI::default()
                    }])
                    .unwrap();
                sorts = heap
                    .allocate_model_storage(
                        [1_i16, 0_i16]
                            .into_iter()
                            .map(|ord_number| crate::source_types::INCHI_SORT {
                                pINChI: [non, SourceMutPointer::null()],
                                ord_number,
                                ..crate::source_types::INCHI_SORT::default()
                            })
                            .collect(),
                    )
                    .unwrap();
                sorts2 = heap
                    .allocate_model_storage(
                        [0_i16, 1_i16]
                            .into_iter()
                            .map(|ord_number| crate::source_types::INCHI_SORT {
                                pINChI: [SourceMutPointer::null(), taut],
                                ord_number,
                                ..crate::source_types::INCHI_SORT::default()
                            })
                            .collect(),
                    )
                    .unwrap();
                components = 2;
            }
            let string_storage = heap.allocate_model_storage(vec![0_i8; 512]).unwrap();
            let mut string_buffer = INCHI_IOS_STRING {
                pStr: string_storage,
                nAllocatedLength: 512,
                nUsedLength: 0,
                nPtr: 32,
            };
            let output_storage = heap.allocate_model_storage(vec![0_i8; 2048]).unwrap();
            let mut output = INCHI_IOSTREAM {
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                s: INCHI_IOS_STRING {
                    pStr: output_storage,
                    nAllocatedLength: 2048,
                    nUsedLength: 0,
                    nPtr: 32,
                },
                ..INCHI_IOSTREAM::default()
            };
            let line_feed = tag_text(&mut heap, "\n");
            let stdout = heap
                .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
                .unwrap();
            let flags = heap.allocate_model_storage(vec![0_i32]).unwrap();
            let mut io = INCHI_OUT_CTL {
                TAUT_MODE: 0,
                pSortPrintINChIFlags: flags,
                bOverflow: overflow,
                bOutputType: if transposition {
                    crate::source_types::OUT_TN as i32
                } else {
                    0
                },
                bOutType: if transposition {
                    crate::source_types::local_ichiprt1::OUT_NONTAUT as i32
                } else {
                    crate::source_types::OUT_T1 as i32
                },
                bPlainTextTags: plain_text_tags,
                bSecondNonTautPass: i32::from(transposition),
                iCurTautMode: crate::source_types::TAUT_YES as i32,
                num_components: components,
                num_iso_H: [1, 2, 3],
                nTag: 2,
                pINChISort: sorts,
                pINChISort2: sorts2,
                ..INCHI_OUT_CTL::default()
            };
            io.bIsotopicRelativeStereo[crate::source_types::TAUT_YES as usize] = relative;
            io.bIsotopicRacemicStereo[crate::source_types::TAUT_YES as usize] = racemic;
            for (segment, value) in [
                (
                    crate::source_types::tagDiffINChISegments_DIFS_i_IATOMS,
                    actions[0],
                ),
                (
                    crate::source_types::tagDiffINChISegments_DIFS_h_H_ATOMS,
                    actions[1],
                ),
                (
                    crate::source_types::tagDiffINChISegments_DIFS_b_SBONDS,
                    actions[2],
                ),
                (
                    crate::source_types::tagDiffINChISegments_DIFS_t_SATOMS,
                    actions[3],
                ),
                (
                    crate::source_types::tagDiffINChISegments_DIFS_m_SP3INV,
                    actions[4],
                ),
                (
                    crate::source_types::tagDiffINChISegments_DIFS_s_STYPE,
                    actions[5],
                ),
            ] {
                io.sDifSegs[0][segment as usize] = value;
            }
            if transposition {
                io.sDifSegs[crate::source_types::tagDiffINChILayers_DIFL_F as usize]
                    [crate::source_types::tagDiffINChISegments_DIFS_o_TRANSP as usize] = 11;
            }
            let result = OutputINCHI_IsotopicLayer(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut string_buffer,
                &0,
                &mut io,
                line_feed,
                line_feed,
                stdout,
            )
            .unwrap();
            let bytes = tag_result(&heap, output.s.pStr);
            let text = String::from_utf8(
                bytes[..bytes.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap();
            let flags = heap.slice(flags.as_const()).unwrap()[0];
            (result, text, io, flags)
        }

        let fill = 11;
        let empty = 4;
        let omit = 0;
        assert_eq!(run_case([omit; 6], 1, 0, 0, 0, false).1, "/////");
        assert_eq!(
            run_case([empty, omit, omit, omit, omit, omit], 1, 0, 0, 0, false).1,
            "/i\n////"
        );
        let filled = run_case([fill; 6], 1, 1, 0, 0, false);
        assert_eq!(filled.0, 0);
        assert!(filled.1.contains("/hT3D2H"));
        assert!(filled.1.contains("/s2\n/"));
        assert_eq!(filled.2.bRelRac, 1);
        assert_eq!(run_case([fill; 6], 1, 0, 1, 0, false).2.bRelRac, 1);
        assert_eq!(run_case([fill; 6], 1, 0, 0, 1, false).0, 1);

        let transposed = run_case([omit; 6], 1, 0, 0, 0, true);
        assert_eq!(transposed.0, 0);
        assert_eq!(
            transposed.3,
            crate::source_types::FLAG_SORT_PRINT_TRANSPOS_BAS as i32
        );
        assert!(transposed.1.ends_with("/o(1,2)\n"));
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_fixedhlayerwithsublayers__line_3750() {
        struct CaseResult {
            status: i32,
            output: String,
            io: INCHI_OUT_CTL,
            flags: i32,
            repeat: i32,
            fixed_sorts: SourceMutPointer<crate::source_types::INCHI_SORT>,
            taut_sorts: SourceMutPointer<crate::source_types::INCHI_SORT>,
        }

        fn run_case(
            basic_or_reconnected: i32,
            out_type: i32,
            output_type: i32,
            identical: i32,
            tautomeric: i32,
            non_tautomeric: i32,
            second_pass: i32,
            formula_difference: i8,
            fixed_h_difference: i8,
            non_taut_non_iso_not_empty: i32,
            non_taut_iso_not_empty: i32,
            initial_overflow: i32,
            n_pzz: i32,
            n_zy: i32,
            fail_fixed_h_allocation: bool,
        ) -> CaseResult {
            let mut heap = SourceHeap::default();
            let fixed_formula = tag_text(&mut heap, "C2H6");
            let taut_formula = tag_text(&mut heap, "CH4");
            let fixed_h = heap
                .allocate_model_storage(vec![
                    if fail_fixed_h_allocation {
                        127_i8
                    } else {
                        1_i8
                    },
                    0,
                ])
                .unwrap();
            let fixed = heap
                .allocate_model_storage(vec![crate::source_types::INChI {
                    nNumberOfAtoms: 2,
                    lenTautomer: 0,
                    szHillFormula: fixed_formula.as_mut(),
                    nNum_H_fixed: fixed_h,
                    ..crate::source_types::INChI::default()
                }])
                .unwrap();
            let taut = heap
                .allocate_model_storage(vec![crate::source_types::INChI {
                    nNumberOfAtoms: 1,
                    lenTautomer: 1,
                    szHillFormula: taut_formula.as_mut(),
                    ..crate::source_types::INChI::default()
                }])
                .unwrap();
            let fixed_sorts = heap
                .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                    pINChI: [fixed, SourceMutPointer::null()],
                    ..crate::source_types::INCHI_SORT::default()
                }])
                .unwrap();
            let taut_sorts = heap
                .allocate_model_storage(vec![crate::source_types::INCHI_SORT {
                    pINChI: [SourceMutPointer::null(), taut],
                    ..crate::source_types::INCHI_SORT::default()
                }])
                .unwrap();
            let sort_sets = heap
                .allocate_model_storage(vec![fixed_sorts, taut_sorts])
                .unwrap();
            let flags_pointer = heap.allocate_model_storage(vec![0_i32]).unwrap();
            let string_storage = heap.allocate_model_storage(vec![0_i8; 256]).unwrap();
            let mut string_buffer = INCHI_IOS_STRING {
                pStr: string_storage,
                nAllocatedLength: 256,
                nUsedLength: 0,
                nPtr: 32,
            };
            let output_storage = heap.allocate_model_storage(vec![0_i8; 1024]).unwrap();
            let mut output = INCHI_IOSTREAM {
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                s: INCHI_IOS_STRING {
                    pStr: output_storage,
                    nAllocatedLength: 1024,
                    nUsedLength: 0,
                    nPtr: 32,
                },
                ..INCHI_IOSTREAM::default()
            };
            let line_feed = tag_text(&mut heap, "\n");
            let stdout = heap
                .allocate_model_storage(vec![crate::source_types::SourceFile::default()])
                .unwrap();
            let mut io = INCHI_OUT_CTL {
                ATOM_MODE: if fail_fixed_h_allocation {
                    crate::source_types::CT_MODE_EQL_H_TOGETHER as i32
                } else {
                    0
                },
                pSortPrintINChIFlags: flags_pointer,
                bOverflow: initial_overflow,
                bOutputType: output_type,
                bOutType: out_type,
                bPlainTextTags: 1,
                bUseMulipliers: 1,
                bNonTautNonIsoIdentifierNotEmpty: non_taut_non_iso_not_empty,
                bNonTautIsoIdentifierNotEmpty: non_taut_iso_not_empty,
                bSecondNonTautPass: second_pass,
                bTautomeric: tautomeric,
                bNonTautomeric: non_tautomeric,
                bNonTautIsIdenticalToTaut: identical,
                bFhTag: 91,
                iCurTautMode: crate::source_types::TAUT_YES as i32,
                num_components: 7,
                nTag: 2,
                num_comp: [1, 3],
                n_pzz,
                n_zy,
                pINChISortTautAndNonTaut: sort_sets,
                pINChISort: if second_pass != 0 {
                    fixed_sorts
                } else {
                    taut_sorts
                },
                pINChISort2: taut_sorts,
                ..INCHI_OUT_CTL::default()
            };
            io.sDifSegs[crate::source_types::tagDiffINChILayers_DIFL_F as usize]
                [crate::source_types::tagDiffINChISegments_DIFS_f_FORMULA as usize] =
                formula_difference;
            io.sDifSegs[crate::source_types::tagDiffINChILayers_DIFL_F as usize]
                [crate::source_types::tagDiffINChISegments_DIFS_h_H_ATOMS as usize] =
                fixed_h_difference;
            let mut repeat = 99;
            if fail_fixed_h_allocation {
                heap.fail_after_allocations(0);
            }
            let status = OutputINCHI_FixedHLayerWithSublayers(
                &mut heap,
                SourceMutPointer::null(),
                &mut output,
                &mut string_buffer,
                &basic_or_reconnected,
                &mut io,
                line_feed,
                line_feed,
                &mut repeat,
                stdout,
            )
            .unwrap();
            let bytes = tag_result(&heap, output.s.pStr);
            let output = String::from_utf8(
                bytes[..bytes.len() - 1]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap();
            let flags = heap.slice(flags_pointer.as_const()).unwrap()[0];
            CaseResult {
                status,
                output,
                io,
                flags,
                repeat,
                fixed_sorts,
                taut_sorts,
            }
        }

        let out_tn = crate::source_types::OUT_TN as i32;
        let out_non_taut = crate::source_types::local_ichiprt1::OUT_NONTAUT as i32;
        let fill = 11_i8;
        let omit = 0_i8;

        let identical_basic = run_case(
            crate::source_types::INCHI_BAS as i32,
            out_tn,
            out_tn,
            1,
            1,
            1,
            0,
            omit,
            omit,
            0,
            0,
            0,
            0,
            0,
            false,
        );
        assert_eq!(identical_basic.status, 0);
        assert_eq!(identical_basic.repeat, 0);
        assert_eq!(identical_basic.output, "");
        assert_eq!(
            identical_basic.flags,
            (crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_BAS
                | crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_BAS) as i32
        );
        let identical_reconnected = run_case(
            crate::source_types::INCHI_REC as i32,
            out_tn,
            out_tn,
            1,
            1,
            1,
            0,
            omit,
            omit,
            0,
            0,
            0,
            0,
            0,
            false,
        );
        assert_eq!(
            identical_reconnected.flags,
            (crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_REC
                | crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_REC) as i32
        );
        assert_eq!(
            run_case(
                0, out_tn, out_tn, 1, 0, 1, 0, omit, omit, 0, 0, 0, 0, 0, false,
            )
            .flags,
            0
        );

        let filled = run_case(
            crate::source_types::INCHI_BAS as i32,
            out_tn,
            out_tn,
            0,
            0,
            0,
            0,
            fill,
            fill,
            4,
            0,
            0,
            0,
            0,
            false,
        );
        assert_eq!(filled.status, 0);
        assert_eq!(filled.repeat, 1);
        assert_eq!(filled.output, "/fC2H6\n/h1H\n");
        assert_eq!(filled.io.bOutType, out_non_taut);
        assert_eq!(filled.io.iCurTautMode, crate::source_types::TAUT_NON as i32);
        assert_eq!(filled.io.pINChISort, filled.fixed_sorts);
        assert_eq!(filled.io.bSecondNonTautPass, 1);
        assert_eq!(
            filled.io.nCurINChISegment,
            crate::source_types::tagDiffINChILayers_DIFL_F as i32
        );
        assert_eq!(filled.io.num_components, 1);
        assert_eq!(
            filled.io.bFhTag,
            crate::source_types::local_ichiprt1::tagIdentLblBit_IL_FIXH as i32
        );
        assert_eq!(filled.io.bNonTautNonIsoIdentifierNotEmpty, 6);
        assert_eq!(filled.io.tot_len, 2);
        assert_eq!(
            filled.io.nSegmAction,
            crate::source_types::tagINChISegmAction_INCHI_SEGM_FILL as i32
        );

        let empty_actions = run_case(0, out_tn, out_tn, 0, 0, 0, 0, 4, omit, 0, 0, 0, 0, 0, false);
        assert_eq!(empty_actions.status, 0);
        assert_eq!(empty_actions.repeat, 1);
        assert_eq!(empty_actions.output, "/f\n");
        assert_eq!(empty_actions.io.tot_len, 0);
        assert_eq!(empty_actions.io.bNonTautNonIsoIdentifierNotEmpty, 0);

        let formula_overflow = run_case(
            0, out_tn, out_tn, 0, 0, 0, 0, fill, fill, 0, 0, 7, 0, 0, false,
        );
        assert_eq!(formula_overflow.status, 1);
        assert_eq!(formula_overflow.repeat, 0);
        assert_eq!(formula_overflow.io.bSecondNonTautPass, 1);
        assert_eq!(formula_overflow.io.bNonTautNonIsoIdentifierNotEmpty, 1);

        let fixed_h_overflow = run_case(
            0, out_tn, out_tn, 0, 0, 0, 0, fill, fill, 0, 0, 0, 0, 0, true,
        );
        assert_eq!(fixed_h_overflow.status, 2);
        assert_eq!(fixed_h_overflow.repeat, 0);
        assert_eq!(fixed_h_overflow.output, "/fC2H6\n");
        assert_eq!(fixed_h_overflow.io.bOverflow, 1);
        assert_eq!(fixed_h_overflow.io.bNonTautNonIsoIdentifierNotEmpty, 1);

        let restored_basic = run_case(
            crate::source_types::INCHI_BAS as i32,
            out_non_taut,
            out_tn,
            0,
            0,
            0,
            1,
            omit,
            omit,
            0,
            0,
            0,
            0,
            0,
            false,
        );
        assert_eq!(restored_basic.status, 0);
        assert_eq!(restored_basic.repeat, 0);
        assert_eq!(restored_basic.io.bOutType, out_tn);
        assert_eq!(
            restored_basic.io.iCurTautMode,
            crate::source_types::TAUT_YES as i32
        );
        assert_eq!(restored_basic.io.pINChISort, restored_basic.taut_sorts);
        assert_eq!(restored_basic.io.bSecondNonTautPass, 0);
        assert_eq!(restored_basic.io.num_components, 3);
        assert_eq!(restored_basic.io.bFhTag, 0);
        assert_eq!(
            restored_basic.flags,
            (crate::source_types::FLAG_SORT_PRINT_NO_NFIX_H_BAS
                | crate::source_types::FLAG_SORT_PRINT_NO_IFIX_H_BAS) as i32
        );

        let restored_nonempty = run_case(
            crate::source_types::INCHI_REC as i32,
            out_non_taut,
            out_tn,
            0,
            0,
            0,
            1,
            omit,
            omit,
            1,
            1,
            0,
            0,
            0,
            false,
        );
        assert_eq!(restored_nonempty.flags, 0);
    }

    #[test]
    fn source_port__ichiprt1__isbondatomnumslesser__line_5007() {
        for (first, second, expected) in [
            ([1, 2], [1, 2], 0),
            ([2, 1], [1, 2], 0),
            ([0, 9], [1, 2], -1),
            ([1, 2], [0, 9], 1),
            ([1, 2], [1, 3], -1),
            ([1, 3], [1, 2], 1),
            ([-5, -1], [-4, i32::MIN], 1),
            ([i32::MIN, i32::MAX], [i32::MIN, i32::MAX - 1], 1),
            ([i32::MIN, i32::MIN], [i32::MIN, i32::MAX], -1),
            ([i32::MAX, i32::MAX], [i32::MAX, i32::MAX], 0),
        ] {
            assert_eq!(IsBondAtomNumsLesser(&first, &second), expected);
            assert_eq!(IsBondAtomNumsLesser(&second, &first), -expected);
            assert_eq!(
                IsBondAtomNumsLesser(&[first[1], first[0]], &second),
                expected
            );
        }
    }

    #[test]
    fn source_port__ichiprt1__editinchi_hidepolymerzz__line_5039() {
        fn stream(heap: &mut SourceHeap, text: &[u8], extra: usize) -> INCHI_IOSTREAM {
            let mut storage = text.iter().map(|byte| *byte as i8).collect::<Vec<_>>();
            storage.push(0);
            storage.resize(storage.len() + extra, 0x55_i8);
            let allocated = storage.len() as i32;
            INCHI_IOSTREAM {
                s: INCHI_IOS_STRING {
                    pStr: heap.allocate_model_storage(storage).unwrap(),
                    nAllocatedLength: allocated,
                    nUsedLength: text.len() as i32,
                    nPtr: 37,
                },
                type_: crate::source_types::INCHI_IOS_TYPE_STRING as i32,
                ..INCHI_IOSTREAM::default()
            }
        }

        fn bytes(heap: &SourceHeap, output: &INCHI_IOSTREAM) -> Vec<u8> {
            heap.slice(output.s.pStr.as_const()).unwrap()[..output.s.nUsedLength as usize]
                .iter()
                .map(|byte| *byte as u8)
                .collect()
        }

        let mut heap = SourceHeap::default();
        let mut null_output = INCHI_IOSTREAM::default();
        assert_eq!(
            EditINCHI_HidePolymerZz(&mut heap, &mut null_output, 0, 1),
            Ok(())
        );
        assert_eq!(
            EditINCHI_HidePolymerZz(&mut heap, &mut null_output, 1, 1),
            Ok(())
        );
        assert_eq!(
            EditINCHI_HidePolymerZz(&mut heap, &mut null_output, 2, 0),
            Err(SourceHeapError::NullPointer)
        );

        for text in [b"InChI=1B/CH4/h1H4".as_slice(), b"prefix/z1-2".as_slice()] {
            let mut output = stream(&mut heap, text, 8);
            let before = output.clone();
            let allocations = heap.live_allocation_count();
            heap.trace_source_allocations();
            assert_eq!(
                EditINCHI_HidePolymerZz(&mut heap, &mut output, 2, 0),
                Ok(())
            );
            assert_eq!(output, before);
            assert_eq!(heap.source_allocation_calls(), 0);
            assert_eq!(heap.live_allocation_count(), allocations);
        }

        const SOURCE_EXAMPLE: &[u8] = b"InChI=1B/C4H4N4.2Zz/c1-5-2-7-4-8-3-6-1;;/h1-4H;;/z101-1-8(9,10-8,3,1,6,2,5,2,7,3,6,1,5,4,7,4,8)/b5-1-,5-2+,6-1+,6-3-,7-2+,7-4+,8-3+,8-4+;;";
        const HIDDEN_EXAMPLE: &[u8] = b"InChI=1B/C4H4N4/c1-5-2-7-4-8-3-6-1/h1-4H/z101-1-8(8,3,1,6,2,5,2,7,3,6,1,5,4,7,4,8)/b5-1-,5-2+,6-1+,6-3-,7-2+,7-4+,8-3+,8-4+";
        for newline in [false, true] {
            let mut input = SOURCE_EXAMPLE.to_vec();
            if newline {
                input.push(b'\n');
            }
            let mut output = stream(&mut heap, &input, 32);
            let pointer = output.s.pStr;
            heap.trace_source_allocations();
            assert_eq!(
                EditINCHI_HidePolymerZz(&mut heap, &mut output, 2, 0),
                Ok(())
            );
            let mut expected = HIDDEN_EXAMPLE.to_vec();
            if newline {
                expected.push(b'\n');
            }
            assert_eq!(bytes(&heap, &output), expected);
            assert_eq!(output.s.nUsedLength, expected.len() as i32);
            assert_eq!(output.s.nPtr, 37);
            assert_eq!(output.s.pStr, pointer);
            assert_eq!(heap.source_allocation_calls(), 1);
            let storage = heap.slice(output.s.pStr.as_const()).unwrap();
            assert_eq!(storage[expected.len()], 0);
            assert!(storage[expected.len() + 1..].contains(&0x55_i8));
        }

        let retained_pattern = b"InChI=1B/C.Zz/c1;;/z1(2-3,4)/h1H;;";
        let retained_expected = b"InChI=1B/C/c1;/z1(2-3,4)/h1H;";
        let mut retained = stream(&mut heap, retained_pattern, 16);
        assert_eq!(
            EditINCHI_HidePolymerZz(&mut heap, &mut retained, 1, 0),
            Ok(())
        );
        assert_eq!(bytes(&heap, &retained), retained_expected);

        let mut allocation_failure = stream(&mut heap, SOURCE_EXAMPLE, 8);
        let before = allocation_failure.clone();
        let before_bytes = heap
            .slice(allocation_failure.s.pStr.as_const())
            .unwrap()
            .to_vec();
        let allocations = heap.live_allocation_count();
        heap.fail_after_allocations(0);
        assert_eq!(
            EditINCHI_HidePolymerZz(&mut heap, &mut allocation_failure, 2, 0),
            Ok(())
        );
        assert_eq!(allocation_failure, before);
        assert_eq!(
            heap.slice(allocation_failure.s.pStr.as_const()).unwrap(),
            before_bytes
        );
        assert_eq!(heap.live_allocation_count(), allocations);
        assert_eq!(heap.source_allocation_calls(), 1);

        let mut print_failure_input = SOURCE_EXAMPLE.to_vec();
        print_failure_input.push(b'\n');
        let mut print_failure = stream(&mut heap, &print_failure_input, 0);
        print_failure.s.nAllocatedLength = 0;
        let original_pointer = print_failure.s.pStr;
        let allocations = heap.live_allocation_count();
        heap.fail_after_allocations(1);
        assert_eq!(
            EditINCHI_HidePolymerZz(&mut heap, &mut print_failure, 2, 0),
            Ok(())
        );
        assert_eq!(print_failure.s.pStr, original_pointer);
        assert_eq!(print_failure.s.nUsedLength, 0);
        assert_eq!(print_failure.s.nAllocatedLength, 0);
        assert_eq!(print_failure.s.nPtr, 37);
        assert_eq!(
            heap.slice(original_pointer.as_const()).unwrap()[SOURCE_EXAMPLE.len()],
            0
        );
        assert_eq!(heap.live_allocation_count(), allocations);
        assert_eq!(heap.source_allocation_calls(), 2);
    }

    #[test]
    fn source_port__ichiprt1__inchi_sort_int_pair_ascending__line_5525() {
        for (first, second, expected) in [
            (1, 2, (1, 2)),
            (2, 1, (1, 2)),
            (0, 0, (0, 0)),
            (-1, -2, (-2, -1)),
            (i32::MIN, i32::MAX, (i32::MIN, i32::MAX)),
            (i32::MAX, i32::MIN, (i32::MIN, i32::MAX)),
        ] {
            let mut a = first;
            let mut b = second;
            inchi_sort_int_pair_ascending(&mut a, &mut b);
            assert_eq!((a, b), expected);
        }
    }

    #[test]
    fn source_port__ichiprt1__getsaveoptletters__line_3100() {
        for bits in u8::MIN..=u8::MAX {
            let mut first = 0_i8;
            let mut second = 0_i8;
            GetSaveOptLetters(bits, &mut first, &mut second);
            assert_eq!(first, (b'A' + (bits & 0x0f)) as i8, "bits={bits:#04x}");
            assert_eq!(
                second,
                (b'A' + ((bits & 0x30) >> 4)) as i8,
                "bits={bits:#04x}"
            );
        }
    }

    #[test]
    fn source_port__ichiprt1__set_line_separators__line_3112() {
        let sentinel: &'static [i8] = &[99, 0];
        for options in [
            0,
            INCHI_OUT_PLAIN_TEXT_COMMENTS as i32,
            -1,
            i32::MIN,
            i32::MAX,
        ] {
            let mut line_feed = sentinel;
            let mut tab = sentinel;
            set_line_separators(options, &mut line_feed, &mut tab);
            let expected_line_feed: &[i8] = if options & INCHI_OUT_PLAIN_TEXT_COMMENTS as i32 != 0 {
                &[b'\n' as i8, 0]
            } else {
                &[0]
            };
            assert_eq!(line_feed, expected_line_feed);
            assert_eq!(tab, &[b'\n' as i8, 0]);
        }
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_versionandkind__line_3137() {
        assert_eq!(
            version_case(0, 0, true),
            (
                b"\nInChI=\n1\n".iter().map(|byte| *byte as i8).collect(),
                vec![b'1' as i8]
            )
        );
        assert_eq!(
            version_case(INCHI_OUT_STDINCHI as i32, 0, true),
            (
                b"\nInChI=\n1S\n".iter().map(|byte| *byte as i8).collect(),
                vec![b'1' as i8, b'S' as i8],
            )
        );
        assert_eq!(
            version_case(INCHI_OUT_STDINCHI as i32, -1, false),
            (
                b"InChI=1B".iter().map(|byte| *byte as i8).collect(),
                vec![b'1' as i8, b'B' as i8],
            )
        );
    }

    #[test]
    fn source_port__ichiprt1__outputinchi_polymerlayer_singleunit__line_4111() {
        use crate::source_types::{
            CLOSING_SRU_DIRADICAL, CLOSING_SRU_HIGHER_ORDER_BOND, CLOSING_SRU_RING, OAD_AtProps,
            OAD_Polymer, OAD_PolymerUnit, ORIG_ATOM_DATA, ORIG_STRUCT, POLYMERS_LEGACY,
            POLYMERS_MODERN, inp_ATOM,
        };

        fn buffer_text(heap: &SourceHeap, buffer: &INCHI_IOS_STRING) -> String {
            let bytes = heap.slice(buffer.pStr.as_const()).unwrap();
            String::from_utf8(
                bytes[..buffer.nUsedLength as usize]
                    .iter()
                    .map(|byte| *byte as u8)
                    .collect(),
            )
            .unwrap()
        }

        fn context(
            heap: &mut SourceHeap,
            treat: i32,
            frame_shift_scheme: i32,
            atom_text: &[u8],
            num_atoms: i32,
        ) -> (ORIG_ATOM_DATA, ORIG_STRUCT) {
            let polymer = heap
                .allocate(vec![OAD_Polymer {
                    treat,
                    frame_shift_scheme,
                    ..OAD_Polymer::default()
                }])
                .unwrap();
            let atoms = heap.allocate(vec![inp_ATOM::default(); 12]).unwrap();
            let atom_text = heap
                .allocate_model_storage(atom_text.iter().map(|byte| *byte as i8).collect())
                .unwrap();
            (
                ORIG_ATOM_DATA {
                    polymer,
                    at: atoms,
                    ..ORIG_ATOM_DATA::default()
                },
                ORIG_STRUCT {
                    num_atoms,
                    szAtoms: atom_text,
                    ..ORIG_STRUCT::default()
                },
            )
        }

        fn properties() -> Vec<OAD_AtProps> {
            let mut values = vec![
                OAD_AtProps {
                    erank: 2,
                    ..OAD_AtProps::default()
                };
                12
            ];
            values[4].erank = 8;
            values
        }

        let mut heap = SourceHeap::default();
        let (original_data, original_structure) =
            context(&mut heap, POLYMERS_LEGACY as i32, 1, b"C\0", 8);
        let alist = heap.allocate(vec![1_i32, 2, 3]).unwrap();
        let mut unsupported = OAD_PolymerUnit {
            type_: 1,
            subtype: 2,
            conn: 3,
            na: 3,
            nb: 3,
            alist,
            ..OAD_PolymerUnit::default()
        };
        let mut used = 0;
        let mut buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut unsupported,
                POLYMERS_LEGACY as i32,
                0,
                &mut used,
                &properties(),
                &[],
                &original_data,
                &original_structure,
                &mut buffer,
            ),
            Ok(12)
        );
        assert_eq!(buffer_text(&heap, &buffer), "123-1-3");

        let crossing = heap.allocate(vec![1_i32, 4, 3, 5]).unwrap();
        let mut crossing_unit = OAD_PolymerUnit {
            type_: 1,
            subtype: 2,
            conn: 3,
            na: 3,
            nb: 2,
            alist,
            blist: crossing,
            ..OAD_PolymerUnit::default()
        };
        let mut legacy_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut crossing_unit,
                POLYMERS_LEGACY as i32,
                0,
                &mut used,
                &properties(),
                &[],
                &original_data,
                &original_structure,
                &mut legacy_buffer,
            ),
            Ok(0)
        );
        assert_eq!(buffer_text(&heap, &legacy_buffer), "123-1-3(4-1,5-3)");
        let mut modern_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut crossing_unit,
                POLYMERS_MODERN as i32,
                0,
                &mut used,
                &properties(),
                &[],
                &original_data,
                &original_structure,
                &mut modern_buffer,
            ),
            Ok(0)
        );
        assert_eq!(buffer_text(&heap, &modern_buffer), "123-1-3(5-3,4-1)");

        let (star_data, star_structure) =
            context(&mut heap, POLYMERS_LEGACY as i32, 1, b"CHHO\0", 10);
        let one_atom = heap.allocate(vec![1_i32]).unwrap();
        let mut diradical = OAD_PolymerUnit {
            na: 1,
            alist: one_atom,
            cap1: 1,
            cap2: 1,
            cyclizable: CLOSING_SRU_DIRADICAL as i32,
            end_atom1: 7,
            ..OAD_PolymerUnit::default()
        };
        used = 0;
        let mut diradical_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut diradical,
                POLYMERS_LEGACY as i32,
                2,
                &mut used,
                &properties(),
                &[],
                &star_data,
                &star_structure,
                &mut diradical_buffer,
            ),
            Ok(0)
        );
        assert_eq!(used, 2);
        assert_eq!(buffer_text(&heap, &diradical_buffer), "000-1(7,8-7)");

        used = 10;
        let mut error_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut diradical,
                POLYMERS_LEGACY as i32,
                2,
                &mut used,
                &properties(),
                &[],
                &star_data,
                &star_structure,
                &mut error_buffer,
            ),
            Ok(11)
        );
        assert_eq!(used, 10);
        assert_eq!(buffer_text(&heap, &error_buffer), "000-1");

        let mut higher_order = OAD_PolymerUnit {
            na: 1,
            alist: one_atom,
            cyclizable: CLOSING_SRU_HIGHER_ORDER_BOND as i32,
            end_atom1: 9,
            end_atom2: 3,
            ..OAD_PolymerUnit::default()
        };
        let mut higher_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut higher_order,
                POLYMERS_LEGACY as i32,
                0,
                &mut used,
                &properties(),
                &[],
                &original_data,
                &original_structure,
                &mut higher_buffer,
            ),
            Ok(0)
        );
        assert_eq!(buffer_text(&heap, &higher_buffer), "000-1(0,0-3.9)");

        let mut ring_without_bonds = OAD_PolymerUnit {
            na: 1,
            alist: one_atom,
            cyclizable: CLOSING_SRU_RING as i32,
            end_atom1: 9,
            end_atom2: 3,
            ..OAD_PolymerUnit::default()
        };
        let mut ring_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut ring_without_bonds,
                POLYMERS_LEGACY as i32,
                0,
                &mut used,
                &properties(),
                &[],
                &original_data,
                &original_structure,
                &mut ring_buffer,
            ),
            Ok(0)
        );
        assert_eq!(buffer_text(&heap, &ring_buffer), "000-1(0,0-3,9)");

        let row0 = heap.allocate(vec![5_i32, 4]).unwrap();
        let row1 = heap.allocate(vec![3_i32, 8]).unwrap();
        let row2 = heap.allocate(vec![2_i32, 7]).unwrap();
        let rows = heap.allocate(vec![row0, row1, row2]).unwrap();
        let mut ring_bonds = OAD_PolymerUnit {
            na: 1,
            alist: one_atom,
            cyclizable: CLOSING_SRU_RING as i32,
            nbkbonds: 3,
            bkbonds: rows,
            ..OAD_PolymerUnit::default()
        };
        let mut bonds_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut ring_bonds,
                POLYMERS_LEGACY as i32,
                0,
                &mut used,
                &properties(),
                &[],
                &original_data,
                &original_structure,
                &mut bonds_buffer,
            ),
            Ok(0)
        );
        assert_eq!(buffer_text(&heap, &bonds_buffer), "000-1(0,0-2,7,3,8,5,4)");

        let (modern_data, modern_structure) =
            context(&mut heap, POLYMERS_MODERN as i32, 0, b"C\0", 8);
        let modern_row0 = heap.allocate(vec![1_i32, 2]).unwrap();
        let modern_row1 = heap.allocate(vec![3_i32, 5]).unwrap();
        let modern_rows = heap.allocate(vec![modern_row0, modern_row1]).unwrap();
        let mut modern_ring = OAD_PolymerUnit {
            na: 1,
            alist: one_atom,
            cap1: 1,
            cap2: 1,
            cyclizable: CLOSING_SRU_RING as i32,
            nbkbonds: 2,
            bkbonds: modern_rows,
            ..OAD_PolymerUnit::default()
        };
        used = 0;
        let mut senior_buffer = merge_hill_buffer(&mut heap, "");
        assert_eq!(
            OutputINCHI_PolymerLayer_SingleUnit(
                &mut heap,
                &mut modern_ring,
                POLYMERS_MODERN as i32,
                2,
                &mut used,
                &properties(),
                &[],
                &modern_data,
                &modern_structure,
                &mut senior_buffer,
            ),
            Ok(0)
        );
        assert_eq!(buffer_text(&heap, &senior_buffer), "000-1(7,8-5,3,1,2)");
        assert_eq!(&heap.slice(modern_row0.as_const()).unwrap()[..2], &[5, 3]);
        assert_eq!(&heap.slice(modern_row1.as_const()).unwrap()[..2], &[1, 2]);
        assert_eq!((modern_ring.end_atom1, modern_ring.end_atom2), (5, 3));
    }
}
use crate::source::base::ichi_bns::nBondsValenceInpAt;
use crate::source::base::ichi_io::{
    inchi_ios_print, inchi_ios_print_nodisplay, inchi_strbuf_printf, inchi_strbuf_reset,
    inchi_strbuf_update,
};
use crate::source::base::ichierr::AddErrorMessage;
use crate::source::base::ichiparm::source_strtod_with_end;
use crate::source::base::ichiprt2::{MakeDelim, MakeIsoHString, inchi_strtol};
use crate::source::base::ichiprt3::{
    bin_AuxTautTrans, str_AuxEqu, str_AuxNumb, str_AuxTautTrans, str_FixedH_atoms, str_HillFormula,
    str_HillFormula2, str_IsoAtoms, str_IsoSp2, str_IsoSp3, str_IsoStereoAbsInv, str_Sp2, str_Sp3,
    str_StereoAbsInv,
};
use crate::source::base::ichisort::{insertions_sort, insertions_sort_AT_RANK};
use crate::source::base::ichister::get_opposite_sb_atom_slice;
use crate::source::base::runichi3::{
    OAD_Polymer_PrepareWorkingSet, OAD_Polymer_SetAtProps, OAD_PolymerUnit_CreateCopy,
    OAD_PolymerUnit_Free,
};
use crate::source::base::util::{
    get_atomic_mass_from_elnum, inchi_calloc, inchi_free, inchi_malloc, is_el_a_metal,
    is_in_the_list, lrtrim, needed_unusual_el_valence,
};
use crate::source_types::{
    AB_MAX_ILL_DEFINED_PARITY, AB_MAX_WELL_DEFINED_PARITY, AB_MIN_ILL_DEFINED_PARITY,
    AB_MIN_WELL_DEFINED_PARITY, AB_PARITY_EVEN, AB_PARITY_IISO, AB_PARITY_ODD, AB_PARITY_UNDF,
    AB_PARITY_UNKN, AT_RANK, BOND_TYPE_ALTERN, BOND_TYPE_DOUBLE, BOND_TYPE_SINGLE,
    BOND_TYPE_TRIPLE, CANON_GLOBALS, FILE, FIX_EITHER_STEREO_IN_AUX_INFO, FLAG_INP_AT_CHIRAL,
    FLAG_INP_AT_NONCHIRAL, INCHI_IOS_STRING, INCHI_IOSTREAM, INCHI_OUT_CTL,
    INCHI_OUT_PLAIN_TEXT_COMMENTS, INCHI_OUT_STDINCHI, INCHI_REC, LEN_COORD,
    MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BONDS, MAXVAL, MOL_COORD, NUM_COORD, ORIG_ATOM_DATA,
    ORIG_STRUCT, RADICAL_DOUBLET, SB_PARITY_MASK, SB_PARITY_SHFT, STEREO_DBLE_EITHER,
    STEREO_SNGL_DOWN, STEREO_SNGL_EITHER, STEREO_SNGL_UP, STRUCT_DATA, SourceConstPointer,
    SourceFormatArgument, SourceHeap, SourceHeapError, SourceMutPointer, SourceVaList, inp_ATOM,
    local_ichiprt1::{INCHI_TAG, ORIG_STR_BUFLEN},
};
