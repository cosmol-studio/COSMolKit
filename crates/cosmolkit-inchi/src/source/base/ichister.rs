use crate::source::base::ichiring::{
    QueueCreate, QueueDelete, is_atom_in_3memb_ring, is_bond_in_Nmax_memb_ring,
};
use crate::source::base::ichisort::{inchi_swap, insertions_sort};
use crate::source::base::util::{get_endpoint_valence, inchi_calloc, inchi_free, is_el_a_metal};
use crate::source_types::local_ichister::{
    MAX_EDGE_RATIO, MAX_SINE, MIN_ANGLE, MIN_ANGLE_DBOND, MIN_ANGLE_RELAXED, MIN_BOND_LEN,
    MIN_LEN_STRAIGHT, MIN_SINE, MIN_SINE_EDGE, MIN_SINE_OUTSIDE, MIN_SINE_RELAXED, MIN_SINE_SQUARE,
    T2D_OKAY, T2D_UNDF, T2D_WARN, ZERO_ANGLE, ZERO_LENGTH, ZTYPE_3D, ZTYPE_DOWN, ZTYPE_EITHER,
    ZTYPE_NONE, ZTYPE_UP,
};
use crate::source_types::{
    AB_MAX_WELL_DEFINED_PARITY, AB_MIN_WELL_DEFINED_PARITY, AB_PARITY_0D, AB_PARITY_CALC,
    AB_PARITY_EVEN, AB_PARITY_IISO, AB_PARITY_NONE, AB_PARITY_ODD, AB_PARITY_UNDF, AB_PARITY_UNKN,
    AMBIGUOUS_STEREO, AMBIGUOUS_STEREO_ERROR, AT_FLAG_ISO_H_POINT, AT_NUMB, AT_RANK, BITS_PARITY,
    BOND_ALT12NS, BOND_ALTERN, BOND_DOUBLE, BOND_MARK_ALL, BOND_SINGLE, BOND_TAUTOM,
    BOND_TYPE_DOUBLE, CANON_GLOBALS, CMODE_NO_ALT_SBONDS, CT_CALC_STEREO_ERR, CT_ERR_MAX,
    CT_ERR_MIN, CT_ISO_H_ERR, CT_OUT_OF_RAM, CT_STEREOBOND_ERROR, FlagSB_0D, FlagSC_0D, INCHI_MODE,
    MAX_CUMULENE_LEN, MAX_NUM_STEREO_ATOM_NEIGH, MAX_NUM_STEREO_BOND_NEIGH, MAX_NUM_STEREO_BONDS,
    MIN_DOT_PROD, MIN_NUM_STEREO_BOND_NEIGH, MULT_STEREOBOND, NUM_H_ISOTOPES,
    PES_BIT_ARSINE_STEREO, PES_BIT_FIX_SP3_BUG, PES_BIT_PHOSPHINE_STEREO,
    PES_BIT_POINT_EDGE_STEREO, QUEUE, RADICAL_SINGLET, REQ_MODE_MIN_SB_RING_MASK,
    REQ_MODE_MIN_SB_RING_SHFT, S_CHAR, SB_PARITY_FLAG, SB_PARITY_MASK, SB_PARITY_SHFT,
    STEREO_DBLE_EITHER, STEREO_SNGL_DOWN, STEREO_SNGL_EITHER, STEREO_SNGL_UP, SourceConstPointer,
    SourceHeap, SourceHeapError, SourceMutPointer, inp_ATOM, qInt, sp_ATOM,
};

#[allow(non_snake_case)]
pub(crate) fn comp_AT_NUMB(a1: AT_NUMB, a2: AT_NUMB) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:140 comp_AT_NUMB
    // INCHI✔️✔️: complete active source frame follows verbatim; typed value access has equivalent constant cost.
    /*
    int comp_AT_NUMB( const void* a1, const void* a2, void *p )
    {
        return (int)(*(const AT_NUMB*)a1) - (int)(*(const AT_NUMB*)a2);
    }
    */
    // END INCHI C FUNCTION: comp_AT_NUMB

    i32::from(a1) - i32::from(a2)
}

#[allow(non_snake_case)]
pub(crate) fn CompDble(a1: i32, a2: i32, pDoubleForSort: &[f64]) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:379 CompDble
    // INCHI✔️❌: complete active source frame follows verbatim; checked slice access adds overhead.
    /*
    int CompDble( const void *a1, const void *a2, void *p )
    {
        double *pDoubleForSort = (double *) p;
        double diff = pDoubleForSort[*(const int*) a1] - pDoubleForSort[*(const int*) a2];
        if (diff > 0.0)
        {
            return 1;
        }
        if (diff < 0.0)
        {
            return -1;
        }
        return 0;
    }
    */
    // END INCHI C FUNCTION: CompDble

    let first_index = usize::try_from(a1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let second_index = usize::try_from(a2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let difference = pDoubleForSort
        .get(first_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        - pDoubleForSort
            .get(second_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if difference > 0.0 {
        return Ok(1);
    }
    if difference < 0.0 {
        return Ok(-1);
    }
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn Get2DTetrahedralAmbiguity(
    _pCG: &mut CANON_GLOBALS,
    at_coord: &[[f64; 3]; MAX_NUM_STEREO_ATOM_NEIGH as usize],
    bAddExplicitNeighbor: i32,
    bFix2DstereoBorderCase: i32,
    vMinAngle: f64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:402 Get2DTetrahedralAmbiguity
    // INCHI✔️❌: complete source frame follows verbatim; checked sort adapters add overhead.
    /*
    int Get2DTetrahedralAmbiguity( CANON_GLOBALS *pCG,
                                   double at_coord[][3],
                                   int bAddExplicitNeighbor,
                                   int bFix2DstereoBorderCase,
                                   double vMinAngle )
    {
        /*    const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
        const double one_pi = 3.14159265358979323846; /* M_PI */
        const double two_pi = 2.0*one_pi;
        const double dAngleAndPiMaxDiff = 2.0*atan2( 1.0, sqrt( 7.0 ) ); /*  min sine between 2 InPlane bonds */
        double *pDoubleForSort;
        int    nBondType[MAX_NUM_STEREO_ATOM_NEIGH], nBondOrder[MAX_NUM_STEREO_ATOM_NEIGH];
        double dBondDirection[MAX_NUM_STEREO_ATOM_NEIGH];
        volatile double dAngle, dAlpha, dLimit, dBisector;
        /* 2010-02-10  added 'volatile': workaround ensuring proper behavior for gcc 32-bit */
        /* cml-enabled compiles at >=O1 for SID484922 and alike (both lin&win had problems) */
        int  nNumNeigh = MAX_NUM_STEREO_ATOM_NEIGH - ( bAddExplicitNeighbor != 0 );
        int  i, num_Up, num_Dn, bPrev_Up, cur_len_Up, cur_first_Up, len_Up, first_Up = 0; /* djb-rwth: initialisation required to avoid garbage values */
        int  ret = 0;

        for (i = 0, num_Up = num_Dn = 0; i < nNumNeigh; i++)
        {
            dAngle = atan2( at_coord[i][1], at_coord[i][0] ); /*  range from -pi to +pi */
            if (dAngle < 0.0)
            {
                dAngle += two_pi;
            }
            dBondDirection[i] = dAngle;
            nBondType[i] = ( at_coord[i][2] > 0.0 ) ? 1 : ( at_coord[i][2] < 0.0 ) ? -1 : 0; /* z-coord sign */
            if (nBondType[i] > 0)
            {
                num_Up++;
            }
            else
            {
                if (nBondType[i] < 0)
                {
                    num_Dn++;
                }
            }
            nBondOrder[i] = i;
        }
        if (num_Up < num_Dn)
        {
            for (i = 0; i < nNumNeigh; i++)
            {
                nBondType[i] = -nBondType[i];
            }
            inchi_swap( (char*) &num_Dn, (char*) &num_Up, sizeof( num_Dn ) );
        }
        if (!num_Up)
        {
            return T2D_UNDF;
        }

        /*  Sort according to the bond orientations */
        pDoubleForSort = dBondDirection;
        insertions_sort( pDoubleForSort, nBondOrder, nNumNeigh,
                         sizeof( nBondOrder[0] ), CompDble );

        /*  Find the longest contiguous sequence of Up bonds */
        if (num_Up == nNumNeigh)
        {
            /*  all bonds are Up */
            len_Up = cur_len_Up = nNumNeigh; /* added cur_len_Up initialization 1/8/2002 */
            first_Up = 0;
        }
        else
        {
            /*  at least one bond is not Up */
            cur_len_Up = len_Up = bPrev_Up = 0;
            /* prev. cycle header version ---
            for ( i = 0; 1; i ++ )
            {
                if ( i >= nNumNeigh && !bPrev_Up )
                {
                    break;
                }
            ----------
            }
            */

            /* look at all bonds and continue (circle therough the beginning) as long as the current bond is Up */
            for (i = 0; i < nNumNeigh || bPrev_Up; i++)
            {
                if (nBondType[nBondOrder[i % nNumNeigh]] > 0)
                {
                    if (bPrev_Up)
                    {
                        cur_len_Up++; /* uncrement number of Up bonds in current contiguous sequence of them */
                    }
                    else
                    {
                        bPrev_Up = 1; /* start new contiguous sequence of Up bonds */
                        cur_len_Up = 1;
                        cur_first_Up = i % nNumNeigh;
                    }
                }
                else
                {
                    if (bPrev_Up)
                    { /* end of contiguous sequence of Up bonds */
                        if (cur_len_Up > len_Up)
                        {
                            first_Up = cur_first_Up; /* store the sequence because it is longer than the ptrvious one */
                            len_Up = cur_len_Up;
                        }
                        bPrev_Up = 0;
                    }
                }
            }
        }
    #if ( FIX_2D_STEREO_BORDER_CASE == 1 )
        /* check if the bonds with ordering numbers first_Up+len_Up and first_Up+len_Up+1 */
        /* have identical angles. In this case switch their order to enlarge the Up sequence */
    #define ZERO_ANGLE  0.000001
        if (nNumNeigh - len_Up >= 2)
        {
            int next1, next2;
            for (i = 1; i < nNumNeigh - len_Up; i++)
            {
                next2 = ( first_Up + len_Up + i ) % nNumNeigh; /* the 2nd after Up sequence */
                if (nBondType[nBondOrder[next2]] > 0)
                {
                    next1 = ( first_Up + len_Up ) % nNumNeigh; /* the 1st after Up sequence */
                    dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                    if (fabs( dAngle ) < ZERO_ANGLE)
                    {
                        inchi_swap( (char*) &nBondOrder[next1], (char*) &nBondOrder[next2], sizeof( nBondOrder[0] ) );
                        len_Up++;
                        break;
                    }
                }
            }
        }
        /* Check whether the not-Up bond (located before the found first-Up) has */
        /* same angle as the Up bond that precedes this not-Up bond */
        if (nNumNeigh - len_Up >= 2)
        {
            int next1, next2;
            for (i = 1; i < nNumNeigh - len_Up; i++)
            {
                next2 = ( first_Up + nNumNeigh - i - 1 ) % nNumNeigh; /* the 2nd before Up sequence */
                if (nBondType[nBondOrder[next2]] > 0)
                {
                    next1 = ( first_Up + nNumNeigh - 1 ) % nNumNeigh; /* the 1st before Up sequence */
                    dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                    if (fabs( dAngle ) < ZERO_ANGLE)
                    {
                        inchi_swap( (char*) &nBondOrder[next1], (char*) &nBondOrder[next2], sizeof( nBondOrder[0] ) );
                        first_Up = next1;
                        len_Up++;
                        break;
                    }
                }
            }
        }
    #else
        if (bFix2DstereoBorderCase)
        {
            /* Check if the bonds with ordering numbers first_Up+len_Up and first_Up+len_Up+1 */
            /* have identical angles. In this case switch their order to enlarge the Up sequence */
    #define ZERO_ANGLE  0.000001
            if (nNumNeigh - len_Up >= 2)
            {
                int next1, next2;
                for (i = 1; i < nNumNeigh - len_Up; i++)
                {
                    next2 = ( first_Up + len_Up + i ) % nNumNeigh; /* the 2nd after Up sequence */
                    if (nBondType[nBondOrder[next2]] > 0)
                    {
                        next1 = ( first_Up + len_Up ) % nNumNeigh; /* the 1st after Up sequence */
                        dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                        if (fabs( dAngle ) < ZERO_ANGLE)
                        {
                            inchi_swap( (char*) &nBondOrder[next1], (char*) &nBondOrder[next2], sizeof( nBondOrder[0] ) );
                            len_Up++;
                            break;
                        }
                    }
                }
            }

            /* Check whether the not-Up bond (located before the found first-Up) has */
            /* same angle as the Up bond that precedes this not-Up bond */
            if (nNumNeigh - len_Up >= 2)
            {
                int next1, next2;
                for (i = 1; i < nNumNeigh - len_Up; i++)
                {
                    next2 = ( first_Up + nNumNeigh - i - 1 ) % nNumNeigh; /* the 2nd before Up sequence */
                    if (nBondType[nBondOrder[next2]] > 0)
                    {
                        next1 = ( first_Up + nNumNeigh - 1 ) % nNumNeigh; /* the 1st before Up sequence */
                        dAngle = dBondDirection[nBondOrder[next1]] - dBondDirection[nBondOrder[next2]];
                        if (fabs( dAngle ) < ZERO_ANGLE)
                        {
                            inchi_swap( (char*) &nBondOrder[next1], (char*) &nBondOrder[next2], sizeof( nBondOrder[0] ) );
                            first_Up = next1;
                            len_Up++;
                            break;
                        }
                    }
                }
            }
        }
    #endif
        /*  Turn all the bonds around the center so that */
        /*  the 1st Up bond has zero radian direction */
        dAlpha = dBondDirection[nBondOrder[first_Up]];
        for (i = 0; i < nNumNeigh; i++)
        {
            if (i == nBondOrder[first_Up])
            {
                dBondDirection[i] = 0.0;
            }
            else
            {
                dAngle = dBondDirection[i] - dAlpha;
                if (dAngle < 0.0)
                {
                    dAngle += two_pi;
                }
                dBondDirection[i] = dAngle;
            }
        }

        /********************************************************
         * Process particular cases
         ********************************************************/

        if (nNumNeigh == 3)
        {
            /************************ 3 bonds ************************/
            switch (num_Up)
            {

                case 0:    /* 0 Up */
                    return T2D_UNDF;

                case 1:    /* 1 Up */
                    if (num_Dn)
                    {
    #ifdef _DEBUG
                        if (num_Dn != 1)  /*  debug only */
                        {
                            return -1;
                        }
    #endif
                        ret = ( T2D_UNDF | T2D_WARN );
                    }
                    else
                    {
                        dAngle = dBondDirection[nBondOrder[( first_Up + 2 ) % nNumNeigh]] -
                            dBondDirection[nBondOrder[( first_Up + 1 ) % nNumNeigh]];

                        if (dAngle < 0.0)
                        {
                            dAngle += two_pi;
                        }
                        if (dAngle - one_pi < -vMinAngle || dAngle - one_pi > vMinAngle)
                        {
                            ret = T2D_OKAY;
                        }
                        else
                        {
                            ret = ( T2D_UNDF | T2D_WARN );
                        }
                    }
                    break;

                case 2:    /* 2 Up */
                    if (num_Dn)
                    {
                        dAlpha = dBondDirection[nBondOrder[( first_Up + 1 ) % nNumNeigh]] -
                            dBondDirection[nBondOrder[( first_Up ) % nNumNeigh]];

                        if (dAlpha < 0.0)
                        {
                            dAlpha += two_pi;
                        }

                        if (dAlpha > one_pi - vMinAngle)
                        {
                            ret = T2D_OKAY;
                        }
                        else if (dAlpha < two_pi / 3.0 - vMinAngle)
                        {
                            ret = ( T2D_UNDF | T2D_WARN );
                        }
                        else
                        {
                            /*  angle between 2 Up bonds is between 120 and 180 degrees */
                            /*  direction of the (Alpha angle bisector) + 180 degrees    */
                            dBisector = dBondDirection[nBondOrder[( first_Up ) % nNumNeigh]];
                            dBisector += dBondDirection[nBondOrder[( first_Up + 1 ) % nNumNeigh]];
                            dBisector /= 2.0;
                            dBisector -= one_pi;
                            if (dBisector < 0.0)
                            {
                                dBisector += two_pi;
                            }
                            if (dAlpha < two_pi / 3.0 + vMinAngle)
                            {
                                /*  dAlpha is inside ( 2pi/3 - eps, 2pi/3 + eps ) interval */
                                dLimit = vMinAngle * 3.0 / 2.0;
                            }
                            else
                            {
                                dLimit = dAlpha * 3.0 / 2.0 - one_pi;
                            }

                            dAngle = dBondDirection[nBondOrder[( first_Up + 2 ) % nNumNeigh]];

                            if (dBisector - dAngle < -dLimit ||
                                  dBisector - dAngle >  dLimit)
                            {
                                ret = ( T2D_UNDF | T2D_WARN );
                            }
                            else
                            {
                                ret = T2D_OKAY;
                            }
                        }
                    } /* if ( num_Dn )  */
                    else
                    {
                        ret = T2D_OKAY;
                    }
                    break;

                case 3:    /* 3 Up */
                    ret = T2D_OKAY;
                    break;

                default:/* other Up */
                    return -1;
            } /* eof switch( num_Up ) at  nNumNeigh == 3 */
        }

        else if (nNumNeigh == 4)
        {
            /******************************* 4 bonds ********************/
            switch (num_Up)
            {

                case 0:    /* 0 Up */
                    return T2D_UNDF;

                case 1:    /* 1 Up */
                    if (num_Dn)
                    {
                        if (nBondType[nBondOrder[( first_Up + 2 ) % nNumNeigh]] < 0)
                        {
                            /*
                            * Up, In Plane, Dn, In Plane. Undefined if angle between
                            * two In Plane bonds is wuthin pi +/- 2*arcsine(1/sqrt(8)) interval
                            * That is, 138.5 to 221.4 degrees; for certainty the interval is
                            * increased by 5.7 degrees at each end to
                            * 134.8 to 227.1 degrees
                            */
                            dAngle = dBondDirection[nBondOrder[( first_Up + 3 ) % nNumNeigh]] -
                                dBondDirection[nBondOrder[( first_Up + 1 ) % nNumNeigh]];
                            if (dAngle < 0.0)
                            {
                                dAngle += two_pi;
                            }
                            if (fabs( dAngle - one_pi ) < dAngleAndPiMaxDiff + vMinAngle)
                            {
                                ret = ( T2D_UNDF | T2D_WARN );
                            }
                            else
                            {
                                ret = T2D_OKAY;
                            }
                        }
                        else
                        {
                            ret = T2D_OKAY;
                        }
    #ifdef _DEBUG
                        if (num_Dn != 1)  /*  debug only */
                        {
                            return -1;
                        }
    #endif
                    }
                    else
                    {
                        ret = T2D_OKAY;
                        dAngle = dBondDirection[nBondOrder[( first_Up + 3 ) % nNumNeigh]] -
                            dBondDirection[nBondOrder[( first_Up + 1 ) % nNumNeigh]];
                        if (dAngle < 0.0)
                        {
                            dAngle += two_pi;
                        }
                        if (dAngle < one_pi - vMinAngle)
                        {
                            ret |= T2D_WARN;
                        }
                    }
                    break;

                case 2:    /* 2 Up */
    #if ( FIX_2D_STEREO_BORDER_CASE == 1 )
                    if (len_Up == 1)
                    {
                        ret = T2D_OKAY;
                    }
                    else
                    {
                        dAngle = dBondDirection[nBondOrder[( first_Up + 3 ) % nNumNeigh]] -
                            dBondDirection[nBondOrder[( first_Up + 0 ) % nNumNeigh]];
                        dAngle = fabs( two_pi - dAngle );
                        dAlpha = dBondDirection[nBondOrder[( first_Up + 2 ) % nNumNeigh]] -
                            dBondDirection[nBondOrder[( first_Up + 1 ) % nNumNeigh]];
                        dAlpha = fabs( dAlpha );
                        if (dAngle < 2.0 * ZERO_ANGLE && dAlpha > vMinAngle ||
                             dAlpha < 2.0 * ZERO_ANGLE && dAngle > vMinAngle)
                        {
                            ret = ( T2D_OKAY | T2D_WARN );
                        }
                        else
                        {
                            ret = ( T2D_UNDF | T2D_WARN );
                        }
                    }
    #else
                    if (bFix2DstereoBorderCase)
                    {
                        /* bug fix */
                        if (len_Up == 1)
                        {
                            ret = T2D_OKAY;
                        }
                        else
                        {
                            dAngle = dBondDirection[nBondOrder[( first_Up + 3 ) % nNumNeigh]] -
                                dBondDirection[nBondOrder[( first_Up + 0 ) % nNumNeigh]];
                            dAngle = fabs( two_pi - dAngle );
                            dAlpha = dBondDirection[nBondOrder[( first_Up + 2 ) % nNumNeigh]] -
                                dBondDirection[nBondOrder[( first_Up + 1 ) % nNumNeigh]];
                            dAlpha = fabs( dAlpha );
                            if ((dAngle < 2.0 * ZERO_ANGLE && dAlpha > vMinAngle) ||
                                 (dAlpha < 2.0 * ZERO_ANGLE && dAngle > vMinAngle)) /* djb-rwth: addressing LLVM warnings */
                            {
                                ret = ( T2D_OKAY | T2D_WARN );
                            }
                            else
                            {
                                ret = ( T2D_UNDF | T2D_WARN );
                            }
                        }
                    }
                    else
                    {
                        /* original InChI v. 1 bug */
                        if (cur_len_Up == 1)
                        {
                            ret = T2D_OKAY;
                        }
                        else
                        {
                            ret = ( T2D_UNDF | T2D_WARN );
                        }
                    }
    #endif
                    break;

                case 3:    /* 3 Up */
                    ret = T2D_OKAY;
                    dAngle = dBondDirection[nBondOrder[( first_Up + 2 ) % nNumNeigh]] -
                        dBondDirection[nBondOrder[( first_Up + 0 ) % nNumNeigh]];
                    if (dAngle < 0.0)
                    {
                        dAngle += two_pi;
                    }
                    if (dAngle < one_pi - vMinAngle)
                    {
                        ret |= T2D_WARN;
                    }
                    break;

                case 4:    /* 4 Up */
                    ret = ( T2D_UNDF | T2D_WARN );
                    break;

                default:/* other Up */
                    return -1; /*  program error */
            } /* eof switch( num_Up ) at  nNumNeigh == 4 */


            if (ret == T2D_OKAY)
            {
                /*  Check whether all bonds are inside a less than 180 degrees sector */
                for (i = 0; i < nNumNeigh; i++)
                {
                    dAngle = dBondDirection[nBondOrder[( i + nNumNeigh - 1 ) % nNumNeigh]] -
                        dBondDirection[nBondOrder[i % nNumNeigh]];
                    if (dAngle < 0.0)
                    {
                        dAngle += two_pi;
                    }
                    if (dAngle < one_pi - vMinAngle)
                    {
                        ret |= T2D_WARN;
                        break;
                    }
                }
            }
        } /* eof nNumNeigh == 4 */

        else
        {
            /*************************** number of bonds != 3 or 4 ******************/

            return -1; /*  error */
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: Get2DTetrahedralAmbiguity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: Get2DTetrahedralAmbiguity
    // INCHI✔️❌: FIX_2D_STEREO_BORDER_CASE == 0, so the runtime bFix2DstereoBorderCase branch is active.
    // INCHI✔️❌: _DEBUG is not defined for the selected GCC/Linux TARGET_API_LIB production target.
    // INCHI✔️❌: T2D_OKAY == 1, T2D_WARN == 2, T2D_UNDF == 4, and ZERO_ANGLE == 0.000001.
    // END INCHI ACTIVE MACRO CONFIGURATION: Get2DTetrahedralAmbiguity

    fn swap_i32(values: &mut [i32], first: usize, second: usize) -> Result<(), SourceHeapError> {
        let width = std::mem::size_of::<i32>();
        inchi_swap(
            bytemuck::cast_slice_mut::<i32, u8>(values),
            first
                .checked_mul(width)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?,
            second
                .checked_mul(width)
                .ok_or(SourceHeapError::AllocationSizeOverflow)?,
            width,
        )
    }

    let one_pi = 3.14159265358979323846_f64;
    let two_pi = 2.0 * one_pi;
    let angle_and_pi_max_difference = 2.0 * 1.0_f64.atan2(7.0_f64.sqrt());
    let number_of_neighbors =
        MAX_NUM_STEREO_ATOM_NEIGH as usize - usize::from(bAddExplicitNeighbor != 0);
    let mut bond_type = [0_i32; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut bond_order = [0_i32; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut bond_direction = [0.0_f64; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut number_up = 0_i32;
    let mut number_down = 0_i32;

    for i in 0..number_of_neighbors {
        let mut angle = at_coord[i][1].atan2(at_coord[i][0]);
        if angle < 0.0 {
            angle += two_pi;
        }
        bond_direction[i] = angle;
        bond_type[i] = if at_coord[i][2] > 0.0 {
            1
        } else if at_coord[i][2] < 0.0 {
            -1
        } else {
            0
        };
        if bond_type[i] > 0 {
            number_up += 1;
        } else if bond_type[i] < 0 {
            number_down += 1;
        }
        bond_order[i] = i as i32;
    }
    if number_up < number_down {
        for value in bond_type.iter_mut().take(number_of_neighbors) {
            *value = -*value;
        }
        let mut counts = [number_down, number_up];
        swap_i32(&mut counts, 0, 1)?;
        number_down = counts[0];
        number_up = counts[1];
    }
    if number_up == 0 {
        return Ok(T2D_UNDF as i32);
    }

    insertions_sort(
        bytemuck::cast_slice_mut::<i32, u8>(&mut bond_order),
        number_of_neighbors,
        std::mem::size_of::<i32>(),
        &mut |left, right| {
            CompDble(
                i32::from_ne_bytes([left[0], left[1], left[2], left[3]]),
                i32::from_ne_bytes([right[0], right[1], right[2], right[3]]),
                &bond_direction,
            )
        },
    )?;

    let mut first_up = 0_usize;
    let mut current_length_up: i32;
    let length_up: i32;
    if number_up == number_of_neighbors as i32 {
        current_length_up = number_of_neighbors as i32;
        length_up = current_length_up;
    } else {
        current_length_up = 0;
        let mut longest_length_up = 0_i32;
        let mut previous_up = false;
        let mut current_first_up = 0_usize;
        let mut i = 0_usize;
        while i < number_of_neighbors || previous_up {
            let position = i % number_of_neighbors;
            let original = bond_order[position] as usize;
            if bond_type[original] > 0 {
                if previous_up {
                    current_length_up += 1;
                } else {
                    previous_up = true;
                    current_length_up = 1;
                    current_first_up = position;
                }
            } else if previous_up {
                if current_length_up > longest_length_up {
                    first_up = current_first_up;
                    longest_length_up = current_length_up;
                }
                previous_up = false;
            }
            i += 1;
        }
        length_up = longest_length_up;
    }
    let mut length_up = length_up;

    if bFix2DstereoBorderCase != 0 {
        if number_of_neighbors as i32 - length_up >= 2 {
            for i in 1..(number_of_neighbors as i32 - length_up) {
                let next2 = (first_up + length_up as usize + i as usize) % number_of_neighbors;
                if bond_type[bond_order[next2] as usize] > 0 {
                    let next1 = (first_up + length_up as usize) % number_of_neighbors;
                    let angle = bond_direction[bond_order[next1] as usize]
                        - bond_direction[bond_order[next2] as usize];
                    if angle.abs() < ZERO_ANGLE {
                        swap_i32(&mut bond_order, next1, next2)?;
                        length_up += 1;
                        break;
                    }
                }
            }
        }
        if number_of_neighbors as i32 - length_up >= 2 {
            for i in 1..(number_of_neighbors as i32 - length_up) {
                let next2 = (first_up + number_of_neighbors - i as usize - 1) % number_of_neighbors;
                if bond_type[bond_order[next2] as usize] > 0 {
                    let next1 = (first_up + number_of_neighbors - 1) % number_of_neighbors;
                    let angle = bond_direction[bond_order[next1] as usize]
                        - bond_direction[bond_order[next2] as usize];
                    if angle.abs() < ZERO_ANGLE {
                        swap_i32(&mut bond_order, next1, next2)?;
                        first_up = next1;
                        length_up += 1;
                        break;
                    }
                }
            }
        }
    }

    let alpha = bond_direction[bond_order[first_up] as usize];
    for i in 0..number_of_neighbors {
        if i == bond_order[first_up] as usize {
            bond_direction[i] = 0.0;
        } else {
            let mut angle = bond_direction[i] - alpha;
            if angle < 0.0 {
                angle += two_pi;
            }
            bond_direction[i] = angle;
        }
    }

    let mut result = 0_i32;
    if number_of_neighbors == 3 {
        match number_up {
            0 => return Ok(T2D_UNDF as i32),
            1 => {
                if number_down != 0 {
                    result = (T2D_UNDF | T2D_WARN) as i32;
                } else {
                    let mut angle = bond_direction
                        [bond_order[(first_up + 2) % number_of_neighbors] as usize]
                        - bond_direction[bond_order[(first_up + 1) % number_of_neighbors] as usize];
                    if angle < 0.0 {
                        angle += two_pi;
                    }
                    if angle - one_pi < -vMinAngle || angle - one_pi > vMinAngle {
                        result = T2D_OKAY as i32;
                    } else {
                        result = (T2D_UNDF | T2D_WARN) as i32;
                    }
                }
            }
            2 => {
                if number_down != 0 {
                    let mut alpha = bond_direction
                        [bond_order[(first_up + 1) % number_of_neighbors] as usize]
                        - bond_direction[bond_order[first_up % number_of_neighbors] as usize];
                    if alpha < 0.0 {
                        alpha += two_pi;
                    }
                    if alpha > one_pi - vMinAngle {
                        result = T2D_OKAY as i32;
                    } else if alpha < two_pi / 3.0 - vMinAngle {
                        result = (T2D_UNDF | T2D_WARN) as i32;
                    } else {
                        let mut bisector = bond_direction
                            [bond_order[first_up % number_of_neighbors] as usize]
                            + bond_direction
                                [bond_order[(first_up + 1) % number_of_neighbors] as usize];
                        bisector /= 2.0;
                        bisector -= one_pi;
                        if bisector < 0.0 {
                            bisector += two_pi;
                        }
                        let limit = if alpha < two_pi / 3.0 + vMinAngle {
                            vMinAngle * 3.0 / 2.0
                        } else {
                            alpha * 3.0 / 2.0 - one_pi
                        };
                        let angle = bond_direction
                            [bond_order[(first_up + 2) % number_of_neighbors] as usize];
                        if bisector - angle < -limit || bisector - angle > limit {
                            result = (T2D_UNDF | T2D_WARN) as i32;
                        } else {
                            result = T2D_OKAY as i32;
                        }
                    }
                } else {
                    result = T2D_OKAY as i32;
                }
            }
            3 => result = T2D_OKAY as i32,
            _ => return Ok(-1),
        }
    } else if number_of_neighbors == 4 {
        match number_up {
            0 => return Ok(T2D_UNDF as i32),
            1 => {
                if number_down != 0 {
                    if bond_type[bond_order[(first_up + 2) % number_of_neighbors] as usize] < 0 {
                        let mut angle = bond_direction
                            [bond_order[(first_up + 3) % number_of_neighbors] as usize]
                            - bond_direction
                                [bond_order[(first_up + 1) % number_of_neighbors] as usize];
                        if angle < 0.0 {
                            angle += two_pi;
                        }
                        if (angle - one_pi).abs() < angle_and_pi_max_difference + vMinAngle {
                            result = (T2D_UNDF | T2D_WARN) as i32;
                        } else {
                            result = T2D_OKAY as i32;
                        }
                    } else {
                        result = T2D_OKAY as i32;
                    }
                } else {
                    result = T2D_OKAY as i32;
                    let mut angle = bond_direction
                        [bond_order[(first_up + 3) % number_of_neighbors] as usize]
                        - bond_direction[bond_order[(first_up + 1) % number_of_neighbors] as usize];
                    if angle < 0.0 {
                        angle += two_pi;
                    }
                    if angle < one_pi - vMinAngle {
                        result |= T2D_WARN as i32;
                    }
                }
            }
            2 => {
                if bFix2DstereoBorderCase != 0 {
                    if length_up == 1 {
                        result = T2D_OKAY as i32;
                    } else {
                        let angle = (two_pi
                            - (bond_direction
                                [bond_order[(first_up + 3) % number_of_neighbors] as usize]
                                - bond_direction
                                    [bond_order[first_up % number_of_neighbors] as usize]))
                            .abs();
                        let alpha = (bond_direction
                            [bond_order[(first_up + 2) % number_of_neighbors] as usize]
                            - bond_direction
                                [bond_order[(first_up + 1) % number_of_neighbors] as usize])
                            .abs();
                        if (angle < 2.0 * ZERO_ANGLE && alpha > vMinAngle)
                            || (alpha < 2.0 * ZERO_ANGLE && angle > vMinAngle)
                        {
                            result = (T2D_OKAY | T2D_WARN) as i32;
                        } else {
                            result = (T2D_UNDF | T2D_WARN) as i32;
                        }
                    }
                } else if current_length_up == 1 {
                    result = T2D_OKAY as i32;
                } else {
                    result = (T2D_UNDF | T2D_WARN) as i32;
                }
            }
            3 => {
                result = T2D_OKAY as i32;
                let mut angle = bond_direction
                    [bond_order[(first_up + 2) % number_of_neighbors] as usize]
                    - bond_direction[bond_order[first_up % number_of_neighbors] as usize];
                if angle < 0.0 {
                    angle += two_pi;
                }
                if angle < one_pi - vMinAngle {
                    result |= T2D_WARN as i32;
                }
            }
            4 => result = (T2D_UNDF | T2D_WARN) as i32,
            _ => return Ok(-1),
        }

        if result == T2D_OKAY as i32 {
            for i in 0..number_of_neighbors {
                let mut angle = bond_direction
                    [bond_order[(i + number_of_neighbors - 1) % number_of_neighbors] as usize]
                    - bond_direction[bond_order[i % number_of_neighbors] as usize];
                if angle < 0.0 {
                    angle += two_pi;
                }
                if angle < one_pi - vMinAngle {
                    result |= T2D_WARN as i32;
                    break;
                }
            }
        }
    } else {
        return Ok(-1);
    }

    Ok(result)
}

pub(crate) fn get_z_coord(
    heap: &SourceHeap,
    at: SourceConstPointer<inp_ATOM>,
    cur_atom: i32,
    neigh_no: i32,
    nType: &mut i32,
    bPointedEdgeStereo: i32,
) -> Result<f64, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:147 get_z_coord
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    double get_z_coord( inp_ATOM* at,
                        int cur_atom,
                        int neigh_no,
                        int *nType,
                        int bPointedEdgeStereo )
    {
        int stereo_value = at[cur_atom].bond_stereo[neigh_no];
        int stereo_type = abs( stereo_value );
        int neigh = (int) at[cur_atom].neighbor[neigh_no];
        double z = at[neigh].z - at[cur_atom].z;
        int    bFlat;

        if ((bFlat = ( fabs( z ) < ZERO_LENGTH ))) /* djb-rwth: addressing LLVM warning */
        {
            int i;
            for (i = 0; i < at[cur_atom].valence; i++)
            {
                if (fabs( at[cur_atom].z - at[(int) at[cur_atom].neighbor[i]].z ) > ZERO_LENGTH)
                {
                    bFlat = 0;
                    break;
                }
            }
        }

        if (bFlat)
        {
            if (!bPointedEdgeStereo || bPointedEdgeStereo * stereo_value >= 0)
            {
                /* bPointedEdgeStereo > 0: define stereo from pointed end of the stereo bond only */
                /* bPointedEdgeStereo < 0: define stereo from wide end of the stereo bond only (case of removed H) */
                switch (stereo_type)
                {
                    /*  1=Up (solid triangle), 6=Down (Dashed triangle), 4=Either (zigzag triangle) */
                    case 0: /*  No stereo */
                        *nType = ZTYPE_NONE;
                        break;
                    case STEREO_SNGL_UP: /*  1= Up */
                        *nType = ZTYPE_UP;
                        break;
                    case STEREO_SNGL_EITHER: /*  4 = Either */
                        *nType = ZTYPE_EITHER;
                        break;
                    case STEREO_SNGL_DOWN: /*  6 = Down */
                        *nType = ZTYPE_DOWN;
                        break;
                    default:
                        *nType = ZTYPE_NONE; /*  ignore unexpected values */
                }
                if (stereo_value < 0 && ( *nType == ZTYPE_DOWN || *nType == ZTYPE_UP ))
                    *nType = -*nType;
            }
            else
            {
                *nType = ZTYPE_NONE; /* no stereo */
            }
        }
        else
        {
            if (stereo_type == STEREO_SNGL_EITHER &&
                ( !bPointedEdgeStereo || bPointedEdgeStereo * stereo_value >= 0 ))
            {
                *nType = ZTYPE_EITHER;
            }
            else
            {
                *nType = ZTYPE_3D;
            }
        }

        return z;
    }
    */
    // END INCHI C FUNCTION: get_z_coord

    let atoms = heap.slice(at)?;
    let current_index =
        usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let order = usize::try_from(neigh_no).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = atoms
        .get(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let stereo_value = i32::from(
        *current
            .bond_stereo
            .get(order)
            .ok_or(SourceHeapError::PointerOutOfBounds)?,
    );
    let stereo_type = stereo_value.wrapping_abs();
    let neighbor = atoms
        .get(usize::from(
            *current
                .neighbor
                .get(order)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let z = neighbor.z - current.z;
    let mut flat = z.abs() < ZERO_LENGTH;
    if flat {
        for i in 0..i32::from(current.valence) {
            let neighbor = atoms
                .get(usize::from(current.neighbor[i as usize]))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if (current.z - neighbor.z).abs() > ZERO_LENGTH {
                flat = false;
                break;
            }
        }
    }
    let accepted = bPointedEdgeStereo == 0 || bPointedEdgeStereo.wrapping_mul(stereo_value) >= 0;
    if flat {
        if accepted {
            *nType = match stereo_type as u32 {
                0 => ZTYPE_NONE as i32,
                STEREO_SNGL_UP => ZTYPE_UP as i32,
                STEREO_SNGL_EITHER => ZTYPE_EITHER as i32,
                STEREO_SNGL_DOWN => ZTYPE_DOWN,
                _ => ZTYPE_NONE as i32,
            };
            if stereo_value < 0 && (*nType == ZTYPE_DOWN || *nType == ZTYPE_UP as i32) {
                *nType = -*nType;
            }
        } else {
            *nType = ZTYPE_NONE as i32;
        }
    } else if stereo_type as u32 == STEREO_SNGL_EITHER && accepted {
        *nType = ZTYPE_EITHER as i32;
    } else {
        *nType = ZTYPE_3D as i32;
    }
    Ok(z)
}

pub(crate) fn len3(c: &[f64; 3]) -> f64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:221 len3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size copy and constant-work arithmetic match.
    /*
    double len3( const double c[] ) /* djb-rwth: avoiding uninitialised values */
    {
        double tmpar[ARR_DIM] = { 0.0 };
    #if USE_BCF
        memcpy_s(tmpar, ARR_DIM * sizeof(c[0]), c, ARR_DIM * sizeof(c[0]));
    #else
        memcpy(tmpar, c, ARR_DIM * sizeof(c[0]));
    #endif
        return sqrt(pow(tmpar[0],2.0) + pow(tmpar[1],2.0) + pow(tmpar[2],2.0) );
    }
    */
    // END INCHI C FUNCTION: len3
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: len3
    // INCHI✔️✔️: ARR_DIM == 3.
    // INCHI✔️✔️: USE_BCF == 0 selects memcpy for the GCC/Linux production target.
    // END INCHI ACTIVE MACRO CONFIGURATION: len3

    let tmpar = *c;
    (tmpar[0].powf(2.0) + tmpar[1].powf(2.0) + tmpar[2].powf(2.0)).sqrt()
}

pub(crate) fn len2(c: &[f64; 2]) -> f64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:234 len2
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size copy and constant-work arithmetic match.
    /*
    double len2( const double c[] ) /* djb-rwth: avoiding uninitialised values */
    {
        double tmpar[ARR_DIM - 1] = { 0,0 };
    #if USE_BCF
        memcpy_s(tmpar, (ARR_DIM - 1) * sizeof(c[0]), c, (ARR_DIM - 1) * sizeof(c[0]));
    #else
        memcpy(tmpar, c, (ARR_DIM - 1) * sizeof(c[0]));
    #endif
        return sqrt(pow(tmpar[0],2.0) + pow(tmpar[1],2.0));
    }
    */
    // END INCHI C FUNCTION: len2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: len2
    // INCHI✔️✔️: ARR_DIM == 3, so ARR_DIM - 1 copies exactly two doubles.
    // INCHI✔️✔️: USE_BCF == 0 selects memcpy for the GCC/Linux production target.
    // END INCHI ACTIVE MACRO CONFIGURATION: len2

    let tmpar = *c;
    (tmpar[0].powf(2.0) + tmpar[1].powf(2.0)).sqrt()
}

pub(crate) fn diff3<'a>(a: &[f64; 3], b: &[f64; 3], result: &'a mut [f64; 3]) -> &'a mut [f64; 3] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:247 diff3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size writes and returned destination identity match.
    /*
    void* diff3( const double a[], const double b[], double result[] ) /* djb-rwth: changed function type */
    {

        result[0] = a[0] - b[0];
        result[1] = a[1] - b[1];
        result[2] = a[2] - b[2];

        return result;
    }
    */
    // END INCHI C FUNCTION: diff3

    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
    result
}

pub(crate) fn add3(a: &[f64; 3], b: &[f64; 3], result: &mut [f64; 3]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:259 add3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size writes and constant-work arithmetic match.
    /*
    void add3( const double a[], const double b[], double result[] ) /* djb-rwth: changed function type */
    {
        result[0] = a[0] + b[0];
        result[1] = a[1] + b[1];
        result[2] = a[2] + b[2];

        /* return result; */
    }
    */
    // END INCHI C FUNCTION: add3

    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
}

pub(crate) fn mult3(a: &[f64; 3], b: f64, result: &mut [f64; 3]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:270 mult3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size writes and constant-work arithmetic match.
    /*
    void mult3( const double a[], double b, double result[] ) /* djb-rwth: changed function type */
    {
        result[0] = a[0] * b;
        result[1] = a[1] * b;
        result[2] = a[2] * b;

        /* return result; */
    }
    */
    // END INCHI C FUNCTION: mult3

    result[0] = a[0] * b;
    result[1] = a[1] * b;
    result[2] = a[2] * b;
}

pub(crate) fn change_sign3(a: &[f64; 3], result: &mut [f64; 3]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:292 change_sign3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size writes and allocation behavior match.
    /*
    void change_sign3( const double a[], double result[] ) /* djb-rwth: changed function type */
    {
        result[0] = -a[0];
        result[1] = -a[1];
        result[2] = -a[2];

        /* return result; */
    }
    */
    // END INCHI C FUNCTION: change_sign3

    result[0] = -a[0];
    result[1] = -a[1];
    result[2] = -a[2];
}

pub(crate) fn dot_prod3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:303 dot_prod3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size arithmetic and evaluation order match.
    /*
    double dot_prod3( const double a[], const double b[] )
    {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    */
    // END INCHI C FUNCTION: dot_prod3

    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn dot_prodchar3(a: &[S_CHAR; 3], b: &[S_CHAR; 3]) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:310 dot_prodchar3
    // INCHI✔️✔️: complete active source frame follows verbatim; scalar arithmetic and allocation behavior match.
    /*
    int dot_prodchar3( const S_CHAR a[], const S_CHAR b[] )
    {
        int prod = ( (int)a[0]*(int)b[0] + (int)a[1]*(int)b[1] + (int)a[2]*(int)b[2] ) / 100;
        if (prod > 100)
        {
            prod = 100;
        }
        else
        {
            if (prod < -100)
            {
                prod = -100;
            }
        }

        return prod;
    }
    */
    // END INCHI C FUNCTION: dot_prodchar3

    let mut product = (i32::from(a[0]) * i32::from(b[0])
        + i32::from(a[1]) * i32::from(b[1])
        + i32::from(a[2]) * i32::from(b[2]))
        / 100;
    if product > 100 {
        product = 100;
    } else if product < -100 {
        product = -100;
    }
    product
}

pub(crate) fn cross_prod3<'a>(
    a: &[f64; 3],
    b: &[f64; 3],
    result: &'a mut [f64; 3],
) -> &'a mut [f64; 3] {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:330 cross_prod3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size temporary and write order match.
    /*
    void* cross_prod3( const double a[], const double b[], double result[] ) /* djb-rwth: changed function type */
    {
        double tmp[3];

        tmp[0] = ( a[1] * b[2] - a[2] * b[1] );
        tmp[1] = -( a[0] * b[2] - a[2] * b[0] );
        tmp[2] = ( a[0] * b[1] - a[1] * b[0] );

        result[0] = tmp[0];
        result[1] = tmp[1];
        result[2] = tmp[2];

        return result;
    }
    */
    // END INCHI C FUNCTION: cross_prod3

    let mut tmp = [0.0; 3];
    tmp[0] = a[1] * b[2] - a[2] * b[1];
    tmp[1] = -(a[0] * b[2] - a[2] * b[0]);
    tmp[2] = a[0] * b[1] - a[1] * b[0];
    result[0] = tmp[0];
    result[1] = tmp[1];
    result[2] = tmp[2];
    result
}

pub(crate) fn triple_prod(
    a: &[f64; 3],
    b: &[f64; 3],
    c: &[f64; 3],
    sine_value: Option<&mut f64>,
) -> f64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:347 triple_prod
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size callee composition and branching match.
    /*
    double triple_prod( double a[], double b[], double c[], double *sine_value )
    {
        double ab[3], dot_prod_ab_c, abs_c, abs_ab;
        cross_prod3( a, b, ab );
        /* ab[0] =  (a[1]*b[2]-a[2]*b[1]); */
        /* ab[1] = -(a[0]*b[2]-a[2]*b[0]); */
        /* ab[2] =  (a[0]*b[1]-a[1]*b[0]); */
        dot_prod_ab_c = dot_prod3( ab, c );
        /* dot_prod_ab_c   =  ab[0]*c[0] + ab[1]*c[1] + ab[2]*c[2]; */
        if (sine_value)
        {
            abs_c = len3( c );
            /* abs_c  = sqrt( c[0]*c[0]   + c[1]*c[1]   + c[2]*c[2] ); */
            abs_ab = len3( ab );
            /* abs_ab = sqrt( ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2] ); */

            if (abs_c > 1.e-7 /* otherwise c has zero length */ && abs_ab > 1.e-7 /* otherwise a is parallel to b*/)
            {
                *sine_value = MPY_SINE * dot_prod_ab_c / ( abs_c * abs_ab );
                /*  *sine_value = dot_prod_ab_c / ( abs_c * abs_ab); */
            }
            else
            {
                *sine_value = 0.0;
            }
        }

        return dot_prod_ab_c;
    }
    */
    // END INCHI C FUNCTION: triple_prod
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: triple_prod
    // INCHI✔️✔️: STEREO_CENTER_BONDS_NORM == 1 selects MPY_SINE == 1.00.
    // END INCHI ACTIVE MACRO CONFIGURATION: triple_prod

    let mut ab = [0.0; 3];
    cross_prod3(a, b, &mut ab);
    let dot_prod_ab_c = dot_prod3(&ab, c);
    if let Some(sine_value) = sine_value {
        let abs_c = len3(c);
        let abs_ab = len3(&ab);
        if abs_c > 1.0e-7 && abs_ab > 1.0e-7 {
            *sine_value = 1.0 * dot_prod_ab_c / (abs_c * abs_ab);
        } else {
            *sine_value = 0.0;
        }
    }
    dot_prod_ab_c
}

pub(crate) fn triple_prod_and_min_abs_sine2(
    at_coord: &[[f64; 3]; 3],
    central_at_coord: &[f64; 3],
    bAddedExplicitNeighbor: i32,
    mut min_sine: Option<&mut f64>,
    mut bAmbiguous: Option<&mut i32>,
    vMinSine: f64,
) -> f64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:926 triple_prod_and_min_abs_sine2
    // INCHI✔️✔️: complete source frame follows verbatim; fixed-size storage, callee composition, and branch costs match.
    /*
    double triple_prod_and_min_abs_sine2( double at_coord[][3],
                                          double central_at_coord[],
                                          int bAddedExplicitNeighbor,
                                          double *min_sine,
                                          int *bAmbiguous,
                                          double vMinSine )
    {
        double min_sine_value = 9999.0, sine_value, min_edge_len, max_edge_len, min_edge_len_NoExplNeigh, max_edge_len_NoExplNeigh;
        double s0, s1, s2, s3, e01, e02, e03, e12, e13, e23, tmp[3], e[3][3];
        double prod, ret, central_prod[4];
        int    bLongEdges;

        if (!min_sine)
        {
            return triple_prod( at_coord[0], at_coord[1], at_coord[2], NULL );
        }

        ret = triple_prod( at_coord[0], at_coord[1], at_coord[2], &sine_value );
        sine_value = MPY_SINE * fabs( sine_value );

        diff3( at_coord[1], at_coord[0], e[2] );
        diff3( at_coord[0], at_coord[2], e[1] );
        diff3( at_coord[2], at_coord[1], e[0] );

        /*  lengths of the 6 edges of the tetrahedra */
        e03 = len3( at_coord[0] ); /* 1 */
        e13 = len3( at_coord[1] );
        e23 = len3( at_coord[2] ); /* includes added neighbor if bAddedExplicitNeighbor*/
        e02 = len3( e[1] );        /* includes added neighbor if bAddedExplicitNeighbor*/
        e12 = len3( e[0] );        /* includes added neighbor if bAddedExplicitNeighbor*/
        e01 = len3( e[2] );

        /*  min & max edge length */
        max_edge_len =
            min_edge_len = e03;

        if (min_edge_len > e13)
        {
            min_edge_len = e13;
        }
        if (min_edge_len > e01)
        {
            min_edge_len = e01;
        }
        min_edge_len_NoExplNeigh = min_edge_len;

        if (min_edge_len > e23)
        {
            min_edge_len = e23;
        }
        if (min_edge_len > e02)
        {
            min_edge_len = e02;
        }
        if (min_edge_len > e12)
        {
            min_edge_len = e12;
        }

        if (max_edge_len < e13)
        {
            max_edge_len = e13;
        }
        if (max_edge_len < e01)
        {
            max_edge_len = e01;
        }
        max_edge_len_NoExplNeigh = max_edge_len;

        if (max_edge_len < e23)
        {
            max_edge_len = e23;
        }
        if (max_edge_len < e02)
        {
            max_edge_len = e02;
        }
        if (max_edge_len < e12)
        {
            max_edge_len = e12;
        }

        if (!bAddedExplicitNeighbor)
        {
            min_edge_len_NoExplNeigh = min_edge_len;
            max_edge_len_NoExplNeigh = max_edge_len;
        }

        bLongEdges = bAddedExplicitNeighbor
            ? ( max_edge_len_NoExplNeigh < MAX_EDGE_RATIO * min_edge_len_NoExplNeigh )
            : ( max_edge_len < MAX_EDGE_RATIO * min_edge_len );

        if (sine_value > vMinSine && ( min_sine || bAmbiguous )) /* djb-rwth: fixing coverity ID #499548 -- unresolved issue -- revision required */
        {
            if (min_sine)
            {
                prod = fabs( ret );
                /*  tetrahedra height = volume(prod) / area of a plane(cross_prod) */
                /*  (instead of a tetrahedra calculate parallelogram/parallelepiped area/volume) */

                /*  4 heights from each of the 4 vertices to the opposite plane */
                s0 = prod / len3((double *)cross_prod3( at_coord[1], at_coord[2], tmp ) ); /* djb-rwth: cast operator added for compatibility */
                s1 = prod / len3((double *)cross_prod3( at_coord[0], at_coord[2], tmp ) ); /* djb-rwth: cast operator added for compatibility */
                s2 = prod / len3((double *)cross_prod3( at_coord[0], at_coord[1], tmp ) ); /* djb-rwth: cast operator added for compatibility */
                s3 = prod / len3((double *)cross_prod3( e[0], e[1], tmp ) ); /* djb-rwth: cast operator added for compatibility */
                /*  abs. value of a sine of an angle between each tetrahedra edge and plane */
                /*  sine = height / edge length */
                if (( sine_value = s0 / e01 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s0 / e02 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s0 / e03 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }

                if (( sine_value = s1 / e01 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s1 / e12 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s1 / e13 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }

                if (( sine_value = s2 / e02 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s2 / e12 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s2 / e23 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }

                if (( sine_value = s3 / e03 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s3 / e13 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                if (( sine_value = s3 / e23 ) < min_sine_value)
                {
                    min_sine_value = sine_value;
                }
                /*  actually use triple sine */
                *min_sine = sine_value = MPY_SINE * min_sine_value;
            }

            if (bAmbiguous && sine_value >= vMinSine)
            {
                /*  Check whether the central atom is outside the tetrahedra (0,0,0), at_coord[0,1,2] */
                /*  compare the tetrahedra volume and the volume of a tetrahedra having central_at_coord[] vertex */
                int i;
                diff3( central_at_coord, at_coord[0], tmp );
                central_prod[0] = triple_prod( at_coord[0], at_coord[1], central_at_coord, NULL );
                central_prod[1] = triple_prod( at_coord[1], at_coord[2], central_at_coord, NULL );
                central_prod[2] = triple_prod( at_coord[2], at_coord[0], central_at_coord, NULL );
                central_prod[3] = triple_prod( e[2], e[1], tmp, NULL );
                for (i = 0; i <= 3; i++)
                {
                    if (central_prod[i] / ret < -MIN_SINE_OUTSIDE)
                    {
                        *bAmbiguous |= AMBIGUOUS_STEREO;
                        break;
                    }
                }
            }
    #if ( STEREO_CENTER_BONDS_NORM == 1 )

            if (bLongEdges && !bAddedExplicitNeighbor && max_edge_len >= MIN_LEN_STRAIGHT)
            {
                /*  possible planar tetragon */
                if (sine_value < MIN_SINE_SQUARE)
                {
                    *min_sine = vMinSine / 2.0; /*  force parity to be undefined */
                    if (bAmbiguous && !*bAmbiguous)
                    {
                        *bAmbiguous |= AMBIGUOUS_STEREO;
                    }
                }
            }

            if (bLongEdges && sine_value < MIN_SINE_SQUARE && sine_value < MIN_SINE_EDGE * min_edge_len_NoExplNeigh)
            {
                *min_sine = vMinSine / 2.0; /*  force parity to be undefined */
                if (bAmbiguous && !*bAmbiguous)
                {
                    *bAmbiguous |= AMBIGUOUS_STEREO;
                }
            }
    #endif
        }
        else
        {
            if (min_sine)
            {
                *min_sine = sine_value;
            }
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: triple_prod_and_min_abs_sine2
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: triple_prod_and_min_abs_sine2
    // INCHI✔️✔️: STEREO_CENTER_BONDS_NORM == 1; MPY_SINE == 1.00 and both long-edge normalization branches are active.
    // END INCHI ACTIVE MACRO CONFIGURATION: triple_prod_and_min_abs_sine2

    if min_sine.is_none() {
        return triple_prod(&at_coord[0], &at_coord[1], &at_coord[2], None);
    }

    let mut sine_value = 0.0;
    let ret = triple_prod(
        &at_coord[0],
        &at_coord[1],
        &at_coord[2],
        Some(&mut sine_value),
    );
    sine_value = sine_value.abs();

    let mut edges = [[0.0; 3]; 3];
    diff3(&at_coord[1], &at_coord[0], &mut edges[2]);
    diff3(&at_coord[0], &at_coord[2], &mut edges[1]);
    diff3(&at_coord[2], &at_coord[1], &mut edges[0]);

    let e03 = len3(&at_coord[0]);
    let e13 = len3(&at_coord[1]);
    let e23 = len3(&at_coord[2]);
    let e02 = len3(&edges[1]);
    let e12 = len3(&edges[0]);
    let e01 = len3(&edges[2]);

    let mut minimum_edge_length = e03;
    let mut maximum_edge_length = e03;
    if minimum_edge_length > e13 {
        minimum_edge_length = e13;
    }
    if minimum_edge_length > e01 {
        minimum_edge_length = e01;
    }
    let mut minimum_edge_length_without_explicit_neighbor = minimum_edge_length;
    if minimum_edge_length > e23 {
        minimum_edge_length = e23;
    }
    if minimum_edge_length > e02 {
        minimum_edge_length = e02;
    }
    if minimum_edge_length > e12 {
        minimum_edge_length = e12;
    }

    if maximum_edge_length < e13 {
        maximum_edge_length = e13;
    }
    if maximum_edge_length < e01 {
        maximum_edge_length = e01;
    }
    let mut maximum_edge_length_without_explicit_neighbor = maximum_edge_length;
    if maximum_edge_length < e23 {
        maximum_edge_length = e23;
    }
    if maximum_edge_length < e02 {
        maximum_edge_length = e02;
    }
    if maximum_edge_length < e12 {
        maximum_edge_length = e12;
    }

    if bAddedExplicitNeighbor == 0 {
        minimum_edge_length_without_explicit_neighbor = minimum_edge_length;
        maximum_edge_length_without_explicit_neighbor = maximum_edge_length;
    }

    let long_edges = if bAddedExplicitNeighbor != 0 {
        maximum_edge_length_without_explicit_neighbor
            < MAX_EDGE_RATIO * minimum_edge_length_without_explicit_neighbor
    } else {
        maximum_edge_length < MAX_EDGE_RATIO * minimum_edge_length
    };

    if sine_value > vMinSine && (min_sine.is_some() || bAmbiguous.is_some()) {
        let product = ret.abs();
        let mut temporary = [0.0; 3];
        let s0 = product / len3(cross_prod3(&at_coord[1], &at_coord[2], &mut temporary));
        let s1 = product / len3(cross_prod3(&at_coord[0], &at_coord[2], &mut temporary));
        let s2 = product / len3(cross_prod3(&at_coord[0], &at_coord[1], &mut temporary));
        let s3 = product / len3(cross_prod3(&edges[0], &edges[1], &mut temporary));

        let mut minimum_sine_value = 9999.0_f64;
        for candidate in [
            s0 / e01,
            s0 / e02,
            s0 / e03,
            s1 / e01,
            s1 / e12,
            s1 / e13,
            s2 / e02,
            s2 / e12,
            s2 / e23,
            s3 / e03,
            s3 / e13,
            s3 / e23,
        ] {
            if candidate < minimum_sine_value {
                minimum_sine_value = candidate;
            }
        }
        sine_value = minimum_sine_value;
        if let Some(output) = min_sine.as_deref_mut() {
            *output = sine_value;
        }

        if let Some(ambiguous) = bAmbiguous.as_deref_mut()
            && sine_value >= vMinSine
        {
            let mut temporary = [0.0; 3];
            diff3(central_at_coord, &at_coord[0], &mut temporary);
            let central_products = [
                triple_prod(&at_coord[0], &at_coord[1], central_at_coord, None),
                triple_prod(&at_coord[1], &at_coord[2], central_at_coord, None),
                triple_prod(&at_coord[2], &at_coord[0], central_at_coord, None),
                triple_prod(&edges[2], &edges[1], &temporary, None),
            ];
            for central_product in central_products {
                if central_product / ret < -MIN_SINE_OUTSIDE {
                    *ambiguous |= AMBIGUOUS_STEREO as i32;
                    break;
                }
            }
        }

        if long_edges
            && bAddedExplicitNeighbor == 0
            && maximum_edge_length >= MIN_LEN_STRAIGHT
            && sine_value < MIN_SINE_SQUARE
        {
            if let Some(output) = min_sine.as_deref_mut() {
                *output = vMinSine / 2.0;
            }
            if let Some(ambiguous) = bAmbiguous.as_deref_mut()
                && *ambiguous == 0
            {
                *ambiguous |= AMBIGUOUS_STEREO as i32;
            }
        }

        if long_edges
            && sine_value < MIN_SINE_SQUARE
            && sine_value < MIN_SINE_EDGE * minimum_edge_length_without_explicit_neighbor
        {
            if let Some(output) = min_sine.as_deref_mut() {
                *output = vMinSine / 2.0;
            }
            if let Some(ambiguous) = bAmbiguous.as_deref_mut()
                && *ambiguous == 0
            {
                *ambiguous |= AMBIGUOUS_STEREO as i32;
            }
        }
    } else if let Some(output) = min_sine.as_deref_mut() {
        *output = sine_value;
    }

    ret
}

pub(crate) fn triple_prod_and_min_abs_sine(
    at_coord: &[[f64; 3]; 3],
    min_sine: Option<&mut f64>,
) -> f64 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1145 triple_prod_and_min_abs_sine
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size callee composition and branching match.
    /*
    double triple_prod_and_min_abs_sine( double at_coord[][3], double *min_sine )
    {
        double min_sine_value = 9999.0, sine_value;
        double prod = 0.0;

        if (!min_sine)
        {
            return triple_prod( at_coord[0], at_coord[1], at_coord[2], NULL );
        }

        prod = triple_prod( at_coord[0], at_coord[1], at_coord[2], &sine_value ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        sine_value = fabs( sine_value );
        min_sine_value = inchi_min( min_sine_value, sine_value );

        prod = triple_prod( at_coord[1], at_coord[2], at_coord[0], &sine_value ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        sine_value = fabs( sine_value );
        min_sine_value = inchi_min( min_sine_value, sine_value );

        prod = triple_prod( at_coord[2], at_coord[0], at_coord[1], &sine_value ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        sine_value = fabs( sine_value );
        min_sine_value = inchi_min( min_sine_value, sine_value );

        *min_sine = min_sine_value;

        return prod;
    }
    */
    // END INCHI C FUNCTION: triple_prod_and_min_abs_sine

    let Some(min_sine) = min_sine else {
        return triple_prod(&at_coord[0], &at_coord[1], &at_coord[2], None);
    };
    let source_min = |first: f64, second: f64| if first < second { first } else { second };
    let mut min_sine_value = 9999.0_f64;
    let mut sine_value = 0.0_f64;
    let _ = triple_prod(
        &at_coord[0],
        &at_coord[1],
        &at_coord[2],
        Some(&mut sine_value),
    );
    min_sine_value = source_min(min_sine_value, sine_value.abs());
    let _ = triple_prod(
        &at_coord[1],
        &at_coord[2],
        &at_coord[0],
        Some(&mut sine_value),
    );
    min_sine_value = source_min(min_sine_value, sine_value.abs());
    let product = triple_prod(
        &at_coord[2],
        &at_coord[0],
        &at_coord[1],
        Some(&mut sine_value),
    );
    min_sine_value = source_min(min_sine_value, sine_value.abs());
    *min_sine = min_sine_value;
    product
}

pub(crate) fn are_3_vect_in_one_plane(at_coord: &[[f64; 3]; 3], min_sine: f64) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1176 are_3_vect_in_one_plane
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size callee composition and comparison match.
    /*
    int are_3_vect_in_one_plane( double at_coord[][3], double min_sine )
    {
        double actual_min_sine;
        double prod; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

        prod = triple_prod_and_min_abs_sine( at_coord, &actual_min_sine ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */

        return actual_min_sine <= min_sine;
    }
    */
    // END INCHI C FUNCTION: are_3_vect_in_one_plane

    let mut actual_min_sine = 0.0_f64;
    let _product = triple_prod_and_min_abs_sine(at_coord, Some(&mut actual_min_sine));
    i32::from(actual_min_sine <= min_sine)
}

pub(crate) fn are_4at_in_one_plane(at_coord: &[[f64; 3]; 4], min_sine: f64) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1190 are_4at_in_one_plane
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-size loops, callee composition, and comparison match.
    /*
    int are_4at_in_one_plane( double at_coord[][3], double min_sine )
    {
        double actual_min_sine, min_actual_min_sine;
        double coord[3][3], prod; /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        int i, k, j;
        for (k = 0; k < 4; k++)
        {
            for (i = j = 0; i < 4; i++)
            {
                if (i != k)
                {
                    diff3( at_coord[i], at_coord[k], coord[j] );
                    j++;
                }
            }

            prod = triple_prod_and_min_abs_sine( coord, &actual_min_sine ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
            if (!k || actual_min_sine < min_actual_min_sine)
            {
                min_actual_min_sine = actual_min_sine;
            }
        }

        return min_actual_min_sine <= min_sine;
    }
    */
    // END INCHI C FUNCTION: are_4at_in_one_plane

    let mut minimum_actual_sine = 0.0_f64;
    for k in 0..4 {
        let mut coordinates = [[0.0_f64; 3]; 3];
        let mut j = 0_usize;
        for i in 0..4 {
            if i != k {
                diff3(&at_coord[i], &at_coord[k], &mut coordinates[j]);
                j += 1;
            }
        }
        let mut actual_minimum_sine = 0.0_f64;
        let _product = triple_prod_and_min_abs_sine(&coordinates, Some(&mut actual_minimum_sine));
        if k == 0 || actual_minimum_sine < minimum_actual_sine {
            minimum_actual_sine = actual_minimum_sine;
        }
    }
    i32::from(minimum_actual_sine <= min_sine)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn triple_prod_char(
    heap: &SourceHeap,
    at: SourceConstPointer<inp_ATOM>,
    at_1: i32,
    i_next_at_1: i32,
    z_dir1: &[S_CHAR; 3],
    at_2: i32,
    i_next_at_2: i32,
    z_dir2: &[S_CHAR; 3],
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1218 triple_prod_char
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and modeled-array copies add overhead.
    /*
    int triple_prod_char( inp_ATOM *at,
                          int at_1,
                          int i_next_at_1,
                          S_CHAR *z_dir1,
                          int at_2,
                          int i_next_at_2,
                          S_CHAR *z_dir2 )
    {
        inp_ATOM *at1, *at2;
        double    pnt[3][3], len;
        int       i;
        int       ret = 0;

        at1 = at + at_1;
        at2 = at + at[at_1].neighbor[i_next_at_1];

        pnt[0][0] = at2->x - at1->x;
        pnt[0][1] = at2->y - at1->y;
        pnt[0][2] = at2->z - at1->z;

        at2 = at + at_2;
        at1 = at + at[at_2].neighbor[i_next_at_2];

        pnt[1][0] = at2->x - at1->x;
        pnt[1][1] = at2->y - at1->y;
        pnt[1][2] = at2->z - at1->z;
        /*
         *  resultant pnt vector directions:
         *
         *         pnt[0]              pnt[1]
         *
         *   [at_1]---->[...]    [...]---->[at_2]
         *
         *
         *  add3 below: (pnt[0] + pnt[1]) -> pnt[1]
         */
        add3( pnt[0], pnt[1], pnt[1] );

        for (i = 0; i < 3; i++)
        {
            pnt[0][i] = (double) z_dir1[i];
            pnt[2][i] = (double) z_dir2[i];
        }
        for (i = 0; i < 3; i++)
        {
            len = len3( pnt[i] );
            if (len < MIN_BOND_LEN)
            {
                if (i == 1 && ( at[at_1].bUsed0DParity || at[at_2].bUsed0DParity ))
                {
                    pnt[i][0] = 0.0;
                    pnt[i][1] = 1.0;
                    pnt[i][2] = 0.0;
                    len = 1.0; /* standard at_1-->at_2 vector coordinates in case of 0D allene */
                }
                else
                {
                    goto exit_function; /*  too short bond */
                }
            }
            mult3( pnt[i], 1.0 / len, pnt[i] );
        }

        len = 100.0*triple_prod( pnt[0], pnt[1], pnt[2], NULL );

        /*
         *   ^ pnt[0]
         *   |                         The orientation on this diagram
         *   |                         produces len = -100
         *  [at_1]------>[at_2]
         *        pnt[1]    /
         *                 /
         *                / pnt[2]  (up from the plane)
         *               v
         *
         * Note: len is invariant upon at_1 <--> at_2 transposition because
         *       triple product changes sign upon pnt[0]<-->pnt[2] transposition and
         *       triple product changes sign upon pnt[1]--> -pnt[1] change of direction:
         *
         * triple_prod(pnt[0],  pnt[1], pnt[2], NULL ) =
         * triple_prod(pnt[2], -pnt[1], pnt[0], NULL )
         *
         */

        ret = len >= 0.0 ? (int) ( floor( len + 0.5 ) ) : -(int) ( floor( 0.5 - len ) );

    exit_function:

        return ret;
    }
    */
    // END INCHI C FUNCTION: triple_prod_char
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: triple_prod_char
    // INCHI✔️❌: MIN_BOND_LEN == 0.000001.
    // END INCHI ACTIVE MACRO CONFIGURATION: triple_prod_char

    const MIN_BOND_LEN: f64 = 0.000001;
    let atoms = heap.slice(at)?;
    let at_1_index = usize::try_from(at_1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let at_2_index = usize::try_from(at_2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let next_1_index =
        usize::try_from(i_next_at_1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let next_2_index =
        usize::try_from(i_next_at_2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let endpoint_1 = atoms
        .get(at_1_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let adjacent_1 = atoms
        .get(usize::from(
            *endpoint_1
                .neighbor
                .get(next_1_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let endpoint_2 = atoms
        .get(at_2_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    let adjacent_2 = atoms
        .get(usize::from(
            *endpoint_2
                .neighbor
                .get(next_2_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        ))
        .ok_or(SourceHeapError::PointerOutOfBounds)?;

    let mut pnt = [[0.0; 3]; 3];
    pnt[0][0] = adjacent_1.x - endpoint_1.x;
    pnt[0][1] = adjacent_1.y - endpoint_1.y;
    pnt[0][2] = adjacent_1.z - endpoint_1.z;
    pnt[1][0] = endpoint_2.x - adjacent_2.x;
    pnt[1][1] = endpoint_2.y - adjacent_2.y;
    pnt[1][2] = endpoint_2.z - adjacent_2.z;
    let first_segment = pnt[0];
    let second_segment = pnt[1];
    add3(&first_segment, &second_segment, &mut pnt[1]);

    for i in 0..3 {
        pnt[0][i] = f64::from(z_dir1[i]);
        pnt[2][i] = f64::from(z_dir2[i]);
    }
    for i in 0..3 {
        let mut length = len3(&pnt[i]);
        if length < MIN_BOND_LEN {
            if i == 1 && (endpoint_1.bUsed0DParity != 0 || endpoint_2.bUsed0DParity != 0) {
                pnt[i] = [0.0, 1.0, 0.0];
                length = 1.0;
            } else {
                return Ok(0);
            }
        }
        let input = pnt[i];
        mult3(&input, 1.0 / length, &mut pnt[i]);
    }

    let length = 100.0 * triple_prod(&pnt[0], &pnt[1], &pnt[2], None);
    Ok(if length >= 0.0 {
        (length + 0.5).floor() as i32
    } else {
        -((0.5 - length).floor() as i32)
    })
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetHalfStereobond0DParity(
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
    nSbNeighOrigAtNumb: &[AT_NUMB],
    nNumExplictAttachments: i32,
    mut bond_parity: i32,
    nFlag: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1725 GetHalfStereobond0DParity
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int GetHalfStereobond0DParity( inp_ATOM *at,
                                   int cur_at,
                                   AT_NUMB nSbNeighOrigAtNumb[],
                                   int nNumExplictAttachments,
                                   int bond_parity,
                                   int nFlag )
    {
        int m, last_parity, cur_parity;
        int i, icur2nxt, icur2neigh, cur_order_parity, nxt_at;
        AT_NUMB nNextSbAtOrigNumb;
        /* find atom parities for all valid streobonds incident to at[cur_at] */
        for (m = 0, last_parity = 0; m < MAX_NUM_STEREO_BONDS && at[cur_at].sb_parity[m]; m++)
        {
            icur2nxt = icur2neigh = -1; /* ordering number of neighbors in nSbNeighOrigAtNumb[] */
            cur_parity = 0;             /* parity for mth stereobond incident to the cur_at */
            nxt_at = at[cur_at].neighbor[(int)at[cur_at].sb_ord[m]]; /* djb-rwth: fixing coverity ID #499549 */
            if (0 <= at[cur_at].sb_ord[m] && at[cur_at].sb_ord[m] < at[cur_at].valence &&
                 at[nxt_at].valence <= MAX_NUM_STEREO_BONDS && /* make sure it is a valid stereobond */
                 ( nNextSbAtOrigNumb = at[nxt_at].orig_at_number ))
            {
                /* since at[cur_at].sn_ord[m] = -1 for explicit H use at[cur_at].sn_orig_at_num[m] */
                for (i = 0; i < nNumExplictAttachments; i++)
                {
                    if (at[cur_at].sn_orig_at_num[m] == nSbNeighOrigAtNumb[i])
                    {
                        icur2neigh = i; /* neighbor */
                    }
                    else
                    {
                        if (nNextSbAtOrigNumb == nSbNeighOrigAtNumb[i])
                        {
                            icur2nxt = i; /* atom connected by a stereobond */
                        }
                    }
                }
                if (icur2neigh >= 0 && icur2nxt >= 0)
                {
                    if (ATOM_PARITY_WELL_DEF( at[cur_at].sb_parity[m] ))
                    {
                        /* parity of at[cur_atom] neighbor permutation to reach this order: { next_atom, neigh_atom, ...} */
                        cur_order_parity = ( icur2nxt + icur2neigh + ( icur2nxt > icur2neigh ) - 1 ) % 2;
                        cur_parity = 2 - ( cur_order_parity + at[cur_at].sb_parity[m] ) % 2;
                    }
                    else
                    {
                                 /* unknowm/undef parities do not depend on the neighbor order */
                        cur_parity = at[cur_at].sb_parity[m];
                    }
                }
            }
            else
            {
                continue;
            }
            /* use a well-known parity if available; if not then use preferably the unknown */
            if (!last_parity)
            {
                last_parity = cur_parity;
            }
            else
            {
                if (last_parity != cur_parity && cur_parity)
                {
                    if (ATOM_PARITY_WELL_DEF( last_parity ))
                    {
                        if (ATOM_PARITY_WELL_DEF( cur_parity ))
                        {
                            last_parity = 0; /* error: all well-defined parities should be same */
                            break;
                        }
                    }
                    else
                    {
                        if (ATOM_PARITY_WELL_DEF( cur_parity ))
                        {
                            /* replace unknown/undefined parity with well-known */
                            last_parity = cur_parity;
                        }
                        else
                        {
                            /* select min unknown/undefined parity (out of AB_PARITY_UNKN and AB_PARITY_UNDF) */
                            last_parity = inchi_min( cur_parity, last_parity );
                        }
                    }
                }
            }
        }

        if (last_parity)
        {
            bond_parity = last_parity;
            at[cur_at].bUsed0DParity |= nFlag; /* set flag: used stereobond 0D parity */
        }

        return bond_parity;
    }
    */
    // END INCHI C FUNCTION: GetHalfStereobond0DParity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetHalfStereobond0DParity
    // INCHI✔️❌: ATOM_PARITY_WELL_DEF(p) is AB_MIN_WELL_DEFINED_PARITY <= p <= AB_MAX_WELL_DEFINED_PARITY.
    // INCHI✔️❌: inchi_min(a, b) selects a when a < b and otherwise selects b.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetHalfStereobond0DParity

    let atoms = heap.slice_mut(at)?;
    let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    if current_index >= atoms.len() {
        return Err(SourceHeapError::PointerOutOfBounds);
    }
    let well_defined = |parity: i32| {
        (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32).contains(&parity)
    };
    let mut last_parity = 0_i32;

    for m in 0..MAX_NUM_STEREO_BONDS as usize {
        let descriptor_parity = i32::from(atoms[current_index].sb_parity[m]);
        if descriptor_parity == 0 {
            break;
        }
        let mut current_to_next = -1_i32;
        let mut current_to_neighbor = -1_i32;
        let mut current_parity = 0_i32;
        let order = i32::from(atoms[current_index].sb_ord[m]);
        let order_index =
            usize::try_from(order).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next_index = usize::from(
            *atoms[current_index]
                .neighbor
                .get(order_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );

        if order >= 0 && order < i32::from(atoms[current_index].valence) {
            let next_atom = atoms
                .get(next_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let next_original_number = next_atom.orig_at_number;
            if i32::from(next_atom.valence) <= MAX_NUM_STEREO_BONDS as i32
                && next_original_number != 0
            {
                for i in 0..nNumExplictAttachments {
                    let attachment = *nSbNeighOrigAtNumb
                        .get(i as usize)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if atoms[current_index].sn_orig_at_num[m] == attachment {
                        current_to_neighbor = i;
                    } else if next_original_number == attachment {
                        current_to_next = i;
                    }
                }
                if current_to_neighbor >= 0 && current_to_next >= 0 {
                    if well_defined(descriptor_parity) {
                        let order_parity = (current_to_next
                            + current_to_neighbor
                            + i32::from(current_to_next > current_to_neighbor)
                            - 1)
                            % 2;
                        current_parity = 2 - (order_parity + descriptor_parity) % 2;
                    } else {
                        current_parity = descriptor_parity;
                    }
                }
            } else {
                continue;
            }
        } else {
            continue;
        }

        if last_parity == 0 {
            last_parity = current_parity;
        } else if last_parity != current_parity && current_parity != 0 {
            if well_defined(last_parity) {
                if well_defined(current_parity) {
                    last_parity = 0;
                    break;
                }
            } else if well_defined(current_parity) {
                last_parity = current_parity;
            } else {
                last_parity = current_parity.min(last_parity);
            }
        }
    }

    if last_parity != 0 {
        bond_parity = last_parity;
        atoms[current_index].bUsed0DParity =
            (i32::from(atoms[current_index].bUsed0DParity) | nFlag) as i8;
    }
    Ok(bond_parity)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn FixSb0DParities(
    heap: &SourceHeap,
    at: SourceConstPointer<inp_ATOM>,
    chain_length: i32,
    at_1: i32,
    i_next_at_1: i32,
    z_dir1: &mut [S_CHAR; 3],
    at_2: i32,
    i_next_at_2: i32,
    z_dir2: &mut [S_CHAR; 3],
    pparity1: &mut i32,
    pparity2: &mut i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1824 FixSb0DParities
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and modeled-array copies add overhead.
    /*
    int FixSb0DParities( inp_ATOM *at,
                         /* inp_ATOM *at_removed_H,
                         int num_removed_H,*/
                         int chain_length,
                         int at_1,
                         int i_next_at_1,
                         S_CHAR z_dir1[],
                         int at_2,
                         int i_next_at_2,
                         S_CHAR z_dir2[],
                         int *pparity1,
                         int *pparity2 )
    {
        int k, parity1, parity2, abs_parity1, abs_parity2;
        int j1, j2, parity_sign;
        /*
        AT_NUMB nSbNeighOrigAtNumb1[MAX_NUM_STEREO_BOND_NEIGH], nSbNeighOrigAtNumb2[MAX_NUM_STEREO_BOND_NEIGH];
        int     nNumExplictAttachments1, nNumExplictAttachments2;
        */
        parity1 = parity2 = AB_PARITY_NONE;
        j1 = j2 = -1;
        parity_sign = ( *pparity1 < 0 || *pparity2 < 0 ) ? -1 : 1;

        abs_parity1 = abs( *pparity1 );
        abs_parity2 = abs( *pparity2 );

        for (k = 0; k < MAX_NUM_STEREO_BONDS && at[at_1].sb_parity[k]; k++)
        {
            if (at[at_1].sb_ord[k] == i_next_at_1)
            {
                parity1 = at[at_1].sb_parity[k];
                j1 = k;
            }
        }
        for (k = 0; k < MAX_NUM_STEREO_BONDS && at[at_2].sb_parity[k]; k++)
        {
            if (at[at_2].sb_ord[k] == i_next_at_2)
            {
                parity2 = at[at_2].sb_parity[k];
                j2 = k;
            }
        }
        switch (( j1 >= 0 ) + 2 * ( j2 >= 0 ))
        {
            case 0:
                /* the bond has no 0D parity */
                *pparity1 = *pparity2 = parity_sign * AB_PARITY_UNDF;
                return 0;
            case 1:
            case 2:
                /* 0D parity data error */
                *pparity1 = *pparity2 = AB_PARITY_NONE;
                return -1;
            case 3:
                /* the bond has 0D parity */
                switch (!( ATOM_PARITY_WELL_DEF( abs_parity1 ) && ATOM_PARITY_WELL_DEF( parity1 ) ) +
                         2 * !( ATOM_PARITY_WELL_DEF( abs_parity2 ) && ATOM_PARITY_WELL_DEF( parity2 ) ))
                {
                    case 0:
                        /* both parities are well-defined; continue */
                        break;
                    case 1:
                        /* 0D parity not well-defined for at_1 */
                        *pparity1 = parity_sign * ( ATOM_PARITY_WELL_DEF( parity1 ) ? abs_parity1 :
                                                   ATOM_PARITY_WELL_DEF( abs_parity1 ) ? parity1 :
                                                   inchi_min( abs_parity1, parity1 ) );
                        *pparity2 = parity_sign * abs_parity2;
                        return -1;
                    case 2:
                        /* 0D parity not well-defined for at_2 */
                        *pparity1 = parity_sign * abs_parity1;
                        *pparity2 = parity_sign * ( ATOM_PARITY_WELL_DEF( parity2 ) ? abs_parity2 :
                                                   ATOM_PARITY_WELL_DEF( abs_parity2 ) ? parity2 :
                                                   inchi_min( abs_parity2, parity2 ) );
                        return -1;
                    case 3:
                        abs_parity1 = ( ATOM_PARITY_WELL_DEF( parity1 ) ? abs_parity1 :
                                        ATOM_PARITY_WELL_DEF( abs_parity1 ) ? parity1 :
                                        inchi_min( abs_parity1, parity1 ) );
                        abs_parity2 = ( ATOM_PARITY_WELL_DEF( parity2 ) ? abs_parity2 :
                                        ATOM_PARITY_WELL_DEF( abs_parity2 ) ? parity2 :
                                        inchi_min( abs_parity2, parity2 ) );
                        *pparity1 = *pparity2 = parity_sign * inchi_min( abs_parity1, abs_parity2 );
                        /*return (parity1 == parity2)? 0 : -1;*/
                        return -1;
                }
                break;
        }
        /* we are here if both end-atoms of the bond have well-defined 0D parities */
        /*
        nNumExplictAttachments1 = GetSbNeighOrigAtNumb( at, at_1, at_removed_H, num_removed_H, nSbNeighOrigAtNumb1 );
        nNumExplictAttachments2 = GetSbNeighOrigAtNumb( at, at_2, at_removed_H, num_removed_H, nSbNeighOrigAtNumb2 );
        parity1 = GetHalfStereobond0DParity( at, at_1, nSbNeighOrigAtNumb1, nNumExplictAttachments1, *pparity1, 0 );
        parity2 = GetHalfStereobond0DParity( at, at_2, nSbNeighOrigAtNumb2, nNumExplictAttachments2, *pparity2, 0 );
        */
        *pparity1 = parity_sign * abs_parity1;
        *pparity2 = parity_sign * abs_parity2;

        if (chain_length % 2)
        {
            /* allene; chain_length = (number of double bonds) - 1 */
            /*
            int zer1 = ( !z_dir1[0] && !z_dir1[1] && !z_dir1[2] );
            int zer2 = ( !z_dir2[0] && !z_dir2[1] && !z_dir2[2] );
            */
            int bWrong_z_dir1 = ( 0 != ( at[at_1].bUsed0DParity & FlagSB_0D ) );
            int bWrong_z_dir2 = ( 0 != ( at[at_2].bUsed0DParity & FlagSB_0D ) );

            if (bWrong_z_dir1 && bWrong_z_dir2)
            {
                goto set_default;
            }
            else
            {
                if (bWrong_z_dir1 || bWrong_z_dir2)
                {
                    double r12[3], zi1[3], zi2[3], abs_r12, abs_zi2;
                    int    at_i1, at_i2, j; /* djb-rwth: ignoring LLVM warning: variables used */
                    S_CHAR   z_dir[3];
                    r12[0] = at[at_2].x - at[at_1].x;
                    r12[1] = at[at_2].y - at[at_1].y;
                    r12[2] = at[at_2].z - at[at_1].z;
                    abs_r12 = len3( r12 );
                    if (abs_r12 < MIN_BOND_LEN)
                    {
                        goto set_default;
                    }
                    /* make r12[] point to the atom with 'good' z_dir[] */
                    if (bWrong_z_dir1)
                    {
                        at_i1 = at_2; /* has good z_dir2[] */ /* djb-rwth: ignoring LLVM warning: variable used */
                        at_i2 = at_1; /* has bad  z_dir1[] */ /* djb-rwth: ignoring LLVM warning: variable used */
                        zi1[0] = z_dir2[0];
                        zi1[1] = z_dir2[1];
                        zi1[2] = z_dir2[2];
                        mult3( r12, 1.0 / abs_r12, r12 ); /* make length = 1 */
                    }
                    else
                    {
                        at_i1 = at_1; /* has good z_dir1[] */ /* djb-rwth: ignoring LLVM warning: variable used */
                        at_i2 = at_2; /* has bad  z_dir2[] */ /* djb-rwth: ignoring LLVM warning: variable used */
                        zi1[0] = z_dir1[0];
                        zi1[1] = z_dir1[1];
                        zi1[2] = z_dir1[2];
                        mult3( r12, -1.0 / abs_r12, r12 ); /* make length = 1 */
                    }
                    cross_prod3( r12, zi1, zi2 );
                    abs_zi2 = len3( zi2 );
                    mult3( zi2, 100.0 / abs_zi2, zi2 ); /* make length = 100 */
                    for (j = 0; j < 3; j++)
                    {
                        z_dir[j] = (S_CHAR) ( zi2[j] >= 0.0 ? floor( 0.5 + zi2[j] ) :
                                                           -floor( 0.5 - zi2[j] ) ); /*  abs(z_dir) = 100 */
                    }
                    if (bWrong_z_dir1)
                    {
                        memcpy(z_dir1, z_dir, sizeof(z_dir));
                    }
                    else
                    {
                        memcpy(z_dir2, z_dir, sizeof(z_dir));
                    }
                }
            }
            return 0;

        set_default:
            /* z_dir1[] = x-direction; z_dir2[] = z-direction; r12[] = y-direction */
            z_dir1[0] = 100;
            z_dir1[1] = z_dir1[2] = 0;
            z_dir2[0] = z_dir2[1] = 0;
            z_dir2[2] = 100;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: FixSb0DParities
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: FixSb0DParities
    // INCHI✔️❌: ATOM_PARITY_WELL_DEF(p) is AB_MIN_WELL_DEFINED_PARITY <= p <= AB_MAX_WELL_DEFINED_PARITY.
    // INCHI✔️❌: MIN_BOND_LEN == 0.000001 and FlagSB_0D == 2.
    // END INCHI ACTIVE MACRO CONFIGURATION: FixSb0DParities

    let atoms = heap.slice(at)?;
    let i1 = usize::try_from(at_1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let i2 = usize::try_from(at_2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let a1 = atoms.get(i1).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let a2 = atoms.get(i2).ok_or(SourceHeapError::PointerOutOfBounds)?;
    let well = |p: i32| {
        (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32).contains(&p)
    };
    let parity_sign = if *pparity1 < 0 || *pparity2 < 0 {
        -1
    } else {
        1
    };
    let mut abs1 = pparity1.wrapping_abs();
    let mut abs2 = pparity2.wrapping_abs();
    let (mut parity1, mut parity2) = (AB_PARITY_NONE as i32, AB_PARITY_NONE as i32);
    let (mut j1, mut j2) = (-1, -1);
    for k in 0..MAX_NUM_STEREO_BONDS as usize {
        if a1.sb_parity[k] == 0 {
            break;
        }
        if i32::from(a1.sb_ord[k]) == i_next_at_1 {
            parity1 = i32::from(a1.sb_parity[k]);
            j1 = k as i32;
        }
    }
    for k in 0..MAX_NUM_STEREO_BONDS as usize {
        if a2.sb_parity[k] == 0 {
            break;
        }
        if i32::from(a2.sb_ord[k]) == i_next_at_2 {
            parity2 = i32::from(a2.sb_parity[k]);
            j2 = k as i32;
        }
    }
    match i32::from(j1 >= 0) + 2 * i32::from(j2 >= 0) {
        0 => {
            let p = parity_sign * AB_PARITY_UNDF as i32;
            *pparity2 = p;
            *pparity1 = p;
            return Ok(0);
        }
        1 | 2 => {
            *pparity2 = AB_PARITY_NONE as i32;
            *pparity1 = AB_PARITY_NONE as i32;
            return Ok(-1);
        }
        _ => {}
    }
    match i32::from(!(well(abs1) && well(parity1))) + 2 * i32::from(!(well(abs2) && well(parity2)))
    {
        0 => {}
        1 => {
            *pparity1 = parity_sign
                * if well(parity1) {
                    abs1
                } else if well(abs1) {
                    parity1
                } else {
                    abs1.min(parity1)
                };
            *pparity2 = parity_sign * abs2;
            return Ok(-1);
        }
        2 => {
            *pparity1 = parity_sign * abs1;
            *pparity2 = parity_sign
                * if well(parity2) {
                    abs2
                } else if well(abs2) {
                    parity2
                } else {
                    abs2.min(parity2)
                };
            return Ok(-1);
        }
        _ => {
            abs1 = if well(parity1) {
                abs1
            } else if well(abs1) {
                parity1
            } else {
                abs1.min(parity1)
            };
            abs2 = if well(parity2) {
                abs2
            } else if well(abs2) {
                parity2
            } else {
                abs2.min(parity2)
            };
            let p = parity_sign * abs1.min(abs2);
            *pparity2 = p;
            *pparity1 = p;
            return Ok(-1);
        }
    }
    *pparity1 = parity_sign * abs1;
    *pparity2 = parity_sign * abs2;
    if chain_length % 2 != 0 {
        let wrong1 = (u32::from(a1.bUsed0DParity as u8) & FlagSB_0D) != 0;
        let wrong2 = (u32::from(a2.bUsed0DParity as u8) & FlagSB_0D) != 0;
        let mut use_default = wrong1 && wrong2;
        if !use_default && (wrong1 || wrong2) {
            let mut r12 = [a2.x - a1.x, a2.y - a1.y, a2.z - a1.z];
            let abs_r12 = len3(&r12);
            if abs_r12 < 0.000001 {
                use_default = true;
            } else {
                let zi1 = if wrong1 {
                    z_dir2.map(f64::from)
                } else {
                    z_dir1.map(f64::from)
                };
                let input = r12;
                mult3(
                    &input,
                    if wrong1 {
                        1.0 / abs_r12
                    } else {
                        -1.0 / abs_r12
                    },
                    &mut r12,
                );
                let mut zi2 = [0.0; 3];
                cross_prod3(&r12, &zi1, &mut zi2);
                let abs_zi2 = len3(&zi2);
                let input = zi2;
                mult3(&input, 100.0 / abs_zi2, &mut zi2);
                let z_dir = zi2.map(|v| {
                    (if v >= 0.0 {
                        (0.5 + v).floor()
                    } else {
                        -(0.5 - v).floor()
                    }) as i8
                });
                if wrong1 {
                    *z_dir1 = z_dir;
                } else {
                    *z_dir2 = z_dir;
                }
            }
        }
        if use_default {
            *z_dir1 = [100, 0, 0];
            *z_dir2 = [0, 0, 100];
        }
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn half_stereo_bond_parity(
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
    at_removed_H: SourceConstPointer<inp_ATOM>,
    num_removed_H: i32,
    mut z_dir: Option<&mut [S_CHAR; 3]>,
    bPointedEdgeStereo: i32,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2121 half_stereo_bond_parity
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and local snapshots add overhead.
    /*
    int half_stereo_bond_parity( inp_ATOM *at,
                                 int cur_at,
                                 inp_ATOM *at_removed_H,
                                 int num_removed_H,
                                 S_CHAR *z_dir,
                                 int bPointedEdgeStereo,
                                 int vABParityUnknown )
    {
        double at_coord[MAX_NUM_STEREO_BOND_NEIGH][3], c, s, tmp[3], tmp1, tmp2, min_tmp, max_tmp, z;
        double temp[3], pnt[3][3];
        int j, k, p0, p1, p2, next, bValence3 = 0, num_z, nType, num_either_single; /* djb-rwth: ignoring LLVM warning: variable used in switch statement; removing redundant variables */
        int nNumExplictAttachments;
        int bond_parity = AB_PARITY_UNDF;
        int    num_H = 0, num_iH, num_eH = 0, num_nH = 0 /* = num_iso_H[0] */;
        int    num_iso_H[NUM_H_ISOTOPES + 1];
        int    index_H[5]; /*  cannot have more than 4 elements: 1 H, 1 1H, 1 D, 1 T atom(s) */
        /*    const double one_pi = 2.0*atan2(1.0 , 0.0 ); */
        const double one_pi = 3.14159265358979323846; /* M_PI */
        const double two_pi = 2.0*one_pi;
        int    bIgnoreIsotopicH = ( 0 != ( at[cur_at].cFlags & AT_FLAG_ISO_H_POINT ) );
        AT_NUMB nSbNeighOrigAtNumb[MAX_NUM_STEREO_BOND_NEIGH];


        if (z_dir && !z_dir[0] && !z_dir[1] && !z_dir[2])
        {
            z_dir[2] = 100;
        }

        num_H = at[cur_at].num_H;
        if (num_H > NUM_H_ISOTOPES)
        {
            return 0; /*  at least 2 H atoms are isotopically identical */
        }

        if (MAX_NUM_STEREO_BOND_NEIGH < at[cur_at].valence + num_H ||
             MIN_NUM_STEREO_BOND_NEIGH > at[cur_at].valence + num_H)
        {
            return 0;
        }

        if (!bCanAtomHaveAStereoBond( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical ))
        {
            return 0;
        }

        if (!bIgnoreIsotopicH)
        {
            for (j = 0, num_nH = num_H; j < NUM_H_ISOTOPES; j++)
            {
                if (( k = (int) at[cur_at].num_iso_H[j] ) > 1)
                {
                    return AB_PARITY_IISO;  /*  two or more identical isotopic H atoms */
                }
                num_nH -= k;
            }
        }

        /*  at this point num_nH = number of non-isotopic H atoms */
        if (num_nH > 1)
        {
            return AB_PARITY_IISO; /*  two or more identical non-isotopic H atoms */
        }

        if (num_nH < 0)
        {
            return CT_ISO_H_ERR;  /*  program error */ /*   <BRKPT> */
        }

        for (j = 0; j < MAX_NUM_STEREO_BOND_NEIGH; j++) /* djb-rwth: avoiding garbage values with proper initialisation */
        {
            for (k = 0; k < 3; k++)
            {
                at_coord[j][k] = 0;
            }
        }

        /********************************************************************
         * Note. At this point all (implicit and explicit) isotopic
         * terminal H neighbors are either different or not present.
         ********************************************************************/

        /*  locate explicit hydrogen atoms */
        /*  (at_removed_H are sorted in ascending isotopic H mass order, non-isotopic first) */
        memset( num_iso_H, 0, sizeof( num_iso_H ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        if (at_removed_H && num_removed_H > 0)
        {
            for (j = 0; j < num_removed_H; j++)
            {
                if (at_removed_H[j].neighbor[0] == cur_at)
                {
                    k = bIgnoreIsotopicH ? 0 : at_removed_H[j].iso_atw_diff;
                    if (0 <= k && k <= NUM_H_ISOTOPES)
                    {
                        if (++num_iso_H[k] > 1)   /*  num_iso_H[0] = number of non-isotopic H atoms */
                            return CT_ISO_H_ERR;    /*  program error in counting hydrogens */ /*   <BRKPT> */
                        index_H[num_eH++] = j;
                    }
                    else
                    {
                        return CT_ISO_H_ERR;  /*  program error */ /*   <BRKPT> */
                    }
                }
            }
            num_iH = num_H - num_eH; /*  number of implicit non-isotopic and isotopic H atoms */
            if (num_iH > 1)
            {
                /*  more than one implicit H: cannot reconstruct the geometry */
                bond_parity = -AB_PARITY_UNDF;
                goto exit_function;
            }
        }
        /* djb-rwth: removing redundant code */
        /*  at this point num_iH = number of implicit non-isotopic and isotopic H atoms */
        if (at[cur_at].valence + num_eH < MIN_NUM_STEREO_BOND_NEIGH)
        {
            /*  =NH or =CHD when no explicit H is present */
            return num_H == 1 ? AB_PARITY_UNDF : -AB_PARITY_UNDF;
        }

        bValence3 = bAtomHasValence3( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical );
        /*
         * Can one explicit hydrogen be added to make asymmetric configuration?
         * For now we can add 1 H atom in case of an appropriate geometry if:
         * (a) one non-isotopic H (even if explicit isotopic H atoms are present), or
         * (b) one isotopic or non-isotopic H if NO explicit isotopic or non-isotopic H atom is present
         * This makes sense only in case chem. valence = 4. In case of chem. valence = 3, do not check.
         */
        if (at[cur_at].valence + num_eH == MIN_NUM_STEREO_BOND_NEIGH && !bValence3 &&
             !(/*(a)*/ (1 == num_nH && !num_iso_H[0]) ||
                 /*(b)*/ (1 == num_H && !num_eH )) /* djb-rwth: addressing LLVM warnings */
           )
        {
            goto exit_function;
            /* return num_H == 1? AB_PARITY_UNDF : -AB_PARITY_UNDF; */
        }

        /*  store neighbors coordinates */
        num_z = num_either_single = 0; /* djb-rwth: ignoring LLVM warning: variable used for switch statement; removing redundant code */
        for (k = nNumExplictAttachments = 0; k < 2; k++)
        {
            switch (k)
            {
                case 0:
                    for (j = 0; j < num_eH; j++, nNumExplictAttachments++)
                    {
                        next = index_H[j];
                        at_coord[nNumExplictAttachments][0] = at_removed_H[next].x - at[cur_at].x;
                        at_coord[nNumExplictAttachments][1] = at_removed_H[next].y - at[cur_at].y;
                        nSbNeighOrigAtNumb[nNumExplictAttachments] = at_removed_H[next].orig_at_number;
                        /* use the fact that (at_removed_H - at) = (number of atoms except removed explicit H) */
                        z = -get_z_coord( at, (int) ( at_removed_H - at ) + next,
                                          0 /*neighbor #*/,
                                          &nType,
                                          -( bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO ) );
                        switch (nType)
                        {
                            case ZTYPE_EITHER:
                                num_either_single++; /*  bond in "Either" direction. */
                                break;
                            case ZTYPE_UP:
                            case ZTYPE_DOWN:
                                nType = -nType; /*  at_removed_H[] contains bonds TO the center, not from */
                                z = len2( at_coord[nNumExplictAttachments] );
                                /*
                                z = sqrt( at_coord[nNumExplictAttachments][0]*at_coord[nNumExplictAttachments][0]
                                        + at_coord[nNumExplictAttachments][1]*at_coord[nNumExplictAttachments][1] );
                                */
                                if (nType == ZTYPE_DOWN)
                                {
                                    z = -z;
                                }
                                /*  no break; here */
                            case ZTYPE_3D:
                                num_z++;
                        }
                        at_coord[nNumExplictAttachments][2] = z;
                    }
                    break;
                case 1:
                    for (j = 0; j < at[cur_at].valence; j++, nNumExplictAttachments++)
                    {
                        next = at[cur_at].neighbor[j];
                        at_coord[nNumExplictAttachments][0] = at[next].x - at[cur_at].x;
                        at_coord[nNumExplictAttachments][1] = at[next].y - at[cur_at].y;
                        nSbNeighOrigAtNumb[nNumExplictAttachments] = at[next].orig_at_number;

                        z = get_z_coord( at, cur_at, j /*neighbor #*/, &nType, ( bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO ) );
                        switch (nType)
                        {
                            case ZTYPE_EITHER:
                                num_either_single++; /*  bond in "Either" direction. */
                                break;
                            case ZTYPE_UP:
                            case ZTYPE_DOWN:
                                z = len2( at_coord[nNumExplictAttachments] );
                                /*
                                z = sqrt( at_coord[nNumExplictAttachments][0]*at_coord[nNumExplictAttachments][0]
                                        + at_coord[nNumExplictAttachments][1]*at_coord[nNumExplictAttachments][1] );
                                */
                                if (nType == ZTYPE_DOWN)
                                {
                                    z = -z;
                                }
                                /*  no break; here */
                            case ZTYPE_3D:
                                num_z++;
                        }
                        at_coord[nNumExplictAttachments][2] = z;
                    }
                    break;
            }
        }

        if (num_either_single)
        {
            bond_parity = vABParityUnknown /*AB_PARITY_UNKN*/;  /*  single bond is 'unknown' */
            goto exit_function;
        }

        /* nNumExplictAttachments is a total number of attachments, including removed explicit terminal hydrogens */
        if (nNumExplictAttachments == 2)
        {
            /*  create coordinates of the implicit hydrogen (or a fictitious atom in case of ==N-X ), */
            /*  coord[2][], attached to the cur_at. */
            for (j = 0; j < 3; j++)
            {
                at_coord[2][j] = -( at_coord[0][j] + at_coord[1][j] );
            }
            nSbNeighOrigAtNumb[nNumExplictAttachments] = 0; /* implicit H or lone pair */
        }
        for (j = 0; j < 3; j++)
        {
            tmp[j] = len3( at_coord[j] );
        }
        min_tmp = inchi_min( tmp[0], inchi_min( tmp[1], tmp[2] ) );
        max_tmp = inchi_max( tmp[0], inchi_max( tmp[1], tmp[2] ) );
        if (min_tmp < MIN_BOND_LEN || min_tmp < MIN_SINE*max_tmp)
        {
            /*  all bonds or some of bonds are too short */
            if (at[cur_at].sb_parity[0])
            {
                /* use bond psrity; the reconciliation in ReconcileAllCmlBondParities()
                * has made all ways to calculate parity produce same result
                */
                bond_parity = GetHalfStereobond0DParity( at, cur_at, nSbNeighOrigAtNumb,
                                                         nNumExplictAttachments, bond_parity, FlagSB_0D );
            }

            goto exit_function;
        }
        /*  normalize lengths to 1 */
        for (j = 0; j < 3; j++)
        {
            mult3( at_coord[j], 1.0 / tmp[j], at_coord[j] );
        }

        /*  find projections of at_coord vector differences on the plane containing their arrowhead ends */
        for (j = 0; j < 3; j++)
        {
            /*  pnt[0..2] = {0-1, 1-2, 2-0} */
            tmp[j] = len3((double *)diff3( at_coord[j], at_coord[( j + 1 ) % 3], pnt[j] ) ); /* djb-rwth: cast operator added for compatibility */
            if (tmp[j] < MIN_SINE)
            {
                goto exit_function; /*  angle #i-cur_at-#j is too small */
            }
            mult3( pnt[j], 1.0 / tmp[j], pnt[j] ); /* 2003-10-06 */
        }
        /*  find pnt[p2], a vector perpendicular to the plane, and its length tmp[p2] */
        /*  replace previous pnt[p2], tmp[p2] with new values; the old values do not have any additional */
        /*  information because pnt[p0]+pnt[p1]+pnt[p2]=0 */
        /*  10-6-2003: a cross-product of one pair pnt[j], pnt[(j+1)%3] can be very small. Find the larges one */
        tmp1 = len3((double *)cross_prod3( pnt[0], pnt[1], temp ) ); /* djb-rwth: cast operator added for compatibility */
        for (j = 1, k = 0; j < 3; j++)
        {
            tmp2 = len3((double *)cross_prod3( pnt[j], pnt[( j + 1 ) % 3], temp ) ); /* djb-rwth: cast operator added for compatibility */
            if (tmp2 > tmp1)
            {
                tmp1 = tmp2;
                k = j;
            }
        }
        /* previously p0=0, p1=1, p2=2 */
        p0 = k;
        p1 = ( k + 1 ) % 3;
        p2 = ( k + 2 ) % 3;
        tmp[p2] = len3((double *)cross_prod3( pnt[p0], pnt[p1], pnt[p2] ) ); /* djb-rwth: cast operator added for compatibility */
        if (tmp[p2] < MIN_SINE*tmp[p0] * tmp[p1])
        {
            goto exit_function; /*  pnt[p0] is almost colinear to pnt[p1] */
        }
        /*  new basis: pnt[p0], pnt[p1], pnt[p2]; set z-coord sign and make abs(pnt[p2]) = 1 */
        mult3( pnt[p2], ( pnt[p2][2] > 0.0 ? 1.0 : -1.0 ) / tmp[p2], pnt[p2] ); /*  unit vector in the new z-axis direction */

        min_tmp = dot_prod3( at_coord[0], pnt[p2] ); /*  non-planarity measure (sine): hight of at_coord[] pyramid */
        mult3( pnt[p2], min_tmp, pnt[p0] ); /*  vector height of the pyramid, ideally 0 */
        /*  find new pnt[p0] = projection of at_coord[p0] on plane orthogonal to pnt[p2] */
        tmp[p0] = len3((double *)diff3( at_coord[0], pnt[p0], pnt[p0] ) ); /* djb-rwth: cast operator added for compatibility */
        mult3( pnt[p0], 1.0 / tmp[p0], pnt[p0] );  /*  new x axis basis vector */
        cross_prod3( pnt[p2], pnt[p0], pnt[p1] ); /*  new y axis basis vector */
        /*  find at_coord in the new basis of {pnt[p0], pnt[p1], pnt[p2]} */
        for (j = 0; j < 3; j++)
        {
            /* copy3(at_coord[j], temp);-- djb-rwth: removing copy3 function */
            memcpy(temp, at_coord[j], sizeof(at_coord[j]));
            for (k = 0; k < 3; k++)
            {
                at_coord[j][k] = dot_prod3( temp, pnt[( k + p0 ) % 3] );
            }
            /*  new xy plane projection length */
            tmp[j] = sqrt( at_coord[j][0] * at_coord[j][0] + at_coord[j][1] * at_coord[j][1] );
            /*  make new xy plane projection length = 1 */
            mult3( at_coord[j], 1.0 / tmp[j], at_coord[j] );
        }

        s = fabs( at_coord[1][0] * at_coord[2][1] - at_coord[1][1] * at_coord[2][0] ); /*  1-2 sine */
        c = at_coord[1][0] * at_coord[2][0] + at_coord[1][1] * at_coord[2][1];   /*  1-2 cosine */
        if (s < MIN_SINE && c > 0.5)
        {
            goto exit_function; /*  bonds to neigh. 1 and 2 have almost same direction; relative angles are undefined */
        }
        c = at_coord[0][0]; /*  cosine of the angle between new Ox axis and a bond to the neighbor 0. Should be 1 */
        s = at_coord[0][1]; /*  sine. Should be 0 */
        /*  turn vectors so that vector #1 (at_coord[0]) becomes {1, 0} */
        for (j = 0; j < MAX_NUM_STEREO_BOND_NEIGH; j++)
        {
            tmp1 = c*at_coord[j][0] + s*at_coord[j][1];
            tmp2 = -s*at_coord[j][0] + c*at_coord[j][1];
            at_coord[j][0] = tmp1;
            at_coord[j][1] = tmp2;
        }
        /*  counterclockwise angles from the direction to neigh 0 to to directions to neighbors 1 and 2: */
        tmp1 = atan2( at_coord[1][1], at_coord[1][0] ); /*  range -pi and +pi */
        tmp2 = atan2( at_coord[2][1], at_coord[2][0] );
        if (tmp1 < 0.0)
        {
            tmp1 += two_pi; /*  range 0 to 2*pi */
        }
        if (tmp2 < 0.0)
        {
            tmp2 += two_pi;
        }
        /*-----------------------------------
                            Example
          1 \               case tmp1 < tmp2
             \              parity is odd
              \             (counterclockwise)
               A------- 0
              /
             /
          2 /

        ------------------------------------*/
        bond_parity = 2 - ( tmp1 < tmp2 );
        for (j = 0; j < 3; j++)
        {
            if (z_dir)
                z_dir[j] = (S_CHAR) ( pnt[p2][j] >= 0.0 ? floor( 0.5 + 100.0 * pnt[p2][j] ) :
                                                   -floor( 0.5 - 100.0 * pnt[p2][j] ) ); /*  abs(z_dir) = 100 */
        }
        /*  check for ambiguity */
        if (nNumExplictAttachments > 2)
        {
            min_tmp = inchi_min( tmp1, tmp2 );
            max_tmp = inchi_max( tmp1, tmp2 );
            if (min_tmp > one_pi - MIN_SINE || max_tmp < one_pi + MIN_SINE || max_tmp - min_tmp > one_pi - MIN_SINE)
            {
                at[cur_at].bAmbiguousStereo |= AMBIGUOUS_STEREO;
            }
            else /* 3D ambiguity 8-28-2002 */
            {
                if (fabs( at_coord[0][2] ) > MAX_SINE)
                { /*  all fabs(at_coord[j][2] (j=0..2) must be equal */
                    at[cur_at].bAmbiguousStereo |= AMBIGUOUS_STEREO;
                }
            }
        }
        else
        {
            if (nNumExplictAttachments == 2)
            {  /* 10-6-2003: added */
                min_tmp = fabs( tmp1 - one_pi );
                if (min_tmp < MIN_SINE)
                {
                    bond_parity = AB_PARITY_UNDF; /* consider as undefined 10-6-2003 */
                }
                else
                {
                    if (min_tmp < MIN_ANGLE_DBOND)
                    {
                        at[cur_at].bAmbiguousStereo |= AMBIGUOUS_STEREO;
                    }
                }
            }
        }

        /*  for 3 neighbors moving implicit H to the index=0 from index=2 position */
        /*  can be done in 2 transpositions and does not change atom's parity */

    exit_function:

        if (num_H > 1 && bond_parity > 0 && !( bond_parity & AB_PARITY_0D ) /*&& PARITY_WELL_DEF(bond_parity)*/)
        {
            /*
             * stereo only if isotopes are counted.             Do not inverse
             * Examples:                                        sign for this:
             *     H                            D
             *    /                            /                    H
             * ==C                      or  ==CH                   /
             *    \                                             ==N  (bValence3=1)
             *     D
             * two explicit         one explicit H isotope (D),
             * isotopic H atoms     one implicit H
             */
            bond_parity = -bond_parity; /*  refers to isotopically substituted structure only */
        }

        return bond_parity;
    }
    */
    // END INCHI C FUNCTION: half_stereo_bond_parity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: half_stereo_bond_parity
    // INCHI✔️❌: MAX_NUM_STEREO_BOND_NEIGH == 3, MIN_NUM_STEREO_BOND_NEIGH == 2, and NUM_H_ISOTOPES == 3.
    // INCHI✔️❌: ATOM_PARITY_WELL_DEF, inchi_min, and inchi_max use the configured source macro behavior.
    // INCHI✔️❌: PES_BIT_POINT_EDGE_STEREO == 1, AT_FLAG_ISO_H_POINT == 1, and FlagSB_0D == 2.
    // END INCHI ACTIVE MACRO CONFIGURATION: half_stereo_bond_parity

    let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(at.as_const())?
        .get(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if let Some(direction) = z_dir.as_deref_mut()
        && direction == &[0, 0, 0]
    {
        direction[2] = 100;
    }

    let num_h = i32::from(current.num_H);
    if num_h > NUM_H_ISOTOPES as i32 {
        return Ok(0);
    }
    let attachment_count = i32::from(current.valence) + num_h;
    if attachment_count > MAX_NUM_STEREO_BOND_NEIGH as i32
        || attachment_count < MIN_NUM_STEREO_BOND_NEIGH as i32
    {
        return Ok(0);
    }
    if bCanAtomHaveAStereoBond(&current.elname, current.charge, current.radical)? == 0 {
        return Ok(0);
    }

    let ignore_isotopic_h = (i32::from(current.cFlags) & AT_FLAG_ISO_H_POINT as i32) != 0;
    let mut num_nonisotopic_h = 0_i32;
    if !ignore_isotopic_h {
        num_nonisotopic_h = num_h;
        for &count in &current.num_iso_H {
            let count = i32::from(count);
            if count > 1 {
                return Ok(AB_PARITY_IISO as i32);
            }
            num_nonisotopic_h -= count;
        }
    }
    if num_nonisotopic_h > 1 {
        return Ok(AB_PARITY_IISO as i32);
    }
    if num_nonisotopic_h < 0 {
        return Ok(CT_ISO_H_ERR);
    }

    let mut at_coord = [[0.0_f64; 3]; MAX_NUM_STEREO_BOND_NEIGH as usize];
    let mut neighbor_original_numbers = [0_u16; MAX_NUM_STEREO_BOND_NEIGH as usize];
    let mut explicit_isotope_counts = [0_i32; NUM_H_ISOTOPES as usize + 1];
    let mut hydrogen_indices = [0_usize; 5];
    let mut num_explicit_h = 0_i32;
    let mut bond_parity = AB_PARITY_UNDF as i32;

    if !at_removed_H.is_null() && num_removed_H > 0 {
        let removed = heap.slice(at_removed_H)?;
        for j in 0..num_removed_H {
            let removed_atom = removed
                .get(j as usize)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if i32::from(removed_atom.neighbor[0]) == cur_at {
                let isotope = if ignore_isotopic_h {
                    0
                } else {
                    i32::from(removed_atom.iso_atw_diff)
                };
                if (0..=NUM_H_ISOTOPES as i32).contains(&isotope) {
                    explicit_isotope_counts[isotope as usize] += 1;
                    if explicit_isotope_counts[isotope as usize] > 1 {
                        return Ok(CT_ISO_H_ERR);
                    }
                    hydrogen_indices[num_explicit_h as usize] = j as usize;
                    num_explicit_h += 1;
                } else {
                    return Ok(CT_ISO_H_ERR);
                }
            }
        }
        let num_implicit_h = num_h - num_explicit_h;
        if num_implicit_h > 1 {
            bond_parity = -(AB_PARITY_UNDF as i32);
            return Ok(finalize_half_stereo_parity(num_h, bond_parity));
        }
    }

    if i32::from(current.valence) + num_explicit_h < MIN_NUM_STEREO_BOND_NEIGH as i32 {
        return Ok(if num_h == 1 {
            AB_PARITY_UNDF as i32
        } else {
            -(AB_PARITY_UNDF as i32)
        });
    }
    let valence_three = bAtomHasValence3(&current.elname, current.charge, current.radical)?;
    if i32::from(current.valence) + num_explicit_h == MIN_NUM_STEREO_BOND_NEIGH as i32
        && valence_three == 0
        && !((num_nonisotopic_h == 1 && explicit_isotope_counts[0] == 0)
            || (num_h == 1 && num_explicit_h == 0))
    {
        return Ok(finalize_half_stereo_parity(num_h, bond_parity));
    }

    let mut number_of_z = 0_i32;
    let mut number_of_either = 0_i32;
    let mut number_of_explicit_attachments = 0_usize;
    if num_explicit_h > 0 {
        let removed = heap.slice(at_removed_H)?.to_vec();
        let removed_offset = at_removed_H.as_mut().difference(at)?;
        for j in 0..num_explicit_h as usize {
            let next = hydrogen_indices[j];
            let removed_atom = removed
                .get(next)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            at_coord[number_of_explicit_attachments][0] = removed_atom.x - current.x;
            at_coord[number_of_explicit_attachments][1] = removed_atom.y - current.y;
            neighbor_original_numbers[number_of_explicit_attachments] = removed_atom.orig_at_number;
            let removed_index = removed_offset
                .checked_add(next as i64)
                .and_then(|value| i32::try_from(value).ok())
                .ok_or(SourceHeapError::PointerDifferenceOverflow)?;
            let mut coordinate_type = ZTYPE_NONE as i32;
            let mut z = -get_z_coord(
                heap,
                at.as_const(),
                removed_index,
                0,
                &mut coordinate_type,
                -(bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO as i32),
            )?;
            match coordinate_type {
                value if value == ZTYPE_EITHER as i32 => number_of_either += 1,
                value if value == ZTYPE_UP as i32 || value == ZTYPE_DOWN => {
                    coordinate_type = -coordinate_type;
                    z = len2(&[
                        at_coord[number_of_explicit_attachments][0],
                        at_coord[number_of_explicit_attachments][1],
                    ]);
                    if coordinate_type == ZTYPE_DOWN {
                        z = -z;
                    }
                    number_of_z += 1;
                }
                value if value == ZTYPE_3D as i32 => number_of_z += 1,
                _ => {}
            }
            at_coord[number_of_explicit_attachments][2] = z;
            number_of_explicit_attachments += 1;
        }
    }
    for j in 0..usize::try_from(current.valence).map_err(|_| SourceHeapError::PointerOutOfBounds)? {
        let next = usize::from(current.neighbor[j]);
        let neighbor = heap
            .slice(at.as_const())?
            .get(next)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .clone();
        at_coord[number_of_explicit_attachments][0] = neighbor.x - current.x;
        at_coord[number_of_explicit_attachments][1] = neighbor.y - current.y;
        neighbor_original_numbers[number_of_explicit_attachments] = neighbor.orig_at_number;
        let mut coordinate_type = ZTYPE_NONE as i32;
        let mut z = get_z_coord(
            heap,
            at.as_const(),
            cur_at,
            j as i32,
            &mut coordinate_type,
            bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO as i32,
        )?;
        match coordinate_type {
            value if value == ZTYPE_EITHER as i32 => number_of_either += 1,
            value if value == ZTYPE_UP as i32 || value == ZTYPE_DOWN => {
                z = len2(&[
                    at_coord[number_of_explicit_attachments][0],
                    at_coord[number_of_explicit_attachments][1],
                ]);
                if coordinate_type == ZTYPE_DOWN {
                    z = -z;
                }
                number_of_z += 1;
            }
            value if value == ZTYPE_3D as i32 => number_of_z += 1,
            _ => {}
        }
        at_coord[number_of_explicit_attachments][2] = z;
        number_of_explicit_attachments += 1;
    }
    let _ = number_of_z;

    if number_of_either != 0 {
        bond_parity = vABParityUnknown;
        return Ok(finalize_half_stereo_parity(num_h, bond_parity));
    }
    if number_of_explicit_attachments == 2 {
        for j in 0..3 {
            at_coord[2][j] = -(at_coord[0][j] + at_coord[1][j]);
        }
        neighbor_original_numbers[2] = 0;
    }

    let mut lengths = at_coord.map(|coordinate| len3(&coordinate));
    let source_min = |a: f64, b: f64| if a < b { a } else { b };
    let source_max = |a: f64, b: f64| if a > b { a } else { b };
    let minimum_length = source_min(lengths[0], source_min(lengths[1], lengths[2]));
    let maximum_length = source_max(lengths[0], source_max(lengths[1], lengths[2]));
    if minimum_length < MIN_BOND_LEN || minimum_length < MIN_SINE * maximum_length {
        if current.sb_parity[0] != 0 {
            bond_parity = GetHalfStereobond0DParity(
                heap,
                at,
                cur_at,
                &neighbor_original_numbers,
                number_of_explicit_attachments as i32,
                bond_parity,
                FlagSB_0D as i32,
            )?;
        }
        return Ok(finalize_half_stereo_parity(num_h, bond_parity));
    }
    for j in 0..3 {
        let input = at_coord[j];
        mult3(&input, 1.0 / lengths[j], &mut at_coord[j]);
    }

    let mut points = [[0.0_f64; 3]; 3];
    for j in 0..3 {
        diff3(&at_coord[j], &at_coord[(j + 1) % 3], &mut points[j]);
        lengths[j] = len3(&points[j]);
        if lengths[j] < MIN_SINE {
            return Ok(finalize_half_stereo_parity(num_h, bond_parity));
        }
        let input = points[j];
        mult3(&input, 1.0 / lengths[j], &mut points[j]);
    }

    let mut temporary = [0.0_f64; 3];
    cross_prod3(&points[0], &points[1], &mut temporary);
    let mut largest_cross_length = len3(&temporary);
    let mut largest_cross_index = 0_usize;
    for j in 1..3 {
        cross_prod3(&points[j], &points[(j + 1) % 3], &mut temporary);
        let candidate = len3(&temporary);
        if candidate > largest_cross_length {
            largest_cross_length = candidate;
            largest_cross_index = j;
        }
    }
    let p0 = largest_cross_index;
    let p1 = (largest_cross_index + 1) % 3;
    let p2 = (largest_cross_index + 2) % 3;
    let first = points[p0];
    let second = points[p1];
    cross_prod3(&first, &second, &mut points[p2]);
    lengths[p2] = len3(&points[p2]);
    if lengths[p2] < MIN_SINE * lengths[p0] * lengths[p1] {
        return Ok(finalize_half_stereo_parity(num_h, bond_parity));
    }
    let input = points[p2];
    mult3(
        &input,
        (if points[p2][2] > 0.0 { 1.0 } else { -1.0 }) / lengths[p2],
        &mut points[p2],
    );

    let pyramid_height = dot_prod3(&at_coord[0], &points[p2]);
    let input = points[p2];
    mult3(&input, pyramid_height, &mut points[p0]);
    let subtraction = points[p0];
    diff3(&at_coord[0], &subtraction, &mut points[p0]);
    lengths[p0] = len3(&points[p0]);
    let input = points[p0];
    mult3(&input, 1.0 / lengths[p0], &mut points[p0]);
    let z_axis = points[p2];
    let x_axis = points[p0];
    cross_prod3(&z_axis, &x_axis, &mut points[p1]);

    for coordinate in &mut at_coord {
        let original = *coordinate;
        for k in 0..3 {
            coordinate[k] = dot_prod3(&original, &points[(k + p0) % 3]);
        }
        let projection_length =
            (coordinate[0] * coordinate[0] + coordinate[1] * coordinate[1]).sqrt();
        let input = *coordinate;
        mult3(&input, 1.0 / projection_length, coordinate);
    }

    let sine_12 = (at_coord[1][0] * at_coord[2][1] - at_coord[1][1] * at_coord[2][0]).abs();
    let cosine_12 = at_coord[1][0] * at_coord[2][0] + at_coord[1][1] * at_coord[2][1];
    if sine_12 < MIN_SINE && cosine_12 > 0.5 {
        return Ok(finalize_half_stereo_parity(num_h, bond_parity));
    }
    let cosine = at_coord[0][0];
    let sine = at_coord[0][1];
    for coordinate in &mut at_coord {
        let x = cosine * coordinate[0] + sine * coordinate[1];
        let y = -sine * coordinate[0] + cosine * coordinate[1];
        coordinate[0] = x;
        coordinate[1] = y;
    }
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut angle1 = at_coord[1][1].atan2(at_coord[1][0]);
    let mut angle2 = at_coord[2][1].atan2(at_coord[2][0]);
    if angle1 < 0.0 {
        angle1 += two_pi;
    }
    if angle2 < 0.0 {
        angle2 += two_pi;
    }
    bond_parity = 2 - i32::from(angle1 < angle2);
    if let Some(direction) = z_dir.as_deref_mut() {
        for j in 0..3 {
            direction[j] = (if points[p2][j] >= 0.0 {
                (0.5 + 100.0 * points[p2][j]).floor()
            } else {
                -(0.5 - 100.0 * points[p2][j]).floor()
            }) as i8;
        }
    }

    if number_of_explicit_attachments > 2 {
        let minimum_angle = source_min(angle1, angle2);
        let maximum_angle = source_max(angle1, angle2);
        if minimum_angle > std::f64::consts::PI - MIN_SINE
            || maximum_angle < std::f64::consts::PI + MIN_SINE
            || maximum_angle - minimum_angle > std::f64::consts::PI - MIN_SINE
            || at_coord[0][2].abs() > MAX_SINE
        {
            set_half_stereo_ambiguous(heap, at, current_index)?;
        }
    } else if number_of_explicit_attachments == 2 {
        let distance_from_pi = (angle1 - std::f64::consts::PI).abs();
        if distance_from_pi < MIN_SINE {
            bond_parity = AB_PARITY_UNDF as i32;
        } else if distance_from_pi < MIN_ANGLE_DBOND {
            set_half_stereo_ambiguous(heap, at, current_index)?;
        }
    }

    Ok(finalize_half_stereo_parity(num_h, bond_parity))
}

fn finalize_half_stereo_parity(num_h: i32, bond_parity: i32) -> i32 {
    if num_h > 1 && bond_parity > 0 && (bond_parity & AB_PARITY_0D as i32) == 0 {
        -bond_parity
    } else {
        bond_parity
    }
}

fn set_half_stereo_ambiguous(
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    current_index: usize,
) -> Result<(), SourceHeapError> {
    let current = heap
        .slice_mut(at)?
        .get_mut(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    current.bAmbiguousStereo =
        (i32::from(current.bAmbiguousStereo) | AMBIGUOUS_STEREO as i32) as i8;
    Ok(())
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn save_a_stereo_bond(
    z_prod: i32,
    result_action: i32,
    at1: i32,
    ord1: i32,
    stereo_bond_neighbor1: &mut [AT_NUMB; MAX_NUM_STEREO_BONDS as usize],
    stereo_bond_ord1: &mut [S_CHAR; MAX_NUM_STEREO_BONDS as usize],
    stereo_bond_z_prod1: &mut [S_CHAR; MAX_NUM_STEREO_BONDS as usize],
    stereo_bond_parity1: &mut [S_CHAR; MAX_NUM_STEREO_BONDS as usize],
    at2: i32,
    ord2: i32,
    stereo_bond_neighbor2: &mut [AT_NUMB; MAX_NUM_STEREO_BONDS as usize],
    stereo_bond_ord2: &mut [S_CHAR; MAX_NUM_STEREO_BONDS as usize],
    stereo_bond_z_prod2: &mut [S_CHAR; MAX_NUM_STEREO_BONDS as usize],
    stereo_bond_parity2: &mut [S_CHAR; MAX_NUM_STEREO_BONDS as usize],
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2542 save_a_stereo_bond
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-array scans, writes, and allocation behavior match.
    /*
    int save_a_stereo_bond( int z_prod,
                            int result_action,
                            int at1,
                            int ord1,
                            AT_NUMB *stereo_bond_neighbor1,
                            S_CHAR *stereo_bond_ord1,
                            S_CHAR *stereo_bond_z_prod1,
                            S_CHAR *stereo_bond_parity1,
                            int at2,
                            int ord2,
                            AT_NUMB *stereo_bond_neighbor2,
                            S_CHAR *stereo_bond_ord2,
                            S_CHAR *stereo_bond_z_prod2,
                            S_CHAR *stereo_bond_parity2 )
    {
        int i1, i2;

        for (i1 = 0; i1 < MAX_NUM_STEREO_BONDS && stereo_bond_neighbor1[i1]; i1++)
        {
            ;
        }
        for (i2 = 0; i2 < MAX_NUM_STEREO_BONDS && stereo_bond_neighbor2[i2]; i2++)
        {
            ;
        }

        if (i1 == MAX_NUM_STEREO_BONDS || i2 == MAX_NUM_STEREO_BONDS)
        {
            return 0;
        }

        stereo_bond_parity1[i1] =
            stereo_bond_parity2[i2] = result_action;

        stereo_bond_neighbor1[i1] = (AT_NUMB) ( at2 + 1 );
        stereo_bond_ord1[i1] = (S_CHAR) ord1;
        stereo_bond_neighbor2[i2] = (AT_NUMB) ( at1 + 1 );
        stereo_bond_ord2[i2] = (S_CHAR) ord2;
        stereo_bond_z_prod1[i1] =
            stereo_bond_z_prod2[i2] = (S_CHAR) z_prod;

        return 1;
    }
    */
    // END INCHI C FUNCTION: save_a_stereo_bond

    let first_slot = stereo_bond_neighbor1
        .iter()
        .position(|&neighbor| neighbor == 0);
    let second_slot = stereo_bond_neighbor2
        .iter()
        .position(|&neighbor| neighbor == 0);
    let (Some(i1), Some(i2)) = (first_slot, second_slot) else {
        return 0;
    };

    stereo_bond_parity2[i2] = result_action as S_CHAR;
    stereo_bond_parity1[i1] = result_action as S_CHAR;
    stereo_bond_neighbor1[i1] = at2.wrapping_add(1) as AT_NUMB;
    stereo_bond_ord1[i1] = ord1 as S_CHAR;
    stereo_bond_neighbor2[i2] = at1.wrapping_add(1) as AT_NUMB;
    stereo_bond_ord2[i2] = ord2 as S_CHAR;
    stereo_bond_z_prod2[i2] = z_prod as S_CHAR;
    stereo_bond_z_prod1[i1] = z_prod as S_CHAR;
    1
}

#[allow(non_snake_case)]
pub(crate) fn half_stereo_bond_action(
    mut nParity: i32,
    bUnknown: i32,
    bIsotopic: i32,
    vABParityUnknown: i32,
) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2901 half_stereo_bond_action
    // INCHI✔️✔️: complete active source frame follows verbatim; scalar control flow and allocation behavior match.
    /*
    int half_stereo_bond_action( int nParity,
                                 int bUnknown,
                                 int bIsotopic,
                                 int vABParityUnknown )
    {
    #define AB_NEGATIVE 0x10
    #define AB_UNKNOWN  0x20
        int nAction;

        if (nParity == AB_PARITY_NONE)
        {
            return AB_PARITY_NONE;
        }

        /*  Unknown (type 1) in the parity value may come from the 'Either' single bond only */
        /*  Treat it as a known single bond geometry and unknown (Either) double bond */
        if (nParity == vABParityUnknown /*AB_PARITY_UNKN*/)
        {
            nParity = AB_PARITY_ODD | AB_UNKNOWN;
        }
        if (nParity == -vABParityUnknown /*AB_PARITY_UNKN*/)
        {
            nParity = AB_PARITY_ODD | AB_UNKNOWN | AB_NEGATIVE;
        }

        /*  make positive, replace AB_PARITY_EVEN with AB_PARITY_ODD  */
        if (nParity < 0)
        {
            nParity = ( ( nParity == -AB_PARITY_EVEN ) ? AB_PARITY_ODD : ( -nParity ) ) | AB_NEGATIVE;
        }
        else
        {
            if (nParity == AB_PARITY_EVEN)
            {
                nParity = AB_PARITY_ODD;
            }
        }

        /*  Unknown (type 2): was detected in the double bond attribute */
        /*  (this 'unknown' came from 'Either' double bond) */
        /*  Treat both unknowns in the same way */
        if (bUnknown)
        {
            nParity |= AB_UNKNOWN;
        }

        if (bIsotopic)
        {
            switch (nParity)
            {
                case AB_PARITY_ODD:
                case AB_PARITY_ODD | AB_NEGATIVE:
                    nAction = AB_PARITY_CALC;
                    break;
                case AB_PARITY_ODD | AB_UNKNOWN:
                case AB_PARITY_UNDF | AB_UNKNOWN:
                case AB_PARITY_ODD | AB_UNKNOWN | AB_NEGATIVE:
                case AB_PARITY_UNDF | AB_UNKNOWN | AB_NEGATIVE:
                    nAction = vABParityUnknown /*AB_PARITY_UNKN*/;
                    break;
                case AB_PARITY_IISO:
                case AB_PARITY_IISO | AB_UNKNOWN:
                    nAction = AB_PARITY_NONE;
                    break;
                case AB_PARITY_UNDF:
                case AB_PARITY_UNDF | AB_NEGATIVE:
                    nAction = AB_PARITY_UNDF;
                    break;
                default:
                    nAction = -1; /*  program error */
            }
        }
        else
        {
         /*  Non-isotopic */
            switch (nParity)
            {
                case AB_PARITY_ODD:
                    nAction = AB_PARITY_CALC;
                    break;
                case AB_PARITY_ODD | AB_UNKNOWN:
                case AB_PARITY_UNDF | AB_UNKNOWN:
                    nAction = vABParityUnknown /*AB_PARITY_UNKN*/;
                    break;
                /* case AB_PARITY_ODD  | AB_UNKNOWN | AB_NEGATIVE: */
                case AB_PARITY_UNDF:
                    nAction = AB_PARITY_UNDF;
                    break;
                case AB_PARITY_ODD | AB_UNKNOWN | AB_NEGATIVE:
                case AB_PARITY_ODD | AB_NEGATIVE:
                case AB_PARITY_IISO:
                case AB_PARITY_IISO | AB_UNKNOWN:
                case AB_PARITY_UNDF | AB_NEGATIVE:
                case AB_PARITY_UNDF | AB_UNKNOWN | AB_NEGATIVE:
                    nAction = AB_PARITY_NONE;
                    break;
                default:
                    nAction = -1; /*  program error */
            }
        }

        return nAction;
    #undef AB_NEGATIVE
    #undef AB_UNKNOWN
    }
    */
    // END INCHI C FUNCTION: half_stereo_bond_action
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: half_stereo_bond_action
    // INCHI✔️✔️: function-local AB_NEGATIVE == 0x10 and AB_UNKNOWN == 0x20.
    // END INCHI ACTIVE MACRO CONFIGURATION: half_stereo_bond_action

    const AB_NEGATIVE: i32 = 0x10;
    const AB_UNKNOWN: i32 = 0x20;

    if nParity == AB_PARITY_NONE as i32 {
        return AB_PARITY_NONE as i32;
    }
    if nParity == vABParityUnknown {
        nParity = AB_PARITY_ODD as i32 | AB_UNKNOWN;
    }
    if nParity == -vABParityUnknown {
        nParity = AB_PARITY_ODD as i32 | AB_UNKNOWN | AB_NEGATIVE;
    }
    if nParity < 0 {
        nParity = (if nParity == -(AB_PARITY_EVEN as i32) {
            AB_PARITY_ODD as i32
        } else {
            nParity.wrapping_neg()
        }) | AB_NEGATIVE;
    } else if nParity == AB_PARITY_EVEN as i32 {
        nParity = AB_PARITY_ODD as i32;
    }
    if bUnknown != 0 {
        nParity |= AB_UNKNOWN;
    }

    if bIsotopic != 0 {
        match nParity {
            value
                if value == AB_PARITY_ODD as i32
                    || value == (AB_PARITY_ODD as i32 | AB_NEGATIVE) =>
            {
                AB_PARITY_CALC as i32
            }
            value
                if value == (AB_PARITY_ODD as i32 | AB_UNKNOWN)
                    || value == (AB_PARITY_UNDF as i32 | AB_UNKNOWN)
                    || value == (AB_PARITY_ODD as i32 | AB_UNKNOWN | AB_NEGATIVE)
                    || value == (AB_PARITY_UNDF as i32 | AB_UNKNOWN | AB_NEGATIVE) =>
            {
                vABParityUnknown
            }
            value
                if value == AB_PARITY_IISO as i32
                    || value == (AB_PARITY_IISO as i32 | AB_UNKNOWN) =>
            {
                AB_PARITY_NONE as i32
            }
            value
                if value == AB_PARITY_UNDF as i32
                    || value == (AB_PARITY_UNDF as i32 | AB_NEGATIVE) =>
            {
                AB_PARITY_UNDF as i32
            }
            _ => -1,
        }
    } else {
        match nParity {
            value if value == AB_PARITY_ODD as i32 => AB_PARITY_CALC as i32,
            value
                if value == (AB_PARITY_ODD as i32 | AB_UNKNOWN)
                    || value == (AB_PARITY_UNDF as i32 | AB_UNKNOWN) =>
            {
                vABParityUnknown
            }
            value if value == AB_PARITY_UNDF as i32 => AB_PARITY_UNDF as i32,
            value
                if value == (AB_PARITY_ODD as i32 | AB_UNKNOWN | AB_NEGATIVE)
                    || value == (AB_PARITY_ODD as i32 | AB_NEGATIVE)
                    || value == AB_PARITY_IISO as i32
                    || value == (AB_PARITY_IISO as i32 | AB_UNKNOWN)
                    || value == (AB_PARITY_UNDF as i32 | AB_NEGATIVE)
                    || value == (AB_PARITY_UNDF as i32 | AB_UNKNOWN | AB_NEGATIVE) =>
            {
                AB_PARITY_NONE as i32
            }
            _ => -1,
        }
    }
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn set_stereo_bonds_parity(
    heap: &mut SourceHeap,
    out_at: SourceMutPointer<sp_ATOM>,
    at: SourceMutPointer<inp_ATOM>,
    at_1: i32,
    at_removed_H: SourceConstPointer<inp_ATOM>,
    num_removed_H: i32,
    nMode: INCHI_MODE,
    q: SourceMutPointer<QUEUE>,
    nAtomLevel: SourceMutPointer<AT_RANK>,
    cSource: SourceMutPointer<S_CHAR>,
    min_sb_ring_size: AT_RANK,
    bPointedEdgeStereo: i32,
    vABParityUnknown: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3009 set_stereo_bonds_parity
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and snapshots add overhead.
    /*
    int set_stereo_bonds_parity( sp_ATOM *out_at,
                                 inp_ATOM *at,
                                 int at_1,
                                 inp_ATOM *at_removed_H,
                                 int num_removed_H,
                                 INCHI_MODE nMode, QUEUE *q,
                                 AT_RANK *nAtomLevel,
                                 S_CHAR *cSource,
                                 AT_RANK min_sb_ring_size,
                                 int bPointedEdgeStereo,
                                 int vABParityUnknown )
    {
        int j, k, i_next_at_1, i_next_at_2, at_2, next_at_2, num_stereo_bonds, bFound, bAllene; /* djb-rwth: removing redundant variables */
        int bond_type, num_2s_1, num_alt_1;
        int num_2s_2, num_alt_2;
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
        int num_wrong_bonds_1, num_wrong_bonds_2;
    #endif
    #if ( N_V_STEREOBONDS == 1 )
        int n2sh, num_2s_hetero[2], num_2s_hetero_next[2], next_next_at, type_N, type_N_next;
    #endif
        int num_stored_isotopic_stereo_bonds; /* djb-rwth: removing redundant variables/code */
        int chain_length, num_chains, cur_chain_length;
        int all_at_2[MAX_NUM_STEREO_BONDS];
        int all_pos_1[MAX_NUM_STEREO_BONDS], all_pos_2[MAX_NUM_STEREO_BONDS];
        S_CHAR all_unkn[MAX_NUM_STEREO_BONDS];
        int /*at_1_parity, at_2_parity,*/ nUnknown, stop = 0;

        /* at_1_parity = AB_PARITY_NONE; */ /*  do not know */

        /*  Check valence */
        if (MAX_NUM_STEREO_BOND_NEIGH < at[at_1].valence + at[at_1].num_H ||
             MIN_NUM_STEREO_BOND_NEIGH > at[at_1].valence + at[at_1].num_H)
        {
            return 0;
        }
        if (!bCanAtomHaveAStereoBond( at[at_1].elname, at[at_1].charge, at[at_1].radical ))
        {
            return 0;
        }
        if (at[at_1].c_point)
        {
            return 0; /* rejects atoms that can lose or gain a (positive) charge. 01-24-2003 */
        }

        /*  middle cumulene atoms, for example, =C=, should be ignored here */
        /*  only atoms at the ends of cumulene chains are considered. */
        if (!at[at_1].num_H && 2 == at[at_1].valence &&
             BOND_DOUBLE == get_allowed_stereo_bond_type( (int) at[at_1].bond_type[0] ) &&
             BOND_DOUBLE == get_allowed_stereo_bond_type( (int) at[at_1].bond_type[1] ))
        {
            return 0;
        }

        /*  count bonds and find the second atom on the stereo bond */
        num_2s_1 = num_alt_1 = 0;
        chain_length = 0;
        num_chains = 0;
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
        num_wrong_bonds_1 = 0;
    #endif
    #if ( N_V_STEREOBONDS == 1 )
        num_2s_hetero[0] = num_2s_hetero[1] = type_N = 0;
        if (0 == at[at_1].num_H && 0 == at[at_1].charge && 0 == at[at_1].radical &&
             3 == get_endpoint_valence( at[at_1].el_number ))
        {
            if (2 == at[at_1].valence && 3 == at[at_1].chem_bonds_valence)
            {
                type_N = 1;
            }
            else
            {
                if (3 == at[at_1].valence && 5 == at[at_1].chem_bonds_valence)
                {
                    type_N = 2; /* unfortunately includes >N# */
                }
            }
        }
    #endif
        for (i_next_at_1 = 0, num_stereo_bonds = 0; i_next_at_1 < at[at_1].valence; i_next_at_1++)
        {
            nUnknown = ( at[at_1].bond_stereo[i_next_at_1] == STEREO_DBLE_EITHER );
            bond_type = get_allowed_stereo_bond_type( (int) at[at_1].bond_type[i_next_at_1] );
            at_2 = -1; /* not found */
            if (bond_type == BOND_ALTERN ||
                 bond_type == BOND_DOUBLE)
            {
                at_2 = at[at_1].neighbor[i_next_at_1]; /* djb-rwth: removing redundant code */
                next_at_2 = at_1;
            }
            switch (bond_type)
            {
                case BOND_ALTERN:
                    num_alt_1++;
    #if ( FIND_RING_SYSTEMS == 1 )
                    if (at[at_1].nRingSystem != at[at_2].nRingSystem)
                    {
                        continue; /* reject alt. bond connecting different ring systems */
                    }
    #endif
                    if (( nMode & CMODE_NO_ALT_SBONDS ) ||
                         !bCanAtomHaveAStereoBond( at[at_2].elname, at[at_2].charge, at[at_2].radical ))
                    {
                        continue; /*  reject non-stereogenic bond to neighbor ord. #i_next_at_1 */
                    }
                    break;
                case BOND_DOUBLE:
                    /*  check for cumulene/allene */
                    num_2s_1++;
                    cur_chain_length = 0;
                    if (bCanAtomBeTerminalAllene( at[at_1].elname, at[at_1].charge, at[at_1].radical ))
                    {
                        /*
                         * Example of cumulene
                         * chain length = 2:     >X=C=C=Y<
                         *                        | | | |
                         *  1st cumulene atom= at_1 | | at_2 =last cumlene chain atom
                         *  next to at_1=   next_at_1 next_at_2  =previous to at_2
                         *
                         *  chain length odd:  stereocenter on the middle atom ( 1=> allene )
                         *  chain length even: "long stereogenic bond"
                         */
                        while (( bAllene =
                            !at[at_2].num_H && at[at_2].valence == 2 &&
                            BOND_DOUBLE == get_allowed_stereo_bond_type( (int) at[at_2].bond_type[0] ) &&
                            BOND_DOUBLE == get_allowed_stereo_bond_type( (int) at[at_2].bond_type[1] ) ) &&
                                bCanAtomBeMiddleAllene( at[at_2].elname, at[at_2].charge, at[at_2].radical ))
                        {
                            k = ( (int) at[at_2].neighbor[0] == next_at_2 ); /*  opposite neighbor position */
                            next_at_2 = at_2;
                            nUnknown += ( at[at_2].bond_stereo[k] == STEREO_DBLE_EITHER );
                            at_2 = (int) at[at_2].neighbor[k];
                            cur_chain_length++;  /*  count =C= atoms */
                        }
                        if (cur_chain_length)
                        {
                            num_chains++;
                            if (bAllene /* at the end of the chain atom Y is =Y=, not =Y< or =Y- */ ||
                                 !bCanAtomBeTerminalAllene( at[at_2].elname, at[at_2].charge, at[at_2].radical ))
                            {
                                cur_chain_length = 0; /* djb-rwth: ignoring LLVM warning: value used */
                                continue; /*  ignore: does not fit cumulene description; go to check next at_1 neighbor */
                            }
                            chain_length = cur_chain_length; /*  accept a stereogenic cumulele */
                        }
                    }
    #if ( N_V_STEREOBONDS == 1 )
                    if (!cur_chain_length &&
                         0 <= ( n2sh = bIsSuitableHeteroInpAtom( at + at_2 ) ))
                    {
                        num_2s_hetero[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
                    }
    #endif
                    if (!cur_chain_length &&
                         !bCanAtomHaveAStereoBond( at[at_2].elname, at[at_2].charge, at[at_2].radical ))
                    {
                        continue; /*  reject non-stereogenic bond to neighbor #i_next_at_1 */
                    }

                    break;

                case BOND_SINGLE:
                case BOND_TAUTOM:
                    continue; /*  reject non-stereogenic bond to neighbor #i_next_at_1 */
                default:
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
                    num_wrong_bonds_1++;
                    continue;
    #else
                    return 0; /*  wrong bond type; */
    #endif
            }

            /*  Check atom at the opposite end of possibly stereogenic bond */

            bFound = ( at_2 >= 0 && at_1 > at_2 ); /*  i_next_at_1 = at_1 stereogenic bond neighbor attachment number */

            if (bFound)
            {
                /*  Check "at_2" atom on the opposite side of the bond or cumulene chain */
                if (MAX_NUM_STEREO_BOND_NEIGH < at[at_2].valence + at[at_2].num_H ||
                     MIN_NUM_STEREO_BOND_NEIGH > at[at_2].valence + at[at_2].num_H)
                    continue;

                /*  Check at_2 neighbors and bonds */
                num_2s_2 = num_alt_2 = 0;
    #if ( N_V_STEREOBONDS == 1 )
                num_2s_hetero_next[0] = num_2s_hetero_next[1] = type_N_next = 0;
                if (0 == at[at_2].num_H && 0 == at[at_2].charge && 0 == at[at_2].radical &&
                     3 == get_endpoint_valence( at[at_2].el_number ))
                {
                    if (2 == at[at_2].valence && 3 == at[at_2].chem_bonds_valence)
                    {
                        type_N_next = 1; /* -N= */
                    }
                    else
                    {
                        if (3 == at[at_2].valence && 5 == at[at_2].chem_bonds_valence)
                        {
                            type_N_next = 2; /* unfortunately includes >N# */
                        }
                    }
                }
    #endif
                i_next_at_2 = -1;  /*  unassigned mark */
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
                num_wrong_bonds_2 = 0;
    #endif
                for (j = 0; j < at[at_2].valence; j++)
                {
                    bond_type = get_allowed_stereo_bond_type( (int) at[at_2].bond_type[j] );
                    if (!bond_type)
                    {
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
                        num_wrong_bonds_2++;
                        continue;  /*  this bond type is not allowed to be adjacent to a stereo bond */
    #else
                        break;
    #endif
                    }
                    if (bond_type == BOND_DOUBLE)
                    {
                        num_2s_2++;
    #if ( N_V_STEREOBONDS == 1 )
                        next_next_at = at[at_2].neighbor[j];
                        if (0 <= ( n2sh = bIsSuitableHeteroInpAtom( at + next_next_at ) ))
                        {
                            num_2s_hetero_next[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
                        }
    #endif
                    }
                    else
                    {
                        num_alt_2 += ( bond_type == BOND_ALTERN );
                    }
                    if ((int) at[at_2].neighbor[j] == next_at_2)
                    {
                        i_next_at_2 = j; /*  assigned */
                    }
                }
                if (
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
                     num_wrong_bonds_2 > 1 || (num_wrong_bonds_2 && 2 >= at[at_2].valence) || /* djb-rwth: addressing LLVM warning */
    #else
                     j < at[at_2].valence /* "next" has a wrong bond type*/ ||
    #endif
    ( num_alt_2 > 0 ) + ( num_2s_2 > 0 ) != 1 || /* all double XOR all alt bonds only */
     /* num_2s_2 > 1  ||*/ /* only one double bond permitted */
                      i_next_at_2 < 0 /* atom next to the opposite atom not found */)
                {
                    bFound = 0;
                }
                else
                    if (at[at_2].c_point)
                    {
                        bFound = 0; /* rejects atoms that can lose or gain a (positive) charge. 01-24-2003 */
                    }
                    else
                        if (num_2s_2 > 2)
                        {
                            bFound = 0;
                        }
                        else
    #if ( N_V_STEREOBONDS == 1 )
                            if ( 3==( type_N | type_N_next ) && ( (2==type_N && !bIsOxide( at, at_1 )) ||
                                 (2==type_N_next && !bIsOxide( at, at_2 )) )) /* djb-rwth: addressing LLVM warnings */
                            {
                                bFound = 0;
                            }
                            else
    #endif
                            {
                                if (2 == num_2s_2)
                                {
    #if ( N_V_STEREOBONDS == 1 )
                                    if (!chain_length &&
                                         1 == ( num_2s_hetero_next[0] | num_2s_hetero_next[1] ) &&
                                         3 == at[at_2].valence + at[at_2].num_H &&
                                         5 == at[at_2].chem_bonds_valence + at[at_2].num_H &&
                                         3 == get_endpoint_valence( at[at_2].el_number ) &&
                                         ( !type_N || bIsOxide( at, at_2 ) ))
                                    {
                                        /*
                                         *   found:
                                         *
                                         *    \      /    \      /    \      /
                                         *     \    /      \    /      \    /
                                         *      N==C   or   N==C   or   N==N
                                         *    //    \     //    \     //    \
                                         *   O  ^    \   N  ^    \   O  ^    \
                                         *      |           |           |
                                         *      |           |           |
                                         *      at[at_2]    at[at_2]    at[at_2]
                                         */
                                        ;
                                    }
                                    else
                                    {
                                        bFound = 0;
                                    }
    #else
                                    bFound = 0;
    #endif
                                    }
                                }

                if (chain_length && num_alt_2)
                {
                    return 0; /*  allow no alt bonds in cumulenes */
                }
            }

            if (bFound)
            {
                all_pos_1[num_stereo_bonds] = i_next_at_1; /* neighbor to at_1 position */
                all_pos_2[num_stereo_bonds] = i_next_at_2; /* neighbor to at_2 position */
                all_at_2[num_stereo_bonds] = at_2;        /* at_2 */
                all_unkn[num_stereo_bonds] = nUnknown;    /* stereogenic bond has Unknown configuration */
                /*
                if ( (at[at_1].bUsed0DParity & 2) || (at[at_2].bUsed0DParity & 2) ) {
                    for ( k = 0; k < MAX_NUM_STEREO_BONDS && at[at_1].sb_parity[k]; k ++ ) {
                        if ( at[at_1].sb_neigh[k] == i_next_at_1 ) {
                            if ( at[at_1].sb_parity[k] == AB_PARITY_UNKN && !nUnknown ) {
                                all_unkn[num_stereo_bonds] = 1;
                            }
                            break;
                        }
                    }
                }
                */
                num_stereo_bonds++;
            }
        }
        if (num_chains > 1)
        {
            return 0; /*  cannot be more than 1 cumulene chain. */
        }
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
        if (num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 && 2 >= at[at_1].valence)) /* djb-rwth: addressing LLVM warning */
        {
            return 0; /* wrong bond type */
        }
    #endif
        /*  Accept only short chains for now */
        /*  chain_length=1: >C=C=C<      tetrahedral center, allene */
        /*  chain_length=2: >C=C=C=C<    stereogenic bond, cumulene */
        if (chain_length && ( num_stereo_bonds != 1 || num_alt_1 || chain_length > MAX_CUMULENE_LEN ))
        {
            return 0;
        }

        /*  We need 1 double bond/chain XOR up to 3 arom. bonds */
        /*  to have a stereogenic bond */
        if (( num_alt_1 > 0 ) + ( num_2s_1 > 0 ) != 1 || !num_stereo_bonds /*|| num_2s_1 > 1*/)
            return 0;

        if (num_2s_1 > 1)
        {
    #if ( N_V_STEREOBONDS == 1 )
            if (2 == num_2s_1 &&
                 2 == type_N &&
                 1 == ( num_2s_hetero[0] | num_2s_hetero[1] ) &&
                 3 == at[at_1].valence + at[at_1].num_H &&
                 5 == at[at_1].chem_bonds_valence + at[at_1].num_H &&
                 3 == get_endpoint_valence( at[at_1].el_number ))
            {
                ;
            }
            else
            {
                return 0;
            }
    #else
            return 0;
    #endif
        }

        /* ================== Calculate parities ====================== */


        /*  Find possibly stereo bonds and save them */
        num_stored_isotopic_stereo_bonds = 0;
        /* djb-rwth: removing redundant code */
        for (k = 0; k < num_stereo_bonds; k++)
        {

            int cur_parity, next_parity, abs_cur_parity, abs_next_parity, dot_prod_z;
            S_CHAR z_dir1[3], z_dir2[3]; /*  3D vectors for half stereo bond parity direction */
            int  chain_len_bits = MAKE_BITS_CUMULENE_LEN( chain_length );
            int  cur_parity_defined, next_parity_defined;
            int  cur_action, next_action, result_action;

            at_2 = all_at_2[k];
            i_next_at_1 = all_pos_1[k];

    #if ( MIN_SB_RING_SIZE > 0 )
            if (at[at_1].nRingSystem == at[at_2].nRingSystem)
            {
                /*  check min. ring size only if both double bond/cumulene */
                /*  ending atoms belong to the same ring system */
                j = is_bond_in_Nmax_memb_ring( at, at_1, i_next_at_1, q, nAtomLevel, cSource, min_sb_ring_size );
                if (j > 0)
                {
                    continue;
                }
                else
                {
                    if (j < 0)
                    {
                        return CT_STEREOBOND_ERROR;
                    }
                }
            }
    #endif

            i_next_at_2 = all_pos_2[k];
            nUnknown = all_unkn[k];
            memset( z_dir1, 0, sizeof( z_dir1 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( z_dir2, 0, sizeof( z_dir2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */

            /********************************************************************************
             * Find atom parities (negative means parity due to H-isotopes only)
             * and half stereo bond parity directions z_dir1, z_dir2.
             *
             * Bond can have unknown or undefined parity or no parity because of:
             * 1. Geometry (poorly defined, cannot calculate, for example linear =C-F
             *    or =CHD with no geometry) -- Undefined parity
             *                                                              H
             * 2. Identical H atoms (no parity in principle, for example =C<  )
             *    -- No parity                                              H
             *
             * 3. The user said double bond stereo is unknown
             *    or at least one of single bonds is in unknown direction
             *    -- Unknown parity
             *
             * These 3 cases (see above) are referred below as 1, 2, 3.
             * Each of the cases may be present or not (2 possibilities)
             * Total number of combination is 2*2*2=8
             *
             * Since a case when all 3 are not present is a well-defined parity,
             * we do not consider this case here. Then 2*2*2-1=7 cases are left.
             *
             * If several cases are present, list them below separated by "+".
             * For example, 1+2 means (1) undefined geometry and (2) no parity
             * is possible because of identical H atoms.
             *
             * N) Decision table, Non-isotopic, 2*2*2-1=7 cases:
             * =================================================
             * none     : 2+any: 1+2(e.g.=CH2); 1+2+3; 2; 2+3  AB_PARITY_NONE=0
             * undefined: 1                                    AB_PARITY_UNDF
             * unknown  : 1+3; 3                               AB_PARITY_UNKN
             *
             * I) Decision table, Isotopic, 2*2*2-1=7 cases:
             * =============================================
             * none     : none
             * undefined: 1; 1+2; 1+2+3; 2; 2+3
             * unknown  : 1+3; 3
             *
             * Note: When defining identical atoms H atoms in case 2,
             *       Isotopic and Non-isotopic cases are different:
             *  N: do NOT take into account the isotopic composition of H atoms
             *  I: DO take into account the isotopic composition of H atoms
             *     (it is assumed that H isotopes are always different)
             *
             * half_stereo_bond_parity() returns:
             * ==================================
             * Note: half_stereo_bond_parity() is unaware of case 3.
             *
             * can't be a half of a stereo bond    AB_PARITY_NONE
             * 1, isotopic & non-isotopic:         AB_PARITY_UNDF
             * 1, isotopic only                   -AB_PARITY_UNDF
             * 2, no parity: identical H isotopes  AB_PARITY_IISO
             * 3, 'Either' single bond(s)          AB_PARITY_UNKN ???
             * 3, 'Either' single bond(s), iso H  -AB_PARITY_UNKN ???
             * defined parity                      AB_PARITY_ODD,  AB_PARITY_EVEN
             * defined parity for isotopic only:  -AB_PARITY_ODD, -AB_PARITY_EVEN
             *
             * Resultant value for the stereo bond parity
             * ---+-------------------+-------+--------+----------------+
             * 3? | half_stereo_bond_ | N or I| case 1,| bond parity    |
             *    |  parity()=        |       | 2 or 3 |                |
             * ---+-------------------+-------+--------+----------------+
             *   ( AB_PARITY_ODD/EVEN) => N&I: -       => AB_PARITY_CALC (=6, calc.later)
             * 3+( AB_PARITY_ODD/EVEN) => N&I: 3       => AB_PARITY_UNKN (=3)
             *   (-AB_PARITY_ODD/EVEN) => N:   2       => AB_PARITY_NONE (=0)
             *   (-AB_PARITY_ODD/EVEN) => I:   -       => AB_PARITY_CALC
             * 3+(-AB_PARITY_ODD/EVEN) => N:   2+3     => AB_PARITY_UNDF (=4)
             * 3+(-AB_PARITY_ODD/EVEN) => I:   3       => AB_PARITY_UNKN
             *   ( AB_PARITY_IISO )    => N:   1+2, 2  => AB_PARITY_NONE (=0)
             *   ( AB_PARITY_IISO )    => I:   1+2, 2  => AB_PARITY_UNDF
             * 3+( AB_PARITY_IISO )    => N:  1+2+3,2+3=> AB_PARITY_NONE
             * 3+( AB_PARITY_IISO )    => I:  1+2+3,2+3=> AB_PARITY_UNDF
             *   ( AB_PARITY_UNDF )    => N&I: 1       => AB_PARITY_UNDF
             * 3+( AB_PARITY_UNDF )    => N&I: 1+3     => AB_PARITY_UNKN
             *   (-AB_PARITY_UNDF )    => N:   1+2     => AB_PARITY_NONE
             *   (-AB_PARITY_UNDF )    => I:   1       => AB_PARITY_UNDF
             * 3+(-AB_PARITY_UNDF )    => N:   1+2+3   => AB_PARITY_NONE
             * 3+(-AB_PARITY_UNDF )    => I:   1+3     => AB_PARITY_UNKN
             * ---+-------------------+-------+--------+----------------+

             * If bond parity is undefined because abs(dot_prod_z) < MIN_DOT_PROD
             * then replace: AB_PARITY_CALC
             *         with: AB_PARITY_UNDF
             * Joining two half_bond_parity() results:
             *
             *
             * atom1 \ atom2   | AB_PARITY_NONE  AB_PARITY_UNKN  AB_PARITY_UNDF  AB_PARITY_CALC
             * ----------------+---------------------------------------------------------------
             *0=AB_PARITY_NONE | AB_PARITY_NONE  AB_PARITY_NONE  AB_PARITY_NONE  AB_PARITY_NONE
             *3=AB_PARITY_UNKN |                 AB_PARITY_UNKN  AB_PARITY_UNKN  AB_PARITY_UNKN
             *4=AB_PARITY_UNDF |                                 AB_PARITY_UNDF  AB_PARITY_UNDF
             *6=AB_PARITY_CALC |                                                 AB_PARITY_CALC
             *
             * that is, take min out of the two
             *********************************************************************************/

            cur_parity = half_stereo_bond_parity( at, at_1, at_removed_H, num_removed_H,
                                                    z_dir1, bPointedEdgeStereo, vABParityUnknown );
            next_parity = half_stereo_bond_parity( at, at_2, at_removed_H, num_removed_H,
                                                    z_dir2, bPointedEdgeStereo, vABParityUnknown );

            if (RETURNED_ERROR( cur_parity ) || RETURNED_ERROR( next_parity ))
            {
                return CT_CALC_STEREO_ERR;
            }
            if (( at[at_1].bUsed0DParity & FlagSB_0D ) || ( at[at_2].bUsed0DParity & FlagSB_0D )) /* djb-rwth: condition corrected */
            {
                FixSb0DParities( at, /* at_removed_H, num_removed_H,*/ chain_length,
                                 at_1, i_next_at_1, z_dir1,
                                 at_2, i_next_at_2, z_dir2, &cur_parity, &next_parity );
            }

            if (cur_parity == AB_PARITY_NONE || abs( cur_parity ) == AB_PARITY_IISO)
            {
                continue;
            }
            if (next_parity == AB_PARITY_NONE || abs( next_parity ) == AB_PARITY_IISO)
            {
                continue;
            }

            cur_action = half_stereo_bond_action( cur_parity, nUnknown, 0, vABParityUnknown ); /*  -1 => program error */
            next_action = half_stereo_bond_action( next_parity, nUnknown, 0, vABParityUnknown );
            result_action = inchi_min( cur_action, next_action );

            if (result_action == -1)
            {
                stop = 1; /*  program error <BRKPT> */
            }

            abs_cur_parity = abs( cur_parity );
            abs_next_parity = abs( next_parity );
            cur_parity_defined = ATOM_PARITY_WELL_DEF( abs_cur_parity );
            next_parity_defined = ATOM_PARITY_WELL_DEF( abs_next_parity );

            if (cur_parity_defined && next_parity_defined)
            {
                /*  find how the whole bond parity depend on geometry */
                /*  if dot_prod_z < 0 then bond_parity := 3-bond_parity */
                /*  can be done only for a well-defined geometry */
                /*
                dot_prod_z  = (chain_len_bits & BIT_CUMULENE_CHI)?
                               triple_prod_char( at, at_1, i_next_at_1, z_dir1, at_2, i_next_at_2, z_dir2 ) :
                               dot_prodchar3(z_dir1, z_dir2);
                */

                dot_prod_z = ( chain_len_bits && BOND_CHAIN_LEN( chain_len_bits ) % 2 )
                    ?  triple_prod_char( at, at_1, i_next_at_1, z_dir1, at_2, i_next_at_2, z_dir2 )
                    :  dot_prodchar3( z_dir1, z_dir2 );

                if (abs( dot_prod_z ) < MIN_DOT_PROD)
                {
                    /*  The geometry is not well-defined. Eliminate AB_PARITY_CALC */
                    result_action = inchi_min( result_action, AB_PARITY_UNDF );
                }
            }
            else
            {
                dot_prod_z = 0;
            }

            if (result_action != AB_PARITY_NONE && result_action != -1)
            {
                /*  Stereo, no isotopes (only positive) */
                if (cur_parity > 0 && next_parity > 0)
                {
                    if (save_a_stereo_bond( dot_prod_z, result_action | chain_len_bits,
                        at_1, i_next_at_1, out_at[at_1].stereo_bond_neighbor,
                        out_at[at_1].stereo_bond_ord, out_at[at_1].stereo_bond_z_prod,
                        out_at[at_1].stereo_bond_parity,
                        at_2, i_next_at_2, out_at[at_2].stereo_bond_neighbor,
                        out_at[at_2].stereo_bond_ord, out_at[at_2].stereo_bond_z_prod,
                        out_at[at_2].stereo_bond_parity ))
                    {
                        if (!out_at[at_1].parity ||
                             (cur_parity_defined && !ATOM_PARITY_WELL_DEF( abs( out_at[at_1].parity )) )) /* djb-rwth: addressing LLVM warning */
                        {
                            out_at[at_1].parity = cur_parity;
                            memcpy(out_at[at_1].z_dir, z_dir1, sizeof(out_at[0].z_dir));
                        }
                        if (!out_at[at_2].parity ||
                             (next_parity_defined && !ATOM_PARITY_WELL_DEF( abs( out_at[at_2].parity )) )) /* djb-rwth: addressing LLVM warning */
                        {
                            out_at[at_2].parity = next_parity;
                            memcpy(out_at[at_2].z_dir, z_dir2, sizeof(out_at[0].z_dir));
                        }
                        out_at[at_1].bAmbiguousStereo |= at[at_1].bAmbiguousStereo;
                        out_at[at_2].bAmbiguousStereo |= at[at_2].bAmbiguousStereo;
                        /* djb-rwth: removing redundant code */
                    }
                }
            }

            /*  Stereo + isotopic (all non-zero) */
            cur_action = half_stereo_bond_action( cur_parity, nUnknown, 1, vABParityUnknown ); /*  -1 => program error */
            next_action = half_stereo_bond_action( next_parity, nUnknown, 1, vABParityUnknown );
            result_action = inchi_min( cur_action, next_action );
            cur_parity = abs_cur_parity;
            next_parity = abs_next_parity;
            if (result_action != AB_PARITY_NONE && result_action != -1)
            {
                /*  Stereo, isotopic */
                if (cur_parity > 0 && next_parity > 0)
                {
                    if (save_a_stereo_bond( dot_prod_z, result_action | chain_len_bits,
                        at_1, i_next_at_1, out_at[at_1].stereo_bond_neighbor2,
                        out_at[at_1].stereo_bond_ord2, out_at[at_1].stereo_bond_z_prod2,
                        out_at[at_1].stereo_bond_parity2,
                        at_2, i_next_at_2, out_at[at_2].stereo_bond_neighbor2,
                        out_at[at_2].stereo_bond_ord2, out_at[at_2].stereo_bond_z_prod2,
                        out_at[at_2].stereo_bond_parity2 ))
                    {
                        if (!out_at[at_1].parity2 ||
                             (cur_parity_defined && !ATOM_PARITY_WELL_DEF( abs( out_at[at_1].parity2 )) )) /* djb-rwth: addressing LLVM warning */
                        {
                            out_at[at_1].parity2 = cur_parity /*| chain_len_bits*/;
                            if (!out_at[at_1].parity)
                            {
                                memcpy(out_at[at_1].z_dir, z_dir1, sizeof(out_at[0].z_dir));
                            }
                        }
                        if (!out_at[at_2].parity2 || /* next line changed from abs(out_at[at_2].parity) 2006-03-05 */
                             (next_parity_defined && !ATOM_PARITY_WELL_DEF( abs( out_at[at_2].parity2 )) )) /* djb-rwth: addressing LLVM warning */
                        {
                            out_at[at_2].parity2 = next_parity /*| chain_len_bits*/;
                            if (!out_at[at_2].parity)
                            {
                                memcpy(out_at[at_2].z_dir, z_dir2, sizeof(out_at[0].z_dir));
                            }
                        }
                        out_at[at_1].bAmbiguousStereo |= at[at_1].bAmbiguousStereo;
                        out_at[at_2].bAmbiguousStereo |= at[at_2].bAmbiguousStereo;
                        num_stored_isotopic_stereo_bonds++;
                    }
                }
            }
            else
            {
                if (result_action == -1)
                {
                    stop = 1; /*  program error? <BRKPT> */
                }
            }
        }

        if (stop)
        {
            return CT_CALC_STEREO_ERR;
        }

        return /*num_stored_stereo_bonds+*/ num_stored_isotopic_stereo_bonds;
    }
    */
    // END INCHI C FUNCTION: set_stereo_bonds_parity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: set_stereo_bonds_parity
    // INCHI✔️❌: ONE_BAD_SB_NEIGHBOR == 1, N_V_STEREOBONDS == 1, FIND_RING_SYSTEMS == 1.
    // INCHI✔️❌: MIN_SB_RING_SIZE == 8 enables is_bond_in_Nmax_memb_ring.
    // INCHI✔️❌: MAX_CUMULENE_LEN == 2, MULT_STEREOBOND == 8, MIN_DOT_PROD == 50.
    // END INCHI ACTIVE MACRO CONFIGURATION: set_stereo_bonds_parity

    let at_1_index = usize::try_from(at_1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(at.as_const())?
        .get(at_1_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let neighbor_count = i32::from(current.valence) + i32::from(current.num_H);
    if neighbor_count > MAX_NUM_STEREO_BOND_NEIGH as i32
        || neighbor_count < MIN_NUM_STEREO_BOND_NEIGH as i32
    {
        return Ok(0);
    }
    if bCanAtomHaveAStereoBond(&current.elname, current.charge, current.radical)? == 0
        || current.c_point != 0
    {
        return Ok(0);
    }
    if current.num_H == 0
        && current.valence == 2
        && get_allowed_stereo_bond_type(i32::from(current.bond_type[0])) == BOND_DOUBLE as i32
        && get_allowed_stereo_bond_type(i32::from(current.bond_type[1])) == BOND_DOUBLE as i32
    {
        return Ok(0);
    }

    let mut num_2s_1 = 0_i32;
    let mut num_alt_1 = 0_i32;
    let mut chain_length = 0_i32;
    let mut num_chains = 0_i32;
    let mut num_wrong_bonds_1 = 0_i32;
    let mut num_2s_hetero = [0_i32; 2];
    let mut type_n = 0_i32;
    if current.num_H == 0
        && current.charge == 0
        && current.radical == 0
        && get_endpoint_valence(current.el_number) == 3
    {
        if current.valence == 2 && current.chem_bonds_valence == 3 {
            type_n = 1;
        } else if current.valence == 3 && current.chem_bonds_valence == 5 {
            type_n = 2;
        }
    }

    let mut all_at_2 = [0_i32; MAX_NUM_STEREO_BONDS as usize];
    let mut all_pos_1 = [0_i32; MAX_NUM_STEREO_BONDS as usize];
    let mut all_pos_2 = [0_i32; MAX_NUM_STEREO_BONDS as usize];
    let mut all_unkn = [0_i8; MAX_NUM_STEREO_BONDS as usize];
    let mut num_stereo_bonds = 0_i32;

    for i_next_at_1 in 0..i32::from(current.valence) {
        let first_order =
            usize::try_from(i_next_at_1).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let mut n_unknown = i32::from(current.bond_stereo[first_order] == STEREO_DBLE_EITHER as i8);
        let bond_type = get_allowed_stereo_bond_type(i32::from(current.bond_type[first_order]));
        let mut at_2 = -1_i32;
        let mut next_at_2 = 0_i32;
        if bond_type == BOND_ALTERN as i32 || bond_type == BOND_DOUBLE as i32 {
            at_2 = i32::from(current.neighbor[first_order]);
            next_at_2 = at_1;
        }
        let mut cur_chain_length = 0_i32;
        let mut reject_current_chain = false;
        match bond_type as u32 {
            BOND_ALTERN => {
                num_alt_1 += 1;
                let opposite = heap
                    .slice(at.as_const())?
                    .get(usize::try_from(at_2).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                if current.nRingSystem != opposite.nRingSystem
                    || nMode & CMODE_NO_ALT_SBONDS as INCHI_MODE != 0
                    || bCanAtomHaveAStereoBond(&opposite.elname, opposite.charge, opposite.radical)?
                        == 0
                {
                    continue;
                }
            }
            BOND_DOUBLE => {
                num_2s_1 += 1;
                if bCanAtomBeTerminalAllene(&current.elname, current.charge, current.radical)? != 0
                {
                    loop {
                        let opposite_index = usize::try_from(at_2)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                        let opposite = heap
                            .slice(at.as_const())?
                            .get(opposite_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .clone();
                        let is_allene = opposite.num_H == 0
                            && opposite.valence == 2
                            && get_allowed_stereo_bond_type(i32::from(opposite.bond_type[0]))
                                == BOND_DOUBLE as i32
                            && get_allowed_stereo_bond_type(i32::from(opposite.bond_type[1]))
                                == BOND_DOUBLE as i32;
                        if !is_allene
                            || bCanAtomBeMiddleAllene(
                                &opposite.elname,
                                opposite.charge,
                                opposite.radical,
                            )? == 0
                        {
                            if cur_chain_length != 0 {
                                num_chains += 1;
                                if is_allene
                                    || bCanAtomBeTerminalAllene(
                                        &opposite.elname,
                                        opposite.charge,
                                        opposite.radical,
                                    )? == 0
                                {
                                    cur_chain_length = 0;
                                    reject_current_chain = true;
                                } else {
                                    chain_length = cur_chain_length;
                                }
                            }
                            break;
                        }
                        let k = usize::from(i32::from(opposite.neighbor[0]) == next_at_2);
                        next_at_2 = at_2;
                        n_unknown += i32::from(opposite.bond_stereo[k] == STEREO_DBLE_EITHER as i8);
                        at_2 = i32::from(opposite.neighbor[k]);
                        cur_chain_length += 1;
                    }
                    if reject_current_chain {
                        continue;
                    }
                }
                if cur_chain_length == 0 {
                    let suitable =
                        bIsSuitableHeteroInpAtom(heap, at.as_const().offset(i64::from(at_2))?)?;
                    if suitable >= 0 {
                        num_2s_hetero[usize::try_from(suitable)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?] += 1;
                    }
                    let opposite = heap
                        .slice(at.as_const())?
                        .get(
                            usize::try_from(at_2)
                                .map_err(|_| SourceHeapError::PointerOutOfBounds)?,
                        )
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    if bCanAtomHaveAStereoBond(&opposite.elname, opposite.charge, opposite.radical)?
                        == 0
                    {
                        continue;
                    }
                }
            }
            BOND_SINGLE | BOND_TAUTOM => continue,
            _ => {
                num_wrong_bonds_1 += 1;
                continue;
            }
        }

        let mut found = at_2 >= 0 && at_1 > at_2;
        let mut i_next_at_2 = -1_i32;
        if found {
            let opposite_index =
                usize::try_from(at_2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let opposite = heap
                .slice(at.as_const())?
                .get(opposite_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let opposite_neighbor_count = i32::from(opposite.valence) + i32::from(opposite.num_H);
            if opposite_neighbor_count > MAX_NUM_STEREO_BOND_NEIGH as i32
                || opposite_neighbor_count < MIN_NUM_STEREO_BOND_NEIGH as i32
            {
                continue;
            }

            let mut num_2s_2 = 0_i32;
            let mut num_alt_2 = 0_i32;
            let mut num_2s_hetero_next = [0_i32; 2];
            let mut type_n_next = 0_i32;
            if opposite.num_H == 0
                && opposite.charge == 0
                && opposite.radical == 0
                && get_endpoint_valence(opposite.el_number) == 3
            {
                if opposite.valence == 2 && opposite.chem_bonds_valence == 3 {
                    type_n_next = 1;
                } else if opposite.valence == 3 && opposite.chem_bonds_valence == 5 {
                    type_n_next = 2;
                }
            }
            let mut num_wrong_bonds_2 = 0_i32;
            for j in 0..i32::from(opposite.valence) {
                let order = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let opposite_bond_type =
                    get_allowed_stereo_bond_type(i32::from(opposite.bond_type[order]));
                if opposite_bond_type == 0 {
                    num_wrong_bonds_2 += 1;
                    continue;
                }
                if opposite_bond_type == BOND_DOUBLE as i32 {
                    num_2s_2 += 1;
                    let next_next_at = i32::from(opposite.neighbor[order]);
                    let suitable = bIsSuitableHeteroInpAtom(
                        heap,
                        at.as_const().offset(i64::from(next_next_at))?,
                    )?;
                    if suitable >= 0 {
                        num_2s_hetero_next[usize::try_from(suitable)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?] += 1;
                    }
                } else {
                    num_alt_2 += i32::from(opposite_bond_type == BOND_ALTERN as i32);
                }
                if i32::from(opposite.neighbor[order]) == next_at_2 {
                    i_next_at_2 = j;
                }
            }
            if num_wrong_bonds_2 > 1
                || (num_wrong_bonds_2 != 0 && opposite.valence <= 2)
                || i32::from(num_alt_2 > 0) + i32::from(num_2s_2 > 0) != 1
                || i_next_at_2 < 0
                || opposite.c_point != 0
                || num_2s_2 > 2
            {
                found = false;
            } else if (type_n | type_n_next) == 3
                && ((type_n == 2 && bIsOxide(heap, at, at_1)? == 0)
                    || (type_n_next == 2 && bIsOxide(heap, at, at_2)? == 0))
            {
                found = false;
            } else if num_2s_2 == 2
                && !(chain_length == 0
                    && (num_2s_hetero_next[0] | num_2s_hetero_next[1]) == 1
                    && i32::from(opposite.valence) + i32::from(opposite.num_H) == 3
                    && i32::from(opposite.chem_bonds_valence) + i32::from(opposite.num_H) == 5
                    && get_endpoint_valence(opposite.el_number) == 3
                    && (type_n == 0 || bIsOxide(heap, at, at_2)? != 0))
            {
                found = false;
            }
            if chain_length != 0 && num_alt_2 != 0 {
                return Ok(0);
            }
        }

        if found {
            let slot = usize::try_from(num_stereo_bonds)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            all_pos_1[slot] = i_next_at_1;
            all_pos_2[slot] = i_next_at_2;
            all_at_2[slot] = at_2;
            all_unkn[slot] = n_unknown as S_CHAR;
            num_stereo_bonds += 1;
        }
    }

    if num_chains > 1 || num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 != 0 && current.valence <= 2) {
        return Ok(0);
    }
    if chain_length != 0
        && (num_stereo_bonds != 1 || num_alt_1 != 0 || chain_length > MAX_CUMULENE_LEN as i32)
    {
        return Ok(0);
    }
    if i32::from(num_alt_1 > 0) + i32::from(num_2s_1 > 0) != 1 || num_stereo_bonds == 0 {
        return Ok(0);
    }
    if num_2s_1 > 1
        && !(num_2s_1 == 2
            && type_n == 2
            && (num_2s_hetero[0] | num_2s_hetero[1]) == 1
            && i32::from(current.valence) + i32::from(current.num_H) == 3
            && i32::from(current.chem_bonds_valence) + i32::from(current.num_H) == 5
            && get_endpoint_valence(current.el_number) == 3)
    {
        return Ok(0);
    }

    let mut num_stored_isotopic_stereo_bonds = 0_i32;
    let mut stop = false;
    for k in
        0..usize::try_from(num_stereo_bonds).map_err(|_| SourceHeapError::PointerOutOfBounds)?
    {
        let at_2 = all_at_2[k];
        let i_next_at_1 = all_pos_1[k];
        let opposite_index =
            usize::try_from(at_2).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let ring_system_1 = heap.slice(at.as_const())?[at_1_index].nRingSystem;
        let ring_system_2 = heap
            .slice(at.as_const())?
            .get(opposite_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?
            .nRingSystem;
        if ring_system_1 == ring_system_2 {
            let ring = is_bond_in_Nmax_memb_ring(
                heap,
                at,
                at_1,
                i_next_at_1,
                q,
                nAtomLevel,
                cSource,
                min_sb_ring_size,
            )?;
            if ring > 0 {
                continue;
            }
            if ring < 0 {
                return Ok(CT_STEREOBOND_ERROR);
            }
        }

        let i_next_at_2 = all_pos_2[k];
        let n_unknown = i32::from(all_unkn[k]);
        let mut z_dir1 = [0_i8; 3];
        let mut z_dir2 = [0_i8; 3];
        let mut cur_parity = half_stereo_bond_parity(
            heap,
            at,
            at_1,
            at_removed_H,
            num_removed_H,
            Some(&mut z_dir1),
            bPointedEdgeStereo,
            vABParityUnknown,
        )?;
        let mut next_parity = half_stereo_bond_parity(
            heap,
            at,
            at_2,
            at_removed_H,
            num_removed_H,
            Some(&mut z_dir2),
            bPointedEdgeStereo,
            vABParityUnknown,
        )?;
        if (CT_ERR_MIN..=CT_ERR_MAX).contains(&cur_parity)
            || (CT_ERR_MIN..=CT_ERR_MAX).contains(&next_parity)
        {
            return Ok(CT_CALC_STEREO_ERR);
        }
        let used_0d_1 =
            i32::from(heap.slice(at.as_const())?[at_1_index].bUsed0DParity) & FlagSB_0D as i32;
        let used_0d_2 = i32::from(
            heap.slice(at.as_const())?
                .get(opposite_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .bUsed0DParity,
        ) & FlagSB_0D as i32;
        if used_0d_1 != 0 || used_0d_2 != 0 {
            let _ = FixSb0DParities(
                heap,
                at.as_const(),
                chain_length,
                at_1,
                i_next_at_1,
                &mut z_dir1,
                at_2,
                i_next_at_2,
                &mut z_dir2,
                &mut cur_parity,
                &mut next_parity,
            )?;
        }
        if cur_parity == AB_PARITY_NONE as i32
            || cur_parity.wrapping_abs() == AB_PARITY_IISO as i32
            || next_parity == AB_PARITY_NONE as i32
            || next_parity.wrapping_abs() == AB_PARITY_IISO as i32
        {
            continue;
        }

        let chain_len_bits = chain_length.wrapping_mul(MULT_STEREOBOND as i32);
        let abs_cur_parity = cur_parity.wrapping_abs();
        let abs_next_parity = next_parity.wrapping_abs();
        let well_defined = |parity: i32| {
            (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
                .contains(&parity)
        };
        let cur_parity_defined = well_defined(abs_cur_parity);
        let next_parity_defined = well_defined(abs_next_parity);
        let cur_action = half_stereo_bond_action(cur_parity, n_unknown, 0, vABParityUnknown);
        let next_action = half_stereo_bond_action(next_parity, n_unknown, 0, vABParityUnknown);
        let mut result_action = cur_action.min(next_action);
        if result_action == -1 {
            stop = true;
        }
        let dot_prod_z = if cur_parity_defined && next_parity_defined {
            let product =
                if chain_len_bits != 0 && (chain_len_bits / MULT_STEREOBOND as i32) % 2 != 0 {
                    triple_prod_char(
                        heap,
                        at.as_const(),
                        at_1,
                        i_next_at_1,
                        &z_dir1,
                        at_2,
                        i_next_at_2,
                        &z_dir2,
                    )?
                } else {
                    dot_prodchar3(&z_dir1, &z_dir2)
                };
            if product.wrapping_abs() < MIN_DOT_PROD as i32 {
                result_action = result_action.min(AB_PARITY_UNDF as i32);
            }
            product
        } else {
            0
        };

        if result_action != AB_PARITY_NONE as i32
            && result_action != -1
            && cur_parity > 0
            && next_parity > 0
        {
            let mut first = heap
                .slice(out_at.as_const())?
                .get(at_1_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let mut second = heap
                .slice(out_at.as_const())?
                .get(opposite_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let saved = save_a_stereo_bond(
                dot_prod_z,
                result_action | chain_len_bits,
                at_1,
                i_next_at_1,
                &mut first.stereo_bond_neighbor,
                &mut first.stereo_bond_ord,
                &mut first.stereo_bond_z_prod,
                &mut first.stereo_bond_parity,
                at_2,
                i_next_at_2,
                &mut second.stereo_bond_neighbor,
                &mut second.stereo_bond_ord,
                &mut second.stereo_bond_z_prod,
                &mut second.stereo_bond_parity,
            );
            if saved != 0 {
                if first.parity == 0
                    || (cur_parity_defined && !well_defined(i32::from(first.parity).wrapping_abs()))
                {
                    first.parity = cur_parity as S_CHAR;
                    first.z_dir = z_dir1;
                }
                if second.parity == 0
                    || (next_parity_defined
                        && !well_defined(i32::from(second.parity).wrapping_abs()))
                {
                    second.parity = next_parity as S_CHAR;
                    second.z_dir = z_dir2;
                }
                first.bAmbiguousStereo = (i32::from(first.bAmbiguousStereo)
                    | i32::from(heap.slice(at.as_const())?[at_1_index].bAmbiguousStereo))
                    as S_CHAR;
                second.bAmbiguousStereo = (i32::from(second.bAmbiguousStereo)
                    | i32::from(
                        heap.slice(at.as_const())?
                            .get(opposite_index)
                            .ok_or(SourceHeapError::PointerOutOfBounds)?
                            .bAmbiguousStereo,
                    )) as S_CHAR;
                heap.slice_mut(out_at)?[at_1_index] = first;
                heap.slice_mut(out_at)?[opposite_index] = second;
            }
        }

        let cur_action = half_stereo_bond_action(cur_parity, n_unknown, 1, vABParityUnknown);
        let next_action = half_stereo_bond_action(next_parity, n_unknown, 1, vABParityUnknown);
        result_action = cur_action.min(next_action);
        cur_parity = abs_cur_parity;
        next_parity = abs_next_parity;
        if result_action != AB_PARITY_NONE as i32 && result_action != -1 {
            if cur_parity > 0 && next_parity > 0 {
                let mut first = heap
                    .slice(out_at.as_const())?
                    .get(at_1_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let mut second = heap
                    .slice(out_at.as_const())?
                    .get(opposite_index)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                let saved = save_a_stereo_bond(
                    dot_prod_z,
                    result_action | chain_len_bits,
                    at_1,
                    i_next_at_1,
                    &mut first.stereo_bond_neighbor2,
                    &mut first.stereo_bond_ord2,
                    &mut first.stereo_bond_z_prod2,
                    &mut first.stereo_bond_parity2,
                    at_2,
                    i_next_at_2,
                    &mut second.stereo_bond_neighbor2,
                    &mut second.stereo_bond_ord2,
                    &mut second.stereo_bond_z_prod2,
                    &mut second.stereo_bond_parity2,
                );
                if saved != 0 {
                    if first.parity2 == 0
                        || (cur_parity_defined
                            && !well_defined(i32::from(first.parity2).wrapping_abs()))
                    {
                        first.parity2 = cur_parity as S_CHAR;
                        if first.parity == 0 {
                            first.z_dir = z_dir1;
                        }
                    }
                    if second.parity2 == 0
                        || (next_parity_defined
                            && !well_defined(i32::from(second.parity2).wrapping_abs()))
                    {
                        second.parity2 = next_parity as S_CHAR;
                        if second.parity == 0 {
                            second.z_dir = z_dir2;
                        }
                    }
                    first.bAmbiguousStereo = (i32::from(first.bAmbiguousStereo)
                        | i32::from(heap.slice(at.as_const())?[at_1_index].bAmbiguousStereo))
                        as S_CHAR;
                    second.bAmbiguousStereo = (i32::from(second.bAmbiguousStereo)
                        | i32::from(
                            heap.slice(at.as_const())?
                                .get(opposite_index)
                                .ok_or(SourceHeapError::PointerOutOfBounds)?
                                .bAmbiguousStereo,
                        )) as S_CHAR;
                    heap.slice_mut(out_at)?[at_1_index] = first;
                    heap.slice_mut(out_at)?[opposite_index] = second;
                    num_stored_isotopic_stereo_bonds += 1;
                }
            }
        } else if result_action == -1 {
            stop = true;
        }
    }

    if stop {
        Ok(CT_CALC_STEREO_ERR)
    } else {
        Ok(num_stored_isotopic_stereo_bonds)
    }
}

#[allow(non_snake_case)]
pub(crate) fn bInpAtomHasRequirdNeigh(
    heap: &SourceHeap,
    at: SourceConstPointer<inp_ATOM>,
    cur_at: i32,
    RequirdNeighType: i32,
    NumDbleBonds: i32,
    bStereoAtZz: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1317 bInpAtomHasRequirdNeigh
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and snapshots add overhead.
    /*
    int bInpAtomHasRequirdNeigh( inp_ATOM *at,
                                 int cur_at,
                                 int RequirdNeighType,
                                 int NumDbleBonds,
                                 int bStereoAtZz )
    {
        /* RequirdNeighType:
              reqired neighbor types (bitmap):
                   0 => any neighbors
                   1 => no terminal hydrogen atom neighbors
                   2 => no terminal -X and -XH together (don't care about -X, -XH bond type, charge, radical)
                        (X = tautomeric endpoint atom)
           NumDbleBonds:
              if non-zero then allow double, alternating and tautomeric bonds
        */
        int i, j, ni, nj, bond_type, num_1s, num_mult, num_other;

        if (at[cur_at].endpoint)
        {
            /*  tautomeric endpoint cannot be a stereo center */
            return 0;
        }

        if (( 1 & RequirdNeighType ) && at[cur_at].num_H)
        {
            return 0;
        }

        if (2 & RequirdNeighType)
        {
            for (i = 0; i < at[cur_at].valence; i++)
            {
                ni = (int) at[cur_at].neighbor[i];
                if (at[ni].valence != 1 ||
                     !get_endpoint_valence( at[ni].el_number ))
                {
                    continue;
                }
                for (j = i + 1; j < at[cur_at].valence; j++)
                {
                    nj = (int) at[cur_at].neighbor[j];
                    if (at[nj].valence != 1 ||
                         at[ni].el_number != at[nj].el_number ||
                         !get_endpoint_valence( at[nj].el_number ))
                    {
                        continue;
                    }
                    /*
                     * if (at[ni].num_H != at[nj].num_H) then the atoms (neighbors of at[cur_at]
                     * are tautomeric endpoints and are indistinguishable => cur_at is not stereogenic
                     * if (at[ni].num_H == at[nj].num_H) then the neighbors are indistinguishable
                     * and cur_at will be found non-sterogenic later
                     * get_endpoint_valence() check will not allow the neighbors to be carbons
                     * Therefore the following "if" is not needed; we may just return 0.
                     */
                    if (at[ni].num_H != at[nj].num_H && strcmp( at[ni].elname, "C" ))
                    {
                        return 0; /*  found -X and -XH neighbors */
                    }
                }
            }
        }

        if (0==bStereoAtZz)
        {
            /* No stereo allowed at centers connected to pseudoelement atoms "Zz","Zy" */
            int nbr, inbr;
            for (inbr = 0; inbr < at[cur_at].valence; inbr++)
            {
                nbr = (int) at[cur_at].neighbor[inbr];
                if (!strcmp( at[nbr].elname, "Zz" ) || !strcmp( at[nbr].elname, "Zy" ))
                {
                    return 0;
                }
            }
        }

        num_1s = num_mult = num_other = 0;

        for (i = 0; i < at[cur_at].valence; i++)
        {
            bond_type = ( at[cur_at].bond_type[i] & ~BOND_MARK_ALL );
            switch (bond_type)
            {
                case BOND_SINGLE:
                    num_1s++;
                    break;
                case BOND_DOUBLE:
                case BOND_ALTERN:
                case BOND_TAUTOM:
                case BOND_ALT12NS:
                    num_mult++;
                    break;
                default:
                    num_other++;
                    break;
            }
        }

        if (num_other)
        {
            return 0;
        }

        if ((NumDbleBonds && NumDbleBonds > num_mult) ||
             (!NumDbleBonds && at[cur_at].valence != num_1s)) /* djb-rwth: addressing LLVM warning */
        {
            return 0;
        }

        return 1;
    }
    */
    // END INCHI C FUNCTION: bInpAtomHasRequirdNeigh

    let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(at)?
        .get(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    if current.endpoint != 0 || (RequirdNeighType & 1 != 0 && current.num_H != 0) {
        return Ok(0);
    }

    if RequirdNeighType & 2 != 0 {
        for i in 0..i32::from(current.valence) {
            let first = heap
                .slice(at)?
                .get(usize::from(current.neighbor[i as usize]))
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            if first.valence != 1 || get_endpoint_valence(first.el_number) == 0 {
                continue;
            }
            for j in i + 1..i32::from(current.valence) {
                let second = heap
                    .slice(at)?
                    .get(usize::from(current.neighbor[j as usize]))
                    .ok_or(SourceHeapError::PointerOutOfBounds)?
                    .clone();
                if second.valence != 1
                    || first.el_number != second.el_number
                    || get_endpoint_valence(second.el_number) == 0
                {
                    continue;
                }
                let first_name_length = first
                    .elname
                    .iter()
                    .position(|&byte| byte == 0)
                    .ok_or(SourceHeapError::MissingNulTerminator)?;
                if first.num_H != second.num_H && first.elname[..first_name_length] != [b'C' as i8]
                {
                    return Ok(0);
                }
            }
        }
    }

    if bStereoAtZz == 0 {
        for i in 0..i32::from(current.valence) {
            let neighbor = heap
                .slice(at)?
                .get(usize::from(current.neighbor[i as usize]))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let length = neighbor
                .elname
                .iter()
                .position(|&byte| byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            if neighbor.elname[..length] == [b'Z' as i8, b'z' as i8]
                || neighbor.elname[..length] == [b'Z' as i8, b'y' as i8]
            {
                return Ok(0);
            }
        }
    }

    let mut num_single = 0_i32;
    let mut num_multiple = 0_i32;
    let mut num_other = 0_i32;
    for i in 0..i32::from(current.valence) {
        let bond_type = u32::from(current.bond_type[i as usize]) & !BOND_MARK_ALL;
        match bond_type {
            BOND_SINGLE => num_single += 1,
            BOND_DOUBLE | BOND_ALTERN | BOND_TAUTOM | BOND_ALT12NS => num_multiple += 1,
            _ => num_other += 1,
        }
    }
    if num_other != 0
        || (NumDbleBonds != 0 && NumDbleBonds > num_multiple)
        || (NumDbleBonds == 0 && i32::from(current.valence) != num_single)
    {
        Ok(0)
    } else {
        Ok(1)
    }
}

#[allow(non_snake_case)]
pub(crate) fn bCanInpAtomBeAStereoCenter(
    heap: &SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
    bPointedEdgeStereo: i32,
    bStereoAtZz: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1432 bCanInpAtomBeAStereoCenter
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and snapshots add overhead.
    /*
    int bCanInpAtomBeAStereoCenter( inp_ATOM *at,
                                    int cur_at,
                                    int bPointedEdgeStereo,
                                    int bStereoAtZz )
    {
    /*************************************************************************************
     * current version
     *************************************************************************************
     *       Use #define to split the stereocenter description table into parts
     *       to make it easier to read
     *
     *                      --------- 4 single bonds stereocenters -------
     *                       0       1       2       3      4       5
     *
     *                       |       |       |       |      |       |
     *                      -C-     -Si-    -Ge-    -Sn-   >As[+]  >B[-]
     *                       |       |       |       |      |       |
     */
    #define SZELEM1         "C\000","Si",   "Ge",   "Sn",   "As", "B\000",
    #define CCHARGE1         0,      0,      0,      0,      1,   -1,
    #define CNUMBONDSANDH1   4,      4,      4,      4,      4,    4,
    #define CCHEMVALENCEH1   4,      4,      4,      4,      4,    4,
    #define CHAS3MEMBRING1   0,      0,      0,      0,      0,    0,
    #define CREQUIRDNEIGH1   0,      0,      0,      0,      3,    0,
    /*
     *                      --------------- S, Se stereocenters ----------
     *                       6       7       8       9      10     11    12    13
     *
     *                               |       |      ||             |     |     ||
     *                      -S=     =S=     -S[+]   >S[+]   -Se=  =Se=  -Se[+] >Se[+]
     *                       |       |       |       |       |     |     |      |
     */
    #define SZELEM2         "S\000","S\000","S\000","S\000","Se", "Se", "Se",  "Se",
    #define CCHARGE2         0,      0,      1,      1,      0,    0,    1,     1,
    #define CNUMBONDSANDH2   3,      4,      3,      4,      3,    4,    3,     4,
    #define CCHEMVALENCEH2   4,      6,      3,      5,      4,    6,    3,     5,
    #define CHAS3MEMBRING2   0,      0,      0,      0,      0,    0,    0,     0,
    #define CREQUIRDNEIGH2   3,      3,      3,      3,      3,    3,    3,     3,
    /*
     *                      ------------------ N, P stereocenters -----------------
     *                        14     15       16     17     18       19       20
     *
     *                                                             Phosphine Arsine
     *                                      X---Y
     *                        |      |       \ /     |       |       \ /      \ /
     *                       =N-    >N[+]     N     >P[+]   =P-       P        As
     *                        |      |        |      |       |        |        |
     */
    #define SZELEM3         "N\000","N\000","N\000","P\000","P\000","P\000", "As",
    #define CCHARGE3         0,      1,      0,      1,      0,      0,       0,
    #define CNUMBONDSANDH3   4,      4,      3,      4,      4,      3,       3,
    #define CCHEMVALENCEH3   5,      4,      3,      4,      5,      3,       3,
    #define CHAS3MEMBRING3   0,      0,      1,      0,      0,      0,       0,
    #define CREQUIRDNEIGH3   3,      3,      1,      3,      3,      2,       2,

    #define PHOSPHINE_STEREO  19  /* the number must match Phosphine number in the comments, see above */
    #define ARSINE_STEREO     20  /* the number must match Arsine number in the comments, see above */

        static const char        szElem[][3] = { SZELEM1         SZELEM2         SZELEM3 };
        static const S_CHAR        cCharge[] = { CCHARGE1        CCHARGE2        CCHARGE3 };
        static const S_CHAR  cNumBondsAndH[] = { CNUMBONDSANDH1  CNUMBONDSANDH2  CNUMBONDSANDH3 };
        static const S_CHAR  cChemValenceH[] = { CCHEMVALENCEH1  CCHEMVALENCEH2  CCHEMVALENCEH3 };
        static const S_CHAR  cHas3MembRing[] = { CHAS3MEMBRING1  CHAS3MEMBRING2  CHAS3MEMBRING3 };
        static const S_CHAR  cRequirdNeigh[] = { CREQUIRDNEIGH1  CREQUIRDNEIGH2  CREQUIRDNEIGH3 };

        static const int n = sizeof( szElem ) / sizeof( szElem[0] );
        /* reqired neighbor types (bitmap):
           0 => check bonds only
           1 => no terminal hydrogen atom neighbors
           2 => no terminal -X and -XH together (don't care the bond type, charge, radical)
                (X = tautomeric endpoint atom)
           Note: whenever cChemValenceH[] > cNumBondsAndH[]
                 the tautomeric and/or alternating bonds
                 are permitted

        */
        int i, ret = 0;
        for (i = 0; i < n; i++)
        {
            if (!strcmp( at[cur_at].elname, szElem[i] ) &&
                ( at[cur_at].charge == cCharge[i]
    #ifdef ALLOW_NO_CHARGE_ON_STEREO_CENTERS
                  || at[cur_at].charge == 0
    #endif
                  ) &&
                  ( !at[cur_at].radical || at[cur_at].radical == 1 ) &&
                 at[cur_at].valence + at[cur_at].num_H == cNumBondsAndH[i] &&
                 at[cur_at].chem_bonds_valence + at[cur_at].num_H == cChemValenceH[i] &&
                 ( cHas3MembRing[i] ? is_atom_in_3memb_ring( at, cur_at ) : 1 ) &&
                 bInpAtomHasRequirdNeigh( at, cur_at, cRequirdNeigh[i], cChemValenceH[i] - cNumBondsAndH[i], bStereoAtZz ))
                /*
                if (!strcmp( at[cur_at].elname, szElem[i] ) &&
                 at[cur_at].charge == cCharge[i] &&
                 ( !at[cur_at].radical || at[cur_at].radical == 1 ) &&
                 at[cur_at].valence + at[cur_at].num_H == cNumBondsAndH[i] &&
                 at[cur_at].chem_bonds_valence + at[cur_at].num_H == cChemValenceH[i] &&
                 ( cHas3MembRing[i] ? is_atom_in_3memb_ring( at, cur_at ) : 1 ) &&
                 bInpAtomHasRequirdNeigh( at, cur_at, cRequirdNeigh[i], cChemValenceH[i] - cNumBondsAndH[i], bStereoAtZz ))*/
            {
                ret = cNumBondsAndH[i];
                break;
            }
        }

        if (i == PHOSPHINE_STEREO && !( bPointedEdgeStereo & PES_BIT_PHOSPHINE_STEREO ))
        {
            ret = 0;
        }
        if (i == ARSINE_STEREO && !( bPointedEdgeStereo & PES_BIT_ARSINE_STEREO ))
        {
            ret = 0;
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: bCanInpAtomBeAStereoCenter
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: bCanInpAtomBeAStereoCenter
    // INCHI✔️❌: NEW_STEREOCENTER_CHECK == 1 selects this table-driven implementation.
    // INCHI✔️❌: ALLOW_NO_CHARGE_ON_STEREO_CENTERS is undefined.
    // INCHI✔️❌: PES_BIT_PHOSPHINE_STEREO == 2 and PES_BIT_ARSINE_STEREO == 4.
    // END INCHI ACTIVE MACRO CONFIGURATION: bCanInpAtomBeAStereoCenter

    const ELEMENTS: [&[i8]; 21] = [
        &[b'C' as i8],
        &[b'S' as i8, b'i' as i8],
        &[b'G' as i8, b'e' as i8],
        &[b'S' as i8, b'n' as i8],
        &[b'A' as i8, b's' as i8],
        &[b'B' as i8],
        &[b'S' as i8],
        &[b'S' as i8],
        &[b'S' as i8],
        &[b'S' as i8],
        &[b'S' as i8, b'e' as i8],
        &[b'S' as i8, b'e' as i8],
        &[b'S' as i8, b'e' as i8],
        &[b'S' as i8, b'e' as i8],
        &[b'N' as i8],
        &[b'N' as i8],
        &[b'N' as i8],
        &[b'P' as i8],
        &[b'P' as i8],
        &[b'P' as i8],
        &[b'A' as i8, b's' as i8],
    ];
    const CHARGES: [i8; 21] = [
        0, 0, 0, 0, 1, -1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
    ];
    const NUM_BONDS_AND_H: [i8; 21] = [
        4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 3, 3,
    ];
    const CHEM_VALENCE_H: [i8; 21] = [
        4, 4, 4, 4, 4, 4, 4, 6, 3, 5, 4, 6, 3, 5, 5, 4, 3, 4, 5, 3, 3,
    ];
    const HAS_3_MEMBER_RING: [i8; 21] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    ];
    const REQUIRED_NEIGHBOR: [i8; 21] = [
        0, 0, 0, 0, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 2, 2,
    ];
    const PHOSPHINE_STEREO: usize = 19;
    const ARSINE_STEREO: usize = 20;

    let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(at.as_const())?
        .get(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let element_length = current
        .elname
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let element = &current.elname[..element_length];
    let mut matched_index = ELEMENTS.len();
    let mut result = 0_i32;
    for i in 0..ELEMENTS.len() {
        if element == ELEMENTS[i]
            && current.charge == CHARGES[i]
            && (current.radical == 0 || current.radical == 1)
            && i32::from(current.valence) + i32::from(current.num_H)
                == i32::from(NUM_BONDS_AND_H[i])
            && i32::from(current.chem_bonds_valence) + i32::from(current.num_H)
                == i32::from(CHEM_VALENCE_H[i])
            && (HAS_3_MEMBER_RING[i] == 0 || is_atom_in_3memb_ring(heap, at, cur_at)? != 0)
            && bInpAtomHasRequirdNeigh(
                heap,
                at.as_const(),
                cur_at,
                i32::from(REQUIRED_NEIGHBOR[i]),
                i32::from(CHEM_VALENCE_H[i]) - i32::from(NUM_BONDS_AND_H[i]),
                bStereoAtZz,
            )? != 0
        {
            matched_index = i;
            result = i32::from(NUM_BONDS_AND_H[i]);
            break;
        }
    }
    if matched_index == PHOSPHINE_STEREO
        && bPointedEdgeStereo & PES_BIT_PHOSPHINE_STEREO as i32 == 0
    {
        result = 0;
    }
    if matched_index == ARSINE_STEREO && bPointedEdgeStereo & PES_BIT_ARSINE_STEREO as i32 == 0 {
        result = 0;
    }
    Ok(result)
}

#[allow(non_snake_case)]
pub(crate) fn can_be_a_stereo_atom_with_isotopic_H(
    heap: &SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
    bPointedEdgeStereo: i32,
    bStereoAtZz: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3690 can_be_a_stereo_atom_with_isotopic_H
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int can_be_a_stereo_atom_with_isotopic_H( inp_ATOM *at,
                                              int cur_at,
                                              int bPointedEdgeStereo,
                                              int bStereoAtZz )
    {
        int nNumNeigh;
        if (( nNumNeigh = bCanInpAtomBeAStereoCenter( at, cur_at, bPointedEdgeStereo, bStereoAtZz ) ) &&
             at[cur_at].valence + at[cur_at].num_H == nNumNeigh &&
             at[cur_at].num_H <= NUM_H_ISOTOPES
           )
        {
            return 1;
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: can_be_a_stereo_atom_with_isotopic_H
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: can_be_a_stereo_atom_with_isotopic_H
    // INCHI✔️❌: NEW_STEREOCENTER_CHECK == 1 and NUM_H_ISOTOPES == 3.
    // END INCHI ACTIVE MACRO CONFIGURATION: can_be_a_stereo_atom_with_isotopic_H

    let number_of_neighbors =
        bCanInpAtomBeAStereoCenter(heap, at, cur_at, bPointedEdgeStereo, bStereoAtZz)?;
    if number_of_neighbors != 0 {
        let current = heap
            .slice(at.as_const())?
            .get(usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if i32::from(current.valence) + i32::from(current.num_H) == number_of_neighbors
            && i32::from(current.num_H) <= NUM_H_ISOTOPES as i32
        {
            return Ok(1);
        }
    }
    Ok(0)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn GetStereocenter0DParity(
    _pCG: &mut CANON_GLOBALS,
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
    j1: i32,
    nSbNeighOrigAtNumb: &mut [AT_NUMB],
    nFlag: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3734 GetStereocenter0DParity
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap and slice access add overhead.
    /*
    int GetStereocenter0DParity( CANON_GLOBALS *pCG,
                                 inp_ATOM *at,
                                 int cur_at,
                                 int j1,
                                 AT_NUMB nSbNeighOrigAtNumb[],
                                 int nFlag )
    {
        int parity = AB_PARITY_NONE;
        if (at[cur_at].p_parity && ( j1 == MAX_NUM_STEREO_ATOM_NEIGH - 1 || j1 == MAX_NUM_STEREO_ATOM_NEIGH ))
        {
            int i, num_trans_inp, num_trans_neigh;
            AT_NUMB nInpNeighOrigAtNumb[MAX_NUM_STEREO_ATOM_NEIGH];
            for (i = 0; i < MAX_NUM_STEREO_ATOM_NEIGH; i++)
            {
                nInpNeighOrigAtNumb[i] = at[cur_at].p_orig_at_num[i];
                if (nInpNeighOrigAtNumb[i] == at[cur_at].orig_at_number)
                {
                    nInpNeighOrigAtNumb[i] = 0; /* lone pair or explicit H */
                }
            }

            num_trans_inp = insertions_sort( pCG, nInpNeighOrigAtNumb, MAX_NUM_STEREO_ATOM_NEIGH, sizeof( nInpNeighOrigAtNumb[0] ), comp_AT_NUMB );
            num_trans_neigh = insertions_sort( pCG, nSbNeighOrigAtNumb, j1, sizeof( nSbNeighOrigAtNumb[0] ), comp_AT_NUMB );

            if (j1 == MAX_NUM_STEREO_ATOM_NEIGH - 1)
            {
                ;
                /*num_trans_neigh += j1;*/
                /* the lone pair or implicit H is implicitly at the top of the list */
            }
            if (!memcmp( nInpNeighOrigAtNumb + MAX_NUM_STEREO_ATOM_NEIGH - j1, nSbNeighOrigAtNumb, j1 * sizeof( AT_NUMB ) ))
            {
                if (ATOM_PARITY_WELL_DEF( at[cur_at].p_parity ))
                {
                    parity = 2 - ( num_trans_inp + num_trans_neigh + at[cur_at].p_parity ) % 2;
                }
                else
                {
                    parity = at[cur_at].p_parity;
                }
                at[cur_at].bUsed0DParity |= nFlag; /* 0D parity used for streocenter parity */
            }
        }

        return parity;
    }
    */
    // END INCHI C FUNCTION: GetStereocenter0DParity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: GetStereocenter0DParity
    // INCHI✔️❌: MAX_NUM_STEREO_ATOM_NEIGH == 4; AT_NUMB is uint16_t.
    // INCHI✔️❌: ATOM_PARITY_WELL_DEF(p) is AB_MIN_WELL_DEFINED_PARITY <= p <= AB_MAX_WELL_DEFINED_PARITY.
    // INCHI✔️❌: memcmp compares the j1 sorted AT_NUMB values byte-for-byte; equality is identical to typed slice equality.
    // END INCHI ACTIVE MACRO CONFIGURATION: GetStereocenter0DParity

    let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let (input_parity, original_number, mut input_neighbors) = {
        let current = heap
            .slice(at.as_const())?
            .get(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        (
            i32::from(current.p_parity),
            current.orig_at_number,
            current.p_orig_at_num,
        )
    };
    let number_of_neighbors = match j1 {
        value if value == MAX_NUM_STEREO_ATOM_NEIGH as i32 - 1 => 3_usize,
        value if value == MAX_NUM_STEREO_ATOM_NEIGH as i32 => 4_usize,
        _ => return Ok(AB_PARITY_NONE as i32),
    };
    if input_parity == 0 {
        return Ok(AB_PARITY_NONE as i32);
    }
    if nSbNeighOrigAtNumb.len() < number_of_neighbors {
        return Err(SourceHeapError::PointerOutOfBounds);
    }

    for neighbor in &mut input_neighbors {
        if *neighbor == original_number {
            *neighbor = 0;
        }
    }
    let input_transactions = insertions_sort(
        bytemuck::cast_slice_mut::<AT_NUMB, u8>(&mut input_neighbors),
        MAX_NUM_STEREO_ATOM_NEIGH as usize,
        std::mem::size_of::<AT_NUMB>(),
        &mut |left, right| {
            Ok(comp_AT_NUMB(
                AT_NUMB::from_ne_bytes([left[0], left[1]]),
                AT_NUMB::from_ne_bytes([right[0], right[1]]),
            ))
        },
    )?;
    let neighbor_transactions = insertions_sort(
        bytemuck::cast_slice_mut::<AT_NUMB, u8>(nSbNeighOrigAtNumb),
        number_of_neighbors,
        std::mem::size_of::<AT_NUMB>(),
        &mut |left, right| {
            Ok(comp_AT_NUMB(
                AT_NUMB::from_ne_bytes([left[0], left[1]]),
                AT_NUMB::from_ne_bytes([right[0], right[1]]),
            ))
        },
    )?;

    if input_neighbors[MAX_NUM_STEREO_ATOM_NEIGH as usize - number_of_neighbors..]
        != nSbNeighOrigAtNumb[..number_of_neighbors]
    {
        return Ok(AB_PARITY_NONE as i32);
    }

    let parity = if (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
        .contains(&input_parity)
    {
        2_i32.wrapping_sub(
            input_transactions
                .wrapping_add(neighbor_transactions)
                .wrapping_add(input_parity)
                .wrapping_rem(2),
        )
    } else {
        input_parity
    };
    let atoms = heap.slice_mut(at)?;
    let current = atoms
        .get_mut(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    current.bUsed0DParity = (i32::from(current.bUsed0DParity) | nFlag) as i8;
    Ok(parity)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn set_stereo_atom_parity(
    pCG: &mut CANON_GLOBALS,
    heap: &mut SourceHeap,
    out_at: SourceMutPointer<sp_ATOM>,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
    at_removed_H: SourceConstPointer<inp_ATOM>,
    num_removed_H: i32,
    bPointedEdgeStereo: i32,
    vABParityUnknown: i32,
    LooseTSACheck: i32,
    bStereoAtZz: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3790 set_stereo_atom_parity
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and atom snapshots add overhead.
    /*
    int set_stereo_atom_parity( CANON_GLOBALS *pCG,
                                sp_ATOM *out_at,
                                inp_ATOM *at,
                                int cur_at,
                                inp_ATOM *at_removed_H,
                                int num_removed_H,
                                int bPointedEdgeStereo,
                                int vABParityUnknown,
                                int LooseTSACheck,
                                int bStereoAtZz )
    {
        int    j, k, next_at, num_z, j1, nType, num_explicit_H, tot_num_iso_H, nMustHaveNumNeigh;
        int    num_explicit_iso_H[NUM_H_ISOTOPES + 1];
                        /*  numbers of removed hydrogen atoms   */
        int    index_H[MAX_NUM_STEREO_ATOM_NEIGH];
                        /*  cannot have more than 4 elements: 1 H, 1 D, 1 T atom(s) */
        double z, sum_xyz[3], min_sine, triple_product;
        double at_coord[MAX_NUM_STEREO_ATOM_NEIGH][3];
        double bond_len_xy[4], rmax = 0.0, rmin = 0.0;
        double at_coord_center[3];
        int    parity, bAmbiguous = 0, bAddExplicitNeighbor = 0, b2D = 0, n2DTetrahedralAmbiguity = 0;
        int    bIgnoreIsotopicH = ( 0 != ( at[cur_at].cFlags & AT_FLAG_ISO_H_POINT ) );
        AT_NUMB nSbNeighOrigAtNumb[MAX_NUM_STEREO_ATOM_NEIGH];

        double vMinAngle = MIN_ANGLE;
        double vMinSine = MIN_SINE;
        if (LooseTSACheck)
        {
            /* Relax check for principal 3-atomic angles close to 180 deg, for 2D tetrahedron, to account for           */
            /* large cycles drawn as appear after applying some structure 'cleaning' algorithms/softwares               */

            if (at[cur_at].nNumAtInRingSystem >= 3) /* Ensure that central atom is in ring (ideally, one should check   */
                                                    /* also that neighbors are in-ring, and for large enough ring)      */
            {
                vMinAngle = MIN_ANGLE_RELAXED;
                vMinSine = MIN_SINE_RELAXED;
            }
        }

        out_at[cur_at].parity =
            out_at[cur_at].parity2 =
            out_at[cur_at].stereo_atom_parity =
            out_at[cur_at].stereo_atom_parity2 =
            AB_PARITY_NONE;
        parity = AB_PARITY_NONE;

        memset( num_explicit_iso_H, 0, sizeof( num_explicit_iso_H ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        num_explicit_H = 0;

    #if ( NEW_STEREOCENTER_CHECK == 1 )
        if (!( nMustHaveNumNeigh = bCanInpAtomBeAStereoCenter( at, cur_at, bPointedEdgeStereo, bStereoAtZz ) ) ||
             at[cur_at].num_H > NUM_H_ISOTOPES)
        {
            goto exit_function;
        }
    #else
        nMustHaveNumNeigh = MAX_NUM_STEREO_ATOM_NEIGH;
        if (!bCanAtomBeAStereoCenter( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical ) ||
             at[cur_at].valence + at[cur_at].num_H != nMustHaveNumNeigh ||
             at[cur_at].num_H > NUM_H_ISOTOPES
           )
        {
            goto exit_function;
        }
        for (j = 0; j < at[cur_at].valence; j++)
        {
            if (( at[cur_at].bond_type[j] & ~BOND_MARK_ALL ) != BOND_SINGLE)
            {
                goto exit_function;
            }
        }
    #endif

        /*  numbers of isotopic H atoms */
        for (j = 0, tot_num_iso_H = 0; j < NUM_H_ISOTOPES; j++)
        {
            if (at[cur_at].num_iso_H[j] > 1)
            {
                goto exit_function; /*  two or more identical hydrogen isotopic neighbors */
            }
            tot_num_iso_H += at[cur_at].num_iso_H[j];
        }
        if (bIgnoreIsotopicH)
        {
            tot_num_iso_H = 0; /* isotopic H considered subject to exchange => ignore isotopic */
        }
        /*  number of non-isotopic H atoms */
        if (at[cur_at].num_H - tot_num_iso_H > 1)
        {
            goto exit_function; /*  two or more identical hydrogen non-isotopic neighbors */
        }

        /*  count removed explicit terminal hydrogens attached to at[cur_at]. */
        /*  the result is num_explicit_H. */
        /*  Removed hydrogens are sorted in increasing isotopic shift order */
        if (at_removed_H && num_removed_H > 0)
        {
            for (j = 0; j < num_removed_H; j++)
            {
                if (at_removed_H[j].neighbor[0] == cur_at)
                {
                    k = at_removed_H[j].iso_atw_diff;
                    /*  iso_atw_diff values: H=>0, 1H=>1, D=2H=>2, T=3H=>3 */
                    if (k < 0 || k > NUM_H_ISOTOPES || bIgnoreIsotopicH)
                    {
                        k = 0; /*  treat wrong H isotopes as non-isotopic H */
                    }
                    num_explicit_iso_H[k] ++;
                    index_H[num_explicit_H++] = j;
                }
            }
        }

        /*  coordinates initialization */
        num_z = 0;
        sum_xyz[0] = sum_xyz[1] = sum_xyz[2] = 0.0;

        at_coord_center[0]      =
            at_coord_center[1]  =
            at_coord_center[2]  =
            0.0;

            /*  fill out stereo center neighbors coordinates */
            /*  and obtain the parity from the geometry */

        for (k = 0, j1 = 0; k < 2; k++)
        {
            switch (k)
            {

                case 0:
                    /*   add coordinates of removed hydrogens */
                    for (j = 0; j < num_explicit_H; j++, j1++)
                    {
                        next_at = index_H[j];
                        /*  Use bond description located at removed_H atom */
                        /*  minus sign at get_z_coord: at_removed_H[] contains bonds TO at[cur_at], not FROM it. */
                        /*  Note: &at[(at_removed_H-at)+ next_at] == &at_removed_H[next_at] */
                        z = -get_z_coord( at,
                                          (int) ( at_removed_H - at ) + next_at,
                                          0 /*neighbor #*/,
                                          &nType,
                                          -( bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO ) );
                        switch (nType)
                        {
                            case ZTYPE_EITHER:
                                parity = vABParityUnknown /*AB_PARITY_UNKN*/; /*  no parity: bond in "Either" direction. */
                                goto exit_function;
                            case ZTYPE_UP:
                            case ZTYPE_DOWN:
                                nType = -nType; /*  at_removed_H[] contains bonds TO the center, not from */
                                b2D++;
                                /*  no break; here */
                            case ZTYPE_3D:
                                num_z++;
                        }

                        if (j1 < MAX_NUM_STEREO_ATOM_NEIGH) /* djb-rwth: fixing oss-fuzz issue #71142 */
                        {
                            nSbNeighOrigAtNumb[j1] = at_removed_H[next_at].orig_at_number;
                            at_coord[j1][0] = at_removed_H[next_at].x - at[cur_at].x;
                            at_coord[j1][1] = at_removed_H[next_at].y - at[cur_at].y;
                            bond_len_xy[j1] = len2(at_coord[j1]);
                            /* bond_len_xy[j1] = sqrt(at_coord[j1][0]*at_coord[j1][0]+at_coord[j1][1]*at_coord[j1][1]); */
                            at_coord[j1][2] = (nType == ZTYPE_3D ? z : nType == ZTYPE_UP
                                ? bond_len_xy[j1] : nType == ZTYPE_DOWN
                                ? -bond_len_xy[j1] : 0.0);
                        }
                        else
                        {
                            break;
                        }
                    }
                    break;
                case 1:
                    /*  add all coordinates of other neighboring atoms */
                    for (j = 0; j < at[cur_at].valence; j++, j1++)
                    {
                        next_at = at[cur_at].neighbor[j];
                        z = get_z_coord( at, cur_at, j, &nType, ( bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO ) );
                        switch (nType)
                        {
                            case ZTYPE_EITHER:
                                parity = vABParityUnknown /*AB_PARITY_UNKN*/; /*  unknown parity: bond in "Either" direction. */
                                goto exit_function;
                            case ZTYPE_UP:
                            case ZTYPE_DOWN:
                                b2D++;
                            case ZTYPE_3D:
                                num_z++;
                        }

                        if (j1 < MAX_NUM_STEREO_ATOM_NEIGH) /* djb-rwth: fixing oss-fuzz issue #71142 */
                        {
                            nSbNeighOrigAtNumb[j1] = at[next_at].orig_at_number;
                            at_coord[j1][0] = at[next_at].x - at[cur_at].x;
                            at_coord[j1][1] = at[next_at].y - at[cur_at].y;
                            bond_len_xy[j1] = len2(at_coord[j1]);
                            /* bond_len_xy[j1] = sqrt(at_coord[j1][0]*at_coord[j1][0]+at_coord[j1][1]*at_coord[j1][1]); */
                            at_coord[j1][2] = (nType == ZTYPE_3D ? z :
                                nType == ZTYPE_UP ? bond_len_xy[j1] :
                                nType == ZTYPE_DOWN ? -bond_len_xy[j1] : 0.0);
                        }
                        else
                        {
                            break;
                        }
                    }
                    break;
            }
        }
        /* j1 is the number of explicit neighbors (that is, all neighbors except implicit H) */
        b2D = ( b2D == num_z && num_z );  /*  1 => two-dimensional */

        if (MAX_NUM_STEREO_ATOM_NEIGH != at[cur_at].valence + num_explicit_H &&
             MAX_NUM_STEREO_ATOM_NEIGH - 1 != at[cur_at].valence + num_explicit_H)
        {
            /*  not enough geometry data to find the central atom parity */
            if (nMustHaveNumNeigh == at[cur_at].valence + at[cur_at].num_H &&
                at[cur_at].num_H > 1)
            {
                /*  only isotopic parity is possible; no non-isotopic parity */
                if (parity == vABParityUnknown /*AB_PARITY_UNKN*/)
                {
                    parity = -vABParityUnknown /*AB_PARITY_UNKN*/; /*  the user marked the center as "unknown" */
                }
                else
                {
                    parity = -AB_PARITY_UNDF; /*  not enough geometry; only isotopic parity is possible */
                }
            }
            else
            {
                parity = AB_PARITY_NONE;      /*  not a stereocenter at all */
            }
            goto exit_function;
        }
        /*  make all vector lengths equal to 1; exit if too short. 9-10-2002 */
        for (j = 0; j < j1; j++)
        {
            z = len3( at_coord[j] );
            if (z < MIN_BOND_LEN)
            {
                /* bond length is too small: use 0D parities */
                if (AB_PARITY_NONE == ( parity = GetStereocenter0DParity( pCG, at, cur_at,
                                                                          j1,nSbNeighOrigAtNumb,
                                                                          FlagSC_0D ) ))
                {
                    parity = AB_PARITY_UNDF;
                }
                goto exit_function;
            }
    #if ( STEREO_CENTER_BONDS_NORM == 1 )
            else
            {
                mult3( at_coord[j], 1.0 / z, at_coord[j] );
            }
    #endif
            rmax = j ? inchi_max( rmax, z ) : z;
            rmin = j ? inchi_min( rmin, z ) : z;
        }
        if (rmin / rmax < MIN_SINE)
        {
            /* bond ratio is too small: use 0D parities */
            if (AB_PARITY_NONE == ( parity = GetStereocenter0DParity( pCG, at, cur_at, j1,
                nSbNeighOrigAtNumb,
                FlagSC_0D ) ))
            {
                parity = AB_PARITY_UNDF;
            }
            goto exit_function;
        }
        for (j = 0; j < j1; j++)
        {
            add3( sum_xyz, at_coord[j], sum_xyz );
        }

        /*  Here j1 is a number of neighbors including explicit terminal isotopic H */
        /*  num_explicit_iso_H[0] = number of explicit non-isotopic hydrogen atom neighbors */
        j = j1;
        /*  Add Explicit Neighbor */
        if (j1 == MAX_NUM_STEREO_ATOM_NEIGH - 1)
        {
            /*  Add an explicit neighbor if possible */
            if (nMustHaveNumNeigh == MAX_NUM_STEREO_ATOM_NEIGH - 1)
            {
                bAddExplicitNeighbor = ADD_EXPLICIT_LONE_PAIR_NEIGH;
            }
            else
            {
                if (nMustHaveNumNeigh == MAX_NUM_STEREO_ATOM_NEIGH)
                {
                    /*  Check whether an explicit non-isotopic hydrogen can be added */
                    /*  to an atom that is a stereogenic atom */
                    if (1 == at[cur_at].num_H - num_explicit_H &&     /*  the atom has only one one implicit hydrogen */
                        1 == at[cur_at].num_H - tot_num_iso_H)
                    {
                        /*  this hydrogen is non-isotopic */
                        bAddExplicitNeighbor = ADD_EXPLICIT_HYDROGEN_NEIGH;
                    }
                }
            }
        }

        if (bAddExplicitNeighbor)
        {
            /***********************************************************
             * May happen only if (j1 == MAX_NUM_STEREO_ATOM_NEIGH-1)
             * 3 neighbors only, no H-neighbors. Create and add coordinates of an implicit H
             * or a fake 4th neighbor, that is, a lone pair
             */
            if (parity == vABParityUnknown /*AB_PARITY_UNKN*/)
            {
                goto exit_function;  /*  the user insists the parity is unknown and the isotopic */
                                     /*  composition of the neighbors does not contradict */
            }
            else
            {
                if (num_z == 0 || are_3_vect_in_one_plane( at_coord, vMinSine ))
                {
                    /*  "hydrogen down" rule is needed to resolve an ambiguity */
                    if (num_z > 0)
                    {
                        bAmbiguous |= AMBIGUOUS_STEREO;
                    }
    #if ( APPLY_IMPLICIT_H_DOWN_RULE == 1 )  /* { */
                    /*  Although H should be at the top of the list, add it to the bottom. */
                    /*  This will be taken care of later by inverting parity 1<->2 */
                    at_coord[j][0] = 0.0;
                    at_coord[j][1] = 0.0;
    #if ( STEREO_CENTER_BONDS_NORM == 1 )
                    at_coord[j][2] = -1.0;
    #else
                    at_coord[j][2] = -( bond_len_xy[0] + bond_len_xy[1] + bond_len_xy[2] ) / 3.0;
    #endif
    #else /* } APPLY_IMPLICIT_H_DOWN_RULE { */
    #if (ALWAYS_SET_STEREO_PARITY == 1)
                    parity = AB_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
    #else
                    /* All 3 bonds are in one plane: try to get 0D parities */
                    if (AB_PARITY_NONE == ( parity = GetStereocenter0DParity( pCG, at, cur_at,
                                                                              j1,
                                                                              nSbNeighOrigAtNumb,
                                                                              FlagSC_0D ) ))
                    {
                        parity = AB_PARITY_UNDF;
                    }
                    /*parity =  AB_PARITY_UNDF;*/ /*  no parity can be calculated found */
    #endif
                    goto exit_function;
    #endif /* } APPLY_IMPLICIT_H_DOWN_RULE */
                }
                else
                {
                    /*  we have enough information to find implicit hydrogen coordinates */
                    /*
                    at_coord[j][0] = -sum_x;
                    at_coord[j][1] = -sum_y;
                    at_coord[j][2] = -sum_z;
                    */
                    /* copy3(sum_xyz, at_coord[j]); -- djb-rwth: removing copy3 function */
                    memcpy(at_coord[j], sum_xyz, sizeof(sum_xyz));
                    change_sign3( at_coord[j], at_coord[j] );
                    z = len3( at_coord[j] );
    #if ( FIX_STEREO_SCALING_BUG == 1 )
                    if (z > 1.0)
                    {
                        rmax *= z;
                    }
                    else
                    {
                        rmin *= z;
                    }
    #else
                    /* Comparing the original bond lengths to lenghts derived from normalized to 1 */
                    /* This bug leads to pronouncing legitimate stereogenic atoms */
                    /* connected by 3 bonds "undefined" if in a nicely drawn 2D structure */
                    /* bond lengths are about 20 or greater. Reported by Reinhard Dunkel 2005-08-05 */
                    if (bPointedEdgeStereo & PES_BIT_FIX_SP3_BUG)
                    {
                        /* coordinate scaling bug fixed here */
                        if (z > 1.0)
                        {
                            rmax *= z;
                        }
                        else
                        {
                            rmin *= z;
                        }
                    }
                    else
                    {
                        /* original InChI v.1 bug */
                        rmax = inchi_max( rmax, z );
                        rmin = inchi_min( rmin, z );
                    }
    #endif
                    if (z < MIN_BOND_LEN || rmin / rmax < MIN_SINE)
                    {
                        /* the new 4th bond is too short: try to get 0D parities */
                        if (AB_PARITY_NONE == ( parity = GetStereocenter0DParity( pCG, at, cur_at, j1, nSbNeighOrigAtNumb, FlagSC_0D ) ))
                        {
                            parity = AB_PARITY_UNDF;
                        }
                        goto exit_function;
                    }
    #if ( STEREO_CENTER_BOND4_NORM == 1 )
                    else
                    {
                        mult3( at_coord[j], 1.0 / z, at_coord[j] );
                    }
    #endif
                }
            }
        }

        else if (j1 != MAX_NUM_STEREO_ATOM_NEIGH)
        {
            if (parity == vABParityUnknown /*AB_PARITY_UNKN*/)
            {
                parity = -AB_PARITY_UNDF; /*  isotopic composition of H-neighbors contradicts 'unknown' */
            }
            goto exit_function;
        }

        else /*  j1 == MAX_NUM_STEREO_ATOM_NEIGH */
        {
            if (num_z == 0 || are_4at_in_one_plane( at_coord, vMinSine ))
            {
                /*  all four neighours in xy plane: undefined geometry. */
                if (num_z > 0)
                {
                    bAmbiguous |= AMBIGUOUS_STEREO;
                }
                if (parity != vABParityUnknown /*AB_PARITY_UNKN*/)
                {
    #if (ALWAYS_SET_STEREO_PARITY == 1)
                    parity = AB_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
    #else
                    /* All 4 bonds are in one plane: try to get 0D parities */
                    if (AB_PARITY_NONE == ( parity = GetStereocenter0DParity( pCG, at, cur_at, j1, nSbNeighOrigAtNumb, FlagSC_0D ) ))
                    {
                        parity = AB_PARITY_UNDF;
                    }
                    else if (ATOM_PARITY_WELL_DEF( parity ))
                    {
                        bAmbiguous &= ~AMBIGUOUS_STEREO; /* 0D parity has resolved the ambiguity */
                    }
    #endif
                }
                goto exit_function;
            }
        }

        /***********************************************************
        * At this point we have 4 neighboring atoms.
        * check for tetrahedral ambiguity in 2D case
        */
        if (b2D)
        {

            n2DTetrahedralAmbiguity = Get2DTetrahedralAmbiguity( pCG,at_coord,
                                                                 bAddExplicitNeighbor,
                                                                 ( bPointedEdgeStereo & PES_BIT_FIX_SP3_BUG ),
                                                                 vMinAngle );

            if (0 < n2DTetrahedralAmbiguity)
            {
                if (T2D_WARN & n2DTetrahedralAmbiguity)
                {
                    bAmbiguous |= AMBIGUOUS_STEREO;
                }
                if (T2D_UNDF & n2DTetrahedralAmbiguity)
                {
                    if (parity != vABParityUnknown /*AB_PARITY_UNKN*/)
                    {
    #if (ALWAYS_SET_STEREO_PARITY == 1)
                        parity = AB_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
    #else
                        parity = AB_PARITY_UNDF; /*  no parity */
    #endif
                    }
                    goto exit_function;
                }
            }
            else if (n2DTetrahedralAmbiguity < 0)
            {
                bAmbiguous |= AMBIGUOUS_STEREO_ERROR; /*  error */
                parity = AB_PARITY_UNDF;
                goto exit_function;
            }
        }

        /************************************************************/
        /*  Move coordinates origin to the neighbor #0 */
        for (j = 1; j < MAX_NUM_STEREO_ATOM_NEIGH; j++)
        {
            diff3( at_coord[j], at_coord[0], at_coord[j] );
        }
        diff3( at_coord_center, at_coord[0], at_coord_center );

        /*
        for ( k = 0; k < 3; k++ )
        {
            for ( j = 1; j < MAX_NUM_STEREO_ATOM_NEIGH; j ++ )
            {
                at_coord[j][k] -= at_coord[0][k];
            }
            at_coord_center[k] -= at_coord[0][k];
        }
        */


        /********************************************************
         *  Find the central (cur_at) atom's parity
         *  (orientation of atoms #1-3 when looking from #0)
         ********************************************************/
        triple_product = triple_prod_and_min_abs_sine2( &at_coord[1], at_coord_center, bAddExplicitNeighbor, &min_sine, &bAmbiguous, vMinSine );

        /*
         *  Check for tetrahedral ambiguity -- leave it out for now
         */
        if (fabs( triple_product ) > ZERO_FLOAT && ( min_sine > vMinSine || (fabs( min_sine ) > ZERO_FLOAT && ( n2DTetrahedralAmbiguity & T2D_OKAY )) )) /* djb-rwth: addressing LLVM warning */
        {
             /* Even => sorted in correct order, Odd=>transposed */
            parity = triple_product > 0.0 ? AB_PARITY_EVEN : AB_PARITY_ODD;
            /* if ( num_explicit_H && at[cur_at].removed_H_parity % 2 )  */
                  /* odd transposition of the removed implicit H */
            /*     out_at[cur_at].parity = 3 - out_at[cur_at].parity; */

            /*  moved; see below */
            /* out_at[cur_at].bAmbiguousStereo |= bAmbiguous; */
            /* at[cur_at].bAmbiguousStereo |= bAmbiguous; */

            /*  for 4 attached atoms, moving the implicit H from index=3 to index=0 */
            /*  can be done in odd number (3) transpositions: (23)(12)(01), which inverts the parity */
            if (j1 == MAX_NUM_STEREO_ATOM_NEIGH - 1)
            {
                parity = 3 - parity;
            }
        }
        else
        {
    #if (ALWAYS_SET_STEREO_PARITY == 1)
            parity = AT_PARITY_EVEN; /*  suppose atoms are pre-sorted (testing) */
    #else
            if (num_z > 0)
            {
                bAmbiguous |= AMBIGUOUS_STEREO;
            }
            parity = AB_PARITY_UNDF; /*  no parity: 4 bonds are in one plane. */
    #endif
        }

    exit_function:

        if (parity)
        {
            out_at[cur_at].bAmbiguousStereo |= bAmbiguous;
            at[cur_at].bAmbiguousStereo |= bAmbiguous;
        }

        /*  Non-isotopic parity */
        if (at[cur_at].num_H > 1 || parity <= 0)
        {
            ; /*  no non-isotopic parity */
        }
        else
        {
            out_at[cur_at].parity = parity;
        }

        /*  isotopic parity */
        if (parity == -AB_PARITY_UNDF || parity == -vABParityUnknown /*AB_PARITY_UNKN*/)
        {
            parity = -parity;
        }
        if (parity < 0)
        {
            parity = AB_PARITY_NONE;
        }
        out_at[cur_at].parity2 = parity;

        parity = PARITY_VAL( out_at[cur_at].parity );
        out_at[cur_at].stereo_atom_parity = ATOM_PARITY_WELL_DEF( parity ) ? AB_PARITY_CALC : parity;
        parity = PARITY_VAL( out_at[cur_at].parity2 );
        out_at[cur_at].stereo_atom_parity2 = ATOM_PARITY_WELL_DEF( parity ) ? AB_PARITY_CALC : parity;
        /*
        out_at[cur_at].parity2 = out_at[cur_at].parity; // save for stereo + isotopic canon.
        if ( out_at[cur_at].parity ) {
            if ( num_explicit_H > 1 || j1 == MAX_NUM_STEREO_ATOM_NEIGH-1 && num_explicit_H ) {
                //              X   H      X
                // for example,  >C<   or   >C-D
                //              Y   D      Y
                // parity exists for stereo + isotopic atoms canonicalization only
                out_at[cur_at].parity  = 0;
            }
        }
        returning 0 means this can be an adjacent to a stereogenic bond atom
        */

        return (int) out_at[cur_at].parity2;
    }
    #undef ADD_EXPLICIT_HYDROGEN_NEIGH
    #undef ADD_EXPLICIT_LONE_PAIR_NEIGH


        */
    // END INCHI C FUNCTION: set_stereo_atom_parity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: set_stereo_atom_parity
    // INCHI✔️❌: NEW_STEREOCENTER_CHECK == 1; the legacy bCanAtomBeAStereoCenter branch is inactive.
    // INCHI✔️❌: STEREO_CENTER_BONDS_NORM == 1 and STEREO_CENTER_BOND4_NORM == 0.
    // INCHI✔️❌: APPLY_IMPLICIT_H_DOWN_RULE == 0 and ALWAYS_SET_STEREO_PARITY == 0.
    // INCHI✔️❌: FIX_STEREO_SCALING_BUG == 0 selects the runtime PES_BIT_FIX_SP3_BUG branch.
    // END INCHI ACTIVE MACRO CONFIGURATION: set_stereo_atom_parity

    let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(at.as_const())?
        .get(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    {
        let output = heap
            .slice_mut(out_at)?
            .get_mut(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        output.parity = AB_PARITY_NONE as i8;
        output.parity2 = AB_PARITY_NONE as i8;
        output.stereo_atom_parity = AB_PARITY_NONE as i8;
        output.stereo_atom_parity2 = AB_PARITY_NONE as i8;
    }

    let mut parity = AB_PARITY_NONE as i32;
    let mut ambiguous = 0_i32;
    let mut number_of_explicit_hydrogens = 0_i32;
    let mut number_of_explicit_isotopic_hydrogens = [0_i32; NUM_H_ISOTOPES as usize + 1];
    let mut hydrogen_indices = [0_usize; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut neighbor_original_numbers = [0_u16; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut coordinates = [[0.0_f64; 3]; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut bond_lengths_xy = [0.0_f64; MAX_NUM_STEREO_ATOM_NEIGH as usize];
    let mut coordinate_center = [0.0_f64; 3];
    let mut sum_xyz = [0.0_f64; 3];
    let mut number_of_z = 0_i32;
    let mut is_2d = 0_i32;
    let mut add_explicit_neighbor = 0_i32;
    let mut tetrahedral_2d_ambiguity = 0_i32;
    let mut minimum_sine = MIN_SINE;
    let mut minimum_angle = MIN_ANGLE;
    let ignore_isotopic_hydrogen = (i32::from(current.cFlags) & AT_FLAG_ISO_H_POINT as i32) != 0;
    let source_min = |left: f64, right: f64| if left < right { left } else { right };
    let source_max = |left: f64, right: f64| if left > right { left } else { right };

    if LooseTSACheck != 0 && current.nNumAtInRingSystem >= 3 {
        minimum_angle = MIN_ANGLE_RELAXED;
        minimum_sine = MIN_SINE_RELAXED;
    }

    'calculate: {
        let required_number_of_neighbors =
            bCanInpAtomBeAStereoCenter(heap, at, cur_at, bPointedEdgeStereo, bStereoAtZz)?;
        if required_number_of_neighbors == 0 || i32::from(current.num_H) > NUM_H_ISOTOPES as i32 {
            break 'calculate;
        }

        let mut total_number_of_isotopic_hydrogens = 0_i32;
        for isotope in 0..NUM_H_ISOTOPES as usize {
            if current.num_iso_H[isotope] > 1 {
                break 'calculate;
            }
            total_number_of_isotopic_hydrogens += i32::from(current.num_iso_H[isotope]);
        }
        if ignore_isotopic_hydrogen {
            total_number_of_isotopic_hydrogens = 0;
        }
        if i32::from(current.num_H) - total_number_of_isotopic_hydrogens > 1 {
            break 'calculate;
        }

        if !at_removed_H.is_null() && num_removed_H > 0 {
            let removed_atoms = heap.slice(at_removed_H)?;
            for j in 0..num_removed_H {
                let removed = removed_atoms
                    .get(usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                if i32::from(removed.neighbor[0]) == cur_at {
                    let mut isotope = i32::from(removed.iso_atw_diff);
                    if !(0..=NUM_H_ISOTOPES as i32).contains(&isotope) || ignore_isotopic_hydrogen {
                        isotope = 0;
                    }
                    number_of_explicit_isotopic_hydrogens[isotope as usize] += 1;
                    let slot = usize::try_from(number_of_explicit_hydrogens)
                        .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    let hydrogen_slot = hydrogen_indices
                        .get_mut(slot)
                        .ok_or(SourceHeapError::PointerOutOfBounds)?;
                    *hydrogen_slot =
                        usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                    number_of_explicit_hydrogens += 1;
                }
            }
        }

        let mut explicit_neighbor_count = 0_usize;
        if number_of_explicit_hydrogens > 0 {
            let removed_atoms = heap.slice(at_removed_H)?.to_vec();
            let removed_offset = at_removed_H.as_mut().difference(at)?;
            for j in 0..usize::try_from(number_of_explicit_hydrogens)
                .map_err(|_| SourceHeapError::PointerOutOfBounds)?
            {
                let next = hydrogen_indices[j];
                let removed = removed_atoms
                    .get(next)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let removed_index = removed_offset
                    .checked_add(
                        i64::try_from(next)
                            .map_err(|_| SourceHeapError::PointerDifferenceOverflow)?,
                    )
                    .and_then(|value| i32::try_from(value).ok())
                    .ok_or(SourceHeapError::PointerDifferenceOverflow)?;
                let mut coordinate_type = ZTYPE_NONE as i32;
                let z = -get_z_coord(
                    heap,
                    at.as_const(),
                    removed_index,
                    0,
                    &mut coordinate_type,
                    -(bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO as i32),
                )?;
                match coordinate_type {
                    value if value == ZTYPE_EITHER as i32 => {
                        parity = vABParityUnknown;
                        break 'calculate;
                    }
                    value if value == ZTYPE_UP as i32 || value == ZTYPE_DOWN => {
                        coordinate_type = -coordinate_type;
                        is_2d += 1;
                        number_of_z += 1;
                    }
                    value if value == ZTYPE_3D as i32 => number_of_z += 1,
                    _ => {}
                }
                if explicit_neighbor_count >= MAX_NUM_STEREO_ATOM_NEIGH as usize {
                    break;
                }
                neighbor_original_numbers[explicit_neighbor_count] = removed.orig_at_number;
                coordinates[explicit_neighbor_count][0] = removed.x - current.x;
                coordinates[explicit_neighbor_count][1] = removed.y - current.y;
                bond_lengths_xy[explicit_neighbor_count] = len2(&[
                    coordinates[explicit_neighbor_count][0],
                    coordinates[explicit_neighbor_count][1],
                ]);
                coordinates[explicit_neighbor_count][2] = if coordinate_type == ZTYPE_3D as i32 {
                    z
                } else if coordinate_type == ZTYPE_UP as i32 {
                    bond_lengths_xy[explicit_neighbor_count]
                } else if coordinate_type == ZTYPE_DOWN {
                    -bond_lengths_xy[explicit_neighbor_count]
                } else {
                    0.0
                };
                explicit_neighbor_count += 1;
            }
        }

        for j in
            0..usize::try_from(current.valence).map_err(|_| SourceHeapError::PointerOutOfBounds)?
        {
            let next = usize::from(current.neighbor[j]);
            let neighbor = heap
                .slice(at.as_const())?
                .get(next)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let mut coordinate_type = ZTYPE_NONE as i32;
            let z = get_z_coord(
                heap,
                at.as_const(),
                cur_at,
                i32::try_from(j).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                &mut coordinate_type,
                bPointedEdgeStereo & PES_BIT_POINT_EDGE_STEREO as i32,
            )?;
            match coordinate_type {
                value if value == ZTYPE_EITHER as i32 => {
                    parity = vABParityUnknown;
                    break 'calculate;
                }
                value if value == ZTYPE_UP as i32 || value == ZTYPE_DOWN => {
                    is_2d += 1;
                    number_of_z += 1;
                }
                value if value == ZTYPE_3D as i32 => number_of_z += 1,
                _ => {}
            }
            if explicit_neighbor_count >= MAX_NUM_STEREO_ATOM_NEIGH as usize {
                break;
            }
            neighbor_original_numbers[explicit_neighbor_count] = neighbor.orig_at_number;
            coordinates[explicit_neighbor_count][0] = neighbor.x - current.x;
            coordinates[explicit_neighbor_count][1] = neighbor.y - current.y;
            bond_lengths_xy[explicit_neighbor_count] = len2(&[
                coordinates[explicit_neighbor_count][0],
                coordinates[explicit_neighbor_count][1],
            ]);
            coordinates[explicit_neighbor_count][2] = if coordinate_type == ZTYPE_3D as i32 {
                z
            } else if coordinate_type == ZTYPE_UP as i32 {
                bond_lengths_xy[explicit_neighbor_count]
            } else if coordinate_type == ZTYPE_DOWN {
                -bond_lengths_xy[explicit_neighbor_count]
            } else {
                0.0
            };
            explicit_neighbor_count += 1;
        }

        is_2d = i32::from(is_2d == number_of_z && number_of_z != 0);
        let explicit_neighbor_count_i32 = i32::try_from(explicit_neighbor_count)
            .map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        if MAX_NUM_STEREO_ATOM_NEIGH as i32
            != i32::from(current.valence) + number_of_explicit_hydrogens
            && MAX_NUM_STEREO_ATOM_NEIGH as i32 - 1
                != i32::from(current.valence) + number_of_explicit_hydrogens
        {
            if required_number_of_neighbors == i32::from(current.valence) + i32::from(current.num_H)
                && current.num_H > 1
            {
                parity = if parity == vABParityUnknown {
                    -vABParityUnknown
                } else {
                    -(AB_PARITY_UNDF as i32)
                };
            } else {
                parity = AB_PARITY_NONE as i32;
            }
            break 'calculate;
        }

        let mut maximum_radius = 0.0_f64;
        let mut minimum_radius = 0.0_f64;
        for j in 0..explicit_neighbor_count {
            let length = len3(&coordinates[j]);
            if length < MIN_BOND_LEN {
                parity = GetStereocenter0DParity(
                    pCG,
                    heap,
                    at,
                    cur_at,
                    explicit_neighbor_count_i32,
                    &mut neighbor_original_numbers,
                    FlagSC_0D as i32,
                )?;
                if parity == AB_PARITY_NONE as i32 {
                    parity = AB_PARITY_UNDF as i32;
                }
                break 'calculate;
            }
            let coordinate = coordinates[j];
            mult3(&coordinate, 1.0 / length, &mut coordinates[j]);
            maximum_radius = if j != 0 {
                source_max(maximum_radius, length)
            } else {
                length
            };
            minimum_radius = if j != 0 {
                source_min(minimum_radius, length)
            } else {
                length
            };
        }
        if minimum_radius / maximum_radius < MIN_SINE {
            parity = GetStereocenter0DParity(
                pCG,
                heap,
                at,
                cur_at,
                explicit_neighbor_count_i32,
                &mut neighbor_original_numbers,
                FlagSC_0D as i32,
            )?;
            if parity == AB_PARITY_NONE as i32 {
                parity = AB_PARITY_UNDF as i32;
            }
            break 'calculate;
        }
        for coordinate in coordinates.iter().take(explicit_neighbor_count) {
            let previous_sum = sum_xyz;
            add3(&previous_sum, coordinate, &mut sum_xyz);
        }

        if explicit_neighbor_count == MAX_NUM_STEREO_ATOM_NEIGH as usize - 1 {
            if required_number_of_neighbors == MAX_NUM_STEREO_ATOM_NEIGH as i32 - 1 {
                add_explicit_neighbor = 2;
            } else if required_number_of_neighbors == MAX_NUM_STEREO_ATOM_NEIGH as i32
                && i32::from(current.num_H) - number_of_explicit_hydrogens == 1
                && i32::from(current.num_H) - total_number_of_isotopic_hydrogens == 1
            {
                add_explicit_neighbor = 1;
            }
        }

        if add_explicit_neighbor != 0 {
            if parity == vABParityUnknown {
                break 'calculate;
            }
            if number_of_z == 0
                || are_3_vect_in_one_plane(
                    &[coordinates[0], coordinates[1], coordinates[2]],
                    minimum_sine,
                ) != 0
            {
                if number_of_z > 0 {
                    ambiguous |= AMBIGUOUS_STEREO as i32;
                }
                parity = GetStereocenter0DParity(
                    pCG,
                    heap,
                    at,
                    cur_at,
                    explicit_neighbor_count_i32,
                    &mut neighbor_original_numbers,
                    FlagSC_0D as i32,
                )?;
                if parity == AB_PARITY_NONE as i32 {
                    parity = AB_PARITY_UNDF as i32;
                }
                break 'calculate;
            }

            coordinates[explicit_neighbor_count] = sum_xyz;
            let coordinate = coordinates[explicit_neighbor_count];
            change_sign3(&coordinate, &mut coordinates[explicit_neighbor_count]);
            let length = len3(&coordinates[explicit_neighbor_count]);
            if bPointedEdgeStereo & PES_BIT_FIX_SP3_BUG as i32 != 0 {
                if length > 1.0 {
                    maximum_radius *= length;
                } else {
                    minimum_radius *= length;
                }
            } else {
                maximum_radius = source_max(maximum_radius, length);
                minimum_radius = source_min(minimum_radius, length);
            }
            if length < MIN_BOND_LEN || minimum_radius / maximum_radius < MIN_SINE {
                parity = GetStereocenter0DParity(
                    pCG,
                    heap,
                    at,
                    cur_at,
                    explicit_neighbor_count_i32,
                    &mut neighbor_original_numbers,
                    FlagSC_0D as i32,
                )?;
                if parity == AB_PARITY_NONE as i32 {
                    parity = AB_PARITY_UNDF as i32;
                }
                break 'calculate;
            }
        } else if explicit_neighbor_count != MAX_NUM_STEREO_ATOM_NEIGH as usize {
            if parity == vABParityUnknown {
                parity = -(AB_PARITY_UNDF as i32);
            }
            break 'calculate;
        } else if number_of_z == 0 || are_4at_in_one_plane(&coordinates, minimum_sine) != 0 {
            if number_of_z > 0 {
                ambiguous |= AMBIGUOUS_STEREO as i32;
            }
            if parity != vABParityUnknown {
                parity = GetStereocenter0DParity(
                    pCG,
                    heap,
                    at,
                    cur_at,
                    explicit_neighbor_count_i32,
                    &mut neighbor_original_numbers,
                    FlagSC_0D as i32,
                )?;
                if parity == AB_PARITY_NONE as i32 {
                    parity = AB_PARITY_UNDF as i32;
                } else if (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
                    .contains(&parity)
                {
                    ambiguous &= !(AMBIGUOUS_STEREO as i32);
                }
            }
            break 'calculate;
        }

        if is_2d != 0 {
            tetrahedral_2d_ambiguity = Get2DTetrahedralAmbiguity(
                pCG,
                &coordinates,
                add_explicit_neighbor,
                bPointedEdgeStereo & PES_BIT_FIX_SP3_BUG as i32,
                minimum_angle,
            )?;
            if tetrahedral_2d_ambiguity > 0 {
                if tetrahedral_2d_ambiguity & T2D_WARN as i32 != 0 {
                    ambiguous |= AMBIGUOUS_STEREO as i32;
                }
                if tetrahedral_2d_ambiguity & T2D_UNDF as i32 != 0 {
                    if parity != vABParityUnknown {
                        parity = AB_PARITY_UNDF as i32;
                    }
                    break 'calculate;
                }
            } else if tetrahedral_2d_ambiguity < 0 {
                ambiguous |= AMBIGUOUS_STEREO_ERROR as i32;
                parity = AB_PARITY_UNDF as i32;
                break 'calculate;
            }
        }

        let coordinate_origin = coordinates[0];
        for j in 1..MAX_NUM_STEREO_ATOM_NEIGH as usize {
            let coordinate = coordinates[j];
            diff3(&coordinate, &coordinate_origin, &mut coordinates[j]);
        }
        let center = coordinate_center;
        diff3(&center, &coordinate_origin, &mut coordinate_center);

        let tetrahedral_coordinates = [coordinates[1], coordinates[2], coordinates[3]];
        let mut actual_minimum_sine = 0.0_f64;
        let triple_product = triple_prod_and_min_abs_sine2(
            &tetrahedral_coordinates,
            &coordinate_center,
            add_explicit_neighbor,
            Some(&mut actual_minimum_sine),
            Some(&mut ambiguous),
            minimum_sine,
        );
        if triple_product.abs() > 1.0e-12
            && (actual_minimum_sine > minimum_sine
                || (actual_minimum_sine.abs() > 1.0e-12
                    && tetrahedral_2d_ambiguity & T2D_OKAY as i32 != 0))
        {
            parity = if triple_product > 0.0 {
                AB_PARITY_EVEN as i32
            } else {
                AB_PARITY_ODD as i32
            };
            if explicit_neighbor_count == MAX_NUM_STEREO_ATOM_NEIGH as usize - 1 {
                parity = 3_i32.wrapping_sub(parity);
            }
        } else {
            if number_of_z > 0 {
                ambiguous |= AMBIGUOUS_STEREO as i32;
            }
            parity = AB_PARITY_UNDF as i32;
        }
    }

    if parity != 0 {
        {
            let output = heap
                .slice_mut(out_at)?
                .get_mut(current_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            output.bAmbiguousStereo = (i32::from(output.bAmbiguousStereo) | ambiguous) as i8;
        }
        let atoms = heap.slice_mut(at)?;
        let current_atom = atoms
            .get_mut(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        current_atom.bAmbiguousStereo =
            (i32::from(current_atom.bAmbiguousStereo) | ambiguous) as i8;
    }

    {
        let output = heap
            .slice_mut(out_at)?
            .get_mut(current_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        if current.num_H <= 1 && parity > 0 {
            output.parity = parity as i8;
        }
        if parity == -(AB_PARITY_UNDF as i32) || parity == -vABParityUnknown {
            parity = -parity;
        }
        if parity < 0 {
            parity = AB_PARITY_NONE as i32;
        }
        output.parity2 = parity as i8;

        let nonisotopic_parity = i32::from(output.parity) & BITS_PARITY as i32;
        output.stereo_atom_parity = if (AB_MIN_WELL_DEFINED_PARITY as i32
            ..=AB_MAX_WELL_DEFINED_PARITY as i32)
            .contains(&nonisotopic_parity)
        {
            AB_PARITY_CALC as i8
        } else {
            nonisotopic_parity as i8
        };
        let isotopic_parity = i32::from(output.parity2) & BITS_PARITY as i32;
        output.stereo_atom_parity2 = if (AB_MIN_WELL_DEFINED_PARITY as i32
            ..=AB_MAX_WELL_DEFINED_PARITY as i32)
            .contains(&isotopic_parity)
        {
            AB_PARITY_CALC as i8
        } else {
            isotopic_parity as i8
        };
        Ok(i32::from(output.parity2))
    }
}

#[allow(non_snake_case)]
pub(crate) fn bAtomHasValence3(
    elname: &[i8],
    charge: i8,
    radical: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1574 bAtomHasValence3
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-table scan and allocation behavior match.
    /*
    int bAtomHasValence3( char *elname, S_CHAR charge, S_CHAR radical )
    {
        static const char   szElem[][3] = { "N\000" };
        static const S_CHAR   cCharge[] = { 0, };
        int i, ret = 0;
        for (i = 0; i < (int) ( sizeof( szElem ) / sizeof( szElem[0] ) ); i++)
        {
            if (!strcmp( elname, szElem[i] ) && ( charge == cCharge[i] ))
            {
                ret = ( !radical || radical == RADICAL_SINGLET );
                break;
            }
        }
        return ret;
    }
    */
    // END INCHI C FUNCTION: bAtomHasValence3

    let length = elname
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let element = &elname[..length];
    let mut ret = 0_i32;
    if element == [b'N' as i8] && charge == 0 {
        ret = i32::from(radical == 0 || radical == RADICAL_SINGLET as i8);
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn bCanAtomHaveAStereoBond(
    elname: &[i8],
    charge: i8,
    radical: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1594 bCanAtomHaveAStereoBond
    // INCHI✔️✔️: complete source frame follows verbatim; fixed-table scan and allocation behavior match.
    /*
    int bCanAtomHaveAStereoBond( char *elname, S_CHAR charge, S_CHAR radical )
    {
        static const char   szElem[][3] = { "C\000", "Si", "Ge", "N\000", "N\000" };
        static const S_CHAR   cCharge[] = { 0,        0,    0,   0,       1, };
        static const int       n = sizeof( szElem ) / sizeof( szElem[0] );
        int i, ret = 0;
        for (i = 0; i < n; i++)
        {
            if (!strcmp( elname, szElem[i] ) && ( charge == cCharge[i] ))
            {
                ret = ( !radical || radical == RADICAL_SINGLET );
                break;
            }
        }
        return ret;
    }
    */
    // END INCHI C FUNCTION: bCanAtomHaveAStereoBond

    let length = elname
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let element = &elname[..length];
    let elements: [(&[i8], i8); 5] = [
        (&[b'C' as i8], 0),
        (&[b'S' as i8, b'i' as i8], 0),
        (&[b'G' as i8, b'e' as i8], 0),
        (&[b'N' as i8], 0),
        (&[b'N' as i8], 1),
    ];
    let mut ret = 0_i32;
    for (expected_element, expected_charge) in elements {
        if element == expected_element && charge == expected_charge {
            ret = i32::from(radical == 0 || radical == RADICAL_SINGLET as i8);
            break;
        }
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn bCanAtomBeMiddleAllene(
    elname: &[i8],
    charge: i8,
    radical: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1615 bCanAtomBeMiddleAllene
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-table scan and allocation behavior match.
    /*
    int bCanAtomBeMiddleAllene( char *elname, S_CHAR charge, S_CHAR radical )
    {
        static const char   szElem[][3] = { "C\000", "Si", "Ge", };
        static const S_CHAR   cCharge[] = { 0,        0,    0, };
        static const int       n = sizeof( szElem ) / sizeof( szElem[0] );
        int i, ret = 0;
        for (i = 0; i < n; i++)
        {
            if (!strcmp( elname, szElem[i] ) && ( charge == cCharge[i] ))
            {
                ret = ( !radical || radical == RADICAL_SINGLET );
                break;
            }
        }
        return ret;
    }
    */
    // END INCHI C FUNCTION: bCanAtomBeMiddleAllene

    let length = elname
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let element = &elname[..length];
    let elements: [&[i8]; 3] = [
        &[b'C' as i8],
        &[b'S' as i8, b'i' as i8],
        &[b'G' as i8, b'e' as i8],
    ];
    let mut ret = 0_i32;
    for expected_element in elements {
        if element == expected_element && charge == 0 {
            ret = i32::from(radical == 0 || radical == RADICAL_SINGLET as i8);
            break;
        }
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn bCanAtomBeTerminalAllene(
    elname: &[i8],
    charge: i8,
    radical: i8,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1705 bCanAtomBeTerminalAllene
    // INCHI✔️✔️: complete active source frame follows verbatim; fixed-table scan and allocation behavior match.
    /*
    int bCanAtomBeTerminalAllene( char *elname, S_CHAR charge, S_CHAR radical )
    {
        static const char   szElem[][3] = { "C\000", "Si", "Ge", };
        static const S_CHAR   cCharge[] = { 0,        0,    0, };
        static const int       n = sizeof( szElem ) / sizeof( szElem[0] );
        int i, ret = 0;
        for (i = 0; i < n; i++)
        {
            if (!strcmp( elname, szElem[i] ) && ( charge == cCharge[i] ))
            {
                ret = ( !radical || radical == RADICAL_SINGLET );
                break;
            }
        }

        return ret;
    }
    */
    // END INCHI C FUNCTION: bCanAtomBeTerminalAllene

    let length = elname
        .iter()
        .position(|&byte| byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let element = &elname[..length];
    let elements: [&[i8]; 3] = [
        &[b'C' as i8],
        &[b'S' as i8, b'i' as i8],
        &[b'G' as i8, b'e' as i8],
    ];
    let mut ret = 0_i32;
    for expected_element in elements {
        if element == expected_element && charge == 0 {
            ret = i32::from(radical == 0 || radical == RADICAL_SINGLET as i8);
            break;
        }
    }
    Ok(ret)
}

#[allow(non_snake_case)]
pub(crate) fn bIsSuitableHeteroInpAtom(
    heap: &SourceHeap,
    at: SourceConstPointer<inp_ATOM>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1634 bIsSuitableHeteroInpAtom
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int bIsSuitableHeteroInpAtom( inp_ATOM  *at )
    {
        int val, num_H;
        if (0 == at->charge &&
            ( !at->radical || RADICAL_SINGLET == at->radical ) &&
             0 < ( val = get_endpoint_valence( at->el_number ) ))
        {
            num_H = at->num_H;
            if (val == at->chem_bonds_valence + num_H)
            {
                switch (val)
                {
                    case 2: /* O */
                        if (!num_H && 1 == at->valence)
                        {
                            return 0; /* =O */
                        }
                        break;        /* not found */
                    case 3: /* N */
                        if ((1 == at->valence && 1 == num_H) ||
                             (2 == at->valence && 0 == num_H)) /* djb-rwth: addressing LLVM warnings */
                        {
                            return 1; /* =N- or =NH */
                        }
                        break;        /* not found */
                }
            }
        }
        return -1;
    }
    */
    // END INCHI C FUNCTION: bIsSuitableHeteroInpAtom

    let atom = heap
        .slice(at)?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if atom.charge == 0 && (atom.radical == 0 || atom.radical == RADICAL_SINGLET as i8) {
        let val = get_endpoint_valence(atom.el_number);
        if val > 0 {
            let num_h = i32::from(atom.num_H);
            if val == i32::from(atom.chem_bonds_valence) + num_h {
                match val {
                    2 if num_h == 0 && atom.valence == 1 => return Ok(0),
                    3 if (atom.valence == 1 && num_h == 1) || (atom.valence == 2 && num_h == 0) => {
                        return Ok(1);
                    }
                    _ => {}
                }
            }
        }
    }
    Ok(-1)
}

#[allow(non_snake_case)]
pub(crate) fn bIsOxide(
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1667 bIsOxide
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int bIsOxide( inp_ATOM  *at, int cur_at )
    {
        int i, bond_type;
        inp_ATOM  *a = at + cur_at, *an;
        for (i = 0; i < a->valence; i++)
        {
            bond_type = ( a->bond_type[i] &= ~BOND_MARK_ALL );
            if (bond_type == BOND_DOUBLE)
            {
                an = at + (int) a->neighbor[i];
                if (1 == an->valence &&
                     !an->charge && !an->num_H && !an->radical &&
                     2 == get_endpoint_valence( an->el_number ))
                {
                    return 1;
                }
            }
            else
            {
                if (bond_type == BOND_TAUTOM || bond_type == BOND_ALT12NS)
                {
                    an = at + (int) a->neighbor[i];
                    if (1 == an->valence &&
                         2 == get_endpoint_valence( an->el_number ))
                    {
                        return 1;
                    }
                }
            }
        }

        return 0;
    }
    */
    // END INCHI C FUNCTION: bIsOxide

    let current_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let valence = heap
        .slice(at.as_const())?
        .get(current_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .valence;
    for i in 0..i32::from(valence) {
        let bond_index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let (bond_type, neighbor) = {
            let atom = heap
                .slice_mut(at)?
                .get_mut(current_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let bond_type = atom
                .bond_type
                .get_mut(bond_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            *bond_type &= !(BOND_MARK_ALL as u8);
            (*bond_type, atom.neighbor[bond_index])
        };
        if bond_type == BOND_TYPE_DOUBLE as u8 {
            let neighbor = heap
                .slice(at.as_const())?
                .get(usize::from(neighbor))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.valence == 1
                && neighbor.charge == 0
                && neighbor.num_H == 0
                && neighbor.radical == 0
                && get_endpoint_valence(neighbor.el_number) == 2
            {
                return Ok(1);
            }
        } else if bond_type == BOND_TAUTOM as u8 || bond_type == BOND_ALT12NS as u8 {
            let neighbor = heap
                .slice(at.as_const())?
                .get(usize::from(neighbor))
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if neighbor.valence == 1 && get_endpoint_valence(neighbor.el_number) == 2 {
                return Ok(1);
            }
        }
    }

    Ok(0)
}

pub(crate) fn get_allowed_stereo_bond_type(mut bond_type: i32) -> i32 {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2588 get_allowed_stereo_bond_type
    // INCHI✔️✔️: complete source frame follows verbatim; active mode.h macro branches are reproduced below.
    /*
    int get_allowed_stereo_bond_type( int bond_type )
    {
    #if (ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS == 0 )
        if (( bond_type & ~BOND_MARK_ALL ) == BOND_TAUTOM)
        {
            return 0; /*  no tautomer bonds allowed */
        }
        else
    #endif
    #if ( EXCL_ALL_AROM_BOND_PARITY  == 1 )  /* { */
        /*  a stereo bond cannot belong to an aromatic atom */
            if (( bond_type &= ~BOND_MARK_ALL ) == BOND_ALTERN)
            {
                return 0;
            }
    #else  /* } { */
    #if ( ADD_6MEMB_AROM_BOND_PARITY == 1 )
        /*  accept any aromatic bond as a stereo bond */
        if (( bond_type &= ~BOND_MARK_ALL ) == BOND_ALTERN)
    #else
        /*  accept only aromatic bonds in non-6-member rings */
        if (( bond_type &= ~BOND_MARK_ALL ) == BOND_ALTERN) )
    #endif
        {
            return BOND_ALTERN;
        }
    #endif  /* } */
        else
        {
            /*  at this point BOND_MARK_ALL bits have been removed from bond_type */
            if (bond_type == BOND_DOUBLE || bond_type == BOND_SINGLE)
            {
                return bond_type;
            }
    #if (ALLOW_TAUT_ATTACHMENTS_TO_STEREO_BONDS == 1 )
            else
            {
                if (bond_type == BOND_TAUTOM)
                {
                    return BOND_TAUTOM;
                }
            }
    #endif
        }

        return 0; /*  wrong bond type */
    }
    */
    // END INCHI C FUNCTION: get_allowed_stereo_bond_type

    bond_type &= !(BOND_MARK_ALL as i32);
    if bond_type == BOND_ALTERN as i32 {
        BOND_ALTERN as i32
    } else if bond_type == BOND_TYPE_DOUBLE as i32 || bond_type == BOND_SINGLE as i32 {
        bond_type
    } else if bond_type == BOND_TAUTOM as i32 {
        BOND_TAUTOM as i32
    } else {
        0
    }
}

pub(crate) fn can_be_a_stereo_bond_with_isotopic_H(
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    cur_at: i32,
    nMode: INCHI_MODE,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2638 can_be_a_stereo_bond_with_isotopic_H
    // INCHI✔️❌: complete source frame follows verbatim; checked SourceHeap access adds overhead.
    /*
    int can_be_a_stereo_bond_with_isotopic_H( inp_ATOM *at,
                                              int cur_at,
                                              INCHI_MODE nMode )
    {
        int i, j, next_at, num_stereo_bonds, bFound;
        int bond_type, num_2s, num_alt;
        int num_2s_next, num_alt_next, num_wrong_bonds_1, num_wrong_bonds_2;
    #if ( N_V_STEREOBONDS == 1 )
        int n2sh, num_2s_hetero[2], num_2s_hetero_next[2], next_next_at, type_N, type_N_next;
    #endif

        if (MAX_NUM_STEREO_BOND_NEIGH < at[cur_at].valence + at[cur_at].num_H ||
             MIN_NUM_STEREO_BOND_NEIGH > at[cur_at].valence + at[cur_at].num_H)
        {
            return 0;
        }
        if (!bCanAtomHaveAStereoBond( at[cur_at].elname, at[cur_at].charge, at[cur_at].radical ))
        {
            return 0;
        }

        /*  Count bonds and find the second atom on the stereo bond */
        num_2s = num_alt = num_wrong_bonds_1 = 0;
    #if ( N_V_STEREOBONDS == 1 )
        num_2s_hetero[0] = num_2s_hetero[1] = type_N = 0;
        if (0 == at[cur_at].num_H && 0 == at[cur_at].charge && 0 == at[cur_at].radical &&
             3 == get_endpoint_valence( at[cur_at].el_number ))
        {
            if (2 == at[cur_at].valence && 3 == at[cur_at].chem_bonds_valence)
            {
                type_N = 1;
            }
            else
            {
                if (3 == at[cur_at].valence && 5 == at[cur_at].chem_bonds_valence)
                {
                    type_N = 2; /* unfortunately includes >N# */
                }
            }
        }
    #endif
        for (i = 0, num_stereo_bonds = 0; i < at[cur_at].valence; i++)
        {
            bFound = 0;
            next_at = at[cur_at].neighbor[i];
            bond_type = get_allowed_stereo_bond_type( (int) at[cur_at].bond_type[i] );
            if (bond_type == BOND_ALTERN)
            {
                num_alt++;
                if (cur_at > next_at && !( nMode & CMODE_NO_ALT_SBONDS ))
                {
                    bFound = 1;
                }
            }
            else
            {
                if (bond_type == BOND_DOUBLE)
                {
                    num_2s++;
    #if ( N_V_STEREOBONDS == 1 )
                    if (0 <= ( n2sh = bIsSuitableHeteroInpAtom( at + next_at ) ))
                    {
                        num_2s_hetero[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
                    }
    #endif
                    if (cur_at > next_at)
                        bFound = 1;
                }
                else
                {
                    if (bond_type != BOND_SINGLE && bond_type != BOND_TAUTOM)
                    {
                        num_wrong_bonds_1++;
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
                        if (num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 && 2 >= at[cur_at].valence)) /* djb-rwth: addressing LLVM warning */
                        {
                            return 0; /* wrong bond type */
                        }
                        else
                        {
                            continue;
                        }
    #else
                        return 0; /*  wrong bond type */
    #endif
                    }
                }
            }

            if (bFound)
            {
                /*  check "next_at" atom on the opposite side of the bond */
                if (MAX_NUM_STEREO_BOND_NEIGH < at[next_at].valence + at[next_at].num_H ||
                     MIN_NUM_STEREO_BOND_NEIGH > at[next_at].valence + at[next_at].num_H)
                {
                    continue;
                }
                if (!bCanAtomHaveAStereoBond( at[next_at].elname, at[next_at].charge, at[next_at].radical ))
                {
                    continue;
                }
                /*  next atom neighbors */
                num_2s_next = num_alt_next = num_wrong_bonds_2 = 0;
    #if ( N_V_STEREOBONDS == 1 )
                num_2s_hetero_next[0] = num_2s_hetero_next[1] = type_N_next = 0;
                if (0 == at[next_at].num_H && 0 == at[next_at].charge && 0 == at[next_at].radical &&
                     3 == get_endpoint_valence( at[next_at].el_number ))
                {
                    if (2 == at[next_at].valence && 3 == at[next_at].chem_bonds_valence)
                    {
                        type_N_next = 1; /* -N= */
                    }
                    else
                    {
                        if (3 == at[next_at].valence && 5 == at[next_at].chem_bonds_valence)
                        {
                            type_N_next = 2; /* unfortunately includes >N# */
                        }
                    }
                }
    #endif
                for (j = 0; j < at[next_at].valence; j++)
                {
                    bond_type = get_allowed_stereo_bond_type( (int) at[next_at].bond_type[j] );
                    if (bond_type == BOND_ALTERN)
                    {
                        num_alt_next++;
                    }
                    else
                    {
                        if (bond_type == BOND_DOUBLE)
                        {
                            num_2s_next++;
    #if ( N_V_STEREOBONDS == 1 )
                            next_next_at = at[next_at].neighbor[j];
                            if (0 <= ( n2sh = bIsSuitableHeteroInpAtom( at + next_next_at ) ))
                            {
                                num_2s_hetero_next[n2sh] ++; /* n2sh=0 -> =N- or =NH; n2sh=1 -> =O */
                            }
    #endif
                        }
                        else
                        {
                            if (bond_type != BOND_SINGLE && bond_type != BOND_TAUTOM)
                            {
                                num_wrong_bonds_2++;
    #if ( ONE_BAD_SB_NEIGHBOR == 1 )
                                if (num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 && 2 >= at[cur_at].valence)) /* djb-rwth: addressing LLVM warning */
                                {
                                    break; /* wrong bond type */
                                }
                                else
                                {
                                    continue;
                                }
    #else
                                break; /*  wrong bond type */
    #endif
                            }
                        }
                    }
                }
                /* figure out whether the at[cur_at]--at[next_at] bond may not be stereogenic */

    #if ( N_V_STEREOBONDS == 1 )
                if (3 == ( type_N | type_N_next ) &&
                    ( (2 == type_N && !bIsOxide( at, cur_at )) ||
                        (2 == type_N_next && !bIsOxide( at, next_at )) )) /* djb-rwth: addressing LLVM warnings */
                {
                    bFound = 0;
                }
                else
    #endif
                {
                    if (j < at[next_at].valence ||                  /* at[next_at] has a wrong bond type*/
                        ( num_alt_next > 0 ) + ( num_2s_next > 0 ) != 1     /* only one type of stereogenic bond permitted */
                           )
                    {
                        bFound = 0;
                    }
                    else
                    {
                        if (2 < num_2s_next)
                        {
                            bFound = 0;
                        }
                        else
                        {
                            if (2 == num_2s_next)
                            {
                                if (2 == at[next_at].valence)
                                {
                                    ; /* only one double bond permitted except cumulenes */
    #if ( N_V_STEREOBONDS == 1 )
                                }
                                else
                                {
                                    if (1 == ( num_2s_hetero_next[0] | num_2s_hetero_next[1] ) &&
                                            3 == at[next_at].valence + at[next_at].num_H &&
                                            5 == at[next_at].chem_bonds_valence + at[next_at].num_H &&
                                            3 == get_endpoint_valence( at[next_at].el_number ) &&
                                            ( !type_N || bIsOxide( at, next_at ) ))
                                    {
                                        ; /*
                                           *   found:
                                           *
                                           *    \      /    \      /    \      /
                                           *     \    /      \    /      \    /
                                           *      N==C   or   N==C   or   N==N
                                           *    //    \     //    \     //    \
                                           *   O  ^    \   N  ^    \   O  ^    \
                                           *      |           |           |
                                           *      |           |           |
                                           *      at[next_at] at[next_at] at[next_at]
                                           */
    #endif
                                    }
                                    else
                                    {
                                        bFound = 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (bFound)
            {
                num_stereo_bonds++;
            }
        }

        if (( num_alt > 0 ) + ( num_2s > 0 ) != 1 || !num_stereo_bonds)
        {
            return 0;
        }

        if (num_2s > 1)
        {
    #if ( N_V_STEREOBONDS == 1 )
            if (2 == num_2s &&
                 1 == ( num_2s_hetero[0] | num_2s_hetero[1] ) &&
                 3 == at[cur_at].valence + at[cur_at].num_H &&
                 5 == at[cur_at].chem_bonds_valence + at[cur_at].num_H &&
                 3 == get_endpoint_valence( at[cur_at].el_number ))
            {
                ;
            }
            else
            {
                return 0;
            }
    #else
            return 0;
    #endif
        }

        return num_stereo_bonds;
    }
    */
    // END INCHI C FUNCTION: can_be_a_stereo_bond_with_isotopic_H

    let cur_index = usize::try_from(cur_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(at.as_const())?
        .get(cur_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?
        .clone();
    let current_neighbors = i32::from(current.valence) + i32::from(current.num_H);
    if !(2..=3).contains(&current_neighbors)
        || bCanAtomHaveAStereoBond(&current.elname, current.charge, current.radical)? == 0
    {
        return Ok(0);
    }

    let mut num_2s = 0_i32;
    let mut num_alt = 0_i32;
    let mut num_wrong_bonds_1 = 0_i32;
    let mut num_2s_hetero = [0_i32; 2];
    let mut type_n = 0_i32;
    if current.num_H == 0
        && current.charge == 0
        && current.radical == 0
        && get_endpoint_valence(current.el_number) == 3
    {
        if current.valence == 2 && current.chem_bonds_valence == 3 {
            type_n = 1;
        } else if current.valence == 3 && current.chem_bonds_valence == 5 {
            type_n = 2;
        }
    }

    let mut num_stereo_bonds = 0_i32;
    for i in 0..i32::from(current.valence) {
        let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next_at = i32::from(current.neighbor[index]);
        let bond_type = get_allowed_stereo_bond_type(i32::from(current.bond_type[index]));
        let mut found = false;
        if bond_type == BOND_ALTERN as i32 {
            num_alt += 1;
            if cur_at > next_at && nMode & CMODE_NO_ALT_SBONDS as INCHI_MODE == 0 {
                found = true;
            }
        } else if bond_type == BOND_TYPE_DOUBLE as i32 {
            num_2s += 1;
            let suitable =
                bIsSuitableHeteroInpAtom(heap, at.as_const().offset(i64::from(next_at))?)?;
            if suitable >= 0 {
                num_2s_hetero[usize::try_from(suitable)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?] += 1;
            }
            if cur_at > next_at {
                found = true;
            }
        } else if bond_type != BOND_SINGLE as i32 && bond_type != BOND_TAUTOM as i32 {
            num_wrong_bonds_1 += 1;
            if num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 != 0 && current.valence <= 2) {
                return Ok(0);
            }
            continue;
        }

        if found {
            let next_index =
                usize::try_from(next_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let next = heap
                .slice(at.as_const())?
                .get(next_index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?
                .clone();
            let next_neighbors = i32::from(next.valence) + i32::from(next.num_H);
            if !(2..=3).contains(&next_neighbors)
                || bCanAtomHaveAStereoBond(&next.elname, next.charge, next.radical)? == 0
            {
                continue;
            }

            let mut num_2s_next = 0_i32;
            let mut num_alt_next = 0_i32;
            let mut _num_wrong_bonds_2 = 0_i32;
            let mut num_2s_hetero_next = [0_i32; 2];
            let mut type_n_next = 0_i32;
            if next.num_H == 0
                && next.charge == 0
                && next.radical == 0
                && get_endpoint_valence(next.el_number) == 3
            {
                if next.valence == 2 && next.chem_bonds_valence == 3 {
                    type_n_next = 1;
                } else if next.valence == 3 && next.chem_bonds_valence == 5 {
                    type_n_next = 2;
                }
            }

            let mut j = 0_i32;
            while j < i32::from(next.valence) {
                let next_bond_index =
                    usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let next_bond_type =
                    get_allowed_stereo_bond_type(i32::from(next.bond_type[next_bond_index]));
                if next_bond_type == BOND_ALTERN as i32 {
                    num_alt_next += 1;
                } else if next_bond_type == BOND_TYPE_DOUBLE as i32 {
                    num_2s_next += 1;
                    let next_next_at = i32::from(next.neighbor[next_bond_index]);
                    let suitable = bIsSuitableHeteroInpAtom(
                        heap,
                        at.as_const().offset(i64::from(next_next_at))?,
                    )?;
                    if suitable >= 0 {
                        num_2s_hetero_next[usize::try_from(suitable)
                            .map_err(|_| SourceHeapError::PointerOutOfBounds)?] += 1;
                    }
                } else if next_bond_type != BOND_SINGLE as i32
                    && next_bond_type != BOND_TAUTOM as i32
                {
                    _num_wrong_bonds_2 += 1;
                    if num_wrong_bonds_1 > 1 || (num_wrong_bonds_1 != 0 && current.valence <= 2) {
                        break;
                    }
                    j += 1;
                    continue;
                }
                j += 1;
            }

            if (type_n | type_n_next) == 3
                && ((type_n == 2 && bIsOxide(heap, at, cur_at)? == 0)
                    || (type_n_next == 2 && bIsOxide(heap, at, next_at)? == 0))
            {
                found = false;
            } else if j < i32::from(next.valence)
                || i32::from(num_alt_next > 0) + i32::from(num_2s_next > 0) != 1
                || num_2s_next > 2
            {
                found = false;
            } else if num_2s_next == 2 && next.valence != 2 {
                if (num_2s_hetero_next[0] | num_2s_hetero_next[1]) != 1
                    || i32::from(next.valence) + i32::from(next.num_H) != 3
                    || i32::from(next.chem_bonds_valence) + i32::from(next.num_H) != 5
                    || get_endpoint_valence(next.el_number) != 3
                    || (type_n != 0 && bIsOxide(heap, at, next_at)? == 0)
                {
                    found = false;
                }
            }
        }
        if found {
            num_stereo_bonds += 1;
        }
    }

    if i32::from(num_alt > 0) + i32::from(num_2s > 0) != 1 || num_stereo_bonds == 0 {
        return Ok(0);
    }
    if num_2s > 1
        && !(num_2s == 2
            && (num_2s_hetero[0] | num_2s_hetero[1]) == 1
            && i32::from(current.valence) + i32::from(current.num_H) == 3
            && i32::from(current.chem_bonds_valence) + i32::from(current.num_H) == 5
            && get_endpoint_valence(current.el_number) == 3)
    {
        return Ok(0);
    }

    Ok(num_stereo_bonds)
}

#[allow(non_snake_case, clippy::too_many_arguments)]
pub(crate) fn set_stereo_parity(
    pCG: &mut CANON_GLOBALS,
    heap: &mut SourceHeap,
    at: SourceMutPointer<inp_ATOM>,
    at_output: SourceMutPointer<sp_ATOM>,
    num_at: i32,
    num_removed_H: i32,
    mut nMaxNumStereoAtoms: Option<&mut i32>,
    mut nMaxNumStereoBonds: Option<&mut i32>,
    nMode: INCHI_MODE,
    bPointedEdgeStereo: i32,
    vABParityUnknown: i32,
    bLooseTSACheck: i32,
    bStereoAtZz: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4398 set_stereo_parity
    // INCHI✔️❌: complete active source frame follows verbatim; checked SourceHeap access and modeled allocations add overhead.
    /*
    int set_stereo_parity( CANON_GLOBALS *pCG,
                           inp_ATOM* at,
                           sp_ATOM* at_output,
                           int num_at,
                           int num_removed_H,
                           int *nMaxNumStereoAtoms,
                           int *nMaxNumStereoBonds,
                           INCHI_MODE nMode,
                           int bPointedEdgeStereo,
                           int vABParityUnknown,
                           int bLooseTSACheck,
                           int bStereoAtZz )
    {
        int num_3D_stereo_atoms = 0;
        int num_stereo_bonds = 0; /* added to fix allene stereo bug reported for FClC=C=CFCl by Burt Leland - 2009-02-05 DT */

        int i, is_stereo, num_stereo, max_stereo_atoms = 0, max_stereo_bonds = 0;
        QUEUE *q = NULL;
        AT_RANK *nAtomLevel = NULL;
        S_CHAR  *cSource = NULL;
        AT_RANK min_sb_ring_size = 0;

        /**********************************************************
         *
         * Note: this parity reflects only relative positions of
         *       the atoms-neighbors and their ordering in the
         *       lists of neighbors.
         *
         *       To obtain the actual parity, the parity of a number
         *       of neighbors transpositions (to obtain a sorted
         *       list of numbers assigned to the atoms) should be
         *       added.
         *
         **********************************************************/

        /*********************************************************************************

         An example of parity=1 for stereogenic center, tetrahedral asymmetric atom



                      (1)
                       |
                       |
                   [C] |
                       |
             (2)------(0)
                      /
                    /
                  /
                /
             (3)


         Notation: (n) is a tetrahedral atom neighbor; n is an index of a neighbor in
         the central_at->neighbor[] array : neighbor atom number is central_at->neighbor[n].

         (0)-(1), (0)-(2), (0)-(3) are lines connecting atom [C] neighbors to neighbor (0)
         (0), (1) and (2) are in the plane
         (0)-(3) is directed from the plain to the viewer
         [C] is somewhere between (0), (1), (2), (3)
         Since (1)-(2)-(3)  are in a clockwise order when looking from (0), parity is 2, or even;
         otherwise parity would be 1, or odd.

        **********************************************************************************

          Examples of a stereogenic bond.

          Notation:   [atom number], (index of a neighbor):
                      [1] and [2] are atoms connected by the stereogenic bond
                      numbers in () are indexes of neighbors of [1] or [2].
                      (12 x 16)z = z-component of [1]-[2] and [1]-[6] cross-product

                                         atom [1]                     atom [2]
         [8]                    [4]      prod01 = (12 x 16)z < 0      prod01 = (21 x 24)z < 0
           \                    /        prod02 = (12 x 18)z > 0      prod02 = (21 x 25)z > 0
            (2)               (1)        0 transpositions because     0 transpositions because
              \              /           double bond is in 0 posit.   double bond is in 0 position
              [1]==(0)(0)==[2]           0 = (prod01 > prod02)        0 = (prod01 > prod02)
              /              \
            (1)               (2)        result: parity = 2, even     result: parity=2, even
           /                    \
         [6]                    [5]



                                         atom [1]                     atom [2]
         [8]                    [5]      prod01 = (12 x 18)z > 0      prod01 = (21 x 24)z > 0
           \                    /        prod02 = (12 x 16)z < 0      prod02 = (21 x 25)z < 0
            (0)               (2)        2 transpositions to move     1 transposition to move
              \              /           at [2] from 2 to 0 pos.      at [1] from 1 to 0 position
              [1]==(2)(1)==[2]           1 = (prod01 > prod02)        1 = (prod01 > prod02)
              /              \
            (1)               (0)        result: parity = (1+2)       result: parity=(1+1)
           /                    \        2-(1+2)%2 = 1, odd           2-(1+1)%2 = 2, even
         [6]                    [4]


        ***********************************************************************************
        Note: atoms' numbers [1], [2], [4],... are not used to calculate parity at this
              point. They will be used for each numbering in the canonicalization.
        Note: parity=3 for a stereo atom means entered undefined bond direction
              parity=4 for an atom means parity cannot be determined from the given geometry
        ***********************************************************************************/

        if (!at_output || !at)
        {
            return -1;
        }

        /*  Clear stereo descriptors */
        for (i = 0; i < num_at; i++)
        {
            at_output[i].parity = 0;
            at_output[i].parity2 = 0;
            memset( &at_output[i].stereo_bond_neighbor[0], 0, sizeof( at_output[0].stereo_bond_neighbor ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( &at_output[i].stereo_bond_neighbor2[0], 0, sizeof( at_output[0].stereo_bond_neighbor2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( &at_output[i].stereo_bond_ord[0], 0, sizeof( at_output[0].stereo_bond_ord ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( &at_output[i].stereo_bond_ord2[0], 0, sizeof( at_output[0].stereo_bond_ord2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( &at_output[i].stereo_bond_z_prod[0], 0, sizeof( at_output[0].stereo_bond_z_prod ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( &at_output[i].stereo_bond_z_prod2[0], 0, sizeof( at_output[0].stereo_bond_z_prod2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( &at_output[i].stereo_bond_parity[0], 0, sizeof( at_output[0].stereo_bond_parity ) ); /* djb-rwth: memset_s C11/Annex K variant? */
            memset( &at_output[i].stereo_bond_parity2[0], 0, sizeof( at_output[0].stereo_bond_parity2 ) ); /* djb-rwth: memset_s C11/Annex K variant? */
        }

        /*  Estimate max numbers of stereo atoms and bonds if isotopic H are added */
        if (nMaxNumStereoAtoms || nMaxNumStereoBonds)
        {
            for (i = 0; i < num_at; i++) /* djb-rwth: removing redundant code */
            {
                int num;
                num = can_be_a_stereo_atom_with_isotopic_H( at, i, bPointedEdgeStereo, bStereoAtZz );
                if (num)
                {
                    max_stereo_atoms += num;
                }
                else
                {
                    if (( num = can_be_a_stereo_bond_with_isotopic_H( at, i, nMode ) ))
                    {
                        /*  accept cumulenes */
                        max_stereo_bonds += num;
                    }
                }
            }
            if (nMaxNumStereoAtoms)
            {
                *nMaxNumStereoAtoms = max_stereo_atoms;
            }
            if (nMaxNumStereoBonds)
            {
                *nMaxNumStereoBonds = max_stereo_bonds;
            }
        }

        /*  Calculate stereo descriptors */
    #if ( MIN_SB_RING_SIZE > 0 )
        min_sb_ring_size = (AT_RANK) ( ( ( nMode & REQ_MODE_MIN_SB_RING_MASK ) >> REQ_MODE_MIN_SB_RING_SHFT ) & AT_RANK_MASK );
        if (min_sb_ring_size >= 3)
        {
            /* Create BFS data structure for finding for each stereo bond its min. ring sizes */
            q = QueueCreate( num_at + 1, sizeof( qInt ) );
            nAtomLevel = (AT_RANK*) inchi_calloc( num_at, sizeof( nAtomLevel[0] ) );
            cSource = (S_CHAR *) inchi_calloc( num_at, sizeof( cSource[0] ) );
            if (!q || !cSource || !nAtomLevel)
            {
                q = QueueDelete(q); /* djb-rwth: fixing coverity ID #499562 */
                inchi_free(nAtomLevel); /* djb-rwth: avoiding memory leak */
                inchi_free(cSource); /* djb-rwth: avoiding memory leak */
                num_3D_stereo_atoms = CT_OUT_OF_RAM;
                goto exit_function;
            }
        }
        else
        {
            min_sb_ring_size = 2;
        }
    #endif

        /* Main cycle: set stereo parities */
        for (i = 0, num_stereo = 0; i < num_at; i++)
        {
            is_stereo = set_stereo_atom_parity( pCG, at_output, at, i, at + num_at, num_removed_H,
                                                bPointedEdgeStereo, vABParityUnknown, bLooseTSACheck, bStereoAtZz );
            if (is_stereo)
            {
                num_3D_stereo_atoms += ATOM_PARITY_WELL_DEF( is_stereo );
            }
            else
            {
                is_stereo = set_stereo_bonds_parity( at_output, at, i, at + num_at,
                                                     num_removed_H, nMode,q,
                                                     nAtomLevel, cSource,
                                                     min_sb_ring_size,
                                                     bPointedEdgeStereo,
                                                     vABParityUnknown );
                if (RETURNED_ERROR( is_stereo ))
                {
                    num_3D_stereo_atoms = is_stereo;
                    break;
                }
                num_stereo_bonds += ( is_stereo != 0 ); /* added to fix bug reported by Burt Leland - 2009-02-05 DT */
            }
            num_stereo += ( is_stereo != 0 );
            /* djb-rwth: removing redundant code */
        }

        /* Added to fix bug reported by Burt Leland - 2009-02-05 DT */
        if (max_stereo_atoms < num_3D_stereo_atoms && nMaxNumStereoAtoms)
        {
            *nMaxNumStereoAtoms = num_3D_stereo_atoms;
        }
        if (max_stereo_bonds < num_stereo_bonds && nMaxNumStereoBonds)
        {
            *nMaxNumStereoBonds = num_stereo_bonds;
        }

        /*
        if ( (nMode & REQ_MODE_SC_IGN_ALL_UU )
        REQ_MODE_SC_IGN_ALL_UU
        REQ_MODE_SB_IGN_ALL_UU
        */

    #if ( MIN_SB_RING_SIZE > 0 )
        if (q)
        {
            q = QueueDelete( q ); /* djb-rwth: ignoring LLVM warning: variable used to store function return value */
        }
        if (nAtomLevel)
        {
            inchi_free( nAtomLevel );
        }
        if (cSource)
        {
            inchi_free( cSource );
        }

    exit_function:
    #endif

        return num_3D_stereo_atoms;
    }
        */
    // END INCHI C FUNCTION: set_stereo_parity
    // BEGIN INCHI ACTIVE MACRO CONFIGURATION: set_stereo_parity
    // INCHI✔️❌: MIN_SB_RING_SIZE == 8, AT_RANK == uint16_t, and qInt == AT_RANK.
    // INCHI✔️❌: inchi_calloc/inchi_free select the GCC/Linux libc macro behavior.
    // END INCHI ACTIVE MACRO CONFIGURATION: set_stereo_parity

    if at_output.is_null() || at.is_null() {
        return Ok(-1);
    }

    if num_at > 0 {
        let output = heap.slice_mut(at_output)?;
        let count = usize::try_from(num_at).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if output.len() < count {
            return Err(SourceHeapError::PointerOutOfBounds);
        }
        for atom in &mut output[..count] {
            atom.parity = 0;
            atom.parity2 = 0;
            atom.stereo_bond_neighbor.fill(0);
            atom.stereo_bond_neighbor2.fill(0);
            atom.stereo_bond_ord.fill(0);
            atom.stereo_bond_ord2.fill(0);
            atom.stereo_bond_z_prod.fill(0);
            atom.stereo_bond_z_prod2.fill(0);
            atom.stereo_bond_parity.fill(0);
            atom.stereo_bond_parity2.fill(0);
        }
    }

    let mut max_stereo_atoms = 0_i32;
    let mut max_stereo_bonds = 0_i32;
    if nMaxNumStereoAtoms.is_some() || nMaxNumStereoBonds.is_some() {
        for i in 0..num_at {
            let number =
                can_be_a_stereo_atom_with_isotopic_H(heap, at, i, bPointedEdgeStereo, bStereoAtZz)?;
            if number != 0 {
                max_stereo_atoms = max_stereo_atoms.wrapping_add(number);
            } else {
                let number = can_be_a_stereo_bond_with_isotopic_H(heap, at, i, nMode)?;
                if number != 0 {
                    max_stereo_bonds = max_stereo_bonds.wrapping_add(number);
                }
            }
        }
        if let Some(output) = nMaxNumStereoAtoms.as_deref_mut() {
            *output = max_stereo_atoms;
        }
        if let Some(output) = nMaxNumStereoBonds.as_deref_mut() {
            *output = max_stereo_bonds;
        }
    }

    let min_sb_ring_size = (((nMode & u64::from(REQ_MODE_MIN_SB_RING_MASK))
        >> REQ_MODE_MIN_SB_RING_SHFT)
        & 0xffff) as AT_RANK;
    let mut q = SourceMutPointer::<QUEUE>::null();
    let mut nAtomLevel = SourceMutPointer::<AT_RANK>::null();
    let mut cSource = SourceMutPointer::<S_CHAR>::null();

    if min_sb_ring_size >= 3 {
        q = QueueCreate(
            heap,
            num_at.wrapping_add(1),
            std::mem::size_of::<qInt>() as i32,
        )?;
        nAtomLevel = match inchi_calloc::<AT_RANK>(
            heap,
            num_at as u64,
            std::mem::size_of::<AT_RANK>() as u64,
        ) {
            Ok(pointer) => pointer,
            Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
            Err(error) => {
                QueueDelete(heap, q)?;
                return Err(error);
            }
        };
        cSource =
            match inchi_calloc::<S_CHAR>(heap, num_at as u64, std::mem::size_of::<S_CHAR>() as u64)
            {
                Ok(pointer) => pointer,
                Err(SourceHeapError::AllocationFailed) => SourceMutPointer::null(),
                Err(error) => {
                    QueueDelete(heap, q)?;
                    inchi_free(heap, nAtomLevel)?;
                    return Err(error);
                }
            };
        if q.is_null() || cSource.is_null() || nAtomLevel.is_null() {
            QueueDelete(heap, q)?;
            inchi_free(heap, nAtomLevel)?;
            inchi_free(heap, cSource)?;
            return Ok(CT_OUT_OF_RAM);
        }
    }
    let effective_min_sb_ring_size = if min_sb_ring_size >= 3 {
        min_sb_ring_size
    } else {
        2
    };

    let result = (|| -> Result<i32, SourceHeapError> {
        let mut num_3D_stereo_atoms = 0_i32;
        let mut num_stereo_bonds = 0_i32;
        for i in 0..num_at {
            let removed_h = at.as_const().offset(i64::from(num_at))?;
            let mut is_stereo = set_stereo_atom_parity(
                pCG,
                heap,
                at_output,
                at,
                i,
                removed_h,
                num_removed_H,
                bPointedEdgeStereo,
                vABParityUnknown,
                bLooseTSACheck,
                bStereoAtZz,
            )?;
            if is_stereo != 0 {
                num_3D_stereo_atoms = num_3D_stereo_atoms.wrapping_add(i32::from(
                    (AB_MIN_WELL_DEFINED_PARITY as i32..=AB_MAX_WELL_DEFINED_PARITY as i32)
                        .contains(&is_stereo.abs()),
                ));
            } else {
                is_stereo = set_stereo_bonds_parity(
                    heap,
                    at_output,
                    at,
                    i,
                    removed_h,
                    num_removed_H,
                    nMode,
                    q,
                    nAtomLevel,
                    cSource,
                    effective_min_sb_ring_size,
                    bPointedEdgeStereo,
                    vABParityUnknown,
                )?;
                if (CT_ERR_MIN..=CT_ERR_MAX).contains(&is_stereo) {
                    num_3D_stereo_atoms = is_stereo;
                    break;
                }
                num_stereo_bonds = num_stereo_bonds.wrapping_add(i32::from(is_stereo != 0));
            }
        }

        if max_stereo_atoms < num_3D_stereo_atoms {
            if let Some(output) = nMaxNumStereoAtoms.as_deref_mut() {
                *output = num_3D_stereo_atoms;
            }
        }
        if max_stereo_bonds < num_stereo_bonds {
            if let Some(output) = nMaxNumStereoBonds.as_deref_mut() {
                *output = num_stereo_bonds;
            }
        }
        Ok(num_3D_stereo_atoms)
    })();

    QueueDelete(heap, q)?;
    inchi_free(heap, nAtomLevel)?;
    inchi_free(heap, cSource)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn FixUnkn0DStereoBonds(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_at: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2006 FixUnkn0DStereoBonds
    // INCHI✔️❌: int FixUnkn0DStereoBonds( inp_ATOM *at, int num_at )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, m, num = 0;
    // INCHI✔️❌:
    // INCHI✔️❌:     /* Add usual Unknown stereobond descriptors to each Unknown bond */
    // INCHI✔️❌:     for (i = 0; i < num_at; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         for (m = 0; m < MAX_NUM_STEREO_BONDS && at[i].sb_parity[m]; m++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (AB_PARITY_UNKN == at[i].sb_parity[m])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 at[i].bond_stereo[(int) at[i].sb_ord[m]] = STEREO_DBLE_EITHER;
    // INCHI✔️❌:                 num++;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return num;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: FixUnkn0DStereoBonds

    let mut num = 0_i32;
    for i in 0..num_at {
        let i = usize::try_from(i).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
        for m in 0..MAX_NUM_STEREO_BONDS as usize {
            let (parity, bond_ordinal) = {
                let atom = heap
                    .slice(atoms.as_const())?
                    .get(i)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                (atom.sb_parity[m], atom.sb_ord[m])
            };
            if parity == 0 {
                break;
            }
            if parity == AB_PARITY_UNKN as i8 {
                let bond_ordinal = usize::try_from(bond_ordinal)
                    .map_err(|_| SourceHeapError::PointerOutOfBounds)?;
                let atom = heap
                    .slice_mut(atoms)?
                    .get_mut(i)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                let stereo = atom
                    .bond_stereo
                    .get_mut(bond_ordinal)
                    .ok_or(SourceHeapError::PointerOutOfBounds)?;
                *stereo = STEREO_DBLE_EITHER as i8;
                num = num
                    .checked_add(1)
                    .ok_or(SourceHeapError::SourceIntegerOverflow)?;
            }
        }
    }
    Ok(num)
}

#[allow(non_snake_case)]
pub(crate) fn ReconcileAllCmlBondParities(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    num_atoms: i32,
    b_disconnected: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4663 ReconcileAllCmlBondParities
    // INCHI✔️❌: int ReconcileAllCmlBondParities( inp_ATOM *at,
    // INCHI✔️❌:                                  int num_atoms,
    // INCHI✔️❌:                                  int bDisconnected )
    // INCHI✔️❌: {
    // INCHI✔️❌:     int i, ret = 0;
    // INCHI✔️❌:     S_CHAR *visited = (S_CHAR*) inchi_calloc( num_atoms, sizeof( *visited ) );
    // INCHI✔️❌:     if (!visited)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return -1; /* out of RAM */
    // INCHI✔️❌:     }
    // INCHI✔️❌:     for (i = 0; i < num_atoms; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         if (at[i].sb_parity[0] && !visited[i] && !( bDisconnected && is_el_a_metal( at[i].el_number ) ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if ((ret = ReconcileCmlIncidentBondParities( at, i, -1, visited, bDisconnected ))) /* djb-rwth: addressing LLVM warning */
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break; /* error */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     inchi_free( visited );
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ReconcileAllCmlBondParities

    let visited = match inchi_calloc::<i8>(heap, num_atoms as u64, 1) {
        Ok(visited) => visited,
        Err(_) => return Ok(-1),
    };
    let result = (|| {
        let mut ret = 0_i32;
        for i in 0..num_atoms {
            let index = usize::try_from(i).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let atom = heap
                .slice(atoms.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            let was_visited = *heap
                .slice(visited.as_const())?
                .get(index)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if atom.sb_parity[0] != 0
                && was_visited == 0
                && !(b_disconnected != 0 && is_el_a_metal(i32::from(atom.el_number))? != 0)
            {
                ret =
                    ReconcileCmlIncidentBondParities(heap, atoms, i, -1, visited, b_disconnected)?;
                if ret != 0 {
                    break;
                }
            }
        }
        Ok(ret)
    })();
    inchi_free(heap, visited)?;
    result
}

#[allow(non_snake_case)]
pub(crate) fn ReconcileCmlIncidentBondParities(
    heap: &mut SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    cur_atom: i32,
    prev_atom: i32,
    visited: SourceMutPointer<i8>,
    b_disconnected: i32,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4690 ReconcileCmlIncidentBondParities
    // INCHI✔️❌: int ReconcileCmlIncidentBondParities( inp_ATOM *at,
    // INCHI✔️❌:                                       int cur_atom,
    // INCHI✔️❌:                                       int prev_atom,
    // INCHI✔️❌:                                       S_CHAR *visited,
    // INCHI✔️❌:                                       int bDisconnected )
    // INCHI✔️❌: {
    // INCHI✔️❌:     /* visited = 0 or parity => atom has not been visited
    // INCHI✔️❌:                  10 + parity => currently is on the stack + its final parity
    // INCHI✔️❌:                  20 + parity => has been visited; is not on the stack anymore + its final parity */
    // INCHI✔️❌:     int i, j, nxt_atom, ret = 0, len;
    // INCHI✔️❌:     int icur2nxt, icur2neigh;   /* cur atom neighbors */
    // INCHI✔️❌:     int inxt2cur, inxt2neigh;   /* next atom neighbors */
    // INCHI✔️❌:     int cur_parity, nxt_parity;
    // INCHI✔️❌:     int cur_order_parity, nxt_order_parity, cur_sb_parity, nxt_sb_parity, bCurMask, bNxtMask;
    // INCHI✔️❌:     /* !(bDisconnected && is_el_a_metal(at[i].el_number) */
    // INCHI✔️❌:
    // INCHI✔️❌:     if (at[cur_atom].valence > MAX_NUM_STEREO_BONDS)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 0; /* ignore */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (!at[cur_atom].sb_parity[0])
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 1; /* wrong call */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (visited[cur_atom] >= 10)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         return 2; /* program error */
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     cur_parity = visited[cur_atom] % 10;
    // INCHI✔️❌:
    // INCHI✔️❌:     visited[cur_atom] += 10;
    // INCHI✔️❌:
    // INCHI✔️❌:     for (i = 0; i < MAX_NUM_STEREO_BONDS && at[cur_atom].sb_parity[i]; i++)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         icur2nxt = (int) at[cur_atom].sb_ord[i];
    // INCHI✔️❌:         len = get_opposite_sb_atom( at, cur_atom, icur2nxt, &nxt_atom, &inxt2cur, &j );
    // INCHI✔️❌:         if (!len)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 4; /* could not find the opposite atom: bond parity data error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (nxt_atom == prev_atom)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (visited[nxt_atom] >= 20)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* back edge, second visit: ignore */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[nxt_atom].valence > MAX_NUM_STEREO_BONDS)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             continue; /* may be treated only after metal disconnection */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (bDisconnected && ( at[cur_atom].sb_parity[i] & SB_PARITY_FLAG ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cur_sb_parity = ( at[cur_atom].sb_parity[i] >> SB_PARITY_SHFT );
    // INCHI✔️❌:             bCurMask = 3 << SB_PARITY_SHFT;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cur_sb_parity = ( at[cur_atom].sb_parity[i] & SB_PARITY_MASK );
    // INCHI✔️❌:             bCurMask = 3;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (bDisconnected && ( at[nxt_atom].sb_parity[j] & SB_PARITY_FLAG ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nxt_sb_parity = ( at[nxt_atom].sb_parity[j] >> SB_PARITY_SHFT );
    // INCHI✔️❌:             bNxtMask = 3 << SB_PARITY_SHFT;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nxt_sb_parity = ( at[nxt_atom].sb_parity[j] & SB_PARITY_MASK );
    // INCHI✔️❌:             bNxtMask = 3;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!ATOM_PARITY_WELL_DEF( cur_sb_parity ) ||
    // INCHI✔️❌:              !ATOM_PARITY_WELL_DEF( nxt_sb_parity ))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_sb_parity == nxt_sb_parity)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 continue;
    // INCHI✔️❌:                 /* goto move_forward;       */
    // INCHI✔️❌:                 /* bypass unknown/undefined */
    // INCHI✔️❌:             }
    // INCHI✔️❌:             return 3; /* sb parities do not match: bond parity data error */
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         icur2neigh = (int) at[cur_atom].sn_ord[i];
    // INCHI✔️❌:         inxt2neigh = (int) at[nxt_atom].sn_ord[j];
    // INCHI✔️❌:         /* Parity of at[cur_atom].neighbor[] premutation to reach this order: { next_atom, neigh_atom, ...} */
    // INCHI✔️❌:
    // INCHI✔️❌:         /* 1. move next_atom  from position=icur2nxt to position=0 =>
    // INCHI✔️❌:          *         icur2nxt permutations
    // INCHI✔️❌:          * 2. move neigh_atom from position=inxt2neigh+(inxt2cur > inxt2neigh) to position=1 =>
    // INCHI✔️❌:          *         inxt2neigh+(inxt2cur > inxt2neigh)-1 permutations.
    // INCHI✔️❌:          * Note if (inxt2cur > inxt2neigh) then move #1 increments neigh_atom position
    // INCHI✔️❌:          * Note add 4 because icur2neigh may be negative due to isotopic H removal
    // INCHI✔️❌:          */
    // INCHI✔️❌:         cur_order_parity = ( 4 + icur2nxt + icur2neigh + ( icur2neigh > icur2nxt ) ) % 2;
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Same for next atom: */
    // INCHI✔️❌:         /* parity of at[nxt_atom].neighbor[] premutation to reach this order: { cur_atom, neigh_atom, ...} */
    // INCHI✔️❌:         nxt_order_parity = ( 4 + inxt2cur + inxt2neigh + ( inxt2neigh > inxt2cur ) ) % 2;
    // INCHI✔️❌:         nxt_parity = visited[nxt_atom] % 10;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!cur_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             cur_parity = 2 - ( cur_order_parity + cur_sb_parity ) % 2;
    // INCHI✔️❌:             visited[cur_atom] += cur_parity;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_parity != 2 - ( cur_order_parity + cur_sb_parity ) % 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /***** Reconcile bond parities *****/
    // INCHI✔️❌:
    // INCHI✔️❌:                 /* Each bond parity is split into two values located at the end atoms.
    // INCHI✔️❌:                    For T (trans) the values are (1,1) or (2,2)
    // INCHI✔️❌:                    For C (cis)   the values are (1,2) or (2,1)
    // INCHI✔️❌:                    The fact that one pair = another with inverted parities, namely
    // INCHI✔️❌:                    Inv(1,1) = (2,2) and Inv(1,2) = (2,1), allows to
    // INCHI✔️❌:                    simultaneouly invert parities of the current bond end atoms
    // INCHI✔️❌:                    (at[cur_atom].sb_parity[i], at[nxt_atom].sb_parity[j])
    // INCHI✔️❌:                    so that the final current atom parity cur_parity
    // INCHI✔️❌:                    calculated later in stereochemical canonicalization for
    // INCHI✔️❌:                    each stereobond incident with the current atomis same.
    // INCHI✔️❌:                    Achieving this is called here RECONCILIATION.
    // INCHI✔️❌:                    If at the closure of an aromatic circuit the parities of
    // INCHI✔️❌:                    next atom cannot be reconciled with already calculated then
    // INCHI✔️❌:                    this function returns 5 (error).
    // INCHI✔️❌:                 */
    // INCHI✔️❌:                 at[cur_atom].sb_parity[i] ^= bCurMask;
    // INCHI✔️❌:                 at[nxt_atom].sb_parity[j] ^= bNxtMask;
    // INCHI✔️❌:                 /* djb-rwth: removing redundant code */
    // INCHI✔️❌:                 nxt_sb_parity ^= 3;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!nxt_parity)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             nxt_parity = 2 - ( nxt_order_parity + nxt_sb_parity ) % 2;
    // INCHI✔️❌:             visited[nxt_atom] += nxt_parity;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (nxt_parity != 2 - ( nxt_order_parity + nxt_sb_parity ) % 2)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 return 5; /* algorithm does not work for Mebius-like structures */
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         /* Move_forward: */
    // INCHI✔️❌:         if (visited[nxt_atom] < 10)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             ret = ReconcileCmlIncidentBondParities( at, nxt_atom, cur_atom, visited, bDisconnected );
    // INCHI✔️❌:             if (ret)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 break;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:     visited[cur_atom] += 10; /* all bonds incident to the current atom have
    // INCHI✔️❌:                                 been processed or an error occurred. */
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: ReconcileCmlIncidentBondParities

    let cur_index = usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
    let current = heap
        .slice(atoms.as_const())?
        .get(cur_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if i32::from(current.valence) > MAX_NUM_STEREO_BONDS as i32 {
        return Ok(0);
    }
    if current.sb_parity[0] == 0 {
        return Ok(1);
    }
    let visited_value = *heap
        .slice(visited.as_const())?
        .get(cur_index)
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if visited_value >= 10 {
        return Ok(2);
    }
    let mut cur_parity = i32::from(visited_value % 10);
    heap.slice_mut(visited)?[cur_index] = visited_value.wrapping_add(10);

    let mut ret = 0_i32;
    for i in 0..MAX_NUM_STEREO_BONDS as usize {
        let current_sb_parity = heap.slice(atoms.as_const())?[cur_index].sb_parity[i];
        if current_sb_parity == 0 {
            break;
        }
        let current_to_next = i32::from(heap.slice(atoms.as_const())?[cur_index].sb_ord[i]);
        let (mut next_atom, mut next_to_current, mut j) = (0_i32, 0_i32, 0_i32);
        let len = get_opposite_sb_atom(
            heap,
            atoms,
            cur_atom,
            current_to_next,
            Some(&mut next_atom),
            Some(&mut next_to_current),
            Some(&mut j),
        )?;
        if len == 0 {
            return Ok(4);
        }
        if next_atom == prev_atom {
            continue;
        }
        let next_index =
            usize::try_from(next_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        if heap.slice(visited.as_const())?[next_index] >= 20 {
            continue;
        }
        if i32::from(heap.slice(atoms.as_const())?[next_index].valence)
            > MAX_NUM_STEREO_BONDS as i32
        {
            continue;
        }
        let j = usize::try_from(j).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let current_parity_bits = i32::from(heap.slice(atoms.as_const())?[cur_index].sb_parity[i]);
        let next_parity_bits = i32::from(heap.slice(atoms.as_const())?[next_index].sb_parity[j]);
        let (cur_sb_parity, current_mask) =
            if b_disconnected != 0 && current_parity_bits & SB_PARITY_FLAG as i32 != 0 {
                (
                    current_parity_bits >> SB_PARITY_SHFT,
                    3_i32 << SB_PARITY_SHFT,
                )
            } else {
                (current_parity_bits & SB_PARITY_MASK as i32, 3)
            };
        let (mut next_sb_parity, next_mask) =
            if b_disconnected != 0 && next_parity_bits & SB_PARITY_FLAG as i32 != 0 {
                (next_parity_bits >> SB_PARITY_SHFT, 3_i32 << SB_PARITY_SHFT)
            } else {
                (next_parity_bits & SB_PARITY_MASK as i32, 3)
            };
        let well_defined = |parity: i32| {
            AB_MIN_WELL_DEFINED_PARITY as i32 <= parity
                && parity <= AB_MAX_WELL_DEFINED_PARITY as i32
        };
        if !well_defined(cur_sb_parity) || !well_defined(next_sb_parity) {
            if cur_sb_parity == next_sb_parity {
                continue;
            }
            return Ok(3);
        }

        let current_to_neighbor = i32::from(heap.slice(atoms.as_const())?[cur_index].sn_ord[i]);
        let next_to_neighbor = i32::from(heap.slice(atoms.as_const())?[next_index].sn_ord[j]);
        let current_order_parity = (4
            + current_to_next
            + current_to_neighbor
            + i32::from(current_to_neighbor > current_to_next))
            % 2;
        let next_order_parity = (4
            + next_to_current
            + next_to_neighbor
            + i32::from(next_to_neighbor > next_to_current))
            % 2;
        let next_visited = heap.slice(visited.as_const())?[next_index];
        let next_parity = i32::from(next_visited % 10);

        if cur_parity == 0 {
            cur_parity = 2 - (current_order_parity + cur_sb_parity) % 2;
            let current_visited = heap.slice(visited.as_const())?[cur_index];
            heap.slice_mut(visited)?[cur_index] = current_visited.wrapping_add(cur_parity as i8);
        } else if cur_parity != 2 - (current_order_parity + cur_sb_parity) % 2 {
            heap.slice_mut(atoms)?[cur_index].sb_parity[i] ^= current_mask as i8;
            heap.slice_mut(atoms)?[next_index].sb_parity[j] ^= next_mask as i8;
            next_sb_parity ^= 3;
        }

        if next_parity == 0 {
            let calculated = 2 - (next_order_parity + next_sb_parity) % 2;
            heap.slice_mut(visited)?[next_index] = next_visited.wrapping_add(calculated as i8);
        } else if next_parity != 2 - (next_order_parity + next_sb_parity) % 2 {
            return Ok(5);
        }

        if heap.slice(visited.as_const())?[next_index] < 10 {
            ret = ReconcileCmlIncidentBondParities(
                heap,
                atoms,
                next_atom,
                cur_atom,
                visited,
                b_disconnected,
            )?;
            if ret != 0 {
                break;
            }
        }
    }
    let value = heap.slice(visited.as_const())?[cur_index];
    heap.slice_mut(visited)?[cur_index] = value.wrapping_add(10);
    Ok(ret)
}

pub(crate) fn get_opposite_sb_atom(
    heap: &SourceHeap,
    atoms: SourceMutPointer<inp_ATOM>,
    cur_atom: i32,
    icur2nxt: i32,
    next_atom_out: Option<&mut i32>,
    next_to_current_out: Option<&mut i32>,
    next_parity_ordinal_out: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4861 get_opposite_sb_atom
    // INCHI✔️❌: int get_opposite_sb_atom( inp_ATOM *at,
    // INCHI✔️❌:                           int cur_atom,
    // INCHI✔️❌:                           int icur2nxt,
    // INCHI✔️❌:                           int *pnxt_atom,
    // INCHI✔️❌:                           int *pinxt2cur,
    // INCHI✔️❌:                           int *pinxt_sb_parity_ord )
    // INCHI✔️❌: {
    // INCHI✔️❌:     AT_NUMB nxt_atom;
    // INCHI✔️❌:     int     j, len;
    // INCHI✔️❌:
    // INCHI✔️❌:     len = 0;
    // INCHI✔️❌:     while (len++ < 20)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         /* Arbitrarily set cumulene length limit to avoid infinite loop */
    // INCHI✔️❌:         nxt_atom = at[cur_atom].neighbor[icur2nxt];
    // INCHI✔️❌:         for (j = 0; j < MAX_NUM_STEREO_BONDS && at[nxt_atom].sb_parity[j]; j++)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             if (cur_atom == at[nxt_atom].neighbor[(int) at[nxt_atom].sb_ord[j]])
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 /* Found the opposite atom */
    // INCHI✔️❌:                 *pnxt_atom = nxt_atom;
    // INCHI✔️❌:                 *pinxt2cur = at[nxt_atom].sb_ord[j];
    // INCHI✔️❌:                 *pinxt_sb_parity_ord = j;
    // INCHI✔️❌:                 return len;
    // INCHI✔️❌:             }
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (j)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /* reached atom(s) with stereobond (sb) parity, the opposite atom has not been found */
    // INCHI✔️❌:         }
    // INCHI✔️❌:         if (at[nxt_atom].valence == 2 && 2 * BOND_TYPE_DOUBLE == at[nxt_atom].chem_bonds_valence)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             /* Follow cumulene =X= path */
    // INCHI✔️❌:             icur2nxt = ( at[nxt_atom].neighbor[0] == cur_atom );
    // INCHI✔️❌:             cur_atom = nxt_atom;
    // INCHI✔️❌:         }
    // INCHI✔️❌:         else
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 0; /* neither atom with a sb parity not middle cumulene could be reached */
    // INCHI✔️❌:         }
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return 0; /* too long chain of cumulene was found */
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: get_opposite_sb_atom

    get_opposite_sb_atom_slice(
        heap.slice(atoms.as_const())?,
        cur_atom,
        icur2nxt,
        next_atom_out,
        next_to_current_out,
        next_parity_ordinal_out,
    )
}

pub(crate) fn get_opposite_sb_atom_slice(
    atoms: &[inp_ATOM],
    mut cur_atom: i32,
    mut icur2nxt: i32,
    mut next_atom_out: Option<&mut i32>,
    mut next_to_current_out: Option<&mut i32>,
    mut next_parity_ordinal_out: Option<&mut i32>,
) -> Result<i32, SourceHeapError> {
    for len in 1..=20_i32 {
        let current_index =
            usize::try_from(cur_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next_ordinal =
            usize::try_from(icur2nxt).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next_atom = i32::from(
            *atoms
                .get(current_index)
                .and_then(|atom| atom.neighbor.get(next_ordinal))
                .ok_or(SourceHeapError::PointerOutOfBounds)?,
        );
        let next_index =
            usize::try_from(next_atom).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
        let next = atoms
            .get(next_index)
            .ok_or(SourceHeapError::PointerOutOfBounds)?;
        let mut j = 0_usize;
        while j < MAX_NUM_STEREO_BONDS as usize && next.sb_parity[j] != 0 {
            let reverse_ordinal =
                usize::try_from(next.sb_ord[j]).map_err(|_| SourceHeapError::PointerOutOfBounds)?;
            let reverse_atom = *next
                .neighbor
                .get(reverse_ordinal)
                .ok_or(SourceHeapError::PointerOutOfBounds)?;
            if cur_atom == i32::from(reverse_atom) {
                *next_atom_out
                    .as_deref_mut()
                    .ok_or(SourceHeapError::NullPointer)? = next_atom;
                *next_to_current_out
                    .as_deref_mut()
                    .ok_or(SourceHeapError::NullPointer)? = i32::from(next.sb_ord[j]);
                *next_parity_ordinal_out
                    .as_deref_mut()
                    .ok_or(SourceHeapError::NullPointer)? = j as i32;
                return Ok(len);
            }
            j += 1;
        }
        if j != 0 {
            return Ok(0);
        }
        if next.valence == 2 && (2 * BOND_TYPE_DOUBLE as i32) == i32::from(next.chem_bonds_valence)
        {
            icur2nxt = i32::from(next.neighbor[0] == cur_atom as u16);
            cur_atom = next_atom;
        } else {
            return Ok(0);
        }
    }
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_port__ichister__comp_at_numb__line_140() {
        assert_eq!(comp_AT_NUMB(0, 0), 0);
        assert_eq!(comp_AT_NUMB(1, 0), 1);
        assert_eq!(comp_AT_NUMB(0, 1), -1);
        assert_eq!(comp_AT_NUMB(AT_NUMB::MAX, 0), i32::from(AT_NUMB::MAX));
        assert_eq!(comp_AT_NUMB(0, AT_NUMB::MAX), -i32::from(AT_NUMB::MAX));
        assert_eq!(comp_AT_NUMB(AT_NUMB::MAX, AT_NUMB::MAX), 0);
    }

    #[test]
    fn source_port__ichister__compdble__line_379() {
        let values = [
            f64::NEG_INFINITY,
            -f64::MAX,
            -0.0,
            0.0,
            f64::MAX,
            f64::INFINITY,
            f64::NAN,
        ];
        assert_eq!(CompDble(4, 1, &values), Ok(1));
        assert_eq!(CompDble(1, 4, &values), Ok(-1));
        assert_eq!(CompDble(2, 3, &values), Ok(0));
        assert_eq!(CompDble(3, 2, &values), Ok(0));
        assert_eq!(CompDble(5, 5, &values), Ok(0));
        assert_eq!(CompDble(6, 4, &values), Ok(0));
        assert_eq!(CompDble(4, 6, &values), Ok(0));
        assert_eq!(CompDble(0, 5, &values), Ok(-1));
        assert_eq!(CompDble(5, 0, &values), Ok(1));
        assert_eq!(
            CompDble(-1, 0, &values),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompDble(0, -1, &values),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            CompDble(values.len() as i32, 0, &values),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }
    use crate::source_types::{AB_PARITY_EVEN, AB_PARITY_ODD};

    #[test]
    fn source_port__ichister__get_z_coord__line_147() {
        fn evaluate(
            stereo: i8,
            current_z: f64,
            selected_z: f64,
            pointed: i32,
            other_z: Option<f64>,
        ) -> (f64, i32) {
            let mut current = inp_ATOM {
                z: current_z,
                valence: if other_z.is_some() { 2 } else { 1 },
                ..inp_ATOM::default()
            };
            current.neighbor[0] = 1;
            current.bond_stereo[0] = stereo;

            let mut atoms = vec![
                current,
                inp_ATOM {
                    z: selected_z,
                    ..inp_ATOM::default()
                },
            ];
            if let Some(z) = other_z {
                atoms[0].neighbor[1] = 2;
                atoms.push(inp_ATOM {
                    z,
                    ..inp_ATOM::default()
                });
            }

            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let mut coordinate_type = i32::MIN;
            let z =
                get_z_coord(&heap, atoms.as_const(), 0, 0, &mut coordinate_type, pointed).unwrap();
            (z, coordinate_type)
        }

        for (stereo, expected_type) in [
            (0, ZTYPE_NONE as i32),
            (STEREO_SNGL_UP as i8, ZTYPE_UP as i32),
            (STEREO_SNGL_EITHER as i8, ZTYPE_EITHER as i32),
            (STEREO_SNGL_DOWN as i8, ZTYPE_DOWN),
            (2, ZTYPE_NONE as i32),
            (i8::MIN, ZTYPE_NONE as i32),
        ] {
            let (z, coordinate_type) = evaluate(stereo, 4.0, 4.0, 0, None);
            assert_eq!(z.to_bits(), 0.0_f64.to_bits());
            assert_eq!(coordinate_type, expected_type);
        }

        assert_eq!(
            evaluate(-(STEREO_SNGL_UP as i8), 0.0, 0.0, 0, None).1,
            -(ZTYPE_UP as i32)
        );
        assert_eq!(
            evaluate(-(STEREO_SNGL_DOWN as i8), 0.0, 0.0, 0, None).1,
            -ZTYPE_DOWN
        );

        for (stereo, pointed, expected_type) in [
            (STEREO_SNGL_UP as i8, 1, ZTYPE_UP as i32),
            (STEREO_SNGL_UP as i8, -1, ZTYPE_NONE as i32),
            (-(STEREO_SNGL_UP as i8), 1, ZTYPE_NONE as i32),
            (-(STEREO_SNGL_UP as i8), -1, -(ZTYPE_UP as i32)),
        ] {
            assert_eq!(evaluate(stereo, 0.0, 0.0, pointed, None).1, expected_type);
        }

        let below = f64::from_bits(ZERO_LENGTH.to_bits() - 1);
        let above = f64::from_bits(ZERO_LENGTH.to_bits() + 1);
        assert_eq!(
            evaluate(STEREO_SNGL_UP as i8, 0.0, below, 0, None),
            (below, ZTYPE_UP as i32)
        );
        assert_eq!(
            evaluate(STEREO_SNGL_UP as i8, 0.0, ZERO_LENGTH, 0, None),
            (ZERO_LENGTH, ZTYPE_3D as i32)
        );
        assert_eq!(
            evaluate(STEREO_SNGL_UP as i8, 0.0, above, 0, None),
            (above, ZTYPE_3D as i32)
        );

        assert_eq!(
            evaluate(STEREO_SNGL_UP as i8, 0.0, below, 0, Some(above)),
            (below, ZTYPE_3D as i32)
        );
        assert_eq!(
            evaluate(STEREO_SNGL_EITHER as i8, 1.0, 2.0, 0, None),
            (1.0, ZTYPE_EITHER as i32)
        );
        assert_eq!(
            evaluate(-(STEREO_SNGL_EITHER as i8), 1.0, 2.0, 1, None),
            (1.0, ZTYPE_3D as i32)
        );
        assert_eq!(
            evaluate(STEREO_SNGL_UP as i8, -2.5, 7.25, 0, None),
            (9.75, ZTYPE_3D as i32)
        );

        let (nan_z, coordinate_type) = evaluate(STEREO_SNGL_EITHER as i8, 0.0, f64::NAN, 0, None);
        assert!(nan_z.is_nan());
        assert_eq!(coordinate_type, ZTYPE_EITHER as i32);

        let (negative_zero, coordinate_type) = evaluate(STEREO_SNGL_UP as i8, 0.0, -0.0, 0, None);
        assert_eq!(negative_zero.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(coordinate_type, ZTYPE_UP as i32);
    }

    #[test]
    fn source_port__ichister__len3__line_221() {
        assert_eq!(len3(&[3.0, 4.0, 12.0]).to_bits(), 13.0_f64.to_bits());
        assert_eq!(len3(&[0.0, -0.0, 0.0]).to_bits(), 0.0_f64.to_bits());
        assert_eq!(len3(&[f64::MIN_POSITIVE, 0.0, 0.0]), 0.0);
        assert_eq!(len3(&[f64::from_bits(1), 0.0, 0.0]), 0.0);
        assert!(len3(&[f64::MAX, 0.0, 0.0]).is_infinite());
        assert!(len3(&[f64::INFINITY, 1.0, 2.0]).is_infinite());
        assert!(len3(&[f64::NEG_INFINITY, 1.0, 2.0]).is_infinite());
        assert!(len3(&[f64::NAN, 1.0, 2.0]).is_nan());
        assert!(len3(&[1.0, f64::NAN, 2.0]).is_nan());
        assert!(len3(&[1.0, 2.0, f64::NAN]).is_nan());
    }

    #[test]
    fn source_port__ichister__len2__line_234() {
        assert_eq!(len2(&[3.0, 4.0]).to_bits(), 5.0_f64.to_bits());
        assert_eq!(len2(&[-3.0, -4.0]).to_bits(), 5.0_f64.to_bits());
        assert_eq!(len2(&[0.0, -0.0]).to_bits(), 0.0_f64.to_bits());
        assert_eq!(len2(&[f64::MIN_POSITIVE, 0.0]), 0.0);
        assert_eq!(len2(&[f64::from_bits(1), 0.0]), 0.0);
        assert!(len2(&[f64::MAX, 0.0]).is_infinite());
        assert!(len2(&[f64::INFINITY, 1.0]).is_infinite());
        assert!(len2(&[f64::NEG_INFINITY, 1.0]).is_infinite());
        assert!(len2(&[f64::NAN, 1.0]).is_nan());
        assert!(len2(&[1.0, f64::NAN]).is_nan());
    }

    #[test]
    fn source_port__ichister__diff3__line_247() {
        let mut result = [91.0, 92.0, 93.0];
        let result_pointer = result.as_mut_ptr();
        let returned = diff3(&[5.0, -2.0, 3.5], &[1.5, 4.0, -0.5], &mut result);
        assert_eq!(returned.as_mut_ptr(), result_pointer);
        assert_eq!(*returned, [3.5, -6.0, 4.0]);

        diff3(&[0.0, -0.0, 0.0], &[0.0, 0.0, -0.0], &mut result);
        assert_eq!(result[0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(result[1].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(result[2].to_bits(), 0.0_f64.to_bits());

        diff3(
            &[f64::INFINITY, f64::INFINITY, f64::NAN],
            &[f64::INFINITY, f64::NEG_INFINITY, 1.0],
            &mut result,
        );
        assert!(result[0].is_nan());
        assert_eq!(result[1], f64::INFINITY);
        assert!(result[2].is_nan());

        diff3(
            &[f64::MAX, f64::MIN, f64::from_bits(2)],
            &[f64::MIN, f64::MAX, f64::from_bits(1)],
            &mut result,
        );
        assert_eq!(result[0], f64::INFINITY);
        assert_eq!(result[1], f64::NEG_INFINITY);
        assert_eq!(result[2].to_bits(), 1);
    }

    #[test]
    fn source_port__ichister__add3__line_259() {
        let mut result = [91.0, 92.0, 93.0];
        add3(&[1.5, -2.0, 3.0], &[2.5, 5.0, -7.0], &mut result);
        assert_eq!(result, [4.0, 3.0, -4.0]);

        add3(
            &[-0.0, 0.0, f64::INFINITY],
            &[-0.0, -0.0, f64::NEG_INFINITY],
            &mut result,
        );
        assert_eq!(result[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(result[1].to_bits(), 0.0_f64.to_bits());
        assert!(result[2].is_nan());

        add3(
            &[f64::MAX, f64::MIN, f64::NAN],
            &[f64::MAX, f64::MIN, 1.0],
            &mut result,
        );
        assert_eq!(result[0], f64::INFINITY);
        assert_eq!(result[1], f64::NEG_INFINITY);
        assert!(result[2].is_nan());
    }

    #[test]
    fn source_port__ichister__mult3__line_270() {
        let mut result = [91.0, 92.0, 93.0];
        mult3(&[1.5, -2.0, 3.0], 4.0, &mut result);
        assert_eq!(result, [6.0, -8.0, 12.0]);

        mult3(&[0.0, -0.0, f64::INFINITY], -2.0, &mut result);
        assert_eq!(result[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(result[1].to_bits(), 0.0_f64.to_bits());
        assert_eq!(result[2], f64::NEG_INFINITY);

        mult3(&[0.0, f64::MAX, f64::NAN], f64::INFINITY, &mut result);
        assert!(result[0].is_nan());
        assert_eq!(result[1], f64::INFINITY);
        assert!(result[2].is_nan());

        mult3(&[f64::MAX, f64::MIN, f64::from_bits(1)], 2.0, &mut result);
        assert_eq!(result[0], f64::INFINITY);
        assert_eq!(result[1], f64::NEG_INFINITY);
        assert_eq!(result[2].to_bits(), 2);
    }

    #[test]
    fn source_port__ichister__change_sign3__line_292() {
        let values = [0.0, -0.0, f64::from_bits(1)];
        let mut result = [91.0, 92.0, 93.0];
        change_sign3(&values, &mut result);
        assert_eq!(result[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(result[1].to_bits(), 0.0_f64.to_bits());
        assert_eq!(result[2].to_bits(), (-f64::from_bits(1)).to_bits());

        change_sign3(&[f64::MAX, f64::MIN, f64::NAN], &mut result);
        assert_eq!(result[0], -f64::MAX);
        assert_eq!(result[1], f64::MAX);
        assert!(result[2].is_nan());

        change_sign3(&[f64::INFINITY, f64::NEG_INFINITY, -1.25], &mut result);
        assert_eq!(result, [f64::NEG_INFINITY, f64::INFINITY, 1.25]);
    }

    #[test]
    fn source_port__ichister__dot_prod3__line_303() {
        assert_eq!(dot_prod3(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]), 32.0);
        assert_eq!(dot_prod3(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0]), 0.0);
        assert_eq!(
            dot_prod3(&[1.0, 1.0, 1.0], &[f64::MAX, -f64::MAX, 1.0]),
            1.0
        );
        assert_eq!(
            dot_prod3(&[-0.0, -0.0, -0.0], &[1.0, 1.0, 1.0]).to_bits(),
            (-0.0_f64).to_bits()
        );
        assert!(dot_prod3(&[0.0, 1.0, 2.0], &[f64::INFINITY, 1.0, 2.0]).is_nan());
        assert!(dot_prod3(&[f64::NAN, 1.0, 2.0], &[1.0, 2.0, 3.0]).is_nan());
        assert_eq!(
            dot_prod3(&[f64::MAX, 0.0, 0.0], &[2.0, 0.0, 0.0]),
            f64::INFINITY
        );
    }

    #[test]
    fn source_port__ichister__dot_prodchar3__line_310() {
        for first in i8::MIN..=i8::MAX {
            for second in i8::MIN..=i8::MAX {
                let raw = i32::from(first) * i32::from(second) / 100;
                let expected = raw.clamp(-100, 100);
                assert_eq!(
                    dot_prodchar3(&[first, 0, 0], &[second, 0, 0]),
                    expected,
                    "first={first} second={second}"
                );
            }
        }

        for (a, b, expected) in [
            ([100, 0, 0], [100, 0, 0], 100),
            ([101, 0, 0], [100, 0, 0], 100),
            ([-100, 0, 0], [100, 0, 0], -100),
            ([-101, 0, 0], [100, 0, 0], -100),
            ([99, 0, 0], [100, 0, 0], 99),
            ([-99, 0, 0], [100, 0, 0], -99),
            ([99, 1, 0], [100, 99, 0], 99),
            ([-99, -1, 0], [100, 99, 0], -99),
            ([127, 127, 127], [127, 127, 127], 100),
            ([-128, -128, -128], [127, 127, 127], -100),
            ([-128, -128, -128], [-128, -128, -128], 100),
            ([1, 1, 1], [33, 33, 34], 1),
            ([-1, -1, -1], [33, 33, 34], -1),
            ([1, 1, 1], [33, 33, 33], 0),
            ([-1, -1, -1], [33, 33, 33], 0),
        ] {
            assert_eq!(dot_prodchar3(&a, &b), expected, "a={a:?} b={b:?}");
        }
    }

    #[test]
    fn source_port__ichister__cross_prod3__line_330() {
        let mut result = [91.0, 92.0, 93.0];
        assert_eq!(
            cross_prod3(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], &mut result),
            &[0.0, 0.0, 1.0]
        );
        assert_eq!(
            cross_prod3(&[0.0, 1.0, 0.0], &[1.0, 0.0, 0.0], &mut result),
            &[0.0, 0.0, -1.0]
        );
        assert_eq!(
            cross_prod3(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0], &mut result),
            &[0.0, 0.0, 0.0]
        );
        assert_eq!(
            cross_prod3(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], &mut result),
            &[-3.0, 6.0, -3.0]
        );

        cross_prod3(&[-0.0, -0.0, -0.0], &[1.0, 1.0, 1.0], &mut result);
        assert_eq!(result[0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(result[1].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(result[2].to_bits(), 0.0_f64.to_bits());

        cross_prod3(&[f64::INFINITY, 0.0, 0.0], &[0.0, 1.0, 0.0], &mut result);
        assert_eq!(result[0].to_bits(), 0.0_f64.to_bits());
        assert!(result[1].is_nan());
        assert_eq!(result[2], f64::INFINITY);

        cross_prod3(&[f64::NAN, 1.0, 2.0], &[3.0, 4.0, 5.0], &mut result);
        assert_eq!(result[0], -3.0);
        assert!(result[1].is_nan());
        assert!(result[2].is_nan());
    }

    #[test]
    fn source_port__ichister__triple_prod__line_347() {
        let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let mut sine = 91.0;
        assert_eq!(triple_prod(&x, &y, &z, Some(&mut sine)), 1.0);
        assert_eq!(sine, 1.0);
        assert_eq!(triple_prod(&y, &x, &z, Some(&mut sine)), -1.0);
        assert_eq!(sine, -1.0);

        sine = 77.0;
        assert_eq!(triple_prod(&x, &y, &z, None), 1.0);
        assert_eq!(sine, 77.0);

        for (length, expected_sine) in [
            (f64::from_bits(1.0e-7_f64.to_bits() - 1), 0.0),
            (1.0e-7, 0.0),
            (f64::from_bits(1.0e-7_f64.to_bits() + 1), 1.0),
        ] {
            sine = 91.0;
            let c = [0.0, 0.0, length];
            assert_eq!(triple_prod(&x, &y, &c, Some(&mut sine)), length);
            assert_eq!(sine, expected_sine);
        }

        sine = 91.0;
        assert_eq!(triple_prod(&x, &x, &z, Some(&mut sine)), 0.0);
        assert_eq!(sine.to_bits(), 0.0_f64.to_bits());

        sine = 91.0;
        assert!(triple_prod(&x, &y, &[f64::NAN, 0.0, 1.0], Some(&mut sine)).is_nan());
        assert_eq!(sine.to_bits(), 0.0_f64.to_bits());

        sine = 91.0;
        assert!(triple_prod(&x, &y, &[0.0, 0.0, f64::INFINITY], Some(&mut sine)).is_infinite());
        assert!(sine.is_nan());
    }

    #[test]
    fn source_port__ichister__triple_prod_and_min_abs_sine__line_1145() {
        let orthogonal = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let direct = triple_prod(&orthogonal[0], &orthogonal[1], &orthogonal[2], None);
        assert_eq!(
            triple_prod_and_min_abs_sine(&orthogonal, None).to_bits(),
            direct.to_bits()
        );
        let mut minimum = -91.0;
        assert_eq!(
            triple_prod_and_min_abs_sine(&orthogonal, Some(&mut minimum)),
            1.0
        );
        assert_eq!(minimum, 1.0);

        let negative = [[1.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 3.0, 0.0]];
        assert_eq!(
            triple_prod_and_min_abs_sine(&negative, Some(&mut minimum)),
            -6.0
        );
        assert_eq!(minimum, 1.0);

        let parallel = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert_eq!(
            triple_prod_and_min_abs_sine(&parallel, Some(&mut minimum)),
            0.0
        );
        assert_eq!(minimum, 0.0);

        let tiny = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0e-8]];
        assert_eq!(
            triple_prod_and_min_abs_sine(&tiny, Some(&mut minimum)),
            1.0e-8
        );
        assert_eq!(minimum, 0.0);

        let nan = [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(triple_prod_and_min_abs_sine(&nan, Some(&mut minimum)).is_nan());
        assert_eq!(minimum.to_bits(), 0.0_f64.to_bits());

        let infinite = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, f64::INFINITY]];
        assert!(triple_prod_and_min_abs_sine(&infinite, Some(&mut minimum)).is_nan());
        assert_eq!(minimum.to_bits(), 0.0_f64.to_bits());

        let zero = [[0.0; 3]; 3];
        let product = triple_prod_and_min_abs_sine(&zero, Some(&mut minimum));
        assert_eq!(product.to_bits(), 0.0_f64.to_bits());
        assert_eq!(minimum.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn source_port__ichister__are_3_vect_in_one_plane__line_1176() {
        let orthogonal = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_eq!(
            are_3_vect_in_one_plane(
                &orthogonal,
                f64::from_bits(1.0_f64.to_bits().wrapping_sub(1)),
            ),
            0
        );
        assert_eq!(are_3_vect_in_one_plane(&orthogonal, 1.0), 1);
        assert_eq!(
            are_3_vect_in_one_plane(
                &orthogonal,
                f64::from_bits(1.0_f64.to_bits().wrapping_add(1)),
            ),
            1
        );

        let planar = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        assert_eq!(are_3_vect_in_one_plane(&planar, -f64::from_bits(1)), 0);
        assert_eq!(are_3_vect_in_one_plane(&planar, -0.0), 1);
        assert_eq!(are_3_vect_in_one_plane(&planar, 0.0), 1);
        assert_eq!(are_3_vect_in_one_plane(&orthogonal, f64::NAN), 0);
    }

    #[test]
    fn source_port__ichister__are_4at_in_one_plane__line_1190() {
        let square = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        assert_eq!(are_4at_in_one_plane(&square, -f64::from_bits(1)), 0);
        assert_eq!(are_4at_in_one_plane(&square, -0.0), 1);
        assert_eq!(are_4at_in_one_plane(&square, 0.0), 1);

        let tetrahedron = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        assert_eq!(are_4at_in_one_plane(&tetrahedron, 0.0), 0);
        assert_eq!(are_4at_in_one_plane(&tetrahedron, 1.0), 1);
        assert_eq!(are_4at_in_one_plane(&tetrahedron, f64::NAN), 0);

        let mut nan_first = tetrahedron;
        nan_first[1][0] = f64::NAN;
        assert_eq!(are_4at_in_one_plane(&nan_first, f64::INFINITY), 1);
    }

    #[test]
    fn source_port__ichister__triple_prod_char__line_1218() {
        fn atom(x: f64, y: f64, z: f64, neighbor: u16, used_0d: i8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                x,
                y,
                z,
                bUsed0DParity: used_0d,
                ..inp_ATOM::default()
            };
            atom.neighbor[0] = neighbor;
            atom
        }

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![
                atom(0.0, 0.0, 0.0, 1, 0),
                atom(1.0, 0.0, 0.0, 0, 0),
                atom(1.0, 1.0, 0.0, 3, 0),
                atom(1.0, 0.0, 0.0, 2, 0),
            ])
            .unwrap();
        assert_eq!(
            triple_prod_char(
                &heap,
                atoms.as_const(),
                0,
                0,
                &[0, 0, 100],
                2,
                0,
                &[100, 0, 0]
            ),
            Ok(-71)
        );
        assert_eq!(
            triple_prod_char(
                &heap,
                atoms.as_const(),
                0,
                0,
                &[0, 0, 100],
                2,
                0,
                &[-100, 0, 0]
            ),
            Ok(71)
        );
        assert_eq!(
            triple_prod_char(
                &heap,
                atoms.as_const(),
                0,
                0,
                &[0, 0, 0],
                2,
                0,
                &[100, 0, 0]
            ),
            Ok(0)
        );

        let short_middle = heap
            .allocate_model_storage(vec![
                atom(0.0, 0.0, 0.0, 1, 0),
                atom(1.0, 0.0, 0.0, 0, 0),
                atom(0.0, 0.0, 0.0, 3, 0),
                atom(1.0, 0.0, 0.0, 2, 0),
            ])
            .unwrap();
        assert_eq!(
            triple_prod_char(
                &heap,
                short_middle.as_const(),
                0,
                0,
                &[0, 0, 100],
                2,
                0,
                &[100, 0, 0],
            ),
            Ok(0)
        );
        heap.slice_mut(short_middle).unwrap()[0].bUsed0DParity = 1;
        assert_eq!(
            triple_prod_char(
                &heap,
                short_middle.as_const(),
                0,
                0,
                &[0, 0, 100],
                2,
                0,
                &[100, 0, 0],
            ),
            Ok(-100)
        );

        assert_eq!(
            triple_prod_char(
                &heap,
                atoms.as_const(),
                0,
                0,
                &[i8::MIN, i8::MAX, 1],
                2,
                0,
                &[i8::MAX, i8::MIN, -1],
            ),
            Ok(0)
        );
        assert_eq!(
            triple_prod_char(&heap, atoms.as_const(), -1, 0, &[1, 0, 0], 2, 0, &[0, 1, 0],),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__bcanatombemiddleallene__line_1615() {
        for element in [
            &[b'C' as i8, 0][..],
            &[b'S' as i8, b'i' as i8, 0][..],
            &[b'G' as i8, b'e' as i8, 0][..],
        ] {
            assert_eq!(bCanAtomBeMiddleAllene(element, 0, 0), Ok(1));
            assert_eq!(
                bCanAtomBeMiddleAllene(element, 0, RADICAL_SINGLET as i8),
                Ok(1)
            );
            assert_eq!(bCanAtomBeMiddleAllene(element, 0, i8::MIN), Ok(0));
            assert_eq!(bCanAtomBeMiddleAllene(element, 0, i8::MAX), Ok(0));
            assert_eq!(bCanAtomBeMiddleAllene(element, 1, 0), Ok(0));
            assert_eq!(bCanAtomBeMiddleAllene(element, -1, 0), Ok(0));
        }
        assert_eq!(bCanAtomBeMiddleAllene(&[b'N' as i8, 0], 0, 0), Ok(0));
        assert_eq!(
            bCanAtomBeMiddleAllene(&[b'C' as i8, b'l' as i8, 0], 0, 0),
            Ok(0)
        );
        assert_eq!(bCanAtomBeMiddleAllene(&[0], 0, 0), Ok(0));
        assert_eq!(
            bCanAtomBeMiddleAllene(&[b'C' as i8], 0, 0),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichister__bcanatombeterminalallene__line_1705() {
        for element in [
            &[b'C' as i8, 0][..],
            &[b'S' as i8, b'i' as i8, 0][..],
            &[b'G' as i8, b'e' as i8, 0][..],
        ] {
            assert_eq!(bCanAtomBeTerminalAllene(element, 0, 0), Ok(1));
            assert_eq!(
                bCanAtomBeTerminalAllene(element, 0, RADICAL_SINGLET as i8),
                Ok(1)
            );
            assert_eq!(bCanAtomBeTerminalAllene(element, 0, i8::MIN), Ok(0));
            assert_eq!(bCanAtomBeTerminalAllene(element, 0, i8::MAX), Ok(0));
            assert_eq!(bCanAtomBeTerminalAllene(element, 1, 0), Ok(0));
            assert_eq!(bCanAtomBeTerminalAllene(element, -1, 0), Ok(0));
        }
        assert_eq!(bCanAtomBeTerminalAllene(&[b'N' as i8, 0], 0, 0), Ok(0));
        assert_eq!(
            bCanAtomBeTerminalAllene(&[b'C' as i8, b'l' as i8, 0], 0, 0),
            Ok(0)
        );
        assert_eq!(bCanAtomBeTerminalAllene(&[0], 0, 0), Ok(0));
        assert_eq!(
            bCanAtomBeTerminalAllene(&[b'C' as i8], 0, 0),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichister__gethalfstereobond0dparity__line_1725() {
        fn fixture(parities: &[i8], reverse_first: bool) -> (Vec<inp_ATOM>, Vec<AT_NUMB>) {
            let mut current = inp_ATOM {
                valence: parities.len() as i8,
                bUsed0DParity: 4,
                ..inp_ATOM::default()
            };
            let mut atoms = vec![current.clone()];
            let mut attachments = Vec::with_capacity(2 * parities.len());
            for (index, &parity) in parities.iter().enumerate() {
                let next_original = 200 + index as u16;
                let neighbor_original = 100 + index as u16;
                current.sb_parity[index] = parity;
                current.sb_ord[index] = index as i8;
                current.sn_orig_at_num[index] = neighbor_original;
                current.neighbor[index] = (index + 1) as u16;
                atoms.push(inp_ATOM {
                    valence: 1,
                    orig_at_number: next_original,
                    ..inp_ATOM::default()
                });
                if reverse_first && index == 0 {
                    attachments.extend([neighbor_original, next_original]);
                } else {
                    attachments.extend([next_original, neighbor_original]);
                }
            }
            atoms[0] = current;
            (atoms, attachments)
        }

        fn evaluate(
            atoms: Vec<inp_ATOM>,
            attachments: &[AT_NUMB],
            count: i32,
            bond_parity: i32,
            flag: i32,
        ) -> (Result<i32, SourceHeapError>, i8) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let result = GetHalfStereobond0DParity(
                &mut heap,
                atoms,
                0,
                attachments,
                count,
                bond_parity,
                flag,
            );
            let used = heap.slice(atoms.as_const()).unwrap()[0].bUsed0DParity;
            (result, used)
        }

        let (atoms, attachments) = fixture(&[AB_PARITY_ODD as i8], false);
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(AB_PARITY_ODD as i32), 6)
        );
        let (atoms, attachments) = fixture(&[AB_PARITY_ODD as i8], true);
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(AB_PARITY_EVEN as i32), 6)
        );
        let (atoms, attachments) = fixture(&[AB_PARITY_EVEN as i8], false);
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(AB_PARITY_EVEN as i32), 6)
        );
        let (atoms, attachments) = fixture(&[AB_PARITY_EVEN as i8], true);
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(AB_PARITY_ODD as i32), 6)
        );

        for parity in [AB_PARITY_UNKN as i8, AB_PARITY_UNDF as i8] {
            let (atoms, attachments) = fixture(&[parity], true);
            assert_eq!(
                evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
                (Ok(parity as i32), 6)
            );
        }

        assert_eq!(
            evaluate(vec![inp_ATOM::default()], &[], 0, 77, 2),
            (Ok(77), 0)
        );
        let (atoms, _) = fixture(&[AB_PARITY_ODD as i8], false);
        assert_eq!(evaluate(atoms, &[999], 1, 77, 2), (Ok(77), 4));
        let (atoms, attachments) = fixture(&[AB_PARITY_ODD as i8], false);
        assert_eq!(evaluate(atoms, &attachments, -1, 77, 2), (Ok(77), 4));

        let (mut atoms, attachments) = fixture(&[AB_PARITY_ODD as i8], false);
        atoms[1].valence = (MAX_NUM_STEREO_BONDS + 1) as i8;
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(77), 4)
        );
        let (mut atoms, attachments) = fixture(&[AB_PARITY_ODD as i8], false);
        atoms[1].orig_at_number = 0;
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(77), 4)
        );

        for (parities, expected) in [
            ([AB_PARITY_UNDF as i8, AB_PARITY_UNKN as i8], AB_PARITY_UNKN),
            ([AB_PARITY_UNDF as i8, AB_PARITY_ODD as i8], AB_PARITY_ODD),
            ([AB_PARITY_ODD as i8, AB_PARITY_UNDF as i8], AB_PARITY_ODD),
            ([AB_PARITY_ODD as i8, AB_PARITY_ODD as i8], AB_PARITY_ODD),
        ] {
            let (atoms, attachments) = fixture(&parities, false);
            assert_eq!(
                evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
                (Ok(expected as i32), 6)
            );
        }
        let (atoms, attachments) = fixture(&[AB_PARITY_ODD as i8, AB_PARITY_EVEN as i8], false);
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(77), 4)
        );

        let (mut atoms, attachments) = fixture(
            &[
                AB_PARITY_UNKN as i8,
                AB_PARITY_ODD as i8,
                AB_PARITY_EVEN as i8,
            ],
            false,
        );
        atoms[0].sb_parity[1] = 0;
        assert_eq!(
            evaluate(atoms, &attachments, attachments.len() as i32, 77, 2),
            (Ok(AB_PARITY_UNKN as i32), 6)
        );

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        assert_eq!(
            GetHalfStereobond0DParity(&mut heap, atoms, -1, &[], 0, 77, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            GetHalfStereobond0DParity(&mut heap, atoms, 1, &[], 0, 77, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let (mut malformed, _) = fixture(&[AB_PARITY_ODD as i8], false);
        malformed[0].sb_ord[0] = -1;
        let malformed = heap.allocate_model_storage(malformed).unwrap();
        assert_eq!(
            GetHalfStereobond0DParity(&mut heap, malformed, 0, &[], 0, 77, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let (atoms, _) = fixture(&[AB_PARITY_ODD as i8], false);
        let atoms = heap.allocate_model_storage(atoms).unwrap();
        assert_eq!(
            GetHalfStereobond0DParity(&mut heap, atoms, 0, &[], 1, 77, 2),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__half_stereo_bond_parity__line_2121() {
        fn atom(element: &[u8]) -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            for (destination, &source) in atom.elname.iter_mut().zip(element) {
                *destination = source as i8;
            }
            atom.elname[element.len()] = 0;
            atom
        }

        fn geometry(coordinates: &[[f64; 3]], num_h: i8) -> Vec<inp_ATOM> {
            let mut center = atom(b"C");
            center.valence = coordinates.len() as i8;
            center.num_H = num_h;
            let mut atoms = vec![center];
            for (index, &coordinate) in coordinates.iter().enumerate() {
                atoms[0].neighbor[index] = (index + 1) as u16;
                atoms.push(inp_ATOM {
                    x: coordinate[0],
                    y: coordinate[1],
                    z: coordinate[2],
                    orig_at_number: (index + 10) as u16,
                    ..inp_ATOM::default()
                });
            }
            atoms
        }

        fn evaluate(
            atoms: Vec<inp_ATOM>,
            cur_at: i32,
            removed_start: Option<usize>,
            num_removed_h: i32,
            initial_direction: [i8; 3],
            pointed: i32,
            unknown: i32,
        ) -> (Result<i32, SourceHeapError>, [i8; 3], Vec<inp_ATOM>) {
            let mut heap = SourceHeap::default();
            let atoms_pointer = heap.allocate_model_storage(atoms).unwrap();
            let removed_pointer = match removed_start {
                Some(index) => atoms_pointer.offset(index as i64).unwrap().as_const(),
                None => SourceConstPointer::null(),
            };
            let mut direction = initial_direction;
            let result = half_stereo_bond_parity(
                &mut heap,
                atoms_pointer,
                cur_at,
                removed_pointer,
                num_removed_h,
                Some(&mut direction),
                pointed,
                unknown,
            );
            let atoms = heap.slice(atoms_pointer.as_const()).unwrap().to_vec();
            (result, direction, atoms)
        }

        let (result, direction, _) = evaluate(
            vec![inp_ATOM::default()],
            0,
            None,
            0,
            [0, 0, 0],
            0,
            AB_PARITY_UNKN as i32,
        );
        assert_eq!(result, Ok(0));
        assert_eq!(direction, [0, 0, 100]);

        let mut too_many_h = atom(b"C");
        too_many_h.num_H = (NUM_H_ISOTOPES + 1) as i8;
        assert_eq!(
            evaluate(
                vec![too_many_h],
                0,
                None,
                0,
                [7, 8, 9],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(0)
        );
        for valence in [1, (MAX_NUM_STEREO_BOND_NEIGH + 1) as i8] {
            let mut center = atom(b"C");
            center.valence = valence;
            assert_eq!(
                evaluate(
                    vec![center],
                    0,
                    None,
                    0,
                    [1, 2, 3],
                    0,
                    AB_PARITY_UNKN as i32,
                )
                .0,
                Ok(0)
            );
        }
        for mut center in [atom(b"O"), atom(b"C")] {
            center.valence = 2;
            if center.elname[0] == b'C' as i8 {
                center.radical = 2;
            }
            assert_eq!(
                evaluate(
                    vec![center],
                    0,
                    None,
                    0,
                    [1, 2, 3],
                    0,
                    AB_PARITY_UNKN as i32,
                )
                .0,
                Ok(0)
            );
        }

        let mut duplicate_isotope = atom(b"C");
        duplicate_isotope.valence = 1;
        duplicate_isotope.num_H = 2;
        duplicate_isotope.num_iso_H[0] = 2;
        assert_eq!(
            evaluate(
                vec![duplicate_isotope],
                0,
                None,
                0,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(AB_PARITY_IISO as i32)
        );
        let mut duplicate_ordinary = atom(b"C");
        duplicate_ordinary.valence = 1;
        duplicate_ordinary.num_H = 2;
        assert_eq!(
            evaluate(
                vec![duplicate_ordinary],
                0,
                None,
                0,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(AB_PARITY_IISO as i32)
        );
        let mut negative_ordinary = atom(b"C");
        negative_ordinary.valence = 1;
        negative_ordinary.num_H = 1;
        negative_ordinary.num_iso_H[..2].copy_from_slice(&[1, 1]);
        assert_eq!(
            evaluate(
                vec![negative_ordinary],
                0,
                None,
                0,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(CT_ISO_H_ERR)
        );

        let mut center = atom(b"C");
        center.valence = 1;
        center.num_H = 2;
        center.num_iso_H[..2].copy_from_slice(&[1, 1]);
        let mut first_removed = inp_ATOM::default();
        first_removed.neighbor[0] = 0;
        first_removed.iso_atw_diff = 1;
        let second_removed = first_removed.clone();
        assert_eq!(
            evaluate(
                vec![center.clone(), first_removed, second_removed],
                0,
                Some(1),
                2,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(CT_ISO_H_ERR)
        );
        let mut invalid_removed = inp_ATOM::default();
        invalid_removed.neighbor[0] = 0;
        invalid_removed.iso_atw_diff = -1;
        assert_eq!(
            evaluate(
                vec![center, invalid_removed],
                0,
                Some(1),
                1,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(CT_ISO_H_ERR)
        );

        let mut implicit_pair = atom(b"C");
        implicit_pair.num_H = 3;
        implicit_pair.num_iso_H.copy_from_slice(&[1, 1, 1]);
        let mut removed = inp_ATOM::default();
        removed.neighbor[0] = 0;
        removed.iso_atw_diff = 1;
        assert_eq!(
            evaluate(
                vec![implicit_pair, removed],
                0,
                Some(1),
                1,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(-(AB_PARITY_UNDF as i32))
        );

        let mut ignored_isotopes = atom(b"C");
        ignored_isotopes.valence = 1;
        ignored_isotopes.num_H = 2;
        ignored_isotopes.cFlags = AT_FLAG_ISO_H_POINT as i8;
        assert_eq!(
            evaluate(
                vec![ignored_isotopes],
                0,
                None,
                0,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(-(AB_PARITY_UNDF as i32))
        );

        let mut one_attachment = atom(b"C");
        one_attachment.valence = 1;
        one_attachment.num_H = 1;
        assert_eq!(
            evaluate(
                vec![one_attachment],
                0,
                None,
                0,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(AB_PARITY_UNDF as i32)
        );
        let mut no_attachment = atom(b"C");
        no_attachment.num_H = 2;
        no_attachment.num_iso_H[..2].copy_from_slice(&[1, 1]);
        assert_eq!(
            evaluate(
                vec![no_attachment],
                0,
                None,
                0,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(-(AB_PARITY_UNDF as i32))
        );

        let planar = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, -1.0, 0.0]];
        let (odd, odd_direction, odd_atoms) = evaluate(
            geometry(&planar, 0),
            0,
            None,
            0,
            [0, 0, 0],
            0,
            AB_PARITY_UNKN as i32,
        );
        assert_eq!(odd, Ok(AB_PARITY_ODD as i32));
        assert_eq!(odd_direction, [0, 0, 100]);
        assert_eq!(odd_atoms[0].bAmbiguousStereo, 0);
        let reversed = [planar[0], planar[2], planar[1]];
        assert_eq!(
            evaluate(
                geometry(&reversed, 0),
                0,
                None,
                0,
                [0, 0, 0],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(AB_PARITY_EVEN as i32)
        );

        let mut either = geometry(&planar, 0);
        either[0].bond_stereo[0] = STEREO_SNGL_EITHER as i8;
        assert_eq!(evaluate(either, 0, None, 0, [1, 2, 3], 0, 77).0, Ok(77));
        let mut pointed = geometry(&planar, 0);
        pointed[0].bond_stereo[0] = STEREO_SNGL_UP as i8;
        let (pointed_parity, pointed_direction, _) =
            evaluate(pointed, 0, None, 0, [0, 0, 0], 1, 77);
        assert_eq!(pointed_parity, Ok(AB_PARITY_ODD as i32));
        assert_ne!(pointed_direction, [0, 0, 100]);

        let three_dimensional = [[1.0, 0.0, 0.25], [0.0, 1.0, -0.5], [-1.0, -1.0, 0.75]];
        let (parity_3d, direction_3d, _) = evaluate(
            geometry(&three_dimensional, 0),
            0,
            None,
            0,
            [0, 0, 0],
            0,
            77,
        );
        assert_eq!(parity_3d, Ok(AB_PARITY_ODD as i32));
        assert_ne!(direction_3d, [0, 0, 100]);

        let mut short = geometry(&[[0.0, 0.0, 0.0]; 3], 0);
        short[0].sb_parity[0] = AB_PARITY_ODD as i8;
        short[0].sb_ord[0] = 0;
        short[0].sn_orig_at_num[0] = short[2].orig_at_number;
        let (short_result, _, short_atoms) =
            evaluate(short, 0, None, 0, [1, 2, 3], 0, AB_PARITY_UNKN as i32);
        assert_eq!(short_result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(short_atoms[0].bUsed0DParity, FlagSB_0D as i8);

        let narrow_half_plane = [
            [1.0, 0.0, 0.0],
            [0.5, 0.8660254037844386, 0.0],
            [-0.5, 0.8660254037844386, 0.0],
        ];
        let (ambiguous_result, _, ambiguous_atoms) = evaluate(
            geometry(&narrow_half_plane, 0),
            0,
            None,
            0,
            [1, 2, 3],
            0,
            AB_PARITY_UNKN as i32,
        );
        assert_eq!(ambiguous_result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(ambiguous_atoms[0].bAmbiguousStereo, AMBIGUOUS_STEREO as i8);

        let near_opposite = [[1.0, 0.0, 0.0], [-0.05_f64.cos(), 0.05_f64.sin(), 0.0]];
        let (two_neighbor_result, _, two_neighbor_atoms) = evaluate(
            geometry(&near_opposite, 1),
            0,
            None,
            0,
            [1, 2, 3],
            0,
            AB_PARITY_UNKN as i32,
        );
        assert_eq!(two_neighbor_result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(
            two_neighbor_atoms[0].bAmbiguousStereo,
            AMBIGUOUS_STEREO as i8
        );
        let too_close_to_opposite = [
            [1.0, 0.0, 0.0],
            [-(MIN_SINE / 2.0).cos(), (MIN_SINE / 2.0).sin(), 0.0],
        ];
        let (undefined_result, _, undefined_atoms) = evaluate(
            geometry(&too_close_to_opposite, 1),
            0,
            None,
            0,
            [1, 2, 3],
            0,
            AB_PARITY_UNKN as i32,
        );
        assert_eq!(undefined_result, Ok(AB_PARITY_UNDF as i32));
        assert_eq!(undefined_atoms[0].bAmbiguousStereo, 0);
        let outside_ambiguous_angle = [
            [1.0, 0.0, 0.0],
            [
                -(MIN_ANGLE_DBOND + 0.01).cos(),
                (MIN_ANGLE_DBOND + 0.01).sin(),
                0.0,
            ],
        ];
        let (defined_result, _, defined_atoms) = evaluate(
            geometry(&outside_ambiguous_angle, 1),
            0,
            None,
            0,
            [1, 2, 3],
            0,
            AB_PARITY_UNKN as i32,
        );
        assert_eq!(defined_result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(defined_atoms[0].bAmbiguousStereo, 0);

        let mut nitrogen_geometry = geometry(&planar[..2], 0);
        nitrogen_geometry[0] = atom(b"N");
        nitrogen_geometry[0].valence = 2;
        nitrogen_geometry[0].neighbor[..2].copy_from_slice(&[1, 2]);
        assert_eq!(
            evaluate(
                nitrogen_geometry,
                0,
                None,
                0,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(AB_PARITY_ODD as i32)
        );

        let mut removed_geometry = geometry(&planar[..2], 1);
        let mut explicit_h = inp_ATOM {
            x: planar[2][0],
            y: planar[2][1],
            z: planar[2][2],
            orig_at_number: 12,
            ..inp_ATOM::default()
        };
        explicit_h.neighbor[0] = 0;
        removed_geometry.push(explicit_h);
        let removed_index = removed_geometry.len() - 1;
        assert_eq!(
            evaluate(
                removed_geometry,
                0,
                Some(removed_index),
                1,
                [0, 0, 0],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(AB_PARITY_ODD as i32)
        );

        let mut isotope_only_geometry = geometry(&planar[..1], 2);
        isotope_only_geometry[0].num_iso_H[..2].copy_from_slice(&[1, 1]);
        for (coordinate, isotope) in [(planar[1], 1_i8), (planar[2], 2_i8)] {
            let mut explicit_h = inp_ATOM {
                x: coordinate[0],
                y: coordinate[1],
                z: coordinate[2],
                iso_atw_diff: isotope,
                orig_at_number: (20 + isotope) as u16,
                ..inp_ATOM::default()
            };
            explicit_h.neighbor[0] = 0;
            isotope_only_geometry.push(explicit_h);
        }
        assert_eq!(
            evaluate(
                isotope_only_geometry,
                0,
                Some(2),
                2,
                [0, 0, 0],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Ok(-(AB_PARITY_ODD as i32))
        );

        let mut dangling = geometry(&planar, 0);
        dangling[0].neighbor[2] = u16::MAX;
        assert_eq!(
            evaluate(dangling, 0, None, 0, [1, 2, 3], 0, AB_PARITY_UNKN as i32,).0,
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut removed_count_overflow = geometry(&planar[..2], 1);
        let mut only_removed = inp_ATOM::default();
        only_removed.neighbor[0] = 0;
        removed_count_overflow.push(only_removed);
        assert_eq!(
            evaluate(
                removed_count_overflow,
                0,
                Some(3),
                2,
                [1, 2, 3],
                0,
                AB_PARITY_UNKN as i32,
            )
            .0,
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let mut heap = SourceHeap::default();
        let atoms = heap
            .allocate_model_storage(geometry(&planar[..2], 1))
            .unwrap();
        let mut explicit_h = inp_ATOM::default();
        explicit_h.neighbor[0] = 0;
        let removed = heap.allocate_model_storage(vec![explicit_h]).unwrap();
        assert_eq!(
            half_stereo_bond_parity(
                &mut heap,
                atoms,
                0,
                removed.as_const(),
                1,
                None,
                0,
                AB_PARITY_UNKN as i32,
            ),
            Err(SourceHeapError::PointerAllocationMismatch)
        );
        assert_eq!(
            half_stereo_bond_parity(
                &mut heap,
                SourceMutPointer::null(),
                0,
                SourceConstPointer::null(),
                0,
                None,
                0,
                AB_PARITY_UNKN as i32,
            ),
            Err(SourceHeapError::NullPointer)
        );
        assert_eq!(
            half_stereo_bond_parity(
                &mut heap,
                atoms,
                -1,
                SourceConstPointer::null(),
                0,
                None,
                0,
                AB_PARITY_UNKN as i32,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            half_stereo_bond_parity(
                &mut heap,
                atoms,
                3,
                SourceConstPointer::null(),
                0,
                None,
                0,
                AB_PARITY_UNKN as i32,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__save_a_stereo_bond__line_2542() {
        const COUNT: usize = MAX_NUM_STEREO_BONDS as usize;

        let mut neighbor1 = [0_u16; COUNT];
        let mut order1 = [11_i8; COUNT];
        let mut z_product1 = [12_i8; COUNT];
        let mut parity1 = [13_i8; COUNT];
        let mut neighbor2 = [0_u16; COUNT];
        let mut order2 = [21_i8; COUNT];
        let mut z_product2 = [22_i8; COUNT];
        let mut parity2 = [23_i8; COUNT];
        assert_eq!(
            save_a_stereo_bond(
                -7,
                AB_PARITY_EVEN as i32,
                3,
                4,
                &mut neighbor1,
                &mut order1,
                &mut z_product1,
                &mut parity1,
                8,
                9,
                &mut neighbor2,
                &mut order2,
                &mut z_product2,
                &mut parity2,
            ),
            1
        );
        let mut expected_neighbor1 = [0_u16; COUNT];
        expected_neighbor1[0] = 9;
        let mut expected_neighbor2 = [0_u16; COUNT];
        expected_neighbor2[0] = 4;
        let mut expected_order1 = [11_i8; COUNT];
        expected_order1[0] = 4;
        let mut expected_order2 = [21_i8; COUNT];
        expected_order2[0] = 9;
        let mut expected_z1 = [12_i8; COUNT];
        expected_z1[0] = -7;
        let mut expected_z2 = [22_i8; COUNT];
        expected_z2[0] = -7;
        let mut expected_parity1 = [13_i8; COUNT];
        expected_parity1[0] = AB_PARITY_EVEN as i8;
        let mut expected_parity2 = [23_i8; COUNT];
        expected_parity2[0] = AB_PARITY_EVEN as i8;
        assert_eq!(neighbor1, expected_neighbor1);
        assert_eq!(neighbor2, expected_neighbor2);
        assert_eq!(order1, expected_order1);
        assert_eq!(order2, expected_order2);
        assert_eq!(z_product1, expected_z1);
        assert_eq!(z_product2, expected_z2);
        assert_eq!(parity1, expected_parity1);
        assert_eq!(parity2, expected_parity2);

        neighbor1[..2].copy_from_slice(&[5, 6]);
        neighbor2[0] = 7;
        assert_eq!(
            save_a_stereo_bond(
                255,
                257,
                65_535,
                130,
                &mut neighbor1,
                &mut order1,
                &mut z_product1,
                &mut parity1,
                -1,
                -129,
                &mut neighbor2,
                &mut order2,
                &mut z_product2,
                &mut parity2,
            ),
            1
        );
        assert_eq!(neighbor1, [5, 6, 0]);
        assert_eq!(neighbor2, [7, 0, 0]);
        assert_eq!(order1[2], -126);
        assert_eq!(order2[1], 127);
        assert_eq!(z_product1[2], -1);
        assert_eq!(z_product2[1], -1);
        assert_eq!(parity1[2], 1);
        assert_eq!(parity2[1], 1);

        let mut full_neighbor1 = [1_u16; COUNT];
        let mut empty_neighbor2 = [0_u16; COUNT];
        let mut full_order1 = [31_i8; COUNT];
        let mut full_z1 = [32_i8; COUNT];
        let mut full_parity1 = [33_i8; COUNT];
        let mut empty_order2 = [41_i8; COUNT];
        let mut empty_z2 = [42_i8; COUNT];
        let mut empty_parity2 = [43_i8; COUNT];
        let before = (
            full_neighbor1,
            full_order1,
            full_z1,
            full_parity1,
            empty_neighbor2,
            empty_order2,
            empty_z2,
            empty_parity2,
        );
        assert_eq!(
            save_a_stereo_bond(
                1,
                2,
                3,
                4,
                &mut full_neighbor1,
                &mut full_order1,
                &mut full_z1,
                &mut full_parity1,
                5,
                6,
                &mut empty_neighbor2,
                &mut empty_order2,
                &mut empty_z2,
                &mut empty_parity2,
            ),
            0
        );
        assert_eq!(
            (
                full_neighbor1,
                full_order1,
                full_z1,
                full_parity1,
                empty_neighbor2,
                empty_order2,
                empty_z2,
                empty_parity2,
            ),
            before
        );

        let mut empty_neighbor1 = [0_u16; COUNT];
        let mut full_neighbor2 = [1_u16; COUNT];
        assert_eq!(
            save_a_stereo_bond(
                1,
                2,
                3,
                4,
                &mut empty_neighbor1,
                &mut full_order1,
                &mut full_z1,
                &mut full_parity1,
                5,
                6,
                &mut full_neighbor2,
                &mut empty_order2,
                &mut empty_z2,
                &mut empty_parity2,
            ),
            0
        );
        assert_eq!(empty_neighbor1, [0_u16; COUNT]);
        assert_eq!(full_neighbor2, [1_u16; COUNT]);
    }

    #[test]
    fn source_port__ichister__half_stereo_bond_action__line_2901() {
        let unknown = AB_PARITY_UNKN as i32;
        for (parity, bond_unknown, isotopic, expected) in [
            (AB_PARITY_NONE as i32, 0, 0, AB_PARITY_NONE as i32),
            (AB_PARITY_ODD as i32, 0, 0, AB_PARITY_CALC as i32),
            (AB_PARITY_EVEN as i32, 0, 0, AB_PARITY_CALC as i32),
            (-(AB_PARITY_ODD as i32), 0, 0, AB_PARITY_NONE as i32),
            (-(AB_PARITY_EVEN as i32), 0, 0, AB_PARITY_NONE as i32),
            (unknown, 0, 0, unknown),
            (-unknown, 0, 0, AB_PARITY_NONE as i32),
            (AB_PARITY_UNDF as i32, 0, 0, AB_PARITY_UNDF as i32),
            (-(AB_PARITY_UNDF as i32), 0, 0, AB_PARITY_NONE as i32),
            (AB_PARITY_IISO as i32, 0, 0, AB_PARITY_NONE as i32),
            (-(AB_PARITY_IISO as i32), 0, 0, -1),
            (AB_PARITY_ODD as i32, -1, 0, unknown),
            (AB_PARITY_UNDF as i32, 2, 0, unknown),
            (AB_PARITY_IISO as i32, i32::MAX, 0, AB_PARITY_NONE as i32),
            (AB_PARITY_ODD as i32, 0, 7, AB_PARITY_CALC as i32),
            (AB_PARITY_EVEN as i32, 0, -1, AB_PARITY_CALC as i32),
            (-(AB_PARITY_ODD as i32), 0, 1, AB_PARITY_CALC as i32),
            (-(AB_PARITY_EVEN as i32), 0, 1, AB_PARITY_CALC as i32),
            (unknown, 0, 1, unknown),
            (-unknown, 0, 1, unknown),
            (AB_PARITY_UNDF as i32, 0, 1, AB_PARITY_UNDF as i32),
            (-(AB_PARITY_UNDF as i32), 0, 1, AB_PARITY_UNDF as i32),
            (AB_PARITY_IISO as i32, 0, 1, AB_PARITY_NONE as i32),
            (-(AB_PARITY_IISO as i32), 0, 1, -1),
            (AB_PARITY_ODD as i32, 1, 1, unknown),
            (AB_PARITY_UNDF as i32, 1, 1, unknown),
            (AB_PARITY_IISO as i32, 1, 1, AB_PARITY_NONE as i32),
            (7, 0, 0, -1),
            (7, 0, 1, -1),
        ] {
            assert_eq!(
                half_stereo_bond_action(parity, bond_unknown, isotopic, unknown),
                expected,
                "parity={parity} bond_unknown={bond_unknown} isotopic={isotopic}"
            );
        }

        const AB_NEGATIVE: i32 = 0x10;
        const AB_UNKNOWN: i32 = 0x20;
        for (normalized, expected) in [
            (AB_PARITY_ODD as i32 | AB_NEGATIVE, AB_PARITY_CALC as i32),
            (AB_PARITY_ODD as i32 | AB_UNKNOWN, unknown),
            (AB_PARITY_UNDF as i32 | AB_UNKNOWN, unknown),
            (AB_PARITY_ODD as i32 | AB_UNKNOWN | AB_NEGATIVE, unknown),
            (AB_PARITY_UNDF as i32 | AB_UNKNOWN | AB_NEGATIVE, unknown),
            (AB_PARITY_IISO as i32 | AB_UNKNOWN, AB_PARITY_NONE as i32),
            (AB_PARITY_UNDF as i32 | AB_NEGATIVE, AB_PARITY_UNDF as i32),
        ] {
            assert_eq!(half_stereo_bond_action(normalized, 0, 1, unknown), expected);
        }
        for normalized in [
            AB_PARITY_ODD as i32 | AB_UNKNOWN | AB_NEGATIVE,
            AB_PARITY_ODD as i32 | AB_NEGATIVE,
            AB_PARITY_IISO as i32 | AB_UNKNOWN,
            AB_PARITY_UNDF as i32 | AB_NEGATIVE,
            AB_PARITY_UNDF as i32 | AB_UNKNOWN | AB_NEGATIVE,
        ] {
            assert_eq!(
                half_stereo_bond_action(normalized, 0, 0, unknown),
                AB_PARITY_NONE as i32
            );
        }

        assert_eq!(half_stereo_bond_action(77, 0, 0, 77), 77);
        assert_eq!(
            half_stereo_bond_action(-77, 0, 1, 77),
            77,
            "dynamic negative unknown is retained for isotopic stereo"
        );
        assert_eq!(
            half_stereo_bond_action(-77, 0, 0, 77),
            AB_PARITY_NONE as i32
        );
    }

    #[test]
    fn source_port__ichister__set_stereo_bonds_parity__line_3009() {
        fn carbon(x: f64, y: f64, ring: u16) -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: 6,
                x,
                y,
                nRingSystem: ring,
                ..inp_ATOM::default()
            };
            atom.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
            atom
        }

        fn alkene(stereo: i8) -> Vec<inp_ATOM> {
            let mut left = carbon(0.0, 0.0, 1);
            left.valence = 2;
            left.num_H = 1;
            left.chem_bonds_valence = 3;
            left.neighbor[..2].copy_from_slice(&[1, 2]);
            left.bond_type[..2].copy_from_slice(&[BOND_DOUBLE as u8, BOND_SINGLE as u8]);
            left.bond_stereo[0] = stereo;
            let mut right = carbon(1.0, 0.0, 2);
            right.valence = 2;
            right.num_H = 1;
            right.chem_bonds_valence = 3;
            right.neighbor[..2].copy_from_slice(&[0, 3]);
            right.bond_type[..2].copy_from_slice(&[BOND_DOUBLE as u8, BOND_SINGLE as u8]);
            right.bond_stereo[0] = stereo;
            let mut left_substituent = carbon(0.0, 1.0, 3);
            left_substituent.valence = 1;
            left_substituent.neighbor[0] = 0;
            left_substituent.bond_type[0] = BOND_SINGLE as u8;
            let mut right_substituent = carbon(1.0, -1.0, 4);
            right_substituent.valence = 1;
            right_substituent.neighbor[0] = 1;
            right_substituent.bond_type[0] = BOND_SINGLE as u8;
            vec![left, right, left_substituent, right_substituent]
        }

        fn evaluate(
            atoms: Vec<inp_ATOM>,
            output: Vec<sp_ATOM>,
            at_1: i32,
        ) -> (Result<i32, SourceHeapError>, Vec<inp_ATOM>, Vec<sp_ATOM>) {
            let mut heap = SourceHeap::default();
            let atoms_pointer = heap.allocate_model_storage(atoms).unwrap();
            let output_pointer = heap.allocate_model_storage(output).unwrap();
            let removed = heap.allocate_model_storage(Vec::<inp_ATOM>::new()).unwrap();
            let result = set_stereo_bonds_parity(
                &mut heap,
                output_pointer,
                atoms_pointer,
                at_1,
                removed.as_const(),
                0,
                0,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                8,
                0,
                AB_PARITY_UNKN as i32,
            );
            (
                result,
                heap.slice(atoms_pointer.as_const()).unwrap().to_vec(),
                heap.slice(output_pointer.as_const()).unwrap().to_vec(),
            )
        }

        let empty_output = vec![sp_ATOM::default(); 4];
        let (result, input, output) = evaluate(alkene(0), empty_output.clone(), 1);
        assert_eq!(result, Ok(1));
        assert_eq!(output[1].stereo_bond_neighbor, [1, 0, 0]);
        assert_eq!(output[0].stereo_bond_neighbor, [2, 0, 0]);
        assert_eq!(output[1].stereo_bond_neighbor2, [1, 0, 0]);
        assert_eq!(output[0].stereo_bond_neighbor2, [2, 0, 0]);
        assert_eq!(output[1].stereo_bond_ord, [0, 0, 0]);
        assert_eq!(output[0].stereo_bond_ord2, [0, 0, 0]);
        assert_eq!(output[1].stereo_bond_parity, [AB_PARITY_CALC as i8, 0, 0]);
        assert_eq!(output[0].stereo_bond_parity2, [AB_PARITY_CALC as i8, 0, 0]);
        assert_eq!(output[1].stereo_bond_z_prod, output[0].stereo_bond_z_prod);
        assert_eq!(output[1].stereo_bond_z_prod, output[1].stereo_bond_z_prod2);
        assert!(i32::from(output[1].stereo_bond_z_prod[0]).abs() >= MIN_DOT_PROD as i32);
        assert!((1..=2).contains(&output[0].parity));
        assert!((1..=2).contains(&output[1].parity));
        assert_eq!(output[0].parity2, output[0].parity);
        assert_eq!(output[1].parity2, output[1].parity);
        assert_ne!(output[0].z_dir, [0; 3]);
        assert_ne!(output[1].z_dir, [0; 3]);
        assert_eq!(input[0].bAmbiguousStereo, 0);
        assert_eq!(input[1].bAmbiguousStereo, 0);

        let (unknown_result, _, unknown_output) =
            evaluate(alkene(STEREO_DBLE_EITHER as i8), empty_output.clone(), 1);
        assert_eq!(unknown_result, Ok(1));
        assert_eq!(
            unknown_output[1].stereo_bond_parity,
            [AB_PARITY_UNKN as i8, 0, 0]
        );
        assert_eq!(
            unknown_output[1].stereo_bond_parity2,
            [AB_PARITY_UNKN as i8, 0, 0]
        );

        let mut full_output = empty_output.clone();
        for atom in &mut full_output[..2] {
            atom.stereo_bond_neighbor = [9, 10, 11];
            atom.stereo_bond_neighbor2 = [9, 10, 11];
        }
        let (full_result, _, full_after) = evaluate(alkene(0), full_output.clone(), 1);
        assert_eq!(full_result, Ok(0));
        assert_eq!(full_after, full_output);

        let mut rejected = alkene(0);
        rejected[1].c_point = 1;
        assert_eq!(evaluate(rejected, empty_output.clone(), 1).0, Ok(0));
        let mut rejected = alkene(0);
        rejected[1].valence = 1;
        rejected[1].num_H = 0;
        assert_eq!(evaluate(rejected, empty_output.clone(), 1).0, Ok(0));
        let mut rejected = alkene(0);
        rejected[1].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
        assert_eq!(evaluate(rejected, empty_output.clone(), 1).0, Ok(0));
        assert_eq!(
            evaluate(alkene(0), empty_output, -1).0,
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__fixsb0dparities__line_1824() {
        fn endpoint(parity: i8, order: i8, flag: i8, xyz: [f64; 3]) -> inp_ATOM {
            let mut atom = inp_ATOM {
                bUsed0DParity: flag,
                x: xyz[0],
                y: xyz[1],
                z: xyz[2],
                ..inp_ATOM::default()
            };
            atom.sb_parity[0] = parity;
            atom.sb_ord[0] = order;
            atom
        }
        fn run(
            atoms: Vec<inp_ATOM>,
            chain: i32,
            mut z1: [i8; 3],
            mut z2: [i8; 3],
            mut p1: i32,
            mut p2: i32,
        ) -> (i32, [i8; 3], [i8; 3], i32, i32) {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            let ret = FixSb0DParities(
                &heap,
                atoms.as_const(),
                chain,
                0,
                0,
                &mut z1,
                1,
                0,
                &mut z2,
                &mut p1,
                &mut p2,
            )
            .unwrap();
            (ret, z1, z2, p1, p2)
        }

        let no_descriptors = vec![inp_ATOM::default(), inp_ATOM::default()];
        assert_eq!(
            run(no_descriptors, 0, [1, 2, 3], [4, 5, 6], -1, 2),
            (0, [1, 2, 3], [4, 5, 6], -4, -4)
        );
        assert_eq!(
            run(
                vec![endpoint(1, 0, 0, [0.0; 3]), inp_ATOM::default()],
                0,
                [1, 2, 3],
                [4, 5, 6],
                1,
                2,
            ),
            (-1, [1, 2, 3], [4, 5, 6], 0, 0)
        );

        let well_atoms = vec![
            endpoint(1, 0, 0, [0.0, 0.0, 0.0]),
            endpoint(2, 0, 0, [0.0, 1.0, 0.0]),
        ];
        assert_eq!(
            run(well_atoms.clone(), 2, [1, 2, 3], [4, 5, 6], 1, 2),
            (0, [1, 2, 3], [4, 5, 6], 1, 2)
        );
        assert_eq!(
            run(well_atoms.clone(), 2, [1, 2, 3], [4, 5, 6], -1, 2),
            (0, [1, 2, 3], [4, 5, 6], -1, -2)
        );

        assert_eq!(
            run(well_atoms.clone(), 0, [1, 2, 3], [4, 5, 6], 3, 2),
            (-1, [1, 2, 3], [4, 5, 6], 3, 2)
        );
        assert_eq!(
            run(well_atoms.clone(), 0, [1, 2, 3], [4, 5, 6], 1, 4),
            (-1, [1, 2, 3], [4, 5, 6], 1, 4)
        );
        assert_eq!(
            run(
                vec![
                    endpoint(3, 0, 0, [0.0; 3]),
                    endpoint(4, 0, 0, [0.0, 1.0, 0.0]),
                ],
                0,
                [1, 2, 3],
                [4, 5, 6],
                4,
                3,
            ),
            (-1, [1, 2, 3], [4, 5, 6], 3, 3)
        );

        let both_wrong = vec![
            endpoint(1, 0, FlagSB_0D as i8, [0.0, 0.0, 0.0]),
            endpoint(2, 0, FlagSB_0D as i8, [0.0, 1.0, 0.0]),
        ];
        assert_eq!(
            run(both_wrong, 1, [1, 2, 3], [4, 5, 6], 1, 2),
            (0, [100, 0, 0], [0, 0, 100], 1, 2)
        );
        let short = vec![
            endpoint(1, 0, FlagSB_0D as i8, [0.0; 3]),
            endpoint(2, 0, 0, [0.0; 3]),
        ];
        assert_eq!(
            run(short, 1, [1, 2, 3], [0, 0, 100], 1, 2),
            (0, [100, 0, 0], [0, 0, 100], 1, 2)
        );
        let wrong1 = vec![
            endpoint(1, 0, FlagSB_0D as i8, [0.0, 0.0, 0.0]),
            endpoint(2, 0, 0, [0.0, 1.0, 0.0]),
        ];
        assert_eq!(
            run(wrong1, 1, [9, 9, 9], [0, 0, 100], 1, 2),
            (0, [100, 0, 0], [0, 0, 100], 1, 2)
        );
        let wrong2 = vec![
            endpoint(1, 0, 0, [0.0, 0.0, 0.0]),
            endpoint(2, 0, FlagSB_0D as i8, [0.0, 1.0, 0.0]),
        ];
        assert_eq!(
            run(wrong2, 1, [100, 0, 0], [9, 9, 9], 1, 2),
            (0, [100, 0, 0], [0, 0, 100], 1, 2)
        );
        assert_eq!(
            run(well_atoms, 1, [1, 2, 3], [4, 5, 6], 1, 2),
            (0, [1, 2, 3], [4, 5, 6], 1, 2)
        );
    }

    #[test]
    fn source_port__ichister__binpatomhasrequirdneigh__line_1317() {
        fn atom(name: &[u8], element: u8) -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: element,
                ..inp_ATOM::default()
            };
            for (destination, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *destination = source as i8;
            }
            atom.elname[name.len()] = 0;
            atom
        }

        fn evaluate(
            atoms: Vec<inp_ATOM>,
            current: i32,
            required: i32,
            double_bonds: i32,
            stereo_at_zz: i32,
        ) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            bInpAtomHasRequirdNeigh(
                &heap,
                atoms.as_const(),
                current,
                required,
                double_bonds,
                stereo_at_zz,
            )
        }

        let mut single_center = atom(b"C", 6);
        single_center.valence = 2;
        single_center.neighbor[..2].copy_from_slice(&[1, 2]);
        single_center.bond_type[..2]
            .copy_from_slice(&[BOND_SINGLE as u8 | BOND_MARK_ALL as u8, BOND_SINGLE as u8]);
        assert_eq!(
            evaluate(
                vec![single_center.clone(), atom(b"C", 6), atom(b"C", 6)],
                0,
                0,
                0,
                0,
            ),
            Ok(1)
        );
        let mut endpoint = single_center.clone();
        endpoint.endpoint = 1;
        assert_eq!(
            evaluate(vec![endpoint, atom(b"C", 6), atom(b"C", 6)], 0, 0, 0, 0),
            Ok(0)
        );
        let mut hydrogen = single_center.clone();
        hydrogen.num_H = 1;
        assert_eq!(
            evaluate(vec![hydrogen, atom(b"C", 6), atom(b"C", 6)], 0, 1, 0, 0),
            Ok(0)
        );

        let mut nitrogen_center = single_center.clone();
        let mut n = atom(b"N", 7);
        n.valence = 1;
        let mut nh = n.clone();
        nh.num_H = 1;
        assert_eq!(
            evaluate(vec![nitrogen_center.clone(), n, nh], 0, 2, 0, 1),
            Ok(0)
        );
        let mut c = atom(b"C", 6);
        c.valence = 1;
        let mut ch = c.clone();
        ch.num_H = 1;
        assert_eq!(evaluate(vec![nitrogen_center, c, ch], 0, 2, 0, 1), Ok(1));

        for pseudo in [b"Zz".as_slice(), b"Zy".as_slice()] {
            let mut center = atom(b"C", 6);
            center.valence = 1;
            center.neighbor[0] = 1;
            center.bond_type[0] = BOND_SINGLE as u8;
            assert_eq!(
                evaluate(vec![center.clone(), atom(pseudo, 0)], 0, 0, 0, 0),
                Ok(0)
            );
            assert_eq!(evaluate(vec![center, atom(pseudo, 0)], 0, 0, 0, 1), Ok(1));
        }

        for multiple in [BOND_DOUBLE, BOND_ALTERN, BOND_TAUTOM, BOND_ALT12NS] {
            let mut center = atom(b"C", 6);
            center.valence = 1;
            center.neighbor[0] = 1;
            center.bond_type[0] = multiple as u8 | BOND_MARK_ALL as u8;
            assert_eq!(
                evaluate(vec![center.clone(), atom(b"C", 6)], 0, 0, 1, 1),
                Ok(1)
            );
            assert_eq!(
                evaluate(vec![center.clone(), atom(b"C", 6)], 0, 0, 2, 1),
                Ok(0)
            );
            assert_eq!(evaluate(vec![center, atom(b"C", 6)], 0, 0, 0, 1), Ok(0));
        }
        let mut other = atom(b"C", 6);
        other.valence = 1;
        other.neighbor[0] = 1;
        other.bond_type[0] = 0;
        assert_eq!(evaluate(vec![other, atom(b"C", 6)], 0, 0, 1, 1), Ok(0));
        assert_eq!(
            evaluate(vec![inp_ATOM::default()], -1, 0, 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__bcaninpatombeastereocenter__line_1432() {
        const ELEMENTS: [&[u8]; 21] = [
            b"C", b"Si", b"Ge", b"Sn", b"As", b"B", b"S", b"S", b"S", b"S", b"Se", b"Se", b"Se",
            b"Se", b"N", b"N", b"N", b"P", b"P", b"P", b"As",
        ];
        const CHARGES: [i8; 21] = [
            0, 0, 0, 0, 1, -1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
        ];
        const NUM_BONDS: [i8; 21] = [
            4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 3, 3,
        ];
        const CHEM_VALENCE: [i8; 21] = [
            4, 4, 4, 4, 4, 4, 4, 6, 3, 5, 4, 6, 3, 5, 5, 4, 3, 4, 5, 3, 3,
        ];

        fn named(name: &[u8]) -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            for (destination, source) in atom.elname.iter_mut().zip(name.iter().copied()) {
                *destination = source as i8;
            }
            atom.elname[name.len()] = 0;
            atom
        }

        fn fixture(index: usize) -> Vec<inp_ATOM> {
            let mut current = named(ELEMENTS[index]);
            current.charge = CHARGES[index];
            current.valence = NUM_BONDS[index];
            current.chem_bonds_valence = CHEM_VALENCE[index];
            current.nRingSystem = 1;
            current.nNumAtInRingSystem = if index == 16 { 3 } else { 1 };
            let valence = usize::try_from(NUM_BONDS[index]).unwrap();
            let multiple = usize::try_from(CHEM_VALENCE[index] - NUM_BONDS[index]).unwrap();
            for i in 0..valence {
                current.neighbor[i] = (i + 1) as u16;
                current.bond_type[i] = if i < multiple {
                    BOND_DOUBLE as u8
                } else {
                    BOND_SINGLE as u8
                };
            }
            let mut atoms = vec![current];
            for _ in 0..valence {
                let mut neighbor = named(b"C");
                neighbor.nRingSystem = 1;
                atoms.push(neighbor);
            }
            if index == 16 {
                atoms[1].valence = 2;
                atoms[1].neighbor[..2].copy_from_slice(&[0, 2]);
                atoms[2].valence = 2;
                atoms[2].neighbor[..2].copy_from_slice(&[0, 1]);
            }
            atoms
        }

        fn evaluate(
            atoms: Vec<inp_ATOM>,
            current: i32,
            pointed: i32,
            stereo_at_zz: i32,
        ) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            bCanInpAtomBeAStereoCenter(&heap, atoms, current, pointed, stereo_at_zz)
        }

        let all_features = PES_BIT_PHOSPHINE_STEREO as i32 | PES_BIT_ARSINE_STEREO as i32;
        for index in 0..21 {
            assert_eq!(
                evaluate(fixture(index), 0, all_features, 0),
                Ok(i32::from(NUM_BONDS[index])),
                "table index {index}"
            );
        }

        assert_eq!(evaluate(fixture(19), 0, 0, 0), Ok(0));
        assert_eq!(
            evaluate(fixture(19), 0, PES_BIT_PHOSPHINE_STEREO as i32, 0),
            Ok(3)
        );
        assert_eq!(evaluate(fixture(20), 0, 0, 0), Ok(0));
        assert_eq!(
            evaluate(fixture(20), 0, PES_BIT_ARSINE_STEREO as i32, 0),
            Ok(3)
        );

        let mut no_ring = fixture(16);
        no_ring[1].valence = 1;
        no_ring[2].valence = 1;
        assert_eq!(evaluate(no_ring, 0, all_features, 0), Ok(0));
        let mut wrong_charge = fixture(0);
        wrong_charge[0].charge = 1;
        assert_eq!(evaluate(wrong_charge, 0, all_features, 0), Ok(0));
        let mut wrong_radical = fixture(0);
        wrong_radical[0].radical = 2;
        assert_eq!(evaluate(wrong_radical, 0, all_features, 0), Ok(0));
        let mut singlet = fixture(0);
        singlet[0].radical = 1;
        assert_eq!(evaluate(singlet, 0, all_features, 0), Ok(4));
        let mut wrong_valence = fixture(0);
        wrong_valence[0].chem_bonds_valence = 3;
        assert_eq!(evaluate(wrong_valence, 0, all_features, 0), Ok(0));
        let mut malformed = fixture(0);
        malformed[0].elname = [b'C' as i8; 6];
        assert_eq!(
            evaluate(malformed, 0, all_features, 0),
            Err(SourceHeapError::MissingNulTerminator)
        );
        assert_eq!(
            evaluate(vec![inp_ATOM::default()], -1, all_features, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__can_be_a_stereo_atom_with_isotopic_h__line_3690() {
        fn carbon() -> inp_ATOM {
            let mut atom = inp_ATOM::default();
            atom.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
            atom
        }
        fn fixture(number_of_hydrogens: i8) -> Vec<inp_ATOM> {
            let mut current = carbon();
            current.num_H = number_of_hydrogens;
            current.valence = 4_i8.wrapping_sub(number_of_hydrogens);
            current.chem_bonds_valence = current.valence;
            for i in 0..usize::try_from(current.valence).unwrap_or(0) {
                current.neighbor[i] = (i + 1) as u16;
                current.bond_type[i] = BOND_SINGLE as u8;
            }
            let mut atoms = vec![current];
            atoms.extend((0..atoms[0].valence).map(|_| carbon()));
            atoms
        }
        fn evaluate(atoms: Vec<inp_ATOM>, current: i32) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            can_be_a_stereo_atom_with_isotopic_H(&heap, atoms, current, 0, 0)
        }

        for hydrogen_count in 0..=NUM_H_ISOTOPES as i8 {
            assert_eq!(evaluate(fixture(hydrogen_count), 0), Ok(1));
        }
        assert_eq!(evaluate(fixture(4), 0), Ok(0));
        let mut wrong_element = fixture(1);
        wrong_element[0].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
        assert_eq!(evaluate(wrong_element, 0), Ok(0));
        let mut wrong_bond = fixture(1);
        wrong_bond[0].bond_type[0] = BOND_DOUBLE as u8;
        assert_eq!(evaluate(wrong_bond, 0), Ok(0));
        assert_eq!(
            evaluate(vec![inp_ATOM::default()], -1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__getstereocenter0dparity__line_3734() {
        fn evaluate(
            atom: inp_ATOM,
            j1: i32,
            mut neighbors: Vec<AT_NUMB>,
            flag: i32,
        ) -> Result<(i32, Vec<AT_NUMB>, inp_ATOM), SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(vec![atom]).unwrap();
            let parity = GetStereocenter0DParity(
                &mut CANON_GLOBALS::default(),
                &mut heap,
                atoms,
                0,
                j1,
                &mut neighbors,
                flag,
            )?;
            let atom = heap.slice(atoms.as_const())?[0].clone();
            Ok((parity, neighbors, atom))
        }

        let mut four = inp_ATOM {
            p_parity: AB_PARITY_ODD as i8,
            p_orig_at_num: [40, 10, 30, 20],
            bUsed0DParity: 2,
            ..inp_ATOM::default()
        };
        let (parity, neighbors, atom) = evaluate(four.clone(), 4, vec![30, 10, 40, 20], 1).unwrap();
        assert_eq!(parity, AB_PARITY_EVEN as i32);
        assert_eq!(neighbors, [10, 20, 30, 40]);
        assert_eq!(atom.bUsed0DParity, 3);

        four.orig_at_number = 99;
        four.p_parity = AB_PARITY_EVEN as i8;
        four.p_orig_at_num = [30, 99, 10, 20];
        four.bUsed0DParity = 4;
        let (parity, neighbors, atom) =
            evaluate(four.clone(), 3, vec![30, 10, 20, 777], 2).unwrap();
        assert_eq!(parity, AB_PARITY_ODD as i32);
        assert_eq!(neighbors, [10, 20, 30, 777]);
        assert_eq!(atom.bUsed0DParity, 6);

        for ill_defined in [AB_PARITY_UNKN, AB_PARITY_UNDF] {
            four.p_parity = ill_defined as i8;
            four.p_orig_at_num = [4, 3, 2, 1];
            four.orig_at_number = 99;
            four.bUsed0DParity = -128;
            let (parity, neighbors, atom) = evaluate(four.clone(), 4, vec![2, 4, 1, 3], 1).unwrap();
            assert_eq!(parity, ill_defined as i32);
            assert_eq!(neighbors, [1, 2, 3, 4]);
            assert_eq!(atom.bUsed0DParity, -127);
        }

        four.p_parity = AB_PARITY_ODD as i8;
        four.p_orig_at_num = [1, 2, 3, 4];
        four.bUsed0DParity = 8;
        let (parity, neighbors, atom) = evaluate(four.clone(), 4, vec![7, 6, 5, 4], 1).unwrap();
        assert_eq!(parity, AB_PARITY_NONE as i32);
        assert_eq!(neighbors, [4, 5, 6, 7]);
        assert_eq!(atom.bUsed0DParity, 8);

        for invalid_count in [i32::MIN, -1, 0, 2, 5, i32::MAX] {
            let (parity, neighbors, atom) =
                evaluate(four.clone(), invalid_count, vec![9, 8, 7, 6], 1).unwrap();
            assert_eq!(parity, AB_PARITY_NONE as i32);
            assert_eq!(neighbors, [9, 8, 7, 6]);
            assert_eq!(atom.bUsed0DParity, 8);
        }

        four.p_parity = 0;
        let (parity, neighbors, atom) = evaluate(four.clone(), 4, vec![9, 8, 7, 6], 1).unwrap();
        assert_eq!(parity, AB_PARITY_NONE as i32);
        assert_eq!(neighbors, [9, 8, 7, 6]);
        assert_eq!(atom.bUsed0DParity, 8);

        four.p_parity = AB_PARITY_ODD as i8;
        assert_eq!(
            evaluate(four.clone(), 4, vec![1, 2, 3], 1),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        let mut heap = SourceHeap::default();
        let atoms = heap.allocate_model_storage(vec![four]).unwrap();
        assert_eq!(
            GetStereocenter0DParity(
                &mut CANON_GLOBALS::default(),
                &mut heap,
                atoms,
                -1,
                4,
                &mut [1, 2, 3, 4],
                1,
            ),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__get2dtetrahedralambiguity__line_402() {
        fn point(angle: f64, z: f64) -> [f64; 3] {
            [angle.cos(), angle.sin(), z]
        }
        fn evaluate(
            values: &[(f64, f64)],
            add_explicit_neighbor: i32,
            fix_border: i32,
            minimum_angle: f64,
        ) -> i32 {
            let mut coordinates = [[0.0; 3]; MAX_NUM_STEREO_ATOM_NEIGH as usize];
            for (slot, &(angle, z)) in coordinates.iter_mut().zip(values) {
                *slot = point(angle, z);
            }
            Get2DTetrahedralAmbiguity(
                &mut CANON_GLOBALS::default(),
                &coordinates,
                add_explicit_neighbor,
                fix_border,
                minimum_angle,
            )
            .unwrap()
        }

        let pi = std::f64::consts::PI;
        let three_neighbor_cases: &[(&str, &[(f64, f64)], i32)] = &[
            ("all in plane", &[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)], 4),
            ("one up one down", &[(0.0, 1.0), (1.0, -1.0), (2.0, 0.0)], 6),
            (
                "one up asymmetric planes",
                &[(0.0, 1.0), (1.0, 0.0), (2.0, 0.0)],
                1,
            ),
            (
                "one up opposite planes",
                &[(0.0, 1.0), (0.5, 0.0), (0.5 + pi, 0.0)],
                6,
            ),
            ("two up no down", &[(0.0, 1.0), (1.0, 1.0), (2.0, 0.0)], 1),
            ("three up", &[(0.0, 1.0), (1.0, 1.0), (2.0, 1.0)], 1),
            (
                "two up short alpha",
                &[(0.0, 1.0), (1.0, 1.0), (2.0, -1.0)],
                6,
            ),
            (
                "two up long alpha",
                &[(0.0, 1.0), (3.2, 1.0), (4.5, -1.0)],
                1,
            ),
            (
                "two up middle alpha inside",
                &[(0.0, 1.0), (2.1, 1.0), (4.19, -1.0)],
                1,
            ),
            (
                "two up middle alpha outside",
                &[(0.0, 1.0), (2.1, 1.0), (3.0, -1.0)],
                6,
            ),
            (
                "minority up inversion",
                &[(0.0, 1.0), (1.0, -1.0), (2.0, -1.0)],
                6,
            ),
        ];
        for &(name, values, expected) in three_neighbor_cases {
            assert_eq!(evaluate(values, 1, 0, 0.1), expected, "{name}");
        }
        assert_eq!(
            evaluate(three_neighbor_cases[2].1, -7, i32::MIN, 0.1),
            three_neighbor_cases[2].2
        );

        let four_neighbor_cases: &[(&str, &[(f64, f64)], i32, i32)] = &[
            (
                "all in plane",
                &[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)],
                0,
                4,
            ),
            (
                "one up alternating down ambiguous",
                &[(0.0, 1.0), (0.5, 0.0), (2.0, -1.0), (0.5 + pi, 0.0)],
                0,
                6,
            ),
            (
                "one up nonalternating down",
                &[(0.0, 1.0), (1.0, -1.0), (2.0, 0.0), (3.0, 0.0)],
                0,
                3,
            ),
            (
                "one up no down",
                &[(0.0, 1.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)],
                0,
                3,
            ),
            (
                "two adjacent up legacy",
                &[(0.0, 1.0), (1.0, 1.0), (2.0, 0.0), (3.0, 0.0)],
                0,
                6,
            ),
            (
                "two separated up legacy",
                &[(0.0, 1.0), (1.0, 0.0), (2.0, 1.0), (3.0, 0.0)],
                0,
                3,
            ),
            (
                "two separated up fixed",
                &[(0.0, 1.0), (1.0, 0.0), (2.0, 1.0), (3.0, 0.0)],
                1,
                3,
            ),
            (
                "two adjacent up fixed undefined",
                &[(0.0, 1.0), (1.0, 1.0), (2.0, 0.0), (3.0, 0.0)],
                1,
                6,
            ),
            (
                "two adjacent up fixed zero inner gap",
                &[(0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (3.0, 0.0)],
                1,
                3,
            ),
            (
                "border fix joins equal directions",
                &[(0.0, 1.0), (1.0, 0.0), (1.0, 1.0), (2.0, 0.0)],
                1,
                3,
            ),
            (
                "three up",
                &[(0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 0.0)],
                0,
                3,
            ),
            (
                "four up",
                &[(0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 1.0)],
                0,
                6,
            ),
            (
                "minority up inversion",
                &[(0.0, 1.0), (1.0, -1.0), (2.0, -1.0), (3.0, -1.0)],
                0,
                3,
            ),
        ];
        for &(name, values, fix, expected) in four_neighbor_cases {
            assert_eq!(evaluate(values, 0, fix, 0.1), expected, "{name}");
        }
        assert_eq!(
            evaluate(four_neighbor_cases[8].1, 0, -1, 0.1),
            four_neighbor_cases[8].3
        );

        let mut nan_coordinates = [
            point(0.0, f64::NAN),
            point(1.0, 0.0),
            point(2.0, 0.0),
            point(3.0, 0.0),
        ];
        nan_coordinates[1][0] = f64::NAN;
        assert_eq!(
            Get2DTetrahedralAmbiguity(
                &mut CANON_GLOBALS::default(),
                &nan_coordinates,
                0,
                1,
                f64::NAN,
            ),
            Ok(T2D_UNDF as i32)
        );
    }

    #[test]
    fn source_port__ichister__triple_prod_and_min_abs_sine2__line_926() {
        let orthogonal = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut untouched_ambiguity = 7;
        assert_eq!(
            triple_prod_and_min_abs_sine2(
                &orthogonal,
                &[0.25, 0.25, 0.25],
                0,
                None,
                Some(&mut untouched_ambiguity),
                0.1,
            ),
            1.0
        );
        assert_eq!(untouched_ambiguity, 7);

        let mut minimum = -99.0;
        let mut ambiguity = 0;
        assert_eq!(
            triple_prod_and_min_abs_sine2(
                &orthogonal,
                &[0.25, 0.25, 0.25],
                0,
                Some(&mut minimum),
                Some(&mut ambiguity),
                0.1,
            ),
            1.0
        );
        assert_eq!(minimum.to_bits(), (1.0_f64 / 3.0_f64.sqrt()).to_bits());
        assert_eq!(ambiguity, 0);

        minimum = -99.0;
        ambiguity = 0;
        assert_eq!(
            triple_prod_and_min_abs_sine2(
                &orthogonal,
                &[2.0, 2.0, 2.0],
                0,
                Some(&mut minimum),
                Some(&mut ambiguity),
                0.1,
            ),
            1.0
        );
        assert_eq!(ambiguity, AMBIGUOUS_STEREO as i32);

        let planar = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        minimum = -99.0;
        ambiguity = 0;
        assert_eq!(
            triple_prod_and_min_abs_sine2(
                &planar,
                &[0.0; 3],
                0,
                Some(&mut minimum),
                Some(&mut ambiguity),
                0.1,
            ),
            0.0
        );
        assert_eq!(minimum.to_bits(), 0.0_f64.to_bits());
        assert_eq!(ambiguity, 0);

        minimum = -99.0;
        ambiguity = 0;
        assert_eq!(
            triple_prod_and_min_abs_sine2(
                &orthogonal,
                &[2.0, 2.0, 2.0],
                0,
                Some(&mut minimum),
                Some(&mut ambiguity),
                f64::NAN,
            ),
            1.0
        );
        assert_eq!(minimum.to_bits(), 1.0_f64.to_bits());
        assert_eq!(ambiguity, 0);

        let almost_planar = [[1.0, 0.0, 0.0], [-1.0, 0.01, 0.0], [0.0, 1.0, 0.01]];
        for added_explicit_neighbor in [0, 1, -1, i32::MAX] {
            minimum = -99.0;
            ambiguity = 0;
            let product = triple_prod_and_min_abs_sine2(
                &almost_planar,
                &[0.0; 3],
                added_explicit_neighbor,
                Some(&mut minimum),
                Some(&mut ambiguity),
                0.001,
            );
            assert!(product > 0.0);
            assert_eq!(minimum.to_bits(), 0.0005_f64.to_bits());
            assert_eq!(ambiguity, AMBIGUOUS_STEREO as i32);
        }

        minimum = -99.0;
        assert_eq!(
            triple_prod_and_min_abs_sine2(
                &orthogonal,
                &[0.25, 0.25, 0.25],
                0,
                Some(&mut minimum),
                None,
                0.1,
            ),
            1.0
        );
        assert_eq!(minimum.to_bits(), (1.0_f64 / 3.0_f64.sqrt()).to_bits());
    }

    #[test]
    fn source_port__ichister__set_stereo_atom_parity__line_3790() {
        #[derive(Debug)]
        struct Observation {
            result: Result<i32, SourceHeapError>,
            output: sp_ATOM,
            input: inp_ATOM,
        }

        fn carbon_center(coordinates: &[[f64; 3]], number_of_hydrogens: i8) -> Vec<inp_ATOM> {
            let mut center = inp_ATOM::default();
            center.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
            center.orig_at_number = 100;
            center.valence = i8::try_from(coordinates.len()).unwrap();
            center.chem_bonds_valence = center.valence;
            center.num_H = number_of_hydrogens;
            for index in 0..coordinates.len() {
                center.neighbor[index] = (index + 1) as AT_NUMB;
                center.bond_type[index] = BOND_SINGLE as u8;
            }
            let mut atoms = vec![center];
            for (index, coordinate) in coordinates.iter().enumerate() {
                let mut neighbor = inp_ATOM::default();
                neighbor.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
                neighbor.orig_at_number = 101 + index as AT_NUMB;
                [neighbor.x, neighbor.y, neighbor.z] = *coordinate;
                atoms.push(neighbor);
            }
            atoms
        }

        #[allow(clippy::too_many_arguments)]
        fn evaluate(
            atoms: Vec<inp_ATOM>,
            current: i32,
            removed_start: Option<usize>,
            number_of_removed_hydrogens: i32,
            pointed_edge_stereo: i32,
            unknown_parity: i32,
            loose_check: i32,
            stereo_at_zz: i32,
        ) -> Observation {
            let mut heap = SourceHeap::default();
            let atoms_pointer = heap.allocate_model_storage(atoms).unwrap();
            let removed_pointer = match removed_start {
                Some(start) => atoms_pointer
                    .offset(i64::try_from(start).unwrap())
                    .unwrap()
                    .as_const(),
                None => SourceConstPointer::null(),
            };
            let atom_count = heap.slice(atoms_pointer.as_const()).unwrap().len();
            let mut output_atoms = vec![sp_ATOM::default(); atom_count];
            if let Ok(index) = usize::try_from(current)
                && let Some(output) = output_atoms.get_mut(index)
            {
                output.parity = 91;
                output.parity2 = 92;
                output.stereo_atom_parity = 93;
                output.stereo_atom_parity2 = 94;
                output.bAmbiguousStereo = 64;
            }
            let output_pointer = heap.allocate_model_storage(output_atoms).unwrap();
            let result = set_stereo_atom_parity(
                &mut CANON_GLOBALS::default(),
                &mut heap,
                output_pointer,
                atoms_pointer,
                current,
                removed_pointer,
                number_of_removed_hydrogens,
                pointed_edge_stereo,
                unknown_parity,
                loose_check,
                stereo_at_zz,
            );
            let index = usize::try_from(current).unwrap_or(0);
            Observation {
                result,
                output: heap
                    .slice(output_pointer.as_const())
                    .unwrap()
                    .get(index)
                    .cloned()
                    .unwrap_or_default(),
                input: heap
                    .slice(atoms_pointer.as_const())
                    .unwrap()
                    .get(index)
                    .cloned()
                    .unwrap_or_default(),
            }
        }

        let tetrahedron = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let defined = evaluate(
            carbon_center(&tetrahedron, 0),
            0,
            None,
            0,
            0,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        assert_eq!(defined.result, Ok(AB_PARITY_EVEN as i32));
        assert_eq!(
            (
                defined.output.parity,
                defined.output.parity2,
                defined.output.stereo_atom_parity,
                defined.output.stereo_atom_parity2,
                defined.output.bAmbiguousStereo,
                defined.input.bAmbiguousStereo,
            ),
            (
                AB_PARITY_EVEN as i8,
                AB_PARITY_EVEN as i8,
                AB_PARITY_CALC as i8,
                AB_PARITY_CALC as i8,
                64,
                0,
            )
        );

        let mut transposed = carbon_center(&tetrahedron, 0);
        transposed[0].neighbor.swap(0, 1);
        let transposed = evaluate(transposed, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(transposed.result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(transposed.output.parity, AB_PARITY_ODD as i8);

        let implicit_geometry = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let implicit = evaluate(
            carbon_center(&implicit_geometry, 1),
            0,
            None,
            0,
            PES_BIT_FIX_SP3_BUG as i32,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        assert_eq!(implicit.result, Ok(AB_PARITY_EVEN as i32));
        assert_eq!(implicit.output.parity, AB_PARITY_EVEN as i8);

        let mut with_removed =
            carbon_center(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, -1.0, 1.0]], 1);
        let mut removed = inp_ATOM::default();
        removed.elname[..2].copy_from_slice(&[b'H' as i8, 0]);
        removed.valence = 1;
        removed.neighbor[0] = 0;
        removed.orig_at_number = 104;
        removed.z = -1.0;
        with_removed.push(removed);
        let removed = evaluate(with_removed, 0, Some(4), 1, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(removed.result, Ok(AB_PARITY_EVEN as i32));
        assert_eq!(removed.output.parity, AB_PARITY_EVEN as i8);

        let planar_coordinates = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let planar = evaluate(
            carbon_center(&planar_coordinates, 0),
            0,
            None,
            0,
            0,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        assert_eq!(planar.result, Ok(AB_PARITY_UNDF as i32));
        assert_eq!(
            (
                planar.output.parity,
                planar.output.parity2,
                planar.output.stereo_atom_parity,
                planar.output.stereo_atom_parity2,
            ),
            (
                AB_PARITY_UNDF as i8,
                AB_PARITY_UNDF as i8,
                AB_PARITY_UNDF as i8,
                AB_PARITY_UNDF as i8,
            )
        );

        let mut zero_dimensional = carbon_center(&planar_coordinates, 0);
        zero_dimensional[0].p_parity = AB_PARITY_ODD as i8;
        zero_dimensional[0].p_orig_at_num = [101, 102, 103, 104];
        let zero_dimensional =
            evaluate(zero_dimensional, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(zero_dimensional.result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(zero_dimensional.output.parity, AB_PARITY_ODD as i8);
        assert_eq!(
            zero_dimensional.input.bUsed0DParity & FlagSC_0D as i8,
            FlagSC_0D as i8
        );

        let mut either = carbon_center(&planar_coordinates, 0);
        either[0].bond_stereo[0] = STEREO_SNGL_EITHER as i8;
        let either = evaluate(either, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(either.result, Ok(AB_PARITY_UNKN as i32));
        assert_eq!(
            (either.output.parity, either.output.parity2),
            (AB_PARITY_UNKN as i8, AB_PARITY_UNKN as i8)
        );

        let mut wedges = carbon_center(&planar_coordinates, 0);
        wedges[0].bond_stereo[..4].copy_from_slice(&[
            STEREO_SNGL_UP as i8,
            STEREO_SNGL_DOWN as i8,
            STEREO_SNGL_UP as i8,
            STEREO_SNGL_DOWN as i8,
        ]);
        let wedges = evaluate(
            wedges,
            0,
            None,
            0,
            PES_BIT_FIX_SP3_BUG as i32,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        assert_eq!(wedges.result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(
            i32::from(wedges.output.bAmbiguousStereo) & AMBIGUOUS_STEREO as i32,
            0
        );
        assert_eq!(
            i32::from(wedges.input.bAmbiguousStereo) & AMBIGUOUS_STEREO as i32,
            0
        );

        let warning_coordinates =
            [0.0_f64, 1.0, 2.0, 3.0].map(|angle| [angle.cos(), angle.sin(), 0.0]);
        let mut warned_wedges = carbon_center(&warning_coordinates, 0);
        warned_wedges[0].bond_stereo[..4].copy_from_slice(&[
            STEREO_SNGL_UP as i8,
            STEREO_SNGL_DOWN as i8,
            0,
            0,
        ]);
        let warned_wedges = evaluate(warned_wedges, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(warned_wedges.result, Ok(AB_PARITY_ODD as i32));
        assert_eq!(
            i32::from(warned_wedges.output.bAmbiguousStereo) & AMBIGUOUS_STEREO as i32,
            AMBIGUOUS_STEREO as i32
        );
        assert_eq!(
            i32::from(warned_wedges.input.bAmbiguousStereo) & AMBIGUOUS_STEREO as i32,
            AMBIGUOUS_STEREO as i32
        );

        let mut short_bond = carbon_center(&tetrahedron, 0);
        [short_bond[1].x, short_bond[1].y, short_bond[1].z] = [1.0e-8, 0.0, 0.0];
        short_bond[0].p_parity = AB_PARITY_EVEN as i8;
        short_bond[0].p_orig_at_num = [101, 102, 103, 104];
        let short_bond = evaluate(short_bond, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(short_bond.result, Ok(AB_PARITY_EVEN as i32));
        assert_ne!(short_bond.input.bUsed0DParity & FlagSC_0D as i8, 0);

        let mut isotopic_only = carbon_center(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], 2);
        isotopic_only[0].num_iso_H[0] = 1;
        isotopic_only[0].num_iso_H[1] = 1;
        let isotopic_only = evaluate(isotopic_only, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(isotopic_only.result, Ok(AB_PARITY_UNDF as i32));
        assert_eq!(
            (isotopic_only.output.parity, isotopic_only.output.parity2),
            (AB_PARITY_NONE as i8, AB_PARITY_UNDF as i8)
        );

        let mut duplicate_isotope = carbon_center(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], 2);
        duplicate_isotope[0].num_iso_H[0] = 2;
        let duplicate_isotope = evaluate(
            duplicate_isotope,
            0,
            None,
            0,
            0,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        assert_eq!(duplicate_isotope.result, Ok(AB_PARITY_NONE as i32));
        assert_eq!(duplicate_isotope.output.parity2, AB_PARITY_NONE as i8);

        let mut invalid_element = carbon_center(&tetrahedron, 0);
        invalid_element[0].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
        let invalid = evaluate(invalid_element, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(invalid.result, Ok(AB_PARITY_NONE as i32));
        assert_eq!(
            (
                invalid.output.parity,
                invalid.output.parity2,
                invalid.output.stereo_atom_parity,
                invalid.output.stereo_atom_parity2,
                invalid.output.bAmbiguousStereo,
            ),
            (0, 0, 0, 0, 64)
        );

        let near_planar = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.03]];
        let mut strict_atoms = carbon_center(&near_planar, 1);
        strict_atoms[0].nNumAtInRingSystem = 3;
        let strict = evaluate(
            strict_atoms.clone(),
            0,
            None,
            0,
            0,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        let loose = evaluate(strict_atoms, 0, None, 0, 0, AB_PARITY_UNKN as i32, 1, 0);
        assert_eq!(strict.result, Ok(AB_PARITY_UNDF as i32));
        assert_eq!(loose.result, Ok(AB_PARITY_UNDF as i32));

        let mut nan_geometry = carbon_center(&tetrahedron, 0);
        nan_geometry[1].x = f64::NAN;
        let nan_geometry = evaluate(nan_geometry, 0, None, 0, 0, AB_PARITY_UNKN as i32, 0, 0);
        assert_eq!(nan_geometry.result, Ok(AB_PARITY_UNDF as i32));

        let negative_index = evaluate(
            carbon_center(&tetrahedron, 0),
            -1,
            None,
            0,
            0,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        assert_eq!(
            negative_index.result,
            Err(SourceHeapError::PointerOutOfBounds)
        );

        let excessive_removed_count = evaluate(
            carbon_center(&tetrahedron, 0),
            0,
            Some(1),
            i32::MAX,
            0,
            AB_PARITY_UNKN as i32,
            0,
            0,
        );
        assert_eq!(
            excessive_removed_count.result,
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__batomhasvalence3__line_1574() {
        for radical in [0, RADICAL_SINGLET as i8] {
            assert_eq!(bAtomHasValence3(&[b'N' as i8, 0], 0, radical), Ok(1));
            assert_eq!(
                bAtomHasValence3(&[b'N' as i8, 0, b'X' as i8], 0, radical),
                Ok(1)
            );
        }

        for radical in [i8::MIN, -1, 2, i8::MAX] {
            assert_eq!(bAtomHasValence3(&[b'N' as i8, 0], 0, radical), Ok(0));
        }
        for charge in [i8::MIN, -1, 1, i8::MAX] {
            assert_eq!(bAtomHasValence3(&[b'N' as i8, 0], charge, 0), Ok(0));
        }
        for element in [
            &[0_i8][..],
            &[b'n' as i8, 0][..],
            &[b'N' as i8, b'N' as i8, 0][..],
            &[b'C' as i8, 0][..],
        ] {
            assert_eq!(bAtomHasValence3(element, 0, 0), Ok(0));
        }
        assert_eq!(
            bAtomHasValence3(&[b'N' as i8], 0, 0),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichister__bcanatomhaveastereobond__line_1594() {
        let valid = [
            (&b"C\0"[..], 0_i8),
            (&b"Si\0"[..], 0_i8),
            (&b"Ge\0"[..], 0_i8),
            (&b"N\0"[..], 0_i8),
            (&b"N\0"[..], 1_i8),
        ];
        for (element, charge) in valid {
            let element = element.iter().map(|&byte| byte as i8).collect::<Vec<_>>();
            assert_eq!(bCanAtomHaveAStereoBond(&element, charge, 0), Ok(1));
            assert_eq!(
                bCanAtomHaveAStereoBond(&element, charge, RADICAL_SINGLET as i8),
                Ok(1)
            );
            for radical in [i8::MIN, -1, 2, i8::MAX] {
                assert_eq!(bCanAtomHaveAStereoBond(&element, charge, radical), Ok(0));
            }
        }

        for (element, charge) in [
            (&b"C\0"[..], 1_i8),
            (&b"Si\0"[..], -1_i8),
            (&b"Ge\0"[..], i8::MAX),
            (&b"N\0"[..], 2_i8),
            (&b"n\0"[..], 0_i8),
            (&b"Nn\0"[..], 0_i8),
            (&b"\0"[..], 0_i8),
        ] {
            let element = element.iter().map(|&byte| byte as i8).collect::<Vec<_>>();
            assert_eq!(bCanAtomHaveAStereoBond(&element, charge, 0), Ok(0));
        }

        assert_eq!(bCanAtomHaveAStereoBond(&[b'C' as i8, 0, 88], 0, 0), Ok(1));
        assert_eq!(
            bCanAtomHaveAStereoBond(&[b'C' as i8], 0, 0),
            Err(SourceHeapError::MissingNulTerminator)
        );
    }

    #[test]
    fn source_port__ichister__bissuitableheteroinpatom__line_1634() {
        fn evaluate(atom: inp_ATOM) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atom = heap.allocate_model_storage(vec![atom]).unwrap();
            bIsSuitableHeteroInpAtom(&heap, atom.as_const())
        }

        for element in [8_u8, 16, 34, 52] {
            let atom = inp_ATOM {
                el_number: element,
                chem_bonds_valence: 2,
                valence: 1,
                ..inp_ATOM::default()
            };
            assert_eq!(evaluate(atom.clone()), Ok(0));
            assert_eq!(
                evaluate(inp_ATOM {
                    radical: RADICAL_SINGLET as i8,
                    ..atom
                }),
                Ok(0)
            );
        }

        assert_eq!(
            evaluate(inp_ATOM {
                el_number: 7,
                chem_bonds_valence: 2,
                num_H: 1,
                valence: 1,
                ..inp_ATOM::default()
            }),
            Ok(1)
        );
        assert_eq!(
            evaluate(inp_ATOM {
                el_number: 7,
                chem_bonds_valence: 3,
                num_H: 0,
                valence: 2,
                radical: RADICAL_SINGLET as i8,
                ..inp_ATOM::default()
            }),
            Ok(1)
        );

        for atom in [
            inp_ATOM {
                el_number: 8,
                chem_bonds_valence: 2,
                valence: 1,
                charge: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 8,
                chem_bonds_valence: 2,
                valence: 1,
                radical: 2,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 6,
                chem_bonds_valence: 2,
                valence: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 8,
                chem_bonds_valence: 1,
                valence: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 8,
                chem_bonds_valence: 1,
                num_H: 1,
                valence: 2,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 7,
                chem_bonds_valence: 3,
                valence: 1,
                ..inp_ATOM::default()
            },
            inp_ATOM {
                el_number: 7,
                chem_bonds_valence: i8::MIN,
                num_H: i8::MAX,
                valence: 2,
                ..inp_ATOM::default()
            },
        ] {
            assert_eq!(evaluate(atom), Ok(-1));
        }

        assert_eq!(
            bIsSuitableHeteroInpAtom(&SourceHeap::default(), SourceConstPointer::null()),
            Err(SourceHeapError::NullPointer)
        );
    }

    #[test]
    fn source_port__ichister__bisoxide__line_1667() {
        fn evaluate(
            atoms: Vec<inp_ATOM>,
            cur_at: i32,
        ) -> (Result<i32, SourceHeapError>, Vec<inp_ATOM>) {
            let mut heap = SourceHeap::default();
            let atoms_pointer = heap.allocate_model_storage(atoms).unwrap();
            let result = bIsOxide(&mut heap, atoms_pointer, cur_at);
            let atoms = heap.slice(atoms_pointer.as_const()).unwrap().to_vec();
            (result, atoms)
        }

        fn bonded_pair(bond_type: u8, neighbor: inp_ATOM) -> Vec<inp_ATOM> {
            let mut current = inp_ATOM::default();
            current.valence = 1;
            current.neighbor[0] = 1;
            current.bond_type[0] = bond_type | BOND_MARK_ALL as u8;
            vec![current, neighbor]
        }

        let oxide = inp_ATOM {
            el_number: 8,
            valence: 1,
            ..inp_ATOM::default()
        };
        let (result, atoms) = evaluate(bonded_pair(BOND_TYPE_DOUBLE as u8, oxide.clone()), 0);
        assert_eq!(result, Ok(1));
        assert_eq!(atoms[0].bond_type[0], BOND_TYPE_DOUBLE as u8);

        for neighbor in [
            inp_ATOM {
                valence: 2,
                ..oxide.clone()
            },
            inp_ATOM {
                charge: 1,
                ..oxide.clone()
            },
            inp_ATOM {
                num_H: 1,
                ..oxide.clone()
            },
            inp_ATOM {
                radical: RADICAL_SINGLET as i8,
                ..oxide.clone()
            },
            inp_ATOM {
                el_number: 6,
                ..oxide.clone()
            },
        ] {
            let (result, atoms) = evaluate(bonded_pair(BOND_TYPE_DOUBLE as u8, neighbor), 0);
            assert_eq!(result, Ok(0));
            assert_eq!(atoms[0].bond_type[0], BOND_TYPE_DOUBLE as u8);
        }

        let permissive_endpoint = inp_ATOM {
            el_number: 16,
            valence: 1,
            charge: i8::MAX,
            num_H: i8::MIN,
            radical: i8::MAX,
            ..inp_ATOM::default()
        };
        for bond_type in [BOND_TAUTOM as u8, BOND_ALT12NS as u8] {
            let (result, atoms) = evaluate(bonded_pair(bond_type, permissive_endpoint.clone()), 0);
            assert_eq!(result, Ok(1));
            assert_eq!(atoms[0].bond_type[0], bond_type);
        }

        for neighbor in [
            inp_ATOM {
                valence: 0,
                ..oxide.clone()
            },
            inp_ATOM {
                el_number: u8::MAX,
                ..oxide.clone()
            },
        ] {
            let (result, _) = evaluate(bonded_pair(BOND_TAUTOM as u8, neighbor), 0);
            assert_eq!(result, Ok(0));
        }

        let mut current = inp_ATOM::default();
        current.valence = 3;
        current.neighbor[..3].copy_from_slice(&[u16::MAX, 1, u16::MAX]);
        current.bond_type[..3].copy_from_slice(&[
            1 | BOND_MARK_ALL as u8,
            BOND_ALT12NS as u8 | BOND_MARK_ALL as u8,
            1 | BOND_MARK_ALL as u8,
        ]);
        let (result, atoms) = evaluate(vec![current, oxide.clone()], 0);
        assert_eq!(result, Ok(1));
        assert_eq!(
            &atoms[0].bond_type[..3],
            &[1, BOND_ALT12NS as u8, 1 | BOND_MARK_ALL as u8]
        );

        let (result, atoms) = evaluate(bonded_pair(BOND_TYPE_DOUBLE as u8, oxide.clone()), 1);
        assert_eq!(result, Ok(0));
        assert_eq!(
            atoms[0].bond_type[0],
            BOND_TYPE_DOUBLE as u8 | BOND_MARK_ALL as u8
        );

        let mut negative_valence = inp_ATOM::default();
        negative_valence.valence = i8::MIN;
        negative_valence.bond_type[0] = BOND_MARK_ALL as u8 | BOND_TYPE_DOUBLE as u8;
        let (result, atoms) = evaluate(vec![negative_valence], 0);
        assert_eq!(result, Ok(0));
        assert_eq!(
            atoms[0].bond_type[0],
            BOND_MARK_ALL as u8 | BOND_TYPE_DOUBLE as u8
        );

        let (result, _) = evaluate(vec![inp_ATOM::default()], -1);
        assert_eq!(result, Err(SourceHeapError::PointerOutOfBounds));
        let (result, _) = evaluate(vec![inp_ATOM::default()], 1);
        assert_eq!(result, Err(SourceHeapError::PointerOutOfBounds));

        let mut dangling = inp_ATOM::default();
        dangling.valence = 1;
        dangling.neighbor[0] = u16::MAX;
        dangling.bond_type[0] = BOND_TYPE_DOUBLE as u8 | BOND_MARK_ALL as u8;
        let (result, atoms) = evaluate(vec![dangling], 0);
        assert_eq!(result, Err(SourceHeapError::PointerOutOfBounds));
        assert_eq!(atoms[0].bond_type[0], BOND_TYPE_DOUBLE as u8);
    }

    #[test]
    fn source_port__ichister__get_allowed_stereo_bond_type__line_2588() {
        let allowed = [
            BOND_SINGLE as i32,
            BOND_TYPE_DOUBLE as i32,
            BOND_ALTERN as i32,
            BOND_TAUTOM as i32,
        ];
        for bond_type in allowed {
            assert_eq!(get_allowed_stereo_bond_type(bond_type), bond_type);
            for mark in (0..=BOND_MARK_ALL).step_by(16) {
                assert_eq!(
                    get_allowed_stereo_bond_type(bond_type | mark as i32),
                    bond_type
                );
            }
        }

        for bond_type in 0_i32..=u8::MAX as i32 {
            let unmarked = bond_type & !(BOND_MARK_ALL as i32);
            let expected = if allowed.contains(&unmarked) {
                unmarked
            } else {
                0
            };
            assert_eq!(get_allowed_stereo_bond_type(bond_type), expected);
        }

        for bond_type in [
            i32::MIN,
            i32::MIN | BOND_MARK_ALL as i32,
            -1,
            0x100 | BOND_SINGLE as i32,
            0x100 | BOND_TAUTOM as i32 | BOND_MARK_ALL as i32,
            i32::MAX,
        ] {
            assert_eq!(get_allowed_stereo_bond_type(bond_type), 0);
        }
    }

    #[test]
    fn source_port__ichister__can_be_a_stereo_bond_with_isotopic_h__line_2638() {
        fn carbon() -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: 6,
                ..inp_ATOM::default()
            };
            atom.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
            atom
        }

        fn evaluate(
            atoms: Vec<inp_ATOM>,
            cur_at: i32,
            mode: INCHI_MODE,
        ) -> Result<i32, SourceHeapError> {
            let mut heap = SourceHeap::default();
            let atoms = heap.allocate_model_storage(atoms).unwrap();
            can_be_a_stereo_bond_with_isotopic_H(&mut heap, atoms, cur_at, mode)
        }

        fn ethene(mark: u8) -> Vec<inp_ATOM> {
            let mut left = carbon();
            left.valence = 2;
            left.neighbor[..2].copy_from_slice(&[1, 2]);
            left.bond_type[..2]
                .copy_from_slice(&[BOND_TYPE_DOUBLE as u8 | mark, BOND_SINGLE as u8]);
            let mut right = carbon();
            right.valence = 2;
            right.neighbor[..2].copy_from_slice(&[0, 3]);
            right.bond_type[..2]
                .copy_from_slice(&[BOND_TYPE_DOUBLE as u8 | mark, BOND_SINGLE as u8]);
            vec![left, right, carbon(), carbon()]
        }

        assert_eq!(evaluate(ethene(0), 1, 0), Ok(1));
        assert_eq!(evaluate(ethene(BOND_MARK_ALL as u8), 1, 0), Ok(1));
        assert_eq!(evaluate(ethene(0), 0, 0), Ok(0));

        let mut alternating_left = carbon();
        alternating_left.valence = 1;
        alternating_left.num_H = 1;
        alternating_left.neighbor[0] = 1;
        alternating_left.bond_type[0] = BOND_ALTERN as u8;
        let mut alternating_right = alternating_left.clone();
        alternating_right.neighbor[0] = 0;
        let alternating = vec![alternating_left, alternating_right];
        assert_eq!(evaluate(alternating.clone(), 1, 0), Ok(1));
        assert_eq!(
            evaluate(alternating, 1, CMODE_NO_ALT_SBONDS as INCHI_MODE),
            Ok(0)
        );

        let mut invalid = ethene(0);
        invalid[1].valence = 0;
        assert_eq!(evaluate(invalid, 1, 0), Ok(0));
        let mut invalid = ethene(0);
        invalid[1].elname[..2].copy_from_slice(&[b'O' as i8, 0]);
        assert_eq!(evaluate(invalid, 1, 0), Ok(0));
        let mut invalid_opposite = ethene(0);
        invalid_opposite[0].charge = 1;
        assert_eq!(evaluate(invalid_opposite, 1, 0), Ok(0));

        let mut mixed = carbon();
        mixed.valence = 2;
        mixed.neighbor[..2].copy_from_slice(&[0, 1]);
        mixed.bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_ALTERN as u8]);
        let mut double_endpoint = carbon();
        double_endpoint.valence = 1;
        double_endpoint.num_H = 1;
        double_endpoint.neighbor[0] = 2;
        double_endpoint.bond_type[0] = BOND_TYPE_DOUBLE as u8;
        let mut altern_endpoint = double_endpoint.clone();
        altern_endpoint.bond_type[0] = BOND_ALTERN as u8;
        assert_eq!(
            evaluate(vec![double_endpoint, altern_endpoint, mixed], 2, 0),
            Ok(0)
        );

        let mut endpoint = carbon();
        endpoint.valence = 2;
        endpoint.neighbor[..2].copy_from_slice(&[2, 1]);
        endpoint.bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_SINGLE as u8]);
        let mut one_bad = carbon();
        one_bad.valence = 3;
        one_bad.neighbor[..3].copy_from_slice(&[0, 3, 4]);
        one_bad.bond_type[..3].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_SINGLE as u8, 0]);
        assert_eq!(
            evaluate(
                vec![
                    endpoint.clone(),
                    carbon(),
                    one_bad.clone(),
                    carbon(),
                    carbon(),
                ],
                2,
                0,
            ),
            Ok(1)
        );
        one_bad.bond_type[1] = 0;
        assert_eq!(
            evaluate(vec![endpoint, carbon(), one_bad, carbon(), carbon()], 2, 0,),
            Ok(0)
        );

        let mut inner_bad = carbon();
        inner_bad.valence = 3;
        inner_bad.neighbor[..3].copy_from_slice(&[2, 3, 4]);
        inner_bad.bond_type[..3].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, 0, 0]);
        let mut current = carbon();
        current.valence = 2;
        current.neighbor[..2].copy_from_slice(&[0, 1]);
        current.bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_SINGLE as u8]);
        assert_eq!(
            evaluate(vec![inner_bad, carbon(), current, carbon(), carbon()], 2, 0,),
            Ok(1)
        );

        let mut middle = carbon();
        middle.valence = 2;
        middle.neighbor[..2].copy_from_slice(&[2, 0]);
        middle.bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_TYPE_DOUBLE as u8]);
        let mut cumulene_current = carbon();
        cumulene_current.valence = 2;
        cumulene_current.neighbor[..2].copy_from_slice(&[1, 3]);
        cumulene_current.bond_type[..2]
            .copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_SINGLE as u8]);
        assert_eq!(
            evaluate(vec![carbon(), middle, cumulene_current, carbon()], 2, 0,),
            Ok(1)
        );

        let mut dangling = carbon();
        dangling.valence = 2;
        dangling.neighbor[..2].copy_from_slice(&[u16::MAX, 0]);
        dangling.bond_type[..2].copy_from_slice(&[BOND_TYPE_DOUBLE as u8, BOND_SINGLE as u8]);
        assert_eq!(
            evaluate(vec![dangling], 0, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            evaluate(vec![carbon()], -1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
        assert_eq!(
            evaluate(vec![carbon()], 1, 0),
            Err(SourceHeapError::PointerOutOfBounds)
        );
    }

    #[test]
    fn source_port__ichister__fixunkn0dstereobonds__line_2006() {
        let mut heap = SourceHeap::default();
        assert_eq!(
            FixUnkn0DStereoBonds(&mut heap, SourceMutPointer::null(), 0).unwrap(),
            0
        );
        assert_eq!(
            FixUnkn0DStereoBonds(&mut heap, SourceMutPointer::null(), -1).unwrap(),
            0
        );

        let mut first = inp_ATOM::default();
        first.sb_parity = [AB_PARITY_UNKN as i8, 0, AB_PARITY_UNKN as i8];
        first.sb_ord = [2, 3, 4];
        let mut second = inp_ATOM::default();
        second.sb_parity = [1, AB_PARITY_UNKN as i8, AB_PARITY_UNKN as i8];
        second.sb_ord = [0, 19, 7];
        let atoms = heap.allocate(vec![first, second]).unwrap();
        assert_eq!(FixUnkn0DStereoBonds(&mut heap, atoms, 2).unwrap(), 3);

        let values = heap.slice(atoms.as_const()).unwrap();
        assert_eq!(values[0].bond_stereo[2], STEREO_DBLE_EITHER as i8);
        assert_eq!(values[0].bond_stereo[3], 0);
        assert_eq!(values[0].bond_stereo[4], 0);
        assert_eq!(values[1].bond_stereo[0], 0);
        assert_eq!(values[1].bond_stereo[19], STEREO_DBLE_EITHER as i8);
        assert_eq!(values[1].bond_stereo[7], STEREO_DBLE_EITHER as i8);
    }

    fn stereo_pair(current_parity: i8, next_parity: i8) -> Vec<inp_ATOM> {
        let mut current = inp_ATOM::default();
        current.valence = 2;
        current.neighbor[0] = 1;
        current.sb_ord[0] = 0;
        current.sn_ord[0] = 0;
        current.sb_parity[0] = current_parity;
        let mut next = inp_ATOM::default();
        next.valence = 2;
        next.neighbor[0] = 0;
        next.sb_ord[0] = 0;
        next.sn_ord[0] = 0;
        next.sb_parity[0] = next_parity;
        vec![current, next]
    }

    #[test]
    fn source_port__ichister__reconcilecmlincidentbondparities__line_4690() {
        let mut heap = SourceHeap::default();

        let mut ignored = inp_ATOM::default();
        ignored.valence = (MAX_NUM_STEREO_BONDS + 1) as i8;
        ignored.sb_parity[0] = 1;
        let ignored = heap.allocate(vec![ignored]).unwrap();
        let ignored_visited = heap.allocate(vec![0_i8]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, ignored, 0, -1, ignored_visited, 0,)
                .unwrap(),
            0
        );

        let no_parity = heap.allocate(vec![inp_ATOM::default()]).unwrap();
        let no_parity_visited = heap.allocate(vec![0_i8]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, no_parity, 0, -1, no_parity_visited, 0,)
                .unwrap(),
            1
        );

        let already = heap.allocate(stereo_pair(1, 1)).unwrap();
        let already_visited = heap.allocate(vec![10_i8, 0]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, already, 0, -1, already_visited, 0,)
                .unwrap(),
            2
        );

        let mut broken = inp_ATOM::default();
        broken.valence = 1;
        broken.neighbor[0] = 1;
        broken.sb_parity[0] = 1;
        let cannot_find = heap.allocate(vec![broken, inp_ATOM::default()]).unwrap();
        let cannot_find_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, cannot_find, 0, -1, cannot_find_visited, 0,).unwrap(),
            4
        );
        assert_eq!(heap.slice(cannot_find_visited.as_const()).unwrap()[0], 10);

        let mismatch = heap.allocate(stereo_pair(3, 4)).unwrap();
        let mismatch_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, mismatch, 0, -1, mismatch_visited, 0,)
                .unwrap(),
            3
        );

        let bypass = heap.allocate(stereo_pair(3, 3)).unwrap();
        let bypass_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, bypass, 0, -1, bypass_visited, 0,).unwrap(),
            0
        );
        assert_eq!(heap.slice(bypass_visited.as_const()).unwrap(), &[20, 0]);

        let normal = heap.allocate(stereo_pair(1, 1)).unwrap();
        let normal_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, normal, 0, -1, normal_visited, 0,).unwrap(),
            0
        );
        assert_eq!(heap.slice(normal_visited.as_const()).unwrap(), &[21, 21]);

        let reconcile = heap.allocate(stereo_pair(1, 1)).unwrap();
        let reconcile_visited = heap.allocate(vec![2_i8, 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, reconcile, 0, -1, reconcile_visited, 0,)
                .unwrap(),
            0
        );
        assert_eq!(heap.slice(reconcile.as_const()).unwrap()[0].sb_parity[0], 2);
        assert_eq!(heap.slice(reconcile.as_const()).unwrap()[1].sb_parity[0], 2);
        assert_eq!(heap.slice(reconcile_visited.as_const()).unwrap(), &[22, 22]);

        let mobius = heap.allocate(stereo_pair(1, 1)).unwrap();
        let mobius_visited = heap.allocate(vec![0_i8, 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, mobius, 0, -1, mobius_visited, 0,).unwrap(),
            5
        );

        let disconnected = heap
            .allocate(stereo_pair(1 << SB_PARITY_SHFT, 1 << SB_PARITY_SHFT))
            .unwrap();
        let disconnected_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(
                &mut heap,
                disconnected,
                0,
                -1,
                disconnected_visited,
                1,
            )
            .unwrap(),
            0
        );
        assert_eq!(
            heap.slice(disconnected_visited.as_const()).unwrap(),
            &[21, 21]
        );

        let completed = heap.allocate(stereo_pair(1, 2)).unwrap();
        let completed_visited = heap.allocate(vec![0_i8, 20]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(&mut heap, completed, 0, -1, completed_visited, 0,)
                .unwrap(),
            0
        );
        assert_eq!(heap.slice(completed_visited.as_const()).unwrap(), &[20, 20]);

        let mut high_valence_pair = stereo_pair(1, 2);
        high_valence_pair[1].valence = (MAX_NUM_STEREO_BONDS + 1) as i8;
        let high_valence = heap.allocate(high_valence_pair).unwrap();
        let high_valence_visited = heap.allocate(vec![0_i8; 2]).unwrap();
        assert_eq!(
            ReconcileCmlIncidentBondParities(
                &mut heap,
                high_valence,
                0,
                -1,
                high_valence_visited,
                0,
            )
            .unwrap(),
            0
        );
        assert_eq!(
            heap.slice(high_valence_visited.as_const()).unwrap(),
            &[20, 0]
        );
    }

    #[test]
    fn source_port__ichister__reconcileallcmlbondparities__line_4663() {
        let mut heap = SourceHeap::default();
        let empty = heap.allocate(Vec::<inp_ATOM>::new()).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, empty, 0, 0).unwrap(),
            0
        );

        let normal = heap.allocate(stereo_pair(1, 1)).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, normal, 2, 0).unwrap(),
            0
        );

        let mismatch = heap.allocate(stereo_pair(3, 4)).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, mismatch, 2, 0).unwrap(),
            3
        );

        let mut metal_pair = stereo_pair(3, 4);
        metal_pair[0].el_number = 25;
        metal_pair[1].el_number = 25;
        let metal = heap.allocate(metal_pair).unwrap();
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, metal, 2, 1).unwrap(),
            0
        );

        let malformed = {
            let mut atom = inp_ATOM::default();
            atom.valence = 1;
            atom.neighbor[0] = 1;
            atom.sb_parity[0] = 1;
            heap.allocate(vec![atom, inp_ATOM::default()]).unwrap()
        };
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, malformed, 2, 0).unwrap(),
            4
        );

        heap.fail_after_allocations(0);
        assert_eq!(
            ReconcileAllCmlBondParities(&mut heap, normal, 2, 0).unwrap(),
            -1
        );
    }

    #[test]
    fn source_port__ichister__get_opposite_sb_atom__line_4861() {
        let mut heap = SourceHeap::default();
        let mut start = inp_ATOM::default();
        start.neighbor[0] = 1;
        let mut opposite = inp_ATOM::default();
        opposite.sb_parity[0] = 1;
        opposite.sb_ord[0] = 2;
        opposite.neighbor[2] = 0;
        let direct = heap.allocate(vec![start.clone(), opposite]).unwrap();
        let (mut next, mut reverse, mut parity) = (-1, -1, -1);
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                direct,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            1
        );
        assert_eq!((next, reverse, parity), (1, 2, 0));

        let mut wrong = inp_ATOM::default();
        wrong.sb_parity[0] = 1;
        wrong.sb_ord[0] = 0;
        wrong.neighbor[0] = 7;
        let no_match = heap.allocate(vec![start.clone(), wrong]).unwrap();
        let (mut next, mut reverse, mut parity) = (-2, -3, -4);
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                no_match,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            0
        );
        assert_eq!((next, reverse, parity), (-2, -3, -4));

        let mut middle = inp_ATOM::default();
        middle.neighbor[0] = 0;
        middle.neighbor[1] = 2;
        middle.valence = 2;
        middle.chem_bonds_valence = (2 * BOND_TYPE_DOUBLE) as i8;
        let mut far = inp_ATOM::default();
        far.sb_parity[0] = 2;
        far.sb_ord[0] = 1;
        far.neighbor[1] = 1;
        let chain = heap.allocate(vec![start, middle, far]).unwrap();
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                chain,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            2
        );
        assert_eq!((next, reverse, parity), (2, 1, 0));

        let mut loop_atom = inp_ATOM::default();
        loop_atom.neighbor[0] = 1;
        loop_atom.neighbor[1] = 0;
        loop_atom.valence = 2;
        loop_atom.chem_bonds_valence = (2 * BOND_TYPE_DOUBLE) as i8;
        let cycle = heap.allocate(vec![loop_atom.clone(), loop_atom]).unwrap();
        assert_eq!(
            get_opposite_sb_atom(
                &heap,
                cycle,
                0,
                0,
                Some(&mut next),
                Some(&mut reverse),
                Some(&mut parity),
            )
            .unwrap(),
            0
        );
    }

    #[test]
    fn source_port__ichister__set_stereo_parity__line_4398() {
        fn carbon(x: f64, y: f64, z: f64, ring: u16) -> inp_ATOM {
            let mut atom = inp_ATOM {
                el_number: 6,
                x,
                y,
                z,
                nRingSystem: ring,
                ..inp_ATOM::default()
            };
            atom.elname[..2].copy_from_slice(&[b'C' as i8, 0]);
            atom
        }

        fn tetrahedron() -> Vec<inp_ATOM> {
            let coordinates = [
                [1.0, 1.0, 1.0],
                [-1.0, -1.0, 1.0],
                [-1.0, 1.0, -1.0],
                [1.0, -1.0, -1.0],
            ];
            let mut center = carbon(0.0, 0.0, 0.0, 1);
            center.orig_at_number = 100;
            center.valence = 4;
            center.chem_bonds_valence = 4;
            center.neighbor[..4].copy_from_slice(&[1, 2, 3, 4]);
            center.bond_type[..4].fill(BOND_SINGLE as u8);
            let mut atoms = vec![center];
            for (index, coordinate) in coordinates.into_iter().enumerate() {
                let mut neighbor = carbon(
                    coordinate[0],
                    coordinate[1],
                    coordinate[2],
                    (index + 2) as u16,
                );
                neighbor.orig_at_number = 101 + index as u16;
                atoms.push(neighbor);
            }
            atoms
        }

        fn alkene() -> Vec<inp_ATOM> {
            let mut left = carbon(0.0, 0.0, 0.0, 1);
            left.valence = 2;
            left.num_H = 1;
            left.chem_bonds_valence = 3;
            left.neighbor[..2].copy_from_slice(&[1, 2]);
            left.bond_type[..2].copy_from_slice(&[BOND_DOUBLE as u8, BOND_SINGLE as u8]);
            let mut right = carbon(1.0, 0.0, 0.0, 2);
            right.valence = 2;
            right.num_H = 1;
            right.chem_bonds_valence = 3;
            right.neighbor[..2].copy_from_slice(&[0, 3]);
            right.bond_type[..2].copy_from_slice(&[BOND_DOUBLE as u8, BOND_SINGLE as u8]);
            let mut left_substituent = carbon(0.0, 1.0, 0.0, 3);
            left_substituent.valence = 1;
            left_substituent.neighbor[0] = 0;
            left_substituent.bond_type[0] = BOND_SINGLE as u8;
            let mut right_substituent = carbon(1.0, -1.0, 0.0, 4);
            right_substituent.valence = 1;
            right_substituent.neighbor[0] = 1;
            right_substituent.bond_type[0] = BOND_SINGLE as u8;
            vec![left, right, left_substituent, right_substituent]
        }

        let mut empty_heap = SourceHeap::default();
        assert_eq!(
            set_stereo_parity(
                &mut CANON_GLOBALS::default(),
                &mut empty_heap,
                SourceMutPointer::null(),
                SourceMutPointer::null(),
                0,
                0,
                None,
                None,
                0,
                0,
                AB_PARITY_UNKN as i32,
                0,
                0,
            ),
            Ok(-1)
        );
        let one_atom = empty_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        assert_eq!(
            set_stereo_parity(
                &mut CANON_GLOBALS::default(),
                &mut empty_heap,
                one_atom,
                SourceMutPointer::null(),
                1,
                0,
                None,
                None,
                0,
                0,
                AB_PARITY_UNKN as i32,
                0,
                0,
            ),
            Ok(-1)
        );

        let mut clear_heap = SourceHeap::default();
        let clear_input = clear_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let mut dirty = sp_ATOM {
            parity: 11,
            parity2: 12,
            stereo_atom_parity: 13,
            z_dir: [14, 15, 16],
            ..sp_ATOM::default()
        };
        dirty.stereo_bond_neighbor.fill(21);
        dirty.stereo_bond_neighbor2.fill(22);
        dirty.stereo_bond_ord.fill(23);
        dirty.stereo_bond_ord2.fill(24);
        dirty.stereo_bond_z_prod.fill(25);
        dirty.stereo_bond_z_prod2.fill(26);
        dirty.stereo_bond_parity.fill(27);
        dirty.stereo_bond_parity2.fill(28);
        let clear_output = clear_heap.allocate_model_storage(vec![dirty]).unwrap();
        assert_eq!(
            set_stereo_parity(
                &mut CANON_GLOBALS::default(),
                &mut clear_heap,
                clear_input,
                clear_output,
                1,
                0,
                None,
                None,
                0,
                0,
                AB_PARITY_UNKN as i32,
                0,
                0,
            ),
            Ok(0)
        );
        let cleared = &clear_heap.slice(clear_output.as_const()).unwrap()[0];
        assert_eq!(cleared.parity, 0);
        assert_eq!(cleared.parity2, 0);
        assert_eq!(cleared.stereo_bond_neighbor, [0; 3]);
        assert_eq!(cleared.stereo_bond_neighbor2, [0; 3]);
        assert_eq!(cleared.stereo_bond_ord, [0; 3]);
        assert_eq!(cleared.stereo_bond_ord2, [0; 3]);
        assert_eq!(cleared.stereo_bond_z_prod, [0; 3]);
        assert_eq!(cleared.stereo_bond_z_prod2, [0; 3]);
        assert_eq!(cleared.stereo_bond_parity, [0; 3]);
        assert_eq!(cleared.stereo_bond_parity2, [0; 3]);
        assert_eq!(cleared.stereo_atom_parity, 0);
        assert_eq!(cleared.z_dir, [14, 15, 16]);

        let mut tetrahedral_heap = SourceHeap::default();
        let tetrahedral_input = tetrahedral_heap
            .allocate_model_storage(tetrahedron())
            .unwrap();
        let tetrahedral_output = tetrahedral_heap
            .allocate_model_storage(vec![sp_ATOM::default(); 5])
            .unwrap();
        let (mut max_atoms, mut max_bonds) = (-1, -1);
        assert_eq!(
            set_stereo_parity(
                &mut CANON_GLOBALS::default(),
                &mut tetrahedral_heap,
                tetrahedral_input,
                tetrahedral_output,
                5,
                0,
                Some(&mut max_atoms),
                Some(&mut max_bonds),
                0,
                0,
                AB_PARITY_UNKN as i32,
                0,
                0,
            ),
            Ok(1)
        );
        assert_eq!((max_atoms, max_bonds), (1, 0));
        let tetrahedral_center = &tetrahedral_heap
            .slice(tetrahedral_output.as_const())
            .unwrap()[0];
        assert_eq!(tetrahedral_center.parity, AB_PARITY_EVEN as i8);
        assert_eq!(tetrahedral_center.parity2, AB_PARITY_EVEN as i8);

        let mut bond_heap = SourceHeap::default();
        let bond_input = bond_heap.allocate_model_storage(alkene()).unwrap();
        let bond_output = bond_heap
            .allocate_model_storage(vec![sp_ATOM::default(); 4])
            .unwrap();
        let (mut max_atoms, mut max_bonds) = (-1, -1);
        assert_eq!(
            set_stereo_parity(
                &mut CANON_GLOBALS::default(),
                &mut bond_heap,
                bond_input,
                bond_output,
                4,
                0,
                Some(&mut max_atoms),
                Some(&mut max_bonds),
                0,
                0,
                AB_PARITY_UNKN as i32,
                0,
                0,
            ),
            Ok(0)
        );
        assert_eq!((max_atoms, max_bonds), (0, 1));
        let bond_output_values = bond_heap.slice(bond_output.as_const()).unwrap();
        assert_ne!(bond_output_values[0].stereo_bond_neighbor, [0; 3]);
        assert_ne!(bond_output_values[1].stereo_bond_neighbor, [0; 3]);

        for minimum_ring_size in [0_u64, 2] {
            let mut heap = SourceHeap::default();
            let input = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let output = heap
                .allocate_model_storage(vec![sp_ATOM::default()])
                .unwrap();
            heap.trace_source_allocations();
            let mode = minimum_ring_size << REQ_MODE_MIN_SB_RING_SHFT;
            assert_eq!(
                set_stereo_parity(
                    &mut CANON_GLOBALS::default(),
                    &mut heap,
                    input,
                    output,
                    1,
                    0,
                    None,
                    None,
                    mode,
                    0,
                    AB_PARITY_UNKN as i32,
                    0,
                    0,
                ),
                Ok(0)
            );
            assert_eq!(heap.source_allocation_calls(), 0);
            assert_eq!(heap.live_allocation_count(), 2);
        }

        let ring_mode = 8_u64 << REQ_MODE_MIN_SB_RING_SHFT;
        let mut ring_heap = SourceHeap::default();
        let ring_input = ring_heap
            .allocate_model_storage(vec![inp_ATOM::default()])
            .unwrap();
        let ring_output = ring_heap
            .allocate_model_storage(vec![sp_ATOM::default()])
            .unwrap();
        ring_heap.trace_source_allocations();
        assert_eq!(
            set_stereo_parity(
                &mut CANON_GLOBALS::default(),
                &mut ring_heap,
                ring_input,
                ring_output,
                1,
                0,
                None,
                None,
                ring_mode,
                0,
                AB_PARITY_UNKN as i32,
                0,
                0,
            ),
            Ok(0)
        );
        assert_eq!(ring_heap.source_allocation_calls(), 4);
        assert_eq!(ring_heap.live_allocation_count(), 2);

        for successful_allocations in 0..4 {
            let mut heap = SourceHeap::default();
            let input = heap
                .allocate_model_storage(vec![inp_ATOM::default()])
                .unwrap();
            let output = heap
                .allocate_model_storage(vec![sp_ATOM::default()])
                .unwrap();
            heap.fail_after_allocations(successful_allocations);
            assert_eq!(
                set_stereo_parity(
                    &mut CANON_GLOBALS::default(),
                    &mut heap,
                    input,
                    output,
                    1,
                    0,
                    None,
                    None,
                    ring_mode,
                    0,
                    AB_PARITY_UNKN as i32,
                    0,
                    0,
                ),
                Ok(CT_OUT_OF_RAM),
                "allocation failure after {successful_allocations} successful allocations"
            );
            assert_eq!(
                heap.live_allocation_count(),
                2,
                "cleanup after allocation failure {successful_allocations}"
            );
        }
    }
}
