use crate::source::base::ichi_io::inchi_fgetsLf;
use crate::source::base::util::{lrtrim, mystrncpy};
use crate::source_types::{
    FILE, INCHI_IOS_TYPE_FILE, INCHI_IOSTREAM, MOL_FMT_INPLINELEN, SourceHeap, SourceHeapError,
    SourceMutPointer,
};

fn source_fseek(
    heap: &mut SourceHeap,
    file: SourceMutPointer<FILE>,
    position: i64,
) -> Result<i32, SourceHeapError> {
    if file.is_null() || position < 0 {
        return Ok(1);
    }
    let file = heap
        .slice_mut(file)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if file.error {
        return Ok(1);
    }
    file.position = position as u64;
    file.eof = false;
    Ok(0)
}

fn source_ftell(heap: &SourceHeap, file: SourceMutPointer<FILE>) -> Result<i64, SourceHeapError> {
    if file.is_null() {
        return Ok(-1);
    }
    let file = heap
        .slice(file.as_const())?
        .first()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if file.error {
        return Ok(-1);
    }
    i64::try_from(file.position).map_err(|_| SourceHeapError::SourceIntegerOverflow)
}

fn source_fputs(
    heap: &mut SourceHeap,
    string: SourceMutPointer<i8>,
    file: SourceMutPointer<FILE>,
) -> Result<i32, SourceHeapError> {
    if file.is_null() {
        return Ok(-1);
    }
    let input = heap.slice(string.as_const())?;
    let length = input
        .iter()
        .position(|byte| *byte == 0)
        .ok_or(SourceHeapError::MissingNulTerminator)?;
    let input = input[..length]
        .iter()
        .map(|byte| *byte as u8)
        .collect::<Vec<_>>();
    let file = heap
        .slice_mut(file)?
        .first_mut()
        .ok_or(SourceHeapError::PointerOutOfBounds)?;
    if file.error {
        return Ok(-1);
    }
    let position =
        usize::try_from(file.position).map_err(|_| SourceHeapError::SourceIntegerOverflow)?;
    if file.bytes.len() < position {
        file.bytes.resize(position, 0);
    }
    let overwrite = input.len().min(file.bytes.len().saturating_sub(position));
    file.bytes[position..position + overwrite].copy_from_slice(&input[..overwrite]);
    file.bytes.extend_from_slice(&input[overwrite..]);
    file.position = file
        .position
        .checked_add(
            u64::try_from(input.len()).map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
        )
        .ok_or(SourceHeapError::SourceIntegerOverflow)?;
    Ok(0)
}

#[allow(non_snake_case)]
pub(crate) fn MolfileSaveCopy(
    heap: &mut SourceHeap,
    input: Option<&mut INCHI_IOSTREAM>,
    start: i64,
    end: i64,
    output: SourceMutPointer<FILE>,
    number: i64,
) -> Result<i32, SourceHeapError> {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:385 MolfileSaveCopy
    // INCHI✔️❌: int MolfileSaveCopy(INCHI_IOSTREAM *inp_file,
    // INCHI✔️❌:                     long fPtrStart,
    // INCHI✔️❌:                     long fPtrEnd,
    // INCHI✔️❌:                     FILE *outfile,
    // INCHI✔️❌:                     long num)
    // INCHI✔️❌: {
    // INCHI✔️❌:     char line[MOL_FMT_INPLINELEN], *p;
    // INCHI✔️❌:     long fPtr;
    // INCHI✔️❌:     int ret = 1;
    // INCHI✔️❌:     char szNumber[32];
    // INCHI✔️❌:
    // INCHI✔️❌:     if (inp_file->type == INCHI_IOS_TYPE_FILE)
    // INCHI✔️❌:     {
    // INCHI✔️❌:
    // INCHI✔️❌:         FILE *infile = inp_file->f;
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!infile)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (!outfile)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (fPtrStart < 0L && fPtrEnd <= fPtrStart)
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         if (0 != fseek(infile, fPtrStart, SEEK_SET))
    // INCHI✔️❌:         {
    // INCHI✔️❌:             return 1;
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         while (fPtrEnd > (fPtr = ftell(infile)) && fPtr >= 0L
    // INCHI✔️❌:                 && inchi_fgetsLf(line, sizeof(line) - 1, inp_file))
    // INCHI✔️❌:         {
    // INCHI✔️❌:
    // INCHI✔️❌:             line[sizeof(line) - 1] = '\0'; /*  unnecessary extra precaution */
    // INCHI✔️❌:
    // INCHI✔️❌:             if (fPtr == fPtrStart && num)
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 int len;
    // INCHI✔️❌:                 lrtrim(line, &len);
    // INCHI✔️❌:                 len = sprintf(szNumber, "#%ld%s", num, len ? "/" : "");
    // INCHI✔️❌:                 mystrncpy(line + len, line, sizeof(line) - len - 1);
    // INCHI✔️❌:                 memcpy(line, szNumber, len);
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             if (!strchr(line, '\n'))
    // INCHI✔️❌:             {
    // INCHI✔️❌:                 p = line + strlen(line);
    // INCHI✔️❌:                 p[0] = '\n';
    // INCHI✔️❌:                 p[1] = '\0';
    // INCHI✔️❌:             }
    // INCHI✔️❌:
    // INCHI✔️❌:             fputs(line, outfile);
    // INCHI✔️❌:         }
    // INCHI✔️❌:
    // INCHI✔️❌:         ret = fseek(infile, fPtrEnd, SEEK_SET);
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else if (inp_file->type == INCHI_IOS_TYPE_STRING)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:     else
    // INCHI✔️❌:     {
    // INCHI✔️❌:         ;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     return ret;
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: MolfileSaveCopy

    let input = input.ok_or(SourceHeapError::NullPointer)?;
    let mut result = 1;
    if input.type_ == INCHI_IOS_TYPE_FILE as i32 {
        if input.f.is_null() || output.is_null() || (start < 0 && end <= start) {
            return Ok(1);
        }
        if source_fseek(heap, input.f, start)? != 0 {
            return Ok(1);
        }

        let line = heap.allocate_model_storage(vec![0_i8; MOL_FMT_INPLINELEN as usize])?;
        loop {
            let position = source_ftell(heap, input.f)?;
            if end <= position || position < 0 {
                break;
            }
            if inchi_fgetsLf(heap, line, MOL_FMT_INPLINELEN as i32 - 1, Some(input))?.is_null() {
                break;
            }
            heap.slice_mut(line)?[MOL_FMT_INPLINELEN as usize - 1] = 0;
            if position == start && number != 0 {
                let mut trimmed_length = 0;
                lrtrim(heap, line, Some(&mut trimmed_length))?;
                let prefix = format!("#{number}{}", if trimmed_length != 0 { "/" } else { "" });
                let prefix_length = prefix.len();
                let shifted = line.offset(
                    i64::try_from(prefix_length)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                )?;
                mystrncpy(
                    heap,
                    shifted,
                    line.as_const(),
                    u32::try_from(MOL_FMT_INPLINELEN as usize - prefix_length - 1)
                        .map_err(|_| SourceHeapError::SourceIntegerOverflow)?,
                )?;
                for (destination, source) in heap.slice_mut(line)?[..prefix_length]
                    .iter_mut()
                    .zip(prefix.bytes())
                {
                    *destination = source as i8;
                }
            }
            let bytes = heap.slice_mut(line)?;
            let length = bytes
                .iter()
                .position(|byte| *byte == 0)
                .ok_or(SourceHeapError::MissingNulTerminator)?;
            if !bytes[..length].contains(&(b'\n' as i8)) {
                if length + 1 >= bytes.len() {
                    heap.free(line)?;
                    return Err(SourceHeapError::PointerOutOfBounds);
                }
                bytes[length] = b'\n' as i8;
                bytes[length + 1] = 0;
            }
            source_fputs(heap, line, output)?;
        }
        result = source_fseek(heap, input.f, end)?;
        heap.free(line)?;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source_types::{INCHI_IOS_TYPE_STRING, SourceFile};

    #[test]
    fn source_port__mol_fmt2__molfilesavecopy__line_385() {
        let mut heap = SourceHeap::default();
        let input_file = heap
            .allocate(vec![SourceFile {
                bytes: b"skip\n  name  \nline without lf".to_vec(),
                ..SourceFile::default()
            }])
            .unwrap();
        let output_file = heap
            .allocate(vec![SourceFile {
                bytes: b"old".to_vec(),
                position: 1,
                ..SourceFile::default()
            }])
            .unwrap();
        let mut stream = INCHI_IOSTREAM {
            f: input_file,
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut stream), 5, 30, output_file, 7),
            Ok(0)
        );
        let output = &heap.slice(output_file.as_const()).unwrap()[0];
        assert_eq!(output.bytes, b"o#7/name\nline without lf\n".to_vec());
        assert_eq!(output.position, 25);
        assert_eq!(heap.slice(input_file.as_const()).unwrap()[0].position, 30);

        heap.slice_mut(input_file).unwrap()[0].position = 20;
        heap.slice_mut(output_file).unwrap()[0] = SourceFile::default();
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut stream), 0, 5, output_file, 0),
            Ok(0)
        );
        assert_eq!(
            heap.slice(output_file.as_const()).unwrap()[0].bytes,
            b"skip\n".to_vec()
        );

        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut stream), -2, -3, output_file, 0),
            Ok(1)
        );
        let mut string_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_STRING as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut string_stream), 0, 1, output_file, 0),
            Ok(1)
        );
        let mut null_file_stream = INCHI_IOSTREAM {
            type_: INCHI_IOS_TYPE_FILE as i32,
            ..INCHI_IOSTREAM::default()
        };
        assert_eq!(
            MolfileSaveCopy(&mut heap, Some(&mut null_file_stream), 0, 1, output_file, 0),
            Ok(1)
        );
        assert_eq!(
            MolfileSaveCopy(
                &mut heap,
                Some(&mut stream),
                0,
                1,
                SourceMutPointer::null(),
                0
            ),
            Ok(1)
        );
        assert_eq!(
            MolfileSaveCopy(&mut heap, None, 0, 1, output_file, 0),
            Err(SourceHeapError::NullPointer)
        );
    }
}
