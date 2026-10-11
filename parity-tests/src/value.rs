//! Canonical comparison bytes, shared with the tiny Python encoder.
//! Expected results need no decoding on the successful comparison path.
use crate::Result;
use std::io::{Cursor, Read};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    Null,
    Bool(bool),
    Int(i64),
    UInt(u64),
    Float(u64),
    Text(String),
    Bytes(Vec<u8>),
    List(Vec<Value>),
}

impl Value {
    pub fn float(value: f64) -> Self {
        Self::Float(value.to_bits())
    }

    pub fn encode(&self) -> Vec<u8> {
        let mut output = Vec::new();
        self.write(&mut output);
        output
    }

    fn write(&self, output: &mut Vec<u8>) {
        match self {
            Self::Null => output.push(0),
            Self::Bool(value) => output.extend([1, u8::from(*value)]),
            Self::Int(value) => {
                output.push(2);
                output.extend(value.to_le_bytes());
            }
            Self::UInt(value) => {
                output.push(3);
                output.extend(value.to_le_bytes());
            }
            Self::Float(bits) => {
                output.push(4);
                output.extend(bits.to_le_bytes());
            }
            Self::Text(value) => {
                output.push(5);
                length(value.len(), output);
                output.extend(value.as_bytes());
            }
            Self::Bytes(value) => {
                output.push(6);
                length(value.len(), output);
                output.extend(value);
            }
            Self::List(values) => {
                output.push(7);
                length(values.len(), output);
                for value in values {
                    value.write(output);
                }
            }
        }
    }
}

fn length(value: usize, output: &mut Vec<u8>) {
    output.extend(
        u32::try_from(value)
            .expect("value exceeds 4 GiB")
            .to_le_bytes(),
    );
}

/// Only diagnostics/tests decode values; bound recursion and allocation by input size.
pub fn decode(data: &[u8]) -> Result<Value> {
    fn read(cursor: &mut Cursor<&[u8]>, depth: usize) -> Result<Value> {
        if depth > 128 {
            return Err("value nesting exceeds 128".into());
        }
        let mut tag = [0];
        cursor.read_exact(&mut tag)?;
        let mut number = [0; 8];
        if (2..=4).contains(&tag[0]) {
            cursor.read_exact(&mut number)?;
        }
        Ok(match tag[0] {
            0 => Value::Null,
            1 => {
                cursor.read_exact(&mut tag)?;
                match tag[0] {
                    0 => Value::Bool(false),
                    1 => Value::Bool(true),
                    _ => return Err("invalid boolean".into()),
                }
            }
            2 => Value::Int(i64::from_le_bytes(number)),
            3 => Value::UInt(u64::from_le_bytes(number)),
            4 => Value::Float(u64::from_le_bytes(number)),
            tag @ 5..=7 => {
                let mut count = [0; 4];
                cursor.read_exact(&mut count)?;
                let count = u32::from_le_bytes(count) as usize;
                if count > cursor.get_ref().len() - cursor.position() as usize {
                    return Err("value length exceeds remaining bytes".into());
                }
                if tag == 7 {
                    let mut values = Vec::new();
                    for _ in 0..count {
                        values.push(read(cursor, depth + 1)?);
                    }
                    Value::List(values)
                } else {
                    let mut bytes = vec![0; count];
                    cursor.read_exact(&mut bytes)?;
                    if tag == 5 {
                        Value::Text(String::from_utf8(bytes)?)
                    } else {
                        Value::Bytes(bytes)
                    }
                }
            }
            _ => return Err("unknown value tag".into()),
        })
    }
    let mut cursor = Cursor::new(data);
    let value = read(&mut cursor, 0)?;
    if cursor.position() != data.len() as u64 {
        return Err("trailing value bytes".into());
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn floating_comparison_preserves_every_bit() {
        assert_ne!(Value::float(0.0).encode(), Value::float(-0.0).encode());
        assert_ne!(
            Value::float(1.0).encode(),
            Value::Float(1f64.to_bits() + 1).encode()
        );
        for bits in [0, 1, u64::MAX, 0x7ff8_0000_0000_1234, (-0f64).to_bits()] {
            let value = Value::Float(bits);
            assert_eq!(decode(&value.encode()).unwrap(), value);
        }
    }

    #[test]
    fn composite_values_roundtrip_without_json() {
        let value = Value::List(vec![
            Value::Null,
            Value::Int(-9),
            Value::UInt(u64::MAX),
            Value::Bool(true),
            Value::Text("λ".into()),
            Value::Bytes(vec![0, 255]),
        ]);
        assert_eq!(decode(&value.encode()).unwrap(), value);
    }

    #[test]
    fn truncated_or_oversized_values_fail_without_large_allocation() {
        assert!(decode(&[5, 255, 255, 255, 255]).is_err());
        assert!(decode(&[7, 255, 255, 255, 255]).is_err());
        assert!(decode(&[1, 4]).is_err());
        assert!(decode(&[0, 0]).is_err());
    }
}
