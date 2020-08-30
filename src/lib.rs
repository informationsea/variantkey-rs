//! Pure rust implementation of [variantkey](https://github.com/Genomicsplc/variantkey)
//!
//! ## Example
//! ```rust
//! use variantkey::{encode_variant_key, decode_variant_key, Variant};
//! assert_eq!(
//!     Ok(0x0807728e88e80000),
//!     encode_variant_key(b"1", 976157, b"T", b"C"),
//! );
//! assert_eq!(
//!     decode_variant_key(0x0806b567a0fee000),
//!     Ok(Variant {
//!         chrom: b"1".to_vec(),
//!         position: 879311,
//!         reference: Some(b"TTTC".to_vec()),
//!         alternative: Some(b"T".to_vec()),
//!     })
//! );
//! ```

mod hash;

use std::convert::TryInto;
use std::fmt;

/**
 * Variant Key Error
 */
#[derive(Debug, PartialEq, Eq)]
pub enum VariantKeyError {
    InvalidChromosome,
    InvalidPosition,
}

impl std::fmt::Display for VariantKeyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VariantKeyError::InvalidChromosome => write!(f, "Invalid Chromosome"),
            VariantKeyError::InvalidPosition => write!(f, "Invalid Position"),
        }
    }
}

impl std::error::Error for VariantKeyError {}

/**
 * Decoded variant key
 */
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub chrom: Vec<u8>,
    pub position: u64,
    pub reference: Option<Vec<u8>>,
    pub alternative: Option<Vec<u8>>,
}

/**
 * Decode variant key into Variant
 */
pub fn decode_variant_key(variant_key: u64) -> Result<Variant, VariantKeyError> {
    let ref_alt = variant_key & ((1 << 31) - 1);
    let pos = (variant_key >> 31) & ((1 << 28) - 1);
    let chrom = (variant_key >> (31 + 28)) & ((1 << 5) - 1);

    let (reference, alternative) = if (ref_alt & 1) == 1 {
        (None, None)
    } else {
        let ref_len = (ref_alt >> 27) & ((1 << 4) - 1);
        let alt_len = (ref_alt >> 23) & ((1 << 4) - 1);
        let ref_value = ref_alt >> (23 - ref_len * 2);
        let alt_value = ref_alt >> (23 - (ref_len + alt_len) * 2);
        (
            Some(decode_binary_sequence(ref_value, ref_len)),
            Some(decode_binary_sequence(alt_value, alt_len)),
        )
    };

    Ok(Variant {
        chrom: decode_chromosome(chrom as u8)?.to_vec(),
        position: pos,
        reference,
        alternative,
    })
}

/**
 * Encode a variant to variantkey
 */
pub fn encode_variant_key(
    chrom: &[u8],
    position: u64,
    reference: &[u8],
    alternative: &[u8],
) -> Result<u64, VariantKeyError> {
    let chrom_val: u64 = chromosome_to_number(chrom)?.into();
    if position > (1 << 28) {
        return Err(VariantKeyError::InvalidPosition);
    }

    let mut data = chrom_val;
    data <<= 28;
    data |= position;
    data <<= 31;

    let ref_alt = if reference.len() + alternative.len() > 11
        || !is_simple_sequence(reference)
        || !is_simple_sequence(alternative)
    {
        hash::encode_refalt_hash(reference, alternative).into()
    } else {
        let ref_len: u64 = (reference.len() & 0xf).try_into().unwrap();
        let alt_len: u64 = (alternative.len() & 0xf).try_into().unwrap();

        let (ref_data, ref_bit_len) = binary_sequence(reference);
        let (alt_data, alt_bit_len) = binary_sequence(alternative);
        let padding_bit_len = 31 - ref_bit_len - alt_bit_len - 8;

        let mut data: u64 = ref_len << 4 | alt_len;
        data <<= ref_bit_len;
        data |= ref_data;
        data <<= alt_bit_len;
        data |= alt_data;
        data <<= padding_bit_len;
        data
    };
    data |= ref_alt;

    Ok(data)
}

fn decode_binary_sequence(value: u64, len: u64) -> Vec<u8> {
    let mut result = Vec::new();

    for pos in 0..len {
        match (value >> ((len - pos - 1) * 2)) & 0x3 {
            0 => result.push(b'A'),
            1 => result.push(b'C'),
            2 => result.push(b'G'),
            3 => result.push(b'T'),
            _ => unreachable!(),
        }
    }

    result
}

fn binary_sequence(seq: &[u8]) -> (u64, usize) {
    let mut data: u64 = 0;
    for one in seq {
        data <<= 2;
        match one {
            b'A' => data |= 0,
            b'C' => data |= 1,
            b'G' => data |= 2,
            b'T' => data |= 3,
            _ => unreachable!(),
        }
    }
    (data, seq.len() * 2)
}

fn is_simple_sequence(seq: &[u8]) -> bool {
    for one in seq {
        match one {
            b'A' | b'T' | b'C' | b'G' => (),
            _ => return false,
        }
    }
    true
}

/**
 * Convert chromosome number into text
 */
pub fn decode_chromosome(chrom: u8) -> Result<&'static [u8], VariantKeyError> {
    match chrom {
        0 => Ok(b"NA"),
        1 => Ok(b"1"),
        2 => Ok(b"2"),
        3 => Ok(b"3"),
        4 => Ok(b"4"),
        5 => Ok(b"5"),
        6 => Ok(b"6"),
        7 => Ok(b"7"),
        8 => Ok(b"8"),
        9 => Ok(b"9"),
        10 => Ok(b"10"),
        11 => Ok(b"11"),
        12 => Ok(b"12"),
        13 => Ok(b"13"),
        14 => Ok(b"14"),
        15 => Ok(b"15"),
        16 => Ok(b"16"),
        17 => Ok(b"17"),
        18 => Ok(b"18"),
        19 => Ok(b"19"),
        20 => Ok(b"20"),
        21 => Ok(b"21"),
        22 => Ok(b"22"),
        23 => Ok(b"X"),
        24 => Ok(b"Y"),
        25 => Ok(b"MT"),
        _ => Err(VariantKeyError::InvalidChromosome),
    }
}

/**
 * Convert chromosome text into number
 */
pub fn chromosome_to_number(chrom: &[u8]) -> Result<u8, VariantKeyError> {
    let chrom = if chrom.starts_with(b"chr") {
        &chrom[3..]
    } else {
        chrom
    };

    match chrom {
        b"X" => Ok(23),
        b"Y" => Ok(24),
        b"MT" => Ok(25),
        b"M" => Ok(25),
        _ => {
            if let Some(x) = atoi::atoi(chrom) {
                if x < 23 {
                    Ok(x)
                } else {
                    Err(VariantKeyError::InvalidChromosome)
                }
            } else {
                Err(VariantKeyError::InvalidChromosome)
            }
        }
    }
}

#[cfg(test)]
mod tests;
