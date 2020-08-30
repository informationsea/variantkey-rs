use super::*;
use atoi::FromRadix16;
use std::io::{self, prelude::*};

#[test]
fn test_encode_variant_key() {
    assert_eq!(
        Ok(0x0807728e88e80000),
        encode_variant_key(b"1", 976157, b"T", b"C"),
    );
}

#[test]
fn test_decode_variant_key() {
    assert_eq!(
        decode_variant_key(0x0806b567a0fee000),
        Ok(Variant {
            chrom: b"1".to_vec(),
            position: 879311,
            reference: Some(b"TTTC".to_vec()),
            alternative: Some(b"T".to_vec()),
        })
    );
    assert_eq!(
        decode_variant_key(0x0807728e88e80000),
        Ok(Variant {
            chrom: b"1".to_vec(),
            position: 976157,
            reference: Some(b"T".to_vec()),
            alternative: Some(b"C".to_vec()),
        })
    );
}

#[test]
fn test_decode_binary_sequence() {
    assert_eq!(decode_binary_sequence(0, 1), b"A");
    assert_eq!(decode_binary_sequence(1, 1), b"C");
    assert_eq!(decode_binary_sequence(2, 1), b"G");
    assert_eq!(decode_binary_sequence(3, 1), b"T");

    assert_eq!(decode_binary_sequence(0, 2), b"AA");
    assert_eq!(decode_binary_sequence(1 << 2 | 2, 2), b"CG");
}

#[test]
fn test_encode_variant_key_all_clinvar() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = io::BufReader::new(flate2::read::MultiGzDecoder::new(std::fs::File::open(
        "testfiles/vk-test.txt.gz",
    )?));

    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        //println!("{}", line.trim());
        let elements: Vec<_> = line.trim().split(' ').collect();
        let expected = u64::from_radix_16(elements[0].as_bytes()).0;
        let pos = elements[2].parse::<u64>().unwrap();
        //println!("{} {} {} {}", elements[1], pos, elements[3], elements[4]);
        assert_eq!(
            Ok(expected),
            encode_variant_key(
                elements[1].as_bytes(),
                pos,
                elements[3].as_bytes(),
                elements[4].as_bytes()
            )
        );

        if elements[3].len() + elements[4].len() <= 11
            && is_simple_sequence(elements[3].as_bytes())
            && is_simple_sequence(elements[4].as_bytes())
        {
            assert_eq!(
                Ok(Variant {
                    chrom: elements[1].as_bytes().to_vec(),
                    position: pos,
                    reference: Some(elements[3].as_bytes().to_vec()),
                    alternative: Some(elements[4].as_bytes().to_vec())
                }),
                decode_variant_key(expected)
            )
        } else {
            assert_eq!(
                Ok(Variant {
                    chrom: elements[1].as_bytes().to_vec(),
                    position: pos,
                    reference: None,
                    alternative: None,
                }),
                decode_variant_key(expected)
            )
        }

        line.clear();
    }

    Ok(())
}

#[test]
fn test_simple_sequence() {
    assert_eq!(is_simple_sequence(b"AT"), true);
    assert_eq!(is_simple_sequence(b"CG"), true);
    assert_eq!(is_simple_sequence(b"CX"), false);
    assert_eq!(is_simple_sequence(b"G"), true);
}

#[test]
fn test_chromosome_number() {
    assert_eq!(Ok(1), chromosome_to_number(b"chr1"));
    assert_eq!(Ok(1), chromosome_to_number(b"1"));
    assert_eq!(Ok(2), chromosome_to_number(b"chr2"));
    assert_eq!(Ok(2), chromosome_to_number(b"2"));
    assert_eq!(Ok(22), chromosome_to_number(b"chr22"));
    assert_eq!(Ok(22), chromosome_to_number(b"22"));

    assert_eq!(Ok(23), chromosome_to_number(b"chrX"));
    assert_eq!(Ok(23), chromosome_to_number(b"X"));
    assert_eq!(Ok(24), chromosome_to_number(b"chrY"));
    assert_eq!(Ok(24), chromosome_to_number(b"Y"));
    assert_eq!(Ok(25), chromosome_to_number(b"chrM"));
    assert_eq!(Ok(25), chromosome_to_number(b"MT"));

    assert_eq!(
        Err(VariantKeyError::InvalidChromosome),
        chromosome_to_number(b"23")
    );
    assert_eq!(
        Err(VariantKeyError::InvalidChromosome),
        chromosome_to_number(b"chr23")
    );
    assert_eq!(
        Err(VariantKeyError::InvalidChromosome),
        chromosome_to_number(b"Z")
    );
}
