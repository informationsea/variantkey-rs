# variantkey-rs

Pure rust implementation of [variantkey](https://github.com/Genomicsplc/variantkey)

## Example
```rust
use variantkey::{encode_variant_key, decode_variant_key, Variant};
assert_eq!(
    Ok(0x0807728e88e80000),
    encode_variant_key(b"1", 976157, b"T", b"C"),
);
assert_eq!(
    decode_variant_key(0x0806b567a0fee000),
    Ok(Variant {
        chrom: b"1".to_vec(),
        position: 879311,
        reference: Some(b"TTTC".to_vec()),
        alternative: Some(b"T".to_vec()),
    })
);
```