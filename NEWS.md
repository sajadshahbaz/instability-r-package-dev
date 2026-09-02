# instability 0.1.0

- Initial validated R-package release of the UV/DES Instability computation.
- Provides `compute_uv_des_instability(vst, metadata)`.
- Uses already-preprocessed VST expression data and metadata as formal inputs.
- Preserves the validated scientific behavior of the authoritative UV/DES
  implementation.
- Exact validation reproduced all eleven frozen UV/DES scientific reference
  tables byte-for-byte and by parsed semantic comparison.
