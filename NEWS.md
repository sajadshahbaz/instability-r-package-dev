# instability 0.1.1

- Document the complete UV/DES input contract, fixed scientific choices, and
  all returned tables and important columns.
- Add a runnable synthetic own-data example that does not use reference data.
- Add actionable errors for structurally empty inputs and missing required UV
  or DES conditions without changing valid-input scientific calculations.
- Preserve the single public API and the validated scientific behavior.

# instability 0.1.0

- Initial validated R-package release of the UV/DES Instability computation.
- Provides `compute_uv_des_instability(vst, metadata)`.
- Uses already-preprocessed VST expression data and metadata as formal inputs.
- Preserves the validated scientific behavior of the authoritative UV/DES
  implementation.
- Exact validation reproduced all eleven frozen UV/DES scientific reference
  tables byte-for-byte and by parsed semantic comparison.
