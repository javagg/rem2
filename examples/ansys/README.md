# Ansys Sidecar Template (Geometry by rmsh)

This folder provides a minimal sidecar metadata template for `rem --project <aedt> --format ansys`.

## Workflow

1. Use rmsh to generate `from_rmsh.msh` (or your mesh file).
2. Copy and edit `rem_ansys_template.json`.
3. Run conversion:

```powershell
cargo run -p rem-cli -- --project path/to/demo.aedt --format ansys --out_config output/demo.json --out_msh output/demo.msh
```

If `--out_msh` does not already exist, converter copies mesh from sidecar `mesh` field.

## Required fields

- `mesh` (when `--out_msh` does not exist)
- `ports[*].attributes` (non-empty)
- Valid sweep in `solver`: `freq_min > 0`, `freq_max >= freq_min`, `freq_step > 0`

## Notes

- Geometry is intentionally out of scope for this converter and must come from rmsh.
- This sidecar format is a P0 bootstrap and can be extended as importer coverage grows.
