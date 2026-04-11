# REM Release Process

The authoritative REM version lives in the workspace manifest:

- `[workspace.package].version` in `Cargo.toml`

All first-party crates inherit that version through `version.workspace = true`.

## Release checklist

1. Update `Cargo.toml` workspace version.
2. Verify the workspace test suite passes:
   - `cargo test --quiet -p rem-core -p rem-config -p rem-mesh -p rem-materials -p rem-bc -p rem-electrostatic -p rem-magnetostatic -p rem-eigenmode -p rem-driven -p rem-transient -p rem-result -p rem-cli -p rem-parallel -p rem-touchstone -p rem-mom -p rem-bem -p rem-sbr -p rem-febi -p rem-ddm -p rem-planar`
3. Verify wasm build passes:
   - `cargo build --release --target wasm32-unknown-unknown -p rem-wasm`
4. Verify Yew app bundle passes:
   - `cd crates/yew-app && trunk build --release`
5. Update top-level capability or release notes if the public version changed.
6. Tag the release in git after verification.

## Notes

- CI mirrors the same three verification steps above.
- First-party test selection intentionally excludes vendored upstream crates whose own test suites are not release-gating for REM.
- Third-party vendored crates keep their own versions; only REM workspace crates inherit the shared version.