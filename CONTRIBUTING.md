# Contributing to Falak

## Getting Started

```bash
git clone https://github.com/MacCracken/falak.git
cd falak
rustup show          # ensure correct toolchain (see rust-toolchain.toml)
make check           # fmt + clippy + test + audit
```

## Development Process

### P(-1): Scaffold Hardening

Run before any new features to establish a solid baseline.

0. Read roadmap, CHANGELOG, and open issues — know what was intended before auditing what was built
1. Test + benchmark sweep of existing code
2. Cleanliness check (see below)
3. Get baseline benchmarks (`./scripts/bench-history.sh`)
4. Internal deep review — gaps, optimizations, security, logging/errors, docs
5. External research — domain completeness, missing capabilities, best practices, world-class accuracy
6. Cleanliness check — must be clean after review
7. Additional tests/benchmarks from findings
8. Post-review benchmarks — prove the wins
9. Repeat if heavy

### Work Loop

Continuous cycle for all development work.

1. Work phase — new features, roadmap items, bug fixes
2. Cleanliness check
3. Test + benchmark additions for new code
4. Run benchmarks (`./scripts/bench-history.sh`)
5. Internal review — performance, memory, security, throughput, correctness
6. Cleanliness check — must be clean after audit
7. Deeper tests/benchmarks from audit observations
8. Run benchmarks again — prove the wins
9. If audit heavy, return to step 5
10. Documentation — update CHANGELOG, roadmap, docs
11. Version check — VERSION, Cargo.toml, recipe all in sync
12. Return to step 1

### Cleanliness Check

Every cycle requires a clean pass through all of these:

```bash
cargo fmt --check
cargo clippy --all-features --all-targets -- -D warnings
cargo audit
cargo deny check
RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
```

### Task Sizing

- **Low/Medium effort** — batch freely, multiple items per work loop cycle.
- **Large effort** — small bites only. Break into sub-tasks, verify each before moving on. Never batch large items together.
- **If unsure** — treat it as large. Smaller bites are always safer than overcommitting.

### Refactoring

- Refactor when the code tells you to — duplication, unclear boundaries, performance bottlenecks.
- Never refactor speculatively. Wait for the third instance before extracting an abstraction.
- Refactoring is part of the work loop, not a separate phase. If review (step 5) reveals structural issues, refactor before moving to step 6.
- Every refactor must pass the same cleanliness + benchmark gates as new code.

## Key Principles

- **Never skip benchmarks.** Numbers don't lie. The CSV history is the proof.
- **Tests + benchmarks are the way.** Minimum 80%+ coverage target.
- **Own the stack.** If an AGNOS crate wraps an external lib, depend on the AGNOS crate.
- **No magic.** Every operation is measurable, auditable, traceable.
- **`#[non_exhaustive]`** on all public enums.
- **`#[must_use]`** on all pure functions.
- **`#[inline]`** on hot-path functions.
- **`write!` over `format!`** — avoid temporary allocations.
- **Cow over clone** — borrow when you can, allocate only when you must.
- **Vec arena over HashMap** — when indices are known, direct access beats hashing.
- **Feature-gate optional deps** — consumers pull only what they need.
- **tracing on all operations** — structured logging for audit trail.

## DO NOTs

- Do not add unnecessary dependencies — keep it lean.
- Do not `unwrap()` or `panic!()` in library code.
- Do not skip benchmarks before claiming performance improvements.
- Do not commit `target/` or `Cargo.lock` (library crate).
- Zero clippy warnings. Zero `println!`.
- Physics functions must document units in their doc comments.

## Code Style

- Follow `rustfmt` defaults (enforced by CI).
- Public API items must have doc comments.
- Use clear, imperative-mood commit messages:
  ```
  add Hohmann transfer delta-v computation
  fix Kepler equation convergence for high eccentricity
  ```

## Testing

```bash
cargo test --all-features       # all tests
make bench                      # benchmarks
make coverage                   # coverage report
```

- All public API changes must include unit tests.
- Integration tests go in `tests/`.
- Benchmarks go in `benches/` — run before/after performance-sensitive changes.
- Bug fixes require a regression test.
- Tests must validate against known physical values, not arbitrary numbers.

## Pull Request Process

1. **Fork and branch** — create a feature branch from `main`.
2. **Keep commits focused** — one logical change per commit.
3. **Write tests** — new features require tests; bug fixes require a regression test.
4. **Run CI locally** before pushing: `make check`
5. **Open a PR** against `main` with a clear description of the change.
6. **Address review feedback** — maintainers may request changes before merging.

## License

By contributing you agree that your contributions will be licensed under GPL-3.0-only.
