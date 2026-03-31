# Falak — Claude Code Instructions

## Project Identity

**Falak** (Arabic/Persian: sky/celestial sphere) — Orbital mechanics and celestial dynamics for AGNOS

- **Type**: Flat library crate
- **License**: GPL-3.0-only
- **MSRV**: 1.89
- **Version**: SemVer 0.1.0

## Consumers

kiran (orbital mechanics for space sims), joshua (celestial simulation), impetus (physics integration)

## Development Reference

See [CONTRIBUTING.md](CONTRIBUTING.md) for development process, coding standards, and key principles. See `docs/architecture/` for technical details.

## DO NOT

- **Do not commit or push** — the user handles all git operations (commit, push, tag)
- **NEVER use `gh` CLI** — use `curl` to GitHub API only

## Documentation Structure

```
Root files (required):
  README.md          — quick start, features, dependency stack, consumers, license
  CHANGELOG.md       — per-version changes (Added/Changed/Fixed/Removed)
  CLAUDE.md          — this file (Claude Code behavior constraints)
  CONTRIBUTING.md    — development process, principles, PR workflow
  SECURITY.md        — supported versions, scope, reporting
  CODE_OF_CONDUCT.md — Contributor Covenant
  LICENSE            — GPL-3.0

docs/ (required):
  architecture/
    overview.md      — module map, data flow, consumers, dependency stack
  development/
    roadmap.md       — completed items, backlog, future features (demand-gated), v1.0 criteria

docs/ (when earned — not scaffolded empty):
  adr/
    NNN-title.md     — architectural decision records (when non-obvious choices are made)
  development/
    threat-model.md  — attack surface, mitigations (when security-relevant)
    dependency-watch.md — deps to monitor for updates/CVEs
  guides/
    usage.md         — patterns, philosophy, code examples
    testing.md       — test count, coverage, testing patterns

ADR format:
  # NNN — Title
  ## Status: Accepted/Superseded
  ## Context: Why this decision was needed
  ## Decision: What we chose
  ## Consequences: Trade-offs, what changes
```

## CHANGELOG Format

Follow [Keep a Changelog](https://keepachangelog.com/):

```markdown
# Changelog

## [Unreleased]
### Added — new features
### Changed — changes to existing features
### Fixed — bug fixes
### Removed — removed features
### Security — vulnerability fixes
### Performance — benchmark-proven improvements (include numbers)

## [X.Y.Z] - YYYY-MM-DD
### Added
- **module_name** — what was added and why
### Changed
- item: old behavior → new behavior
### Fixed
- issue description (root cause → fix)
### Performance
- benchmark_name: before → after (−XX%)
```

Rules:
- Every PR/commit that changes behavior gets a CHANGELOG entry
- Performance claims MUST include benchmark numbers
- Breaking changes get a **Breaking** section with migration guide
- Group by module when multiple changes in one release
- Link to ADR if a change was driven by an architectural decision
