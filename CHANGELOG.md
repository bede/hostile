# Hostile EIT - Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.0rc2] - 2024-07-24

### Added

- `hostile-eit` entrypoint can be used, not expected to be picked up, added just in case.

### Changed

- Project references to `hostile-eit` to enable pushing to PyPi without conflicting with
  existing naming, `hostile` is still the entrypoint.
- GitHub workflow attempts to use an external re-usable action.

## [2.0.0rc1] - 2024-07-18

### Added

- This change log.
- GitHub workflow to build and push containers to the OCIR registry in the `gpasltd` tenancy.

### Changed

- Updated bucket reference to equivalent hosted in EIT tenancy.
- Formatted `README.md` and updated links and blurb to reflect fork intentions.
- Various updates to EIT Pathogena references.
- Updated Dockerfile to match changes used in GPAS (`hostile-s3`)

## [1.1.0] - 2024-04-10

[Unreleased](https://github.com/EIT-Pathogena/hostile-eit/compare/2.0.0rc1...HEAD)
[2.0.0rc1](https://github.com/EIT-Pathogena/hostile-eit/compare/2.0.0rc1...1.1.0)
[1.1.0](https://github.com/bede/hostile/releases/tag/1.1.0)