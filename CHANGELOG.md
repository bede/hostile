# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).



## [Unreleased]

### Added

- Faster long read decontamination due to automatic caching of MMI indexes
- Support for streaming data input via stdin for both single and interleaved paired FASTQ
- Support for streaming data output via stdout for both single and interleaved paired FASTQ
- Automatic allocation of available cores between alignment and compression tasks
- Support for Illumina CASAVA 1.8+ read headers
- Ability to override remote index repository URL

### Changed

- Fixed bug in inverted mode `--invert` with paired reads
- Reorganised index-related functionality into `hostile index` subcommands
- `--offline` renamed to `--airplane`
- Improved warnings and error messages
- Pins Bowtie2>=2.5.4 and Minimap2>=2.28
- Removed space before `/1` and `/2` when writing paired FASTQ



## [1.1.0] - 2024-04-10

### Added
- Allows cross-platform override of the default index cache directory using the environment variable `HOSTILE_CACHE_DIR`. Previously, overrides using `XDG_DATA_HOME` worked only on Linux. The default remains a platform-specific directory determined by the `platformdirs` library. (#32)

### Changed
- Dockerfile simplification (suggested by @bdklahn). (#33)
- Improved usage examples, prompted by @Ackia. (#34)



## [1.0.0] - 2024-01-22
### Added
- Automatically download and use any standard index by name. Example: `hostile clean --index human-t2t-hla-argos985 --fastq1 reads.fq.gz`. (#28)
- Two new standard indexes:
  - `human-t2t-hla.rs-viral-202401_ml-phage-202401` (RefSeq viral & Millard Lab phage)
  - `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` (ARGOS bacteria, RefSeq viral, & Millard Lab phage)
- Improved masking protocol for custom index creation (`hostile mask`):
  - Uses `minimap2` with adjusted secondary alignment limits and score thresholds.
  - K-merises target genomes into 150mers (10bp step) before alignment. Adds `dnaio` dependency.
- New options for `hostile clean`: `--aligner-args`, `--invert`, and `--offline`.
- Verification of downloaded index checksums in `manifest.json`. (#20)
- Adds `version`, `aligner`, and `options` fields to log output.

### Changed
- Prevents corruption of empty `fastq.gz` files when no reads remain after decontamination. (#24)



## [0.4.0] - 2023-11-24

### Added
- `hostile fetch` subcommand for viewing and downloading prebuilt references.
- Automatic determination of near-ideal alignment thread count based on CPU availability.

### Changed
- Renamed `--sort-by-name` to `--reorder`, providing deterministic read order without necessarily sorting.
- Improved decontamination speed for Linux by using Bowtie2’s `--reorder` flag.



## [0.3.0] - 2023-11-22

### Added
- `--aligner-args` option for arbitrary alignment parameter passthrough (e.g., `--aligner-args="--ignore-quals"`).
- `--sort-by-name` option for deterministic output order when using Bowtie2 with multiple threads.



## [0.2.0] - 2023-11-10

### Changed
- Performance improvements for Illumina reads / Bowtie.
- Added a space before `/1` and `/2` when renaming Illumina reads (with `--rename`), enhancing compatibility.
- Improved logging.



## [0.1.0] - 2023-07-23

### Changed
- Addressed issues #10, #19, #20, #22, #23, and #25.
- Added extensive new tests.
