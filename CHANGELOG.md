# Changelog

## [1.2.0](https://github.com/RIVM-bioinformatics/juno-snp/compare/v1.1.2...v1.2.0) (2025-01-20)


### Features

* allow masking with custom ref ([d436eae](https://github.com/RIVM-bioinformatics/juno-snp/commit/d436eae096c48d3aca0bd973390fa498b27a1cca))


### Bug Fixes

* updated dependency issues in main environment ([7971c7a](https://github.com/RIVM-bioinformatics/juno-snp/commit/7971c7ad1d95a3de273470e65e44c7b2fe20ff4a))

## [1.1.2](https://github.com/RIVM-bioinformatics/juno-snp/compare/v1.1.1...v1.1.2) (2024-04-26)


### Bug Fixes

* Fix type error doing dryrun with fixed reference ([8539c88](https://github.com/RIVM-bioinformatics/juno-snp/commit/8539c88be5e182958ea70e3634fb2328e66fb98d))

## [1.1.1](https://github.com/RIVM-bioinformatics/juno-snp/compare/v1.1.0...v1.1.1) (2024-02-23)


### Bug Fixes

* add hpc resources ([b15eaed](https://github.com/RIVM-bioinformatics/juno-snp/commit/b15eaed3302a660e300d55cd1dc86f8fc20a297e))
* convert hpc config values to int ([55e5656](https://github.com/RIVM-bioinformatics/juno-snp/commit/55e56569029942fdbecfc718a7881f7c40cedff3))
* improve error message for outlier genomes ([2c038b4](https://github.com/RIVM-bioinformatics/juno-snp/commit/2c038b4bd1eec86b13266a0e2b001c570c5adcfc))

## [1.1.0](https://github.com/RIVM-bioinformatics/juno-snp/compare/v1.0.0...v1.1.0) (2024-02-23)


### Features

* add QC per cluster ([65960ef](https://github.com/RIVM-bioinformatics/juno-snp/commit/65960efa3a8d5dd8eb43c5618f61eb636f9758f3))
* allow toggling snippy --report flag ([0949a61](https://github.com/RIVM-bioinformatics/juno-snp/commit/0949a61b522b29a89d33f4be9805ccbc31ccfcae))


### Bug Fixes

* allow gbk and fasta refs ([7dbea61](https://github.com/RIVM-bioinformatics/juno-snp/commit/7dbea6137ddbeabc1d33417592c259b47ab72510))

## [1.0.0](https://github.com/RIVM-bioinformatics/juno-snp/compare/v0.0.3...v1.0.0) (2024-02-21)


### ⚠ BREAKING CHANGES

* update juno library to v2

### Code Refactoring

* update juno library to v2 ([347971f](https://github.com/RIVM-bioinformatics/juno-snp/commit/347971f22d28726b274b296b552128f934320a40))

## [0.0.3](https://github.com/RIVM-bioinformatics/juno-snp/compare/v0.0.2...v0.0.3) (2023-12-05)


### Bug Fixes

* rm unneeded rules ([5df54ff](https://github.com/RIVM-bioinformatics/juno-snp/commit/5df54ff9307d561f3ecbd267eeca69425ab1531c))
* set prefix to str ([89ea436](https://github.com/RIVM-bioinformatics/juno-snp/commit/89ea4361b912ea6d419e3944b151a3e74952dabb))


### Dependencies

* re-add containers ([3930b55](https://github.com/RIVM-bioinformatics/juno-snp/commit/3930b5539edd4824d956770faae94b75236345c3))


### Documentation

* update README ([6dcb6c6](https://github.com/RIVM-bioinformatics/juno-snp/commit/6dcb6c60b7235c041ced7803935e328bc854cf5e))

## [0.0.2](https://github.com/RIVM-bioinformatics/juno-snp/compare/v0.0.1...v0.0.2) (2023-07-24)


### Bug Fixes

* fix environments ([7e2857e](https://github.com/RIVM-bioinformatics/juno-snp/commit/7e2857e2b17e61c276841c79bcd95cbe02aacb28))


### Dependencies

* clean up envs ([4eab33e](https://github.com/RIVM-bioinformatics/juno-snp/commit/4eab33e53c19bfd4cd49b6d118df84ea660c8fbe))
* fix dependencies for multiple environments/tools ([3f8e2e6](https://github.com/RIVM-bioinformatics/juno-snp/commit/3f8e2e641d7ae812dd98338a8d67ef64ea474156))
* remove anaconda and defaults channel and add nodefaults ([84d6132](https://github.com/RIVM-bioinformatics/juno-snp/commit/84d6132fbf8084f4fefc2f28181d834c5204c975))
