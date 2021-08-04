# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2021-08-04
### Added
- ability to manually remove samples
- custom version of janitor::make_clean_names that optionally allows duplicate values

### Changed
- Replaced some hard coded variables with abstraction.
  - Replaced comparing by disease_class with a variable `comparison_grouping_variable`
- Palette generation is a little smarter
- Module data is now split by a custom function (helps with NSE)
- Update C5 MSigDb file

### Removed
- Eliminated the excessive number of times I was mutating that comparison column
  into a factor.

## [1.0.0] - 2021-03-08
### Added
- Started project
- CHANGELOG.md
- README.md

### Changed
- Rearranged directory layout
- Split analysis plan into parts
