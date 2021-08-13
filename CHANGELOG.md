# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2021-08-13
Pulling changes from updates added during BLAST analysis

### Added

  - ability to manually remove samples
  - custom version of janitor::make_clean_names that optionally allows duplicate values

### Changed

  -Replaced some hard coded variables with abstraction.
    - Replaced comparing by disease_class with a variable `comparison_grouping_variable`
  - Palette generation is a little smarter
  - Module data is now split by a custom function (helps with NSE)
  - Update C5 MSigDb file

### Removed

  - Eliminated the excessive number of times I was mutating that comparison
    column into a factor.


## [1.1.0] - 2021-04-22
### Changed
- Switched to using a generic metadata template instead of trying to customize
  the import_metadata function for each new dataset.

### Added
- Metadata template


## [1.0.0] - 2021-03-08
### Added
- Started project
- CHANGELOG.md
- README.md

### Changed
- Rearranged directory layout
- Split analysis plan into parts

[1.1.0]: https://github.com/olivierlacan/keep-a-changelog/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/1.0.0
