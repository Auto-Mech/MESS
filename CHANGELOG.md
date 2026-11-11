# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions are `YYYY.M.N`, where `N` counts the releases within the month.

## [Unreleased]

## [2026.11.0] - 2026-11-11

### Added

- The `mess` conda package: the executables linked against conda-forge's ILP64 OpenBLAS and GSL.
- Releases: a version tag builds and tests both packages, publishes them, and makes a GitHub release with the packages and a tarball of the static executables.

### Changed

- `mess-static` is built from source with CMake, instead of repackaging prebuilt executables. The executables are `mess`, `messpf` and `messabs`, and their output for the examples is unchanged.
- `mess` and `mess-static` share one version.
- The packages are published to the [`automech`](https://prefix.dev/channels/automech) channel on prefix.dev. The `auto-mech` channel on anaconda.org keeps the earlier `mess-static` packages and gets no new ones.
- `examples/mess.inp` sets `WithCRM`, which its radiative oscillators need.

### Removed

- The `gumbo` executable from `mess-static`: its example input (`examples/messsym.inp`) does not run, so it cannot be tested.

## [2026.2.4] - 2026-02-04

- `mess-static` on anaconda.org, with executables built from the source of that date.

## [2025.9.0] - 2025-09-25

- `mess-static` on anaconda.org, packaged with Pixi.

## [2024.5.16] - 2024-05-16

- `mess-static` on anaconda.org, the first conda package of the static executables.
