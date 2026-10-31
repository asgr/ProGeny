# AGENTS.md — ProGeny

R package (LGPL-3) to generate simple stellar population (SSP) libraries by combining IMFs, isochrones and spectral atmospheres.

## Layout
- `R/` source. `IMF.R`/`IMF_evo.R` (IMFs), `utility.R` (bin/interp helpers), `progeny*.R` (isochrone/track/SSP pipeline), `shinyProGeny.R` (app), `zzz.R` (load hooks).
- `man/` Rd files (hand-maintained; update when changing exported signatures). `vignettes/`, `data/`, `inst/`.
- `tests/testthat/` testthat (edition 3) unit tests; entry point `tests/testthat.R`.
- `NAMESPACE` is hand-written: `exportPattern("^[[:alpha:]]+")` exports everything not starting with `.`. Do not run roxygen to overwrite it. Dot-prefixed functions are internal; test them with `ProGeny:::`.

## Commands
- Tests: `Rscript -e 'devtools::test()'` (or `testthat::test_local()`).
- Check: `R CMD build . && R CMD check --no-manual --as-cran ProGeny_*.tar.gz`.
- Coverage: `Rscript -e 'covr::package_coverage()'`.
- Dependencies Rfits, ParmOff, ProPane, ProSpect, Rwcs come from GitHub (`asgr/...`).

## Testing guidelines
- Tests must be fast, deterministic and offline: no downloads (`progenyDownload`), no Shiny, no large data fetches, no plotting.
- Prefer testing numerical invariants (normalisation, continuity, symmetry, monotonic bins) over snapshots.
- Heavier pipeline functions should be tested with tiny synthetic isochrone/spectra inputs; use `skip_on_cran()`/`skip_if_not_installed()` for optional (Suggests) packages.

## Conventions
- Code uses `=` for assignment and base-R style; keep style consistent with surrounding code.
- Bump `Version`/`Date` in `DESCRIPTION` for user-facing changes.
- CI: `.github/workflows/main.yml` runs R CMD check (which runs the tests) on several OSes; `test-coverage.yml` reports coverage.
