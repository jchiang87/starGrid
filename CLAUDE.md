# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Run tests:**
```bash
python -m unittest tests.test_registration
```

**Run a single test:**
```bash
python -m unittest tests.test_registration.StarGridRegistrationTestCase.test_skycatalog_yaml_loading
```

**Install (EUPS-based environment):**
```bash
setup star_grid
```

## Architecture

`starGrid` is a plugin for the [skyCatalogs](https://github.com/LSSTDESC/skyCatalogs) framework that generates synthetic star catalogs arranged in regular RA/Dec grids. It is used for astronomical image simulation (e.g., with imsim).

### Plugin Integration

skyCatalogs discovers object types via YAML config files. Each object type entry must specify `module: star_grid` and `collection_class: StarGridCollection`. When `skyCatalogs.open_catalog()` is called, it invokes `register_objects()` from the `star_grid` module for each matching entry, which calls `StarGridCollection.register()` to register the type with skyCatalogs' internal `cat_cxt`.

When a region query is made, skyCatalogs calls `StarGridCollection.load_collection()`, which reads `num_stars`, `sed_path`, and `magnorm` from the YAML config and instantiates a `StarGridCollection`.

### Key Classes

- **`StarGridCollection`** (`star_grid/starGrid.py`): Generates a meshgrid of RA/Dec points to cover the queried region. Grid dimensions are calculated so that `nra * ndec ≈ num_stars` with spacing adjusted for declination via `cos(dec)`. All stars share a single SED (loaded from `sed_path` and normalized to `magnorm`).

- **`StarGridObject`** (`star_grid/starGrid.py`): Represents one star in the grid. Returns a `galsim.DeltaFunction` as its GSObject and the collection's shared SED. The only valid SED component name is `"this_object"`.

### Configuration

Object types are named arbitrarily (e.g., `star_grid_100`) and can co-exist in one YAML file. The `num_stars` parameter controls grid density — actual star count may differ slightly because grid dimensions are rounded up via `np.ceil`. Object IDs are formatted as `{object_type}_{index}`.

The test in `tests/test_registration.py` loads `examples/skyCatalog.yaml`, which has hardcoded absolute paths for `sed_path` pointing to a local skyCatalogs installation.
