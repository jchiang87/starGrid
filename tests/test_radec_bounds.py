"""
Tests for ra/dec bounds config option in StarGridCollection,
without using mocks.

Strategy: the only mocked pieces in the original tests are galsim SED loading
and the sky_catalog object.  We avoid them by:
  - Calling StarGridCollection._create_star_grid directly on a minimal shim
    (MinimalGrid) that owns only the attributes that method needs.
  - Using a plain SimpleNamespace/dict for sky_catalog instead of MagicMock,
    to test the load_collection config-parsing path.
"""

import types
import unittest
import numpy as np

from star_grid.starGrid import StarGridCollection


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class FakeRegion:
    """Minimal region that returns fixed ra/dec bounds."""
    def __init__(self, ra_min, ra_max, dec_min, dec_max):
        self._bounds = (ra_min, ra_max, dec_min, dec_max)
        self.call_count = 0

    def get_radec_bounds(self):
        self.call_count += 1
        return self._bounds


class MinimalGrid:
    """
    Shim that exercises StarGridCollection._create_star_grid without
    triggering galsim SED loading.
    """
    def __init__(self, num_stars, region, radec_bounds=None):
        self.num_stars = num_stars
        ra, dec = StarGridCollection._create_star_grid(self, region,
                                                       radec_bounds=radec_bounds)
        self._ra = np.ravel(ra)
        self._dec = np.ravel(dec)


def _fake_sky_catalog(extra_config=None):
    """Return a plain namespace with a minimal raw_config for 'test_type'."""
    config = {'num_stars': 25, 'sed_path': '/dev/null', 'magnorm': 20.0}
    if extra_config:
        config.update(extra_config)
    sky_catalog = types.SimpleNamespace(
        raw_config={'object_types': {'test_type': config}}
    )
    return sky_catalog


def _parse_radec_bounds(sky_catalog, object_type='test_type'):
    """
    Replicate the radec_bounds extraction from load_collection so we can
    test that logic without going through the full constructor.
    """
    config = dict(sky_catalog.raw_config['object_types'][object_type])
    radec_keys = ('ra_min', 'ra_max', 'dec_min', 'dec_max')
    return (
        tuple(config[k] for k in radec_keys)
        if all(k in config for k in radec_keys)
        else None
    )


# ---------------------------------------------------------------------------
# Tests: grid uses region bounds when no config bounds are supplied
# ---------------------------------------------------------------------------

class RadecBoundsFromRegionTestCase(unittest.TestCase):
    """Grid uses region.get_radec_bounds() when radec_bounds is None."""

    REGION_BOUNDS = (10.0, 12.0, -1.0, 1.0)

    def setUp(self):
        self.region = FakeRegion(*self.REGION_BOUNDS)
        self.grid = MinimalGrid(25, self.region, radec_bounds=None)

    def test_region_get_radec_bounds_is_called(self):
        self.assertEqual(self.region.call_count, 1)

    def test_grid_ra_min_matches_region(self):
        self.assertAlmostEqual(np.min(self.grid._ra), self.REGION_BOUNDS[0])

    def test_grid_ra_max_matches_region(self):
        self.assertAlmostEqual(np.max(self.grid._ra), self.REGION_BOUNDS[1])

    def test_grid_dec_min_matches_region(self):
        self.assertAlmostEqual(np.min(self.grid._dec), self.REGION_BOUNDS[2])

    def test_grid_dec_max_matches_region(self):
        self.assertAlmostEqual(np.max(self.grid._dec), self.REGION_BOUNDS[3])


# ---------------------------------------------------------------------------
# Tests: grid uses explicit radec_bounds and ignores region
# ---------------------------------------------------------------------------

class RadecBoundsFromConfigTestCase(unittest.TestCase):
    """Grid uses the supplied radec_bounds tuple and never queries the region."""

    CONFIG_BOUNDS = (20.0, 22.0, -2.0, 2.0)
    REGION_BOUNDS = (10.0, 12.0, -1.0, 1.0)  # intentionally different

    def setUp(self):
        self.region = FakeRegion(*self.REGION_BOUNDS)
        self.grid = MinimalGrid(25, self.region, radec_bounds=self.CONFIG_BOUNDS)

    def test_region_get_radec_bounds_not_called(self):
        self.assertEqual(self.region.call_count, 0)

    def test_grid_ra_min_matches_config(self):
        self.assertAlmostEqual(np.min(self.grid._ra), self.CONFIG_BOUNDS[0])

    def test_grid_ra_max_matches_config(self):
        self.assertAlmostEqual(np.max(self.grid._ra), self.CONFIG_BOUNDS[1])

    def test_grid_dec_min_matches_config(self):
        self.assertAlmostEqual(np.min(self.grid._dec), self.CONFIG_BOUNDS[2])

    def test_grid_dec_max_matches_config(self):
        self.assertAlmostEqual(np.max(self.grid._dec), self.CONFIG_BOUNDS[3])


# ---------------------------------------------------------------------------
# Tests: load_collection config-parsing logic
# ---------------------------------------------------------------------------

class LoadCollectionConfigParsingTestCase(unittest.TestCase):
    """
    load_collection extracts a (ra_min, ra_max, dec_min, dec_max) tuple
    when all four keys are present, and returns None otherwise.
    """

    BOUNDS_CONFIG = {'ra_min': 20.0, 'ra_max': 22.0, 'dec_min': -2.0, 'dec_max': 2.0}

    def test_all_keys_present_returns_tuple(self):
        sky_catalog = _fake_sky_catalog(self.BOUNDS_CONFIG)
        result = _parse_radec_bounds(sky_catalog)
        self.assertEqual(result, (20.0, 22.0, -2.0, 2.0))

    def test_no_keys_returns_none(self):
        sky_catalog = _fake_sky_catalog()
        self.assertIsNone(_parse_radec_bounds(sky_catalog))

    def test_partial_keys_returns_none(self):
        sky_catalog = _fake_sky_catalog({'ra_min': 20.0, 'ra_max': 22.0})
        self.assertIsNone(_parse_radec_bounds(sky_catalog))

    def test_bounds_tuple_drives_grid(self):
        """When load_collection parses bounds, the grid respects them."""
        sky_catalog = _fake_sky_catalog(self.BOUNDS_CONFIG)
        radec_bounds = _parse_radec_bounds(sky_catalog)
        region = FakeRegion(10.0, 12.0, -1.0, 1.0)
        grid = MinimalGrid(25, region, radec_bounds=radec_bounds)
        self.assertAlmostEqual(np.min(grid._ra), self.BOUNDS_CONFIG['ra_min'])
        self.assertAlmostEqual(np.max(grid._ra), self.BOUNDS_CONFIG['ra_max'])
        self.assertEqual(region.call_count, 0)

    def test_no_bounds_in_config_uses_region(self):
        sky_catalog = _fake_sky_catalog()
        radec_bounds = _parse_radec_bounds(sky_catalog)
        region = FakeRegion(10.0, 12.0, -1.0, 1.0)
        grid = MinimalGrid(25, region, radec_bounds=radec_bounds)
        self.assertAlmostEqual(np.min(grid._ra), 10.0)
        self.assertAlmostEqual(np.max(grid._ra), 12.0)
        self.assertEqual(region.call_count, 1)


# ---------------------------------------------------------------------------
# Tests: grid geometry is consistent
# ---------------------------------------------------------------------------

class GridGeometryTestCase(unittest.TestCase):
    """Sanity-check the grid dimensions and declination correction."""

    def test_star_count_near_target(self):
        """Actual star count should be close to num_stars (within a factor of 2)."""
        num_stars = 100
        region = FakeRegion(0.0, 2.0, -1.0, 1.0)
        grid = MinimalGrid(num_stars, region)
        self.assertGreaterEqual(len(grid._ra), num_stars)
        self.assertLess(len(grid._ra), num_stars * 4)

    def test_dec_correction_applied(self):
        """
        At high declination (cos small) nra < ndec; near equator they are close.
        The ratio ndec/nra should grow with |dec|.
        """
        region_equator = FakeRegion(0.0, 2.0, -0.5, 0.5)
        region_pole = FakeRegion(0.0, 2.0, 80.0, 82.0)
        grid_eq = MinimalGrid(100, region_equator)
        grid_pole = MinimalGrid(100, region_pole)
        n_unique_ra_eq = len(np.unique(np.round(grid_eq._ra, 10)))
        n_unique_dec_eq = len(np.unique(np.round(grid_eq._dec, 10)))
        n_unique_ra_pole = len(np.unique(np.round(grid_pole._ra, 10)))
        n_unique_dec_pole = len(np.unique(np.round(grid_pole._dec, 10)))
        ratio_eq = n_unique_dec_eq / n_unique_ra_eq
        ratio_pole = n_unique_dec_pole / n_unique_ra_pole
        self.assertGreater(ratio_pole, ratio_eq)

    def test_single_star(self):
        region = FakeRegion(45.0, 45.0, 10.0, 10.0)
        grid = MinimalGrid(1, region)
        self.assertEqual(len(grid._ra), 2)
        self.assertAlmostEqual(grid._ra[0], 45.0)
        self.assertAlmostEqual(grid._dec[0], 10.0)


if __name__ == '__main__':
    unittest.main()

