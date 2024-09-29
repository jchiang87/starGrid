import numpy as np
import galsim
from skycatalogs.objects import BaseObject, ObjectCollection
from skycatalogs.utils import normalize_sed


__all__ = ["StarGridCollection", "StarGridObject"]


class StarGridObject(BaseObject):

    def get_observer_sed_component(self, component, mjd=None):
        if component != "this_object":
            raise RuntimeError("Unknown SED component: %s", component)
        return self.belongs_to.sed

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}


class StarGridCollection(ObjectCollection):
    """
    Arrange stars in a grid in RA, Dec to cover the specified region
    of the sky.
    """
    _object_type = "star_grid"
    def __init__(self, region, sky_catalog, num_stars, sed_path, magnorm,
                 obj_id_offset=0, **kwds):
        # Create a grid of stars that cover the region.
        ra, dec = self._create_star_grid(region, num_stars)

        # Compute the SED that all StarGridObjects will use.
        lut = galsim.LookupTable.from_file(sed_path, interpolant='linear')
        self.sed = normalize_sed(
            galsim.SED(lut, wave_type='angstrom', flux_type='flambda'),
            magnorm)

        # Fill the private attributes required by the base class.
        self._ra = np.ravel(ra)
        self._dec = np.ravel(dec)
        self._id = [str(_) for _ in np.arange(len(self)) + obj_id_offset]
        self._sky_catalog = sky_catalog
        self._object_type_unique = self._object_type
        self._object_class = StarGridObject
        self._uniform_object_type = True

    def _create_star_grid(self, region, num_stars):
        # Create a grid of stars covering the specified region.
        ra_min, ra_max, dec_min, dec_max = region.get_radec_bounds()
        ra0 = (ra_max + ra_min)/2.0
        dec0 = (dec_max + dec_min)/2.0

        # Set the number of grid points in ra and dec so that
        # nra*ndec~num_stars and the spacings in ra and dec are
        # approximately the same.
        ratio = np.cos(np.radians(dec0))
        ndec = int(np.ceil(np.sqrt(num_stars/ratio)))
        nra = int(np.ceil(np.sqrt(num_stars*ratio)))

        # Compute the grid points.
        ra_vals = np.linspace(ra_min, ra_max, nra)
        dec_vals = np.linspace(dec_min, dec_max, ndec)
        return np.meshgrid(ra_vals, dec_vals)

    @property
    def native_columns(self):
        return ()

    def __len__(self):
        return len(self._ra)

    @staticmethod
    def register(sky_catalog):
        sky_catalog.cat_cxt\
                   .register_source_type(StarGridCollection._object_type,
                                         object_class=StarGridObject,
                                         collection_class=StarGridCollection,
                                         custom_load=True)

    @staticmethod
    def load_collection(region, sky_catalog, mjd=None, exposure=None):
        # Get catalog parameters from config file.
        config = dict(sky_catalog.raw_config["object_types"]\
                      [StarGridCollection._object_type])
        return StarGridCollection(region, sky_catalog, **config)
