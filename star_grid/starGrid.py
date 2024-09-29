import numpy as np
import galsim
from skycatalogs.objects import BaseObject, ObjectCollection
from skycatalogs.utils import normalize_sed


__all__ = ["StarGridCollection", "StarGridObject"]


class StarGridObject(BaseObject):

    def __init__(self, ra, dec, obj_id, parent_collection, index):
        super().__init__(ra, dec, obj_id, "star_grid",
                         belongs_to=parent_collection, belongs_index=index)

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
    def __init__(self, region, num_stars, magnorm, sed_path, obj_id_offset,
                 sky_catalog):
        self._sky_catalog = sky_catalog
        self._object_type_unique = "star_grid"
        # Create grid of stars covering the specified region.
        ra_min, ra_max, dec_min, dec_max = region.get_radec_bounds()
        ra0 = (ra_max + ra_min)/2.0
        dec0 = (dec_max + dec_min)/2.0

        # Compute the number of grid points in ra and dec so
        # that nra*ndec ~ num_stars and the spacing in ra and dec
        # are approximately the same.
        ratio = np.cos(np.radians(dec0))
        ndec = int(np.ceil(np.sqrt(num_stars/ratio)))
        nra = int(np.ceil(np.sqrt(num_stars*ratio)))

        # Compute grid points.
        ra_vals = np.linspace(ra_min, ra_max, nra)
        dec_vals = np.linspace(dec_min, dec_max, ndec)
        ra, dec = np.meshgrid(ra_vals, dec_vals)
        self.ra = np.ravel(ra)
        self.dec = np.ravel(dec)
        self.index = np.arange(len(self.ra))
        self.obj_id = [str(_) for _ in self.index + obj_id_offset]
        lut = galsim.LookupTable.from_file(sed_path, interpolant='linear')
        self.sed = normalize_sed(
            galsim.SED(lut, wave_type='angstrom', flux_type='flambda'),
            magnorm)

    @property
    def native_columns(self):
        return ()

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            return StarGridObject(self.ra[key], self.dec[key],
                                  self.obj_id[key], self, key)
        elif isinstance(key, slice):
            return [self.__getitem__(i) for i in self.index[key]]
        raise TypeError(f"Index must be integers or slices, {type(key)}")

    def __len__(self):
        return len(self.index)

    @staticmethod
    def register(sky_catalog):
        sky_catalog.cat_cxt\
                   .register_source_type('star_grid',
                                         object_class=StarGridObject,
                                         collection_class=StarGridCollection,
                                         custom_load=True)

    @staticmethod
    def load_collection(region, skycatalog, mjd=None, exposure=None):
        # Get catalog parameters from config file.  A more general
        # way of obtaining parameters would be nice.
        config = skycatalog.raw_config["object_types"]["star_grid"]
        num_stars = config["num_stars"]
        magnorm = config["magnorm"]
        sed_path = config["sed_path"]
        obj_id_offset = config["obj_id_offset"]

        return StarGridCollection(region, num_stars, magnorm, sed_path,
                                  obj_id_offset, skycatalog)
