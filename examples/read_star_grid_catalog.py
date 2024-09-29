from skycatalogs import skyCatalogs
from skycatalogs.utils import Disk

yaml_file = "./skyCatalog.yaml"
skycatalog_root = "."

sky_cat = skyCatalogs.open_catalog(yaml_file, skycatalog_root=skycatalog_root)

ra0, dec0 = 9.5, -44.0
radius = 0.17*3600.0
region = Disk(ra0, dec0, radius)
objs = sky_cat.get_objects_by_region(region, obj_type_set={'star_grid'})

print(len(objs))
print(objs[0:20])
obj = objs[0]
print(obj.get_LSST_flux('r'))
