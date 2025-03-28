from .utils import object_type_config
from .starGrid import StarGridCollection

COLLECTION_CLASS_MAP = {"StarGridCollection": StarGridCollection}

def register_objects(sky_catalog, object_type):
    config = object_type_config(sky_catalog, object_type)
    if "collection_class" not in config:
        return
    collection_class = COLLECTION_CLASS_MAP[config["collection_class"]]
    collection_class.register(sky_catalog, object_type)



