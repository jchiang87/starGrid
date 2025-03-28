"""
Utilities for interacting with skyCatalogs code.
"""
__all__ = ["object_type_config"]


def object_type_config(sky_catalog, object_type):
    return dict(sky_catalog.raw_config["object_types"][object_type])
