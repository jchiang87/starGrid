"""Unit tests for registering with skyCatalogs"""

import os
from pathlib import Path
import unittest
from skycatalogs import skyCatalogs
from skycatalogs.utils import Disk


PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent)))


class StarGridRegistrationTestCase(unittest.TestCase):
    """
    TestCase class for StarGrid registration.
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_skycatalog_yaml_loading(self):
        """Test loading of object types via skyCatalogs yaml config file."""
        yaml_file = os.path.join(PACKAGE_DIR, "examples", "skyCatalog.yaml")
        skycatalog_root = os.path.dirname(yaml_file)
        skycat = skyCatalogs.open_catalog(yaml_file,
                                          skycatalog_root=skycatalog_root)
        region = Disk(0, 0, 3600.0)
        object_params = {"star_grid_100": 100, "star_grid_25": 25}
        for object_type, num_stars in object_params.items():
            # Check that object types are registered by name.
            self.assertIn(object_type, skycat.raw_config["object_types"])
            self.assertIn(object_type, skycat.cat_cxt._source_type_dict)
            # Check that the object type parameters are set in the collection class.
            obj = skycat.get_object_type_by_region(region, object_type)[0]
            self.assertEqual(obj.belongs_to.num_stars, num_stars)


if __name__ == "__main__":
    unittest.main()

