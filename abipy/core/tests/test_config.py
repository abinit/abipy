"""
Unit tests for config.py module
"""
import os

from abipy.core.testing import AbipyTest
from abipy.core.config import AbipyConfig, get_config


class TestAbipyConfig(AbipyTest):

    def test_base(self):

        assert get_config()

        schema = AbipyConfig.schema()
        assert schema
        #print(schema)

        yaml_string = """
mongo_host: mongo.host.com
mongo_port: 27017
mongo_dbname: my_database
mongo_username: john_doe
mongo_password: secret
worker_scratchdir: $TMPDIR
#hello: 1
"""

        config = AbipyConfig.from_yaml_string(yaml_string)
        assert config.mongo_port == 27017
        assert config.mongo_password == "secret"
        assert config.worker_scratchdir == os.environ["TMPDIR"]

        bad_string = yaml_string.replace("_port", "_prt")
        with self.assertRaises(ValueError):
            AbipyConfig.from_yaml_string(bad_string)

        bad_string = yaml_string.replace("27017", "27017,")
        with self.assertRaises(ValueError):
            AbipyConfig.from_yaml_string(bad_string)

        # Read official configuration file if present.
        config_path = os.path.expanduser(os.path.join("~", ".abinit", "abipy", "config.yml"))
        if os.path.exists(config_path):
            AbipyConfig.from_yaml_file(config_path)

        # This to test the case in which the file does not exist
        # and we have to create it with the json schema
        tmp_path = self.get_tmpname()
        #print("tmp_path", tmp_path)
        os.remove(tmp_path)
        AbipyConfig.from_yaml_file(tmp_path)
        #assert 0


