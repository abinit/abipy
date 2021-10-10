"""
Global configuration options for Abipy stored in the config.yml file
"""
from __future__ import annotations

import os

from typing import Optional, List
from pydantic import BaseSettings, Field
from monty.json import jsanitize
from abipy.tools.iotools import yaml_safe_load


_CONFIG = None


def get_config() -> AbipyConfig:
    """
    Build and return a singleton with the configuration options.
    This is the main entry point for client code.
    """
    global _CONFIG
    if _CONFIG is not None: return _CONFIG
    _CONFIG = AbipyConfig.from_yaml_file()
    return _CONFIG


class AbipyConfig(BaseSettings):
    """
    Pydantic model with the global configuration options.
    Since we inherit from BaseSettings, it is possible to use env variables to specify the values.
    See https://pydantic-docs.helpmanual.io/usage/settings/
    """

    mongo_username: str = Field(None, description="User name for MongDB authentication. Implies password.")

    mongo_password: str = Field(None, description="Password for MongDB authentication")

    mongo_host: str = Field(None, description="Host address e.g. 'localhost'")

    mongo_port: int = Field(27017, description="MongoDB server port e.g. 27017")

    mongo_dbname: str = Field(None, description="MongoDB database")

    worker_scratchdir: str = Field(None, description="Scratch directory used by AbipyWorkers to generate Flows")

    remote_hosts_for_workers: List[str] = Field(None, description="List of remote hosts used to run AbipyWorkers")

    class Config:
        case_sensitive = False
        #env_prefix = 'abipy_'  # defaults to no prefix, i.e. ""


    @classmethod
    def from_yaml_file(cls, filepath: Optional[str] = None) -> AbipyConfig:
        """
        Read the configuration options from a Yaml file.
        Use default configuration file in ~/.abinit/abipy if filepath is None.
        """
        if filepath is None:
            filepath = os.path.expanduser(os.path.join("~", ".abinit", "abipy", "config.yml"))

        if not os.path.exists(filepath):
            # Create file and fill it with commented schema in json format.
            s = cls.schema_json(indent=2)
            lines = [
                "# This is the JSON schema with the variables that can be used in the Abipy configuration file.",
                "# Use the schema as a reference but keep in mind that JSON can be replaced by Yaml syntax.",
                "# For instance, one can use:",
                "#",
                "# mongo_username: John",
                "# mongo_port: 1234",
                "#",
                "#",
            ]
            lines += [f"#{t}" for t in  s.split("\n")]
            s = "\n".join(lines)

            with open(filepath, "wt") as fh:
                fh.write(s)

            return cls()

        with open(filepath, "rt") as fh:
            return cls.from_yaml_string(fh.read())

    @classmethod
    def from_yaml_string(cls, yaml_string: str) -> AbipyConfig:
        """
        Read the configuration options from a Yaml file.
        Supports shell variables e.g.

            home: $HOME
        """
        data = yaml_safe_load(yaml_string)
        # From Yaml to JSON.
        data = jsanitize(data)
        unknown_keys = {k for k in data if k not in cls.__fields__}
        if unknown_keys:
            raise ValueError(f"Found unknown keys in AbiPy configuration file:\n{unknown_keys}")

        # Expand shell variables
        for k, v in data.items():
            if isinstance(v, str) and v.startswith("$"):
                data[k] = os.environ[v[1:].strip()]

        return cls(**data)



if __name__ == "__main__":
    from pprint import pprint
    print(AbipyConfig.schema_json(indent=2))
    #print(AbipyConfig.__fields__)
    #config = get_config()
    #print(config)
