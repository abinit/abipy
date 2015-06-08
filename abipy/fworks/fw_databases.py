# coding: utf-8
"""
Utilities for database insertion
"""

from pymatgen.serializers.json_coders import PMGSONable
import gridfs
import json
import pymongo


class MongoDatabase(PMGSONable):
    """
    MongoDB database class for access, insertion, update, ... in a MongoDB database
    """

    def __init__(self, host, port, database, username, password, collection, gridfs_collection=None):
        self._host = host
        self._port = port
        self._database = database
        self._username = username
        self._password = password
        self._collection = collection
        self._gridfs_collection = gridfs_collection
        self._connect()

    def _connect(self):
        self.server = pymongo.MongoClient(host=self._host, port=self._port)
        self.database = self.server[self._database]
        self.database.authenticate(name=self._username, password=self._password)
        self.collection = self.database[self._collection]
        if self._gridfs_collection is not None:
            self.gridfs = gridfs.GridFS(self.database, collection=self._gridfs_collection)
        else:
            self.gridfs = None

    def insert_entry(self, entry):
        self.collection.insert(entry)

    def update_entry(self, query, entry_update):
        count = self.collection.find(query).count()
        if count != 1:
            raise RuntimeError("Number of entries != 1, found : {:d}".format(count))
        entry = self.collection.find_one(query)
        entry.update(entry_update)
        self.collection.save(entry)

    def as_dict(self):
        """
        Json-serializable dict representation of a MongoDatabase
        """
        dd = {"@module": self.__class__.__module__,
              "@class": self.__class__.__name__,
              "host": self._host,
              "port": self._port,
              "database": self._database,
              "username": self._username,
              "password": self._password,
              "collection": self._collection,
              "gridfs_collection": self._gridfs_collection}
        return dd

    @classmethod
    def from_dict(cls, d):
        return cls(host=d['host'], port=d['port'], database=d['database'],
                   username=d['username'], password=d['password'], collection=d['collection'],
                   gridfs_collection=d['gridfs_collection'])