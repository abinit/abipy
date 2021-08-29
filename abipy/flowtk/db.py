"""
"""
from pydantic import BaseModel, Field


class MongoConnector(BaseModel):

    host: str = Field(..., description="Host address e.g. 0.0.0.0")

    port: int = Field(27017, description="MongoDB server port. Default: 27017")

    db_name: str = Field("abipy", description="Name of the MongoDB database")

    collection: str = Field(..., description="Name of the collection")

    user: str = Field(None, description="User name. Default: None")

    password: str = Field(None, description="password for authentication. Default: None")

    def get_client(self):
        """
        Establish a connection with the database. Return MongoClient
        """
        from pymongo import MongoClient
        if self.host and self.port:
            client = MongoClient(host=self.host, port=self.port)
        else:
            client = MongoClient()

        return client

    def get_collection(self):
        """
        Returns MongoDb collection
        """
        client = self.get_client()
        print(client)

        db = client[self.db_name]
        # Authenticate if needed
        if self.user and self.password:
            db.autenticate(self.user, password=self.password)

        return db[self.collection]


if __name__ == "__main__":
    connector = MongoConnector(host="localhost", collection="foo")
    print(connector)
    col = connector.get_collection()
    print(col)
    r = col.insert_one({"foo": 1})
    print(r.inserted_id)