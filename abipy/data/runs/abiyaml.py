#!/usr/bin/env python
from __future__ import division, print_function

import collections
#import os
#import numpy as np
import yaml

class YamlDoc(object):
    """A string with a YAML tag."""

    def __init__(self, string):
        self.string = string 
        self.tag = None

    #def __repr__(self):
    #    return self.string

    def __str__(self):
        return self.string


class YamlFileReaderError(Exception):
    """Exception raised by `YamlFileReader`."""


class YamlFileReader(collections.Iterator):
    """
    A file locking mechanism that has context-manager support 
    so you can use it in a with statement. 
    """
    Error = YamlFileReaderError

    def __init__(self, filename):
        self.stream = open(filename, "r")

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        self.stream.close()

    def seek(self, offset, whence=0):
        """
        seek(offset[, whence]) -> None.  Move to new file position.

        Argument offset is a byte count.  Optional argument whence defaults to
        0 (offset from start of file, offset should be >= 0); other values are 1
        (move relative to current position, positive or negative), and 2 (move
        relative to end of file, usually negative, although many platforms allow
        seeking beyond the end of a file).  If the file is opened in text mode,
        only offsets returned by tell() are legal.  Use of other offsets causes
        undefined behavior.
        Note that not all file objects are seekable.
        """
        return self.stream.seek(offset, whence)

    def all_yaml_docs(self):
        """
        Returns a list with all the YAML docs found. Seek the stream before returning.
        """
        docs = all_yaml_docs(self.stream)
        self.seek(0)
        return docs

    # Python 3 compatibility
    def __next__(self):
        return self.next()

    def next(self):
        """Return the next YAML document in stream."""
        return next_yaml_doc(self.stream)

    def next_doc_with_tag(self, doc_tag):
        """
        Returns the next document with the specified tag.
        Empty string is no doc is found.
        """
        for doc in self:
            if doc_tag in doc:
                return doc
        else:
            return ""

    def all_docs_with_tag(self, doc_tag, seek=True):
        """
        Returns all the documents with the specified tag.
        """
        docs = []

        while True:
            try:
                doc = self.next_doc_with(doc_tag)
                docs.append(doc)

            except StopIteration():
                break

        if seek:
            self.seek(0)
        
        return docs


def all_yaml_docs(stream):
    """
    Returns a list with all the YAML documents found in stream.

    .. warning:

        Assume that all the YAML docs (with the exception of the last one) 
        are closed explicitely with the sentinel '...'
    """
    docs, in_doc = [], False

    for line in stream:
        if line.startswith("---"):
            doc, in_doc = [], True

        if in_doc:
            doc.append(line)

        if in_doc and line.startswith("..."):
            in_doc = False
            docs.append("".join(doc))
            doc = []

    if doc:
        docs.append("".join(doc))

    return docs


def next_yaml_doc(stream, doc_tag="---"):
    """
    Returns the first YAML document in stream.

    .. warning:

        Assume that the YAML document are closed explicitely with the sentinel '...'
    """
    in_doc, lines = None, []
    for i, line in enumerate(stream):
        if line.startswith(doc_tag):
            in_doc = True

        if in_doc:
            lines.append(line)

        if in_doc and line.startswith("..."):
            break

    if lines:
        return "".join(lines)
    else:
        raise StopIteration()


def yaml_read_kpoints(filename, doc_tag="!Kpoints"):

    with YamlFileReader(filename) as r:
        doc = r.next_doc_with_tag(doc_tag)
        #doc = doc.replace(doc_tag, "")
        d = yaml.load(doc)

        return np.array(d["reduced_coordinates_of_qpoints"])
        #return KpointList(reciprocal_lattice, frac_coords, weights=None, names=None)


def yaml_read_irred_perts(filename, doc_tag="!IrredPerts"):

    with YamlFileReader(filename) as r:
        doc = r.next_doc_with_tag(doc_tag)
        #doc = doc.replace(doc_tag, "")
        print("doc",doc)
        d = yaml.load(doc)

        return d ["irred_perts"]


if __name__ == "__main__":

    string = """
---
none: [~, null]
bool: [true, false, on, off]
int: 42
float: 3.14159
list: [LITE, RES_ACID, SUS_DEXT]
dict: {hp: 13, sp: 5}
...

this is not a YAML document!
and the reader will ignore it

--- !Monster
name: Cave spider
hp: [2,6]    # 2d6
ac: 16
attacks: [BITE, HURT]
...

This is not a proper document since it does not start with --- 
the end tag below is ignored 
...
--- !Monster
name: Dragon
hp: [2,6]    # 2d6
ac: 32
attacks: [BITE, HURT]
...
"""
    print(string)

    filename = "foo.yaml"
    with open(filename, "w") as fh:
        fh.write(string)

    with YamlFileReader(filename) as r:

        # Read all docs present in file.
        all_docs = r.all_yaml_docs()
        print(all_docs)
        assert len(all_docs) == 3

        # We should be at the begining at the file.
        assert all_docs == r.all_yaml_docs()

        # Generate the docs
        r.seek(0)
        for i, doc in enumerate(r):
            print("doc", doc, "all", all_docs[i])
            assert doc == all_docs[i]

        # Find documents by tag.
        r.seek(0)
        monster = r.next_doc_with_tag("!Monster")
        assert monster == all_docs[1]

        monster = r.next_doc_with_tag("!Monster")
        assert monster == all_docs[2]

        monster = r.next_doc_with_tag("!Monster")
        assert monster == ""
    
