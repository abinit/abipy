#!/usr/bin/env python
"""
Script to generate/check the md5 checksum associated to pseudopotential files.

If the file does not contain the md5, a new line of the form `#md5: checksum`
is addded at the end of the file.
If the checksum is already present, a consistency check is performed.
"""
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse
import json
import hashlib


def read_dojoreport(filepath):
    """Read the DojoReport from file."""
    with open(filepath, "rt") as fh:
        lines = fh.readlines()
        start = lines.index("<DOJO_REPORT>\n")
        stop = lines.index("</DOJO_REPORT>\n")

        return json.loads("".join(lines[start+1:stop]))


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    #epilog=str_examples(), 

    parser.add_argument('path', help="Pseudopotential file")
    options = parser.parse_args()
    path = options.path

    # Open file and look for md5 line and beginning of DOJO_REPORT
    dojoreport_start = -1
    with open(path, "rt") as fh:
        text = fh.readlines()
        md5 = None
        for i, line in enumerate(text):
            if line == "<DOJO_REPORT>\n": dojoreport_start = i
            if not line.startswith("#md5:"): continue
            md5 = line[5:].strip()
            break

    if md5 is not None:
        # Recompute hash from the text before `#md5:` 
        # Do not include DOJO_REPORT section if present.
        if dojoreport_start == -1:
            s = "".join(text[:i])
        else:
            s = "".join(text[:dojoreport_start])

        check = hashlib.md5(s.encode("utf-8")).hexdigest()

        # Perform consistency test.
        if md5 != check:
            msg = ("md5 checksum read from file %s does not agree with computed value:\n"
                   "file md5:     %s\n"
                   "computed md5: %s\n" % (path, md5, check))
            print(msg)
            return 1
        else:
            print("md5 %s read from file agrees with recomputed values." % md5)

    else:
        # Compute hash and append new line at the end of the file.
        # Do not include DOJO_REPORT section in the checksum if present.
        if dojoreport_start == -1:
            s = "".join(text)
        else:
            s = "".join(text[:dojoreport_start])

        md5 = hashlib.md5(s.encode("utf-8")).hexdigest()
        md5_line = "#md5: %s" % md5

        text.append(md5_line)
        print("Pseudo does not contain md5 record.") 
        print("Will add", md5_line, "at the end of the file")

        with open(path, "wt") as fh:
            fh.writelines(text)

    # If the pseudos contains a DOJO_REPORT section, 
    # we make sure that the md5 agrees with the one stored in the report.
    if dojoreport_start != -1:
        report = read_dojoreport(path)
        if md5 != report["md5"]:
            msg = ("pseudodojo checksum read from file %s does not agree with computed value:\n"
                   "pseudodojo md5: %s\n"
                   "computed md5:   %s\n" % (path, report["md5"], md5))
            print(msg)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

