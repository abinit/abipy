#!/usr/bin/env python
import sys

from pymatgen.io.abinitio.events import EventsParser

def prof_main():
    import pstats, cProfile
    cProfile.runctx("main()", globals(), locals(), "Profile.prof")
    s = pstats.Stats("Profile.prof")
    s.strip_dirs().sort_stats("time").print_stats()
    return 0

def main():
    parser = EventsParser()
    events = parser.parse(sys.argv[1])
    print(events)

    return 0

if __name__ == "__main__":
    sys.exit(prof_main())
    
    
    
