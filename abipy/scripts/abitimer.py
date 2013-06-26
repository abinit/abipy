#!/usr/bin/env python
"""Script for analyzing the data saved by the ABINT timer in the TIMER section"""
from __future__ import print_function, division
import sys
import matplotlib.pyplot as plt

from abipy.htc.abitimer import AbiTimerParser, build_timer_parser
from abipy.tools.text import pprint_table

__version__ = "0.2"
__author__  = "Matteo Giantomassi"

######################################################################


def main():
    opt_parser = build_timer_parser()

    # Add files to the list of arguments.
    opt_parser.add_argument("files", nargs="+", help="File(s) to analyze")

    options = opt_parser.parse_args()

    parser = AbiTimerParser()
    parser.read(options.files)

    args = options._get_args()
    kwargs = dict(options._get_kwargs())
    print(args, kwargs)

    parser.main(*args, **kwargs)
    return 0

    if options.efficiency:
        parser.show_efficiency()

    if options.pie:
        parser.show_pie()

    if options.table:
        for timer in parser.timers(): 
            table = timer.totable(stop=options.table)
            c = "*"
            print(80*c + "\n" + str(timer) + "\n" + 80*c + "\n")
            pprint_table(timer.totable(stop=options.table))
            print(80*c + "\n")

    if options.csv:
        for timer in parser.timers():
            timer.tocsv()

    if options.histogram:
        for timer in parser.timers():
            timer.hist2()
        #for timer in parser.timers(): timer.cpuwall_histogram(title=timer.fname)

    if options.stacked_histogram:
        parser.show_stacked_hist()

    return 0

###############################################################################

if __name__ == "__main__":
  sys.exit(main())
