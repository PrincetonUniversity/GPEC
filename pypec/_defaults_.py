"""
HARDCODED VALUES

Edit this file to set the defaults for your pypec package.

"""

import os

######################### gpec run default keyword arguments
email = os.getlogin() + "@pppl.gov"
mailon = "ae"

######################### gui defaults

# starting inputs taken from
inputdir = "."
# inputs used
inputs = ["equil", "dcon", "gpec", "coil", "pent"]
# temp file location for brute editing (this is the temp dir in this dir)
tempdir = "~/GPEC/temp"
