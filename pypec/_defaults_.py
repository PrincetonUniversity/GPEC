"""
HARDCODED VALUES

Edit this file to set the defaults for your pypec package.

"""

######################### gpec run default keyword arguments
email    = 'jpark@pppl.gov'
mailon   = 'ae'
rundir   = 'LinuxLahey64'



######################### gui defaults

# starting inputs taken from
inputdir = "."
# inputs used
inputs   = ['equil','dcon','ipec','coil','pent']
# temp file location for brute editing (this is the temp dir in this dir)
import os
from string import join
tempdir  = join(os.path.abspath(__file__).split('/')[:-1],'/')+'/temp'
