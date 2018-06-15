#import asg_tools as asg
from asg_tools import *
import os.path
import sys
import numpy as np
import subprocess
import threading
import random

#THIS IS MEANT TO BE RUN WITH THE "naive" GRID-PACKING!
###################################################################################################################
#SET UP THE RUN

#Get the list of filenames

out_folder = "./naive_grid_timer_scan_output/"
os.system("mkdir -p " + out_folder)

print_file = "dcon_print.out_"

nruns = 40

minIntervals = 300
maxIntervals = 300
stepInterval = 1

minThreads = 28
maxThreads = 28
stepThread = 1

if minIntervals == maxIntervals and minThreads == maxThreads:
    print_each_run = True
else:
    print_each_run = False

for nIntervals in range(minIntervals,maxIntervals+1,stepInterval):
    for nThreads in range(minThreads,maxThreads+1,stepThread):
        os.system("rm -f "+out_folder+"dcon*")
        os.system("rm -f "+out_folder+"*.in")
        os.system("cp -f ./dcon.in " + out_folder + "dcon.in")
        os.system("cp -f ./equil.in " + out_folder + "equil.in")
        os.system("cp -f ./vac.in " + out_folder + "vac.in")

        dcon_cmd = "./dcon "+str(nIntervals)+" "+str(nThreads)+" > "
        file_print_head = out_folder+print_file+str(nIntervals)+"_"+str(nThreads)

        #Run code with output to file
        for irun in range(0, nruns):
            if print_each_run:
                print(irun)
            file_print = file_print_head+"_"+str(irun)
            dcon_file_cmd = dcon_cmd+file_print
            os.system(dcon_file_cmd)

        #Collect singular surface data
        file_print = file_print_head+"_0"
        if os.path.isfile(file_print):
            f = open(file_print, "r")
            lines = f.readlines()
            f.close()

            searchtxt = "Sing# ="
            ising = 0
            for i, line in enumerate(lines):
                line_clean = clean_string(line)
                if left(line_clean,len(searchtxt)) == searchtxt:
                    line_sing = clean_string(lines[i+3])
                    parts = line_sing.split(" ")
                    ising = ising+1
                    print("Sing #"+str(ising)+" = "+str(parts[-1]))

        #Collect interval runtime data
        itime_interval = []
        itime_timer = []
        itime_t0 = []
        itime_t1 = []
        for irun in range(0, nruns):
            file_print = file_print_head+"_"+str(irun)
            if os.path.isfile(file_print):
                f = open(file_print, "r")
                lines = f.readlines()
                f.close()

                for iInterval in range (1,maxIntervals+1):
                    searchtxt1 = str(iInterval)+","
                    for i, line in enumerate(lines):
                        line_clean = clean_string(line)
                        if left(line_clean,len(searchtxt1)) == searchtxt1:
                            parts = line_clean.split(",")
                            theInterval = int(parts[0])
                            theT0 = float(parts[1])
                            theT1 = float(parts[2])
                            theTime = float(parts[3])
                            itime_interval.append(theInterval)
                            itime_timer.append(theTime)
                            itime_t0.append(min(theT0,theT1)) #This handles reversed intervals (with min).
                            itime_t1.append(max(theT0,theT1)) #This handles reversed intervals (with max).

        #Now sort/median interval times...
        itime_interval = np.asarray(itime_interval)
        itime_timer = np.asarray(itime_timer)
        itime_t0 = np.asarray(itime_t0)
        itime_t1 = np.asarray(itime_t1)
        for iInterval in range (1,maxIntervals+1):
            i_t0 = np.median(itime_t0[itime_interval==iInterval])
            i_t1 = np.median(itime_t1[itime_interval==iInterval])
            i_median = np.median(itime_timer[itime_interval==iInterval])
            print(str(iInterval)+","+str(i_t0)+","+str(i_t1)+","+str(i_median))
