#import asg_tools as asg
from asg_tools import *
import os.path
import sys
import numpy as np
import subprocess
import threading
import random

###################################################################################################################
#SET UP THE RUN

#Get the list of filenames

out_folder = "./timer_scan_output/"
os.system("mkdir -p " + out_folder)
os.system("module load intel")
os.system("module load openmpi")
os.system("module load intel-vtune")
os.system("module load intel-advisor")
os.system("module load intel-python")

print_file = "dcon_print.out_"

nruns = 30
output_on_console = True
console_target = "thread_data.out"

minIntervals = 18
maxIntervals = 24
stepInterval = 1

minThreads = 18
maxThreads = 24
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

        for irun in range(0, nruns):
            if print_each_run:
                if output_on_console:
                    print(irun)
                else:
                    fh = open(console_target, "a")
                    fh.write(str(irun)+"\n")
                    fh.close

            file_print = file_print_head+"_"+str(irun)
            dcon_file_cmd = dcon_cmd+file_print
            os.system(dcon_file_cmd)

        for irun in range(0, nruns):
            file_print = file_print_head+"_"+str(irun)
            if os.path.isfile(file_print):
                f = open(file_print, "r")
                lines = f.readlines()
                f.close()

                itime_cat = []
                itime_timer = []

                searchtxt1 = "***"
                searchtxt2 = "time="
                for i, line in enumerate(lines):
                    line_clean = clean_string(line)
                    if left(line_clean,len(searchtxt1)) == searchtxt1:
                        if searchtxt2 in line_clean:
                            parts = line_clean.split(" ")
                            if parts[0] == "***" and len(parts) == 4:
                                itime_cat.append(parts[1])
                                itime_timer.append(float(parts[3]))
                            else:
                                error("Unknown time line in output!")

                data_list = zip(itime_cat,itime_timer)
                data_list = sorted(data_list)
                [itime_cat,itime_timer] = zip(*data_list)
                itime_cat = list(itime_cat)
                itime_cat.insert(0,"nThreads")
                itime_cat.insert(0,"nIntervals")
                itime_timer = list(itime_timer)
                itime_timer = np.insert(itime_timer,0,nThreads)
                itime_timer = np.insert(itime_timer,0,nIntervals)
                if irun == 0:
                    time_cat = itime_cat
                    time_array = np.asarray(itime_timer)
                    if print_each_run:
                        if output_on_console:
                            print(",".join(map(str,itime_cat)))
                            print(",".join(map(str,itime_timer)))
                        else:
                            fh = open(console_target, "a")
                            fh.write(",".join(map(str,itime_cat)))
                            fh.write(",".join(map(str,itime_timer)) + "\n")
                            fh.close
                    elif nIntervals == minIntervals and nThreads == minThreads:
                        if output_on_console:
                            print(",".join(map(str,itime_cat)))
                        else:
                            fh = open(console_target, "a")
                            fh.write(",".join(map(str,itime_cat)) + "\n")
                            fh.close
                elif itime_cat == time_cat:
                    time_array = np.vstack((time_array,np.asarray(itime_timer)))
                    if print_each_run:
                        if output_on_console:
                            print(",".join(map(str,itime_timer)))
                        else:
                            fh = open(console_target, "a")
                            fh.write(",".join(map(str,itime_timer)) + "\n")
                            fh.close
                else:
                    print("Irregular timing pattern for file: "+file_print)
                    sys.exit(0)
            else:
                print("File not found: "+file_print)
                sys.exit(0)

        if not print_each_run:
            if output_on_console:
                print(",".join(map(str,np.median(time_array,axis=0))))
            else:
                fh = open(console_target, "a")
                fh.write(",".join(map(str,np.median(time_array,axis=0))) + "\n")
                fh.close
