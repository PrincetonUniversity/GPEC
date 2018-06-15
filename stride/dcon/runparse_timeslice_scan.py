#import asg_tools as asg
from asg_tools import *
import os.path
import sys
import numpy as np
import subprocess
import threading
import random

###################################################################################################################
###################################################################################################################
#MAIN SETTINGS

#operation_type = "run"
operation_type = "parse"
#operation_type = "spectrum"
#operation_type = "run_thread_scan"
#operation_type = "parse_thread_scan"

#eqtype_str = "CAKE"
#eqtype_str = "EFIT02"
#eqtype_str = "EFIT03"
eqtype_str = "EFIT04"
#eqtype_str = "EFIT02+03"
#eqtype_str = "HANDFIT"
#eqtype_str = "NSTX"

#shot_num = 166738
#shot_num = 163518
#shot_num = 148798
shot_num = 156746
#shot_num = 150312
#shot_num = 204118
#shot_num = 139285

equilibria_folder = "/home/asg5/stride_1.0/equilibria/new_CAKE_tests/"
output_folder = "./parse_NEW_timeslice_scan_output_novac_dp/"

###################################################################################################################
###################################################################################################################
#CHECKS AND TIMESLICE-SETTING
##############################################################
dotstr = ".0"
if shot_num == 166738:
    #date_str = "20180212"
    #date_str = "20180215"
    date_str = "20180214" #*For paper data
    time_range = map(str,range(2000,2505,5))
elif shot_num == 163518:
    ##For CAKE, EFIT03, et. al.
    #date_str = "20180210" #*For paper data
    #time_range = map(str,range(2045,4045,40))

    #For HANDFIT
    date_str = "20171017" #*For HANDFIT paper data
    time_range = map(str,[2045,2245,2250,2255,2305,2390,2850,3120,3455,3460,3470,3500,3510,3525,3640,3775,3780,3825,3905,3990])
elif shot_num == 148798:
    #date_str = "20180213"
    #date_str = "20180215"
    date_str = "20180214" #*For paper data
    time_range = map(str,range(4700,5205,5))
elif shot_num == 156746:
    date_str = "20180210"
    #date_str = "Matthijs"
    time_range = map(str,range(4530,4785,5))
elif shot_num == 149403:
    date_str = "20171001"
    time_range = map(str,range(3175,5175,20))
elif shot_num == 150312:
    date_str = "20170728"
    if operation_type == "run":
        time_range = map(str,range(2620,3220,20))
    elif operation_type == "parse":
        time_range = map(str,range(2700,3220,20))
    else:
        print("Not sure.")
        sys.exit(0)
elif shot_num == 204118:
    date_str = "201604"
    dotstr = "."
    time_range = ['00135', '00146', '00157', '00168', '00179', '00190', '00201', '00212', '00223', '00234', '00245', '00256', '00267', '00278', '00289', '00300', '00311', '00322', '00333', '00344', '00355', '00366', '00377', '00388', '00399', '00410', '00421', '00432', '00443', '00454', '00465', '00476', '00487', '00498', '00509', '00520', '00531', '00542', '00553', '00564', '00575', '00586', '00597', '00608', '00619', '00630', '00641', '00652', '00663', '00674', '00685', '00696', '00707', '00718', '00729', '00740', '00751', '00762', '00773', '00784', '00795', '00806', '00817', '00828', '00839', '00850', '00861', '00872', '00883', '00894', '00905', '00916', '00927', '00938', '00949', '00960', '00971', '00982', '00993', '01004', '01015', '01026', '01037', '01048', '01059', '01070', '01081', '01092', '01103', '01114', '01125', '01136', '01147', '01158', '01169', '01180', '01191', '01202', '01213', '01224', '01235', '01246']
elif shot_num == 139285:
    date_str = "201007"
    dotstr = "."
    time_range = map(str,range(121,1033,8))
    for i,t in enumerate(time_range):
        while len(t) < 5:
            t = "0"+t
        time_range[i] = t
else:
    print("Unknown shot number!")
    sys.exit(0)
##############################################################
if operation_type == "run":
    print("Running files...")
elif operation_type == "parse":
    print(eqtype_str + " #" + str(shot_num))
    array_out = np.chararray((1, 1))
elif operation_type == "spectrum":
    #rand_ifile = random.randrange(0,len(time_range),1)
    #rand_ifile = 10
    rand_ifile = time_range.index("4530")
    print(eqtype_str + " spectrum")
elif operation_type == "run_thread_scan":
    max_thread = 28
    rand_ifile = 10
    print("Scanning threads...")
elif operation_type == "parse_thread_scan":
    max_thread = 28
    rand_ifile = 10
    if eqtype_str != "EFIT03" or shot_num != 163518:
        print("Thread data was only generated for EFIT03 #163518.")
    else:
        print("Parsing threads...")
else:
    print("Unknown operation type!")
    sys.exit(0)
###################################################################################################################
###################################################################################################################
#SET UP THE RUN

#Get the list of filenames
#file_list = [("g" + str(shot_num) + ".0" + t) for t in time_range]
file_list = [("g" + str(shot_num) + dotstr + t) for t in time_range]
file_list = sorted(file_list)

data_folder = eqtype_str+"_"+str(shot_num)+"_"+date_str+"/"
folder_path = equilibria_folder+data_folder
out_folder = output_folder+data_folder
os.system("mkdir -p " + out_folder)
os.system("cp -f ./dcon.in " + out_folder + "dcon.in")
os.system("cp -f ./equil.in " + out_folder + "equil.in")
os.system("cp -f ./vac.in " + out_folder + "vac.in")

out_file = "dcon.out"
deltap_file = "delta_prime.out"
print_file = "dcon_print.out_"
if operation_type == "run_thread_scan":
    dcon_cmd1 = "./dcon 200 "
    dcon_cmd2 = " > "
else:
    dcon_cmd = "./stride 100 40 > "

file_list_len = len(file_list)
for ifile in range(0, file_list_len):
    gfile_str = file_list[ifile]

    gprint_file = print_file+gfile_str

    file_dp = out_folder+deltap_file+"_"+gfile_str+eqtype_str
    file_out = out_folder+out_file+"_"+gfile_str+eqtype_str
    file_print = out_folder+print_file+gfile_str+eqtype_str

    if operation_type == "run":
###################################################################################################################
###################################################################################################################
#RUN FILES

        os.system("cp -f "+folder_path+gfile_str+" g")
        print("Copied "+gfile_str+" to g")

        dcon_file_cmd = dcon_cmd+print_file+gfile_str
        os.system(dcon_file_cmd)

        os.system("cp -f " + gprint_file + " " + file_print)
        os.system("rm -f " + gprint_file)
        os.system("cp -f " + deltap_file + " " + file_dp)
        os.system("rm -f " + deltap_file)
        os.system("cp -f " + out_file + " " + file_out)
        os.system("rm -f " + out_file)

    elif operation_type == "parse":
###################################################################################################################
###################################################################################################################
#PARSE FILES

        #Parse screen-printed (console-printed) DCON output.
        if os.path.isfile(file_print):
            f = open(file_print, "r")
            lines = f.readlines()
            f.close()

            searchtxt0 = "q = "
            searchtxt1 = "q =   1.00000000"
            searchtxt2 = "q =   2.00000000"
            searchtxt3 = "q =   3.00000000"
            searchtxt4 = "q =   4.00000000"
            searchtxt5 = "q =   5.00000000"
            searchtxt6 = "q =   6.00000000"
            searchtxt7 = "q =   7.00000000"
            searchtxt8 = "q =   8.00000000"
            searchtxt9 = "q =   9.00000000"
            for i, line in enumerate(lines):
                line_clean = line.strip(" \t\n\r")
                if left(line_clean,len(searchtxt0)) == searchtxt0:
                    searchtxtB = "q =  "
                    for k in range(1,21):
                        if k < 10:
                            searchtxt = searchtxtB+" "+str(k)+".000000"
                        else:
                            searchtxt = searchtxtB+str(k)+".000000"
                        if left(line_clean,len(searchtxt)) == searchtxt:
                            row1_q = k #1L is the first row
                            break
                    if k == 20:
                        print("Error, q line was: "+line_clean)
                        sys.exit(0)
                    else:
                        break

            searchtxt0 = "*** calc-modes-for-wp time="
            searchtxt1 = "*** wp-calc time="
            for i, line in enumerate(lines):
                line_clean = line.strip(" \t\n\r")
                if left(line_clean,len(searchtxt0)) == searchtxt0:
                    k=1
                    while not left(lines[i+k],2).isspace():
                        k=k+1
                        if i+k == len(lines):
                            print("Could not find first eigenvalue.")
                            sys.exit(0)
                    maxWpSingVal = lines[i+k]
                    maxWpSingVal = maxWpSingVal.rstrip()
                    maxWpSingVal = maxWpSingVal.strip(" \t\n\r")
                elif left(line_clean,len(searchtxt1)) == searchtxt1:
                    k=1
                    while not left(lines[i-k],2).isspace():
                        k=k+1
                        if i-k == 0:
                            print("Could not find last eigenvalue.")
                            sys.exit(0)
                    minWpSingVal = lines[i-k]
                    minWpSingVal = minWpSingVal.rstrip()
                    minWpSingVal = minWpSingVal.strip(" \t\n\r")
                    break
        else:
            row1_q = "NaN"
            minWpSingVal = "NaN"
            maxWpSingVal = "NaN"
#            print("The parse operation was not directed to the correct files.")
#            sys.exit(0)

###################################################################################################################

        #Now parse the delta_prime matrix.
        if os.path.isfile(file_dp):
            f = open(file_dp, "r")
            lines = f.readlines()
            f.close()

            dpDiag = []
            for i, line in enumerate(lines):
                #The D' output file has a text row first, but no leading text column. First row is #0.
                if i >= 1:
                    line_clean = line.strip(" \t\n\r")
                    parts = line_clean.split("|")
                    for ii, part in enumerate(parts):
                        parts[ii] = clean_string(parts[ii])
                    dpDiag.append(parts[i-1]) #Gets the diagonal elements of the D' matrix

        else:
            dpDiag = []

###################################################################################################################

        #Now parse the dcon.out file.
        minPlasmaEnergyEig = "NaN"
        minTotEnergyEig = "NaN"
        if os.path.isfile(file_out):
            f = open(file_out, "r")
            lines = f.readlines()
            f.close()

            searchtxt = "Total Energy Eigenvalues:"
            for i, line in enumerate(lines):
                line_clean = line.strip(" \t\n\r")
                if left(line_clean,len(searchtxt)) == searchtxt:
                    j = i+4
                    while not lines[j].isspace():
                        lines[j] = clean_string(lines[j])
                        data_str = lines[j].split(" ")
                        eP = 1 #1=Plasma,2=Vacuum,3=Total
                        eT = 3 #1=Plasma,2=Vacuum,3=Total
                        data_str[eP] = float(data_str[eP])
                        data_str[eT] = float(data_str[eT])
                        if j == i+4:
                            plasmaEigs = [data_str[eP]]
                            totEigs = [data_str[eT]]
                        else:
                            plasmaEigs = np.vstack((plasmaEigs, data_str[eP]))
                            totEigs = np.vstack((totEigs, data_str[eT]))
                        j+=1

                    minPlasmaEnergyEig = str(min(plasmaEigs)[0])
                    minTotEnergyEig = str(min(totEigs)[0])

        #Add all data to array_out
        shot_time = right(gfile_str,5)
        out_row = np.asarray([shot_time, minTotEnergyEig, minPlasmaEnergyEig, minWpSingVal, maxWpSingVal, str(row1_q)+str(row1_q)+"L"] + dpDiag)
        if ifile == 0:
            array_out = out_row
            print("time,deltaW,deltaW_P,minWpEig,maxWpEig,dpMat_row1_q,DPdiag")
        elif len(out_row) > array_out.shape[1]:
            #The number of columns does not match.
            a0 = nans_str([array_out.shape[0],len(out_row)-array_out.shape[1]])
            array_out = np.concatenate((array_out, a0), axis=1)
        else:
            #The number of columns does not match.
            while len(out_row) < array_out.shape[1]:
                out_row = np.append(out_row,"NaN")
        array_out = np.vstack((array_out, out_row))

    elif operation_type == "spectrum":
###################################################################################################################
###################################################################################################################
#RANDOM SPECTRUM GENERATION

        #Print out spectrum
        eig = []
        if ifile == rand_ifile:
            if os.path.isfile(file_print):
                f = open(file_print, "r")
                lines = f.readlines()
                f.close()

                searchtxt0 = "Computing free"
                searchtxt1 = "vac loop time="
                for i, line in enumerate(lines):
                    line_clean = line.strip(" \t\n\r")
                    if left(line_clean,len(searchtxt0)) == searchtxt0:
                        i += 1
                        line_clean = lines[i].strip(" \t\n\r")
                        k = 0
                        while left(line_clean,len(searchtxt1)) != searchtxt1:
                            eig.append(lines[i+k])
                            eig[k] = eig[k].rstrip()
                            eig[k] = eig[k].strip(" \t\n\r")
                            k += 1
                            line_clean = lines[i+k].strip(" \t\n\r")
                        break
                print("Shot.Slice=" + str(shot_num) + "." + str(time_range[rand_ifile]))
                for e in eig:
                    print e

    elif operation_type == "run_thread_scan":
###################################################################################################################
###################################################################################################################
#SCAN OVER PARALLEL THREADS

        if ifile == rand_ifile:
            os.system("cp -f "+folder_path+gfile_str+" g")
            print("Copied "+gfile_str+" to g")

            for ithread in range(1,max_thread+1):
                dcon_file_cmd = dcon_cmd1+str(ithread)+dcon_cmd2+print_file+gfile_str
                os.system(dcon_file_cmd)

                os.system("cp -f " + gprint_file + " " + file_print + "_thread" + str(ithread))
                os.system("rm -f " + gprint_file)
                os.system("cp -f " + deltap_file + " " + file_dp + "_thread" + str(ithread))
                os.system("rm -f " + deltap_file)
                os.system("cp -f " + out_file + " " + file_out + "_thread" + str(ithread))
                os.system("rm -f " + out_file)

    elif operation_type == "parse_thread_scan":
###################################################################################################################
###################################################################################################################
#PARSE OVER PARALLEL THREADS

        if ifile == rand_ifile:
            for ithread in range(1,max_thread+1):
                file_thread = file_print+ "_thread" + str(ithread)
                if os.path.isfile(file_thread):
                    f = open(file_thread, "r")
                    lines = f.readlines()
                    f.close()

                    searchtxt = "*** full ode time="
                    for i, line in enumerate(lines):
                        line_clean = line.strip(" \t\n\r")
                        if left(line_clean,len(searchtxt)) == searchtxt:
                            line_clean = clean_string(line_clean)
                            data_str = line_clean.split(" ")
                            ode_time = data_str[-1]
                            print(str(ithread) + "," + str(ode_time))

if operation_type == "parse":
    #Cut out 0 from the range here because the first row is 0's
    for i in range(1,len(array_out)):
        print(",".join(array_out[i,:]))
