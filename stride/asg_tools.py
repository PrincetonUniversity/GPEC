__all__ = ['left','right','mid','clean_string','nans','nans_str','replace_CHEASE_field_val','replace_in_field','isclose','gradient_asg','interp_asg','command','field']

import os
import sys
import numpy as np
import subprocess
import threading

def left(s, amount):
    return s[:amount]

def right(s, amount):
    return s[-amount:]

def mid(s, offset, amount):
    return s[offset:offset+amount]

def clean_string(s):
    sclean = s.strip(" \t\n\r")
    sclean = sclean.strip()
    while "  " in sclean:
        sclean = sclean.replace("  "," ")
    return sclean

def nans(shape, dtype=float):
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a

def nans_str(shape, dtype=object):
    a = np.empty(shape, dtype)
    a.fill("NaN")
    return a

def replace_CHEASE_field_val(file_str,field_str,valNew):
    f = open(file_str, "r")
    lines = f.readlines()
    f.close()
    for i, line in enumerate(lines):
        sSpot = line.find(field_str+"=")
        if sSpot != -1:
            cSpot = line.find(",",sSpot,len(line))
            if cSpot != -1:
                valStr = mid(line,sSpot+len(field_str)+1,cSpot-sSpot-len(field_str)-1)
                if valStr.find(".") != -1:
                    valNow = float(valStr)
                else:
                    valNow = int(valStr)
                nowStr = field_str+"="+str(valNow)
                nextStr= field_str+"="+str(valNew)

                with open(file_str, 'r') as file :
                    filedata = file.read()
                    filedata = filedata.replace(nowStr, nextStr)
                with open(file_str, 'w') as file:
                    file.write(filedata)

def replace_in_field(file_str,field_str,valNew):
    f = open(file_str, "r")
    lines = f.readlines()
    f.close()

    found_it = False
    eSpot = -1
    for i, line in enumerate(lines):
        line_clean = clean_string(line)
        if left(line_clean,len(field_str)+2) == field_str+"  =" or left(line_clean,len(field_str)+1) == field_str+"=":
            if found_it:
                print("Found two inputs for one field!")
                exit()
            else:
                found_it=True

            eSpot = line.find("=")
            leadStr = left(line,eSpot+1) #The "+1" includes the "=" sign

            exSpot = line.find("!",eSpot+1,len(line))
            if exSpot != -1:
                valStr = mid(line,eSpot+1,exSpot-1)
                oldStr = leadStr+valStr
            else:
                oldStr = line.rstrip(" \n\r")
            newStr = leadStr+str(valNew)

    with open(file_str, 'r') as file :
        filedata = file.read()
        filedata = filedata.replace(oldStr, newStr)
    with open(file_str, 'w') as file:
        file.write(filedata)

def isclose(a, b, rel_tol=1e-10, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def gradient_asg(x,f):
    #Calculates O(h^2) derivative.
    #################################
    hL = x[:-2] - x[1:-1]   #Must always subtract the point it is calculating!
    hL = np.append(x[1] - x[0], hL)
    hL = np.append(hL, x[-2] - x[-1])
    #################################
    hR = x[2:] - x[1:-1]
    hR = np.append(x[2] - x[0], hR)
    hR = np.append(hR, x[-3] - x[-1])
    #################################
    fL = f[:-2]
    fL = np.append(f[1], fL)
    fL = np.append(fL, f[-2])
    #################################
    fR = f[2:]
    fR = np.append(f[2], fR)
    fR = np.append(fR, f[-3])
    #################################

    f_1 = (1/(hL-hR)) * (  (hL/hR)*(fR-f) - (hR/hL)*(fL-f) )
    return f_1

def interp_asg(x_pt,x,y):
    if not np.all(np.diff(x) > 0):
        if np.all(np.diff(x[1:]) > 0):
            return(np.interp(x_pt,x[1:],y[1:]))
        else:
            print("Fed a non-sorted x vector to interp_asg")
            exit()
    else:
        return(np.interp(x_pt,x,y))

class command(object):
    def __init__(self, cmd, cwd, shell):
        self.cmd = cmd
        self.cwd = cwd
        self.shell = shell
        self.process = None

    def run(self, timeout):
        def target():
            print("Thread started")
            self.process = subprocess.Popen(self.cmd, cwd=self.cwd, shell=self.shell)
            self.process.communicate()
            print("Thread finished")

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print("Terminating process")
            self.process.terminate()
            thread.join()
        print(self.process.returncode)

class field:
    def __init__(self, name=None):
        self.name = name
        self.vec = []
