# wipp.py
#
# Builds and runs the full WIPP code. 
# I can't believe I'm writing this shit yet again
#
# Mmm... spaghetti.
#

import os
import shutil
import datetime
import random
import string
import time
import numpy as np
import itertools

freqs = [200,240,289,347,418,502,603,725,872,1048,1259,1514,1819,2187,2629,3160,3798,4565,5487,6596,7928,9530,11455,13769,16550,19893,23912,28742,34549,41528,49916,60000]
I0    = -100000.0
center_lat = 35.0
L_targ= 3.0

# Working directory

root_dir = '/shared/users/asousa/WIPP/WIPPy/python/c'
ray_dir  = '/shared/users/asousa/WIPP/WIPPy/rays/full_kp0_20sec'
out_dir  = os.path.join(root_dir,'out_35')


# Compile C code:
print "---- COMPILING ----"
os.system("gcc -o calc_scattering calc_scattering.c -lm")

# Move compiled files to out dir

shutil.copy(root_dir + "/calc_scattering", out_dir + "/calc_scattering")
shutil.copy(root_dir + "/run_job.pbs", out_dir + "/run_job.pbs")

print "---- SCATTERING ----"

prefix = "scatter"
qq = "batch"

for i in xrange(len(freqs) - 1):


    f_low = freqs[i]
    f_high= freqs[i+1]

    jobname = prefix + "_%s" % (f_low)
    logfile = root_dir + "/logs/log_scatter_%s.txt" % (f_low)

    
    # The command to be run at each process
    job_cmd = "./calc_scattering %s %s %s %s %s %s" %(ray_dir, I0, center_lat, f_low, f_high, L_targ)

    # The command to start a batch process
    os.system("qsub -N %s -j oe -o %s -l nodes=1:ppn=1 -l walltime=24:00:00 -q %s " % (jobname, logfile, qq) +
              "-v run_path=%s,"                                                     % (out_dir) +
              "cmd=\"%s\" "                                                         % (job_cmd) +
              "%s/run_job.pbs" % out_dir);


