#!/bin/bash

# --- the name of your job
#PBS -N _DIIID_kinetic_example
 
# --- Mail options: execution("b"), termination ("e"), or interruption ("a")
#PBS -m ae
# --- emailhere
  
#PBS -l nodes=1
#PBS -l mem=10000mb
#PBS -l walltime=40:00:00
#PBS -r n
#PBS -q mque

# --- export all my environment variables to the job
#PBS -V

# --- combine both output and error streams into one file with the
#PBS -o /p/gpec/GPEC-1.0/docs/examples/DIIID_kinetic_example.oe
#PBS -e /p/gpec/GPEC-1.0/docs/examples/DIIID_kinetic_example.err
#PBS -j oe

# ------------------------------------------------------------
# Log interesting information
#
echo " "
echo "-------------------"
echo "This is a $PBS_ENVIRONMENT job"
echo "This job was submitted to the queue: $PBS_QUEUE"
echo "The job's id is: $PBS_JOBID"
echo "-------------------"
echo "The master node of this job is: $PBS_O_HOST"

# ------------------------------------------------------------
 
echo -n 'Started job at : ' ; date
 
cd /p/gpec/GPEC-1.0/docs/examples/DIIID_kinetic_example

/p/gpec/GPEC-1.0/bin/dcon 
/p/gpec/GPEC-1.0/bin/gpec 
/p/gpec/GPEC-1.0/bin/pentrc 

 
rm *.dat
rm draw*

echo -n 'Ended job at  : ' ; date
echo " " 
exit
