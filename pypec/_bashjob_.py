bashjob="""#!/bin/bash

#SBATCH --job-name=jobnamehere

#SBATCH --nodes=1
#SBATCH --mem=memhereM
#SBATCH --time=40:00:00
#SBATCH --partition=mque

# --- export all my environment variables to the job
#SBATCH --export=ALL

# --- default combines both output and error streams into one file
#SBATCH --output=runlochere.oe

# --- mail notifications (NONE, BEGIN, END, FAIL, REQUEUE, ALL)
#SBATCH --mail-type=mailtypehere
#SBATCH --mail-user=mailuserhere

# ------------------------------------------------------------

echo "The job's id is: $SLURM_JOBID"
echo "The master node of this job is: $SLURM_SUBMIT_HOST"
echo -n 'Started job at : ' ; date
echo " "

cd runlochere

exelisthere
 
rm *.dat
rm draw*

echo " "
echo -n 'Ended job at  : ' ; date
echo " " 
exit
"""
