bashjob="""#!/bin/bash

#SBATCH --job-name={name:}

#SBATCH --ntasks={ntasks:}
#SBATCH --mem={mem:}M
#SBATCH --time={days:}-{hours:}:00:00
#SBATCH --partition={partition:}

# --- export all my environment variables to the job
#SBATCH --export=ALL

# --- default combines both output and error streams into one file
#SBATCH --output={location:}.oe

# --- mail notifications (NONE, BEGIN, END, FAIL, REQUEUE, ALL)
#SBATCH --mail-type={mailtype:}
#SBATCH --mail-user={mailuser:}

# ------------------------------------------------------------

echo "The job's id is: $SLURM_JOBID"
echo "The job's partition is: $SLURM_JOB_PARTITION"
echo "The master node of this job is: $SLURM_SUBMIT_HOST"
echo -n 'Started job at : ' ; date
echo " "

cd {location:}

module purge
module load gpec/{version:}
{exes:}
 
rm *.dat
rm draw*

echo " "
echo -n 'Ended job at  : ' ; date
echo " " 
exit

"""
