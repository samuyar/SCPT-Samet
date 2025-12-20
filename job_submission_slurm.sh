#!/usr/local_rwth/bin/zsh

#SBATCH -p c18m

### Job name
#SBATCH --job-name=Test00
### Send e-mail when job begins and ends

#SBATCH --mail-user=samet.uyar@rwth-aachen.de
#SBATCH --mail-type=ALL

### File for the output (or error)
#SBATCH --output=MPIJOB_OUTPUT.out.%j
#SBATCH -e MPIJOB_ERROR.err.%j

### Time your job needs to execute, hh:mm:ss
#SBATCH --time=10:00:00

### Memory your job needs per node, e. g. 250 MB = 250M, 1 Gigabyte = 1G
#SBATCH --mem=32G

### Use more than one node for parallel jobs on distributed-memory systems, e. g. 2
##SBATCH --nodes=1
##SBATCH -n 8

### Number of CPUS per task (for distributed-memory parallelisation, use 1)
#SBATCH --cpus-per-task=48

### Disable hyperthreading by setting the tasks per core to 1
##SBATCH --ntasks-per-core=10

### Number of processes per node, e. g. 6 (6 processes on 2 nodes = 12 processes in total)
##SBATCH --ntasks-per-node=12

##SBATCH --ntasks=1

### Load the same Python module
module load Python/3.12.3

### Add your local bin to PATH so pip/pkgs are visible
export PATH=$HOME/.local/bin:$PATH

### Make sure Python knows where to find user packages
export PYTHONUSERBASE=$HOME/.local

### Add the site-packages explicitly (optional but recommended) 
export PYTHONPATH=$HOME/.local/lib/python3.12/site-packages:$PYTHONPATH

### Change to working directory
cd /rwthfs/rz/cluster/home/qp605411/SCPT-00

### Run your parallel application
python3 /rwthfs/rz/cluster/home/qp605411/SCPT-00/bruteforce_local.py