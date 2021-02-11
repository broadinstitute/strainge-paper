#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3
use .openmpi-uge-4.0.1

# strainest environment has python 2
source activate strainest

if [ $# -lt 1 ]; then 
    echo "usage: ${0} sample_list.txt"
    exit 1
fi

SAMPLE=($(sed -n ${SGE_TASK_ID}p "$1"))

set -euxo pipefail

sample_name=${SAMPLE[0]%/*}
sample_name=${sample_name#samples/}
dirname=${SAMPLE[0]%/*}
fname=${SAMPLE[0]##*/}

mkdir -p "log/${sample_name}"

# Redirect all stdout/stderr to a file
exec > log/${sample_name}/sigma.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}
CFG=${R1%.1.fq.gz}.sigma.cfg

TRIMMED_R1=${R1%.1.fq.gz}.trimmed.1.fq
TRIMMED_R2=${R2%.2.fq.gz}.trimmed.2.fq
SINGLES=${R1%.1.fq.gz}.singles.fq

output_dir=sigma/${sample_name}
mkdir -p ${output_dir}

START=$(date --rfc-3339=seconds)

sickle pe -f ${R1} -r ${R2} -q 20 -t sanger \
    -o ${TRIMMED_R1} -p ${TRIMMED_R2} -s ${SINGLES}

./sigma/bin/sigma-align-reads -c ${CFG} -w ${output_dir}
./sigma/bin/sigma -c ${CFG} -w ${output_dir}

PRE_BOOTSTRAP=$(date --rfc-3339=seconds)

# Can't get MPI programs to run
# mpirun -np 1 ./sigma/bin/sigma-bootstrap-mpi -t 2 -c ${CFG} -w ${output_dir}
# mpirun -np 1 ./sigma/bin/sigma-jackknife-mpi -t 2 -c ${CFG} -w ${output_dir}

rm ${TRIMMED_R1} ${TRIMMED_R2} ${SINGLES}

END=$(date --rfc-3339=seconds)
