#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l os=RedHat7
#$ -pe smp 4
#$ -j y
#$ -R y
#$ -q gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda

source /cil/shed/apps/external/seQuoia/miniconda3/bin/activate /cil/shed/apps/internal/seQc_GAEMR/conda_env/

if [ $# -lt 1 ]; then 
    echo "usage: ${0} sample_list.txt"
    exit 1
fi

SAMPLE=($(sed -n ${SGE_TASK_ID}p "$1"))

set -euxo pipefail

sample_name=${SAMPLE[0]%_1.fastq.gz}
sample_name=${sample_name#data/}
assembly_dir=assemblies/${sample_name}

mkdir -p "log/${sample_name}"

# Redirect all stdout/stderr to a file
exec > log/${sample_name}/gaemr.log
exec 2>&1

cd $assembly_dir

/cil/shed/apps/external/seQuoia/seQc_GAEMR/GAEMR/bin/GAEMR.py \
    -t 4 -s assembly.scaffolds.fasta -c assembly.contigs.fasta -a assembly.agp \
    -m spades -l gaemr_config.csv -b dc-megablast --analyze_rna
