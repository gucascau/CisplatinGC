#!/bin/sh
#SBATCH --partition=himem
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --mail-user=xin.wang@nationwidechildrens.org
#SBATCH --job-name=Pipseeker
#SBATCH --time=40:00:00
#SBATCH --ntasks=12

set -e

ml purge

ml STAR/2.7.9a

id=Cisplatin-GC

pipseeker=/home/gdbecknelllab/xxw004/Software/PipSeeker3.3/pipseeker-v3.3.0-linux/pipseeker
starindex=/home/gdbecknelllab/xxw004/Datasets/pipseeker-gex-reference-GRCm39-2022.04
script=/home/gdbecknelllab/xxw004/gdmanlab/Projects/Databases/RawDZscRNAseq/Scripts

in=/home/gdbecknelllab/xxw004/gdmanlab/Projects/Databases/RawDZscRNAseq/${id}
out=/home/gdbecknelllab/xxw004/gdmanlab/Projects/Databases/RawDZscRNAseq/Results/${id}



cd ${in}

#gunzip *.gz

#
perl ${script}/ExtractSequenceCorrectedBaseOnLength.pl -i 54 -f ${id}_S4_L001_R1_001.fastq -r ${id}_S4_L001_R2_001.fastq -o ${id}_Cfiltered_L001

perl ${script}/ExtractSequenceCorrectedBaseOnLength.pl  -i 54 -f ${id}_S4_L002_R1_001.fastq -r ${id}_S4_L002_R2_001.fastq -o ${id}_Cfiltered_L002

#
gzip *.fastq

mkdir -p ${out}
cd ${out}

${pipseeker} full --chemistry v4 --fastq ${in}/${id}_Cfiltered   --output-path ${out} --star-index-path ${starindex} --threads 10