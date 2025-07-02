#!/bin/bash
#SBATCH --job-name=samtools_stats
#SBATCH --output=logs/samtools_stats_%A_%a.out
#SBATCH --error=logs/samtools_stats_%A_%a.err
#SBACTH --partition=amilan128c,csu,amilan
#SBATCH --account=csu96_alpine1
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=0-48

module load samtools

# define input directory
BAM_DIR="results/bqsr-round-0/rmdup"

# get list of BAM files into array
BAM_FILES=($(ls $BAM_DIR/*_rmdup.bam))

# pick the BAM file for this task
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# extract sample name (strip path and suffix)
SAMPLE=$(basename "$BAM_FILE" _rmdup.bam)

# define output directory and file
OUT_DIR="results/bqsr-round-0/qc/samtools_stats"
mkdir -p "$OUT_DIR"
OUT_FILE="$OUT_DIR/${SAMPLE}.txt"

# run samtools stats
samtools stats "$BAM_FILE" > "$OUT_FILE"

