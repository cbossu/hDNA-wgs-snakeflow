#!/bin/bash
#SBATCH --job-name=mapdamage
#SBATCH --output=logs/mapdamage_%A_%a.out
#SBATCH --error=logs/mapdamage_%A_%a.err
#SBATCH --time=04:00:00
#SBACTH --partition=amilan128c,csu,amilan
#SBATCH --account=csu96_alpine1
#SBATCH --cpus-per-task=10
#SBATCH --mem=37G
#SBATCH --array=0-46

# source activate your_mapdamage_env

BAM_FILES=(results/bqsr-round-0/rmdup/*_rmdup.bam)
BAM=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM" _rmdup.bam)

RESCALED_BAM="results/bqsr-round-0/mapdamage2/rescaled_bams/${SAMPLE}.rescaled.bam"

if [[ -f "$RESCALED_BAM" ]]; then
  echo "Rescaled BAM for sample $SAMPLE already exists at $RESCALED_BAM. Skipping."
else
  mamba run -n historical-prep mapDamage \
    --input "$BAM" \
    --reference resources/genome.fasta \
    --folder results/bqsr-round-0/mapdamage2/"$SAMPLE" \
    --rescale-out "$RESCALED_BAM" \
    --rescale
fi
