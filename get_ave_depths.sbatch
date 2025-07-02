#!/bin/bash
#SBATCH --job-name=ave_depths
#SBATCH --output=logs/get_ave_depths.out
#SBATCH --error=logs/get_ave_depths.err
#SBACTH --partition=amilan128c,csu,amilan
#SBATCH --account=csu96_alpine1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Set variables
NUM_BASES=$(cat results/bqsr-round-0/DS_control/genome_length.txt)
OUT_FILE="results/bqsr-round-0/DS_control/sample_info.tsv"
STATS_DIR="results/bqsr-round-0/qc/samtools_stats"

#echo -e "sample\tave_depth" > "$OUT_FILE"

#for STATS_FILE in "$STATS_DIR"/*.txt; do
#    SAMPLE=$(basename "$STATS_FILE" .txt)
#    DEPTH=$(awk -v f="$SAMPLE" -v NumBases="$NUM_BASES" '
#      BEGIN {OFS="\t"}
#      $2 == "bases" && $3 == "mapped" && $4 == "(cigar):" { print f, $5 / NumBases }
#    ' "$STATS_FILE")
#    echo -e "${SAMPLE}\t${DEPTH}" >> "$OUT_FILE"
#done

echo -e "sample\tave_depth" > "$OUT_FILE"

for STATS_FILE in "$STATS_DIR"/*.txt; do
    SAMPLE=$(basename "$STATS_FILE" .txt)
    DEPTH=$(awk -v f="$SAMPLE" -v NumBases="$NUM_BASES" '
      BEGIN {OFS="\t"}
      $2 == "bases" && $3 == "mapped" && $4 == "(cigar):" { print f, $5 / NumBases }
    ' "$STATS_FILE")
    echo -e "$DEPTH" >> "$OUT_FILE"
done

