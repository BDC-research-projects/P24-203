#!/bin/bash


#$ -N get_fastq
#$ -q development.q
#$ -cwd
#$ -o get_fastq.out -e get_fastq.err
#$ -M kristina.benevides@gu.se
#$ -m bea
#$ -pe mpi 8

# Usage: run_fasterq_dump_parallel.sh <SRR_FILE> <OUTPUT_DIR> <NUM_THREADS>

srr_file=$1
out_dir=$2
threads=${3:-4}  # Default to 4 threads if not specified
sif=$4

if [ -z "$srr_file" ] || [ -z "$out_dir" ]; then
    echo "Usage: $0 <SRR_FILE> <OUTPUT_DIR> [NUM_THREADS] [SINGULARITY]"
    exit 1
fi

# Function to process a single SRR ID
process_srr() {
    local sif=$3
    local out_dir=$2
    local srr_id=$1

    echo "Processing $srr_id..."


    # Check if all expected files already exist
    all_files_exist=true
    for suffix in 1 2 3; do
        file="${out_dir}/${srr_id}_${suffix}.fastq.gz"
        if [ ! -f "$file" ]; then
            all_files_exist=false
            break
        fi
    done

    if [ "$all_files_exist" = true ]; then
        echo "All files for $srr_id already exist. Skipping..."
        return
    fi

    # Run fastq-dump into temporary directory
    singularity exec -B /medstore:/medstore "$sif" fasterq-dump --split-files --include-technical --outdir "$out_dir" "$srr_id"


    # Compress with pigz
    for file in "$out_dir/${srr_id}"_*; do
        pigz -p "${threads}" "$file"  # Compress using pigz with specified threads
    done


}

export -f process_srr  # Export function for parallel
export out_dir         # Export output directory for parallel
export sif

# Run the function in parallel
cat "$srr_file" | parallel -j "$threads" process_srr {} "$out_dir" "$sif"
