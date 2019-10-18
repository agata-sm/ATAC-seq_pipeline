#!/bin/bash

# Run snakemake
snakemake \
    --snakefile /proj/uppstore2019111/nobackup/private/nbis4833/tst/pipeline/ATAC-seq_pipeline/atacseq_proc_snakefile \
    --rerun-incomplete \
    --jobs 50 \
    --cluster-config cluster.yml \
    --cluster "sbatch \
                  -A {cluster.account} \
                  -t {cluster.time} \
                  -p {cluster.partition} \
                  -n {cluster.N} \
                  -J {cluster.jobname}"


