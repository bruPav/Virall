# Virall HPC Memory Tuning Notes (2026-03-04)

## Why this happens

The `ASSEMBLE` step (especially SPAdes) is usually the highest-memory stage in Virall.
For large pooled single-cell FASTQ inputs, SPAdes memory usage can increase sharply as
read volume and k-mer complexity increase. If free/allocated RAM is lower than required,
SPAdes exits with a memory error.

## What this is (and is not)

Increasing ASSEMBLE memory is usually **resource tuning**, not a tool update.
In HPC runs, stability depends on matching both:

1. Nextflow process memory request for `ASSEMBLE`
2. Scheduler-level memory allocation (e.g. Slurm `--mem`)

If either is too low, jobs can fail even when pipeline logic is correct.

## Manual runs

This still applies when running SPAdes manually.
You must set both:

1. SPAdes memory (`-m`, GB)
2. Job/session memory in scheduler (`--mem`)

Example:

```bash
spades.py -1 R1.fastq.gz -2 R2.fastq.gz -t 16 -m 80 -o spades_out
```

If this runs inside Slurm, the job should request at least that much memory (typically more headroom).

## When tuning is commonly needed

- Very large read sets (especially pooled scRNA-seq)
- Large host/background content before filtering
- Multi-lane merged FASTQ files
- Running many samples concurrently on one node

## Recommended policy

Short term:

- Document ASSEMBLE memory overrides clearly (README/run params)
- Instruct users to align Slurm `--mem` with ASSEMBLE memory settings

Long term:

- Keep automatic memory scaling based on input size
- Optionally provide an assembly-skip/lightweight path for ultra-large runs
