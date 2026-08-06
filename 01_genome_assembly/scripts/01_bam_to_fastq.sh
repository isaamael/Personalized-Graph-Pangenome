#!/bin/bash
# Convert the PacBio HiFi BAM files to compressed FASTQ for read filtering and assembly.

bam2fastq -o input/TS input/TS.hifi_reads.bam
bam2fastq -o input/MM input/MM.hifi_reads.bam
