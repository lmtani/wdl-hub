touch alignment1.bam alignment2.bam alignment3.bam

miniwdl run --task Merge -i tests/tasks/samtools/samtools_merge/inputs.json tasks/samtools/samtools_merge.wdl
