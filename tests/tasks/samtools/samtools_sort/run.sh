touch alignment.unsorted.bam

miniwdl run --task Sort -i tests/tasks/samtools/samtools_sort/inputs.json tasks/samtools/samtools_sort.wdl
