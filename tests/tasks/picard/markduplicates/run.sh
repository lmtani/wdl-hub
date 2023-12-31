touch alignment1.bam alignment2.bam
miniwdl run --task MarkDuplicates -i tests/tasks/picard/markduplicates/inputs.json tasks/picard/markduplicates.wdl
