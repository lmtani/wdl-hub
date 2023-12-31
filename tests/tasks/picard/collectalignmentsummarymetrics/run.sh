touch alignment.cram alignment.cram.crai reference.fasta reference.fasta.fai

miniwdl run --task CollectAlignmentSummaryMetrics -i tests/tasks/picard/collectalignmentsummarymetrics/inputs.json tasks/picard/collectalignmentsummarymetrics.wdl
