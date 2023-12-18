touch reference.fasta reference.fasta.fai alignment.cram alignment.cram.crai

miniwdl run --task DeepVariant -i tests/tasks/deepvariant/standard/inputs.json tasks/deepvariant/standard.wdl
