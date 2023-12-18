touch reference.fasta reference.fasta.fai alignment.cram alignment.cram.crai a_model.tar

miniwdl run --task DeepVariantCustomModel -i tests/tasks/deepvariant/custom/inputs.json tasks/deepvariant/custom.wdl
