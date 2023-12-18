touch reference.fasta

miniwdl run --task HISAT2Build -i tests/tasks/hisat2/index_build/inputs.json tasks/hisat2/index_build.wdl
