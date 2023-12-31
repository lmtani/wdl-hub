touch reference.fasta

miniwdl run --task MakeBwaIndex -i tests/tasks/bwa/bwaindex/inputs.json tasks/bwa/bwaindex.wdl
