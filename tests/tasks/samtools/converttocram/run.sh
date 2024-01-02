touch alignment.bam reference.fasta reference.fasta.fai

miniwdl run --task ConvertToCram -i tests/tasks/samtools/converttocram/inputs.json tasks/samtools/converttocram.wdl
