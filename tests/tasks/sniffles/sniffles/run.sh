touch alignment.bam alignment.bam.bai reference.fasta reference.fasta.fai

miniwdl run --task RunSniffles -i tests/tasks/sniffles/sniffles/inputs.json tasks/sniffles/sniffles.wdl
