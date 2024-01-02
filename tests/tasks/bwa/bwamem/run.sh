touch reference.fasta \
      reference.fasta.fai \
      reference.fasta.fai \
      reference.fasta.amb \
      reference.fasta.ann \
      reference.fasta.bwt \
      reference.fasta.pac \
      reference.fasta.sa \
      reads_1.fastq.gz \
      reads_2.fastq.gz

miniwdl run --task BwaMem -i tests/tasks/bwa/bwamem/inputs.json tasks/bwa/bwamem.wdl
