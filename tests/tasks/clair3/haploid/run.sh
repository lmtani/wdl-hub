 touch \
    reference.fasta \
    reference.fasta.fai \
    alignment.cram \
    alignment.cram.crai \
    model.tar

 miniwdl run --task Clair3Haploid -i tests/tasks/clair3/haploid/inputs.json tasks/clair3/haploid.wdl
