touch reads_1.fq.gz reads_2.fq.gz reference.tar

miniwdl run --task HISAT2PairedEnd -i tests/tasks/hisat2/paired_end/inputs.json tasks/hisat2/paired_end.wdl
