version 1.0

# Adapted from: https://github.com/HumanCellAtlas/skylab/blob/e31492cd0219ff6f236cd0500401004f16f0fe41/library/tasks/HISAT2.wdl

task HISAT2PairedEnd {
  input {
    File hisat2_ref
    File fastq1
    File fastq2
    String ref_name
    String output_basename
    String sample_name

    Boolean stub = false

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
  Int machine_mem_mb = 16500
  Int cpu = 4
  # Using (fastq1 + fastq2) x 100 gives factor of a few buffer. BAM can be up to ~5 x (fastq1 + fastq2).
  # Need room for unsorted + sorted bam + temp sorting space + zipped and unzipped ref. Add 10 GiB buffer.
  Int disk = ceil((size(fastq1, "GiB") + size(fastq2, "GiB")) * 100 + size(hisat2_ref, "GiB") * 2 + 10)
  Int preemptible = 5
}
  meta {
    description: "HISAT2 alignment task will align paired-end fastq reads to reference genome."
  }

  parameter_meta {
    hisat2_ref: "HISAT2 reference"
    fastq1: "gz forward fastq file"
    fastq2: "gz reverse fastq file"
    ref_name: "the basename of the index for the reference genome"
    output_basename: "basename used for output files"
    sample_name: "sample name of input"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }



  command {
    # Note that files MUST be gzipped or the module will not function properly
    # This will be addressed in the future either by a change in how Hisat2 functions or a more
    # robust test for compression type.

    set -e
    hisat2 --version | grep hisat2 > version.txt

    if [ ~{stub} == "true" ]; then
        touch ~{output_basename}.bam ~{output_basename}.bam.bai ~{output_basename}.log ~{output_basename}.hisat2.met.txt
        exit 0
    fi


    tar --no-same-owner -xvf "~{hisat2_ref}"

    # run HISAT2 to genome reference with dedault parameters
    # --seed to fix pseudo-random number and in order to produce deterministics results
    # --secondary reports secondary alignments for multimapping reads. -k 10
    # searches for up to 10 primary alignments for each read
    hisat2 -t \
      -x ~{ref_name}/~{ref_name} \
      -1 ~{fastq1} \
      -2 ~{fastq2} \
      --rg-id=~{sample_name} --rg SM:~{sample_name} --rg LB:~{sample_name} \
      --rg PL:ILLUMINA --rg PU:~{sample_name} \
      --new-summary --summary-file ~{output_basename}.log \
      --met-file ~{output_basename}.hisat2.met.txt --met 5 \
      --seed 12345 \
      -k 10 \
      --secondary \
      -p ~{cpu} -S >(samtools view -1 -h -o ~{output_basename}_unsorted.bam)

    samtools sort -@ ~{cpu} -O bam -o "~{output_basename}.bam" "~{output_basename}_unsorted.bam"
    samtools index "~{output_basename}.bam"
  }

  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file = "~{output_basename}.log"
    File met_file = "~{output_basename}.hisat2.met.txt"
    File output_bam = "~{output_basename}.bam"
    File bam_index = "~{output_basename}.bam.bai"
    File version = "version.txt"
  }
}
