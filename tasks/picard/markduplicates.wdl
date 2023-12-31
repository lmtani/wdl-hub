version 1.0

task MarkDuplicates {

  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    String gatk_container = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.2.4-1469632282"
    Boolean stub = false
  }

  Int input_size = ceil(2.0 * size(input_bams, 'G') + 10)

  command <<<
    set -e

    java -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics --version 2>&1| sed 's/Version:/Picard.MarkDuplicates: /g' > version.txt
    if [ ~{stub} == "true" ]; then
      echo "Stubbing out MarkDuplicates"
      touch ~{output_bam_basename}.bam
      echo "0.0" > duplicates.txt
      exit 0
    fi


    java -Xmx4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=" INPUT="  input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      CREATE_MD5_FILE=true

    grep -A 1 "^LIBRARY" ~{metrics_filename} | cut -f 9 | tail -n 1 > duplicates.txt
  >>>

  runtime {
    docker: gatk_container
    memory: "7 GB"
    disks: "local-disk " + input_size + " HDD"
  }

  output {
    File output_bam = "${output_bam_basename}.bam"
    Float duplicates = read_float("duplicates.txt")
    File version = "version.txt"
  }
}
