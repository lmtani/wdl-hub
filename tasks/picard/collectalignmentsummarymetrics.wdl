version 1.0

task CollectAlignmentSummaryMetrics {
    input {
        File alignment
        File alignment_index
        File reference_fasta
        File reference_fasta_index
        Boolean stub = false
    }

    String basename = basename(basename(alignment, ".bam"), ".cram")
    Int disk_size = ceil(size(alignment, "GiB") + 10)

    command <<<
        set -e

        java -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics --version 2>&1| sed 's/Version:/Picard.CollectAlignmentSummaryMetrics: /g' > version.txt

        if [ ~{stub} == "true" ]; then
            echo "Stubbing out the command"
            touch ~{basename}_alignment_metrics.txt
            exit 0
        fi

        java -Xms4g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
          -R ~{reference_fasta} \
          -I ~{alignment} \
          -O ~{basename}_alignment_metrics.txt
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.27.1"
        preemptible: 3
        memory: "6 GiB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File metrics = "~{basename}_alignment_metrics.txt"
        File version = "version.txt"
    }
}
