version 1.0

# About minimal coverage:
# https://github.com/fritzsedlazeck/Sniffles/issues/137#issuecomment-479144579
task RunSniffles {
    input {
        File alignment
        File alignment_index
        File reference
        File reference_index

        Int genotype_ploidy
        Boolean stub = false
    }

    String outname = basename(basename(alignment, ".cram"), ".bam")

    command <<<
        set -e
        sniffles --version > version.txt
        if [ ~{stub} == "true" ]; then
            echo "Stubbing out Sniffles"
            touch ~{outname}.sniffles.vcf.gz
            touch ~{outname}.sniffles.vcf.gz.tbi
            exit 0
        fi
        sniffles --input ~{alignment} --genotype-ploidy ~{genotype_ploidy} --reference ~{reference} --vcf ~{outname}.sniffles.vcf.gz
    >>>

    runtime {
        docker: "quay.io/biocontainers/sniffles:2.2--pyhdfd78af_0"
        memory: "8 GB"
        cpu: 4
        preemptible: 3
    }

    output {
        File vcf = "~{outname}.sniffles.vcf.gz"
        File vcf_index = "~{outname}.sniffles.vcf.gz.tbi"
        File version = "version.txt"
    }
}
