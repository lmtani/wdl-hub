version 1.0

task Sort {
    input {
        File vcf
        Boolean stub = false
    }

    String outname = basename(vcf, ".vcf.gz")

    command <<<
        set -e

        bcftools --version | grep bcftools > version.txt
        if [ ~{stub} == "true" ]; then
            echo "Stubbing out the task"
            touch ~{outname}.sorted.vcf.gz
            touch ~{outname}.sorted.vcf.gz.tbi
            exit 0
        fi

        # keep only snps
        bcftools view -v snps ~{vcf} -Oz -o ~{outname}.snps.vcf.gz
        bcftools index --tbi ~{outname}.snps.vcf.gz

        # Sort the VCF file
        bcftools sort -Oz -o ~{outname}.sorted.vcf.gz ~{outname}.snps.vcf.gz
        bcftools index --tbi ~{outname}.sorted.vcf.gz
    >>>

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
    }

    output {
        File sorted_vcf = "~{outname}.sorted.vcf.gz"
        File sorted_vcf_index = "~{outname}.sorted.vcf.gz.tbi"
        File version = "version.txt"
    }
}
