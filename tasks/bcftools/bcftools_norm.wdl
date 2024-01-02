version 1.0

task Norm {
    input{
        File vcf
        String extra_args = "-m-any"
        Boolean stub = false
    }

    String outname = basename(vcf, ".vcf.gz")

    command <<<
        set -e

        bcftools --version | grep bcftools > version.txt

        if [ ~{stub} == "true" ]; then
            echo "Stubbing out Norm"
            touch ~{outname}.norm.vcf.gz
            touch ~{outname}.norm.vcf.gz.tbi
            exit 0
        fi

        bcftools norm ~{extra_args} ~{vcf} -Oz -o ~{outname}.norm.vcf.gz
        bcftools index --tbi ~{outname}.norm.vcf.gz
    >>>

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        memory: "2 GB"
        cpu: 1
        preemptible: 3
        disks: "local-disk 15 SSD"
    }

    output {
        File norm_vcf = "~{outname}.norm.vcf.gz"
        File norm_vcf_index = "~{outname}.norm.vcf.gz.tbi"
        File version = "version.txt"
    }
}
