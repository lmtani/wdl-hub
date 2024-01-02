version 1.0

task ExtractRegionsWithNCoverage {
    input {
        File per_base_coverage
        Int min_coverage
        Boolean stub = false
    }

    String basename = basename(per_base_coverage, ".per-base.bed.gz")

    command <<<
        set -e -o pipefail
        bedtools --version > version.txt
        if [ ~{stub} = true ]; then
            echo "Stubbing out the command"
            touch ~{basename}_min_of_~{min_coverage}_coverage_regions.bed
            exit 0
        fi

        gzip -dc ~{per_base_coverage} \
            | awk -v OFS="\t" -v min_coverage=~{min_coverage} '$4 >= min_coverage { print }' \
            | bedtools merge -d 1 -c 4 -o mean -i - > ~{basename}_min_of_~{min_coverage}_coverage_regions.bed
    >>>

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"
        cpu: 1
        memory: "1G"
        disks: "local-disk 10 HDD"
    }

    output {
        File covered = "~{basename}_min_of_~{min_coverage}_coverage_regions.bed"
        File version = "version.txt"
    }
}
