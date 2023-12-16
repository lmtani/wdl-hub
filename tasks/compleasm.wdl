version 1.0


task Download {
    input {
        String lineage  # BUSCO compatible lineage name
        String lineage_dir = "mb_downloads"
        String container = "quay.io/biocontainers/compleasm:0.2.2--pyh7cba7a3_0"
        Int disk_size = 10
        Boolean stub = false
    }

    command <<<
        if [ ~{stub} == "true" ]; then
            touch ~{lineage}.tar.gz
            exit 0
        fi

        compleasm download --library_path ~{lineage_dir} ~{lineage}
        tar -czvf ~{lineage}.tar.gz ~{lineage_dir}
    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        docker: container
        disk: "local-disk ~{disk_size} HDD"
    }

    output {
        File lineage_tar = "~{lineage}.tar.gz"
    }
}


task Run {
    input {
        File fasta
        File lineage_tar
        String lineage
        String output_directory
        Int threads = 4
        Int memory = 12
        String lineage_dir = "mb_downloads"
        String extra_args = ""
        String container = "quay.io/biocontainers/compleasm:0.2.2--pyh7cba7a3_0"
        Int disk_size = 10
        Boolean stub = false
    }

    command <<<
        set -e
        if [ ~{stub} == "true" ]; then
            mkdir -p ~{output_directory}/~{lineage}
            touch ~{output_directory}/summary.txt
            touch ~{output_directory}/~{lineage}/full_table.tsv
            exit 0
        fi

        tar -xf ~{lineage_tar}
        compleasm run --library_path ~{lineage_dir} \
                    --assembly_path ~{fasta} \
                    --threads ~{threads} \
                    --lineage ~{lineage} \
                    --output_dir ~{output_directory} \
                    ~{extra_args}
    >>>

    runtime {
        cpu: threads
        memory: "~{memory} GB"
        docker: container
        disk: "local-disk ~{disk_size} HDD"
    }

    output {
        File summary = "~{output_directory}/summary.txt"
        File full_table = "~{output_directory}/~{lineage}/full_table.tsv"
    }
}
