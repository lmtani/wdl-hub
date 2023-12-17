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
