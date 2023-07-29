version 1.0


task MinibuscoDownload {
    input {
        String lineage  # BUSCO compatible lineage name
        String lineage_dir = "mb_downloads"
        Boolean stub = false
    }

    command <<<
        if [ ~{stub} == "true" ]; then
            touch ~{lineage}.tar.gz
            exit 0
        fi

        minibusco download --library_path ~{lineage_dir} ~{lineage}
        tar -czvf ~{lineage}.tar.gz ~{lineage_dir}
    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        docker: "quay.io/biocontainers/minibusco:0.2.1--pyh7cba7a3_0"
        disk: "local-disk 10 HDD"
    }

    output {
        File lineage_tar = "~{lineage}.tar.gz"
    }
}


task MinibuscoRun {
    input {
        File fasta
        File lineage_tar
        String lineage
        String output_directory
        Int threads = 4
        Int memory = 12
        String lineage_dir = "mb_downloads"
        Boolean stub = false
    }

    command <<<
        set -e
        if [ ~{stub} == "true" ]; then
            mkdir -p ~{output_directory} && touch ~{output_directory}/summary.txt
            exit 0
        fi

        tar -xf ~{lineage_tar}
        minibusco run --library_path ~{lineage_dir} -a ~{fasta} -t ~{threads} -l ~{lineage} -o ~{output_directory}
    >>>

    runtime {
        cpu: threads
        memory: memory
        docker: "quay.io/biocontainers/minibusco:0.2.1--pyh7cba7a3_0"
        disk: "local-disk 10 HDD"
    }

    output {
        File report = "~{output_directory}/summary.txt"
    }
}
