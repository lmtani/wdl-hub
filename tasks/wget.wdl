version 1.0

task RunWget {
    input {
        String file_remote_url
        Boolean stub = false
    }

    String file_name = basename(file_remote_url)

    command <<<
        if [ ~{stub} == true ]; then
            touch ~{file_name}
            exit 0
        fi

        wget ~{file_remote_url}
    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        docker: "quay.io/biocontainers/wget:1.20.1"
        disk: "local-disk 10 HDD"
    }

    output {
        File downloaded_file = file_name
    }
}
