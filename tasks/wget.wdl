version 1.0

task RunWget {
    input {
        String file_remote_url
        String extra_args = ""

        String container = "quay.io/biocontainers/wget:1.20.1"
        Int disk_size = 10
        Boolean stub = false
    }

    String file_name = basename(file_remote_url)

    command <<<
        if [ ~{stub} == true ]; then
            touch ~{file_name}
            exit 0
        fi

        wget ~{extra_args} ~{file_remote_url}
    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        docker: container
        disk: "local-disk ~{disk_size} HDD"
    }

    output {
        File downloaded_file = file_name
    }
}
