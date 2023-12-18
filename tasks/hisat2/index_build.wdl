version 1.0


task HISAT2Build {
    input {
        String reference_name
        File reference_genome

        Boolean stub = false
    }

    command <<<
        set -e

        hisat2-build --version | grep hisat2-build > version.txt

        if [ ~{stub} == "true" ]; then
            touch ~{reference_name}.tar
            exit 0
        fi

        mkdir -p ~{reference_name}
        hisat2-build ~{reference_genome} ~{reference_name}/~{reference_name}

        tar -cf ~{reference_name}.tar ~{reference_name}/
    >>>

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
        memory: "16 GB"
        disks: "local-disk 10 HDD"
        preemptible: 5
    }

    output {
        File hisat2_index = "~{reference_name}.tar"
        File version = "version.txt"
    }
}
