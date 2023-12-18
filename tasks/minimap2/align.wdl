version 1.0


task Minimap2Align {
    input {
        File reads_tar
        File reference_mmi
        Int cpus = 8

        String sample_name
        String sample_id
        String technology
        String library

        Boolean stub = false
    }

    Int preemptible_tries = 3
    Int memory = cpus * 2
    Int disk_size = ceil(size(reads_tar, "GiB") * 8 + size(reference_mmi, "GiB") + 10)

    String prefix = basename(basename(reads_tar, ".gz"), ".tar")

    command <<<

        minimap2 --version > version.txt

        if [ ~{stub} == "true" ]; then
            touch ~{prefix}.sam
            exit 0
        fi

        mkdir -p reads/
        tar -xf ~{reads_tar} -C reads/
        rm ~{reads_tar}

        minimap2 -R "@RG\tID:~{sample_id}\tSM:~{sample_name}\tPL:~{technology}\tLB:~{library}" \
            -t ~{cpus} -a -Y -o ~{prefix}.sam -x map-ont ~{reference_mmi} reads/*
    >>>

    runtime {
        docker: "quay.io/biocontainers/minimap2:2.26--he4a0461_1"
        disks: "local-disk ~{disk_size} HDD"
        cpu: "~{cpus}"
        memory: "~{memory} GB"
        preemptible: preemptible_tries
    }

    output {
        File alignment = "~{prefix}.sam"
        File version = "version.txt"
    }
}
