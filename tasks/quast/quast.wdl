version 1.0

task Quast {
    input {
        Array[File] contigs
        Int threads = 4
        String extra_args = ""
        String output_dir = "./quast_output"
        File? reference
        String container = "quay.io/biocontainers/quast:5.2.0--py39pl5321h4e691d4_3"
        Boolean stub = false
    }

    Int disk_size = ceil(size(reference, "GiB") + size(contigs, "GiB")) + 10

    command <<<
        if [ ~{stub} == "true" ]; then
            mkdir -p "~{output_dir}"
            touch "~{output_dir}.tar.gz"
            touch "~{output_dir}/report.tsv"
            exit 0
        fi

        quast.py  ~{extra_args} \
            --threads ~{threads} \
            --output-dir ~{output_dir} \
            ~{"-r " + reference}  \
            ~{sep=' ' contigs}

        tar -czvf ~{output_dir}.tar.gz ~{output_dir}
    >>>

    runtime {
        docker: container
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: threads
    }

    output {
        File output_tar_dir = "~{output_dir}.tar.gz"
        File report_tsv = "~{output_dir}/report.tsv"
    }
}
