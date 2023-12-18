version 1.0


task MakeMinimap2Index {
    input {
        File reference_genome
        Boolean stub = false
    }

    String basename = basename(reference_genome, ".fasta")

    command <<<

        minimap2 --version > version.txt

        if [ ~{stub} == "true" ]; then
            touch ~{basename}.mmi
            exit 0
        fi

        minimap2 -d ~{basename}.mmi ~{reference_genome}
    >>>

    runtime {
        docker: "quay.io/biocontainers/minimap2:2.26--he4a0461_1"
    }

    output {
        File index = "~{basename}.mmi"
        File version = "version.txt"
    }
}
