version 1.0


task FetchNCBI {
    input {
        Array[String] accessions
        String your_email = "your.email@domain.com"

        String container = "quay.io/biocontainers/biopython:1.75"
        Int disk_size = 10
        Boolean stub = false
    }

    command <<<
        if [ ~{stub} == "true" ]; then
            echo "a-url-addres-for-ftp-download"
            exit 0
        fi

        python <<CODE
        from Bio import Entrez

        # Always tell NCBI who you are
        Entrez.email = "~{your_email}"


        def get_ftp_url(accession):
            # Fetch the assembly summary
            handle = Entrez.esearch(db="assembly", term=accession)
            id_list = Entrez.read(handle)["IdList"]
            ids = ",".join(id_list)
            ids = id_list[0]
            handle = Entrez.esummary(db="assembly", id=ids)
            record = Entrez.read(handle)
            return record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_Stats_rpt"].replace("_assembly_stats.txt", "_genomic.fna.gz")


        # read lines from accessions.txt (example: "GCA_949128135.1")
        with open("~{write_lines(accessions)}") as f:
            accessions = f.readlines()

        for accession in accessions:
            print(get_ftp_url(accession))
        CODE
    >>>

    runtime {
        cpu: 1
        memory: "2 GB"
        disk: "local-disk ~{disk_size} HDD"
        docker: container
    }

    output {
        Array[String] assemblies = read_lines(stdout())
    }
}
