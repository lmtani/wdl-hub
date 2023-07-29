version 1.0

import "../tasks/compleasm.wdl"
import "../tasks/ncbi.wdl"
import "../tasks/wget.wdl"


workflow MinibuscoWorkflow {
    input {
        Array[String] accessions
        String busco_lineage
        Boolean stub = false
    }

    call compleasm.MinibuscoDownload {
        input:
            lineage = busco_lineage,
            stub = stub
    }

    call ncbi.FetchNCBI {
        input:
            accessions = accessions,
            stub = stub
    }

    scatter (fasta_url in FetchNCBI.assemblies) {
        call wget.RunWget {
            input:
                file_remote_url = fasta_url,
                stub = stub
        }

        call compleasm.MinibuscoRun {
            input:
                fasta = RunWget.downloaded_file,
                lineage_tar = MinibuscoDownload.lineage_tar,
                output_directory = basename(fasta_url, ".fna.gz"),
                lineage = busco_lineage,
                stub = stub
        }
    }

    output {
        Array[File] reports = MinibuscoRun.report
    }
}
