version 1.0

import "../../../tasks/compleasm.wdl"
import "../../../tasks/ncbi.wdl"
import "../../../tasks/wget.wdl"


workflow BenchmarkGenomeAssemblies {
    input {
        # NCBI assembly accessions
        Array[String] accessions

        # COMPLEASM parameters
        String compleasm_lineage
        String compleasm_extra_args

        # Stub for testing
        Boolean stub = false
    }

    call compleasm.CompleasmDownload {
        input:
            lineage = compleasm_lineage,
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

        call compleasm.CompleasmRun {
            input:
                fasta = RunWget.downloaded_file,
                lineage_tar = CompleasmDownload.lineage_tar,
                output_directory = basename(fasta_url, ".fna.gz"),
                lineage = compleasm_lineage,
                extra_args = compleasm_extra_args,
                stub = stub
        }
    }

    output {
        Array[File] reports = CompleasmRun.summary
        Array[File] assemblies = CompleasmRun.full_table
    }
}
