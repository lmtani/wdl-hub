version 1.0

import "../../../tasks/compleasm.wdl"
import "../../../tasks/ncbi.wdl"
import "../../../tasks/wget.wdl"
import "../../../tasks/quast.wdl"


workflow BenchmarkGenomeAssemblies {
    input {
        # NCBI assembly accessions
        Array[String] accessions

        # COMPLEASM parameters
        String compleasm_lineage
        String compleasm_extra_args

        # QUAST reference - either a file or an accession, but not both
        File? quast_reference_file
        String? quast_reference_accession

        # Stub for testing
        Boolean stub = false
    }

    Boolean download_reference = defined(quast_reference_accession)

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


    # Download reference if necessary
    if (download_reference) {

        String reference_accession = select_first([quast_reference_accession])

        call ncbi.FetchNCBI as FetchReference {
            input:
                accessions = [reference_accession],
        }

        String reference_url = FetchReference.assemblies[0]
        call wget.RunWget as RunWgetReference {
            input:
                file_remote_url = reference_url,
                stub = stub
        }
    }

    call quast.RunQuast {
        input:
            reference = select_first([quast_reference_file, RunWgetReference.downloaded_file]),
            contigs=RunWget.downloaded_file,
            extra_args="--glimmer --eukaryote",
            stub = stub
    }

    output {
        Array[File] reports = CompleasmRun.summary
        Array[File] assemblies = CompleasmRun.full_table
        File quast_report = RunQuast.output_tar_dir
    }
}
