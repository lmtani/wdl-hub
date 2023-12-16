version 1.0

import "../../tasks/compleasm.wdl"
import "../../tasks/ncbi.wdl"
import "../../tasks/wget.wdl"
import "../../tasks/quast.wdl"


workflow BenchmarkGenomeAssemblies {
    input {
        # NCBI assembly accessions
        Array[String] accessions

        # QUAST reference - either a file or an accession, but not both
        # In case both are specified, the file will be used
        File? quast_reference_file
        String? quast_reference_accession
        String quast_extra_args = "--glimmer --eukaryote"


        # COMPLEASM parameters
        String compleasm_lineage
        String compleasm_extra_args


        # Stub for testing
        Boolean stub = false
    }

    Boolean available_reference_file = defined(quast_reference_file)

    call compleasm.Download {
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
        call wget.Wget {
            input:
                file_remote_url = fasta_url,
                stub = stub
        }

        call compleasm.Run {
            input:
                fasta = Wget.downloaded_file,
                lineage_tar = Download.lineage_tar,
                output_directory = basename(fasta_url, ".fna.gz"),
                lineage = compleasm_lineage,
                extra_args = compleasm_extra_args,
                stub = stub
        }
    }


    # Download reference if reference_file is not specified
    if (! available_reference_file) {

        String reference_accession = select_first([quast_reference_accession])

        call ncbi.FetchNCBI as FetchReference {
            input:
                accessions = [reference_accession],
        }

        String reference_url = FetchReference.assemblies[0]  # because we only fetch one
        call wget.Wget as RunWgetReference {
            input:
                file_remote_url = reference_url,
                stub = stub
        }
    }

    File reference_sequence = select_first([quast_reference_file, RunWgetReference.downloaded_file])

    call compleasm.Run as reference_compleasm {
        input:
            fasta = reference_sequence,
            lineage_tar = Download.lineage_tar,
            output_directory = "reference",
            lineage = compleasm_lineage,
            extra_args = compleasm_extra_args,
            stub = stub
    }

    call quast.Quast {
        input:
            reference = select_first([quast_reference_file, RunWgetReference.downloaded_file]),
            contigs=Wget.downloaded_file,
            extra_args=quast_extra_args,
            stub = stub
    }

    output {
        Array[File] compleasm_reports = flatten([Run.summary, [reference_compleasm.summary]])
        Array[File] compleasm_full_tables = flatten([Run.full_table, [reference_compleasm.full_table]])
        File quast_report = Quast.output_tar_dir
    }
}
