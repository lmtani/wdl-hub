version 1.0

# TODO: maybe modify it to make generic for any model
task Clair3Haploid {
    input {
        File alignment
        File alignment_index
        File reference
        File reference_index
        String model_name
        File model_tar
        String output_basename

        Int threads = 8
        Boolean stub = false
    }

    command <<<
    set -e

    /opt/bin/run_clair3.sh --version > version.txt

    if [ ~{stub} == "true" ]; then
        mkdir -p ~{output_basename}
        touch "~{output_basename}/full_alignment.vcf.gz" \
                "~{output_basename}/full_alignment.vcf.gz.tbi" \
                "~{output_basename}/pileup.vcf.gz" \
                "~{output_basename}/pileup.vcf.gz.tbi" \
                "~{output_basename}.merged.vcf.gz" \
                "~{output_basename}.merged.vcf.gz.tbi"
        exit 0
    fi

    tar -xf ~{model_tar} -C ./

    /opt/bin/run_clair3.sh \
        --bam_fn=~{alignment} \
        --ref_fn=~{reference} \
        --threads=~{threads} \
        --platform="ont" \
        --model_path="./~{model_name}" \
        --output=~{output_basename} \
        --no_phasing_for_fa \
        --include_all_ctgs \
        --haploid_precise \
        --sample_name=~{output_basename}

    cp ~{output_basename}/merge_output.vcf.gz ./~{output_basename}.merged.vcf.gz
    cp ~{output_basename}/merge_output.vcf.gz.tbi ./~{output_basename}.merged.vcf.gz.tbi
    >>>

    runtime {
        docker: "docker.io/hkubal/clair3:v1.0.4"
    }

    output {
        File full_alignment = "~{output_basename}/full_alignment.vcf.gz"
        File full_alignment_index = "~{output_basename}/full_alignment.vcf.gz.tbi"
        File pileup = "~{output_basename}/pileup.vcf.gz"
        File pileup_index = "~{output_basename}/pileup.vcf.gz.tbi"
        File vcf = "~{output_basename}.merged.vcf.gz"
        File vcf_index = "~{output_basename}.merged.vcf.gz.tbi"
        File version = "version.txt"
    }
}
