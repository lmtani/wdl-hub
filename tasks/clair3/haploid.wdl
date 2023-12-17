version 1.0


task Clair3Haploid {
    input {
        File alignment
        File alignment_index
        File reference
        File reference_index
        String output_basename
        String model_name
        File model_tar
        Int threads = 8
    }

    command <<<
    set -e
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
    }
}
