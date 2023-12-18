version 1.0


task DeepVariant {
    input {
        File reference
        File reference_index
        File alignment
        File alignment_index
        String model_type
        File? regions
        String? make_examples_extra_args
        Boolean stub = false
    }

    String basename = basename(alignment, ".cram")

    command <<<

        run_deepvariant --version > version.txt

        if [ ~{stub} == "true" ]; then
            touch "~{basename}.deepvariant.vcf.gz" \
                  "~{basename}.deepvariant.vcf.gz.tbi" \
                  "~{basename}.deepvariant.visual_report.html"
            exit 0
        fi


        run_deepvariant \
            --model_type=~{model_type} \
            --ref=~{reference} \
            --reads=~{alignment} \
            --output_vcf=~{basename}.deepvariant.vcf.gz \
            --num_shards="$(nproc)" \
            ~{"--regions=" + regions} \
            ~{"--make_examples_extra_args=" + make_examples_extra_args} \
            --intermediate_results_dir output/intermediate_results_dir
    >>>

    runtime {
        docker: "google/deepvariant:1.6.0"
    }

    output {
        File output_vcf = "~{basename}.deepvariant.vcf.gz"
        File output_vcf_index = "~{basename}.deepvariant.vcf.gz.tbi"
        File html = "~{basename}.deepvariant.visual_report.html"
        File version = "version.txt"
    }
}
