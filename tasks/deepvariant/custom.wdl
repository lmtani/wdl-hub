version 1.0


task DeepVariantCustomModel {
    input {
        File reference
        File reference_index
        File alignment
        File alignment_index
        String model_type
        File model
        File? regions
        String? make_examples_extra_args

        Boolean stub = false
    }

    String basename = basename(alignment, ".cram")
    String model_name = basename(model, ".tar")

    command <<<
        set -e

        run_deepvariant --version > version.txt

        if [ ~{stub} == "true" ]; then
            touch "~{basename}.deepvariant.vcf.gz" \
                  "~{basename}.deepvariant.vcf.gz.tbi" \
                  "~{basename}.deepvariant.visual_report.html"
            exit 0
        fi

        tar -xf ~{model}

        run_deepvariant \
            --model_type=~{model_type} \
            --customized_model=~{model_name}/model.ckpt \
            --ref=~{reference} \
            --reads=~{alignment} \
            --output_vcf=~{basename}.deepvariant.vcf.gz \
            --num_shards="$(nproc)" \
            ~{"--regions=" + regions} \
            ~{"--make_examples_extra_args=" + make_examples_extra_args} \
            --intermediate_results_dir output/intermediate_results_dir
    >>>

    runtime {
        docker: "google/deepvariant:1.4.0"
    }

    output {
        File output_vcf = "~{basename}.deepvariant.vcf.gz"
        File output_vcf_index = "~{basename}.deepvariant.vcf.gz.tbi"
        File html = "~{basename}.deepvariant.visual_report.html"
        File version = "version.txt"
    }
}
