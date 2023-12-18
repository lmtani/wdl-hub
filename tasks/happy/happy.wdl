version 1.0

task HAPPY {
    input {
        File truth_vcf
        File query_vcf
        String prefix
        File fasta
        File fasta_fai
        File? regions_bed
        File? targets_bed
        File? false_positives_bed
        File? stratification_tsv

        Int cpus
        Boolean stub = false
    }


    command <<<
        hap.py -v > version.txt  # It's not working as expected. Just prints Hap.py

        if [ ~{stub} == "true" ]; then
            touch ~{prefix}.summary.csv \
                  ~{prefix}.roc.all.csv.gz \
                  ~{prefix}.roc.Locations.INDEL.csv.gz \
                  ~{prefix}.roc.Locations.INDEL.PASS.csv.gz \
                  ~{prefix}.roc.Locations.SNP.csv.gz \
                  ~{prefix}.roc.Locations.SNP.PASS.csv.gz \
                  ~{prefix}.extended.csv \
                  ~{prefix}.runinfo.json \
                  ~{prefix}.metrics.json.gz \
                  ~{prefix}.vcf.gz \
                  ~{prefix}.vcf.gz.tbi
            exit 0
        fi

        hap.py \
            ~{truth_vcf} \
            ~{query_vcf} \
            ~{"--reference " + fasta} \
            ~{"--threads " + cpus} \
            ~{"--R " + regions_bed } \
            ~{"--T " + targets_bed } \
            ~{"--false-positives " + false_positives_bed} \
            ~{"--stratification " + stratification_tsv} \
            -o ~{prefix}
    >>>

    output {
        File summary_csv = "~{prefix}.summary.csv"
        File roc_all_csv = "~{prefix}.roc.all.csv.gz"
        File roc_indel_locations_csv = "~{prefix}.roc.Locations.INDEL.csv.gz"
        File roc_indel_locations_pass_csv = "~{prefix}.roc.Locations.INDEL.PASS.csv.gz"
        File roc_snp_locations_csv = "~{prefix}.roc.Locations.SNP.csv.gz"
        File roc_snp_locations_pass_csv = "~{prefix}.roc.Locations.SNP.PASS.csv.gz"
        File extended_csv = "~{prefix}.extended.csv"
        File runinfo = "~{prefix}.runinfo.json"
        File metrics_json = "~{prefix}.metrics.json.gz"
        File vcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
        File version = "version.txt"
    }

    runtime {
        docker: "quay.io/biocontainers/hap.py:0.3.14--py27h5c5a3ab_0"
        cpu: cpus
    }
}
