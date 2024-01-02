touch variants.vcf.gz

miniwdl run --task Norm -i tests/tasks/bcftools/bcftools_norm/inputs.json tasks/bcftools/bcftools_norm.wdl
