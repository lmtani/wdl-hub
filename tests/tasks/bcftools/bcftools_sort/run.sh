touch variants.vcf.gz

miniwdl run --task Sort -i tests/tasks/bcftools/bcftools_sort/inputs.json tasks/bcftools/bcftools_sort.wdl
