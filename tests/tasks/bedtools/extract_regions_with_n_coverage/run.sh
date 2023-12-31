touch sample.per-base.bed.gz

miniwdl run --task ExtractRegionsWithNCoverage -i tests/tasks/bedtools/extract_regions_with_n_coverage/inputs.json tasks/bedtools/extract_regions_with_n_coverage.wdl
