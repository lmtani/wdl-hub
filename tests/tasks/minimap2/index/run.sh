touch reference.fasta

miniwdl run --task MakeMinimap2Index -i tests/tasks/minimap2/index/inputs.json tasks/minimap2/index.wdl
