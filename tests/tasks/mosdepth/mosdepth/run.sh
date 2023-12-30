touch reference.fasta reference.fasta.fai alignment.cram alignment.cram.crai regions.bed

miniwdl run --task RunMosDepth -i $1 tasks/mosdepth/mosdepth.wdl
