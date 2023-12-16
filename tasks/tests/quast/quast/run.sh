set -e
touch assembly.fasta
miniwdl run --task Quast -i tasks/tests/quast/quast/inputs.json tasks/quast.wdl
