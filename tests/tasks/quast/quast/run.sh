set -e
touch assembly.fasta
miniwdl run --task Quast -i tests/tasks/quast/quast/inputs.json tasks/quast/quast.wdl
