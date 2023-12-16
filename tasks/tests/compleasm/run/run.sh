set -e

touch assembly.fasta lineage.tar

miniwdl run --task Run -i tasks/tests/compleasm/run/inputs.json tasks/compleasm.wdl
