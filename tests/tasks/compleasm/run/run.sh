set -e

touch assembly.fasta lineage.tar

miniwdl run --task Run -i tests/tasks/compleasm/run/inputs.json tasks/compleasm/run.wdl
