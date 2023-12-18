set -e

touch truth.vcf.gz query.vcf.gz reference.fasta reference.fasta.fai

miniwdl run --task HAPPY -i tests/tasks/happy/happy/inputs.json tasks/happy/happy.wdl
