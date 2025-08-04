Coverage analysis of the transcriptome and genomes of the different species of CWR of wheat. In order to determine if the data is usable for detecting selection.

Use Bams to analyse the coverage.
index the Bams
```bash
module load bioinfo-cirad
module load samtools/1.14-bin

for bam in *.bam; do
    samtools index "$bam"
done
```
