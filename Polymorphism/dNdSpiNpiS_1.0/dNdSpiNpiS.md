**dNdSpiNpiS Analysis Pipeline**
Calculates dN/dS (non-synonymous/synonymous substitution rates), piN/piS (non-synonymous/synonymous polymorphism
ratios), and evolutionary distances between species pairs.

## Requirements
- SLURM workload manager
- Compiled `dNdSpiNpiS_1.0` executable in $PATH
- FASTA alignment file

## Usage
```bash
sbatch dNdSpiNpiS.sbatch <alignment_file> <output_dir> <ingroup> <outgroup>
```
Example:
```bash
sbatch dNdSpiNpiS.sbatch \
    /data/alignments/merged_cat.fasta \
    /results/dNdS_analysis \
    Ae_speltoides \
    Vu_hordeum
```
