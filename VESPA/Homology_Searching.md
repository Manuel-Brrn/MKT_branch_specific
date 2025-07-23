----------------------------------------------------------------------------------
|||||                       Phase 2: Homology Searching                      |||||
----------------------------------------------------------------------------------

**Orthofinder:** 
OrthoFinder is a fast, accurate and comprehensive platform for comparative genomics. It finds orthogroups and orthologs, infers rooted gene trees for all orthogroups and identifies all of the gene duplication events in those gene trees. It also infers a rooted species tree for the species being analysed and maps the gene duplication events from the gene trees to branches in the species tree. OrthoFinder also provides comprehensive statistics for comparative genomic analyses. 

```bash
mkdir orthofinder_input
```
**Installation:** 
```bash
module load bioinfo-cirad
module load python/3.9.13
wget https://github.com/davidemms/OrthoFinder/archive/refs/tags/3.0.1b1.tar.gz
tar xzf 3.0.1b1.tar.gz
pip install --upgrade pip
pip install scipy
pip install scikit-learn
pip install bio
python3
```

**Run Orthofinder:** 
The directory must point where the species fasta files of proteins sequences are present
```bash
python3 orthofinder.py -f /home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/dn_ds_pipeline/VESPA/cds_sequences/database/orthofinder_input
```

**Orthologues:**
Present in orthofinder_input/OrthoFinder/Results_Jul23/Single_Copy_Orthologue_Sequences
