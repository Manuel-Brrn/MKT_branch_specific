# Vespa

*Repository of scripts to run the VESPA pipeline, in order to run codeml*

##  **Goal**  
- Run codeml easily (selection analysis)
- Clean data
- Obtain orthologs
- Phylogeny reconstruction
- Obtain dn/ds on a specific branch of a phylogeny
- Avoid outgroup biaise

##  **Installation**  
```bash
wget https://github.com/aewebb80/VESPA/archive/refs/tags/1.0.1.tar.gz
tar -xzf 1.0.1.tar.gz
cd VESPA
chmod +x vespa.py
mv vespa.py /usr/local/bin
```

## **Dependencies**  
```bash
module load bioinfo-cirad
module load python/2.7.18
```

----------------------------------------------------------------------------------
|||||                        **Phase 1: Data Preparation**                       |||||
----------------------------------------------------------------------------------
This phase included for users new to bioinformatics. The phase prepares downloaded
genomes for homology searching using the two VESPA supported homology search tools:
BLAST and HMMER.

Commands: clean, ensemble_clean, translate, create_database, gene_selection
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||                       Phase 2: Homology Searching                      |||||
----------------------------------------------------------------------------------
Details: This phase is concerned with identifying groups of similar sequences from
either BLAST or HMMER homology searches.

Commands: similarity_groups, reciprocal_groups, best_reciprocal_groups
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||        Phase 3: Alignment Assessment & Phylogeny Reconstruction        |||||
----------------------------------------------------------------------------------
Details: This phase combines multiple third-party programs to automate the
assessment, selection, and phylogenetic reconstruction of protein MSAs.

Commands: metal_compare, prottest_setup, prottest_reader, mrbayes_setup
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||                 Phase 4: Selection Analysis Preparation                |||||
----------------------------------------------------------------------------------
Details: This phase automates large-scale selective pressure analysis using codeML
from the PAML package.

Commands: map_alignments, infer_genetree, mrbayes_reader, link_input, codeml_setup,
          create_subtrees, create_branch
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||                 Phase 5: Selection Analysis Assessment                 |||||
----------------------------------------------------------------------------------
Details: This phase automatically parses the codeML directory structure
and create simplified summary files

Command: codeml_reader
----------------------------------------------------------------------------------
