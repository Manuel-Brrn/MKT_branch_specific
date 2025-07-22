# Branch-Specific Corrected MK Test

## Overview  
The McDonald-Kreitman Test (MKT) analyzes:  
- **Input**: Alignment of focal species sequences + two or more outgroup (for branch-specific substitutions).  

---

## Divergence Analysis  

### **1. Quality Control & Orthology**  
- **VESPA Pipeline**:  
  - Validates CDS sequences (length must be divisible by 3).  
  - Retains only protein-coding sequences.  
- **Outgroup Strategy**:  
  - Two outgroups per focal species:  
    1. **Closely related sister species**.  
    2. **Divergent wild barley** (*Hordeum spontaneum*).  
- **Ortholog Identification**:  
  - **Reciprocal Best BLAST Hit (RBB)** with thresholds:  
    - E-value < 0.001  
    - Alignment length ≥ 75 amino acids  
    - Identity ≥ 75%  

### **2. Alignment & Phylogeny**  
- **Protein MSAs**:  
  - Generated with `MUSCLE`.  
  - Cleaned with `HMMcleaner` to remove alignment errors.  
- **Evolutionary Model Selection**:  
  - Best-fit model chosen via `ProtTest3` (integrated in VESPA).  
- **Gene Trees**:  
  - Reconstructed with `MrBayes` (200,000 MCMC generations, 25% burn-in).  

### **3. Selection Analysis**  
- **Codon Alignments**:  
  - Protein-guided alignments generated with `map_alignments`.  
- **Branch-Site Models**:  
  - Run via `codeML` (PAML) to estimate **dN/dS** for the focal branch.  

---

## Polymorphism Data  

### **1. Sequence Processing**  
- **Reference Transcriptomes**:  
From Ensembl BioMart or generated from genome assembly using `AGAT`.  
- **Read Mapping**:  
  - Performed with `Gecko` workflows.  
  - **Filtering**: Removed soft-clipped reads (potential structural variants).  
- **Coverage Criteria**:  
  - Sites retained if ≥10 reads in ≥7 individuals.  
  - Genes classified as "high quality" if ≥50% positions met coverage.  

### **2. SNP Calling & Alignment**  
- **Consensus Sequences**:  
  - Generated with `Reads2snp_2.0` (min frequency threshold applied).  
- **CDS Extraction**:  
  - Used `ORF_extractor.pl` (PopPhyl) + GFF annotations.  
- **Codon-Aware Alignment**:  
  - Aligned with `MACSE`, cleaned with `HMMcleaner`.  
- **Polymorphism Filtering**:  
  - Excluded genes with <5 polymorphic sites or poor alignment quality.  

---

## Key Tools  
| **Step**               | **Tool**          |  
|------------------------|-------------------|  
| Orthology              | BLAST (RBB)       |  
| Protein Alignment      | MUSCLE, MACSE     |  
| Alignment Cleaning     | HMMcleaner        |  
| Phylogenetics          | MrBayes           |  
| Selection Analysis     | codeML (PAML)     |  
| SNP Calling            | Reads2snp_2.0      |  

For details, refer to [VESPA pipeline](https://example.com) and [PopPhyl](http://kimura.univ-montp2.fr/PopPhyl/).  
