# GenePrediction-TP: Gene Prediction in Listeria monocytogenes

**Anaïs DELASSUS** – anais.delassus@etu.u-paris.fr  
Université Paris Cité

## Description du projet
This project implements a prokaryotic gene prediction pipeline for the genome of Listeria monocytogenes using only raw FASTA sequences.
Les étapes principales sont :

The main steps of the prediction workflow are:
	
1.	Reading the bacterial genome from a FASTA file (.fna)
   
2.	Searching for start codons and stop codons
	
3.	Detecting Shine–Dalgarno motifs upstream of candidate genes
	
4.	Filtering genes based on:
- Minimum gene length (min_gene_len, default: 50 bp)
- Maximum Shine–Dalgarno distance (max_shine_dalgarno_distance, default: 16 bp)
- Minimum gap between consecutive genes (min_gap, default: 40 bp)
	
5.	Predicting genes on both DNA strands:
- Forward strand (5’ → 3’)
- Reverse-complement strand (3’ → 5’)
	
6.	Writing predicted gene positions to a CSV file (predict_genes.csv)
	
7.	Writing predicted gene sequences to a FASTA file (genes.fna)
   
---

## Installation

### Clone the repository
```bash
git clone https://github.com/anaisdlss/geneprediction-tp.git
cd geneprediction-tp
```
### Create and activate the Conda environment
```bash
conda env create -f environment.yml
conda activate gpred
```

---

## Running the Gene Prediction Script
```bash
python3 gpred/gpred.py -i data/listeria.fna
```
Optional parameters:
- `-g` : Minimum gene length (défaut 50)
- `-s` : Maximum Shine–Dalgarno distance (défaut 16)
- `-d` : Minimum gap between two genes (défaut 40)
- `-p` : Output CSV file for gene positions (predict_genes.csv)
- `-o` : Output FASTA file for gene sequences (genes.fna)


---

## Output Files

- `predict_genes.csv` : List of all predicted gene coordinates
- `genes.fna` : FASTA sequences of all predicted genes

---

## Comparing Predicted Genes with Reference Data

After generating your predicted gene list using gpred.py, you can evaluate the accuracy by comparing:
- your predictions (predict_genes.csv)
- Prodigal predictions (prodigal.csv)
- experimentally validated positions (position.csv)

### Run the comparison script

```bash
python3 gpred/comparaison.py
```
