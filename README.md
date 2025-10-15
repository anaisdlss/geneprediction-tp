# GenePrediction-TP: Prédiction de gènes chez Listeria monocytogenes

**Anaïs DELASSUS** – anais.delassus@etu.u-paris.fr  
Université Paris Cité

## Description du projet
Ce projet permet de **prédire les gènes** dans un génome bactérien à partir de séquences FASTA.  
Les étapes principales sont :

1. Lecture du fichier FASTA du génome (`.fna`)  
2. Recherche des codons de départ (start codons) et d’arrêt (stop codons)  
3. Identification des motifs Shine-Dalgarno en amont des gènes  
4. Filtrage des gènes selon :  
   - Longueur minimale (`min_gene_len`, défaut 50)  
   - Distance maximale Shine-Dalgarno (`max_shine_dalgarno_distance`, défaut 16)  
   - Gap minimal entre deux gènes (`min_gap`, défaut 40)  
5. Prédiction des gènes sur le brin direct (5'→3') et le brin complémentaire (3'→5')  
6. Écriture des positions des gènes dans un fichier CSV (`predict_genes.csv`)  
7. Écriture des séquences des gènes dans un fichier FASTA (`genes.fna`)  

---

## Installation

### Cloner le projet
```bash
git clone https://github.com/anaisdlss/geneprediction-tp.git
cd geneprediction-tp
```
### Créer et activer l'environnement Conda
```bash
conda env create -f environment.yml
conda activate gpred
```

---

## Execution du script de prédiction de gènes
```bash
python3 gpred/gpred.py -i data/listeria.fna
```
Les paramètres optionnels sont :
- `-g` : longueur minimale du gène (par défaut 50)
- `-s` : distance maximale Shine-Dalgarno (par défaut 16)
- `-d` : gap minimal entre deux gènes (par défaut 40)
- `-p` : fichier CSV de sortie des positions des gènes
- `-o` : fichier FASTA de sortie des séquences des gènes

---

## Résultats
- `predict_genes.csv` : positions des gènes prédites par votre script
- `genes.fna` : séquences des gènes prédites
