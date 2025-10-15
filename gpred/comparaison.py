import pandas as pd

# Charger les fichiers CSV
predicted = pd.read_csv("predict_genes.csv")
prodigal = pd.read_csv("data/prodigal.csv")
experimental = pd.read_csv("data/position.csv")

prodigal.columns = ["Start", "Stop"]
experimental.columns = ["Start", "Stop"]

# Comparer les positions prédictes avec Prodigal
merged_prod = pd.merge(predicted, prodigal,
                       on=["Start", "Stop"],
                       how="outer",
                       indicator=True)
print("Comparaison avec Prodigal:")
print(merged_prod)

# Comparer les positions prédictes avec les positions expérimentales
merged_exp = pd.merge(predicted, experimental,
                      on=["Start", "Stop"], how="outer", indicator=True)
print("Comparaison avec positions expérimentales:")
print(merged_exp)

# Compter combien de gènes sont identiques ou différents
print("\nGènes identiques à Prodigal:",
      (merged_prod["_merge"] == "both").sum())
print("Gènes identiques aux positions expérimentales:",
      (merged_exp["_merge"] == "both").sum())
