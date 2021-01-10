# Calcul des OTU

## Utilisation

* Pour générer une nouvelle liste d'OTU:

```sh
python agc/agc.py -i data/amplicon.fasta.gz -o data/otu_result.fna
```

* Pour tester la qualité du résultat fourni à l’aide de vsearch et des références 16S de nos séquences d’entrée: mock_16S.fasta.

```sh
vsearch --usearch_global data/otu_result.fna --db data/mock_16S.fasta --id 0.9 --alnout quality.txt
```

Le compte-rendu du test de qualité se trouve dans le fichier `quality.txt` à la racine du projet.

*Note: Notre programme trouve 3 OTU avec `OTU_3` très similaire à `OTU_2` (94.1% de similarité).*

## Tests unitaires

Cette version du programme de calcul des OTU ajoute deux fonctions supplémentaires `compute_id_matrix` et `detect_chimera` avec des tests associés.

Par ailleurs, la fonction `write_OTU` ne valide pas les tests car la fonction de hashage ne donne pas le résultat attendu malgrès le fait que le fichier contenant les OTU soit dans le format attendu:

```txt
>OTU_{numéro partant de 1} occurrence:{nombre d’occurrence à la déréplication}
{séquence au format fasta}
```
