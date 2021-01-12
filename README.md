# Achaz-Projet

Ce projet est réalisé par Lou Moizo (IPFB), ainsi que par Margerit William et Ferdinand Petit (BI).

Ce logiciel permet d'effectuer des alignements de reads sur une séquence d'intéret en utilisant la transformées de Burrows-Wheeller. 
Des options pour prendre en compte un certains nombre d'insertions, de délétions et de substitutions pour l'allignements des reads sont disponibles, voir partie Usage. 

## Requirements

ce logiciel a besoins de :
	- Python 3

## Usage : 

Pour lancer le logiciel :

```
python main.py -h
```

Pour les options de lancement :

usage: main.py [-h] [-f [FASTA FILE]] [-r [READ FASTQ FILE]] [-o CSV FILE] [-m M]
               [-i I] [-d D] [-s S] [-l L]

optional arguments:
  -h, --help            show this help message and exit
  -f [FASTA FILE]       fichier fasta du génome analysé
  -r [READ FASTQ FILE]  fichier fastq des séquences read recherchées
  -o CSV FILE           fichier d'output
  -m M                  Recherche read sans mutation, mettre à 0 pour désactiver l'option.
  -i I                  Recherche read avec insertion, mettre à 1 pour activer l'option.
  -d D                  Recherche read avec deletion, mettre à 1 pour activer l'option.
  -s S                  Recherche read avec substitution, mettre à 1 pour activer l'option.
  -l L                  Taille des insertion/deletion/substitution

## Test run

```
python main.py -f Hu-1.fasta -i 1 -d 1 -s 1 -r ../data/read.fastq -o output.csv
```

