# Plan — Application de classification de séquences sans alignement (démo publique)

## Contexte

Reprise de l'application Streamlit développée en stage chez Affilogic pour la classification
de Nanofitines (ligands d'affinité) sans alignement préalable. **Les séquences de Nanofitines
sont la propriété d'Affilogic** — à vérifier auprès de l'entreprise (contrat de stage /
accord de confidentialité) avant toute publication. Ce dépôt démontre la même méthodologie
sur des séquences protéiques publiques du même type fonctionnel (petits domaines de liaison).

## Objectif

Classifier et regrouper des séquences protéiques sans alignement de référence, via
vectorisation par graphe de séquence (SGT), réduction dimensionnelle et visualisation
interactive, dans une application web utilisable par un non-bioinformaticien.

## Données

| Source | Description | Accès |
|---|---|---|
| SAbDab (Structural Antibody Database) | Séquences de nanobodies/domaines VHH, publiques | `https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/` |
| UniProt | Familles de domaines de liaison, alternative/complément | `https://www.uniprot.org` |

## Méthodologie

1. Récupération d'un jeu de séquences publiques (FASTA) de domaines de liaison
2. Vectorisation par **Sequence Graph Transform (SGT)** — sans alignement préalable
3. Réduction dimensionnelle (PCA, t-SNE)
4. Clustering et visualisation interactive 2D/3D (Plotly)
5. Génération de logos de séquence par cluster (Logomaker) pour identifier positions
   conservées/variables
6. Interface Streamlit : upload FASTA, paramètres ajustables (perplexité t-SNE, nombre de
   clusters), export des résultats

## Stack technique

`Python` · `Streamlit` · `scikit-learn` · `sgt` · `Plotly` · `Logomaker` · `Docker`

## Structure de dépôt prévue

```
sequence-classifier-app/
├── README.md                # précise la substitution Nanofitines → séquences publiques
├── app.py
├── src/
│   ├── vectorization.py
│   ├── dimensionality_reduction.py
│   └── clustering.py
├── data/
│   └── example_sequences.fasta
├── Dockerfile
└── requirements.txt
```

## Livrables

- Application Streamlit fonctionnelle, conteneurisée (Docker)
- Démo publique en ligne si possible (Streamlit Community Cloud, gratuit)
- README expliquant clairement la substitution de données et le contexte du stage d'origine

## Écart au plan initial : globines plutôt que SAbDab

Le plan prévoyait des séquences de nanobodies (SAbDab/UniProt). En pratique, UniProt ne
catalogue quasiment pas de nanobodies nommément (ce sont des constructions issues de
recherche/brevets, pas des entrées naturelles curées). Le dépôt utilise donc des séquences
de **globines humaines** (hémoglobines, myoglobine — UniProt/Swiss-Prot, `reviewed:true`),
facilement récupérables via l'API REST UniProt et suffisamment diverses fonctionnellement
pour démontrer la méthode (vectorisation, clustering, logos). Le code n'a pas été
réutilisé du stage (écrit intégralement pour ce dépôt) : la question de confidentialité
Affilogic ne s'applique donc qu'aux données, pas au code, et ne bloque pas la publication.

## Statut

- [x] Récupération du jeu de séquences publiques (100 globines, UniProt/Swiss-Prot)
- [x] Implémentation de l'app (vectorisation SGT, PCA/t-SNE, K-means, logos)
- [x] Test du pipeline en local (hors Docker) — fonctionnel
- [ ] Déploiement démo publique (Streamlit Community Cloud)
