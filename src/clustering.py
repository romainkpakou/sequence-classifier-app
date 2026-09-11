"""Clustering des séquences à partir de leur embedding SGT/PCA."""
from __future__ import annotations

import pandas as pd
from sklearn.cluster import KMeans


def cluster_sequences(embedding: pd.DataFrame, n_clusters: int = 4, random_state: int = 42) -> pd.Series:
    """
    Clustering K-means sur les features SGT standardisées (avant réduction
    t-SNE, pour clusterer sur l'information complète plutôt que sur la seule
    projection 2D).
    """
    from sklearn.preprocessing import StandardScaler

    n_clusters = min(n_clusters, embedding.shape[0])
    scaled = StandardScaler().fit_transform(embedding.values)
    km = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
    labels = km.fit_predict(scaled)
    return pd.Series(labels, index=embedding.index, name="cluster")


def cluster_medoid(sequences: dict[str, str], cluster_ids: list[str]) -> str:
    """
    Retourne l'ID de la séquence la plus "centrale" d'un cluster (distance
    d'édition totale minimale aux autres membres) -- sert de référence pour
    l'alignement guidé utilisé par le logo de séquence.
    """
    import difflib

    if len(cluster_ids) == 1:
        return cluster_ids[0]

    best_id, best_score = None, float("inf")
    for candidate in cluster_ids:
        total = sum(
            1 - difflib.SequenceMatcher(None, sequences[candidate], sequences[other]).ratio()
            for other in cluster_ids if other != candidate
        )
        if total < best_score:
            best_score, best_id = total, candidate
    return best_id
