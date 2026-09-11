"""Réduction dimensionnelle des embeddings SGT : PCA puis t-SNE."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler


def reduce_dimensions(
    embedding: pd.DataFrame,
    n_pca_components: int = 10,
    perplexity: float = 15.0,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Standardise les features SGT, réduit via PCA (dénoise avant t-SNE), puis
    projette en 2D avec t-SNE pour la visualisation interactive.
    """
    n_samples = embedding.shape[0]
    n_pca = min(n_pca_components, n_samples - 1, embedding.shape[1])
    perplexity = min(perplexity, max(5, (n_samples - 1) / 3))

    scaled = StandardScaler().fit_transform(embedding.values)
    pca = PCA(n_components=n_pca, random_state=random_state)
    pca_coords = pca.fit_transform(scaled)

    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        random_state=random_state,
        init="pca",
    )
    tsne_coords = tsne.fit_transform(pca_coords)

    result = pd.DataFrame(
        tsne_coords, columns=["tsne_1", "tsne_2"], index=embedding.index
    )
    result["pca_variance_explained"] = pca.explained_variance_ratio_.sum()
    return result
