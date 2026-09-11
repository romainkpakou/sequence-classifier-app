"""
Vectorisation de séquences sans alignement préalable, via Sequence Graph
Transform (SGT). Chaque séquence protéique est transformée en un vecteur de
taille fixe qui capture ses motifs de proximité entre acides aminés (courte
et longue portée), sans nécessiter d'alignement de référence.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sgt import SGT


def sequences_to_sgt_matrix(sequences: dict[str, str], kappa: float = 5.0) -> pd.DataFrame:
    """
    Convertit un dict {id: séquence} en une matrice de features SGT
    (une ligne par séquence).

    kappa : paramètre de portée du SGT (plus il est grand, plus les
    interactions longue-distance entre positions pèsent dans l'embedding).
    """
    corpus = pd.DataFrame({
        "id": list(sequences.keys()),
        "sequence": [list(seq) for seq in sequences.values()],
    })

    sgt = SGT(kappa=kappa, lengthsensitive=False)
    embedding = sgt.fit_transform(corpus)
    embedding = embedding.set_index("id")
    return embedding


def read_fasta(path_or_buffer) -> dict[str, str]:
    """Parseur FASTA minimal, sans dépendance à Biopython pour cette étape."""
    sequences: dict[str, str] = {}
    current_id = None
    current_seq: list[str] = []

    lines = path_or_buffer.read().decode("utf-8").splitlines() \
        if hasattr(path_or_buffer, "read") else open(path_or_buffer).read().splitlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                sequences[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences
