"""
Application de classification de séquences protéiques sans alignement
préalable. Démonstration publique de la méthodologie développée en stage
(Affilogic, classification de Nanofitines) -- ici appliquée à des séquences
publiques de globines (UniProt/Swiss-Prot), les séquences de Nanofitines
étant propriétaires. Voir PLAN.md et README.md pour le détail de cette
substitution.
"""
from __future__ import annotations

import io

import logomaker
import pandas as pd
import plotly.express as px
import streamlit as st
from Bio.Align import PairwiseAligner, substitution_matrices

from src.clustering import cluster_medoid, cluster_sequences
from src.dimensionality_reduction import reduce_dimensions
from src.vectorization import read_fasta, sequences_to_sgt_matrix

st.set_page_config(page_title="Classification de séquences sans alignement", layout="wide")

st.title("Classification de séquences protéiques sans alignement")
st.caption(
    "Vectorisation SGT (Sequence Graph Transform) → PCA → t-SNE → clustering. "
    "Démo sur séquences publiques de globines (UniProt) -- méthodologie développée "
    "initialement en stage sur des séquences propriétaires (Nanofitines, Affilogic)."
)

with st.sidebar:
    st.header("Données")
    uploaded = st.file_uploader("Fichier FASTA", type=["fasta", "fa", "txt"])
    use_example = st.checkbox("Utiliser le jeu d'exemple (globines, UniProt)", value=uploaded is None)

    st.header("Paramètres")
    kappa = st.slider("Kappa (portée SGT)", 1.0, 10.0, 5.0, 0.5)
    n_clusters = st.slider("Nombre de clusters (K-means)", 2, 10, 4)
    perplexity = st.slider("Perplexité t-SNE", 5, 50, 15)

if use_example:
    sequences = read_fasta("data/example_sequences.fasta")
elif uploaded is not None:
    sequences = read_fasta(uploaded)
else:
    st.info("Charge un fichier FASTA ou coche la case pour utiliser le jeu d'exemple.")
    st.stop()

if len(sequences) < 4:
    st.error("Il faut au moins 4 séquences pour une classification informative.")
    st.stop()

st.write(f"**{len(sequences)} séquences** chargées.")

with st.spinner("Vectorisation SGT en cours..."):
    embedding = sequences_to_sgt_matrix(sequences, kappa=kappa)

with st.spinner("Réduction dimensionnelle (PCA → t-SNE)..."):
    coords = reduce_dimensions(embedding, perplexity=perplexity)

labels = cluster_sequences(embedding, n_clusters=n_clusters)

plot_df = coords.join(labels).reset_index()
plot_df.columns = ["id", "tsne_1", "tsne_2", "pca_variance_explained", "cluster"]
plot_df["cluster"] = plot_df["cluster"].astype(str)
plot_df["length"] = plot_df["id"].map(lambda i: len(sequences[i]))

col1, col2 = st.columns([2, 1])

with col1:
    fig = px.scatter(
        plot_df, x="tsne_1", y="tsne_2", color="cluster", hover_name="id",
        hover_data={"length": True, "tsne_1": False, "tsne_2": False},
        title="Projection t-SNE des séquences (colorées par cluster)",
    )
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color="white")))
    st.plotly_chart(fig, use_container_width=True)

with col2:
    st.metric("Variance expliquée (PCA)", f"{plot_df['pca_variance_explained'].iloc[0]:.1%}")
    st.dataframe(
        plot_df[["id", "cluster", "length"]].sort_values("cluster"),
        use_container_width=True, height=400,
    )

st.divider()
st.subheader("Logo de séquence par cluster")
st.caption(
    "Alignement guidé par référence (séquence médoïde du cluster) via alignement "
    "global Biopython/BLOSUM62 -- une approximation légère, pas une MSA rigoureuse "
    "(type MUSCLE/Clustal), suffisante pour visualiser les positions conservées au "
    "sein d'un cluster de séquences homologues."
)

selected_cluster = st.selectbox("Cluster à visualiser", sorted(plot_df["cluster"].unique()))
cluster_ids = plot_df.loc[plot_df["cluster"] == selected_cluster, "id"].tolist()

if len(cluster_ids) < 2:
    st.warning("Ce cluster ne contient qu'une seule séquence -- pas de logo possible.")
else:
    ref_id = cluster_medoid(sequences, cluster_ids)
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    ref_seq = sequences[ref_id]
    aligned_rows = []
    for sid in cluster_ids:
        alignment = aligner.align(ref_seq, sequences[sid])[0]
        aligned_ref, aligned_query = alignment[0], alignment[1]
        row = "".join(q for r, q in zip(aligned_ref, aligned_query) if r != "-")
        aligned_rows.append(row[: len(ref_seq)].ljust(len(ref_seq), "-"))

    counts_df = logomaker.alignment_to_matrix(aligned_rows)
    fig_logo, ax = __import__("matplotlib.pyplot", fromlist=["subplots"]).subplots(figsize=(12, 3))
    logomaker.Logo(counts_df, ax=ax, color_scheme="chemistry")
    ax.set_title(f"Cluster {selected_cluster} — {len(cluster_ids)} séquences (référence : {ref_id})")
    st.pyplot(fig_logo)

st.divider()
csv_buffer = io.StringIO()
plot_df.to_csv(csv_buffer, index=False)
st.download_button("Exporter les résultats (CSV)", csv_buffer.getvalue(), "classification_resultats.csv", "text/csv")
