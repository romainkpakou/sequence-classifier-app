#!/usr/bin/env python3
"""
Génère un rapport HTML statique (results/report.html) à partir du jeu
d'exemple, pour consultation directe sans avoir à lancer l'application
Streamlit -- même pipeline que app.py, exécuté en mode non-interactif.
"""
import base64
import io
import sys
from pathlib import Path

import logomaker
import matplotlib
import plotly.express as px
from Bio.Align import PairwiseAligner, substitution_matrices

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.clustering import cluster_medoid, cluster_sequences
from src.dimensionality_reduction import reduce_dimensions
from src.vectorization import read_fasta, sequences_to_sgt_matrix

KAPPA = 5.0
N_CLUSTERS = 4
PERPLEXITY = 15


def fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def main():
    sequences = read_fasta(ROOT / "data" / "example_sequences.fasta")
    embedding = sequences_to_sgt_matrix(sequences, kappa=KAPPA)
    coords = reduce_dimensions(embedding, perplexity=PERPLEXITY)
    labels = cluster_sequences(embedding, n_clusters=N_CLUSTERS)

    plot_df = coords.join(labels).reset_index()
    plot_df.columns = ["id", "tsne_1", "tsne_2", "pca_variance_explained", "cluster"]
    plot_df["cluster"] = plot_df["cluster"].astype(str)
    plot_df["length"] = plot_df["id"].map(lambda i: len(sequences[i]))

    fig = px.scatter(
        plot_df, x="tsne_1", y="tsne_2", color="cluster", hover_name="id",
        hover_data={"length": True, "tsne_1": False, "tsne_2": False},
        title="Projection t-SNE des 100 globines (colorées par cluster K-means)",
    )
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color="white")))
    scatter_html = fig.to_html(full_html=False, include_plotlyjs="cdn")

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    logo_images = {}
    for cluster_id in sorted(plot_df["cluster"].unique()):
        cluster_ids = plot_df.loc[plot_df["cluster"] == cluster_id, "id"].tolist()
        if len(cluster_ids) < 2:
            continue
        ref_id = cluster_medoid(sequences, cluster_ids)
        ref_seq = sequences[ref_id]
        aligned_rows = []
        for sid in cluster_ids:
            alignment = aligner.align(ref_seq, sequences[sid])[0]
            aligned_ref, aligned_query = alignment[0], alignment[1]
            row = "".join(q for r, q in zip(aligned_ref, aligned_query) if r != "-")
            aligned_rows.append(row[: len(ref_seq)].ljust(len(ref_seq), "-"))

        counts_df = logomaker.alignment_to_matrix(aligned_rows)
        fig_logo, ax = plt.subplots(figsize=(12, 2.5))
        logomaker.Logo(counts_df, ax=ax, color_scheme="chemistry")
        ax.set_title(f"Cluster {cluster_id} — {len(cluster_ids)} séquences (référence : {ref_id})")
        logo_images[cluster_id] = fig_to_base64(fig_logo)

    table_rows = "\n".join(
        f"<tr><td>{r.id}</td><td>{r.cluster}</td><td>{r.length}</td></tr>"
        for r in plot_df.sort_values("cluster").itertuples()
    )
    logos_html = "\n".join(
        f'<h3>Cluster {cid}</h3><img src="data:image/png;base64,{img}" style="max-width:100%">'
        for cid, img in logo_images.items()
    )

    html = f"""<!DOCTYPE html>
<html lang="fr"><head><meta charset="utf-8">
<title>Classification de séquences sans alignement — rapport d'exemple</title>
<style>
body {{ font-family: -apple-system, sans-serif; max-width: 960px; margin: 40px auto; padding: 0 20px; color: #1a1a1a; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; margin: 20px 0; }}
th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: left; }}
th {{ background: #2A2A2A; color: #F0EBE0; }}
tr:nth-child(even) {{ background: #f7f5f2; }}
</style></head>
<body>
<h1>Classification de séquences sans alignement — rapport d'exemple</h1>
<p>Jeu de données : 100 séquences de globines humaines (UniProt/Swiss-Prot).
Voir <a href="../README.md">README</a> et <a href="../PLAN.md">PLAN.md</a> pour la méthodologie
complète. Ce rapport est un instantané statique — l'application Streamlit interactive
permet d'ajuster tous les paramètres (kappa SGT, nombre de clusters, perplexité).</p>

<h2>Projection t-SNE</h2>
{scatter_html}

<h2>Table des résultats</h2>
<table><tr><th>ID</th><th>Cluster</th><th>Longueur</th></tr>
{table_rows}
</table>

<h2>Logos de séquence par cluster</h2>
<p>Alignement guidé par référence (médoïde du cluster, BLOSUM62) — approximation légère,
pas une MSA rigoureuse.</p>
{logos_html}

</body></html>"""

    out_path = ROOT / "results" / "report.html"
    out_path.write_text(html, encoding="utf-8")
    print(f"Rapport généré : {out_path}")


if __name__ == "__main__":
    main()
