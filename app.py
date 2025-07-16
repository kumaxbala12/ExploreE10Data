
import streamlit as st
import scanpy as sc
import pandas as pd
import plotly.express as px
import numpy as np
import requests
import tempfile
import matplotlib.pyplot as plt
import gdown
from io import BytesIO

st.set_page_config(page_title="Single Cell Explorer", layout="wide")

@st.cache_resource

@st.cache_resource
def load_from_google_drive():
    url = "https://drive.google.com/uc?id=1drzl3mGt_3nKnIDRqs8cRthV6HnoXuwY"
    output = "/tmp/temp_data.h5ad"
    gdown.download(url, output, quiet=False)
    return sc.read_h5ad(output)

if "adata" not in st.session_state:
    st.session_state.adata = load_from_google_drive()

adata = st.session_state.adata
st.sidebar.title("🧬 Navigation")
page = st.sidebar.radio("Go to", ["📁 Upload", "🔬 Preprocess", "🌐 Visualize", "🧱 Cluster", "📊 DEGs", "⏳ Pseudotime"])

if page == "📁 Upload":
    st.title("📁 Upload Data")
    st.success("Dataset loaded from Google Drive!")
    st.write(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    st.dataframe(adata.obs.head())

if page == "🔬 Preprocess":
    st.title("🔬 Preprocess Data")
    st.write("Filtering, normalization, log1p, scaling, and HVG selection.")
    if st.button("Run Preprocessing"):
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata)
        sc.pp.scale(adata, max_value=10)
        st.success("Preprocessing complete!")

if page == "🌐 Visualize":
    st.title("🌐 Dimensionality Reduction")
    if st.button("Run PCA + UMAP"):
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        st.success("UMAP ready!")
    if "X_umap" in adata.obsm:
        df = pd.DataFrame(adata.obsm['X_umap'], columns=["UMAP1", "UMAP2"])
        df['cluster'] = adata.obs['leiden'] if 'leiden' in adata.obs else "NA"
        fig = px.scatter(df, x="UMAP1", y="UMAP2", color='cluster')
        st.plotly_chart(fig)

if page == "🧱 Cluster":
    st.title("🧱 Leiden Clustering")
    res = st.slider("Resolution", 0.1, 2.0, 0.5, 0.1)
    if st.button("Run Clustering"):
        sc.tl.leiden(adata, resolution=res)
        st.success("Clustering done!")

if page == "📊 DEGs":
    st.title("📊 Differential Expression")
    if 'leiden' in adata.obs:
        cluster = st.selectbox("Select cluster", sorted(adata.obs['leiden'].unique()))
        if st.button("Run DEG"):
            sc.tl.rank_genes_groups(adata, 'leiden', groups=[cluster], reference='rest')
            result = sc.get.rank_genes_groups_df(adata, group=cluster)
            st.dataframe(result.head(10))
            fig = px.scatter(result, x="logfoldchanges", y=-np.log10(result["pvals_adj"]), hover_name="names")
            st.plotly_chart(fig)

if page == "⏳ Pseudotime":
    st.title("⏳ Pseudotime Inference")
    if st.button("Run Pseudotime (DPT)"):
        sc.tl.paga(adata)
        sc.tl.dpt(adata)
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)
        sc.pl.umap(adata, color='dpt_pseudotime', show=False)
        buf = BytesIO()
        plt.savefig(buf, format="png", bbox_inches='tight')
        st.image(buf.getvalue())
        st.success("Pseudotime complete!")
