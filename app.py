import streamlit as st
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import gdown
import os

# ---- CONFIG ----
st.set_page_config(page_title="scRNA-seq Viewer", layout="wide")
st.title("🔬 scRNA-seq Data Explorer")

# ---- DOWNLOAD FROM GOOGLE DRIVE ----
DATA_URL = "https://drive.google.com/uc?id=1drzl3mGt_3nKnIDRqs8cRthV6HnoXuwY"
LOCAL_H5AD = "data.h5ad"

if not os.path.exists(LOCAL_H5AD):
    with st.spinner("Downloading data from Google Drive..."):
        gdown.download(DATA_URL, LOCAL_H5AD, quiet=False)

# ---- LOAD DATA ----
@st.cache_resource
def load_data():
    return sc.read_h5ad(LOCAL_H5AD)

adata = load_data()

# ---- UMAP CHECK ----
with st.spinner("Checking and computing UMAP if needed..."):
    if "X_umap" not in adata.obsm:
        st.warning("UMAP not found. Running preprocessing and UMAP...")
        if "X_pca" not in adata.obsm:
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            adata = adata[:, adata.var.highly_variable]
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

# ---- UMAP PLOT ----
st.subheader("UMAP")
color_by = st.selectbox("Color UMAP by:", list(adata.obs.columns), index=adata.obs.columns.get_loc("seurat_clusters") if "seurat_clusters" in adata.obs else 0)
fig, ax = plt.subplots()
sc.pl.umap(adata, color=color_by, ax=ax, show=False)
st.pyplot(fig)
plt.clf()

# ---- CLUSTER FILTERING ----
cluster_key = "seurat_clusters"
if cluster_key in adata.obs.columns:
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)
    clusters = sorted(adata.obs[cluster_key].unique())
    selected = st.multiselect("Select clusters", clusters, default=clusters)
    filtered = adata[adata.obs[cluster_key].isin(selected)]
else:
    st.warning(f"'{cluster_key}' not found in metadata. Showing full dataset.")
    filtered = adata

# ---- GENE EXPRESSION PLOTS ----
st.subheader("Gene Expression")

gene = st.text_input("Enter a gene name for plots (case-sensitive):")

col1, col2 = st.columns(2)
if gene:
    try:
        with col1:
            st.markdown("#### Violin Plot")
            fig, ax = plt.subplots()
            sc.pl.violin(filtered, keys=gene, groupby=cluster_key, ax=ax, show=False)
            st.pyplot(fig)
            plt.clf()

        with col2:
            st.markdown("#### Dot Plot")
            fig = sc.pl.dotplot(filtered, var_names=[gene], groupby=cluster_key, show=False)
            st.pyplot(fig)
            plt.clf()
    except Exception as e:
        st.error(f"Could not plot for gene '{gene}'. Error: {e}")

# ---- METADATA TABLE ----
st.subheader("🔍 Metadata Preview")
st.dataframe(adata.obs.head())
