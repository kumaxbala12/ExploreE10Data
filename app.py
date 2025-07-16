import streamlit as st
import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO
import gdown
import os

# ---- CONFIG ----
st.set_page_config(page_title="scRNA-seq Viewer", layout="wide")
st.title("scRNA-seq Data Explorer")

# ---- DOWNLOAD FROM GOOGLE DRIVE ----
DATA_URL = "https://drive.google.com/uc?id=1drzl3mGt_3nKnIDRqs8cRthV6HnoXuwY"
LOCAL_H5AD = "data.h5ad"

if not os.path.exists(LOCAL_H5AD):
    with st.spinner("Downloading data..."):
        gdown.download(DATA_URL, LOCAL_H5AD, quiet=False)

# ---- LOAD DATA ----
@st.cache_resource
def load_data():
    return sc.read_h5ad(LOCAL_H5AD)

adata = load_data()

st.subheader("UMAP Plot")

if "X_umap" in filtered.obsm:
    fig, ax = plt.subplots()
    sc.pl.umap(filtered, color=cluster_key, ax=ax, show=False)
    st.pyplot(fig)
else:
    st.warning("UMAP coordinates not found in the .h5ad file.")


# ---- FILTERING ----
cluster_key = "seurat_clusters"
adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)
all_clusters = sorted(adata.obs[cluster_key].unique())
selected_clusters = st.multiselect("Select clusters to include", all_clusters, default=all_clusters)

filtered = adata[adata.obs[cluster_key].isin(selected_clusters)]

# ---- PLOTS ----
st.subheader("Plots")

# Gene input
gene = st.text_input("Enter gene name (case-sensitive):")

# Layout: columns
col1, col2 = st.columns(2)

with col1:
    st.markdown("### Violin Plot")
    if gene:
        fig, ax = plt.subplots()
        sc.pl.violin(filtered, keys=gene, groupby=cluster_key, ax=ax, show=False)
        st.pyplot(fig)

with col2:
    st.markdown("### Dot Plot")
    if gene:
        fig = sc.pl.dotplot(filtered, var_names=[gene], groupby=cluster_key, show=False)
        st.pyplot(fig.figure)

# ---- Metadata Preview ----
st.subheader("Metadata")
st.dataframe(adata.obs.head())

