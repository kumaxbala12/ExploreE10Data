import streamlit as st
import scanpy as sc
import gdown
import os
import anndata as ad
import matplotlib.pyplot as plt
from anndata import read as read_h5ad


st.set_page_config(layout="wide")
st.title("E10 scRNA-seq Data Explorer")

DATA_PATH = "data/e10adata.h5ad"
GDRIVE_URL = "https://drive.google.com/uc?id=1drzl3mGt_3nKnIDRqs8cRthV6HnoXuwY"

@st.cache_data
def load_data():
    if not os.path.exists(DATA_PATH):
        os.makedirs("data", exist_ok=True)
        gdown.download(GDRIVE_URL, DATA_PATH, quiet=False)
    return sc.read_h5ad(DATA_PATH)

adata = load_data()

clusters = sorted(adata.obs['leiden'].unique()) if 'leiden' in adata.obs else []
selected_clusters = st.multiselect("Select clusters", clusters, default=clusters)

filtered_data = adata[adata.obs['leiden'].isin(selected_clusters), :] if 'leiden' in adata.obs else adata

st.write(f"Showing {filtered_data.n_obs} cells and {filtered_data.n_vars} genes.")

gene = st.text_input("Enter gene name to visualize (e.g., Sox2)")
if gene and gene in filtered_data.var_names:
    sc.pl.violin(filtered_data, keys=gene, groupby='leiden', show=False)
    st.pyplot(plt.gcf())
    plt.clf()

    sc.pl.dotplot(filtered_data, var_names=[gene], groupby="leiden", show=False)
    st.pyplot(plt.gcf())
    plt.clf()
elif gene:
    st.warning(f"{gene} not found in dataset.")
