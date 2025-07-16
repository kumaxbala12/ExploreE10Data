import streamlit as st
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import pandas as pd

# Title
st.title("🧬 Explore E10 scRNA-seq Data")

# Load data from Google Drive-mounted path
@st.cache_resource
def load_data():
    return sc.read_h5ad("/mount/src/exploree10data/seurat_object_rna.h5ad")

adata = load_data()

# Show metadata options for filtering
if 'seurat_clusters' in adata.obs.columns:
    cluster_options = sorted(adata.obs['seurat_clusters'].unique())
    selected_clusters = st.multiselect("Select clusters to include", cluster_options, default=cluster_options)
    filtered_data = adata[adata.obs['seurat_clusters'].isin(selected_clusters)].copy()
else:
    st.warning("No 'seurat_clusters' found in metadata. Skipping cluster filtering.")
    filtered_data = adata

# UMAP Plot
st.subheader("UMAP")
color_by = st.selectbox("Color UMAP by:", filtered_data.obs.columns)
sc.pl.umap(filtered_data, color=color_by, show=False)
fig = plt.gcf()
st.pyplot(fig)
plt.clf()

# Violin Plot
st.subheader("Violin Plot")
gene = st.text_input("Enter gene name for violin plot:")
if gene:
    try:
        sc.pl.violin(filtered_data, keys=gene, groupby='seurat_clusters', show=False)
        fig = plt.gcf()
        st.pyplot(fig)
        plt.clf()
    except KeyError:
        st.error(f"Gene '{gene}' not found in dataset.")

# Dot Plot
st.subheader("Dot Plot")
genes = st.text_area("Enter gene names for dot plot (comma-separated):")
if genes:
    gene_list = [g.strip() for g in genes.split(",")]
    try:
        sc.pl.dotplot(filtered_data, var_names=gene_list, groupby='seurat_clusters', show=False)
        fig = plt.gcf()
        st.pyplot(fig)
        plt.clf()
    except Exception as e:
        st.error(f"Error in generating dot plot: {e}")
