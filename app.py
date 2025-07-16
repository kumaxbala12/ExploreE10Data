
import streamlit as st
import scanpy as sc
import tempfile
import requests

st.set_page_config(page_title="Single Cell Explorer", layout="wide")

def load_from_google_drive(file_id):
    url = f"https://drive.google.com/uc?id={file_id}"
    r = requests.get(url)
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad")
    temp.write(r.content)
    temp.flush()
    return sc.read_h5ad(temp.name)

st.title("🧬 Single Cell Explorer")

file_id = "1drzl3mGt_3nKnIDRqs8cRthV6HnoXuwY"
st.write("Loading dataset from Google Drive...")
adata = load_from_google_drive(file_id)

st.success(f"Loaded dataset with {adata.n_obs} cells and {adata.n_vars} genes.")

# Example visualization
if 'X_umap' in adata.obsm:
    import plotly.express as px
    import pandas as pd
    df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
    df['cluster'] = adata.obs['leiden'] if 'leiden' in adata.obs else 'NA'
    fig = px.scatter(df, x='UMAP1', y='UMAP2', color='cluster', title='UMAP Plot')
    st.plotly_chart(fig)
else:
    st.warning("UMAP not found in this dataset.")
