# Script Streamlit: Visualização Interativa de Genes por Espécie via Upload de .zip

import streamlit as st
import pandas as pd
import plotly.express as px
import os
import tempfile
import zipfile
from Bio import SeqIO

st.set_page_config(layout="wide")

st.title("Filtragem Interativa de Genes por Cobertura de Espécies (.zip)")

# === Upload do .zip ===
zip_file = st.file_uploader("Faça upload de um arquivo .zip contendo arquivos .fasta/.fas:", type="zip")

extensoes_aceitas = ('.fasta', '.fas')

if zip_file is not None:
    with tempfile.TemporaryDirectory() as tempdir:
        zip_path = os.path.join(tempdir, "dados.zip")
        with open(zip_path, "wb") as f:
            f.write(zip_file.read())

        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(tempdir)

        arquivos = [f for f in os.listdir(tempdir) if f.endswith(extensoes_aceitas)]
        if not arquivos:
            st.error("Nenhum arquivo .fasta ou .fas encontrado no .zip.")
        else:
            genes_especies = {}
            for arquivo in arquivos:
                caminho = os.path.join(tempdir, arquivo)
                nome_gene = os.path.splitext(arquivo)[0]

                try:
                    especies = set()
                    for seq in SeqIO.parse(caminho, "fasta"):
                        especies.add(seq.id.strip())

                    if especies:
                        genes_especies[nome_gene] = especies

                except Exception as e:
                    st.warning(f"Erro no arquivo {arquivo}: {e}")

            if not genes_especies:
                st.error("Nenhum gene foi processado. Verifique os arquivos no .zip.")
            else:
                todas_especies = sorted(set.union(*genes_especies.values()))
                matriz = pd.DataFrame(0, index=sorted(genes_especies.keys()), columns=todas_especies)

                for gene, especies in genes_especies.items():
                    matriz.loc[gene, list(especies)] = 1

                matriz["n_especies"] = matriz.sum(axis=1)

                st.markdown(f"🔢 Total de espécies detectadas na matriz original: **{matriz.shape[1] - 1}**")

                min_especies = st.slider("Mínimo de espécies representadas por gene", min_value=1, max_value=int(matriz.shape[1]-1), value=10)

                matriz_filtrada = matriz[matriz["n_especies"] >= min_especies].drop(columns="n_especies")
                genes_filtrados = matriz_filtrada.index.tolist()
                especies_retidas = matriz_filtrada.columns[(matriz_filtrada.sum(axis=0) > 0)].tolist()

                st.markdown(f"✅ Genes mantidos: **{len(genes_filtrados)}**")
                st.markdown(f"✅ Espécies representadas após filtro: **{len(especies_retidas)}**")

                if len(genes_filtrados) > 0:
                    fig = px.imshow(
                        matriz_filtrada.loc[genes_filtrados, especies_retidas],
                        labels=dict(x="Espécies", y="Genes", color="Presente"),
                        color_continuous_scale="Blues",
                        aspect="auto"
                    )
                    st.plotly_chart(fig, use_container_width=True)

                with st.expander("🔍 Ver tabela de genes mantidos"):
                    st.dataframe(matriz_filtrada)

                with st.expander("🔍 Ver espécies representadas"):
                    st.dataframe(pd.DataFrame(especies_retidas, columns=["Espécie"]))

                csv_path = os.path.join(tempdir, "matriz_genes_especies_filtrada.csv")
                matriz_filtrada.to_csv(csv_path)
                st.download_button("📥 Baixar CSV filtrado", data=open(csv_path, "rb"), file_name="matriz_filtrada.csv")
