# Script Streamlit: Visualização Interativa de Genes por Espécie

import streamlit as st
import pandas as pd
import plotly.express as px
import os
from Bio import SeqIO

st.set_page_config(layout="wide")

st.title("Filtragem Interativa de Genes por Cobertura de Espécies")

# === CONFIGURAÇÃO ===
pasta_alinhamentos = st.text_input("Caminho da pasta com arquivos .fasta/.fas:", value="C:/Users/User/CAMINHO/PARA/SEUS/ARQUIVOS")
extensoes_aceitas = ('.fasta', '.fas')

if pasta_alinhamentos and os.path.isdir(pasta_alinhamentos):
    # === COLETA DE DADOS ===
    genes_especies = {}
    for arquivo in os.listdir(pasta_alinhamentos):
        if not arquivo.endswith(extensoes_aceitas):
            continue

        caminho = os.path.join(pasta_alinhamentos, arquivo)
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
        st.error("Nenhum gene foi processado. Verifique os arquivos na pasta.")
    else:
        todas_especies = sorted(set.union(*genes_especies.values()))
        matriz = pd.DataFrame(0, index=sorted(genes_especies.keys()), columns=todas_especies)

        for gene, especies in genes_especies.items():
            matriz.loc[gene, list(especies)] = 1

        # === Processamento ===
        matriz["n_especies"] = matriz.sum(axis=1)

        # === Slider ===
        min_especies = st.slider("Mínimo de espécies representadas por gene", min_value=1, max_value=int(matriz.shape[1]-1), value=10)

        matriz_filtrada = matriz[matriz["n_especies"] >= min_especies].drop(columns="n_especies")
        genes_filtrados = matriz_filtrada.index.tolist()
        especies_retidas = matriz_filtrada.columns[(matriz_filtrada.sum(axis=0) > 0)].tolist()

        st.markdown(f"✅ Genes mantidos: **{len(genes_filtrados)}**")
        st.markdown(f"✅ Espécies representadas: **{len(especies_retidas)}**")

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

        # === Exportar CSV ===
        matriz_filtrada.to_csv("matriz_genes_especies_filtrada.csv")
        st.success("Arquivo 'matriz_genes_especies_filtrada.csv' salvo com sucesso.")

else:
    st.info("Informe um caminho válido para iniciar a análise.")
