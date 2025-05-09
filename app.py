# Script Streamlit: Visualização Interativa de Genes por Amostra via Upload de .zip

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.io as pio
import os
import tempfile
import zipfile
from Bio import SeqIO
from io import StringIO

st.set_page_config(layout="wide")

st.title("Filtragem Interativa de Genes por Cobertura de Amostras (.zip)")

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
            genes_amostras = {}
            for arquivo in arquivos:
                caminho = os.path.join(tempdir, arquivo)
                nome_gene = os.path.splitext(arquivo)[0]

                try:
                    amostras = set()
                    for seq in SeqIO.parse(caminho, "fasta"):
                        amostras.add(seq.id.strip())

                    if amostras:
                        genes_amostras[nome_gene] = amostras

                except Exception as e:
                    st.warning(f"Erro no arquivo {arquivo}: {e}")

            if not genes_amostras:
                st.error("Nenhum gene foi processado. Verifique os arquivos no .zip.")
            else:
                todas_amostras = sorted(set.union(*genes_amostras.values()))
                matriz = pd.DataFrame(0, index=sorted(genes_amostras.keys()), columns=todas_amostras)

                for gene, amostras in genes_amostras.items():
                    matriz.loc[gene, list(amostras)] = 1

                matriz["n_amostras"] = matriz.sum(axis=1)
                total_amostras = len(todas_amostras)

                st.markdown(f"🔢 Total de amostras detectadas na matriz original: **{total_amostras}**")

                min_val, max_val = st.slider(
                    "Intervalo de amostras representadas por gene",
                    min_value=1,
                    max_value=total_amostras,
                    value=(10, 20),
                    step=1
                )

                matriz_filtrada = matriz[
                    (matriz["n_amostras"] >= min_val) &
                    (matriz["n_amostras"] <= max_val)
                ].drop(columns="n_amostras")

                genes_filtrados = matriz_filtrada.index.tolist()
                amostras_retidas = matriz_filtrada.columns[(matriz_filtrada.sum(axis=0) > 0)].tolist()

                st.markdown(f"✅ Genes mantidos: **{len(genes_filtrados)}**")
                st.markdown(f"✅ Amostras representadas após filtro: **{len(amostras_retidas)}**")

                # Permitir destacar genes com baixa cobertura
                st.markdown("### 🔍 Destacar genes com baixa cobertura")
                cobertura_minima = st.slider(
                    "Cobertura mínima desejada (%)",
                    min_value=0,
                    max_value=100,
                    value=80,
                    step=5
                )
                cobertura_limite = total_amostras * (cobertura_minima / 100)
                genes_destacados = matriz_filtrada.index[matriz_filtrada.sum(axis=1) < cobertura_limite].tolist()

                if len(genes_filtrados) > 0 and len(amostras_retidas) > 0:
                    altura = min(max(300, 20 * len(genes_filtrados)), 1200)
                    largura = min(max(600, 12 * len(amostras_retidas)), 4000)

                    # Ajuste automático de fonte
                    font_size = 10
                    if len(amostras_retidas) > 100:
                        font_size = 8
                    if len(amostras_retidas) > 200:
                        font_size = 6

                    fig = px.imshow(
                        matriz_filtrada.loc[genes_filtrados, amostras_retidas],
                        labels=dict(x="Amostras", y="Genes"),
                        color_continuous_scale=[[0, "white"], [1, "black"]],
                        zmin=0, zmax=1,
                        aspect="auto"
                    )
                    fig.update_layout(
                        height=altura,
                        width=largura,
                        margin=dict(l=150, b=180, t=50, r=20),
                        xaxis_tickangle=45,
                        xaxis_tickfont=dict(size=font_size),
                        yaxis_tickfont=dict(size=10),
                        coloraxis_showscale=False,
                        title=f"Genes com cobertura abaixo de {cobertura_minima}% destacados em vermelho"
                    )

                    # Adicionar destaque
                    for gene in genes_destacados:
                        fig.add_shape(
                            type="rect",
                            x0=-0.5,
                            x1=len(amostras_retidas)-0.5,
                            y0=genes_filtrados.index(gene)-0.5,
                            y1=genes_filtrados.index(gene)+0.5,
                            line=dict(color="red", width=2),
                            fillcolor="rgba(0,0,0,0)",
                            layer="above"
                        )

                    st.plotly_chart(fig, use_container_width=True)

                    # Exportar gráfico como HTML interativo
                    st.markdown("### Exportar gráfico interativo (formato HTML):")
                    html_buffer = StringIO()
                    fig.write_html(html_buffer)
                    st.download_button("📤 Baixar gráfico como HTML", data=html_buffer.getvalue(), file_name="heatmap_interativo.html", mime="text/html")

                with st.expander("🔍 Ver tabela de genes mantidos"):
                    st.dataframe(matriz_filtrada)
                    csv_genes = os.path.join(tempdir, "genes_mantidos.csv")
                    matriz_filtrada.to_csv(csv_genes)
                    st.download_button("⬇️ Baixar genes em CSV", data=open(csv_genes, "rb"), file_name="genes_mantidos.csv")

                with st.expander("🔍 Ver amostras representadas"):
                    amostras_df = pd.DataFrame(amostras_retidas, columns=["Amostra"])
                    st.dataframe(amostras_df)
                    csv_amostras = os.path.join(tempdir, "amostras_representadas.csv")
                    amostras_df.to_csv(csv_amostras, index=False)
                    st.download_button("⬇️ Baixar amostras em CSV", data=open(csv_amostras, "rb"), file_name="amostras_representadas.csv")

                csv_path = os.path.join(tempdir, "matriz_genes_amostras_filtrada.csv")
                matriz_filtrada.to_csv(csv_path)
                st.download_button("📥 Baixar CSV filtrado (matriz)", data=open(csv_path, "rb"), file_name="matriz_filtrada.csv")
