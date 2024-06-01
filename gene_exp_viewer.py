import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from io import BytesIO
import math
import time
import csv

st.set_page_config(layout="wide")

# Custom CSS to hide "Made with Streamlit" and "About" in the hamburger menu
hide_streamlit_style = """
    <style>

    footer {visibility: hidden;}

    .stApp { margin-top: -50px; }  /* Adjust the top margin */
    </style>
    """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

# Load the data only once
@st.cache_data
def load_data(data_path):
    start_time = time.time()
    try:
        reader = pd.read_csv(data_path, encoding='utf-8-sig', low_memory=False)
        print("Data loaded in", time.time() - start_time, "seconds")
        return reader  # Convert list of dictionaries to DataFrame
    except Exception as e:
        print("Error Loading data: ", e)
        return None

# Custom function to handle non-numeric values and calculate median
def custom_median(series):
    numeric_series = pd.to_numeric(series, errors='coerce')
    numeric_series = numeric_series.dropna()
    return numeric_series.median()


@st.cache_data
def gene_input(options, data, check=None):
    start_time = time.time()
    fig = go.Figure()

    for gene in options:
        if gene == 'Mutation':
            continue
        else:
            row_data = data[data['Gene'] == gene]
            stat = row_data.iloc[:, 2:]

            if stat.empty:
                st.warning(f"No data found for gene:{gene}")
                continue

            stat = stat.select_dtypes(include=[np.number])

            if not stat.empty:
                stat_data = stat.T.reset_index()
                stat_data.columns = ['sample', 'FPKM']
                hover_y = 'FPKM'
                if check:
                    hover_y = 'log2(FPKM)'
                    stat_data['FPKM'] = stat_data['FPKM'].apply(lambda x: math.log2(float(x) + 0.001))

                box_trace = go.Box(y=stat_data['FPKM'], boxpoints='all', jitter=0.2, pointpos=0, hoverinfo='text',
                                   hovertext=[f"Sample: {idx}<br>{hover_y}: {val}" for idx, val in
                                              zip(stat_data['sample'], stat_data['FPKM'])],
                                   boxmean=True, name=gene, marker=dict(size=10))
                fig.add_trace(box_trace)

                median = custom_median(stat_data['FPKM'])
                fig.add_annotation(x=gene, y=median, text=f"Median: {median:.2f}", showarrow=False, xshift=5,
                                   font=dict(color='Black'))

    if check:
        title_y = 'log2(FPKM+0.001)'
    else:
        title_y = 'FPKM'

    fig.update_layout(
        yaxis_title=title_y,
        boxmode='group',
        width=2500,
        height=1000,
        font=dict(size=16),
        boxgap=0,
        boxgroupgap=0,
        title='Boxplot -Gene Expression'
    )
    print("Gene input took", time.time() - start_time, "seconds")
    return fig


# Function to create boxplot for single gene
@st.cache_data
def gene_input_single(gene, group_names, data, data_group, check=None):
    start_time = time.time()
    fig = go.Figure()
    if gene == 'Mutation':
        st.warning(f"No data found for gene '{gene}'.")
    else:
        gene_data = data[data['Gene'] == gene]

    if gene_data.empty:
        st.warning(f"No data found for gene '{gene}'.")
        return fig

    for group_name in group_names:
        if group_name not in data_group.columns:
            st.warning(f"Group '{group_name}' not found in data.")
            continue

        non_empty_samples = data_group.loc[data_group[group_name].notna(), 'Sample']
        stat = gene_data.loc[:, non_empty_samples]
        stat_R = []
        stat_R.append(stat)

        if not stat.empty:
            stat_data = stat.T.reset_index()
            stat_data.columns = ['sample', 'FPKM']
            hover_y = 'FPKM'
            if check:
                hover_y = 'log2(FPKM)'
                stat_data['FPKM'] = stat_data['FPKM'].apply(lambda x: math.log2(float(x) + 0.001))
                title = 'log2(FPKM+0.001)'

            box_trace = go.Box(
                y=stat_data['FPKM'],
                name=group_name,
                boxpoints='all',
                jitter=0.3,
                pointpos=0,
                hoverinfo='text',
                hovertext=[f"Sample: {idx}<br>{hover_y}: {val}" for idx, val in
                           zip(stat_data['sample'], stat_data['FPKM'])],
                boxmean=True, marker=dict(size=10)
            )
            fig.add_trace(box_trace)

            median = custom_median(stat_data['FPKM'])
            fig.add_annotation(x=group_name, y=median, text=f"Median: {median:.2f}", showarrow=False, xshift=5,
                               font=dict(color='Black'))

    if check:
        title_y = 'log2(FPKM+0.001)'
    else:
        title_y = 'FPKM'

    fig.update_layout(
        yaxis_title=title_y,
        xaxis=dict(
            title='Groups',
            tickfont=dict(size=16, family='Lexand')
        ),
        yaxis=dict(
            title=title_y,
            tickfont=dict(size=16, family='Lexand')
        ),
        boxmode='group',
        width=2500,
        height=1000,
        font=dict(size=16),
        boxgap=0.5,
        boxgroupgap=0.5,
        title=f'Boxplot - Gene Expression for {gene}'
    )
    print("Gene input single took", time.time() - start_time, "seconds")
    return fig


data = load_data('big_summary_updated_July_2023.csv')
data_group = load_data('gene_expression_viewer_groups.csv')

gene_list = data['Gene'].tolist()
group_list = data_group.columns.tolist()

tab1, tab2 = st.tabs(['**Single Gene Viewer**', '**Multi Gene Viewer**'])

st.sidebar.title("**Gene Expression viewer**")

default_genes = ["DPM1", "SCYL3"]
default_groups = ["CD34", "CD34+MA9", "Model AML"]

options = st.sidebar.multiselect('**Select genes**', gene_list, default=default_genes)
check = st.sidebar.checkbox('Log2 transform Axis')

with tab1:
    col2 = st.container()
    with col2:
        if options:
            fig = gene_input(options, data, check)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.title("Please select one or more genes and groups to display the boxplot.")

with tab2:
    groups = st.multiselect("**Select Groups**", group_list, default=default_groups)

    for option in options:
        col1 = st.container()
        with col1:
            if option and groups:
                figs = gene_input_single(option, groups, data, data_group, check)
                st.plotly_chart(figs, use_container_width=True)
            else:
                st.title("Please select one or more genes to display the boxplot.")
