import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

import cProfile


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
def load_data(data_path,sheet_name):
    return pd.read_excel(data_path, sheet_name=sheet_name)

data_path='Z:\Lab Science\Data_Wilhelm Lab\Lab projects\web page data viewer\dash\\big_summary_updated_July_2023.xlsx'
# Store the loaded data in a global variable





##################################################
@st.cache_data
def gene_input(options,data,highlight_samples=None):
    #global data # Access the global data variable
    fig = go.Figure()

    for gene in options:
        #for _, row in data.iterrows():
            if gene == 'Mutation':
                continue
            else:  
                row_data=data[data['Gene'] == gene]
                #st.write(row_data)
                #stat = row.iloc[2:]
                stat=row_data.iloc[:,2:]
                #st.write(stat)
                if stat.empty:
                     st.warning(f"No data found for gene:{gene}")
                     continue
                #stat = stat[stat.apply(lambda x: isinstance(x, (int, float)))]
                stat=stat.select_dtypes(include=[np.number])

                

                if not stat.empty:
                    stat_data = stat.T.reset_index()  # Convert the series to a data frame
                    stat_data.columns = ['sample', 'FPKM']  # Rename the columns

                    # Calculate quantiles, skipping non-numerical values
                    #quartiles = np.nanquantile(stat_data['FPKM'].values, [0.25, 0.5, 0.75]).tolist()
                    #min_val = stat_data['FPKM'].min()
                    #max_val = stat_data['FPKM'].max()

                    box_trace = go.Box(y=stat_data['FPKM'], boxpoints='all', jitter=0.2, pointpos=0, hoverinfo='text',
                                        hovertext=[f"Sample: {idx}<br>FPKM: {val}" for idx, val in zip(stat_data['sample'], stat_data['FPKM'])],
                                        boxmean=True, name=gene,marker=dict(size=10))  # Show the mean value
                    fig.add_trace(box_trace)

                    # Add annotations for median, max, and min values
                    median = stat_data['FPKM'].median()
                    fig.add_annotation(x=gene, y=median, text=f"Median: {median:.2f}", showarrow=False, xshift=5, font=dict(color='Black'))
                    #fig.add_annotation(x=gene, y=max_val, text=f"Max: {max_val:.2f}", showarrow=False, xshift=5, font=dict(color='Black'))
                    #fig.add_annotation(x=gene, y=min_val, text=f"Min: {min_val:.2f}", showarrow=False, xshift=5, font=dict(color='Black'))
                    #fig.add_annotation(x=gene, y=quartiles[0], text=f"Q1: {quartiles[0]:.2f}", showarrow=False, yshift=0, xshift=10, font=dict(color='Black'))
                    #fig.add_annotation(x=gene, y=quartiles[2], text=f"Q3: {quartiles[2]:.2f}", showarrow=False, yshift=0, xshift=10, font=dict(color='Black'))
                    

    fig.update_layout(
        yaxis_title='FPKM',
        boxmode='group',
        width=2500,  # Increase the width of the figure
        height=1000,  # Increase the height of the figure
        font=dict(size=14),  # Increase the font size
        boxgap=0,  # Increase the gap between boxes
        boxgroupgap=0,
        #dragmode='zoom',
        #scrollZoom=True,
        title='Boxplot -Gene Expression'  # Increase the gap between groups of boxes
    )
    return fig



@st.cache_data
def gene_input_single(gene, group_names,data,data_group):
    #global data_group # Access the global data variable
    #global data # Access the global data variable

    fig = go.Figure()
    if gene =='Mutation':
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
            stat_R=[]
            stat_R.append(stat)

            if not stat.empty:
                stat_data = stat.T.reset_index()
                stat_data.columns = ['sample', 'FPKM']

                box_trace = go.Box(
                    y=stat_data['FPKM'], 
                    name=group_name, 
                    boxpoints='all', 
                    jitter=0.3, 
                    pointpos=0, 
                    hoverinfo='text',
                    hovertext=[f"Sample: {idx}<br>FPKM: {val}" for idx, val in zip(stat_data['sample'], stat_data['FPKM'])],
                    boxmean=True,marker=dict(size=10)
                )
                fig.add_trace(box_trace)
                
                median = stat_data['FPKM'].median()
                #quartiles = np.nanquantile(stat_data['FPKM'].values, [0.25, 0.5, 0.75]).tolist()
                #min_val = stat_data['FPKM'].min()
                #max_val = stat_data['FPKM'].max()

                fig.add_annotation(x=group_name, y=median, text=f"Median: {median:.2f}", showarrow=False, xshift=5, font=dict(color='Black'))
                #fig.add_annotation(x=group_name, y=max_val, text=f"Max: {max_val:.2f}", showarrow=False, xshift=5, font=dict(color='Black'))
                #fig.add_annotation(x=group_name, y=min_val, text=f"Min: {min_val:.2f}", showarrow=False, xshift=5, font=dict(color='Black'))
                #fig.add_annotation(x=group_name, y=quartiles[0], text=f"Q1: {quartiles[0]:.2f}", showarrow=False, yshift=0, xshift=10, font=dict(color='Black'))
                #fig.add_annotation(x=group_name, y=quartiles[2], text=f"Q3: {quartiles[2]:.2f}", showarrow=False, yshift=0, xshift=10, font=dict(color='Black'))
            else:
                st.warning(f"No data found for gene '{gene}' in group '{group_name}'.")

    fig.update_layout(
        yaxis_title='FPKM',
        xaxis=dict(
             title='Gene',
             tickfont=dict(size=16,family='Lexand')
        ),
        yaxis=dict(
             title='FPKM',
             tickfont=dict(size=16,family='Lexand')
        ),
        boxmode='group',
        width=2500,
        height=1000,
        font=dict(size=14),
        boxgap=0.5,
        boxgroupgap=0.5,
        #dragmode='zoom',
        #scrollZoom=True,
        title=f'Boxplot - Gene Expression for {gene}'
    )
    
    return fig


#######################################################


data = load_data(data_path,sheet_name='summary_genes_FPKM.annotated')
data_group=load_data(data_path,sheet_name='Sheet1_new')
        
gene_list = data['Gene'].tolist()
group_list = data_group.columns.tolist()

tab1,tab2=st.tabs(['**Single Gene Viewer**','**Multi Gene Viewer**'])

st.sidebar.title("**Gene Expression viewer**")
        # Create a two-column layout

# Display the boxplot in the second column

# Display the multiselect input in the first column
default_genes=["DPM1", "SCYL3"]
default_groups=[ "CD34", "CD34+MA9", "Model AML"]
#default_groups = [group for group in default_group if group in group_list]
options=st.sidebar.multiselect('**Select genes**',gene_list,default=default_genes)

groups=st.sidebar.multiselect("**Select Groups**",group_list,default=default_groups)





#####################################################




#cProfile.run('gene_input(options,data)')




with tab2:
         
        for option in options: 
            col1 = st.container()
            with col1:
                if option and groups:
                  
                    figs=gene_input_single(option,groups,data,data_group)
                    st.plotly_chart(figs, use_container_width=True)
                else:
                    st.title("Please select one or more genes to display the boxplot.",)


with tab1:
        
        col2 = st.container()
        with col2:
            if options:
                #st.write(highlight)
                
                fig = gene_input(options,data)
                # Create a FigureWidget and enable scroll zoom
                #fig_widget = go.FigureWidget(fig)
                #fig_widget.layout.scrollZoom = True
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.title("Please select one or more genes and groups to display the boxplot.",)
