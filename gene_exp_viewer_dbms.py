import pandas as pd
import plotly.graph_objects as go
from sqlalchemy import create_engine, MetaData, Table, select,inspect,column,text
import streamlit as st
import math
import numpy as np
from authenticator import authenticate
#layout
st.set_page_config(layout="wide")


# Custom CSS to hide "Made with Streamlit" and "About" in the hamburger menu
hide_streamlit_style = """
    <style>
    .stDeployButton {
            visibility: hidden;
        }
    [data-testid="stHeader"] {
        visibility: hidden;
    }
  
    footer {visibility: hidden;}
    
    .stApp { margin-top: -70px; }  /* Adjust the top margin */
    </style>
    """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)


# Define your database connection details
db_url = 'mysql+pymysql://wilhelm2024:XctyzMr3hr2>@mysql-lab.iric.ca/wilhelm2024'
engine = create_engine(db_url)
metadata = MetaData()
big_summary_groups = Table('big_summary_groups', metadata, autoload_with=engine)
big_summary_gene = Table('big_summary_gene', metadata, autoload_with=engine)
target_AML_groups = Table('target_AML_groups', metadata, autoload_with=engine)
target_AML_genes = Table('target_AML_genes', metadata, autoload_with=engine)

connection = engine.connect()






# Function to calculate custom median handling non-numeric values

def custom_median(series):
    numeric_series = pd.to_numeric(series, errors='coerce')
    numeric_series = numeric_series.dropna()
    return numeric_series.median()

# Streamlit cache decorator for data fetching

@st.cache_data()
def fetch_gene_data(gene,group_name,data_file):
    global result2
    if(data_file=='big_summary_updated_June_2024'):
    # Construct the select query to get the sample names for the given group
        query = select(big_summary_groups.c['Sample']).where(big_summary_groups.c[group_name] != '')
        result = connection.execute(query)
        selected_samples = [row[0] for row in result]
        
            #query = select(big_summary_groups.c['Sample']).where(big_summary_groups.c[group_name] IS NOT NULL)
        # Construct the select query to get the values for the selected samples and gene of interest
        selected_columns = [big_summary_gene.c[sample] for sample in selected_samples]
        query2 = select(*selected_columns).where(big_summary_gene.c['Gene'] == gene)
        result2 = connection.execute(query2)

    elif(data_file=='Target AML'):
        query = select(target_AML_groups.c['Sample']).where(target_AML_groups.c['Type'] == group_name)
        result = connection.execute(query)
        #st.write(result.fetchall())
        selected_samples = [row[0] for row in result]
        #st.write(selected_samples)
        selected_columns = [target_AML_genes.c[sample] for sample in selected_samples]
        

       
        query2 = select(*selected_columns).where(target_AML_genes.c['Sample_Id'] == gene)
        result2 = connection.execute(query2)
       
        
        
        
    # Execute the query and fetch results
    
    df = pd.DataFrame(result2.fetchall())
  

    # Close the connection
 

    return df

def gene_input(options,data_file, check=None):
    fig=go.Figure()

    
    for option in options:
        

        if option == 'Mutation':
            continue
        else:
            
            global result_gene
            if data_file=='big_summary_updated_June_2024':

                query_gene=big_summary_gene.select().where(big_summary_gene.c['Gene']==option)
                result_gene=connection.execute(query_gene)
                stat=pd.DataFrame(result_gene.fetchall())
                stat=stat.iloc[:,2:]

            elif data_file=='Target AML':

                query_gene=target_AML_genes.select().where(target_AML_genes.c['Sample_Id']==option)
                result_gene=connection.execute(query_gene)
                stat=pd.DataFrame(result_gene.fetchall())
                stat=stat.iloc[:,1:]

            
           #st.write(stat)
            
            
            if stat.empty:
                st.warning(f"No data found for gene:{option}")
                continue

            #stat = stat.select_dtypes(include=[np.number])

            if not stat.empty:
                stat_data = stat.T.reset_index()
                stat_data.columns = ['sample', 'FPKM']
                stat_data['FPKM'] = pd.to_numeric(stat_data['FPKM'], errors='coerce').fillna(0.0)
                #st.write(stat_data['FPKM'].dtypes)
                
                hover_y = 'FPKM'
                if check:
                    
                    stat_data['FPKM'] = np.log2(stat_data['FPKM'] + 1)
                    hover_y = 'log2(FPKM+0.001)'

                box_trace = go.Box(y=stat_data['FPKM'], boxpoints='all', jitter=0.2, pointpos=0, hoverinfo='text',
                                   hovertext=[f"Sample: {idx}<br>{hover_y}: {val}" for idx, val in
                                              zip(stat_data['sample'], stat_data['FPKM'])],
                                   boxmean=True, name=option, marker=dict(size=10))
                fig.add_trace(box_trace)

                median = custom_median(stat_data['FPKM'])
                fig.add_annotation(x=option, y=median, text=f"Median: {median:.2f}", showarrow=False, xshift=5,
                                   font=dict(color='Black'))

    if check:
        title_y = 'log2(FPKM+1)'
    else:
        title_y = 'FPKM'

    fig.update_layout(
        yaxis_title=dict(text=title_y,font=dict(size=18, family='Bold Lexand')),
        xaxis=dict(
            title=dict(text='Genes',font=dict(size=18, family='Bold Lexand')),
            tickfont=dict(size=16, family='Bold Lexand')
        ),
        boxmode='group',
        width=2500,
        height=1000,
        font=dict(size=16),
        boxgap=0,
        boxgroupgap=0,
        title='Boxplot -Gene Expression'
    )
    #print("Gene input took", time.time() - start_time, "seconds")
    return fig





# Streamlit function to create box plots
def gene_input_single(gene, group_names,data_file,log_scale=None):
    
    fig = go.Figure()
    if(gene=='Fusion'):
        st.warning("Please select a gene from the dropdown menu.")
        return fig
    
    for group_name in group_names:
        # Fetch data using the cached function
        
        data_file=data_file
        df = fetch_gene_data(gene, group_name,data_file)
        size_group = len(df.columns)
        df=df.T.reset_index()
        df.columns=['sample', 'FPKM']
       
        
        if df.empty:
            st.warning(f"No data found for gene '{gene}' in the database.")
            

        stat_data = df
       
        stat_data.columns = ['sample', 'FPKM']
  
        
        hover_y = 'FPKM'
        
        if log_scale:

            hover_y = 'log2(FPKM)'
            stat_data['FPKM'] = stat_data['FPKM'].apply(lambda x: math.log2(float(x) + 1) )






        box_trace = go.Box(
            y=stat_data['FPKM'],
            name=f"{group_name}<br>count={size_group}",
            boxpoints='all',
            jitter=0.3,
            pointpos=0,
            hoverinfo='text',
            hovertext=[f"Sample: {sample}<br>{hover_y}: {fpkm}" for sample, fpkm in zip(stat_data['sample'], stat_data['FPKM'])],
            boxmean=True,
            marker=dict(size=10)
        )
        fig.add_trace(box_trace)

        median = custom_median(stat_data['FPKM'])
        fig.add_annotation(x=f"{group_name}<br>count={size_group}", y=median, text=f"Median: {median:.2f}", showarrow=False, xshift=5,
                            font=dict(color='Black'))

    # Update layout for the plot
    if log_scale:
        title_y = 'log2(FPKM+1)'
    else:
        title_y = 'FPKM'

    fig.update_layout(
        yaxis_title=title_y,
        xaxis=dict(
            title=dict(text='Groups',font=dict(size=18, family='Bold Lexand')),
            tickfont=dict(size=16, family='Bold Lexand')
        ),
        yaxis=dict(
            title=dict(text=title_y,font=dict(size=18, family='Bold Lexand')),
            tickfont=dict(size=18, family='Bold Lexand')
        ),
        boxmode='group',
        width=1200,
        height=800,
        font=dict(size=16),
        boxgap=0.5,
        boxgroupgap=0.5,
        title=f'Boxplot - Gene Expression for {gene}'
    )

    return fig

# Streamlit app
def main():
    # Streamlit UI
    st.sidebar.title("**Gene Expression viewer**")
    #st.sidebar.subheader("**Select Data File**")
    data_files = ["big_summary_updated_June_2024", "Target AML"]
    Data_selected = st.sidebar.radio('**Select data file**', data_files)

    if not Data_selected:
        st.error("Please select a data file.")
        return



    inspector = inspect(engine)
    columns = inspector.get_columns('big_summary_groups')
    

    
    

    if "big_summary_updated_June_2024" in Data_selected:
        res=select(big_summary_gene.c['Gene']).distinct()
        res2=select(big_summary_groups.c).distinct()
        resu=connection.execute(res)
        resu2=connection.execute(res2)
        gene_list = [row[0] for row in resu if (row[0]!='Fusion' and row[0]!='MLL-AF6' and row[0]!='MLL-AF9' and row[0]!='MLL-AF10')]
        column_names = [column['name'] for column in columns if column['name']!='Sample']
        column_names=column_names[:28]
        default_genes = ["DPM1", "SCYL3"]
        default_groups = ["CD34", "CD34+MA9", "Model AML"]
        options = st.sidebar.multiselect('**Select genes**', gene_list, default=default_genes)
        
    elif "Target AML" in Data_selected:
        res=select(target_AML_genes.c['Sample_Id']).distinct()
        res2=select(target_AML_groups.c['Type']).distinct()
        resu=connection.execute(res)
        resu2=connection.execute(res2)
        
        gene_list = [row[0] for row in resu if row[0]!='Fusion']
        column_names = [ro[0] for ro in resu2]
        
        default_genes = ["DPM1", "SCYL3"]
        default_groups = ["Other", "BM"]
        options = st.sidebar.multiselect('**Select genes**', gene_list, default=default_genes)
    

    

    tab1, tab2 = st.tabs(['**Group Viewer**', '**Gene Viewer**'])

    
    

    
    check = st.sidebar.checkbox('Log2 transform Axis', key='log',value=False)



    with tab1:
        
        groups = st.multiselect("**Select Groups**", column_names, default=default_groups)
      
        for option in options:
            col1 = st.container()
            with col1:
                if groups:
                    
                    figs = gene_input_single(option,groups,Data_selected, check)
                    st.plotly_chart(figs, use_container_width=True)
                else:
                    st.title("Please select one or more groups to display the boxplot.")
                    break

    
    with tab2:
        col2 = st.container()
        with col2:
            if options:
                
                
                fig = gene_input(options,Data_selected,check)
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.title("Please select one or more genes to display the boxplot.")

if __name__ == '__main__':
    if authenticate():
        main()
    else:
        st.error("Please login to the application.")
    
