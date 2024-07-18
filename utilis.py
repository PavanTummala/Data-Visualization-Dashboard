
# To compartmentalize your code, you should split it into multiple files based on functionality. For example, you can separate database operations, plotting functions, and the main Streamlit application logic into different modules. Here's how you can do it:

# Database Operations: Handle all database interactions.
# Plotting Functions: Handle all plotting logic.
# Main Streamlit Application: Integrate everything together.
# Step 1: Create the Directory Structure
# plaintext
# Copy code
# my_streamlit_app/
# ├── db_operations.py
# ├── plotting.py
# ├── utils.py
# ├── styles.css
# ├── main.py
# └── requirements.txt
# Step 2: db_operations.py
# This file will handle all database operations.

# python
# Copy code

import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from sqlalchemy import create_engine, MetaData, Table, select,text,Column,Text
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from connection import get_engine, get_metadata, get_tables
import time
# Database connection
engine = get_engine()
metadata = get_metadata(engine)
Drug_response, Cpds_info, drug_struct = get_tables(metadata)

@st.cache_data
def plots(df_result, df_combines):
    
    fig = go.Figure()
    df_result['Sample'] = pd.Categorical(df_result['Sample'], ordered=True)
    df_combines['Sample'] = pd.Categorical(df_combines['Sample'], ordered=True)
    sample_mapping = {sample: i for i, sample in enumerate(df_result['Sample'].unique())}
    df_result['Sample_Num'] = df_result['Sample'].map(sample_mapping)
    #st.write(df_result)
    df_combines['Sample_Num'] = df_combines['Sample'].map(sample_mapping)

    # Add scatter plot for df_result
    #fig = px.scatter(df_result, x='Sample',y='C_Score',height=700,width=1300,hover_data=df_result,title=f'Response for {selected_drug}')
    fig.add_trace(go.Scatter(
        x=df_result['Sample_Num'],
        y=df_result['C_Score'],
        mode='markers',
        marker=dict(size=10, color='Red'),
        name='C_Score',
        hoverinfo='text',
        text=[f"Well: {well}<br>Plate: {plate}<br>Sample: {sample}<br>C_Score: {c_score}" 
          for well, plate, sample, c_score in zip(df_result['Well'], df_result['Plate'], df_result['Sample'], df_result['C_Score'])]
        )
    )

    # Add vertical lines for percentiles
    for _, row in df_combines.iterrows():
        fig.add_trace(go.Box(
        y=[row['25th Percentile'], row['75th Percentile']],
        x=[row['Sample_Num']] * 2,  # Repeating the sample number for each point
        name=f"Sample {row['Sample_Num']}",
        marker=dict(color='light blue'),
        boxpoints=False,  # Hides the original points
        hovertext=f"Q1(DMSO): {row['25th Percentile']:.3f}<br>Q3(DMSO): {row['75th Percentile']:.3f}",
        hoverinfo='text'
    ))
        fig.add_trace(go.Scatter(
        x=[row['Sample_Num'], row['Sample_Num']],
        y=[row['25th Percentile'], row['75th Percentile']],
        mode='markers',
        marker=dict(size=0.1, color='blue', opacity=0.0),
        hoverinfo='text',
        text=f"Q1(DMSO): {row['25th Percentile']:.3f}<br>Q3(DMSO): {row['75th Percentile']:.3f}",
        showlegend=False
    ))


    # Update layout
    fig.update_layout(
        title='',
        xaxis=dict(
            title='Sample',
            tickmode='array',
            tickvals=list(sample_mapping.values()),
            ticktext=list(sample_mapping.keys())
        ),
        yaxis_title='C_Score'
    )

    return fig 

@st.cache_data
def smiles_to_figure(smiles: str):
    # Convert SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    fig=plt.figure(figsize=(3,3))
    if mol is None:
        st.error("Invalid SMILES string.")
        return None

    # Generate figure from the molecule object
    fig=Draw.MolToImage(mol,size=(450,370))
    return fig

@st.cache_data
def get_drug_name(drug_id,retries=3):
    attempt = 0
    while attempt < retries:
        try:
            return get_drug(drug_id)
        except Exception as e:
            print(f"Retry {attempt + 1}/{retries} failed with error: {e}")
            attempt += 1
            time.sleep(2)  # Wait for 2 seconds before retrying
    st.error("Failed to fetch data after several retries.")
    return None,None,None

@st.cache_data
def get_drug(drug_id):
    if(drug_id =='DMSO'):
        return None
    
    with engine.connect() as connection :
        query=select(Cpds_info.c.Object_Drug_Name , Cpds_info.c.Object_Supplier,Cpds_info.c.Object_Name).where(Cpds_info.c.Object_Id==drug_id)
        result=connection.execute(query)
        info=result.fetchall()
        #st.write(info)
        if info is not None and len(info) > 0:
            source=info[0][1]
            drug_name=info[0][0]
            drug_info=info[0][2]

            #st.write(" Source :  ", source)
        
            if source =='ApexBio':
                #st.write(drug_name)
                return drug_info,source,drug_name
            else :
                return drug_info,source,None
            
        else:
            return None,None,None
        
@st.cache_data       
def get_apexbio():
    with engine.connect() as connection:
        query=text("SELECT Object_Id,Object_Drug_Name FROM wilhelm2024.Compounds_Info where Object_Supplier='ApexBio'")
        #.bindparams(ApexBio=ApexBio)
        result=connection.execute(query)
        columns = result.keys()
        data = result.fetchall()
        if data:
            df = pd.DataFrame(data, columns=columns)
            #df = df.drop(columns=['Object_Supplier'])
            return df
        
    return None
        
@st.cache_data
def get_dmso_data():
    with engine.connect() as connection:
        query = Drug_response.select().where(Drug_response.c.ObjectId == 'DMSO')
        result = connection.execute(query)
        
        columns = result.keys()
        data = result.fetchall()
        if data:
            df = pd.DataFrame(data, columns=columns)
            df = df.drop(columns=['ObjectId'])
            melted_df = pd.melt(df, id_vars=['Well', 'Plate'], var_name='Sample', value_name='C_Score')
            return melted_df
    return None

@st.cache_data
def get_distinct_drug_ids():
    with engine.connect() as connection:
        query = select(Cpds_info.c.Object_Id,Cpds_info.c.Object_Drug_Name)
        result = connection.execute(query)
        drug_ids=result.fetchall()
        #drug_ids = [{row[0]:row[1]} for row in result.fetchall() if row[0]!= 'DMSO']
    return drug_ids

@st.cache_data
def get_drug_response_data(drug_id):
    with engine.connect() as connection:
        query = Drug_response.select().where(Drug_response.c.ObjectId == drug_id)
        result = connection.execute(query)
        
        columns = result.keys()
        data = result.fetchall()
        #st.write(data)
        if data:
            df = pd.DataFrame(data, columns=columns)
            df=df.drop(columns=['ObjectId'])
            #st.write(df)
            melted_df = pd.melt(df, id_vars=['Well','Plate'], var_name='Sample',value_name='C_Score')
            #melted_df=melted_df.iloc[:,1:3]
            #st.write(melted_df)
            return melted_df
    return None           