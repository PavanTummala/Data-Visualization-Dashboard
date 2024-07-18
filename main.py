import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from sqlalchemy import create_engine, MetaData, Table, select,text,Column,Text
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from utilis import plots,smiles_to_figure,get_drug_name,get_dmso_data,get_distinct_drug_ids,get_drug_response_data,get_apexbio
from connection import get_engine, get_metadata, get_tables
from authenticator import authenticate

st.set_page_config(layout="wide")
# Database connection
if authenticate():
    engine = get_engine()
    metadata = get_metadata(engine)
    Drug_response, Cpds_info, drug_struct = get_tables(metadata)

    toolbarMode="viewer"
    
    hide_streamlit_style = """
        <style>
        # .stDeployButton {
        #         visibility: hidden;
        #     }
        # [data-testid="stHeader"] {
        #     visibility: hidden;
        # }
        # .css-18e3th9 {padding-top: 0rem;}  /* Adjust padding to remove space at the top */
        # .css-1d391kg {padding-top: 0rem;} 
        # .viewerBadge_container__1QSob {display: none;}
        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        
        .stApp { margin-top: -75px; }  /* Adjust the top margin */
        </style>
        """
    style_column="""
    <style> 
    .column-border-left { 
    border-left: 2px solid black;
    padding-left: 0px; 
    } 
    </style>
    """

    custom_css = """
    <style>
        .stApp {
            background-color: #f0f8ff;
        }
        
        # .st-emotion-cache-1wmy9hl.e1f1d6gn1 {
        # /* Add any general styles for this class here */
        
        #   border: 1.2px solid #4CAF50;
        # }
        div[data-testid="stVerticalBlock"] > div:has(div.stButton) {
            background-color: #e6f3ff;
            border: 2px solid #4CAF50;
            padding: 10px;
            border-radius: 10px;
            align: center;
        }
        .stButton > button {
            
            background-color: #4CAF50;
            color: white;
        }
        div[data-baseweb="select"] > div {
            background-color:#f0f0f0;
            border: 1px solid #000000;
        }
    </style>
    """
    st.markdown("""
    <style>
    .custom-container {
        background-color: #f0f0f0;
        padding: 10px;
        border-radius: 10px;
        border: 1px solid #000000;
    }
    </style>
    """, unsafe_allow_html=True)

    st.markdown(custom_css, unsafe_allow_html=True)
    st.markdown(hide_streamlit_style, unsafe_allow_html=True)


    drug_ids = get_distinct_drug_ids()
    
    drug_show=[f"{drug_ids[id][0]}-{drug_ids[id][1]}" for id in range(0,len(drug_ids))]
    #st.write(drug_ids)

    #samples = [col.name for col in Drug_response.columns]

    colu=st.container()
    colu.markdown('<div class="custom-container">Drug Response Viewer</div>', unsafe_allow_html=True)
    #col[1].write('<style> .column-border-left { border-left: 2px solid black; padding-left: 10px; } </style>', unsafe_allow_html=True)
    with colu :
        #Create a container for the border
        #st.write(get_drug_name(drug_id))  
        # drug_show=[f"{drug_id} - {get_drug_name(drug_id)}" for drug_id in drug_ids]
        # st.write(drug_show)
        # selected_drug = st.selectbox('**Select drugs**', drug_show)
        # selected_drug = selected_drug.split('-')[0]
        #apexbio=st.checkbox('ApexBio', key='apexbio')
        # if apexbio :
        #     drugids=get_apexbio()
        #     #st.write(drugids)
        #     drug_show=[f"{row[1]['Object_Id']} - {row[1]['Object_Drug_Name']}" for row in drugids.iterrows()]

        #     #st.write(drug_show)
        #     select_drug = st.selectbox('**Select drugs**', drug_show)
        #     select_drug = select_drug.split('-')
        #     selected_drug=select_drug[0]
        #     drug_name=select_drug[1]
        #     #st.write(drug_name)
        #     source='ApexBio'
        # else:
        selected_drug = st.selectbox('**Select drugs**', drug_show)
        selected_drug = selected_drug.split('-')[0]
        info,source,drug_name=get_drug_name(selected_drug)

        
        st.write(f"**Source**: {source} &nbsp; &nbsp; **Drug Name** :  ", drug_name)
        #st.write(" Drug Name :  ", drug_name)
        button = st.button('plot', key='plot')

    

        
    col=st.columns([4,1])
    with col[0]:
        if button:
            try:
                #st.markdown(style_column, unsafe_allow_html=True)
                df_result = get_drug_response_data(selected_drug)
                #st.write(df_result)
                if df_result is not None:
                    #sdf_file_path = 'Z:\Lab Science\Data_Wilhelm Lab\Lab projects\web page data viewer\drug_response viewer\L1048-Inhibitor-Library-SDF.sdf'
                    # Plot the molecules
                    with col[1]:
                        info,_,drug_name=get_drug_name(selected_drug)
                        #print(drug_name)
                        #st.write(drug_name)
                        if drug_name is not None and drug_name!='':
                            with engine.connect() as connection :
                                queryD=text("SELECT SMILES FROM Drug_Structure WHERE Item_Name =:drug_name").bindparams(drug_name=drug_name)
                                #queryD = select(Column('SMILES', Text(), table=drug_struct)).where(drug_struct.c.Item_Name == drug_name)
                                resu=connection.execute(queryD)
                                smiles=str(resu.fetchone()[0])
                                #smiles='CC(C)CC(C(=O)NC(CC(=O)O)C(=O)NC(CO)C(=O)O)NC(=O)C(C(C)C)NC(=O)C(CCCCN)NC(=O)C(CCCCN)NC(=O)C(CCSC)NC(=O)C(CCC(=O)N)NC(=O)C(CCSC)N'
                                #st.write("**SMILES string:**",smiles)
                                df=pd.DataFrame({
                                    'SMILES':[smiles],'Drug_Name':[drug_name],'Information':[info]
                                    })
                                st.write(df.T)
                                img=smiles_to_figure(smiles)
                                st.image(img)
                                #st.write(df_result)
                        else:
                            st.write("No data found")
                    df=df_result.iloc[:,1:3]
                    # df=df_result.T.reset_index()
                    df.columns=['Sample', 'Inhibition(%)']
                    #st.write(df_result)
                    dmso_df=get_dmso_data()
                    dmso_df['C_Score'] = pd.to_numeric(dmso_df['C_Score'], errors='coerce')
                    #st.write(dmso_df)
                    #st.write(dmso_df.groupby('Sample')['C_Score'])  
                    percentile_25 = dmso_df.groupby('Sample')['C_Score'].quantile(0.25).reset_index()
                    percentile_75 = dmso_df.groupby('Sample')['C_Score'].quantile(0.75).reset_index()
                    #st.write(percentile_25, percentile_75) 
                    #df_result=df_result.merge(percentiles_25, on='Sample', how='left')
                    #df_result=df_result.merge(percentile_75, on='Sample', how='left')
                    #st.write(df_result)
                    df_combines= pd.concat([percentile_25, percentile_75], axis=1)
                    df_combines.columns = ['Sample', '25th Percentile','sample_2', '75th Percentile']
                    df_combines=df_combines.drop(columns=['sample_2'])
                
                    
                    # fig = px.scatter(df_result, x='Sample',y='C_Score',height=700,width=1300,hover_data=df_result,title=f'Response for {selected_drug}')
                    # fig=go.figure()
                    # fig.update_traces(marker=dict(color='red',size=10,line=dict(color='cyan', width=2)))

                    # for _, row in df_combines.iterrows():
                    #     sample_value = row['Sample_Num']
                    #     fig.add_shape(
                    #         type="line",
                    #         x0=sample_value,
                    #         x1=sample_value,
                    #         y0=row['25th Percentile'],
                    #         y1=row['75th Percentile'],
                    #         line=dict(color="blue", width=2, dash="dash")
                    #         #hoverinfo="25th Percentile: " + str(row['25th Percentile']) + "<br>" + "75th Percentile: " + str(row['75th Percentile'])
                    #     )
                    #     fig.add_trace(fig.Scatter(
                    #     x=[sample_value, sample_value],
                    #     y=[row['25th Percentile'], row['75th Percentile']],
                    #     mode='markers',
                    #     marker=dict(opacity=0),
                    #     hoverinfo='text',
                    #     text=[f'Q1(DMSO): {row["25th Percentile"]:.3f}', f'Q3(DMSO): {row["75th Percentile"]:.3f}']
                    # ))
                    #     fig.update_layout(showlegend=False)
                    #     # max_sample_num = df_result['Sample'].max()

                    #     # # Update the layout to exclude the last 100 x-values
                    #     # fig.update_layout(
                    #     #     xaxis=dict(range=[0, 37])
                    #     # )

                    fig=plots(df_result, df_combines)
                    #fig.update_traces(marker=dict(color='blue',size=10,line=dict( width=2)))
                    fig.update_layout(showlegend=False)
                    st.plotly_chart(fig)
                else:
                    st.write("No data found for the selected drug.")
            except Exception as e:
                st.write(f"An error occurred: {str(e)}")

else:
    st.write("Please Login with correct credentials.")