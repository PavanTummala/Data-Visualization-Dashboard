import streamlit as st

import streamlit.components.v1 as components

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
    
    .stApp { margin-top: -95px; }  /* Adjust the top margin */
    </style>
    """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)




bg= """
    <style>
    body {
        background-color:#5499C7 ;
        
    }
    </style>
    <h1 style="color:white;text-align:center;font-size:50px">IRIC Tools</h1>

    """

components.html(bg,height=100)
 



col1,col2=st.columns([1.5,1])


with col1:
        row1 = st.container()
        row2 = st.columns(1)
        row3 = st.columns(1)
        row4 = st.columns(1)
        row5 = st.columns(1)

      


        buttons = ["MISTIC", "DASH","MICKYVIEW"]
        titles = ["Gene Correlation", "Dose Response","Single Dose Response"]
        links = ["http://misticbw.iric.ca", "https://vdr.wilhelm.iric.ca","http://cluster.iric.ca:43000"]
        contents=["MISTIC integrates direct visualization and comparison of the gene correlation structure between datasets, analysis of the molecular causes underlying co-variations in gene expression and clinical annotation defined by the combined expression of selected biomarkers",
                "DASH allows users to explore drug response profiles across hundreds of cancer cell lines and identify drugs selectively toxic to molecular subtypes of interest",
                "MICKYVIEW is a Drug Response viewer useful for finding C_score values of required gene using scatter plots",
                "GENIE is a genexpression tool useful for finding expression levels of required gene using boxplots and variances"]
        
        with row1:
            data="Cevin and Genie are gene expression tools useful for finding expression levels of required gene using boxplots and variances"
            tit="Gene Expression"
            but1="CEVIN"
            but2="GENIE"
            lin1="https://bioinfo.iric.ca/~wilhelmb/CEVIN"
            lin2="http://cluster.iric.ca:43001"
            html=f"<h1 style='color:#3498DB;'> {tit} </h1>" 
            st.markdown(html, unsafe_allow_html=True)
            
            html_box = f"""
                <div style="border: 1px; padding: 10px; border-radius: 5px;height:80px;width:850px;overflow:auto;">{data}
                </div>
                """
        
            st.markdown(html_box, unsafe_allow_html=True)
            
            
            col1_button, col2_button = st.columns([0.2,1],gap="small")

            with col1_button:
                if st.button(f"Open {but1}", key=f'button_0', type='primary'):
                    js = f"window.open('{lin1}')"
                    html = f"""
                        <script>
                            {js}
                        </script>
                    """
                    components.html(html, height=0, width=0)

            with col2_button:
                if st.button(f"Open {but2}", key=f'button_1', type='primary'):
                    js = f"window.open('{lin2}')"
                    html = f"""
                        <script>
                            {js}
                        </script>
                    """
                    components.html(html, height=0, width=0)

             
        for i, (col, link, title, button,content) in enumerate(zip( row2 + row3 + row4 + row5, links, titles, buttons,contents)):
            with col:
                html=f"<h1 style='color:#3498DB;'> {title} </h1>" 
                st.markdown(html, unsafe_allow_html=True)
                
                html_box = f"""
                    <div style="border: 1px; padding: 10px; border-radius: 5px;height:80px;width:850px;overflow:auto;">{content}
                    </div>
                    """
            
                st.markdown(html_box, unsafe_allow_html=True)
                
                
                if st.button(f"Open {button}", key=f'cevin{i}',type='primary'):
                    js = f"window.open('{link}')"
                    html = f"""
                        <script>
                            {js}
                        </script>
                    """
                    
                    components.html(html, height=0, width=0)

with col2:

    cont=st.container(height=700)
    with cont:
            st.title(''':blue[IRIC]''')
            # ebw="""
            #     <iframe src="https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.cancertodaymag.org%2Fspring2022%2Fcracking-the-code-of-acute-myeloid-leukemia%2F&psig=AOvVaw2yy6XgRAmex6hy5oahkeL7&ust=1720661251975000&source=images&cd=vfe&opi=89978449&ved=0CBQQjRxqFwoTCKD8nOaom4cDFQAAAAAdAAAAABAX" width="100%" height="500px" ></iframe>
            #     """
            st.image('Cracking-Code-AML.jpg',width=600)
            #components.html(ebw,height=600)





