import streamlit as st

import streamlit.components.v1 as components



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
        row1 = st.columns(1)
        row2 = st.columns(1)
        row3 = st.columns(1)




        buttons = ["CEVIN", "MISTIC", "DASH"]
        titles = ["Gene Expression", "Gene Correlation", "Dose Response"]
        links = ["https://bioinfo.iric.ca/~wilhelmb/CEVIN", "http://misticbw.iric.ca", "https://vdr.wilhelm.iric.ca"]
        contents=[" Cevin is a genexpression tool useful for finding expression levels of required gene using boxplots and variances",
                "MiSTIC integrates direct visualization and comparison of the gene correlation structure between datasets, analysis of the molecular causes underlying co-variations in gene expression and clinical annotation defined by the combined expression of selected biomarkers",
                "DASH allows users to explore drug response profiles across hundreds of cancer cell lines and identify drugs selectively toxic to molecular subtypes of interest"]
        for i, (col, link, title, button,content) in enumerate(zip(row1 + row2 + row3, links, titles, buttons,contents)):
            with col:
                html=f"<h1 style='color:#3498DB;'> {title} </h1>" 
                st.markdown(html, unsafe_allow_html=True)
                
                html_box = f"""
                    <div style="border: 1px; padding: 10px; border-radius: 5px;height:100px;width:750px;overflow:auto;">{content}
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

    cont=st.container(height=600)
    with cont:
            st.title(''':blue[IRIC]''')
            ebw="""
                <iframe src="https://www.iric.ca/en" width="100%" height="500px" ></iframe>
                """
            components.html(ebw,height=600)





