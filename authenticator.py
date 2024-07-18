import streamlit as st


def check_login():
    #st.button("Login")
    if st.session_state.username == "admin" and st.session_state.password == "admin":
        st.session_state.authenticated = True
        #st.experimental_rerun()
    else:
        st.session_state.authenticated = False
        st.error("Invalid credentials")

def authenticate():
    if "authenticated" not in st.session_state:
        st.session_state.authenticated = False

    if not st.session_state.authenticated:
        st.text_input(label="Username:",value="",key="username")
        st.text_input(label="Password:",type="password",value="",key="password",on_change=check_login)
       # check=st.button("Login")
        
    return st.session_state.authenticated
    

