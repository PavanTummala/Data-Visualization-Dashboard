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

# import streamlit as st

# def authenticate():
#     """Display a login modal dialog and return authentication status."""
#     # Initialize the login status
#     if "authenticated" not in st.session_state:
#         st.session_state.authenticated = False

#     if not st.session_state.authenticated:
#         with st.expander("Login", expanded=True):
#             st.write("Please log in to continue")
#             username = st.text_input("Username")
#             password = st.text_input("Password", type="password")
#             submit_button = st.button("Login")

#             if submit_button:
#                 if username == "admin" and password == "admin":
#                     st.session_state.authenticated = True
#                     st.success("You have successfully logged in!")
#                 else:
#                     st.error("Invalid credentials")

#     return st.session_state.authenticated


    

