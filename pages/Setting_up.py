import streamlit as st

st.title('Setting up a development environment')
st.subheader('Installing python with Anaconda')
st.write('''The Anaconda package manager is a good choice for 
        scientific computing and it's fast to get started with.''')
st.markdown('''1. Download and install Anaconda with default 
           settings during installation: https://www.anaconda.com/products/distribution''')
st.markdown('''2. Search for 'Anaconda prompt' on your computer and open it.
            ''')
st.markdown('''3. Create an isolated development environment using the
            code below. Chose a name for the environment.
            ''')
code = '''conda create --name WRITE_NAME_HERE python '''
st.code(code, language='python')
st.markdown('''4. Activate the environment with the following line:
            ''')
code_2 = '''conda activate YOUR_ENV_NAME'''
st.code(code_2, language='python')
st.markdown('''5. Install python libraries using the pip command.
            Libraries like Pandas and Jupyter are good starting points:
            ''')
code_3 = '''pip install pandas jupyter '''
st.code(code_3, language='python')
st.markdown('''6. The code installed in this environment is now isolated
            from the rest of your python installation. This is important if you need
            to work with older software or versions incompatible with your python version.
            A jupyter notebook can be used to execute individual lines of code and see thier output.
            To create a notebook use this command in the anaconda prompt:
            ''')
code_4 = '''jupyter notebook'''
st.code(code_4, language='python')

