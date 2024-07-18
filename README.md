# **Ordering of Files on //wilsrv3.iric.ca**
   ```plaintext
Webpage data viewer/
│──bin
|─ drug_response_viewer
    ├── Bin                       #contains all my previous iterations and csv files i used
├   |---response_viewer/
│     ├── main.py                 # Main logic, Streamlit layout, HTML, CSS code
│     ├── utilis.py               # Utility functions
│     ├── connection.py           # MySQL connections
│     ├── authenticator.py        # Adds authentication to the website
│
|-dash
|   ├──bin
└── |--gene_viewer/
|      |--gene_exp_viewer_dbms.py
|      |--authenticator.py
|      |--Requirements.txt
|
|-WebApp_IRIC
|  ├──bin
|__|--GoTo/
|      |--App_all.py
|      |--Cracking-Code-AML.JPG                    #any image file you want to put in the right column(box)
|
        # Other applications
```

---

# **Ordering of Files on Cluster**

```plaintext
pavan_viewer/                    #directory in cluster where all my files are located
│
├── response_viewer/
│   ├── main.py                 # Main logic, Streamlit layout, HTML, CSS code
│   ├── utilis.py               # Utility functions
│   ├── connection.py           # MySQL connections
│   ├── authenticator.py        # Adds authentication to the website
│
└── gene_viewer/
|    |--gene_exp_viewer_dbms.py
|    |--authenticator.py
|    |--Requirements.txt
|
|
|__ GoTo/
|    |--App_all.py
|    |--Cracking-Code-AML.JPG
|
        # Other applications
```

## **Setting Up the Environment**

The project uses Miniconda to manage the virtual environment. Below are the steps to set up and activate the environment.

### 1. Activate the Virtual Environment

Navigate to the directory where to pavan_viewer where script files are located and run the following commands:

```bash
conda activate pavan_viewer  //my virtual environment name is pavan_viewer
```

## How to Run the Application

1. **Navigate to the Project Directory**:
   ```bash
   cd pavan_viewer/response_viewer                                             // on cluster
   Lab Science\Data_Wilhelm Lab\Lab projects\web page data viewer\Webapp_IRIC  // on wilsrv3
   ```
2. **Activate the Virtual Environment** 

Navigate to the directory where script files are located and run the following commands in the command line(if you are connected to cluster via SSH) or in VS code(if your running on your local computer):
- Note : The project uses Miniconda for managing the virtual environment. I assume you know how to set up a conda virtual environment(miniconda or anaconda) if not please refer to the installation process of anaconda/miniconda.

```bash
conda activate pavan_viewer  
```

3. **Run the Application**:
   ```bash
   python main.py            //compile and run python script
   streamlit run main.py    // opens up a new tab in your deafult where you can see the application
   ```

## **File Description of Drug Response Viewer** :

### main.py
This file contains the main logic for the application, including the Streamlit layout, HTML, and CSS code.

### utilis.py
plots()--> Responsible for plots
- This file contains all the utility functions used in the application.

### connection.py
This file contains the MySQL connections to the database. If there are any issues with the SQL connection, start troubleshooting here.

### authenticator.py
This file adds authentication to the website. The default username and password are:
- Username: `admin`
- Password: `admin`
You can always change it by simply changing username and password in the authenticator.py file in that directory.

## **File Description of Genie**

### gene_exp_viewer_dbms.py
This file contains all the functions, layout, and main logic in the same file 

### authenticator.py
This file adds the authentication to the script

## **File Description of SortingHat**

## App_all.py
This file contains the main logic and HTML

## Troubleshooting

If you encounter any issues, especially with the MySQL connection, please check the `connection.py` file first. Ensure that your database credentials are correct and that the database server is running.

## Dependencies

Ensure that you have all the necessary dependencies installed. You can usually install them using:

```bash
pip install -r requirements.txt
```

If `requirements.txt` is not provided, manually install the required packages as per the imports in the scripts.

## Contributing

Contributions are welcome! Please create a pull request with a detailed description of the changes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

By following the steps in this README file, you should be able to set up and run the Response Viewer (MickyView) web application. If you have any questions or need further assistance, feel free to reach out.
