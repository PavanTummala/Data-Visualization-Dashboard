from sqlalchemy import create_engine, MetaData, Table



def get_engine():
    db_url = 'mysql+pymysql://wilhelm2024:XctyzMr3hr2>@mysql-lab.iric.ca/wilhelm2024'
    return create_engine(db_url)

def get_metadata(engine):
    metadata = MetaData()
    metadata.reflect(bind=engine)
    return metadata

def get_tables(metadata):
    engine = get_engine()
    Drug_response = Table("Drug_response", metadata, autoload_with=engine)
    Cpds_info = Table("Compounds_Info", metadata, autoload_with=engine)
    drug_struct = Table("Drug_Structure", metadata, autoload_with=engine)
    return Drug_response, Cpds_info, drug_struct