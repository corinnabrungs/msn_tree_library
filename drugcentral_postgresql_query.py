from configparser import ConfigParser

import logging
from functools import lru_cache

from meta_constants import MetaColumns
from pandas_utils import notnull, notnull_not_empty

import pandas as pd
import psycopg2

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)

DRUGCENTRAL_SQL = """
select structures.id, cas_reg_no, structures.name, structures.inchikey, inchi, smiles, structures.stem,
       mrdef, status,
       string_agg(distinct concat_ws(':', i.id_type, i.identifier), ';') as "database_id",
       string_agg(distinct concat_ws(':', a.type, approval::text), ';') as "date_of_approval",
       string_agg(distinct a.type, ',') as "administration",
       string_agg(distinct orphan::text, ',') as "orphan"
from structures
left join approval a on structures.id = a.struct_id
left join identifier i on structures.id = i.struct_id

-- condition defined in code
where {}
GROUP BY structures.id, cas_reg_no, structures.name, inchi, smiles, structures.stem, mrdef, structures.inchikey, status
limit 1;
"""

DRUGCENTRAL_ADD_SQL = [
    """
select structures.id, 
       string_agg(distinct syn.name, ',') as "synonyms",
       string_agg(distinct pc.type, ',') as "pharma_type",
       string_agg(distinct pc.name, ',') as "pharma_class"
from structures
left join synonyms syn on structures.id = syn.id
left join pharma_class pc on structures.id = pc.struct_id
-- condition defined in code
where structures.id = %s
GROUP BY structures.id
limit 1;
""",
    """
select structures.id, 
       string_agg(distinct concat_ws('', ai.quantity::text, ai.unit), ',') as "dosage",
       string_agg(distinct stem.definition, ',') as "stem_definition"
from structures
left join active_ingredient ai on structures.id = ai.struct_id
left join inn_stem stem on structures.stem = stem.stem
-- condition defined in code
where structures.id = %s
GROUP BY structures.id
limit 1;
""",
    """
select structures.id,
       string_agg(distinct (s2p.parent_id::text), ',') as "parent_id",
       string_agg(distinct str.strength, ';') as "strength"
from structures
left join struct2parent s2p on structures.id = s2p.struct_id
left join struct2obprod str on structures.id = str.struct_id
-- condition defined in code
where structures.id = %s
GROUP BY structures.id
limit 1;
""",
    """
select structures.id,
       string_agg(distinct concat_ws('', atc_ddd.ddd::text, atc_ddd.unit_type), ',') as "who_defined_daily_dose",
       string_agg(distinct s2atc.atc_code, ',') as "atc"
from structures
left join atc_ddd on structures.id = atc_ddd.struct_id
left join struct2atc s2atc on structures.id = s2atc.struct_id
-- condition defined in code
where structures.id = %s
GROUP BY structures.id
limit 1;
""",
    """
select structures.id,
       string_agg(distinct concat_ws(':', o.relationship_name, o.concept_name), ';') as "indication; contraindication; off_label"
from structures
left join omop_relationship o on structures.id = o.struct_id
-- condition defined in code
where structures.id = %s
GROUP BY structures.id
limit 1;
""",
]

DRUGCENTRAL_ADD2_SQL = """
select structures.id,
       string_agg(distinct stem.definition, ',') as "stem_definition",
       string_agg(distinct (s2p.parent_id::text), ',') as "parent_id",
       string_agg(distinct str.strength, ';') as "strength",
       string_agg(distinct concat_ws('', atc_ddd.ddd::text, atc_ddd.unit_type), ',') as "who_defined_daily_dose",
       string_agg(distinct concat_ws(':', o.relationship_name, o.concept_name), ';') as "indication; contraindication; off_label",
       string_agg(distinct v.species, ';') as "animal_species",
       string_agg(distinct concat_ws(':', v.relationship_type, v.concept_name), ';') as "indication"
from structures
left join inn_stem stem on structures.stem = stem.stem
left join struct2parent s2p on structures.id = s2p.struct_id
left join atc_ddd on structures.id = atc_ddd.struct_id
left join omop_relationship o on structures.id = o.struct_id
left join struct2obprod str on structures.id = str.struct_id
left join vetomop v on structures.id = v.struct_id

-- condition defined in code
where structures.id = %s
GROUP BY structures.id
limit 1;
"""

#  dict column name in our dataframe, and the SQL query where condition
EXTERNAL_IDS = {
    "unii": "i.id_type = 'UNII' and i.identifier = %s",
    "drugbank_id": "i.id_type = 'DRUGBANK_ID' and i.identifier = %s",
    "chembl_id": "i.id_type = 'ChEMBL_ID' and i.identifier = %s",
    "pubchem_cid": "i.id_type = 'PUBCHEM_CID' and i.identifier = %s",
    "input_pubchem_cid": "i.id_type = 'PUBCHEM_CID' and i.identifier = %s",
    "inchikey": "structures.inchikey = %s",
    "compound_name": "lower(structures.name) = lower(%s)",
    "split_inchikey": "structures.inchikey ~ %s",
}

conn = None
is_connected = False

INFO = """
1. Download postgresql 15 or later and install. On Windows: Locate pgAdmin or psql and create a new database, e.g.:
"C:/Program Files/PostgreSQL/15/pgAdmin 4/bin/pgAdmin4.exe"
"C:/Program Files/PostgreSQL/15/scripts/runpsql.bat"
2. Windows users might want to add "cmd.exe /c chcp 1252" as the first line in the runpsql.bat script to avoid the warning about characters space not matching.
3. Download the postgresql dump from https://drugcentral.org/download and unzip the .gz file (e.g., by 7-zip). 
4. Run runpsql.bat and connect to the database (with the chosen name), in this console (make sure to use / instead of \\):
5. Load all data by calling the command \i C:/data/drugcentral_dump_xy.sql
6. Connect to the database using this script after changing the user, password, and database name in drugcentral_database.ini 
"""


def config(filename="drugcentral_database.ini", section="postgresql"):
    # create a parser
    parser = ConfigParser()
    # read config file
    parser.read(filename)

    # get section, default to postgresql
    db = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            db[param[0]] = param[1]
    else:
        raise Exception(
            "Section {0} not found in the {1} file".format(section, filename)
        )

    return db


def connect():
    global conn, is_connected
    try:
        # read connection parameters
        params = config()

        # connect to the PostgreSQL server
        print("Connecting to the PostgreSQL database...")
        conn = psycopg2.connect(**params)

        # execute a statement
        print("PostgreSQL database version:")
        with conn.cursor() as cur:
            cur.execute("SELECT version()")
            # display the PostgreSQL database server version
            db_version = cur.fetchone()
            print(db_version)

        is_connected = True
    except (Exception, psycopg2.DatabaseError) as error:
        logging.warning(error)
        logging.warning(INFO)


def deconnect():
    global conn, is_connected
    is_connected = False
    if conn is not None:
        conn.close()
        print("Database connection closed.")


def drugcentral_for_row(row):
    global is_connected, conn
    if not is_connected:
        logging.info("First connect to the DrugCentral database")
        logging.warning(INFO)
        return None, None
    try:
        if row is None:
            raise ValueError("Row needs to be defined")

        if "compound_name" in row:
            logging.info("Running row {}".format(row["compound_name"]))
        for column_name, sql_condition in EXTERNAL_IDS.items():
            try:
                value = row.get(column_name)
                if column_name == "split_inchikey":
                    value = f"^{value}"  # regular expression match: startswith
                if notnull_not_empty(value):
                    with conn.cursor() as cur:
                        try:
                            query = DRUGCENTRAL_SQL.format(sql_condition)
                            # logging.info(query)
                            cur.execute(query, (value,))
                            structure = cur.fetchone()
                            if structure:
                                return cur.description, structure
                        except Exception as err:
                            # pass exception to function
                            logging.exception(
                                f"Error in postgresql drugcentral query for {column_name} with value {value}"
                            )
                            # rollback the previous transaction before starting another
                            conn.rollback()

            except Exception as e:
                logging.exception("Something went wrong while querying drugcentral")

        return None, None
    except (Exception, psycopg2.DatabaseError) as error:
        logging.warning(error)
        return None, None


def drugcentral_additional_query(dc_id, sql_query):
    global is_connected, conn
    if not is_connected:
        logging.info("First connect to the DrugCentral database")
        logging.warning(INFO)
        return None, None
    try:
        try:
            if notnull_not_empty(dc_id):
                with conn.cursor() as cur:
                    try:
                        # query = sql_query.format(str(dc_id))
                        # logging.info(query)
                        cur.execute(sql_query, (dc_id,))
                        structure = cur.fetchone()
                        if structure:
                            return cur.description, structure
                    except Exception as err:
                        # pass exception to function
                        logging.exception("Error in postgresql drugcentral query")
                        # rollback the previous transaction before starting another
                        conn.rollback()

        except Exception as e:
            logging.exception("Something went wrong while querying drugcentral")

        return None, None
    except (Exception, psycopg2.DatabaseError) as error:
        logging.warning(error)
        return None, None


@lru_cache
def drugcentral_postgresql(inchikey=None, split_inchikey=None):
    global is_connected
    if not is_connected:
        logging.info("First connect to the DrugCentral database")
        logging.warning(INFO)
        return None, None
    try:
        structure = None

        if not (inchikey or split_inchikey):
            raise ValueError("At least one structure identifier need to be a value")

        with conn.cursor() as cur:
            if inchikey:
                cur.execute(
                    DRUGCENTRAL_SQL.format(EXTERNAL_IDS["inchikey"].format(inchikey))
                )
                structure = cur.fetchone()
            if not structure and split_inchikey:
                cur.execute(
                    DRUGCENTRAL_SQL.format(
                        EXTERNAL_IDS["split_inchikey"].format(split_inchikey)
                    )
                )
                structure = cur.fetchone()
            if not structure:
                logging.info(
                    "NO Drugcentral match FOR: Inchikey:{} and split_inchikey:{}".format(
                        inchikey, split_inchikey
                    )
                )
                return None, None
            else:
                return cur.description, structure
    except (Exception, psycopg2.DatabaseError) as error:
        logging.warning(error)
        return None, None
