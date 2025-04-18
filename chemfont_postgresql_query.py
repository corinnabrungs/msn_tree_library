from configparser import ConfigParser

import logging
from functools import lru_cache

import pandas as pd
import psycopg2

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)

CHEMFONT_SQL = """
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

conn = None
is_connected = False

INFO = """
1. Download postgresql 15 or later and install. On Windows: Locate pgAdmin or psql and create a new database, e.g.:
"C:/Program Files/PostgreSQL/15/pgAdmin 4/bin/pgAdmin4.exe"
"C:/Program Files/PostgreSQL/15/scripts/runpsql.bat"
2. Windows users might want to add "cmd.exe /c chcp 1252" as the first line in the runpsql.bat script to avoid the warning about characters space not matching.
3. Download the postgresql dump from https://www.chemfont.ca/simple/download and unzip the .zip file (e.g., by 7-zip). 
4. Run runpsql.bat and connect to the database (with the chosen name), in this console (make sure to use / instead of \\):
5. Load all data by calling the command \i C:/data/chemfont.sql
6. Connect to the database using this script after changing the user, password, and database name in chemfont_database.ini 
"""


def config(filename="chemfont_database.ini", section="postgresql"):
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


def chemfont_for_row(row):
    global is_connected, conn
    if not is_connected:
        logging.info("First connect to the Chemfont database")
        logging.warning(INFO)
        return None, None
    try:
        if row is None:
            raise ValueError("Row needs to be defined")

        if "product_name" in row:
            logging.info("Running row {}".format(row["product_name"]))
        for column_name, sql_condition in EXTERNAL_IDS.items():
            try:
                value = row.get(column_name)
                if notnull_not_empty(value):
                    with conn.cursor() as cur:
                        try:
                            query = CHEMFONT_SQL.format(
                                sql_condition.format(str(value))
                            )
                            # logging.info(query)
                            cur.execute(query)
                            structure = cur.fetchone()
                            if structure:
                                return cur.description, structure
                        except Exception as err:
                            # pass exception to function
                            logging.exception("Error in postgresql chemfont query")
                            # rollback the previous transaction before starting another
                            conn.rollback()

            except Exception as e:
                logging.exception("Something went wrong while querying chemfont")

        return None, None
    except (Exception, psycopg2.DatabaseError) as error:
        logging.warning(error)
        return None, None


@lru_cache
def chemfont_postgresql(inchikey=None, split_inchikey=None):
    global is_connected
    if not is_connected:
        logging.info("First connect to the chemfont database")
        logging.warning(INFO)
        return None, None
    try:
        structure = None

        if not (inchikey or split_inchikey):
            raise ValueError("At least one structure identifier need to be a value")

        with conn.cursor() as cur:
            if inchikey:
                cur.execute(
                    CHEMFONT_SQL.format(EXTERNAL_IDS["inchikey"].format(inchikey))
                )
                structure = cur.fetchone()
            if not structure and split_inchikey:
                cur.execute(
                    CHEMFONT_SQL.format(
                        EXTERNAL_IDS["split_inchikey"].format(split_inchikey)
                    )
                )
                structure = cur.fetchone()
            if not structure:
                logging.info(
                    "NO Chemfont match FOR: Inchikey:{} and split_inchikey:{}".format(
                        inchikey, split_inchikey
                    )
                )
                return None, None
            else:
                return cur.description, structure
    except (Exception, psycopg2.DatabaseError) as error:
        logging.warning(error)
        return None, None
