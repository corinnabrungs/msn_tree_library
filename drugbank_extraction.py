import xml.etree.ElementTree as ET
import pandas as pd
import logging


def from_xml_node_path(node, path_list):
    if node is None:
        return None
    if len(path_list) == 1:
        return from_xml_node(node, path_list[0])
    else:
        return from_xml_node_path(node.find(path_list[0]), path_list[1:])


def from_xml_node(node, node_id):
    if node is None:
        return None
    subnode = node.find(node_id)
    return subnode.text if subnode is not None else None


def from_xml_attribute(node, node_id):
    return node.attrib.get(node_id) if node else None



def extraction_file(drugbank_file):
    logging.info("Will run on %s", drugbank_file)

    tree = ET.parse(drugbank_file)
    root = tree.getroot()

    rows = []
    resources = ["PubChem Compound", "ChEMBL"]

    for node in root:
        groups_node = node.find("groups")

        external = []
        atc_node = None
        mode_of_action = None
        food_interaction = None
        targets_node = None

        smiles = None
        inchikey = None

        try:
            calculated_properties_node = node.find("calculated-properties")
            if calculated_properties_node:
                smiles = next((from_xml_node(prop, "value") for prop in calculated_properties_node if
                               from_xml_node(prop, "kind") == "SMILES"), None)
                inchikey = next((from_xml_node(prop, "value") for prop in calculated_properties_node if
                                 from_xml_node(prop, "kind") == "InChIKey"), None)
        except:
            pass


        try:
            atc_node = node.find("atc-codes").find("atc-code")
            mode_of_action = ", ".join([level.text for level in atc_node])
        except:
            pass

        try:
            food_interaction_node = node.find("food-interactions")
            food_interaction = ", ".join([interaction.text for interaction in food_interaction_node])
        except:
            pass
        external_dict = {}
        try:
            external_identifiers_node = node.find("external-identifiers")
            external_dict = {from_xml_node(id_node, "resource"): from_xml_node(id_node, "identifier") for id_node in
                             external_identifiers_node}
        except:
            pass

        targets = None
        try:
            targets_node = node.find("targets")
            targets = [(from_xml_node(target, "name"), from_xml_node(target, "organism"),
                        from_xml_node_path(target, path_list=["actions", "action"])) for target in targets_node]
            targets = ", ".join(["({})".format("; ".join(target)) for target in targets])
        except:
            pass

        rows.append({"drugbank_id": from_xml_node(node, "drugbank-id"),
                     "compound_name": from_xml_node(node, "name"),
                     "chembl_id": external_dict.get("CHEMBL", None),
                     # currently we only extract pubchem and chembl, that's why this works (excluding chembl, only shows pubchem)
                     "pubchem_cid": external_dict.get("PubChem Compound", None),
                     "cas": from_xml_node(node, "cas-number"),
                     "unii": from_xml_node(node, "unii"),
                     "smiles": smiles,
                     "inchikey": inchikey,
                     "type": from_xml_attribute(node, "type"),
                     "description": from_xml_node(node, "description"),
                     "approved": from_xml_node(groups_node, "group"),
                     "indication": from_xml_node(node, "indication"),
                     "pharmacodynamics": from_xml_node(node, "pharmacodynamics"),
                     "mechanism_of_action": from_xml_node(node, "mechanism-of-action"),
                     "toxicity": from_xml_node(node, "toxicity"),
                     "metabolism": from_xml_node(node, "metabolism"),
                     "absorption": from_xml_node(node, "absorption"),
                     "half_life": from_xml_node(node, "half-life"),
                     "route_of_elimination": from_xml_node(node, "route-of-elimination"),
                     "clearance": from_xml_node(node, "clearance"),
                     "atc_code": from_xml_attribute(atc_node, "code"),
                     "mode_of_action": mode_of_action,
                     "food_interaction": food_interaction,
                     "targets": targets,
                     })

    db_df = pd.DataFrame(rows)

    # export metadata file
    db_df.to_csv("data/drugbank.tsv", sep="\t", index=False)

if __name__ == "__main__":
    extraction_file(r"data\drugbank_database.xml")

# download from: https://go.drugbank.com/releases/latest, approved access needed,