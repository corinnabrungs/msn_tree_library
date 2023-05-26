from enum import auto

from strenum import StrEnum


class MetaColumns(StrEnum):
    # run dates
    date_pubchem_cid_search = auto()
    date_pubchem_name_search = auto()
    date_pubchem_structure_search = auto()
    date_unichem_search = auto()
    date_chembl_search = auto()
    date_npclassifier = auto()
    date_classyfire = auto()
    date_npatlas = auto()
    date_broad_drug_list = auto()
    date_drugbank_search = auto()
    date_drugcentral_search = auto()
    date_wikidata_lotus_search = auto()
    #
    structure_source = auto()

    # sample sequence specific
    well_location = auto()
    plate_id = auto()
    # calculated properties
    monoisotopic_mass = auto()
    # compound specific
    compound_name = auto()
    input_name = auto()
    synonyms = auto()
    formula = auto()
    smiles = auto()  # input smiles but should be cleaned to canonical_smiles
    inchi = auto()
    inchikey = auto()
    canonical_smiles = auto()
    isomeric_smiles = auto()
    iupac = auto()
    # drugs
    clinical_phase = auto()
    approval_withdrawn = auto()
    target = auto()
    prodrug = auto()
    administration = auto()
    countries_approved = auto()
    otc_prescription = auto()
    # urls
    unichem_url = auto()
    # ids
    cas = auto()
    unii = auto()
    pubchem_cid = auto()
    pubchem_cid_parent = auto()
    chembl_id = auto()
    chebi_id = auto()
    drugbank_id = auto()
    drugcentral_id = auto()
    unichem_id = auto()
    zinc_id = auto()
    schembl_id = auto()
    nmrshiftdb2_id = auto()
    kegg_ligand_id = auto()
    hmdb_id = auto()
    # classifications
    class_results_npclassifier = auto()
    superclass_results_npclassifier = auto()
    pathway_results_npclassifier = auto()
    isglycoside_npclassifier = auto()
    kingdom_classyfire = auto()
    superclass_classyfire = auto()
    class_classyfire = auto()
    subclass_classyfire = auto()
    intermediate_nodes_classyfire = auto()
    alternative_parents_classyfire = auto()
    direct_parent_classyfire = auto()
    molecular_framework_classyfire = auto()
    substituents_classyfire = auto()
    description_classyfire = auto()
    external_descriptors_classyfire = auto()
    ancestors_classyfire = auto()
    predicted_chebi_terms_classyfire = auto()
    predicted_lipidmaps_terms_classyfire = auto()
    classification_version_classyfire = auto()
