from enum import auto

from strenum import StrEnum


class MetaColumns(StrEnum):
    # run dates
    date_pubchem_parent_cid_search = auto()
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
    filename = auto()  # usually generated in sequence creation
    unique_sample_id = (
        auto()
    )  # generated from plate_id and well_location and a library id
    well_location = auto()
    plate_id = auto()
    library_id = auto()
    # calculated properties
    monoisotopic_mass = auto()
    # compound specific
    compound_name = auto()
    input_name = auto()
    synonyms = auto()
    formula = auto()
    smiles = auto()  # input smiles but should be cleaned to canonical_smiles
    smarts = auto()
    inchi = auto()
    inchikey = auto()
    split_inchikey = auto()
    canonical_smiles = auto()
    isomeric_smiles = auto()
    iupac = auto()
    logp = auto()
    # drugs
    first_approval = auto()
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
    input_pubchem_cid = auto()
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
    npclassifier_class_results = auto()
    npclassifier_superclass_results = auto()
    npclassifier_pathway_results = auto()
    npclassifier_isglycoside = auto()
    classyfire_kingdom = auto()
    classyfire_superclass = auto()
    classyfire_class = auto()
    classyfire_subclass = auto()
    classyfire_intermediate_nodes = auto()
    classyfire_alternative_parents = auto()
    classyfire_direct_parent = auto()
    classyfire_molecular_framework = auto()
    classyfire_substituents = auto()
    classyfire_description = auto()
    classyfire_external_descriptors = auto()
    classyfire_ancestors = auto()
    classyfire_predicted_chebi_terms = auto()
    classyfire_predicted_lipidmaps_terms = auto()
    classyfire_classification_version = auto()
