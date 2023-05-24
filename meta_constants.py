from enum import auto

from strenum import StrEnum


class MetaColumns(StrEnum):
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
