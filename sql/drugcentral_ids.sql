select * from identifier where id_type = 'ChEMBL_ID' and identifier.identifier = 'CHEMBL3289398';
select * from identifier where id_type = 'UNII' and identifier.identifier = 'RIK029813G';
select * from identifier where struct_id = '677' and identifier.identifier = 'CHEMBL278255';
select * from identifier where id_type = 'ChEMBL_ID' and identifier.identifier = 'CHEMBL298734';


select structures.id, cas_reg_no, structures.name, structures.inchikey, inchi, smiles, substance_unii, structures.stem,
       string_agg(distinct stem.definition, ',') as "stem_definition",
       mrdef, status,
       string_agg(distinct concat_ws(':', i.id_type, i.identifier), ';') as "database_id",
       string_agg(distinct concat_ws(':', a.type, approval::text), ';') as "date_of_approval",
       string_agg(distinct a.type, ',') as "administration",
       string_agg(distinct orphan::text, ',') as "orphan",
       string_agg(distinct syn.name, ',') as "synonyms",
       string_agg(distinct (s2p.parent_id::text), ',') as "parent_id",
       string_agg(distinct pc.type, ',') as "pharma_type",
       string_agg(distinct pc.name, ',') as "pharma_class",
       string_agg(distinct s2atc.atc_code, ',') as "atc",
       string_agg(distinct concat_ws('', ai.quantity::text, ai.unit), ',') as "dosage",
       string_agg(distinct str.strength, ';') as "strength",
       string_agg(distinct concat_ws('', atc_ddd.ddd::text, atc_ddd.unit_type), ',') as "who_defined_daily_dose",
       string_agg(distinct concat_ws(':', o.relationship_name, o.concept_name), ';') as "indication; contraindication; off_label",
--        string_agg(distinct faers.meddra_name, ',') as "faers",
       string_agg(distinct v.species, ';') as "animal_species",
       string_agg(distinct concat_ws(':', v.relationship_type, v.concept_name), ';') as "indication"
from structures
left join inn_stem stem on structures.stem = stem.stem
left join approval a on structures.id = a.struct_id
left join synonyms syn on structures.id = syn.id
left join identifier i on structures.id = i.struct_id
left join active_ingredient ai on structures.id = ai.struct_id
left join pharma_class pc on structures.id = pc.struct_id
left join struct2parent s2p on structures.id = s2p.struct_id
left join struct2atc s2atc on structures.id = s2atc.struct_id
left join atc_ddd on structures.id = atc_ddd.struct_id
left join omop_relationship o on structures.id = o.struct_id
left join struct2obprod str on structures.id = str.struct_id
left join vetomop v on a.struct_id = v.struct_id
-- left join faers on structures.id = faers.struct_id
-- condition defined in code
where i.id_type = 'ChEMBL_ID' and i.identifier = 'CHEMBL278255'
GROUP BY structures.id, cas_reg_no, structures.name, inchi, smiles, structures.stem, mrdef, structures.inchikey, status, substance_unii;


--
-- select structures.id, cas_reg_no, structures.name, structures.inchikey, inchi, smiles, substance_unii, structures.stem,
--        string_agg(distinct stem.definition, ',') as "stem_definition",
--        mrdef, status,
--        string_agg(distinct concat_ws(':', i.id_type, i.identifier), ';') as "database_id",
--        string_agg(distinct concat_ws(':', a.type, approval::text), ';') as "date_of_approval",
--        string_agg(distinct a.type, ',') as "administration",
--        string_agg(distinct orphan::text, ',') as "orphan",
--        string_agg(distinct syn.name, ',') as "synonyms",
--        string_agg(distinct (s2p.parent_id::text), ',') as "parent_id",
--        string_agg(distinct pc.type, ',') as "pharma_type",
--        string_agg(distinct pc.name, ',') as "pharma_class",
--        string_agg(distinct s2atc.atc_code, ',') as "atc",
--        string_agg(distinct concat_ws('', ai.quantity::text, ai.unit), ',') as "dosage",
--        string_agg(distinct str.strength, ';') as "strength",
--        string_agg(distinct concat_ws('', atc_ddd.ddd::text, atc_ddd.unit_type), ',') as "who_defined_daily_dose",
--        string_agg(distinct concat_ws(':', o.relationship_name, o.concept_name), ';') as "indication; contraindication; off_label",
-- --        string_agg(distinct faers.meddra_name, ',') as "faers",
--        string_agg(distinct v.species, ';') as "animal_species",
--        string_agg(distinct concat_ws(':', v.relationship_type, v.concept_name), ';') as "indication"
-- from structures
-- left join inn_stem stem on structures.stem = stem.stem
-- left join approval a on structures.id = a.struct_id
-- left join synonyms syn on structures.id = syn.id
-- left join identifier i on structures.id = i.struct_id
-- left join active_ingredient ai on structures.id = ai.struct_id
-- left join pharma_class pc on structures.id = pc.struct_id
-- left join struct2parent s2p on structures.id = s2p.struct_id
-- left join struct2atc s2atc on structures.id = s2atc.struct_id
-- left join atc_ddd on structures.id = atc_ddd.struct_id
-- left join omop_relationship o on structures.id = o.struct_id
-- left join struct2obprod str on structures.id = str.struct_id
-- left join vetomop v on a.struct_id = v.struct_id
-- -- left join faers on structures.id = faers.struct_id
-- where (i.id_type = 'ChEMBL_ID' and i.identifier = 'CHEMBL2109334')
-- -- or (i.id_type = 'UNII' and i.identifier = 'XZA0KB1BDQ')
-- GROUP BY structures.id, cas_reg_no, structures.name, inchi, smiles, structures.stem, mrdef, structures.inchikey, status, substance_unii
