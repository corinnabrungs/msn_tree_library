# Example for exogenous
select *
from compound c
where c.chemfont_id = 'CFc000000116';

## example for endogenous
select *
from compound c
where c.name = 'Aldosterone';

select m.chemfont_id, m.mol_name, oh.html_content
from molecules m
         left outer join ontology_htmls oh on m.mol_id = oh.molecule_id
where m.chemfont_id = 'CF000000178';






