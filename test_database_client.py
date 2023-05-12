from unittest import TestCase

import chembl_client
import pubchem_client


class TestDatabases(TestCase):

    def test_chembl_id(self):
        from chembl_webresource_client.new_client import new_client as chembl
        comp = chembl.molecule.get('CHEMBL461522')
        assert comp

    def test_chembl_otc_drug(self):
        from chembl_webresource_client.new_client import new_client as chembl
        comp = chembl.molecule.get('CHEMBL112')
        assert comp

    def test_chembl_prescription_drug(self):
        from chembl_webresource_client.new_client import new_client as chembl
        comp = chembl.molecule.get('CHEMBL70')
        assert comp

    def test_get_chembl_mol(self):
        comp = chembl_client.get_chembl_mol(inchikey="BGVLELSCIHASRV-QPEQYQDCSA-N")
        assert comp
        assert comp["max_phase"] == 4
        assert comp["molecule_chembl_id"] == 'CHEMBL3301594'

    def test_get_chembl_mol_by_wrong_inchikey(self):
        comp = chembl_client.get_chembl_mol(inchikey="coco")
        assert not comp

    def test_inchikey_list(self):
        inchikeys = """BGVLELSCIHASRV-QPEQYQDCSA-N
SWMDAPWAQQTBOG-UHFFFAOYSA-N
PZUSGRHVYDQLHR-UHFFFAOYSA-N
KPFZCKDPBMGECB-WGDLNXRISA-N
RYFZBPVMVYTEKZ-KBPBESRZSA-N
JEGHXKRHKHPBJD-UHFFFAOYSA-N
PCMHOSYCWRRHTG-UHFFFAOYSA-N
KFRKRECSIYXARE-HMAPJEAMSA-N
LHSBZAWDPSTOEY-UHFFFAOYSA-N
CYVVJSKZRBZHAV-UNZYHPAISA-N
ZGBAJMQHJDFTQJ-DEOSSOPVSA-N
IZAOBRWCUGOKNH-OAHLLOKOSA-N
HPHUVLMMVZITSG-UHFFFAOYSA-N
RDJGLLICXDHJDY-UHFFFAOYSA-N
RCYPVQCPYKNSTG-UHFFFAOYSA-N
SSZHESNDOMBSRV-UHFFFAOYSA-N
HRDQQHUKUIKFHT-UHFFFAOYSA-N
XGOYIMQSIKSOBS-UHFFFAOYSA-N
JNUGFGAVPBYSHF-UHFFFAOYSA-N
XDLYKKIQACFMJG-WKILWMFISA-N
KFRKRECSIYXARE-HYARGMPZSA-N""".split("\n")
        for inchikey in inchikeys:
            comp = chembl_client.get_chembl_mol(inchikey)
            assert comp
            assert comp["molecule_structures"]["standard_inchi_key"] == inchikey


    def test_get_chembl_mol_by(self):
        comp = chembl_client.get_chembl_mol(inchikey="LVWZTYCIRDMTEY-UHFFFAOYSA-N")
        assert comp

    def test_search_pubchem_by_name(self):
        comp = pubchem_client.search_pubchem_by_name("BGVLELSCIHASRV-QPEQYQDCSA-N")
        assert comp

    def test_search_pubchem_by_structure(self):
        comp = pubchem_client.search_pubchem_by_structure(inchikey="BGVLELSCIHASRV-QPEQYQDCSA-N")
        assert comp


    def test_pubchem_cache(self):
        inchikeys = """BGVLELSCIHASRV-QPEQYQDCSA-N
SWMDAPWAQQTBOG-UHFFFAOYSA-N
PZUSGRHVYDQLHR-UHFFFAOYSA-N
KPFZCKDPBMGECB-WGDLNXRISA-N
RYFZBPVMVYTEKZ-KBPBESRZSA-N
JEGHXKRHKHPBJD-UHFFFAOYSA-N
PCMHOSYCWRRHTG-UHFFFAOYSA-N
KFRKRECSIYXARE-HMAPJEAMSA-N
LHSBZAWDPSTOEY-UHFFFAOYSA-N
CYVVJSKZRBZHAV-UNZYHPAISA-N
ZGBAJMQHJDFTQJ-DEOSSOPVSA-N
IZAOBRWCUGOKNH-OAHLLOKOSA-N
HPHUVLMMVZITSG-UHFFFAOYSA-N
RDJGLLICXDHJDY-UHFFFAOYSA-N
RCYPVQCPYKNSTG-UHFFFAOYSA-N
SSZHESNDOMBSRV-UHFFFAOYSA-N
HRDQQHUKUIKFHT-UHFFFAOYSA-N
XGOYIMQSIKSOBS-UHFFFAOYSA-N
JNUGFGAVPBYSHF-UHFFFAOYSA-N
XDLYKKIQACFMJG-WKILWMFISA-N
KFRKRECSIYXARE-HYARGMPZSA-N""".split("\n")
        search_pubchem(inchikeys)
        print("new chance: this should be chached now \n\n\n\n")
        search_pubchem(inchikeys)

    def test_pubchem_to_df(self):
        inchikeys = """BGVLELSCIHASRV-QPEQYQDCSA-N
SWMDAPWAQQTBOG-UHFFFAOYSA-N
PZUSGRHVYDQLHR-UHFFFAOYSA-N
KPFZCKDPBMGECB-WGDLNXRISA-N
RYFZBPVMVYTEKZ-KBPBESRZSA-N
JEGHXKRHKHPBJD-UHFFFAOYSA-N
PCMHOSYCWRRHTG-UHFFFAOYSA-N
KFRKRECSIYXARE-HMAPJEAMSA-N
LHSBZAWDPSTOEY-UHFFFAOYSA-N
CYVVJSKZRBZHAV-UNZYHPAISA-N
ZGBAJMQHJDFTQJ-DEOSSOPVSA-N
IZAOBRWCUGOKNH-OAHLLOKOSA-N
HPHUVLMMVZITSG-UHFFFAOYSA-N
RDJGLLICXDHJDY-UHFFFAOYSA-N
RCYPVQCPYKNSTG-UHFFFAOYSA-N
SSZHESNDOMBSRV-UHFFFAOYSA-N
HRDQQHUKUIKFHT-UHFFFAOYSA-N
XGOYIMQSIKSOBS-UHFFFAOYSA-N
JNUGFGAVPBYSHF-UHFFFAOYSA-N
XDLYKKIQACFMJG-WKILWMFISA-N
KFRKRECSIYXARE-HYARGMPZSA-N""".split("\n")
        import pubchempy
        df = pubchempy.get_compounds(inchikeys, "inchikey", as_dataframe=True)


def search_pubchem(inchikeys):
    for inchikey in inchikeys:
        comp = pubchem_client.search_pubchem_by_structure(inchikey=inchikey)
        # print(comp)
