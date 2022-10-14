from unittest import TestCase

import database_client


class TestDatabases(TestCase):

    def test_get_chembl_mol_by_inchikey(self):
        comp = database_client.get_chembl_mol_by_inchikey("SWMDAPWAQQTBOG-UHFFFAOYSA-N")
        assert comp
        assert comp["max_phase"] == 4
        assert comp["molecule_chembl_id"] == 'CHEMBL3301594'

    def test_get_chembl_mol_by_wrong_inchikey(self):
        comp = database_client.get_chembl_mol_by_inchikey("coco")
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
            comp = database_client.get_chembl_mol_by_inchikey(inchikey)
            assert comp
            assert comp["molecule_structures"]["standard_inchi_key"] == inchikey


    def test_get_chembl_mol_by_inchikey(self):
        comp = database_client.get_chembl_mol_by_inchikey("LVWZTYCIRDMTEY-UHFFFAOYSA-N")
        assert comp

    def test_search_pubchem_by_name(self):
        comp = database_client.search_pubchem_by_name("BGVLELSCIHASRV-QPEQYQDCSA-N")
        assert comp

    def test_search_pubchem_by_structure(self):
        comp = database_client.search_pubchem_by_structure(inchikey="BGVLELSCIHASRV-QPEQYQDCSA-N")
        assert comp

    def test_get_openfda_information(self):
        comp = database_client.get_openfda_information(r"PARACETAMOL")
        assert comp



