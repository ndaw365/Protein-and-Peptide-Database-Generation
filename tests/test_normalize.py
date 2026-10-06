import pytest

from rag_assistant.ingest import parse_variant_header
from rag_assistant.normalize import (chrom, find_protein_changes, find_rsids, protein_change_key,
                                     rsid, uniprot_base)


@pytest.mark.parametrize("text, key", [
    ("p.P750Q", "P750Q"),
    ("P750Q", "P750Q"),
    ("p.Pro750Gln", "P750Q"),
    ("p.(Gly12Cys)", "G12C"),
    ("p.R306Ter", "R306*"),
    ("NM_000546.6(TP53):c.916C>T (p.Arg306Ter)", "R306*"),
    ("p.K267RfsTer9", "K267fs"),
    ("p.Lys267fs", "K267fs"),
    ("p.P1142_P1143delinsQT", None),
    ("", None),
    (None, None),
])
def test_protein_change_key(text, key):
    assert protein_change_key(text) == key


def test_gene_symbols_are_not_protein_changes():
    assert find_protein_changes("MAP2K1 ATP1A2 TP53 CDC45") == []
    assert find_protein_changes("Is MAP2K1 D67N or NRAS p.Gly12Cys pathogenic?") == ["D67N", "G12C"]


@pytest.mark.parametrize("value, expected", [
    ("121913250.0", "rs121913250"),  # DepMap stores dbSNP IDs as floats
    ("rs121913250", "rs121913250"),
    ("121913250", "rs121913250"),
    ("-1", None),                    # ClinVar's "no rsID"
    ("", None),
])
def test_rsid(value, expected):
    assert rsid(value) == expected


def test_identifier_helpers():
    assert find_rsids("see rs121913250 and RS1") == ["rs121913250", "rs1"]
    assert uniprot_base("Q5SV97-1") == "Q5SV97"
    assert chrom("chr1") == "1" and chrom("chrM") == "MT" and chrom("X") == "X"


def test_parse_variant_header():
    assert parse_variant_header("Fwd_spD242N|P62136-1|PPP1CA_D242N OS=Homo sapiens") == ("PPP1CA", "D242N")
    assert parse_variant_header("Fwd_spA292T|Q9UHD8-1|SEPTIN9_A292T") == ("SEPTIN9", "A292T")
    assert parse_variant_header("Rev_spD242N|P62136-1|PPP1CA_D242N") is None  # decoy
