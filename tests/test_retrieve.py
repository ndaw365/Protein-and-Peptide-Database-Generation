def test_lookup_returns_depmap_and_peptide(index):
    (record,) = index.lookup(gene="PERM1", protein_change="P750Q")
    assert record["source_ids"][:2] == ["depmap:row1", "peptide:row1"]
    assert record["peptide"]["peptide_sequence_fwd"] == "RGPVPSFAFSQNDMCLVFVAFATWAVRTSDQHTPDAWKTALLANVGTISAIR"
    assert (record["peptide"]["peptide_start_fwd"], record["peptide"]["peptide_end_fwd"]) == ("720", "771")
    assert record["peptide"]["Peptide_Header_Fwd"].endswith("GN=PERM1")


def test_three_letter_and_lowercase_queries(index):
    assert index.lookup(gene="perm1", protein_change="p.Pro750Gln")[0]["source_ids"][0] == "depmap:row1"


def test_clinvar_match_paths(index):
    nras = index.lookup(gene="NRAS", protein_change="G12C")[0]["clinvar"]
    assert [(c["variation_id"], c["match"]) for c in nras] == [("900001", "same_allele")]  # GRCh37 copy skipped
    tp53 = index.lookup(gene="TP53", protein_change="R306*")[0]["clinvar"]
    assert [(c["variation_id"], c["match"]) for c in tp53] == [("900002", "same_rsid")]
    perm1 = index.lookup(gene="PERM1", protein_change="P750Q")[0]["clinvar"]
    assert [(c["variation_id"], c["match"]) for c in perm1] == [("900003", "same_protein_change")]
    assert index.conn.execute("SELECT COUNT(*) FROM clinvar WHERE gene='NOTAMOLT4GENE'").fetchone()[0] == 0


def test_lookup_by_rsid_and_uniprot(index):
    assert index.lookup(rsid="rs121913250")[0]["gene"] == "NRAS"
    genes = {r["gene"] for r in index.lookup(uniprot="Q5SV97-1")}
    assert genes == {"PERM1"}


def test_ms_detection_flags_whether_peptide_spans_the_variant(index):
    ppp1ca = index.lookup(gene="PPP1CA", protein_change="D242N")[0]["ms_detection"]
    assert ppp1ca[0]["peptide"] == "AHQVVEDGYEFFAK" and ppp1ca[0]["covers_variant"] is False
    evl = index.lookup(gene="EVL", protein_change="D20N")[0]["ms_detection"]
    assert evl[0]["peptide"] == "ASVMVYDNTSK" and evl[0]["covers_variant"] is True


def test_recovered_splice_region_missense_is_in_peptide_db(index):
    # missense_variant&splice_region_variant rows were dropped before the R filter fix
    (record,) = index.lookup(gene="STK11", protein_change="M125I")
    assert record["variant_info"] == "missense_variant&splice_region_variant"
    assert record["peptide"] is not None


def test_dropped_variant_reason(index):
    (record,) = index.lookup(gene="HELZ2", protein_change="R563L")
    assert record["peptide"] is None
    assert record["dropped_reason"].startswith("Reference mismatch")


def test_retrieve_exact_and_search_modes(index):
    exact = index.retrieve("Is TP53 p.Arg306Ter in ClinVar?")
    assert exact["mode"] == "exact"
    assert exact["records"][0]["source_ids"][0] == "depmap:row1535"

    search = index.retrieve("variant peptide detected by mass spectrometry")
    assert search["mode"] == "search"
    assert all(any(d["covers_variant"] for d in r["ms_detection"]) for r in search["records"])


def test_filter_variants(index):
    hotspots = index.filter_variants(flag="Hotspot", limit=100)
    assert hotspots and all(r["depmap"]["Hotspot"] == "True" for r in hotspots)
    stops = index.filter_variants(gene="TP53", variant_type="stop_gained")
    assert [r["protein_change"] for r in stops] == ["p.R306Ter"]
