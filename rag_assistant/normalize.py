"""Normalise the identifiers that link DepMap, ClinVar and the peptide database."""

import re

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Sec": "U", "Pyl": "O", "Ter": "*",
}
_AA3 = "|".join(AA3_TO_1)

# Substitutions, nonsense and frameshift notations in 1- or 3-letter form, with or
# without "p." and ClinVar's parentheses, e.g. P750Q, p.Pro750Gln, p.(Arg306Ter).
# The look-arounds stop gene symbols such as MAP2K1 being read as "P2K".
_CHANGE_RE = re.compile(
    rf"(?<![A-Za-z0-9])(?P<prefix>p\.)?\(?(?P<ref>{_AA3}|[A-Z])(?P<pos>\d+)"
    rf"(?P<alt>(?:{_AA3}|[A-Z])?fs(?:Ter\d+|\*\d+)?|{_AA3}|[A-Z*])\)?(?![A-Za-z0-9])"
)
_RSID_RE = re.compile(r"\brs(\d+)\b", re.IGNORECASE)
_UNIPROT_RE = re.compile(
    r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\b"
)


def _aa(code):
    return AA3_TO_1.get(code, code)


def _key(m):
    ref, pos, alt = _aa(m["ref"]), m["pos"], m["alt"]
    if "fs" in alt:
        return f"{ref}{pos}fs"
    return f"{ref}{pos}{_aa(alt)}"


def find_protein_changes(text):
    """All protein-change keys in free text, e.g. "P750Q", "R306*", "K267fs".

    Frameshifts are reduced to "<ref><pos>fs" because DepMap and ClinVar spell the
    new residue and stop distance differently.
    """
    if not text:
        return []
    return list(dict.fromkeys(_key(m) for m in _CHANGE_RE.finditer(str(text))))


def protein_change_key(text):
    """The canonical key for one protein change, preferring an explicit "p." form."""
    if not text:
        return None
    matches = list(_CHANGE_RE.finditer(str(text)))
    if not matches:
        return None
    explicit = [m for m in matches if m["prefix"]]
    return _key((explicit or matches)[0])


def rsid(value):
    """Normalise dbSNP IDs: DepMap stores them as floats ("121913250.0")."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {"-1", "na", "NA", "nan"}:
        return None
    m = _RSID_RE.search(text)
    if m:
        return f"rs{m[1]}"
    try:
        number = int(float(text))
    except ValueError:
        return None
    return f"rs{number}" if number > 0 else None


def find_rsids(text):
    return [f"rs{n}" for n in dict.fromkeys(_RSID_RE.findall(text))]


def uniprot_base(accession):
    """Strip the isoform suffix: "Q5SV97-1" -> "Q5SV97"."""
    if not accession:
        return None
    return str(accession).strip().split("-")[0] or None


def find_uniprot_ids(text):
    return list(dict.fromkeys(m.group(1) for m in _UNIPROT_RE.finditer(text)))


def chrom(value):
    """"chr1" / "1" / "chrM" / "MT" -> "1" / "1" / "MT" / "MT"."""
    if value is None:
        return None
    text = str(value).strip().removeprefix("chr")
    return "MT" if text in {"M", "MT"} else text or None
