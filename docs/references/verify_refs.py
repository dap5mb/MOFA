"""Check every reference row against PubMed (and Crossref for DOI-only rows).

Usage: python3 verify_refs.py streams/*_refs.csv references.csv ...
Writes <input>.verified.csv next to each input and prints a summary.

Per row:
  id_check    PMID record matches title (similarity >= 0.90), first-author surname and year
              (publication or e-publication year).
  doi_check   DOI given equals the DOI PubMed lists for the PMID (or, for DOI-only rows,
              Crossref resolves it and its title matches).
  quote_check evidence_quote occurs in the PubMed abstract (exact after normalising case,
              whitespace and punctuation; "fuzzy" if a >= 0.90 near-match exists).
  status      OK when id_check and doi_check pass and quote_check is exact or fuzzy;
              REVIEW when only the quote fails or is missing (claim may rest on full text);
              FAIL otherwise. FAIL rows must not be cited.
"""

import csv
import difflib
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

CACHE = Path(__file__).parent / "data" / "pubmed_cache"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def fetch(url, key):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / key
    if path.exists():
        return path.read_text()
    req = urllib.request.Request(url, headers={"User-Agent": "apoe-discussion-verify/1.0"})
    for attempt in range(6):
        time.sleep(0.4 * 2 ** attempt)  # NCBI allows 3 requests/s without an API key; back off on 429
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                text = r.read().decode("utf-8")
            break
        except urllib.error.HTTPError as e:
            if e.code != 429 or attempt == 5:
                raise
    path.write_text(text)
    return text


def norm(s):
    s = unicodedata.normalize("NFKD", s or "").encode("ascii", "ignore").decode()
    s = re.sub(r"<[^>]+>", " ", s)  # PubMed/Crossref titles may carry <i> tags
    s = re.sub(r"[^a-z0-9 ]", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def similar(a, b):
    return difflib.SequenceMatcher(None, norm(a), norm(b)).ratio()


def pubmed_summary(pmid):
    data = json.loads(fetch(f"{EUTILS}esummary.fcgi?db=pubmed&retmode=json&id={pmid}", f"sum_{pmid}.json"))
    rec = data.get("result", {}).get(pmid)
    if not rec or "error" in rec:
        return None
    return rec


def pubmed_abstract(pmid):
    xml = fetch(f"{EUTILS}efetch.fcgi?db=pubmed&retmode=xml&id={pmid}", f"abs_{pmid}.xml")
    root = ET.fromstring(xml)
    return " ".join("".join(el.itertext()) for el in root.iter("AbstractText"))


def crossref(doi):
    key = "cr_" + re.sub(r"[^A-Za-z0-9]", "_", doi) + ".json"
    try:
        data = json.loads(fetch("https://api.crossref.org/works/" + urllib.parse.quote(doi), key))
        return data["message"]
    except Exception:
        return None


def years(rec):
    return {m for d in (rec.get("pubdate", ""), rec.get("epubdate", "")) for m in re.findall(r"\d{4}", d)}


def quote_check(quote, abstract):
    q, a = norm(quote), norm(abstract)
    if not q:
        return "missing"
    if q in a:
        return "exact"
    m = difflib.SequenceMatcher(None, q, a).find_longest_match(0, len(q), 0, len(a))
    window = a[max(0, m.b - m.a): m.b - m.a + len(q) + 20]
    best = max((difflib.SequenceMatcher(None, q, window[i:i + len(q)]).ratio()
                for i in range(0, max(1, len(window) - len(q) + 1))), default=0)
    return "fuzzy" if best >= 0.90 else "absent"


def check(row):
    pmid = (row.get("pmid") or "").strip()
    doi = (row.get("doi") or "").strip().lower().removeprefix("https://doi.org/")
    out = {"id_check": "", "doi_check": "", "quote_check": "", "pubmed_title": "", "detail": ""}
    if pmid:
        rec = pubmed_summary(pmid)
        if rec is None:
            out.update(id_check="pmid_not_found")
            return out
        out["pubmed_title"] = rec.get("title", "")
        problems = []
        if similar(row.get("title"), rec.get("title")) < 0.90:
            problems.append("title")
        first = (rec.get("authors") or [{}])[0].get("name", "")
        if norm(row.get("first_author")) not in norm(first):
            problems.append(f"first_author (PubMed: {first})")
        if str(row.get("year", "")).strip() not in years(rec):
            problems.append(f"year (PubMed: {rec.get('pubdate')})")
        out["id_check"] = "ok" if not problems else "mismatch: " + "; ".join(problems)
        pm_doi = next((a["value"].lower() for a in rec.get("articleids", []) if a.get("idtype") == "doi"), "")
        if not doi:
            out["doi_check"] = "no_doi_given"
        elif not pm_doi:
            out["doi_check"] = "pubmed_has_no_doi"
        else:
            out["doi_check"] = "ok" if doi == pm_doi else f"mismatch (PubMed: {pm_doi})"
        out["quote_check"] = quote_check(row.get("evidence_quote"), pubmed_abstract(pmid))
    elif doi:
        msg = crossref(doi)
        if msg is None:
            out.update(id_check="doi_not_found")
            return out
        title = (msg.get("title") or [""])[0]
        out["pubmed_title"] = title
        out["id_check"] = "ok (crossref)" if similar(row.get("title"), title) >= 0.90 else "mismatch: title"
        out["doi_check"] = "ok"
        out["quote_check"] = quote_check(row.get("evidence_quote"), msg.get("abstract", ""))
    else:
        out["id_check"] = "no_identifier"
    return out


def status(o):
    id_ok = o["id_check"].startswith("ok")
    doi_ok = o["doi_check"] in ("ok", "no_doi_given", "pubmed_has_no_doi")
    if not (id_ok and doi_ok):
        return "FAIL"
    return "OK" if o["quote_check"] in ("exact", "fuzzy") else "REVIEW"


def main(paths):
    for p in map(Path, paths):
        with p.open(newline="") as f:
            rows = list(csv.DictReader(f))
        counts = {}
        for row in rows:
            try:
                o = check(row)
            except Exception as e:  # network or parse error: never silently pass
                o = {"id_check": f"error: {e}", "doi_check": "", "quote_check": "", "pubmed_title": "", "detail": ""}
            o["status"] = status(o)
            row.update(o)
            counts[o["status"]] = counts.get(o["status"], 0) + 1
        out = p.with_suffix(".verified.csv")
        with out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["status"])
            w.writeheader()
            w.writerows(rows)
        print(f"{p.name}: {len(rows)} rows {counts}")
        for row in rows:
            if row["status"] != "OK":
                print(f"  {row['status']:6} {row.get('ref_id')}: id={row['id_check']} doi={row['doi_check']} quote={row['quote_check']}")


if __name__ == "__main__":
    main(sys.argv[1:])
