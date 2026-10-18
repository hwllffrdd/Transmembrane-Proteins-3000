#!/usr/bin/env python3
"""
download_uniprot_snapshot.py
============================
Downloads all reviewed human proteins with transmembrane annotations from the
UniProt REST API and saves them as a local JSON snapshot file.

Run this locally (or on the server) before running update_imps_database.py.
The output file is then uploaded to /srv/app/data/ alongside the other data
files and consumed by the update script.

Usage
-----
  python download_uniprot_snapshot.py
  python download_uniprot_snapshot.py --out /some/other/path/uniprot_snapshot.json
  python download_uniprot_snapshot.py --all   # include proteins with 0 TM helices too

What gets downloaded
--------------------
  Query: reviewed:true AND organism_id:9606
  Fields: accession, id (entry name, e.g. NCAM1_HUMAN), protein_name, gene_names,
          sequence, ft_transmem

By default (without --all) only proteins that UniProt annotates with at least
one helical transmembrane region are stored. This keeps the file small and
matches what the batch tool actually needs.

Output format
-------------
  A JSON file: list of objects, one per protein:
  {
    "accession":       "P12345",
    "entry_name":      "OPSD_HUMAN",
    "protein_name":    "Rhodopsin",
    "gene_name":       "RHO",
    "sequence_length": 348,
    "sequence":        "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIML...",
    "tm_regions":      [[36, 62], [68, 95], ...]   // [start, end] 1-based inclusive
  }

The full amino-acid sequence is included so the web tool can display the actual
residues of each predicted transmembrane region. This makes the snapshot roughly
5x larger (~14 MB instead of ~2.8 MB), which is still small enough to upload
comfortably.

Resumption
----------
  If the download is interrupted, re-running will restart from scratch (the
  UniProt /stream endpoint does not support range requests). For the full human
  SwissProt (~20,000 proteins) the download takes about 30–90 seconds depending
  on your connection.
"""

import argparse
import json
import sys
import time
from datetime import datetime

import requests
from requests.adapters import HTTPAdapter, Retry

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

QUERY  = "reviewed:true AND organism_id:9606"
FIELDS = "accession,id,protein_name,gene_names,sequence,ft_transmem"

DEFAULT_OUT = "uniprot_snapshot.json"

SESSION = requests.Session()
SESSION.mount("https://", HTTPAdapter(max_retries=Retry(
    total=8,
    backoff_factor=2,
    status_forcelist=[429, 500, 502, 503, 504],
)))
SESSION.headers.update({
    "User-Agent": "IMPS/1.0 (transmembrane protein prediction tool; contact via de.profundis.cz)",
    "Accept": "application/json",
})


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_entry(entry, tm_only):
    """
    Extract the fields we care about from one UniProt JSON entry.
    Returns None if tm_only=True and the protein has no helical TM regions.
    """
    accession = entry.get("primaryAccession", "")
    entry_name = entry.get("uniProtkbId", "")   # e.g. "NCAM1_HUMAN"

    protein_name = ""
    pd = entry.get("proteinDescription", {})
    rn = pd.get("recommendedName", {})
    protein_name = rn.get("fullName", {}).get("value", "") if rn else ""
    if not protein_name:
        # Fall back to submitted name
        sn = pd.get("submissionNames", [{}])
        protein_name = sn[0].get("fullName", {}).get("value", "") if sn else ""

    gene_name = ""
    for gene in entry.get("genes", []):
        if "geneName" in gene:
            gene_name = gene["geneName"].get("value", "")
            break

    sequence_obj = entry.get("sequence", {})
    sequence_length = sequence_obj.get("length", 0)
    sequence = sequence_obj.get("value", "")

    tm_regions = []
    for feature in entry.get("features", []):
        if (feature.get("type") == "Transmembrane" and
                feature.get("description", "").startswith("Helical")):
            start = feature["location"]["start"]["value"]
            end   = feature["location"]["end"]["value"]
            tm_regions.append([start, end])

    if tm_only and not tm_regions:
        return None

    return {
        "accession":       accession,
        "entry_name":      entry_name,
        "protein_name":    protein_name,
        "gene_name":       gene_name,
        "sequence_length": sequence_length,
        "sequence":        sequence,
        "tm_regions":      tm_regions,
    }


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download(out_path, tm_only):
    print("=" * 65)
    print("UniProt snapshot downloader")
    print(f"Query : {QUERY}")
    print(f"Output: {out_path}")
    print(f"Mode  : {'TM-annotated proteins only' if tm_only else 'all reviewed human proteins'}")
    print("=" * 65)

    params = {
        "format":   "json",
        "query":    QUERY,
        "fields":   FIELDS,
        "size":     500,      # page size — stream delivers everything in one response
                              # but size controls internal chunking on UniProt's side
    }

    print("Connecting to UniProt…")
    t0 = time.time()

    try:
        response = SESSION.get(UNIPROT_STREAM, params=params, timeout=120, stream=True)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"\n✗ Request failed: {e}")
        sys.exit(1)

    print("Receiving data…")

    try:
        data = response.json()
    except ValueError as e:
        print(f"\n✗ Failed to parse JSON response: {e}")
        sys.exit(1)

    raw_entries = data.get("results", [])
    elapsed_dl = time.time() - t0
    print(f"  Received {len(raw_entries):,} entries in {elapsed_dl:.1f}s")

    # Parse
    print("Parsing entries…")
    proteins = []
    skipped  = 0
    for entry in raw_entries:
        parsed = parse_entry(entry, tm_only)
        if parsed is None:
            skipped += 1
        else:
            proteins.append(parsed)

    tm_annotated = sum(1 for p in proteins if p["tm_regions"])
    with_sequence = sum(1 for p in proteins if p["sequence"])
    with_entry_name = sum(1 for p in proteins if p["entry_name"])
    total_residues = sum(len(p["sequence"]) for p in proteins)
    print(f"  Kept    : {len(proteins):,} proteins")
    if tm_only:
        print(f"  Skipped : {skipped:,} proteins (no helical TM annotation)")
    print(f"  With TM : {tm_annotated:,}")
    print(f"  With seq: {with_sequence:,}  ({total_residues:,} residues total)")
    print(f"  With entry name: {with_entry_name:,}")

    # Write
    snapshot = {
        "_meta": {
            "generated_at":    datetime.utcnow().isoformat() + "Z",
            "query":           QUERY,
            "fields":          FIELDS,
            "total_entries":   len(raw_entries),
            "stored_proteins": len(proteins),
            "tm_only":         tm_only,
        },
        "proteins": proteins,
    }

    print(f"Writing {out_path}…")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(snapshot, fh, ensure_ascii=False, separators=(",", ":"))

    size_mb = __import__("os").path.getsize(out_path) / 1_048_576
    print(f"\n✓ Done.  {len(proteins):,} proteins written to {out_path} ({size_mb:.1f} MB)")
    print(f"  Total time: {time.time() - t0:.1f}s")
    print()
    print("Next step: upload the snapshot to the server and run update_imps_database.py uniprot")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download UniProt TM snapshot for IMPS.")
    parser.add_argument("--out",  default=DEFAULT_OUT,
                        help=f"Output file path (default: {DEFAULT_OUT})")
    parser.add_argument("--all",  action="store_true",
                        help="Include proteins with zero TM helices (larger file)")
    args = parser.parse_args()

    download(out_path=args.out, tm_only=not args.all)
