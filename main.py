from flask import Flask, render_template, request, redirect, url_for, jsonify
import requests
from requests.adapters import HTTPAdapter, Retry
import re
import os
import json
import csv
import uuid
import threading
from io import StringIO
from datetime import datetime

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

BATCH_RESULTS_DIR = os.path.join(BASE_DIR, 'batch_results')
os.makedirs(BATCH_RESULTS_DIR, exist_ok=True)

batch_storage = {}
batch_storage_lock = threading.Lock()


# ── Local data files ──────────────────────────────────────────────────────────
# The hosted web app reads these from a PostgreSQL database (populated from the
# very same files by update_imps_database.py). This standalone version parses
# them directly instead, applying the same parsing rules as that loader.

DATA_FILES = {
    'phobius':   os.path.join(BASE_DIR, 'phobius.csv'),
    'deeptmhmm': os.path.join(BASE_DIR, 'deeptmhmm.tsv'),
    'tmbed':     os.path.join(BASE_DIR, 'tmbed_combined.txt'),
    'topcons':   os.path.join(BASE_DIR, 'topcons_combined.txt'),
    'uniprot':   os.path.join(BASE_DIR, 'uniprot_snapshot.json'),
}

PREDICTION_METHODS = ['phobius', 'deeptmhmm', 'tmbed', 'topcons']


def _parse_phobius(file_path):
    """
    Semicolon-delimited CSV.
    Relevant columns: Accession, TM string
    TM string format:  o6-23i302-325o   (regions separated by i/o)
    """
    results = {}
    with open(file_path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter=';')
        for row in reader:
            accession = (row.get('Accession') or '').strip()
            if not accession:
                continue
            regions = results.setdefault(accession, [])
            tm_string = str(row.get('TM string', '')).strip()
            # Split on topology-orientation characters i/o, keep coordinate parts
            for segment in re.split(r'[io]', tm_string):
                for part in segment.split('/'):
                    part = part.strip()
                    if '-' in part:
                        try:
                            start, end = map(int, part.split('-'))
                            regions.append((start, end))
                        except ValueError:
                            pass
    return results


def _parse_deeptmhmm(file_path):
    """
    Custom TSV.
    Header lines: # sp|<accession>|<name> Number of predicted TMRs: N
    Data lines:   <name>\t<type>\t<start>\t<end>
    End-of-record: //
    """
    results = {}
    current = None
    with open(file_path, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#') and 'Number of predicted TMRs:' in line:
                parts = line.split('|')
                if len(parts) >= 2:
                    current = parts[1].strip()
                    results.setdefault(current, [])
            elif line == '//':
                current = None
            elif current and '\t' in line:
                parts = line.split('\t')
                if len(parts) >= 4 and parts[1].strip() == 'TMhelix':
                    try:
                        results[current].append((int(parts[2]), int(parts[3])))
                    except ValueError:
                        pass
    return results


def _parse_fasta_topology(file_path, tm_chars):
    """
    Shared parser for TMBED and TOPCONS.
    FASTA header:  >sp|<accession>|<name> ...
    Next line:     topology string (one char per residue)
    tm_chars:      set of characters that denote TM helices (e.g. {'H','h'} or {'M'})
    Returns (regions per accession, topology-string length per accession).
    """
    results = {}
    lengths = {}
    current = None
    with open(file_path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('|')
                current = parts[1].strip() if len(parts) >= 2 else None
                continue
            if current and line:
                regions = []
                i = 0
                while i < len(line):
                    if line[i] in tm_chars:
                        start = i + 1
                        while i < len(line) and line[i] in tm_chars:
                            i += 1
                        regions.append((start, i))
                    else:
                        i += 1
                results[current] = regions
                lengths[current] = len(line)
                current = None  # each protein is exactly two lines
    return results, lengths


def _load_uniprot_snapshot(file_path):
    """Load the local UniProt snapshot (see download_uniprot_snapshot.py) into
    an accession -> protein dict."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            snapshot = json.load(f)
        return {p['accession']: p for p in snapshot.get('proteins', []) if p.get('accession')}
    except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
        print(f"Warning: could not load {file_path}: {e}")
        return {}


class TMData:
    """Local-file stand-in for the web app's TMDatabase.

    All five data files are parsed once, on first use, and kept in memory.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._loaded = False
        self.proteins = {}      # accession -> {protein_name, gene_name, sequence_length, sequence, tm_regions}
        self.predictions = {}   # method -> {accession -> [(start, end), ...]}

    def _ensure_loaded(self):
        if self._loaded:
            return
        with self._lock:
            if self._loaded:
                return
            print("Loading local data files…")

            def load(method, parser):
                path = DATA_FILES[method]
                if not os.path.exists(path):
                    print(f"Warning: {path} not found")
                    return {}, {}
                result = parser(path)
                return result if isinstance(result, tuple) else (result, {})

            predictions = {}
            predictions['phobius'], _ = load('phobius', _parse_phobius)
            predictions['deeptmhmm'], _ = load('deeptmhmm', _parse_deeptmhmm)
            predictions['tmbed'], tmbed_lengths = load(
                'tmbed', lambda p: _parse_fasta_topology(p, {'H', 'h'}))
            predictions['topcons'], topcons_lengths = load(
                'topcons', lambda p: _parse_fasta_topology(p, {'M'}))

            proteins = {}
            for acc, entry in _load_uniprot_snapshot(DATA_FILES['uniprot']).items():
                proteins[acc] = {
                    'protein_name':    entry.get('protein_name') or None,
                    'gene_name':       entry.get('gene_name') or None,
                    'sequence_length': entry.get('sequence_length') or None,
                    'sequence':        entry.get('sequence') or None,
                    'tm_regions':      [tuple(r) for r in entry.get('tm_regions', []) if len(r) == 2],
                }

            # Mirror the database's proteins table, which also holds accessions
            # that appear only in a prediction file (with no UniProt metadata).
            for method in PREDICTION_METHODS:
                for acc in predictions[method]:
                    if acc not in proteins:
                        proteins[acc] = {'protein_name': None, 'gene_name': None,
                                         'sequence_length': None, 'sequence': None,
                                         'tm_regions': []}
            for acc, protein in proteins.items():
                if not protein['sequence_length']:
                    protein['sequence_length'] = tmbed_lengths.get(acc) or topcons_lengths.get(acc)

            self.proteins = proteins
            self.predictions = predictions
            self._loaded = True
            print(f"Loaded {len(proteins):,} proteins.")

    def get_protein(self, accession):
        self._ensure_loaded()
        return self.proteins.get(accession)

    def get_predictions_for_accessions(self, accessions):
        if not accessions:
            return {}
        self._ensure_loaded()
        organized = {}
        for accession in accessions:
            protein = self.proteins.get(accession)
            if protein is None:
                continue
            methods = {}
            for method in PREDICTION_METHODS:
                regions = self.predictions[method].get(accession)
                if regions:
                    methods[method] = {'transmembrane': sorted(regions)}
            organized[accession] = {
                'protein_name': protein['protein_name'],
                'gene_name': protein['gene_name'],
                'sequence_length': protein['sequence_length'],
                'methods': methods,
            }
        return organized


tm_db = TMData()

# ── Shared requests session with automatic retry (per UniProt docs) ───────────
_uniprot_session = requests.Session()
_uniprot_session.mount("https://", HTTPAdapter(max_retries=Retry(
    total=5,
    backoff_factor=2,
    status_forcelist=[429, 500, 502, 503, 504],
)))
_uniprot_session.headers.update({
    "User-Agent": "IMPS/1.0 (transmembrane protein prediction tool; contact via de.profundis.cz)",
    "Accept": "application/json",
})

# ── UniProt accession format validation ───────────────────────────────────────
# Covers both 6-char ([OPQ][0-9][A-Z0-9]{3}[0-9], e.g. P12345) and
# 10-char ([A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]) formats.
# Rejects Ensembl (ENSG…, ENSP…), RefSeq (NP_…) and any other namespace.
_UNIPROT_ACC_RE = re.compile(
    r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)

def _is_valid_uniprot_accession(acc):
    return bool(_UNIPROT_ACC_RE.match(acc.upper()))


# ── UniProt helpers ───────────────────────────────────────────────────────────

def _parse_uniprot_entry(result):
    """Extract name / gene / length / TM regions from one UniProt JSON entry."""
    protein_name = "Unknown"
    gene_name = "Unknown"
    primary_accession = result.get("primaryAccession", "N/A")

    if "proteinDescription" in result and "recommendedName" in result["proteinDescription"]:
        fn = result["proteinDescription"]["recommendedName"].get("fullName", {})
        protein_name = fn.get("value", "Unknown")

    sequence_length = result.get("sequence", {}).get("length", 0)

    if "genes" in result and result["genes"]:
        for gene in result["genes"]:
            if "geneName" in gene:
                gene_name = gene["geneName"].get("value", "Unknown")
                break

    transmembrane_helical_features = []
    if "features" in result:
        for feature in result["features"]:
            if (feature.get("type") == "Transmembrane" and
                    feature.get("description", "").startswith("Helical")):
                start = feature["location"]["start"]["value"]
                end = feature["location"]["end"]["value"]
                transmembrane_helical_features.append((start, end))

    return primary_accession, protein_name, gene_name, sequence_length, transmembrane_helical_features


def _build_feature_info(tm_regions, sequence_length, has_features):
    """Build the human-readable feature_info string."""
    if not has_features:
        return "No feature data available"
    if not tm_regions:
        return "0 transmembrane helical regions"
    pct = calculate_tm_percentage(tm_regions, sequence_length)
    info = (f"{len(tm_regions)} transmembrane helical region"
            f"{'s' if len(tm_regions) != 1 else ''}: ")
    info += ", ".join(f"{s}-{e}" for s, e in tm_regions)
    info += f"; comprise {pct:.1f}% of the sequence"
    return info


def get_uniprot_data(query, is_accession=False):
    """Single-query search (used for the interactive single-protein search)."""
    organism_id = "9606"
    reviewed = "true"
    if is_accession:
        query_string = f"(accession:{query}) AND (organism_id:{organism_id}) AND (reviewed:{reviewed})"
    else:
        query_string = f"((organism_id:{organism_id}) AND (reviewed:{reviewed}) AND ({query}))"

    response = _uniprot_session.get(
        "https://rest.uniprot.org/uniprotkb/stream",
        params={"format": "json", "query": query_string},
        timeout=15,
    )
    response.raise_for_status()
    data = response.json()

    results = []
    for entry in data["results"]:
        acc, name, gene, length, tm_regions = _parse_uniprot_entry(entry)
        feature_info = _build_feature_info(tm_regions, length, bool(entry.get("features")))

        if is_accession or (
            query.lower() in name.lower() or
            query.lower() == gene.lower() or
            query.lower() == acc.lower()
        ):
            results.append({
                "name": name,
                "gene_name": gene,
                "accession": acc,
                "feature_info": feature_info,
                "sequence_length": length,
                "transmembrane_regions": tm_regions,
            })
    return results


def get_batch_uniprot_data(accessions):
    """
    Retrieve UniProt annotation data for a list of accessions from the local
    uniprot_snapshot.json snapshot (see download_uniprot_snapshot.py), rather
    than hitting the live API.

    Invalid accessions (Ensembl, RefSeq, etc.) are filtered out up front and
    reported in not_found so the rest of the batch is never affected.
    """
    if not accessions:
        return [], []

    # Pre-filter: reject non-UniProt accession formats immediately
    valid   = [a for a in accessions if _is_valid_uniprot_accession(a)]
    invalid = [a for a in accessions if not _is_valid_uniprot_accession(a)]
    if invalid:
        print(f"Skipping {len(invalid)} non-UniProt accession(s): "
              f"{', '.join(invalid[:10])}{'…' if len(invalid) > 10 else ''}")

    if not valid:
        return [], invalid

    results = []
    for acc in valid:
        protein = tm_db.get_protein(acc)
        if protein is None:
            continue
        tm_regions = sorted(protein['tm_regions'])
        sequence_length = protein['sequence_length'] or 0
        results.append({
            'name':                  protein['protein_name'] or 'Unknown',
            'gene_name':             protein['gene_name'] or 'Unknown',
            'accession':             acc,
            'feature_info':          _build_feature_info(tm_regions, sequence_length, True),
            'sequence_length':       sequence_length,
            'transmembrane_regions': tm_regions,
            'uniprot_count':         len(tm_regions),
            'uniprot_regions':       ', '.join(f"{s}-{e}" for s, e in tm_regions),
        })

    found     = {r['accession'] for r in results}
    not_found = invalid + [a for a in valid if a not in found]

    return results, not_found


def get_database_predictions(accessions):
    try:
        return tm_db.get_predictions_for_accessions(accessions)
    except Exception as e:
        print(f"Data error: {e}")
        return {}


def calculate_tm_percentage(tm_regions, sequence_length):
    if not tm_regions or sequence_length == 0:
        return 0.0
    return sum(end - start + 1 for start, end in tm_regions) / sequence_length * 100


def format_method_data(method_data, method_name, sequence_length=0):
    if not method_data:
        return {'feature_info': '0 transmembrane helical regions', 'count': 0,
                'regions': [], 'regions_str': '', 'percentage': ''}

    tm_regions = method_data.get('transmembrane', [])
    count = len(tm_regions)

    if count == 0:
        return {'feature_info': '0 transmembrane helical regions', 'count': 0,
                'regions': [], 'regions_str': '', 'percentage': ''}

    regions_str = ", ".join(f"{s}-{e}" for s, e in tm_regions)
    feature_info = f"{count} transmembrane helical region{'s' if count > 1 else ''}: {regions_str}"
    percentage = ''
    if sequence_length > 0:
        pct = calculate_tm_percentage(tm_regions, sequence_length)
        feature_info += f"; comprise {pct:.1f}% of the sequence"
        percentage = f"{pct:.1f}"

    return {'feature_info': feature_info, 'count': count,
            'regions': tm_regions, 'regions_str': regions_str, 'percentage': percentage}


def process_batch_query(batch_accessions):
    # Dedupe while preserving the order the user entered accessions in — a
    # plain set() would scramble it. Results are re-ordered to match this
    # list at the end, since the lookups below don't preserve it either.
    seen = set()
    accessions = []
    for acc in batch_accessions:
        acc = acc.strip()
        if acc and acc not in seen:
            seen.add(acc)
            accessions.append(acc)

    uniprot_results, not_found_accessions = get_batch_uniprot_data(accessions)
    db_predictions = get_database_predictions(accessions)

    combined_results = []
    for uniprot_item in uniprot_results:
        accession = uniprot_item['accession']
        sequence_length = uniprot_item.get('sequence_length', 0)
        uniprot_tm_regions = uniprot_item.get('transmembrane_regions', [])
        uniprot_percentage = (f"{calculate_tm_percentage(uniprot_tm_regions, sequence_length):.1f}"
                              if uniprot_tm_regions and sequence_length else "")

        db_data = db_predictions.get(accession, {})
        methods = db_data.get('methods', {})

        phobius   = format_method_data(methods.get('phobius'),   'phobius',   sequence_length)
        deeptmhmm = format_method_data(methods.get('deeptmhmm'), 'deeptmhmm', sequence_length)
        tmbed     = format_method_data(methods.get('tmbed'),     'tmbed',     sequence_length)
        topcons   = format_method_data(methods.get('topcons'),   'topcons',   sequence_length)

        combined_results.append({
            'name': uniprot_item['name'],
            'gene_name': uniprot_item['gene_name'],
            'accession': accession,
            'sequence_length': sequence_length,
            'transmembrane_regions': uniprot_tm_regions,
            'feature_info': uniprot_item.get('feature_info', ''),
            'uniprot_count': uniprot_item.get('uniprot_count', 0),
            'uniprot_regions': uniprot_item.get('uniprot_regions', ''),
            'uniprot_percentage': uniprot_percentage,
            'phobius_feature_info':   phobius['feature_info'],
            'phobius_count':          phobius['count'],
            'phobius_regions':        phobius['regions'],
            'phobius_regions_str':    phobius['regions_str'],
            'phobius_percentage':     phobius['percentage'],
            'deeptmhmm_feature_info': deeptmhmm['feature_info'],
            'deeptmhmm_count':        deeptmhmm['count'],
            'deeptmhmm_regions':      deeptmhmm['regions'],
            'deeptmhmm_regions_str':  deeptmhmm['regions_str'],
            'deeptmhmm_percentage':   deeptmhmm['percentage'],
            'tmbed_feature_info':     tmbed['feature_info'],
            'tmbed_count':            tmbed['count'],
            'tmbed_regions':          tmbed['regions'],
            'tmbed_regions_str':      tmbed['regions_str'],
            'tmbed_percentage':       tmbed['percentage'],
            'topcons_feature_info':   topcons['feature_info'],
            'topcons_count':          topcons['count'],
            'topcons_regions':        topcons['regions'],
            'topcons_regions_str':    topcons['regions_str'],
            'topcons_percentage':     topcons['percentage'],
        })

    # Re-sort combined_results to match the order the user entered
    # accessions in (matched case-insensitively, since the user's casing may
    # differ from what's stored in the data files).
    order_index = {}
    for i, acc in enumerate(accessions):
        order_index.setdefault(acc.upper(), i)
    combined_results.sort(key=lambda r: order_index.get(r['accession'].upper(), len(accessions)))

    return combined_results, not_found_accessions, len(accessions)


def _run_batch_in_background(batch_id, accessions):
    """Worker executed in a daemon thread so the request never blocks."""
    batch_file = os.path.join(BATCH_RESULTS_DIR, f"{batch_id}.json")
    try:
        results, not_found, total = process_batch_query(accessions)

        # Serialise region tuples to lists before JSON dump (tuples are not
        # JSON-serialisable; doing it here avoids a second full-copy pass).
        for result in results:
            for field in ['transmembrane_regions', 'phobius_regions',
                          'deeptmhmm_regions', 'tmbed_regions', 'topcons_regions']:
                if field in result:
                    result[field] = [list(r) for r in result[field]]

        payload = {
            'accessions':    accessions,
            'results':       results,
            'not_found':     not_found,
            'total':         total,
            'processed_at':  datetime.now().isoformat(),
            'status':        'completed',
        }

        with open(batch_file + '.tmp', 'w', encoding='utf-8') as f:
            json.dump(payload, f)
        os.replace(batch_file + '.tmp', batch_file)   # atomic swap

        with batch_storage_lock:
            batch_storage[batch_id].update(payload)

    except Exception as e:
        import traceback
        msg = f"{e}\n{traceback.format_exc()}"
        print(f"Background batch worker error for {batch_id}: {msg}")

        error_payload = {
            'status':        'error',
            'error_message': str(e),
        }
        try:
            with open(batch_file, 'w', encoding='utf-8') as f:
                json.dump(error_payload, f)
        except Exception as write_err:
            print(f"Failed to write error file for {batch_id}: {write_err}")

        with batch_storage_lock:
            batch_storage[batch_id].update(error_payload)


def get_single_search_predictions(uniprot_results):
    if not uniprot_results:
        return []
    accessions = [r['accession'] for r in uniprot_results]
    db_predictions = get_database_predictions(accessions)

    combined_results = []
    for uniprot_item in uniprot_results:
        accession = uniprot_item['accession']
        sequence_length = uniprot_item.get('sequence_length', 0)
        uniprot_tm_regions = uniprot_item.get('transmembrane_regions', [])
        uniprot_percentage = (f"{calculate_tm_percentage(uniprot_tm_regions, sequence_length):.1f}"
                              if uniprot_tm_regions and sequence_length else "")

        db_data = db_predictions.get(accession, {})
        methods = db_data.get('methods', {})

        phobius   = format_method_data(methods.get('phobius'),   'phobius',   sequence_length)
        deeptmhmm = format_method_data(methods.get('deeptmhmm'), 'deeptmhmm', sequence_length)
        tmbed     = format_method_data(methods.get('tmbed'),     'tmbed',     sequence_length)
        topcons   = format_method_data(methods.get('topcons'),   'topcons',   sequence_length)

        combined_item = uniprot_item.copy()
        combined_item.update({
            'uniprot_count':          len(uniprot_tm_regions),
            'uniprot_regions':        ', '.join(f"{s}-{e}" for s, e in uniprot_tm_regions),
            'uniprot_percentage':     uniprot_percentage,
            'phobius_feature_info':   phobius['feature_info'],
            'phobius_count':          phobius['count'],
            'phobius_regions':        phobius['regions'],
            'phobius_regions_str':    phobius['regions_str'],
            'phobius_percentage':     phobius['percentage'],
            'deeptmhmm_feature_info': deeptmhmm['feature_info'],
            'deeptmhmm_count':        deeptmhmm['count'],
            'deeptmhmm_regions':      deeptmhmm['regions'],
            'deeptmhmm_regions_str':  deeptmhmm['regions_str'],
            'deeptmhmm_percentage':   deeptmhmm['percentage'],
            'tmbed_feature_info':     tmbed['feature_info'],
            'tmbed_count':            tmbed['count'],
            'tmbed_regions':          tmbed['regions'],
            'tmbed_regions_str':      tmbed['regions_str'],
            'tmbed_percentage':       tmbed['percentage'],
            'topcons_feature_info':   topcons['feature_info'],
            'topcons_count':          topcons['count'],
            'topcons_regions':        topcons['regions'],
            'topcons_regions_str':    topcons['regions_str'],
            'topcons_percentage':     topcons['percentage'],
        })
        combined_results.append(combined_item)

    return combined_results


_SINGLE_SEARCH_METHOD_FIELDS = [
    'uniprot_count', 'phobius_count', 'deeptmhmm_count', 'tmbed_count', 'topcons_count',
]


def _build_single_search_summary(query, combined_results):
    """One-line summary shown above the single-search results table.

    A syntactically valid accession means the user was looking up one
    specific protein, so the phrasing is singular/definite. Anything else is
    a name/gene search, which can legitimately match several proteins — for
    those, X counts proteins with at least one method predicting membrane
    (a union across methods, so "up to X" rather than "at least X": X is an
    upper bound on how many are actually membrane proteins, since methods
    disagree on individual proteins).
    """
    total = len(combined_results)

    if _is_valid_uniprot_accession(query.strip()):
        if total == 0:
            return 'No protein found for this accession.'
        methods_predicting = sum(
            1 for field in _SINGLE_SEARCH_METHOD_FIELDS if combined_results[0].get(field, 0) > 0
        )
        return (f'Your accession is a membrane protein as predicted by '
                f'{methods_predicting} of {len(_SINGLE_SEARCH_METHOD_FIELDS)} methods.')

    if total == 0:
        return 'Your query found 0 proteins.'
    imp_count = sum(
        1 for r in combined_results
        if any(r.get(field, 0) > 0 for field in _SINGLE_SEARCH_METHOD_FIELDS)
    )
    protein_word = 'protein' if total == 1 else 'proteins'
    membrane_phrase = 'is a membrane protein' if imp_count == 1 else 'are membrane proteins'
    return (f'Your query found {total} {protein_word}, '
            f'up to {imp_count} of them {membrane_phrase}.')


def get_batch_sequences(accessions):
    """Bulk-fetch full amino-acid sequences for a set of accessions, keyed by
    accession. Batch results don't carry full sequences (only sequence_length)
    to keep the per-batch JSON file light — see protein_sequence() — so
    export, which needs the residues for each predicted region, looks them up
    separately.
    """
    sequences = {}
    for acc in {a for a in accessions if a}:
        protein = tm_db.get_protein(acc)
        if protein and protein['sequence']:
            sequences[acc] = protein['sequence']
    return sequences


def _regions_to_sequences(regions, sequence):
    """Comma-separated residue sequence for each (start, end) region, in
    order — used for the per-method 'Sequences' export columns below.
    Malformed entries (not a (start, end) pair) are skipped, matching the
    same defensive filter _fix_region_strings() applies elsewhere."""
    if not sequence or not regions:
        return ''
    return ', '.join(sequence[r[0] - 1:r[1]] for r in regions if len(r) == 2)


def export_to_csv(results, format='csv'):
    delimiter = ',' if format == 'csv' else '\t'
    output = StringIO()

    sequences = get_batch_sequences([item['accession'] for item in results])

    def method_seqs(item, regions_field):
        return _regions_to_sequences(item.get(regions_field, []), sequences.get(item['accession'], ''))

    if format == 'csv':
        writer = csv.writer(output, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_ALL)
        header = [
            'Accession', 'Protein Name', 'Gene Name',
            'UniProt TM Count', 'UniProt TM Regions', 'UniProt TM Sequences', 'UniProt TM %',
            'Phobius TM Count', 'Phobius TM Regions', 'Phobius TM Sequences', 'Phobius TM %',
            'DeepTMHMM TM Count', 'DeepTMHMM TM Regions', 'DeepTMHMM TM Sequences', 'DeepTMHMM TM %',
            'TMBED TM Count', 'TMBED TM Regions', 'TMBED TM Sequences', 'TMBED TM %',
            'TOPCONS TM Count', 'TOPCONS TM Regions', 'TOPCONS TM Sequences', 'TOPCONS TM %',
        ]
        writer.writerow(header)
        for item in results:
            writer.writerow([
                item['accession'], item['name'], item['gene_name'],
                item.get('uniprot_count', 0), item.get('uniprot_regions', ''), method_seqs(item, 'transmembrane_regions'), item.get('uniprot_percentage', ''),
                item.get('phobius_count', 0), item.get('phobius_regions_str', ''), method_seqs(item, 'phobius_regions'), item.get('phobius_percentage', ''),
                item.get('deeptmhmm_count', 0), item.get('deeptmhmm_regions_str', ''), method_seqs(item, 'deeptmhmm_regions'), item.get('deeptmhmm_percentage', ''),
                item.get('tmbed_count', 0), item.get('tmbed_regions_str', ''), method_seqs(item, 'tmbed_regions'), item.get('tmbed_percentage', ''),
                item.get('topcons_count', 0), item.get('topcons_regions_str', ''), method_seqs(item, 'topcons_regions'), item.get('topcons_percentage', ''),
            ])
    else:
        writer = csv.writer(output, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
        header = [
            'Protein Name', 'Gene Name', 'Accession',
            'UniProt TM Count', 'UniProt TM Regions', 'UniProt TM Sequences', 'UniProt TM %',
            'Phobius TM Count', 'Phobius TM Regions', 'Phobius TM Sequences', 'Phobius TM %',
            'DeepTMHMM TM Count', 'DeepTMHMM TM Regions', 'DeepTMHMM TM Sequences', 'DeepTMHMM TM %',
            'TMBED TM Count', 'TMBED TM Regions', 'TMBED TM Sequences', 'TMBED TM %',
            'TOPCONS TM Count', 'TOPCONS TM Regions', 'TOPCONS TM Sequences', 'TOPCONS TM %',
        ]
        writer.writerow(header)
        for item in results:
            writer.writerow([
                item['name'], item['gene_name'], item['accession'],
                item.get('uniprot_count', 0), item.get('uniprot_regions', ''), method_seqs(item, 'transmembrane_regions'), item.get('uniprot_percentage', ''),
                item.get('phobius_count', 0), item.get('phobius_regions_str', ''), method_seqs(item, 'phobius_regions'), item.get('phobius_percentage', ''),
                item.get('deeptmhmm_count', 0), item.get('deeptmhmm_regions_str', ''), method_seqs(item, 'deeptmhmm_regions'), item.get('deeptmhmm_percentage', ''),
                item.get('tmbed_count', 0), item.get('tmbed_regions_str', ''), method_seqs(item, 'tmbed_regions'), item.get('tmbed_percentage', ''),
                item.get('topcons_count', 0), item.get('topcons_regions_str', ''), method_seqs(item, 'topcons_regions'), item.get('topcons_percentage', ''),
            ])

    return output.getvalue()


def _ensure_batch_loaded(batch_id):
    """Return the batch record, reading from disk when needed.

    batch_storage is in-memory only, so after a restart a batch is only known
    through its JSON file. A cached entry still marked 'processing' is
    re-read on every call, in case the worker has finished it since.
    """
    batch_file = os.path.join(BATCH_RESULTS_DIR, f"{batch_id}.json")

    with batch_storage_lock:
        cached = batch_storage.get(batch_id)

        if cached is None or cached.get('status') == 'processing':
            if os.path.exists(batch_file):
                try:
                    with open(batch_file, 'r', encoding='utf-8') as f:
                        batch_storage[batch_id] = json.load(f)
                except (ValueError, OSError):
                    # Most likely a partially-written file. Keep whatever we
                    # had; the next poll retries.
                    return cached
            elif cached is None:
                return None

        return batch_storage[batch_id]


def _fix_region_strings(batch_data):
    """Convert region data to list-of-tuples as expected by templates.
    Regions are stored as JSON lists [[start,end],...] so just convert
    each inner list to a tuple. Also handles the legacy string format in
    case old batch files are still on disk.
    """
    if 'results' not in batch_data:
        return
    for result in batch_data['results']:
        for field in ['transmembrane_regions', 'phobius_regions', 'deeptmhmm_regions',
                      'tmbed_regions', 'topcons_regions']:
            val = result.get(field)
            if val is None:
                result[field] = []
            elif isinstance(val, list):
                # Normal case: [[start,end],...] from JSON
                result[field] = [tuple(r) for r in val if len(r) == 2]
            elif isinstance(val, str):
                # Legacy: stringified Python list "[(36, 62), ...]"
                try:
                    region_str = val.strip('[]')
                    regions = []
                    for t_str in re.findall(r'\(\d+,\s*\d+\)', region_str):
                        nums = re.findall(r'\d+', t_str)
                        if len(nums) == 2:
                            regions.append((int(nums[0]), int(nums[1])))
                    result[field] = regions
                except Exception:
                    result[field] = []


# ── Routes ────────────────────────────────────────────────────────────────────

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        search_type = request.form.get('search_type', 'single')

        if search_type == 'batch':
            batch_query = request.form.get('batch_query', '').strip()
            if not batch_query:
                return redirect(url_for('index'))

            accessions = re.split(r'[\s\n\r]+', batch_query)
            # Some users paste "protein group" accessions in the format
            # ACCESSION1;ACCESSION2;ACCESSION3 (as some proteomics tools
            # export them) rather than one accession per token. We only
            # support single accessions, so keep just the first member of
            # each group instead of rejecting the whole line.
            accessions = [a.split(';', 1)[0] for a in accessions]
            accessions = [a for a in accessions if a]
            batch_id = str(uuid.uuid4())

            initial_state = {
                'accessions': accessions,
                'submitted_at': datetime.now().isoformat(),
                'status': 'processing',
            }

            with batch_storage_lock:
                batch_storage[batch_id] = initial_state

            # Write a 'processing' marker to disk immediately, so the status
            # endpoint can still find the batch if in-memory state is lost.
            try:
                marker = os.path.join(BATCH_RESULTS_DIR, f"{batch_id}.json")
                with open(marker, 'w', encoding='utf-8') as f:
                    json.dump(initial_state, f)
            except OSError as e:
                print(f"Could not write processing marker for {batch_id}: {e}")

            t = threading.Thread(
                target=_run_batch_in_background,
                args=(batch_id, accessions),
                daemon=True,
            )
            t.start()

            return render_template('batch_processing.html',
                                   batch_id=batch_id,
                                   total_accessions=len(accessions))

        query = request.form['query']
        uniprot_results = get_uniprot_data(query)
        combined_results = get_single_search_predictions(uniprot_results)
        query_summary = _build_single_search_summary(query, combined_results)
        return render_template('results.html', results=combined_results, query=query,
                               query_summary=query_summary)

    return render_template('index.html')


@app.route('/status/<batch_id>')
def batch_status(batch_id):
    batch_data = _ensure_batch_loaded(batch_id)
    if batch_data is None:
        return jsonify({'status': 'not_found'}), 404

    status = batch_data.get('status', 'processing')
    response = {'status': status}
    if status == 'completed':
        response['result_count'] = len(batch_data.get('results', []))
    elif status == 'error':
        response['error_message'] = batch_data.get('error_message', 'Unknown error')
    return jsonify(response)


@app.route('/batch/<batch_id>')
def batch_results(batch_id):
    batch_data = _ensure_batch_loaded(batch_id)
    if batch_data is None:
        return redirect(url_for('index'))

    status = batch_data.get('status', 'processing')

    if status == 'processing':
        return render_template('batch_processing.html',
                               batch_id=batch_id,
                               total_accessions=len(batch_data.get('accessions', [])))

    if status == 'error':
        return render_template('batch_processing.html',
                               batch_id=batch_id,
                               total_accessions=len(batch_data.get('accessions', [])),
                               error=batch_data.get('error_message', 'An unknown error occurred.'))

    _fix_region_strings(batch_data)

    return render_template('batch_results.html',
                           batch_id=batch_id,
                           results=batch_data['results'],
                           not_found_accessions=batch_data.get('not_found', []),
                           total_accessions=batch_data.get('total', 0))


@app.route('/export/<batch_id>')
def export_batch(batch_id):
    format = request.args.get('format', 'csv')
    if format not in ['csv', 'tsv']:
        format = 'csv'

    batch_data = _ensure_batch_loaded(batch_id)
    if batch_data is None or batch_data.get('status') != 'completed':
        return redirect(url_for('index'))

    _fix_region_strings(batch_data)

    output_data = export_to_csv(batch_data['results'], format)
    if format == 'csv':
        output_data = '﻿' + output_data

    mimetype = 'text/csv' if format == 'csv' else 'text/tab-separated-values'
    filename = f"transmembrane_domains_{batch_id}.{format}"
    return output_data, 200, {
        'Content-Type': f'{mimetype}; charset=utf-8',
        'Content-Disposition': f'attachment; filename="{filename}"'
    }


@app.route('/sequence/<accession>')
def protein_sequence(accession):
    """
    Return the amino-acid sequence for one accession as JSON.

    Used by the results visualizer to show the residues of each predicted
    transmembrane region on hover. Fetched lazily (only when a user expands a
    row) so that result pages stay light regardless of batch size.
    """
    accession = accession.strip().upper()

    if not _is_valid_uniprot_accession(accession):
        return jsonify({'error': 'invalid accession format'}), 400

    try:
        protein = tm_db.get_protein(accession)
    except Exception as e:
        print(f"Data error in protein_sequence({accession}): {e}")
        return jsonify({'error': 'data error'}), 500

    if protein is None or not protein['sequence']:
        return jsonify({'error': 'sequence not available'}), 404

    return jsonify({
        'accession': accession,
        'sequence':  protein['sequence'],
        'length':    protein['sequence_length'],
    })


@app.route('/stats/<batch_id>')
def batch_stats(batch_id):
    batch_data = _ensure_batch_loaded(batch_id)
    if batch_data is None or batch_data.get('status') != 'completed':
        return redirect(url_for('index'))

    _fix_region_strings(batch_data)
    results_json = json.dumps(batch_data['results'])

    return render_template('batch_stats.html',
                           batch_id=batch_id,
                           total=len(batch_data['results']),
                           results_json=results_json)


if __name__ == '__main__':
    app.run(debug=True)
