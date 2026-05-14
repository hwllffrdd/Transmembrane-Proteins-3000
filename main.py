from flask import Flask, render_template, request, redirect, url_for, send_file
import requests
import pandas as pd
import re
import os
import json
import csv
import uuid
from io import StringIO
from datetime import datetime
from collections import defaultdict

# Create a directory for storing batch results
os.makedirs('batch_results', exist_ok=True)

# Dictionary to store batch results in memory (for development purposes)
batch_storage = {}

def get_uniprot_data(query, is_accession=False):
    organism_id = "9606"
    reviewed = "true"
    
    if is_accession:
        query_string = f"(accession:{query}) AND (organism_id:{organism_id}) AND (reviewed:{reviewed})"
    else:
        query_string = f"((organism_id:{organism_id}) AND (reviewed:{reviewed}) AND ({query}))"
    
    parameters = {
        "format": "json",
        "query": query_string
    }
    response = requests.get("https://rest.uniprot.org/uniprotkb/stream", params=parameters)
    response.raise_for_status()
    data = response.json()

    results = []

    for result in data["results"]:
        protein_name = "Unknown"
        gene_name = "Unknown"
        primary_accession = result.get("primaryAccession", "N/A")
        
        if "proteinDescription" in result and "recommendedName" in result["proteinDescription"]:
            recommended_name = result["proteinDescription"]["recommendedName"]
            if "fullName" in recommended_name and "value" in recommended_name["fullName"]:
                protein_name = recommended_name["fullName"]["value"]
        
        sequence_length = 0
        if "sequence" in result and "length" in result["sequence"]:
            sequence_length = result["sequence"]["length"]

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

        feature_info = ""
        if not result.get("features"):
            feature_info = "No feature data available"
        elif not transmembrane_helical_features:
            feature_info = "0 transmembrane helical regions"
        else:
            tm_percentage = calculate_tm_percentage(transmembrane_helical_features, sequence_length)
            feature_info = f"{len(transmembrane_helical_features)} transmembrane helical region{'s' if len(transmembrane_helical_features) != 1 else ''}: "
            feature_info += ", ".join([f"{start}-{end}" for start, end in transmembrane_helical_features])
            feature_info += f"; comprise {tm_percentage:.1f}% of the sequence"

        if is_accession or (query.lower() in protein_name.lower() or
            query.lower() == gene_name.lower() or
            query.lower() == primary_accession.lower()):
            results.append({
                "name": protein_name,
                "gene_name": gene_name,
                "accession": primary_accession,
                "feature_info": feature_info,
                "sequence_length": sequence_length,
                "transmembrane_regions": transmembrane_helical_features
            })

    return results

def get_batch_uniprot_data(accessions):
    results = []
    not_found = []
    batch_size = 10
    
    for i in range(0, len(accessions), batch_size):
        batch_accessions = accessions[i:i+batch_size]
        
        for accession in batch_accessions:
            try:
                protein_data = get_uniprot_data(accession, is_accession=True)
                
                if protein_data and len(protein_data) > 0:
                    for result in protein_data:
                        transmembrane_regions = result.get("transmembrane_regions", [])
                        if not transmembrane_regions:
                            uniprot_count = 0
                            uniprot_regions = ""
                        else:
                            uniprot_count = len(transmembrane_regions)
                            uniprot_regions = ", ".join([f"{start}-{end}" for start, end in transmembrane_regions])
                        
                        result["uniprot_count"] = uniprot_count
                        result["uniprot_regions"] = uniprot_regions
                        results.append(result)
                else:
                    not_found.append(accession)
            except Exception as e:
                print(f"Error processing accession {accession}: {e}")
                not_found.append(accession)
    
    return results, not_found

def calculate_tm_percentage(tm_regions, sequence_length):
    if not tm_regions or sequence_length == 0:
        return 0.0
    total_tm_residues = sum(end - start + 1 for start, end in tm_regions)
    percentage = (total_tm_residues / sequence_length) * 100
    return percentage

def process_tm_string(tm_string):
    parts = re.split('[io]', tm_string)
    tm_regions = []
    for part in parts:
        if part:
            region = part.split('/')[-1]
            if '-' in region:
                start, end = map(int, region.split('-'))
                tm_regions.append((start, end))
    return tm_regions


def get_csv_data(accessions=None, sequence_length_dict=None):
    df = pd.read_csv('phobius.csv', sep=';')
    
    if accessions:
        df = df[df['Accession'].isin(accessions)]
    
    processed_data = []
    for _, row in df.iterrows():
        accession = row['Accession']
        signal_peptide = 'Yes' if row['Signal peptide'] == 'Y' else 'No'

        tm_regions = process_tm_string(row['TM string'])
        count = len(tm_regions)
        regions_str = ", ".join([f"{start}-{end}" for start, end in tm_regions])

        feature_info = f"{count} transmembrane helical region{'s' if count != 1 else ''}"
        if tm_regions:
            feature_info += ": " + regions_str
            if sequence_length_dict and accession in sequence_length_dict:
                seq_len = sequence_length_dict[accession]
                tm_percentage = calculate_tm_percentage(tm_regions, seq_len)
                feature_info += f"; comprise {tm_percentage:.1f}% of the sequence"

        processed_data.append({
            'accession': accession,
            'feature_info': feature_info,
            'data_source': 'CSV',
            'tm_regions': tm_regions,   # list of (int, int) tuples — used for visualisation
            'count': count,
            'regions': regions_str      # display string — used for table cell
        })

    return processed_data


def parse_tsv_data(file_path, accessions=None):
    results = {}
    current_protein = None
    tm_regions = []
    has_signal = False

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                if 'Number of predicted TMRs:' in line:
                    parts = line.split('|')
                    accession = parts[1]
                    
                    if accessions and accession not in accessions:
                        current_protein = None
                        continue
                        
                    tm_count = int(line.split(':')[-1].strip())
                    current_protein = accession
                    results[current_protein] = {'tm_count': tm_count, 'tm_regions': [], 'has_signal': False, 'tm_regions_parsed': []}
            elif line == '//':
                current_protein = None
                tm_regions = []
                has_signal = False
            elif current_protein and '\t' in line:
                parts = line.split('\t')
                if parts[1] == 'TMhelix':
                    start = int(parts[2])
                    end = int(parts[3])
                    tm_regions.append((start, end))
                    results[current_protein]['tm_regions'].append(f"{start}-{end}")
                    results[current_protein]['tm_regions_parsed'].append((start, end))
                elif parts[1] == 'signal':
                    has_signal = True

            if current_protein:
                results[current_protein]['has_signal'] = has_signal

    return results


def get_tsv_data(accessions=None, sequence_length_dict=None):
    tsv_data = parse_tsv_data('deeptmhmm.tsv', accessions)
    formatted_results, regions_dict = format_tsv_results(tsv_data, sequence_length_dict)
    return formatted_results, regions_dict


def format_tsv_results(tsv_data, sequence_length_dict=None):
    formatted_results = {}
    regions_dict = {}
    for accession, data in tsv_data.items():
        tm_count = data['tm_count']
        tm_regions = data.get('tm_regions', [])
        has_signal = data['has_signal']
        
        parsed_regions = data.get('tm_regions_parsed', [])
        if not parsed_regions and tm_regions:
            parsed_regions = []
            for region in tm_regions:
                if '-' in region:
                    start, end = map(int, region.split('-'))
                    parsed_regions.append((start, end))

        if tm_count == 0:
            feature_info = "0 transmembrane helical regions"
        else:
            feature_info = f"{tm_count} transmembrane helical region{'s' if tm_count > 1 else ''}"
            feature_info += ": " + ", ".join(tm_regions)
            if sequence_length_dict and accession in sequence_length_dict:
                seq_len = sequence_length_dict[accession]
                tm_percentage = calculate_tm_percentage(parsed_regions, seq_len)
                feature_info += f"; comprise {tm_percentage:.1f}% of the sequence"

        formatted_results[accession] = feature_info
        regions_dict[accession] = parsed_regions  # list of (int, int) tuples

    return formatted_results, regions_dict


def parse_tmbed_data(file_path, accessions=None):
    results = {}
    current_accession = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('|')
                if len(parts) >= 2:
                    current_accession = parts[1]
                    if accessions and current_accession not in accessions:
                        current_accession = None
                        continue
                    results[current_accession] = {
                        'tm_regions': [],
                        'signal_regions': [],
                        'sequence_length': 0
                    }
                else:
                    current_accession = None
            elif current_accession and line:
                sequence_length = len(line)
                results[current_accession]['sequence_length'] = sequence_length

                tm_regions = []
                signal_regions = []

                i = 0
                while i < len(line):
                    char = line[i]
                    if char in 'Hh':
                        start = i + 1
                        while i < len(line) and line[i] in 'Hh':
                            i += 1
                        end = i
                        tm_regions.append((start, end))
                    elif char == 'S':
                        start = i + 1
                        while i < len(line) and line[i] == 'S':
                            i += 1
                        end = i
                        signal_regions.append((start, end))
                    else:
                        i += 1

                results[current_accession]['tm_regions'] = tm_regions
                results[current_accession]['signal_regions'] = signal_regions

    return results


def get_tmbed_data(accessions=None):
    tmbed_data = parse_tmbed_data('tmbed_combined.txt', accessions)
    formatted_results = {}
    regions_dict = {}

    for accession, data in tmbed_data.items():
        tm_regions = data['tm_regions']
        tm_count = len(tm_regions)
        sequence_length = data.get('sequence_length', 0)

        if tm_count == 0:
            feature_info = "0 transmembrane helical regions"
        else:
            regions_str = ", ".join([f"{start}-{end}" for start, end in tm_regions])
            feature_info = f"{tm_count} transmembrane helical region{'s' if tm_count > 1 else ''}: {regions_str}"
            if sequence_length > 0:
                tm_percentage = calculate_tm_percentage(tm_regions, sequence_length)
                feature_info += f"; comprise {tm_percentage:.1f}% of the sequence"

        formatted_results[accession] = feature_info
        regions_dict[accession] = tm_regions  # list of (int, int) tuples

    return formatted_results, regions_dict

def parse_topcons_data(file_path, accessions=None):
    results = {}
    current_accession = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                parts = line.split('|')
                if len(parts) >= 2:
                    current_accession = parts[1]
                    if accessions and current_accession not in accessions:
                        current_accession = None
                        continue
                    results[current_accession] = {
                        'tm_regions': [],
                        'signal_regions': [],
                        'sequence_length': 0
                    }
                else:
                    current_accession = None
            elif current_accession and line:
                sequence_length = len(line)
                results[current_accession]['sequence_length'] = sequence_length

                tm_regions = []
                signal_regions = []

                i = 0
                while i < len(line):
                    char = line[i]
                    if char == 'M':
                        start = i + 1
                        while i < len(line) and line[i] == 'M':
                            i += 1
                        end = i
                        tm_regions.append((start, end))
                    elif char == 'S':
                        start = i + 1
                        while i < len(line) and line[i] == 'S':
                            i += 1
                        end = i
                        signal_regions.append((start, end))
                    else:
                        i += 1

                results[current_accession]['tm_regions'] = tm_regions
                results[current_accession]['signal_regions'] = signal_regions

    return results


def get_topcons_data(accessions=None):
    topcons_data = parse_topcons_data('topcons_combined.txt', accessions)
    formatted_results = {}
    regions_dict = {}

    for accession, data in topcons_data.items():
        tm_regions = data['tm_regions']
        tm_count = len(tm_regions)
        sequence_length = data.get('sequence_length', 0)

        if tm_count == 0:
            feature_info = "0 transmembrane helical regions"
        else:
            regions_str = ", ".join([f"{start}-{end}" for start, end in tm_regions])
            feature_info = f"{tm_count} transmembrane helical region{'s' if tm_count > 1 else ''}: {regions_str}"
            if sequence_length > 0:
                tm_percentage = calculate_tm_percentage(tm_regions, sequence_length)
                feature_info += f"; comprise {tm_percentage:.1f}% of the sequence"

        formatted_results[accession] = feature_info
        regions_dict[accession] = tm_regions  # list of (int, int) tuples

    return formatted_results, regions_dict


def process_batch_query(batch_accessions):
    accessions = [acc.strip() for acc in batch_accessions if acc.strip()]
    accessions = list(set(accessions))

    uniprot_results, not_found_accessions = get_batch_uniprot_data(accessions)
    sequence_length_dict = {result['accession']: result['sequence_length'] for result in uniprot_results}

    csv_results = get_csv_data(accessions, sequence_length_dict)
    csv_dict = {item['accession']: item for item in csv_results}

    tsv_results, tsv_regions = get_tsv_data(accessions, sequence_length_dict)
    tmbed_results, tmbed_regions = get_tmbed_data(accessions)
    topcons_results, topcons_regions = get_topcons_data(accessions)

    combined_results = []

    for uniprot_item in uniprot_results:
        accession = uniprot_item['accession']
        seq_len = sequence_length_dict.get(accession, 0)

        # UniProt — tuples already in transmembrane_regions
        uniprot_tm_regions = uniprot_item.get('transmembrane_regions', [])
        uniprot_percentage = f"{calculate_tm_percentage(uniprot_tm_regions, seq_len):.1f}" if uniprot_tm_regions and seq_len else ""

        # Phobius — tm_regions is the list of tuples; regions is the display string
        phobius_item = csv_dict.get(accession)
        phobius_count = 0
        phobius_regions_display = ''
        phobius_regions_list = []
        phobius_percentage = ''
        if phobius_item:
            phobius_count = phobius_item['count']
            phobius_regions_display = phobius_item['regions']   # "23-45, 67-89"
            phobius_regions_list = phobius_item['tm_regions']   # [(23,45), (67,89)]
            if phobius_regions_list and seq_len:
                phobius_percentage = f"{calculate_tm_percentage(phobius_regions_list, seq_len):.1f}"

        # DeepTMHMM
        deeptmhmm_regions_list = tsv_regions.get(accession, [])
        deeptmhmm_count = len(deeptmhmm_regions_list)
        deeptmhmm_regions_display = ', '.join([f"{s}-{e}" for s, e in deeptmhmm_regions_list])
        deeptmhmm_percentage = f"{calculate_tm_percentage(deeptmhmm_regions_list, seq_len):.1f}" if deeptmhmm_regions_list and seq_len else ""

        # TMBED
        tmbed_regions_list = tmbed_regions.get(accession, [])
        tmbed_count = len(tmbed_regions_list)
        tmbed_regions_display = ', '.join([f"{s}-{e}" for s, e in tmbed_regions_list])
        tmbed_percentage = f"{calculate_tm_percentage(tmbed_regions_list, seq_len):.1f}" if tmbed_regions_list and seq_len else ""

        # TOPCONS
        topcons_regions_list = topcons_regions.get(accession, [])
        topcons_count = len(topcons_regions_list)
        topcons_regions_display = ', '.join([f"{s}-{e}" for s, e in topcons_regions_list])
        topcons_percentage = f"{calculate_tm_percentage(topcons_regions_list, seq_len):.1f}" if topcons_regions_list and seq_len else ""

        combined_results.append({
            'name': uniprot_item['name'],
            'gene_name': uniprot_item['gene_name'],
            'accession': accession,
            'sequence_length': seq_len,
            # Display strings for table cells
            'uniprot_count': uniprot_item['uniprot_count'],
            'uniprot_regions': uniprot_item['uniprot_regions'],
            'uniprot_percentage': uniprot_percentage,
            'phobius_count': phobius_count,
            'phobius_regions': phobius_regions_display,
            'phobius_percentage': phobius_percentage,
            'deeptmhmm_count': deeptmhmm_count,
            'deeptmhmm_regions': deeptmhmm_regions_display,
            'deeptmhmm_percentage': deeptmhmm_percentage,
            'tmbed_count': tmbed_count,
            'tmbed_regions': tmbed_regions_display,
            'tmbed_percentage': tmbed_percentage,
            'topcons_count': topcons_count,
            'topcons_regions': topcons_regions_display,
            'topcons_percentage': topcons_percentage,
            # Coordinate lists for visualisation (stored as lists of [start, end])
            'uniprot_regions_list': [list(r) for r in uniprot_tm_regions],
            'phobius_regions_list': [list(r) for r in phobius_regions_list],
            'deeptmhmm_regions_list': [list(r) for r in deeptmhmm_regions_list],
            'tmbed_regions_list': [list(r) for r in tmbed_regions_list],
            'topcons_regions_list': [list(r) for r in topcons_regions_list],
        })

    return combined_results, not_found_accessions, len(accessions)


def export_to_csv(results, format='csv'):
    delimiter = ',' if format == 'csv' else '\t'
    output = StringIO()
    writer = csv.writer(output, delimiter=delimiter)
    
    header = [
        'Protein Name', 'Gene Name', 'Accession',
        'UniProt TM Count', 'UniProt TM Regions', 'UniProt TM %',
        'Phobius TM Count', 'Phobius TM Regions', 'Phobius TM %',
        'DeepTMHMM TM Count', 'DeepTMHMM TM Regions', 'DeepTMHMM TM %',
        'TMBED TM Count', 'TMBED TM Regions', 'TMBED TM %',
        'TOPCONS TM Count', 'TOPCONS TM Regions', 'TOPCONS TM %',
    ]
    writer.writerow(header)
    
    for item in results:
        row = [
            item['name'], item['gene_name'], item['accession'],
            item['uniprot_count'], item['uniprot_regions'], item.get('uniprot_percentage', ''),
            item['phobius_count'], item['phobius_regions'], item.get('phobius_percentage', ''),
            item['deeptmhmm_count'], item['deeptmhmm_regions'], item.get('deeptmhmm_percentage', ''),
            item.get('tmbed_count', 0), item.get('tmbed_regions', ''), item.get('tmbed_percentage', ''),
            item.get('topcons_count', 0), item.get('topcons_regions', ''), item.get('topcons_percentage', ''),
        ]
        writer.writerow(row)
    
    return output.getvalue()

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        search_type = request.form.get('search_type', 'single')
        
        if search_type == 'batch':
            return redirect(url_for('batch'))
        
        query = request.form['query']
        uniprot_results = get_uniprot_data(query)
        
        sequence_length_dict = {result['accession']: result['sequence_length'] for result in uniprot_results}
        
        csv_results = get_csv_data(sequence_length_dict=sequence_length_dict)
        tsv_results, tsv_regions = get_tsv_data(sequence_length_dict=sequence_length_dict)
        tmbed_results, tmbed_regions = get_tmbed_data()
        topcons_results, topcons_regions = get_topcons_data()
        csv_dict = {item['accession']: item for item in csv_results}

        combined_results = []
        for uniprot_item in uniprot_results:
            combined_item = uniprot_item.copy()
            accession = uniprot_item['accession']
            csv_item = csv_dict.get(accession)
            combined_item['phobius_feature_info'] = csv_item['feature_info'] if csv_item else 'No data'
            combined_item['phobius_regions'] = csv_item['tm_regions'] if csv_item else []
            combined_item['deeptmhmm_feature_info'] = tsv_results.get(accession, 'No data')
            combined_item['deeptmhmm_regions'] = tsv_regions.get(accession, [])
            combined_item['tmbed_feature_info'] = tmbed_results.get(accession, 'No data')
            combined_item['tmbed_regions'] = tmbed_regions.get(accession, [])
            combined_item['topcons_feature_info'] = topcons_results.get(accession, 'No data')
            combined_item['topcons_regions'] = topcons_regions.get(accession, [])
            combined_results.append(combined_item)
        return render_template('results.html', results=combined_results, query=query)
    return render_template('index.html')

@app.route('/batch', methods=['POST'])
def batch():
    batch_query = request.form.get('batch_query', '').strip()
    
    if not batch_query:
        return redirect(url_for('index'))
    
    accessions = re.split(r'[\s\n\r]+', batch_query)
    batch_id = str(uuid.uuid4())
    
    batch_storage[batch_id] = {
        'accessions': accessions,
        'submitted_at': datetime.now().isoformat(),
        'status': 'processing'
    }
    
    return render_template('batch_processing.html', batch_id=batch_id, total_accessions=len(accessions))

@app.route('/batch/<batch_id>')
def batch_results(batch_id):
    if batch_id not in batch_storage:
        return redirect(url_for('index'))
    
    batch_data = batch_storage[batch_id]
    
    if batch_data['status'] == 'processing':
        accessions = batch_data['accessions']
        results, not_found, total = process_batch_query(accessions)
        
        batch_data['results'] = results
        batch_data['not_found'] = not_found
        batch_data['total'] = total
        batch_data['processed_at'] = datetime.now().isoformat()
        batch_data['status'] = 'completed'
        
        batch_file = os.path.join('batch_results', f"{batch_id}.json")
        with open(batch_file, 'w') as f:
            json.dump(batch_data, f)
    
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
    
    if batch_id not in batch_storage:
        return redirect(url_for('index'))
    
    batch_data = batch_storage[batch_id]
    
    if 'results' not in batch_data or batch_data['status'] != 'completed':
        return redirect(url_for('batch_results', batch_id=batch_id))
    
    output_data = export_to_csv(batch_data['results'], format)
    mimetype = 'text/csv' if format == 'csv' else 'text/tab-separated-values'
    filename = f"transmembrane_domains_{batch_id}.{format}"
    
    return output_data, 200, {
        'Content-Type': f'{mimetype}; charset=utf-8',
        'Content-Disposition': f'attachment; filename="{filename}"'
    }

@app.route('/stats/<batch_id>')
def batch_stats(batch_id):
    if batch_id not in batch_storage:
        return redirect(url_for('index'))

    batch_data = batch_storage[batch_id]

    if 'results' not in batch_data or batch_data['status'] != 'completed':
        return redirect(url_for('batch_results', batch_id=batch_id))

    results_json = json.dumps(batch_data['results'])

    return render_template('batch_stats.html',
                           batch_id=batch_id,
                           total=len(batch_data['results']),
                           results_json=results_json)

if __name__ == '__main__':
    app.run(debug=True)
