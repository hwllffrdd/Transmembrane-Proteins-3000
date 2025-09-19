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
        # Direct accession query
        query_string = f"(accession:{query}) AND (organism_id:{organism_id}) AND (reviewed:{reviewed})"
    else:
        # General search query
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
        
        # Get protein name
        if "proteinDescription" in result and "recommendedName" in result["proteinDescription"]:
            recommended_name = result["proteinDescription"]["recommendedName"]
            if "fullName" in recommended_name and "value" in recommended_name["fullName"]:
                protein_name = recommended_name["fullName"]["value"]
        
        # Get protein sequence length
        sequence_length = 0
        if "sequence" in result and "length" in result["sequence"]:
            sequence_length = result["sequence"]["length"]

        # Get gene name
        if "genes" in result and result["genes"]:
            for gene in result["genes"]:
                if "geneName" in gene:
                    gene_name = gene["geneName"].get("value", "Unknown")
                    break  # Take the first gene name found

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
            feature_info = f"{len(transmembrane_helical_features)} transmembrane helical region{'s' if len(transmembrane_helical_features) != 1 else ''}: "
            feature_info += ", ".join([f"{start}-{end}" for start, end in transmembrane_helical_features])

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
    """Fetch data for multiple accessions in batches from UniProt"""
    results = []
    not_found = []
    
    # Process in smaller batches to avoid overwhelming the API
    batch_size = 10
    
    for i in range(0, len(accessions), batch_size):
        batch_accessions = accessions[i:i+batch_size]
        
        for accession in batch_accessions:
            try:
                # Get data for this accession
                protein_data = get_uniprot_data(accession, is_accession=True)
                
                if protein_data and len(protein_data) > 0:
                    # Process the protein data
                    for result in protein_data:
                        # Format uniprot feature info
                        transmembrane_regions = result.get("transmembrane_regions", [])
                        if not transmembrane_regions:
                            uniprot_count = 0
                            uniprot_regions = ""
                        else:
                            uniprot_count = len(transmembrane_regions)
                            uniprot_regions = ", ".join([f"{start}-{end}" for start, end in transmembrane_regions])
                        
                        # Add uniprot count and regions to result
                        result["uniprot_count"] = uniprot_count
                        result["uniprot_regions"] = uniprot_regions
                        
                        # Add to results
                        results.append(result)
                else:
                    not_found.append(accession)
            except Exception as e:
                print(f"Error processing accession {accession}: {e}")
                not_found.append(accession)
    
    return results, not_found

def process_tm_string(tm_string):
    # Split the string by 'i' and 'o'
    parts = re.split('[io]', tm_string)
    # Filter out empty strings and process each part
    tm_regions = []
    for part in parts:
        if part:
            # Split by '/' and take the last part
            region = part.split('/')[-1]
            # Only add if it contains a '-'
            if '-' in region:
                start, end = map(int, region.split('-'))
                tm_regions.append((start, end))
    return tm_regions


def get_csv_data(accessions=None):
    # Read the CSV file
    df = pd.read_csv('phobius.csv', sep=';')
    
    # Filter by accessions if provided
    if accessions:
        df = df[df['Accession'].isin(accessions)]
    
    processed_data = []
    for _, row in df.iterrows():
        accession = row['Accession']
        signal_peptide = 'Yes' if row['Signal peptide'] == 'Y' else 'No'

        tm_regions = process_tm_string(row['TM string'])
        
        # Count and format regions
        count = len(tm_regions)
        regions_str = ", ".join([f"{start}-{end}" for start, end in tm_regions])

        feature_info = f"{count} transmembrane helical region{'s' if count != 1 else ''}"
        # uncomment the following if you want to report on signal peptides
        # if signal_peptide == 'Yes':
        #    feature_info += " (with a signal peptide)"
        if tm_regions:
            feature_info += ": " + regions_str

        processed_data.append({
            'accession': accession,
            'feature_info': feature_info,
            'data_source': 'CSV',
            'tm_regions': tm_regions,
            'count': count,
            'regions': regions_str
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
                    
                    # Skip if not in the requested accessions
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


def get_tsv_data(accessions=None):
    tsv_data = parse_tsv_data('deeptmhmm.tsv', accessions)
    formatted_results, regions_dict = format_tsv_results(tsv_data)
    return formatted_results, regions_dict


def format_tsv_results(tsv_data):
    formatted_results = {}
    regions_dict = {}
    for accession, data in tsv_data.items():
        tm_count = data['tm_count']
        tm_regions = data.get('tm_regions', [])
        has_signal = data['has_signal']
        
        # Get the parsed regions (tuples of start-end)
        parsed_regions = data.get('tm_regions_parsed', [])
        if not parsed_regions and tm_regions:
            # Convert string regions to tuples if needed
            parsed_regions = []
            for region in tm_regions:
                if '-' in region:
                    start, end = map(int, region.split('-'))
                    parsed_regions.append((start, end))

        if tm_count == 0:
            feature_info = "0 transmembrane helical regions"
        else:
            feature_info = f"{tm_count} transmembrane helical region{'s' if tm_count > 1 else ''}"
            # if has_signal:
            #    feature_info += " (with a signal peptide)"
            feature_info += ": " + ", ".join(tm_regions)

        formatted_results[accession] = feature_info
        regions_dict[accession] = parsed_regions

    return formatted_results, regions_dict


def parse_tmbed_data(file_path, accessions=None):
    """Parse TMBED prediction file and extract transmembrane regions"""
    results = {}
    current_accession = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # Parse FASTA header to extract accession
                # Format: >sp|ACCESSION|PROTEIN_NAME ...
                parts = line.split('|')
                if len(parts) >= 2:
                    current_accession = parts[1]
                    # Skip if not in requested accessions
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
                # Process prediction string
                sequence_length = len(line)
                results[current_accession]['sequence_length'] = sequence_length

                # Find transmembrane regions (H or h)
                tm_regions = []
                signal_regions = []

                i = 0
                while i < len(line):
                    char = line[i]

                    if char in 'Hh':
                        # Start of transmembrane region
                        start = i + 1  # 1-based indexing
                        while i < len(line) and line[i] in 'Hh':
                            i += 1
                        end = i  # 1-based indexing
                        tm_regions.append((start, end))

                    elif char == 'S':
                        # Start of signal peptide region
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
    """Get TMBED data and format it for display"""
    tmbed_data = parse_tmbed_data('tmbed_combined.txt', accessions)
    formatted_results = {}
    regions_dict = {}

    for accession, data in tmbed_data.items():
        tm_regions = data['tm_regions']
        tm_count = len(tm_regions)

        if tm_count == 0:
            feature_info = "0 transmembrane helical regions"
        else:
            regions_str = ", ".join([f"{start}-{end}" for start, end in tm_regions])
            feature_info = f"{tm_count} transmembrane helical region{'s' if tm_count > 1 else ''}: {regions_str}"

        formatted_results[accession] = feature_info
        regions_dict[accession] = tm_regions

    return formatted_results, regions_dict



def process_batch_query(batch_accessions):
    """Process a batch of accessions and return the combined results"""
    # Clean and deduplicate accessions
    accessions = [acc.strip() for acc in batch_accessions if acc.strip()]
    accessions = list(set(accessions))  # Remove duplicates
    
    # Get UniProt data
    uniprot_results, not_found_accessions = get_batch_uniprot_data(accessions)
    
    # Create map of accession to UniProt result for easy lookup
    uniprot_map = {result['accession']: result for result in uniprot_results}
    
    # Get Phobius data for these accessions
    csv_results = get_csv_data(accessions)
    csv_dict = {item['accession']: item for item in csv_results}
    
    # Get DeepTMHMM data for these accessions
    tsv_results, tsv_regions = get_tsv_data(accessions)

    # Get TMBED data for these accessions
    tmbed_results, tmbed_regions = get_tmbed_data(accessions)

    # Create the combined results
    combined_results = []
    
    # Start with all accessions found in UniProt
    for uniprot_item in uniprot_results:
        accession = uniprot_item['accession']
        
        # Get Phobius data
        phobius_item = csv_dict.get(accession, None)
        phobius_feature_info = 'No data'
        phobius_count = 0
        phobius_regions = ''
        if phobius_item:
            phobius_feature_info = phobius_item['feature_info']
            phobius_count = phobius_item['count']
            phobius_regions = phobius_item['regions']
        
        # Get DeepTMHMM data
        deeptmhmm_feature_info = tsv_results.get(accession, 'No data')
        deeptmhmm_regions = tsv_regions.get(accession, [])
        deeptmhmm_count = len(deeptmhmm_regions)
        deeptmhmm_regions_str = ', '.join([f"{start}-{end}" for start, end in deeptmhmm_regions])

        # Get TMBED data
        tmbed_feature_info = tmbed_results.get(accession, 'No data')
        tmbed_regions_list = tmbed_regions.get(accession, [])
        tmbed_count = len(tmbed_regions_list)
        tmbed_regions_str = ', '.join([f"{start}-{end}" for start, end in tmbed_regions_list])
        
        # Add to combined results
        combined_results.append({
            'name': uniprot_item['name'],
            'gene_name': uniprot_item['gene_name'],
            'accession': accession,
            'uniprot_count': uniprot_item['uniprot_count'],
            'uniprot_regions': uniprot_item['uniprot_regions'],
            'phobius_count': phobius_count,
            'phobius_regions': phobius_regions,
            'deeptmhmm_count': deeptmhmm_count,
            'deeptmhmm_regions': deeptmhmm_regions_str,
            'tmbed_count': tmbed_count,
            'tmbed_regions': tmbed_regions_str
        })
    
    return combined_results, not_found_accessions, len(accessions)


def export_to_csv(results, format='csv'):
    """Export results to CSV or TSV format"""
    delimiter = ',' if format == 'csv' else '\t'
    
    output = StringIO()
    writer = csv.writer(output, delimiter=delimiter)
    
    # Write header
    header = [
        'Protein Name', 'Gene Name', 'Accession',
        'UniProt TM Count', 'UniProt TM Regions',
        'Phobius TM Count', 'Phobius TM Regions',
        'DeepTMHMM TM Count', 'DeepTMHMM TM Regions',
        'TMBED TM Count', 'TMBED TM Regions'
    ]
    writer.writerow(header)
    
    # Write data rows
    for item in results:
        row = [
            item['name'], item['gene_name'], item['accession'],
            item['uniprot_count'], item['uniprot_regions'],
            item['phobius_count'], item['phobius_regions'],
            item['deeptmhmm_count'], item['deeptmhmm_regions'],
            item.get('tmbed_count', 0), item.get('tmbed_regions', '')
        ]
        writer.writerow(row)
    
    return output.getvalue()

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Check if it's a single search or batch search
        search_type = request.form.get('search_type', 'single')
        
        if search_type == 'batch':
            # Redirect to batch processing route
            return redirect(url_for('batch'))
        
        # Regular single search
        query = request.form['query']
        uniprot_results = get_uniprot_data(query)
        csv_results = get_csv_data()
        tsv_results, tsv_regions = get_tsv_data()
        tmbed_results, tmbed_regions = get_tmbed_data()
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
            combined_results.append(combined_item)
        return render_template('results.html', results=combined_results, query=query)
    return render_template('index.html')

@app.route('/batch', methods=['POST'])
def batch():
    # Get batch query input
    batch_query = request.form.get('batch_query', '').strip()
    
    if not batch_query:
        return redirect(url_for('index'))
    
    # Split by space or newline
    accessions = re.split(r'[\s\n\r]+', batch_query)
    
    # Generate a unique batch ID
    batch_id = str(uuid.uuid4())
    
    # Store the accession list
    batch_storage[batch_id] = {
        'accessions': accessions,
        'submitted_at': datetime.now().isoformat(),
        'status': 'processing'
    }
    
    # Redirect to processing page
    return render_template('batch_processing.html', batch_id=batch_id, total_accessions=len(accessions))

@app.route('/batch/<batch_id>')
def batch_results(batch_id):
    # Check if batch exists
    if batch_id not in batch_storage:
        return redirect(url_for('index'))
    
    # Get batch data
    batch_data = batch_storage[batch_id]
    
    # If status is processing, process the batch now
    if batch_data['status'] == 'processing':
        accessions = batch_data['accessions']
        
        # Process the batch
        results, not_found, total = process_batch_query(accessions)
        
        # Update batch storage
        batch_data['results'] = results
        batch_data['not_found'] = not_found
        batch_data['total'] = total
        batch_data['processed_at'] = datetime.now().isoformat()
        batch_data['status'] = 'completed'
        
        # Save to file for persistence
        batch_file = os.path.join('batch_results', f"{batch_id}.json")
        with open(batch_file, 'w') as f:
            json.dump(batch_data, f)
    
    # Return the results page
    return render_template('batch_results.html',
                          batch_id=batch_id,
                          results=batch_data['results'],
                          not_found_accessions=batch_data.get('not_found', []),
                          total_accessions=batch_data.get('total', 0))

@app.route('/export/<batch_id>')
def export_batch(batch_id):
    # Get the export format
    format = request.args.get('format', 'csv')
    if format not in ['csv', 'tsv']:
        format = 'csv'
    
    # Check if batch exists
    if batch_id not in batch_storage:
        return redirect(url_for('index'))
    
    # Get batch data
    batch_data = batch_storage[batch_id]
    
    # Check if results are available
    if 'results' not in batch_data or batch_data['status'] != 'completed':
        return redirect(url_for('batch_results', batch_id=batch_id))
    
    # Generate CSV/TSV data
    output_data = export_to_csv(batch_data['results'], format)
    
    # Create response
    mimetype = 'text/csv' if format == 'csv' else 'text/tab-separated-values'
    filename = f"transmembrane_domains_{batch_id}.{format}"
    
    # Return the file as an attachment
    return output_data, 200, {
        'Content-Type': f'{mimetype}; charset=utf-8',
        'Content-Disposition': f'attachment; filename="{filename}"'
    }

if __name__ == '__main__':
    app.run(debug=True)
