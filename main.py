from flask import Flask, render_template, request
import requests
import pandas as pd
import re


def get_uniprot_data(query):
    organism_id = "9606"
    reviewed = "true"
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
            feature_info = f"{len(transmembrane_helical_features)} transmembrane helical regions: "
            feature_info += ", ".join([f"{start}-{end}" for start, end in transmembrane_helical_features])

        # Check for matches: exact for gene name and accession, partial for protein name
        if (query.lower() in protein_name.lower() or
            query.lower() == gene_name.lower() or
            query.lower() == primary_accession.lower()):
            results.append({
                "name": protein_name,
                "gene_name": gene_name,
                "accession": primary_accession,
                "feature_info": feature_info
            })

    return results

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
                tm_regions.append(region)
    return tm_regions


def get_csv_data():
    # Read the CSV file
    df = pd.read_csv('phobius.csv', sep=';')

    processed_data = []
    for _, row in df.iterrows():
        accession = row['Accession']
        signal_peptide = 'Yes' if row['Signal peptide'] == 'Y' else 'No'

        tm_regions = process_tm_string(row['TM string'])

        feature_info = f"{len(tm_regions)} transmembrane helical region{'s' if len(tm_regions) != 1 else ''}"
        if signal_peptide == 'Yes':
            feature_info += " (with a signal peptide)"
        if tm_regions:
            feature_info += ": " + ", ".join(tm_regions)

        processed_data.append({
            'accession': accession,
            'feature_info': feature_info,
            'data_source': 'CSV'
        })

    return processed_data


def parse_tsv_data(file_path):
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
                    tm_count = int(line.split(':')[-1].strip())
                    current_protein = accession
                    results[current_protein] = {'tm_count': tm_count, 'tm_regions': [], 'has_signal': False}
            elif line == '//':
                current_protein = None
                tm_regions = []
                has_signal = False
            elif current_protein and '\t' in line:
                parts = line.split('\t')
                if parts[1] == 'TMhelix':
                    tm_regions.append(f"{parts[2]}-{parts[3]}")
                elif parts[1] == 'signal':
                    has_signal = True

            if current_protein:
                results[current_protein]['tm_regions'] = tm_regions
                results[current_protein]['has_signal'] = has_signal

    return results


def get_tsv_data():
    tsv_data = parse_tsv_data('deeptmhmm.tsv')
    return format_tsv_results(tsv_data)


def format_tsv_results(tsv_data):
    formatted_results = {}
    for accession, data in tsv_data.items():
        tm_count = data['tm_count']
        tm_regions = data['tm_regions']
        has_signal = data['has_signal']

        if tm_count == 0:
            feature_info = "0 transmembrane helical regions"
        else:
            feature_info = f"{tm_count} transmembrane helical region{'s' if tm_count > 1 else ''}"
            if has_signal:
                feature_info += " (with a signal peptide)"
            feature_info += ": " + ", ".join(tm_regions)

        formatted_results[accession] = feature_info

    return formatted_results


app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        query = request.form['query']
        uniprot_results = get_uniprot_data(query)
        csv_results = get_csv_data()
        tsv_results = get_tsv_data()
        csv_dict = {item['accession']: item for item in csv_results}

        combined_results = []
        for uniprot_item in uniprot_results:
            combined_item = uniprot_item.copy()
            accession = uniprot_item['accession']
            csv_item = csv_dict.get(accession)
            combined_item['csv_feature_info'] = csv_item['feature_info'] if csv_item else 'No data'
            combined_item['tsv_feature_info'] = tsv_results.get(accession, 'No data')
            combined_results.append(combined_item)
        return render_template('results.html', results=combined_results, query=query)
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)