from flask import Flask, render_template, request
import requests
import pandas as pd
import re

def get_uniprot_data(query):
    organism_id = "9606"
    reviewed = "true"
    query_string = f"((organism_id:{organism_id}) AND (reviewed:{reviewed}) AND {query})"
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
        if "proteinDescription" in result and "recommendedName" in result["proteinDescription"]:
            recommended_name = result["proteinDescription"]["recommendedName"]
            if "fullName" in recommended_name and "value" in recommended_name["fullName"]:
                protein_name = recommended_name["fullName"]["value"]

        transmembrane_helical_features = []
        if "features" in result:
            for feature in result["features"]:
                if (feature.get("type") == "Transmembrane" and
                        feature.get("description") == "Helical"):
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

        results.append({
            "name": protein_name,
            "accession": result.get("primaryAccession", "N/A"),
            "feature_info": feature_info
        })

    return results

def process_tm_string(tm_string):
    parts = re.split('[io]', tm_string)
    tm_regions = [part.split('/')[-1] for part in parts if part]
    return tm_regions


def get_csv_data():
    df = pd.read_csv('phobius.csv', sep=';')

    processed_data = []
    for _, row in df.iterrows():
        accession = row['Accession']
        tm_count = row['TM domains']
        signal_peptide = 'Yes' if row['Signal peptide'] == 'Y' else 'No'

        tm_regions = process_tm_string(row['TM string'])

        feature_info = f"{tm_count} transmembrane helical region{'s' if tm_count != 1 else ''}"
        if signal_peptide == 'Yes':
            feature_info += " (with a signal peptide)"
        feature_info += ": " + ", ".join(tm_regions)

        processed_data.append({
            'accession': accession,
            'feature_info': feature_info,
            'data_source': 'CSV'
        })

    return processed_data


app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        query = request.form['query']
        uniprot_results = get_uniprot_data(query)
        csv_results = get_csv_data()
        csv_dict = {item['accession']: item for item in csv_results}
        combined_results = []
        for uniprot_item in uniprot_results:
            combined_item = uniprot_item.copy()
            csv_item = csv_dict.get(uniprot_item['accession'])
            if csv_item:
                combined_item['phobius_feature_info'] = csv_item['feature_info']
            else:
                combined_item['phobius_feature_info'] = 'No data'
            combined_results.append(combined_item)
        return render_template('results.html', results=combined_results, query=query)
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)


