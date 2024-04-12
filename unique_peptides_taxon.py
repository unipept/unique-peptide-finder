import sys
import requests
import re
from tqdm import tqdm

def get_protein_for_taxon(taxon_id):
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A{taxon_id}%29%29"
    response = requests.get(url)
    response.raise_for_status()
    fasta_data = response.text.split('\n>')
    protein_list = [fasta.split('\n', 1)[1].replace('\n', '') for fasta in fasta_data]
    return protein_list

def tryptically_digest_proteins(protein_list):
    peptides = set()
    for protein in protein_list:
        matches = re.split(r'(?<=[KR])(?!P)', protein)
        peptides.update(matches)
    return list(peptides)

def get_unique_peptides_for_taxa(peptides, uniq_taxon):
    unique_peptides = set()
    batch_size=50
    for i in tqdm(range(0, len(peptides), batch_size)):
        batch = peptides[i:i+batch_size]
        url = "http://api.unipept.ugent.be/mpa/pept2data"
        data = {'peptides': batch, 'equate_il': True, 'missed': False}
        response = requests.post(url, json=data)
        response.raise_for_status()
        taxa_data = response.json()
        for item in taxa_data["peptides"]:
            peptide = item["sequence"]
            if item["lca"] == uniq_taxon:
                unique_peptides.add(peptide)
    return unique_peptides

def main():
    if len(sys.argv) < 2:
        print("Please provide a valid NCBI taxon ID as a command line argument.")
        return

    taxon = sys.argv[1]

    proteins = get_protein_for_taxon(taxon)
    peptides = tryptically_digest_proteins(proteins)

    print("Total peptides: ", len(peptides))
    
    uniques = get_unique_peptides_for_taxa(peptides, int(taxon))
    for pep in uniques:
        print(pep)

if __name__ == "__main__":
    main()
