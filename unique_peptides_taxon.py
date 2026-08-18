import sys
import requests
import re
from tqdm import tqdm

# The public Unipept API. `mpa/pept2data` (used previously) is a private endpoint that
# exists for the Unipept web and desktop applications, and it also computes functional
# annotations that we throw away here.
UNIPEPT_URL = "https://api.unipept.ugent.be/api/v2/pept2lca.json"

# Peptides shorter than this are not present in the Unipept index, so there is no point
# in sending them.
MIN_PEPTIDE_LENGTH = 5

# The API accepts request bodies up to 50 MiB, so we can afford large batches.
BATCH_SIZE = 2000

# UniProt proteomes and Unipept lookups both take a while for large taxa.
REQUEST_TIMEOUT = 300

def get_protein_for_taxon(taxon_id):
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A{taxon_id}%29%29"
    response = requests.get(url, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    fasta_data = response.text.split('\n>')
    protein_list = [fasta.split('\n', 1)[1].replace('\n', '') for fasta in fasta_data]
    return protein_list

def tryptically_digest_proteins(protein_list):
    peptides = set()
    for protein in protein_list:
        matches = re.split(r'(?<=[KR])(?!P)', protein)
        peptides.update(matches)
    return [peptide for peptide in peptides if len(peptide) >= MIN_PEPTIDE_LENGTH]

def get_unique_peptides_for_taxa(peptides, uniq_taxon):
    unique_peptides = set()
    for i in tqdm(range(0, len(peptides), BATCH_SIZE)):
        batch = peptides[i:i+BATCH_SIZE]
        data = {'input': batch, 'equate_il': True}
        response = requests.post(UNIPEPT_URL, json=data, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        for item in response.json():
            if item["taxon_id"] == uniq_taxon:
                unique_peptides.add(item["peptide"])
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
