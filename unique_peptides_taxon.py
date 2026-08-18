import argparse
import re
import sys

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

# The public Unipept API. `mpa/pept2data` (used previously) is a private endpoint that
# exists for the Unipept web and desktop applications, and it also computes functional
# annotations that we throw away here.
UNIPEPT_URL = "https://api.unipept.ugent.be/api/v2/pept2lca.json"

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/stream"

# Peptides shorter than this are not present in the Unipept index, so there is no point
# in sending them.
MIN_PEPTIDE_LENGTH = 5

# The API accepts request bodies up to 50 MiB, so we can afford large batches.
BATCH_SIZE = 2000

# UniProt proteomes and Unipept lookups both take a while for large taxa.
REQUEST_TIMEOUT = 300

def create_session():
    """A session that retries on the transient errors a long run is likely to hit.

    Retrying POST is safe here: both endpoints we call are read-only lookups.
    """
    retries = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session

def get_proteins_for_taxon(session, taxon_id):
    params = {"format": "fasta", "query": f"((taxonomy_id:{taxon_id}))"}
    response = session.get(UNIPROT_URL, params=params, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()

    # An unknown taxon, or one without any UniProt entries, returns an empty body.
    if not response.text.strip():
        return []

    fasta_data = response.text.split('\n>')
    return [fasta.split('\n', 1)[1].replace('\n', '') for fasta in fasta_data]

def tryptically_digest_proteins(protein_list):
    peptides = set()
    for protein in protein_list:
        matches = re.split(r'(?<=[KR])(?!P)', protein)
        peptides.update(matches)
    return [peptide for peptide in peptides if len(peptide) >= MIN_PEPTIDE_LENGTH]

def get_unique_peptides_for_taxa(session, peptides, uniq_taxon):
    unique_peptides = set()
    for i in tqdm(range(0, len(peptides), BATCH_SIZE)):
        batch = peptides[i:i+BATCH_SIZE]
        data = {'input': batch, 'equate_il': True}
        response = session.post(UNIPEPT_URL, json=data, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        for item in response.json():
            if item["taxon_id"] == uniq_taxon:
                unique_peptides.add(item["peptide"])
    return unique_peptides

def main():
    parser = argparse.ArgumentParser(
        description="Report all tryptic peptides that are unique to a specific taxon."
    )
    parser.add_argument(
        "taxon_id",
        type=int,
        help="the NCBI taxon ID to find unique peptides for (e.g. 83333)",
    )
    args = parser.parse_args()

    session = create_session()

    try:
        proteins = get_proteins_for_taxon(session, args.taxon_id)
    except requests.RequestException as error:
        print(f"Could not download proteins from UniProt: {error}", file=sys.stderr)
        return 1

    if not proteins:
        print(f"UniProt has no proteins for taxon {args.taxon_id}.", file=sys.stderr)
        return 1

    peptides = tryptically_digest_proteins(proteins)

    # Progress information goes to stderr so that stdout stays a clean list of peptides.
    print(f"Total peptides: {len(peptides)}", file=sys.stderr)

    try:
        uniques = get_unique_peptides_for_taxa(session, peptides, args.taxon_id)
    except requests.RequestException as error:
        print(f"Could not query the Unipept API: {error}", file=sys.stderr)
        return 1

    # Sorting keeps the output stable between runs, which makes it diffable.
    for pep in sorted(uniques):
        print(pep)

    return 0

if __name__ == "__main__":
    sys.exit(main())
