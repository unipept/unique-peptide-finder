# Unique Peptide Finder

Reports all tryptic peptides that are unique to a given NCBI taxon.

A peptide counts as unique to a taxon when its lowest common ancestor, as computed by
[Unipept](https://unipept.ugent.be) over all of UniProtKB, is exactly that taxon. In
other words: every protein in UniProt that contains this peptide belongs to this taxon,
and the taxon is the most specific statement you can make about it.

## Requirements

* Python 3
* The `requests` and `tqdm` packages:

```
pip3 install -r requirements.txt
```

## Usage

```
python3 unique_peptides_taxon.py <NCBI_TAXON_ID>
```

Peptides are written to stdout, one per line; progress goes to stderr. Redirecting gives
you a clean file:

```
python3 unique_peptides_taxon.py 83333 > ecoli_k12.txt
```

## Example

```
> python3 unique_peptides_taxon.py 83333
Total peptides: 83537
100%|███████████████████████████████████████████| 42/42 [00:46<00:00,  1.10s/it]
AAALDELIPGLLSEYNR
AADSGSIVLTHLSK
AAHAPFISHPAEFCHLLVALK
AALADFIVDNR
AALAEMVSGDELVIEFDCTQATEAIPQWAAEEGHAITDYQQIGDAAWSITVQK
AALLSSQDLSVYSMNTPGFIPGIDFSDHLNYWQHDIPAIMITDTAFYR
AANPNGFLVYFSDHGEEVYDTPPHK
AAPNFNIAEDFR
...
```

Taxon 83333 is *Escherichia coli* K-12. It yields 2529 unique peptides out of 83537
tryptic peptides, and takes about a minute end to end.

## How it works

1. Download every UniProtKB protein for the taxon (and its descendants) from the
   [UniProt REST API](https://www.uniprot.org/help/api).
2. Digest those proteins in silico with trypsin: cut after every K or R that is not
   followed by a P.
3. Look up the lowest common ancestor of each peptide through the public Unipept
   [`pept2lca`](https://unipept.ugent.be/apidocs/pept2lca) endpoint, in batches.
4. Keep the peptides whose LCA is the requested taxon.

## Notes

* Isoleucine and leucine are treated as equal, which is the Unipept default. A peptide
  that is unique only because of an I/L difference is not reported.
* Peptides shorter than 5 amino acids are skipped, because the Unipept index does not
  contain them.
* Uniqueness is judged against whatever UniProtKB knows today, so results shift as
  UniProt and the Unipept index are updated.
* Peptides whose LCA is a *descendant* of the requested taxon (a strain below it, say)
  are unique to that taxon too, but are not reported. This only matters for taxa whose
  proteins are annotated at a more specific rank than the one you asked for.
