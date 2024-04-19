# Unique Peptide Finder
This script queries the Unipept API and reports all tryptic peptides that are unique to a specific taxon.

<img width="1250" alt="image" src="https://github.com/unipept/unique-peptide-finder/assets/9608686/ebf60781-83b7-475b-81bc-548e89cb1d62">

## Requirements
* Python3 needs to be installed locally
* `tdqm` module (do `pip3 install tdqm` if this is not yet installed)

## Usage
**Note:** Make sure that you're running these commands from the root of this repository's folder.

```python
python3 unique_peptides_taxon.py <NCBI_ID>
```

## Example
```
> python3 unique_peptides_taxon.py 83333
AAGLATGNVSTAELQDATPAALEAHVTSR
LTTLSHTSEGHR
VANLEAQLAEAQTR
TVLTWTVLP
LIIFFIGTVFVLAALIPMQQVGVEK
IYDYIDIIIEDYENTK
GDFMNIFKPISYIASLAPR
VTNIWHGR
...
```
