# SDfileParser
Pure-Python implementation of SD file parser.

# Usage
## from command line.

```sh
# typical use. Only json mode is implemented.
$ python3 sdfparser.py tojson data/PubChem_compound_text_egfr_records.sdf 
# display options.
$ python3 sdfparser.py tojson -h
# set maximum number of molecule for readling
$ python3 sdfparser.py tojson data/PubChem_compound_text_egfr_records.sdf --maxNumOfMol 10
# skip 10 molecule from the head record.
$ python3 sdfparser.py tojson data/PubChem_compound_text_egfr_records.sdf --numOfSkippedMol 5
```

## As module.
```python
import sdfparser as sdf
import json

file1 = 'data/PubChem_compound_text_egfr_records.sdf'
mols1 = sdf.SDFileParser(file1, maxNumOfMol=10)
print(json.dumps(mols, indent=2))

```

## requirement
tested on python 3.9

## license
MIT
