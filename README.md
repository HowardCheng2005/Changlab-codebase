# Changlab-codebase
Howard's programs for the Chang lab 2026

## Primer Analysis Pipeline
This pipeline uses the IDT Oglioanalyzer API tool to find the best forward and reverse binding sequence based on 
self-dimerization and heterogenous dimerization scores

### Parameters: (requires an API key and secret for your IDT account)
* -f: potential forward binding sequence area
* -r: potential reverse binding sequence area
* -a: forward primer sequence without binding area
* -b: reverse primer sequence without binding area
* -sd: self-dimerization weight for tuning (optional s = 1)
* -h: heterogeneous dimerization weight for tuning (optional h = 1)
* -l: length for the binding sequence itself
* -i: id for API
* -sc: secret for API
* -u: username for IDT
* -p: password for IDT
* -o: output file pathway
