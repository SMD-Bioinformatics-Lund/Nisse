# v1.1.0 - Sync with Tomte v4.1.2
* Run Tomte pipeline initialization step before Tomte itself.
* Take a single csv argument --input and use for both Tomte and Nisse.
* Update `parse_qc_for_cdm.py` such that it can deal with empty trailing fields in the input file. This is common as single-read and paired-read fields are mixed in the same file on different rows.
* README updates.

# v1.0.3
* Changes to parse_qc_for_cdm.py, rename `results` field to `result`

# v1.0.2
* Changes to parse_qc_for_cdm.py, add "id" field to QC entries

# v1.0.1

* Changes to parse_qc_for_cdm.py, certain labels and values has been altered

# v1.0.0

* First full version, going from csv -> Tomte -> postprocesing -> published results
