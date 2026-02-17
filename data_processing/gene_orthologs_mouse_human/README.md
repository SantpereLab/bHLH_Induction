
### GSM to SRX Accession Mapper

This script automates the retrieval of SRA Experiment (SRX) identifiers from GEO Sample (GSM) codes to facilitate raw data downloading.

* **Role**: Links GEO metadata to the Sequence Read Archive (SRA).
* **Function**:
  * Parses `chipseq_datasets.csv` to extract and clean a list of unique GSM identifiers.
  * Uses the `GEOquery` library to fetch online metadata for each sample.
  * Extracts the specific **SRX** accession from the GEO "relation" field.


* **Output**: Generates a mapping file (`srx_accessions_df.txt`) used as input for SRA-Toolkit downloading scripts.
