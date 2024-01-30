# GRIFO
GRIFO is a versatile pipeline for environmental metabarcoding using nanopore sequences with density based clustering for error correction and other utilities.
The name is derived of the mythical creature <em>grifo</em> (Esperanto) which is partly eagle, lion, and deer. The name describes how the pipeline integrates existing tools into a new workflow.

## Usage

GRIFO is under development. There is currently no user guide and source code is unstable and the latest version may not even be functional at times.
Currently, GRIFO runs in a jupyter notebook in Python using wrapper functions for Linux tools.

## Dependencies

GRIFO uses the clustering module of <em>ashure</em>, which itself relies on the following tools:

``` bash
pip install pandas          # for organizing underlying data
pip install scikit-learn    # for clustering
pip install hdbscan         # for clustering
pip install spoa            # for clustering
pip install parasail
```

Other depencies are: qcat, nanofilt, vsearch, minimap2, and crest4

# Citation

GRIFO makes use of several tools and datasets from third parties. Please cite them too when using GRIFO.

- ASHURE - https://github.com/BBaloglu/ASHURE
- BOLD - https://www.boldsystems.org/index.php/datarelease
- BLAST - https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- CREST - https://github.com/lanzen/CREST
- GUPPY - https://nanoporetech.com/
- MINIMAP2 - https://github.com/lh3/minimap2
- NANOFILT - https://github.com/wdecoster/nanofilt/
- PR2 database - https://pr2-database.org/
- QCAT - https://github.com/nanoporetech/qcat
- RDPClassifier - https://github.com/rdpstaff/classifier
- SILVA database - https://www.arb-silva.de/no_cache/download/archive/current/Exports/
- eDNA dataset demonstrator - https://github.com/iobis/dataset-edna
- VSEARCH-2.9.1 - https://github.com/torognes/vsearch/releases/tag/v2.9.1
- BioPython
