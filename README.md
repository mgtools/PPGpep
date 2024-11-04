
- **Reduced misidentification** of human peptides as microbial peptides through a more complete proteome representation.

## Acknowledgements# PPGpep: A Panproteome Graph for Enhanced Peptide and Protein Identification

This repository contains the code and resources for generating the **Panproteome Graph (PPG)** and performing statistical analysis, as described in our paper. The PPG enables improved peptide and protein identification from proteomic data by utilizing the **human pangenome**, providing deeper insights into the shared peptides across the proteome.

## Data Availability
The datasets supporting this work are accessible at [PPGnet](https://omics.informatics.indiana.edu/PPGnet/).

## Overview
In this study, we developed a novel data structure called the **Panproteome Graph (PPG)**, where **nodes represent tryptic peptides**. The PPG was constructed using **47 human proteomes** from the Human Pangenome Reference Consortium (HPRC) and **UniProt human proteins**, resulting in over **4.2 million tryptic peptides**—a **26% increase** compared to using UniProt proteins alone.

The PPG was applied to multiple **human proteomic and metaproteomic datasets**, leading to an **8% increase in identified peptides** across all collections. This approach also addresses the issue of **misidentification of human peptides as microbial peptides**, which had previously been explored primarily through genomic data.


## Repository Structure
```
PPGpep/
│
├── analysis/               # Statistical analysis scripts and outputs
│   ├── analyze/            # Generated plots of distribution
│   ├── filter_graph.py     # Generate filtered disconnected components for analysis
│   └── path_match/         # Match components to proteins
│
├── PanProtGraph/           # Python programs for generating and managing the PPG
│   ├── PanProteomeGraph.py      # Generate PPG from all protein sequences
│   ├── Graph2Pep.py             # Generate cleaved peptides from the PPG
│   ├── UpdatePPG.py             # Update PPG when new genome proteins are available
│   ├── sortpeptides.py          # Sort peptides in lexicographic order
│   ├── EthnicPeptides.py        # Generate tryptic peptides for each ethnic group
│   ├── EthnicPeptideStats.py    # Count tryptic peptides by ethnic group
│   └── CountPTMs.py             # Count PTM sites by ethnic group
│
└── README.md               # Project documentation (this file)
```

## Installation
Clone the repository and install the required dependencies:

```bash
git clone https://github.com/your-username/PPGpep.git
cd PPGpep
pip install networkx matplotlib pandas
```

## Usage

- Use **`Graph2Pep.py`** to generate cleaved peptides from the PPG.
- Run **`UpdatePPG.py`** to incorporate new protein data into the existing PPG.
- **`EthnicPeptides.py`** and **`EthnicPeptideStats.py`** allow generation and analysis of peptides by ethnic groups.

## Results
- **26% more tryptic peptides** identified with the PPG compared to using UniProt alone.
- **8% increase in peptide identification** across human proteomic and metaproteomic datasets.
This project is based on data from the **Human Pangenome Reference Consortium (HPRC)** and **UniProt**.
