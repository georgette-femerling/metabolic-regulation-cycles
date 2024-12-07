# metabolic-regulation-cycles

Overview  
This algorithm is designed to integrate gene regulation and metabolic information to uncover non-obvious relationships between biological functions. This workflow leverages connections from both metaobolic pathways and regulatory networks to identify cycles that link transcription factors (TFs) with their regulated genes.
Key Features

* Integration of Data Sources: Combines regulatory data from RegulonDB with additional regulatory networks to enhance the analysis of gene functions.
* Cycle Discovery: Identifies cycles starting from an initial transcription factor, retrieving functional categories from MultiFun and the genes regulated by that TF.
* Functional Relationship Mapping: Analyzes genes across different functional categories that catalyze reactions involving shared metabolites.
* Direct and Indirect Connections: Closes cycles by establishing either direct connections (genes regulated by the initial TF) or indirect connections (genes regulated by other TFs influenced by the initial TF).

