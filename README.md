# CATCh V2: Enhanced Chimera Detection Tool

## Overview

**CATCh V2** is an advanced chimera detection tool designed to provide robust and accurate detection of chimeric sequences in microbial community studies. It enhances the original CATCh ensemble classifier by integrating multiple chimera detection tools, including both de novo and reference-based approaches, to deliver improved microbial diversity profiling.

## Background

Microbial communities, including bacteria and archaea, play significant roles in ecosystems, human health, and industrial processes. The detection of chimeric sequences, which are artificial sequences generated during PCR amplification, remains a significant challenge in microbial community profiling. Chimeras can lead to inflated diversity estimates and incorrect identification of microbial taxa. 

CATCh V2 leverages multiple detection tools — Perseus, UCHIME1, VSEARCH, ChimeraSlayer — combining their outputs using machine learning classifiers to provide a more comprehensive and reliable chimera detection system.

## Methodology

### Data Collection and Processing

- **Dataset**: The project utilizes 11 mock microbial communities from the Mockrobiota database. Each mock was processed using the DADA2 pipeline in R to filter, trim, and identify amplicon sequence variants (ASVs).
- **Preprocessing**: ASVs were categorized into chimeric and non-chimeric sequences, with outputs used as input features for various chimera detection tools.

### Chimera Detection Tools

- **Tools Used**: Perseus, UCHIME1, VSEARCH, ChimeraSlayer.
- **Integration Strategy**: Outputs from these tools, including scores and decisions, were combined as input features for machine learning classifiers.

### Building the CATCh V2 Classifier

- **Machine Learning Models**: Multiple classifiers, including SVM, Random Forest, Neural Network, LightGBM, and XGBoost, were tested to select the optimal model.
- **Feature Engineering**: Features extracted from different tools were preprocessed, filtered, and combined for model training and testing.
- **Evaluation**: The classifier's performance was evaluated using accuracy, precision, recall, F1 score, and ROC analysis.

### Model Validation

- **Validation Dataset**: Models were tested on independent validation datasets to verify effectiveness.
- **Performance Metrics**: Sensitivity, specificity, and accuracy were calculated, along with ROC curve analysis for robustness.

## Results and Findings

CATCh V2 improves performance over individual chimera detection tools, especially in environments with limited reference databases. By integrating multiple tools, it offers a more accurate and reliable method for microbial diversity profiling.

## Getting Started

### Prerequisites

- **Operating System**: Linux
- **Programming Languages**: Python, R
- **Libraries and Tools**:
  - Python: scikit-learn, LightGBM, XGBoost, matplotlib, pandas, numpy
  - R: DADA2 package
  - Other tools: Perseus, UCHIME1, VSEARCH, ChimeraSlayer

## Usage Instructions

To use **CATCh V2**, follow these steps:

### Set Up Your Directory

1. **Create a main directory** on your computer.
2. Inside the main directory, **create a subdirectory** for storing results (e.g., `output`).
3. Download the following files from the repository and place them in the main directory:
   - Mock FASTA file (e.g., `mock14.fasta`)
   - Mock count table (e.g., `mock14.count_table`)
   - `merge.py` script
   - `scaler_and_model.joblib` file

4. Additionally, download and install **Mothur** and the **SILVA.bacteria** reference files.

### Run the Chimera Detection Tools

1. Open your terminal and navigate to the main directory.
2. Run the following command:

    ```bash
    python merge.py <input_fasta> <input_count_table> <output_folder> scaler_and_model.joblib
    ```

    **Example:**

    ```bash
    python merge.py mock14.fasta mock14.count_table output scaler_and_model.joblib
    ```

    This command will execute all the chimera detection tools on the input files and save the results to the specified output folder.

### Example Output

The output will be saved in the `output` folder. The chimera detection results for the mock dataset will be stored in a CSV file named `result.csv`. This file will include:

- Predicted chimeric status for each sequence.


#### Sample Output Format

| SequenceName | Chimeric       |
|-------------|-----------------|
| Sequence_1      | Non-chimeric    | 
| Sequence_2      | Non-chimeric    | 

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.


## References

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*. [Link](https://pubmed.ncbi.nlm.nih.gov/27214047/)
- Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R. UCHIME improves sensitivity and speed of chimera detection. *Bioinformatics (Oxford, England)*. [Link](https://pubmed.ncbi.nlm.nih.gov/21700674/)
- Edgar RC. UNOISE2: Improved error-correction for Illumina 16S and its amplicon sequencing. *bioRxiv*. [Link](https://www.biorxiv.org/content/10.1101/081257v1)
- Mysara M, Saeys Y, Leys N, Raes J, Monsieurs P. CATCh: An ensemble classifier for chimera detection in 16S rRNA sequencing studies. *Applied and Environmental Microbiology*. [Link](https://pubmed.ncbi.nlm.nih.gov/25527546/)
- Quince C, Lanzen A, Davenport RJ, Turnbaugh PJ. Removing noise from pyrosequenced amplicons. *BMC Bioinformatics*. [Link](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-38)
- Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: A versatile open-source tool for metagenomics. *PeerJ*. [Link](https://peerj.com/articles/2584/)
- Wright ES, Yilmaz LS, Noguera DR. DECIPHER: A search-based approach to chimera identification for 16S rRNA sequences. *Applied and Environmental Microbiology*. [Link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3264099/)

