# CATCh V2: Enhanced Chimera Detection Tool

## Overview

**CATCh V2** is an advanced chimera detection tool designed to provide robust and accurate detection of chimeric sequences in microbial community studies. It enhances the original CATCh ensemble classifier by integrating multiple chimera detection tools, including both de novo and reference-based approaches, to deliver improved microbial diversity profiling.

## Background

Microbial communities, including bacteria and archaea, play significant roles in ecosystems, human health, and industrial processes. The detection of chimeric sequences, which are artificial sequences generated during PCR amplification, remains a significant challenge in microbial community profiling. Chimeras can lead to inflated diversity estimates and incorrect identification of microbial taxa. 

CATCh V2 leverages multiple detection tools — Perseus, UCHIME, Vsearch, ChimeraSlayer, UCHIME3, and DADA2 — combining their outputs using machine learning classifiers to provide a more comprehensive and reliable chimera detection system.

## Methodology

### Data Collection and Processing

- **Dataset**: The project utilizes 11 mock microbial communities from the Mockrobiota database. Each mock was processed using the DADA2 pipeline in R to filter, trim, and identify amplicon sequence variants (ASVs).
- **Preprocessing**: ASVs were categorized into chimeric and non-chimeric sequences, with outputs used as input features for various chimera detection tools.

### Chimera Detection Tools

- **Tools Used**: Perseus, UCHIME, Vsearch, ChimeraSlayer, UCHIME3, and DADA2.
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
  - Other tools: Perseus, UCHIME, Vsearch, ChimeraSlayer

### Usage Instructions

#### Step 1: Data Input
* Prepare your microbial datasets according to the instructions in the methodology section.
* Ensure that your data is correctly formatted for processing.

#### Step 2: Run Chimera Detection Tools
* Execute the following command in your terminal to run the chimera detection tools:
* python merge.py <input_fasta> <input_count_table> <output_folder> scaler_and_model.joblib


#### The `merge.py` script will:
 * 1. Run different chimera detection tools on your specified mock dataset.
 * 2. Combine the outputs from these tools into a unified format.
 * 3. Feed the combined data into CATCh V2 for classification as chimeric or non-chimeric.

#### Step 3: View Results
* The output from the `merge.py` script will indicate the classification results.
* Each sequence in the mock dataset will be labeled as either chimeric or non-chimeric.
