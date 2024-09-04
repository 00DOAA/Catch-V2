import argparse
import os
import csv
import pandas as pd
import numpy as np
from mothur_py import Mothur
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, f1_score, precision_score
import joblib
import xgboost as xgb
import shutil

def run_mothur_commands(fasta_file, count_file, output_dir):
    m = Mothur(mothur_path='mothur/mothur.exe')

   
    before_files = set(os.listdir())

    # Run Mothur commands
    m.chimera.vsearch(fasta=fasta_file, count=count_file)
    m.chimera.perseus(fasta=fasta_file, count=count_file)
    m.chimera.uchime(fasta=fasta_file, count=count_file)

    alignment_output = fasta_file.replace('.fasta', '.align')
    m.align.seqs(fasta=fasta_file, reference='silva.bacteria.fasta')

    m.chimera.slayer(fasta=alignment_output, count=count_file, reference='self')


    after_files = set(os.listdir())

  
    new_files = after_files - before_files

    # Move all new files to the output directory
    for file in new_files:
        if os.path.isfile(file):  
            shutil.move(file, os.path.join(output_dir, file))

    print("Mothur commands executed successfully and files moved.")




def convert_chimera_to_csv(input_file_path, output_file_path):
    with open(input_file_path, 'r') as chimera_file:
        chimera_data = chimera_file.read()
 

    data_rows = [line.split('\t') for line in chimera_data.strip().split('\n')]

    if 'vsearch.chimeras' in os.path.basename(input_file_path):
        header = ['vs1', 'vs2', 'vs3', 'vs4', 'vs5', 'vs6', 'vs7', 'vs8', 'vs9', 'vs10', 'vs11', 'vs12', 'vs13', 'vs14', 'vs15', 'vs16', 'vs17', 'vs18']
        data_rows.insert(0, header)

    with open(output_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data_rows)

    print(f"Data successfully converted to {output_file_path}")

def sort_csv_file(file_path, column_name):
    df = pd.read_csv(file_path)
    df = df.sort_values(by=column_name)
    df.to_csv(file_path, index=False)
    print(f"CSV file {file_path} sorted by {column_name} and saved.")

def merge_csv_files(csv_files, merged_file_path):
    dataframes = [pd.read_csv(file) for file in csv_files]
    merged_df = pd.concat(dataframes, axis=1)
    
    column_order = [
        'Score', 'Query', 'ParentA', 'ParentB', 'IdQM', 'IdQA', 'IdQB', 'IdAB', 'IdQT',
        'LY', 'LN', 'LA', 'RY', 'RN', 'RA', 'Div', 'YN', 'SequenceName', 'Name',
        'DiffsToBestMatch', 'BestMatchIndex', 'BestMatchName', 'DiffstToChimera',
        'IndexofLeftParent', 'IndexOfRightParent', 'NameOfLeftParent', 'NameOfRightParent',
        'DistanceToBestMatch', 'cIndex', '(cIndex - singleDist)', 'loonIndex',
        'MismatchesToChimera', 'MismatchToTrimera', 'ChimeraBreakPoint', 'LogisticProbability',
        'TypeOfSequence', 'LeftParent', 'RightParent', 'DivQLAQRB', 'PerIDQLAQRB', 'BootStrapA',
        'DivQLBQRA', 'PerIDQLBQRA', 'BootStrapB', 'Flag', 'LeftWindow', 'RightWindow', 'vs1',
        'vs2', 'vs3', 'vs4', 'vs5', 'vs6', 'vs7', 'vs8', 'vs9', 'vs10', 'vs11', 'vs12', 'vs13',
        'vs14', 'vs15', 'vs16', 'vs17', 'vs18'
    ]
    
    merged_df = merged_df[[col for col in column_order if col in merged_df.columns]]
    merged_df.to_csv(merged_file_path, index=False)
    print(f"CSV files merged into {merged_file_path}")

def process_csv_file(file_path):
    df = pd.read_csv(file_path)
    df['abundance'] = df['Query'].str.extract(r'\b(?:ab=)(\d+)', expand=False)
    df['Query'] = df['Query'].str.extract(r'^(.*?)\/', expand=False)
    df.to_csv(file_path, index=False)
    print(f"CSV file {file_path} processed and saved.")

def rename(file_path):
    df = pd.read_csv(file_path)
    df = df.rename(columns={'Name': 'SequenceName'})
    df.to_csv(file_path, index=False)
    print(f"CSV file {file_path} column 'Name' renamed to 'SequenceName' and saved.")

def preprocess_dataframe(df):
    df.replace('*', np.nan, inplace=True)
    df['abundance'] = df['Query'].str.extract(r'\b(?:ab=)(\d+)', expand=False)
    
    df_name = df[['SequenceName']]
    
    columns_to_remove = [
        'Name', 'vs2', 'vs3', 'vs4', 'vs5', 'LeftParent', 'RightParent',
        'NameOfLeftParent', 'NameOfRightParent', 'ParentA', 'ParentB',
        'BestMatchName', 'Query', 'SequenceName', 'LeftWindow', 'RightWindow',
        'SequenceIndex','Name.1','YN_'
    ]
    df.drop(columns=columns_to_remove, inplace=True, errors='ignore')
    columns_to_convert = ['IdQM', 'IdQA', 'IdQB', 'IdAB', 'IdQT', 'Div', 'vs6', 'vs7', 'vs8', 'vs9', 'vs10', 'vs17']
    for col in columns_to_convert:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    columns_to_convert_to_int = ['LY', 'LN', 'LA', 'RY', 'RN', 'RA', 'abundance', 'vs11', 'vs12', 'vs13', 'vs14', 'vs15', 'vs16']
    for col in columns_to_convert_to_int:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
    most_frequent_value = df['vs18'].value_counts().idxmax()
    df['vs18'].replace('?', most_frequent_value, inplace=True)
    df['Flag'].fillna('no', inplace=True)
    df['YN'] = df['YN'].fillna('N')
    df['TypeOfSequence'] = df['TypeOfSequence'].apply(lambda x: 0 if x == 'good' else 1)
    df['YN'] = df['YN'].apply(lambda x: 0 if x == 'N' else 1)
    df['Flag'] = df['Flag'].apply(lambda x: 0 if x == 'no' else 1)
    df['vs18'] = df['vs18'].apply(lambda x: 0 if x == 'N' else 1)
    df = df.fillna(0)
    
    return df, df_name

def load_model_and_scaler(model_path):
    scaler, model = joblib.load(model_path)
    return scaler, model

def calculate_metrics(yy_true, yy_pred):
    conf_matrix = confusion_matrix(yy_true, yy_pred)
    sensitivity = conf_matrix[1, 1] / (conf_matrix[1, 0] + conf_matrix[1, 1]) if (conf_matrix[1, 0] + conf_matrix[1, 1]) > 0 else 0
    specificity = conf_matrix[0, 0] / (conf_matrix[0, 0] + conf_matrix[0, 1]) if (conf_matrix[0, 0] + conf_matrix[0, 1]) > 0 else 0
    accuracy = (conf_matrix[0, 0] + conf_matrix[1, 1]) / conf_matrix.sum()
    f1 = f1_score(yy_true, yy_pred)
    precision = precision_score(yy_true, yy_pred)
    return sensitivity, specificity, accuracy, f1, precision

def main():
    parser = argparse.ArgumentParser(description='Process .fasta and .count_table files using Mothur, convert and merge resulting .chimera files, and classify sequences.')
    parser.add_argument('fasta', metavar='fasta_file', type=str, help='The input .fasta file')
    parser.add_argument('count', metavar='count_file', type=str, help='The input .count_table file')
    parser.add_argument('output_dir', metavar='output_dir', type=str,nargs='?',default=os.getcwd(), help='The output directory where the processed files will be saved. This folder should already exist.')
    parser.add_argument('model', metavar='model_file', type=str, help='The path to the trained model file')
    parser.add_argument('--threshold', type=float, default=0.488, help='Probability threshold for classification')


    args = parser.parse_args()
    

    run_mothur_commands(args.fasta, args.count, args.output_dir)

    chimera_files = [os.path.join(args.output_dir, f) for f in os.listdir(args.output_dir) if f.endswith('.chimeras')]
    csv_files = []
    for chimera_file in chimera_files:
        csv_file = chimera_file.replace('.chimeras', '.csv')
        convert_chimera_to_csv(chimera_file, csv_file)
        csv_files.append(csv_file)

    merged_csv = os.path.join(args.output_dir, 'merged_output.csv')
    merge_csv_files(csv_files, merged_csv)

    process_csv_file(merged_csv)
    # sort_csv_file(merged_csv, 'Score')
    rename(merged_csv)

    df = pd.read_csv(merged_csv)
    df, df_name = preprocess_dataframe(df)

    scaler, model = load_model_and_scaler(args.model)

    drop_columns = [
        'numParents', 'mock_id', 'DiffsToBestMatch', 'BestMatchIndex', 'DistanceToBestMatch',
        '(cIndex - singleDist)', 'loonIndex', 'ChimeraBreakPoint', 'TypeOfSequence', 'ischimera', 'uchime3', 'BootStrapA', 'DivQLBQRA', 'BootStrapB', 'Flag'
    ]
    XX = df.drop(columns=drop_columns, errors='ignore')

    for col in XX.columns:
        XX[col] = pd.to_numeric(XX[col], errors='coerce')

    XX_scaled = scaler.transform(XX)

    probabilities = model.predict_proba(XX_scaled)[:, 1]
    predictions = (probabilities >= args.threshold).astype(int)

    df['Tool_result'] = predictions
    df['Chimeric'] = df['Tool_result'].apply(lambda x: 'Chimeric' if x == 1 else 'Non-chimeric')

    output_df = df_name.copy()
    output_df['Chimeric'] = df['Chimeric']
    result_file_path = os.path.join(args.output_dir, 'result.csv')
    output_df.to_csv(result_file_path, index=False)
    print(f"Results saved to {result_file_path}")



if __name__ == '__main__':
    main()
