import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# fasta files in the same folder
def extract_features(folder_name):
    # 20 amino acid
    acid_type = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # 400 dipeptide
    dual_peptide_type = []
    for acid1 in acid_type:
        for acid2 in acid_type:
            dual_peptide_type.append(acid1 + acid2)

    # traverse folder_path files
    folder_path = os.listdir(folder_name)
    cnt = 0  # count the number of fasta file
    feature_data = pd.DataFrame()
    for file_name in folder_path:
        single_data = {'uniprot_id': file_name.split('.')[0]}

        # traverse the content of the fasta
        fasta_file_path = folder_name + r'/' + file_name
        sequence = Seq("")  # Sequence
        for sequence_record in SeqIO.parse(fasta_file_path, "fasta"):
            # single record
            sequence = sequence_record.seq

        # filter unknown amid data
        seq_str = str(sequence)
        if 'X' in seq_str:
            single_features = pd.DataFrame(pd.Series(single_data), columns=[file_name])
            feature_data = pd.concat([feature_data, single_features.T])
            continue
        if 'U' in seq_str:
            single_features = pd.DataFrame(pd.Series(single_data), columns=[file_name])
            feature_data = pd.concat([feature_data, single_features.T])
            continue
        if 'Z' in seq_str:
            single_features = pd.DataFrame(pd.Series(single_data), columns=[file_name])
            feature_data = pd.concat([feature_data, single_features.T])
            continue

        cnt += 1  # count the effective data

        # transfer to the format of ProteinAnalysis
        seq_pro = ProteinAnalysis(seq_str)

        acid_frequency = seq_pro.count_amino_acids()
        # print(acid_frequency)
        # single_data['acid frequencies'] = seq_pro.count_amino_acids()
        for key, value in acid_frequency.items():
            single_data[key] = value

        dual_peptide_frequency = {}
        for Dual_Peptide in dual_peptide_type:
            dual_peptide_frequency[Dual_Peptide] = sequence.count(Dual_Peptide)
        for key, value in dual_peptide_frequency.items():
            single_data[key] = value
        # print(Dual_Peptide_frequency)

        # protein_length =  seq_pro.length
        # print("%d" % protein_length)
        single_data['protein length'] = seq_pro.length

        # molecular_weight = seq_pro.molecular_weight()
        # print("%0.2f" % molecular_weight)
        single_data['molecular weight'] = seq_pro.molecular_weight()

        # aromaticity = seq_pro.aromaticity()
        # print("%0.2f" % aromaticity)
        single_data['aromaticity'] = seq_pro.aromaticity()

        # instability_index = seq_pro.instability_index()
        # print("%0.2f" % instability_index)
        single_data['instability index'] = seq_pro.instability_index()

        # isoelectric_point = seq_pro.isoelectric_point()
        # print("%0.2f" % isoelectric_point)
        single_data['isoelectric point'] = seq_pro.isoelectric_point()

        sec_struc = seq_pro.secondary_structure_fraction()
        # print("%0.2f, %0.2f, %0.2f" % (sec_struc[0], sec_struc[1], sec_struc[2]))
        single_data['helix'] = sec_struc[0]
        single_data['turn'] = sec_struc[1]
        single_data['sheet'] = sec_struc[2]

        single_features = pd.DataFrame(pd.Series(single_data), columns=[file_name.split('.')[0]])
        # print(single_features)
        feature_data = pd.concat([feature_data, single_features.T])

        print(cnt)
        # observe data
        # if cnt >= 100:
        #     break

    return feature_data


if __name__ == '__main__':
    # fasta files in the same folder
    fasta_folder_path = r'fasta_name'
    features = extract_features(fasta_folder_path)  # DataFrame
    print(features)
    organism_data = pd.read_excel('merge_excel_1.xlsx')
    print(organism_data)
    final_dataset = pd.merge(organism_data, features, how='outer')
    print(final_dataset)
    final_dataset.to_excel('final_dataset_2022_11_26_21.xlsx')
