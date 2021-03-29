# -*- coding: utf-8 -*-

import pycountry_convert as pyCountry
import pycountry
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file genomes to be used")
    parser.add_argument("--metadata1", required=True, help="Metadata file from NextStrain")
    parser.add_argument("--metadata2", required=False, help="Custom lab metadata file")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Reformatted, final FASTA file")
    args = parser.parse_args()

    genomes = args.genomes
    metadata1 = args.metadata1
    metadata2 = args.metadata2
    output1 = args.output1
    output2 = args.output2

    # genomes = path + 'pre-analyses/temp_sequences.fasta'
    # metadata1 = path + 'pre-analyses/metadata_nextstrain.tsv'
    # metadata2 = path + 'pre-analyses/SC2_Project2.xlsx'
    # output1 = path + 'pre-analyses/metadata_filtered.tsv'
    # output2 = path + 'pre-analyses/sequences.fasta'

    pd.set_option('max_columns', 100)

    # create a dict of existing sequences
    print('\nProcessing sequence file...\n')
    sequences = {}
    for fasta in SeqIO.parse(open(genomes), 'fasta'):
        id, seq = fasta.description, fasta.seq
        if id not in sequences.keys():
            sequences[id] = str(seq)

    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # nextstrain metadata
    dfN = pd.read_csv(metadata1, encoding='utf-8', sep='\t', dtype='str')

    dfN.insert(5, 'iso', '')
    dfN['category'] = ''
    dfN['batch'] = ''
    dfN['sequencing_date'] = ''

    dfN.fillna('', inplace=True)
    list_columns = dfN.columns.values  # list of column in the original metadata file

    # Lab genomes metadata
    dfL = pd.read_excel(metadata2, index_col=None, header=0, sheet_name=0,
                        # 'sheet_name' must be changed to match the Excel sheet name
                        converters={'sample': str, 'collection-date': str, 'category': str, 'batch': str})  # this need to be tailored to your lab's naming system
    dfL = dfL.rename(columns={'collection-date': 'date', 'lab': 'originating_lab'})

    dfL.fillna('', inplace=True)

    # output dataframe
    outputDF = pd.DataFrame(columns=list_columns)
    found = []
    lab_label = {}

    # process metadata from excel sheet
    for idx, row in dfL.iterrows():
        id = dfL.loc[idx, 'id']
        if id in sequences:
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfL.loc[idx, col]  # add values to dictionary

            if dict_row['location'] in ['', None]:
                dict_row['location'] = dfL.loc[idx, 'location']

            collection_date = ''
            if len(str(dict_row['date'])) > 1:
                collection_date = dict_row['date'].split(' ')[0].replace('.', '-').replace('/', '-')
                dict_row['date'] = collection_date

            dfL.loc[idx, 'country'] = 'USA'
            dfL.loc[idx, 'region'] = 'North America'
            dfL.loc[idx, 'country_exposure'] = ''
            dfL.loc[idx, 'division_exposure'] = ''

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                dict_row[level_exposure] = dfL.loc[idx, level_exposure]
                if dict_row[level_exposure] in ['', None]:
                    if level_exposure == 'country_exposure':
                        dict_row[level_exposure] = dict_row[level]
                    else:
                        if dict_row['country_exposure'] != dfL.loc[idx, 'country']:
                            dict_row[level_exposure] = dict_row['country_exposure']
                        else:
                            dict_row[level_exposure] = dict_row[level]

            if row['state'] == '':
                code = 'un'  # change this line to match the acronym of the most likely state of origin if the 'State' field is unknown
            else:
                code = row['state']
            strain = code + '/' + row['sample'] + '/' + id # new strain name

            dict_row['strain'] = strain
            dict_row['iso'] = get_iso(dict_row['country'])
            dict_row['originating_lab'] = dfL.loc[idx, 'originating_lab']
            dict_row['submitting_lab'] = 'Grubaugh Lab - Yale School of Public Health'
            dict_row['authors'] = 'Fauver et al'
            dict_row['category'] = row['category']
            dict_row['batch'] = 'Batch' + str('0' * (3 - len(row['batch']))) + row['batch']
            dict_row['sequencing_date'] = row['sequencing-collection-date'].strftime('%Y-%m-%d')

            # add lineage
            lineage = ''
            if dfL.loc[idx, 'pango_lineage'] != '':
                lineage = dfL.loc[idx, 'pango_lineage']
            dict_row['pango_lineage'] = lineage

            found.append(strain)
            lab_label[id] = strain

            outputDF = outputDF.append(dict_row, ignore_index=True)

    # process metadata from TSV
    dfN = dfN[dfN['strain'].isin(sequences.keys())]
    for idx, row in dfN.iterrows():
        strain = dfN.loc[idx, 'strain']
        if strain in sequences:
            if strain in outputDF['strain'].to_list():
                continue
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfN.loc[idx, col]

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                if dict_row[level_exposure] in ['', None]:
                    dict_row[level_exposure] = dict_row[level]

            dict_row['iso'] = get_iso(dict_row['country'])
            found.append(strain)

            outputDF = outputDF.append(dict_row, ignore_index=True)


    # write new metadata files
    outputDF = outputDF.drop(columns=['region'])
    outputDF.to_csv(output1, sep='\t', index=False)


    # write sequence file
    exported = []
    with open(output2, 'w') as outfile2:
        # export new metadata lines
        for id, sequence in sequences.items():
            if id in lab_label:
                if lab_label[id] not in exported:
                    entry = '>' + lab_label[id] + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    print('* Exporting newly sequenced genome and metadata for ' + id)
                    exported.append(lab_label[id])
            else:
                if id not in exported:
                    entry = '>' + id + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    exported.append(id)



    # # write fasta file
    # exported = []
    # print('\n### Exporting genomes and metadata')
    # print('\t Exporting all selected lab sequences, publicly available genomes and metadata')
    # with open(output2, 'w') as outfile2:
    #     # export new fasta entries
    #     for id, sequence in sequences.items():
    #         if len(id) < 5:
    #             if lab_label[id] not in exported:
    #                 entry = '>' + lab_label[id] + '\n' + sequence + '\n'
    #                 outfile2.write(entry)
    #                 print('\t\t* Newly sequenced genome and metadata: ' + id)
    #                 exported.append(lab_label[id])
    #         elif len(id) > 4:
    #             if id not in exported:
    #                 entry = '>' + id + '\n' + sequence + '\n'
    #                 outfile2.write(entry)
    #                 exported.append(id)

print('\nMetadata file successfully reformatted and exported!\n')
