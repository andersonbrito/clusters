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

    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_nfl/nextstrain/batch01_20201012e/pre-analyses/"
    # genomes = path + 'temp_sequences.fasta'
    # metadata1 = path + 'metadata_nextstrain.tsv'
    # metadata2 = path + 'SC2_Project2.xlsx'
    # output1 = path + 'metadata_filtered.tsv'
    # output2 = path + 'sequences.fasta'

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
    try:
        dfN = dfN[['strain', 'gisaid_epi_isl', 'genbank_accession', 'date', 'country', 'division', 'location',
                      'region_exposure', 'country_exposure', 'division_exposure', 'originating_lab', 'submitting_lab']]
        dfN.insert(4, 'iso', '')
    except:
        pass
    dfN['category'] = ''
    dfN['batch'] = ''
    dfN['sequencing_date'] = ''

    dfN.fillna('', inplace=True)
    lColumns = dfN.columns.values  # list of column in the original metadata file

    # Lab genomes metadata
    dfL = pd.read_excel(metadata2, index_col=None, header=0, sheet_name=0,
                        # 'sheet_name' must be changed to match the Excel sheet name
                        converters={'sample': str, 'collection-date': str, 'category': str, 'batch': str})  # this need to be tailored to your lab's naming system
    dfL.fillna('', inplace=True)
    dfL.set_index('id', inplace=True)

    dHeaders = {}
    notFound = []
    lstNewMetadata = []
    found = []
    lab_label = {}
    for id in sequences.keys():
        # check nextstrain metadata first
        dRow = {}
        if id in dfN['strain'].to_list() and len(id) > 4:
            fields = {column: '' for column in lColumns}
            row = dfN.loc[lambda dfN: dfN['strain'] == id]

            strain = row.strain.values[0]
            country = row.country.values[0]
            division = row.division.values[0]
            location = row.location.values[0]

            country_exposure = row.country_exposure.values[0].strip()
            # print('\t > ', country_exposure, country)

            # fix exposure
            if country_exposure in ['', None]:  # ignore travel cases
                # print('\t - ', country_exposure, country)
                country_exposure = country

            division_exposure = row.division_exposure.values[0].strip()
            # print('\t - ', division, division_exposure)

            if division_exposure in ['', None]:  # ignore travel cases
                # print('\t - ', division, division_exposure)
                division_exposure = division

            if len(country) < 2:
                row.country.values[0] = ''
            if len(division) < 2:
                row.division.values[0] = country

            iso = get_iso(country_exposure)
            row.iso.values[0] = iso  # needed for exporting a renaming file
            date = row.date.values[0]
            header = '|'.join([strain, iso, division.replace(' ', '-'), date])
            dHeaders[strain] = header

            lValues = row.values[0]
            for field, value in zip(fields.keys(), lValues):
                # print(id)
                if value in ['', np.nan, None]:
                    value = ''

                fields[field] = value
            if country == '':
                continue

            dRow[id] = fields
            found.append(strain)
            # print('Exporting metadata for ' + id)

        # check lab's metadata otherwise
        if id not in dRow.keys():
            # check lab metadata
            if len(id) < 5:
                id = str(id)
                lab_label[id] = ''

                if id in dfL.index:
                    fields = {column: '' for column in lColumns}
                    row = dfL.loc[id]
                    if row['state'] == '':
                        code = 'XX'  # change this line to match the acronym of the most likely state of origin if the 'State' field is unknown
                    else:
                        code = row['state']
                    strain = code + '/' + row['sample'] + '/' + id # new strain name

                    if strain not in found:
                        gisaid_epi_isl = ''
                        genbank_accession = ''
                        if len(str(row['collection-date'])) > 1:
                            date = row['collection-date'].split(' ')[0].replace('.', '-').replace('/', '-')
                        else:
                            date = ''

                        country = row['country']
                        division = row['division']
                        if row['location'] in ['', '?', 'N/A']:
                            location = ''
                        else:
                            location = str(row['location'])

                        region_exposure = 'North America'
                        country_exposure = 'USA'
                        iso = get_iso(country_exposure)
                        division_exposure = division
                        originating_lab = row['lab']
                        submitting_lab = 'Grubaugh Lab - Yale School of Public Health'  # change this line to match you lab's name
                        category = row['category']
                        batch = 'Batch' + str('0' * (3 - len(row['batch']))) + row['batch']
                        sequencing_date = row['sequencing-collection-date'].strftime('%Y-%m-%d')

                        lValues = [strain, gisaid_epi_isl, genbank_accession, date, iso, country, division, location,
                                   region_exposure, country_exposure, division_exposure, originating_lab, submitting_lab,
                                   category, batch, sequencing_date]

                        header = '|'.join([strain, country, division.replace(' ', '-'), date])
                        dHeaders[strain] = header

                        for field, value in zip(fields.keys(), lValues):
                            fields[field] = value
                        dRow[id] = fields
                        found.append(strain)

                        if id in lab_label.keys():
                            lab_label[id] = strain
                    else:
                        continue

            else:  # Assign 'NA' if no metadata is available
                header = '|'.join([id, 'NA', 'NA', 'NA', 'NA'])
                dHeaders[id] = header
                notFound.append(id)
        lstNewMetadata = lstNewMetadata + list(dRow.values())



    # write new metadata files
    outputDF = pd.DataFrame(lstNewMetadata, columns=list(lColumns))
    outputDF.to_csv(output1, sep='\t', index=False)

    if len(notFound) > 0:
        print('\nPlease check for inconsistencies (see above).')

    # write fasta file
    exported = []
    print('\n### Exporting genomes and metadata')
    print('\t Exporting all selected lab sequences, publicly available genomes and metadata')
    with open(output2, 'w') as outfile2:
        # export new fasta entries
        for id, sequence in sequences.items():
            if len(id) < 5:
                if lab_label[id] not in exported:
                    entry = '>' + lab_label[id] + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    print('\t\t* Newly sequenced genome and metadata: ' + id)
                    exported.append(lab_label[id])
            elif len(id) > 4:
                if id not in exported:
                    entry = '>' + id + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    exported.append(id)

print('\nMetadata file successfully reformatted and exported!\n')
