#!/usr/bin/python
import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file with latest genomes from GISAID")
    parser.add_argument("--new-genomes", required=True, help="FASTA file with newly sequenced genomes")
    parser.add_argument("--keep", required=True, help="TXT file with accession number of genomes to be included")
    parser.add_argument("--remove", required=True, help="TXT file with accession number of genomes to be removed")
    parser.add_argument("--output", required=True, help="FASTA file containing filtered sequences")
    args = parser.parse_args()

    genomes = args.genomes
    new_genomes = args.new_genomes
    keep = args.keep
    remove = args.remove
    outfile = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/ncov_20210219_testipipeline/'
    # genomes = path + "pre-analyses/gisaid_hcov-19_short.fasta"
    # new_genomes = path + "pre-analyses/new_genomes.fasta"
    # keep = path + 'config/keep.txt'
    # remove = path + "config/remove.txt"
    # outfile = path + "pre-analyses/temp_sequences.fasta"


    # store only new sequences in a dictionary, ignoring existing ones
    print('\n### Loading new sequences\n')
    newly_sequenced = {}
    for fasta in SeqIO.parse(open(new_genomes),'fasta'):
        id, seq = fasta.description, fasta.seq
        if id not in newly_sequenced.keys(): # avoid potential duplicates
            newly_sequenced[id] = str(seq)
    print('\t- Done!\n')


    # create a list of sequences to be added in all instances
    keep_sequences = []
    for id in open(keep, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            if id not in newly_sequenced:
                if id not in keep_sequences:
                    keep_sequences.append(id)


    # create a list of sequences to be ignored in all instances
    remove_sequences = []
    for id in open(remove, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            remove_sequences.append(id)


    # export only sequences to be used in the nextstrain build
    c = 1
    print('\n### Exporting sequences\n')
    exported = []
    all_sequences = []
    with open(outfile, 'w') as output:
        for fasta in SeqIO.parse(open(genomes), 'fasta'):
            id, seq = fasta.description, fasta.seq
            all_sequences.append(id)

            if id not in remove_sequences:
                if id in keep_sequences: # filter out unwanted sequences
                    entry = ">" + id + "\n" + str(seq).upper() + "\n"
                    exported.append(id)
                    output.write(entry)
                    print(str(c) + '. ' + id)
                    c += 1

        for id, seq in newly_sequenced.items():
                print('* ' + str(c) + '. ' + id)
                entry = ">" + id + "\n" + seq.upper() + "\n"
                exported.append(id)
                output.write(entry)
                c += 1
    print('\n- Done!\n')


    # mismatched sequence headers
    mismatch = [genome for genome in keep_sequences if genome not in all_sequences]
    if len(mismatch) > 0:
        print('\n### Possible sequence header mismatches\n')
        m = 1
        for id in mismatch:
            print(str(m) + '. ' + id)
            m += 1
    else:
        print('\nNo sequence name mismatches found...')


    print('\n\n### Final result\n')

    print('Lab file contains ' + str(len(newly_sequenced)) + ' sequences')
    print('GISAID file contains ' + str(len(all_sequences)) + ' sequences\n')

    print(str(len(mismatch)) + ' genomes in keep.txt were NOT FOUND on GISAID database')
    print(str(len(keep_sequences)) + ' genomes ADDED from GISAID dataset')
    print(str(len(newly_sequenced)) + ' newly sequenced genomes were added')
    print(str(len(remove_sequences)) + ' genomes were REMOVED according to remove.txt\n')
    print(str(len(exported)) + ' genomes included in FINAL dataset\n')