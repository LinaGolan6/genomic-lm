import pandas as pd
from random import randrange
from Bio import SeqIO
import csv
import gzip
import glob

# Function to read the genome file.  It isn't necessary here, but could be helpful
def get_chromo_genome(genome_file):
    identifiers = []
    lengths = []
    seq=[]
    for record in SeqIO.parse(genome_file, "fasta"):
        identifiers.append(record.id)
        lengths.append(len(record))
        seq.append(record.seq)

    percentile_list = pd.DataFrame(
        {'ID': identifiers,
         'lengths': lengths,
         'seq': seq
         })
    return percentile_list


def check_distance(candidate, known_indexes):
    for number in known_indexes:
        if abs(candidate - number) < 600:
            return False
    return True


def generate_negative_indexes(indexes):
    # select random negative TSS within 0-max(known positive index)
    known_indexes = [int(el[1]) for el in indexes]
    negative_indexes = []
    amount = len(known_indexes)
    max_index = max(known_indexes)
    while len(negative_indexes) < amount:
        candidate = randrange(max_index)
        if check_distance(candidate,known_indexes):
            negative_indexes.append(candidate)
            known_indexes.append(candidate)

    return negative_indexes


def read_bed_files(bed_file):
    # read bed file and extract data
    content = []
    with open(bed_file)as f:
        for line in f:
            content.append(line.strip().split())
    # group them into chr
    chr_groups = {}
    for x in content:
        chr_groups.setdefault(x[0], []).append(x)
    neg_data = {}
    for chromo,indexes in chr_groups.items():
        neg_index = generate_negative_indexes(indexes)
        neg_data[chromo] = neg_index

    return neg_data


# Function to read the gzipped FASTA file and return the sequences
def read_fasta(file_path):
    sequences = {}
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = str(record.seq)
    return sequences


# Function to extract subsequence given an index
def extract_subsequence(sequence, index, flanking=300):
    start = max(0, index - flanking)
    end = min(len(sequence), index + flanking + 1)
    return sequence[start:end]


# Main function
def main(fasta_file, indexes_dict, output_csv_path, organism_name):
    # Read the sequences from the gzipped FASTA file
    sequences = read_fasta(fasta_file)

    # Prepare to write to CSV
    with open(output_csv_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write header with the additional label column
        csvwriter.writerow(['organism_name', 'Chromosome', 'Index', 'Subsequence'])

        # Extract subsequences around each index
        for label, indexes in indexes_dict.items():
            if label in sequences:
                sequence = sequences[label]
                for idx in indexes:
                    subseq = extract_subsequence(sequence, idx)
                    csvwriter.writerow([organism_name, label, idx, subseq])


def merge_csv_files(input_folder, output_csv_path):
    # Get list of all CSV files in the input folder
    csv_files = glob.glob(f"{input_folder}/*.csv")

    # List to hold dataframes
    dataframes = []

    # Read each CSV file and append to the list
    for file in csv_files:
        df = pd.read_csv(file)
        dataframes.append(df)

    # Concatenate all dataframes
    merged_df = pd.concat(dataframes, ignore_index=True)

    # Write the merged dataframe to a new CSV file
    merged_df.to_csv(output_csv_path, index=False)


if __name__ == '__main__':
    ### celegans ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\celegans\celegans\celegans.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\celegans\celegans\ce6.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\celegans\celegans\negative_genome.csv"
    # main(fasta_file, indexes_dict, output_csv, 'celegans')

    ### gallus ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\gallus\gallus\gallus.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\gallus\gallus\galGal5.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\gallus\gallus\negative_genome_gallus.csv"
    # main(fasta_file, indexes_dict, output_csv, 'gallus')

    ### human ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\human\human\human.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\human\human\GRCh38.p13.genome.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\human\human\negative_genome_human.csv"
    # main(fasta_file, indexes_dict, output_csv, 'human')

    ### melanogaster ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\melanogaster\melanogaster\melanogaster.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\melanogaster\melanogaster\dm6.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\melanogaster\melanogaster\negative_genome_melanogaster.csv"
    # main(fasta_file, indexes_dict, output_csv, 'melanogaster')

    ### mulatta ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\mulatta\mulatta\mulatta.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\mulatta\mulatta\rheMac8.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\mulatta\mulatta\negative_genome_mulatta.csv"
    # main(fasta_file, indexes_dict, output_csv, 'mulatta')

    ### musculus ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\musculus\musculus\musculus.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\musculus\musculus\mm10.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\musculus\musculus\negative_genome_musculus.csv"
    # main(fasta_file, indexes_dict, output_csv, 'musculus')

    ### norvegicus ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\norvegicus\norvegicus\norvegicus.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\norvegicus\norvegicus\rn6.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\norvegicus\norvegicus\negative_genome_norvegicus.csv"
    # main(fasta_file, indexes_dict, output_csv, 'norvegicus')

    ### rerio ###
    # neg_data = read_bed_files(bed_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\rerio\rerio\rerio.bed")
    # fasta_file = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\rerio\rerio\danRer7.fa.gz"
    # indexes_dict = neg_data
    # output_csv = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\rerio\rerio\negative_genome_rerio.csv"
    # main(fasta_file, indexes_dict, output_csv, 'rerio')

    # Example usage
    input_folder = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\negative data"
    output_csv_path = r"C:\Users\User\OneDrive\Documents\university\4 year\פרויקט\דאטה גנום\negative data/merged_negative_genome.csv"

    merge_csv_files(input_folder, output_csv_path)