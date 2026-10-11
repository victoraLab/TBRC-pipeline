import argparse
import subprocess
import os
import matplotlib.pyplot as plt

def get_fasta_ids(filename):
    with open(filename, 'r') as f:
        ids = [line.strip()[1:] for line in f if line.startswith('>')]
    return ids

def filter_and_get_lengths_from_fastq(fastq_file, ids):
    with open("temp_ids.txt", "w") as f:
        for id in ids:
            f.write(id + "\n")

    # Filter the FASTQ file using seqtk
    filtered_fastq = "filtered_" + os.path.basename(fastq_file)
    cmd = ["seqtk", "subseq", fastq_file, "temp_ids.txt"]
    with open(filtered_fastq, "w") as out_file:
        result = subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise ValueError(f"Error filtering FASTQ with seqtk: {result.stderr}")

    # Now, let's extract the lengths from the filtered FASTQ file
    lengths = []
    with open(filtered_fastq, "r") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # This is a sequence line
                lengths.append(len(line.strip()))

    # Cleanup
    #os.remove("temp_ids.txt")
    #os.remove(filtered_fastq)

    return lengths


def calculate_subtracted_lengths(fasta_lengths, r1_lengths, r2_lengths):
    if len(fasta_lengths) != len(r1_lengths) or len(fasta_lengths) != len(r2_lengths):
        raise ValueError("Mismatch in number of sequences among the provided files.")

    return [fasta - r1 - r2 for fasta, r1, r2 in zip(fasta_lengths, r1_lengths, r2_lengths)]

def plot_histogram(data, output_prefix, title):
    plt.hist(data, bins=50, color='blue', alpha=0.7)
    plt.title(title)
    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.savefig(output_prefix + '_' + title.replace(' ', '_').lower() + '.png')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate histograms from FASTA and FASTQ data.")
    parser.add_argument('-f', '--fasta', required=True, help="Input FASTA file.")
    parser.add_argument('--r1', required=True, help="Input FASTQ R1 file.")
    parser.add_argument('--r2', required=True, help="Input FASTQ R2 file.")
    parser.add_argument('-o', '--output', required=True, help="Output prefix for the histogram plots.")

    args = parser.parse_args()

    ids = get_fasta_ids(args.fasta)

    fasta_lengths = filter_and_get_lengths_from_fastq(args.fasta, ids)
    r1_lengths = filter_and_get_lengths_from_fastq(args.r1, ids)
    r2_lengths = filter_and_get_lengths_from_fastq(args.r2, ids)

    subtracted_lengths = calculate_subtracted_lengths(fasta_lengths, r1_lengths, r2_lengths)

    plot_histogram(fasta_lengths, args.output, "FASTA Lengths")
    plot_histogram(subtracted_lengths, args.output, "Subtracted Lengths")

if __name__ == "__main__":
    main()

