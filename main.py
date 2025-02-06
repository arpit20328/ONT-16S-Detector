from collections import Counter
import gzip

def classify_reads(fastq_file, length_threshold=350, uniformity_percent=90):
    read_lengths = []
    
    # Open gzipped or normal FASTQ file
    open_func = gzip.open if fastq_file.endswith(".gz") else open
    with open_func(fastq_file, 'rt') as fq:
        for i, line in enumerate(fq):
            if i % 4 == 1:  # FASTQ sequence line
                read_lengths.append(len(line.strip()))

    if not read_lengths:
        return "Error: No reads found."

    # Count occurrences of each read length
    length_counts = Counter(read_lengths)
    total_reads = len(read_lengths)

    # Get the most common read length and its percentage
    most_common_length, most_common_count = length_counts.most_common(1)[0]
    most_common_percentage = (most_common_count / total_reads) * 100

    # Classification Logic
    if most_common_length <= length_threshold and most_common_percentage >= uniformity_percent:
        return f"Short-read data (Illumina) - {most_common_percentage:.2f}% reads are {most_common_length} bp"
    elif max(read_lengths) > 1000:
        return f"Long-read data (Nanopore/PacBio) - Max read length: {max(read_lengths)} bp"
    else:
        return "Uncertain classification, mixed-length data."

# Example usage
fastq_file = "your_reads.fastq.gz"  # Change this to your file path
print(classify_reads(fastq_file))
