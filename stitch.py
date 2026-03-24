import os
from Bio import SeqIO

GENE_FOLDER = r"D:\Bioinformatics_Work\Core_Genes_Ready_For_Tree"
OUTPUT_FILE = r"D:\Bioinformatics_Work\supermatrix.fasta"

def create_supermatrix():
    gene_files = [f for f in os.listdir(GENE_FOLDER) if f.endswith(".fasta")]
    
    # 1. First pass: Find total target length
    # We need to know what the total length should be by adding up max lengths of genes
    species_seqs = {}
    
    # Let's track the total length we are aiming for
    total_expected_length = 0
    
    # We need to know the length of each gene segment
    gene_lengths = []
    
    print("Calculating total matrix dimensions...")
    for f in gene_files:
        lengths = [len(rec.seq) for rec in SeqIO.parse(os.path.join(GENE_FOLDER, f), "fasta")]
        max_len = max(lengths)
        gene_lengths.append(max_len)
        total_expected_length += max_len

    print(f"Targeting total alignment length: {total_expected_length} bp")

    # 2. Second pass: Build the matrix with explicit padding
    all_species = set()
    for f in gene_files:
        for record in SeqIO.parse(os.path.join(GENE_FOLDER, f), "fasta"):
            all_species.add(record.id)
            
    final_matrix = {sp: "" for sp in all_species}
    
    for i, f in enumerate(gene_files):
        gene_len = gene_lengths[i]
        current_data = {rec.id: str(rec.seq) for rec in SeqIO.parse(os.path.join(GENE_FOLDER, f), "fasta")}
        
        for sp in all_species:
            if sp in current_data:
                seq = current_data[sp]
                # Pad/Truncate specifically to gene_len
                if len(seq) < gene_len:
                    seq = seq.ljust(gene_len, '-')
                elif len(seq) > gene_len:
                    seq = seq[:gene_len]
                final_matrix[sp] += seq
            else:
                final_matrix[sp] += "-" * gene_len

    # 3. Save
    with open(OUTPUT_FILE, "w") as f_out:
        for sp, full_seq in final_matrix.items():
            f_out.write(f">{sp}\n{full_seq}\n")
            
    print("SUCCESS: supermatrix.fasta created with uniform length.")

if __name__ == "__main__":
    create_supermatrix()