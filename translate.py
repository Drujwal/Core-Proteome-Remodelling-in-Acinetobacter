import os

# 1. Set your directories
path = r"D:\Statstics"
dna_folder = os.path.join(path, "Final_Core_With_Outgroup")
prot_folder = os.path.join(path, "174_Aligned_Proteins")

if not os.path.exists(prot_folder):
    os.makedirs(prot_folder)

# Bacterial Genetic Code (Table 11)
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate(dna):
    protein = ""
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3].upper()
        if len(codon) == 3:
            protein += codon_table.get(codon, 'X') # 'X' for unknown/gaps
    return protein

# 2. Translate all 174 files
dna_files = [f for f in os.listdir(dna_folder) if f.endswith(".fasta")]
print(f"Translating {len(dna_files)} files...")

for file in dna_files:
    with open(os.path.join(dna_folder, file), 'r') as f_in:
        with open(os.path.join(prot_folder, file.replace(".fasta", ".faa")), 'w') as f_out:
            lines = f_in.read().split('>')
            for entry in lines:
                if not entry.strip(): continue
                parts = entry.split('\n')
                header = parts[0]
                dna_seq = "".join(parts[1:]).strip().replace("-", "") # Remove gaps for translation
                prot_seq = translate(dna_seq)
                f_out.write(f">{header}\n{prot_seq}\n")

print(f"Done! 174 Protein files created in {prot_folder}")