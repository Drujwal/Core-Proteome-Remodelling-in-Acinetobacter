import os, subprocess
from Bio import SeqIO
import multiprocessing as mp
from functools import partial

# --- CONFIGURATION ---
BASE_DIR = r"D:\Bioinformatics_Work"
GENOME_DIR = os.path.join(BASE_DIR, "Blast_Temp")
RESULTS_DIR = os.path.join(BASE_DIR, "Core_Genes_Ready_For_Tree")
MAP_DIR = os.path.join(RESULTS_DIR, "Gene_Maps")
DB_DIR = os.path.join(BASE_DIR, "databases")
REF_FILE = os.path.join(GENOME_DIR, "reference.fa")
FINAL_MATRIX = os.path.join(BASE_DIR, "supermatrix.fasta")
BLAST_BIN = r"C:\Program Files\NCBI\blast-2.17.0+\bin"
THREADS = 4

# Filters
OCCUPANCY_THRESHOLD = 40  
IDENTITY_THRESHOLD = 80.0
LENGTH_THRESHOLD = 0.80

# Setup
for d in [RESULTS_DIR, MAP_DIR, DB_DIR]:
    if not os.path.exists(d): os.makedirs(d)

def get_cmd(cmd): return os.path.join(BLAST_BIN, cmd + ".exe")

def process_gene(ref_rec, all_genomes):
    # FIXED: Sanitized the full unique ID for filename safety
    safe_id = ref_rec.id.replace(':', '_').replace('|', '_').replace(' ', '_')
    gene_id = f"Gene_{safe_id}"
    tmp_q = os.path.join(DB_DIR, f"temp_{gene_id}_{os.getpid()}.fa")
    species_found = {}
    
    try:
        with open(tmp_q, "w") as f: f.write(f">query\n{str(ref_rec.seq)}")
        for g_file in all_genomes:
            db_name = os.path.join(DB_DIR, f"{os.path.splitext(g_file)[0]}_db")
            cmd = [get_cmd("blastn"), "-query", tmp_q, "-db", db_name, 
                   "-outfmt", "6 pident length qlen sseq sseqid", "-max_target_seqs", "1"]
            out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode().strip()
            if out:
                pident, aln, qlen, sseq, sseqid = out.split('\t', 4)
                if float(pident) >= IDENTITY_THRESHOLD and (int(aln)/int(qlen)) >= LENGTH_THRESHOLD:
                    species_found[os.path.splitext(g_file)[0]] = (sseq, sseqid)
        
        if len(species_found) >= OCCUPANCY_THRESHOLD:
            with open(os.path.join(RESULTS_DIR, f"{gene_id}.fasta"), "w") as f:
                f.write(f">reference\n{str(ref_rec.seq)}\n")
                for sp, (seq, _) in species_found.items(): f.write(f">{sp}\n{seq}\n")
            return gene_id
        return None
    finally:
        if os.path.exists(tmp_q): os.remove(tmp_q)

def concatenate_genes():
    print("\nConcatenating genes into supermatrix...")
    species_data = {}
    fasta_files = [f for f in os.listdir(RESULTS_DIR) if f.endswith(".fasta")]
    for f in fasta_files:
        for rec in SeqIO.parse(os.path.join(RESULTS_DIR, f), "fasta"):
            if rec.id not in species_data: species_data[rec.id] = ""
            species_data[rec.id] += str(rec.seq)
    with open(FINAL_MATRIX, "w") as out:
        for sp, seq in species_data.items():
            out.write(f">{sp}\n{seq}\n")
    print(f"DONE! Supermatrix saved to: {FINAL_MATRIX}")

if __name__ == "__main__":
    ref_genes = list(SeqIO.parse(REF_FILE, "fasta"))
    genomes = [f for f in os.listdir(GENOME_DIR) if f.lower().endswith(".fa") and f != "reference.fa"]
    
    print("Ensuring BLAST databases exist...")
    for g in genomes:
        db_name = os.path.join(DB_DIR, f"{os.path.splitext(g)[0]}_db")
        if not os.path.exists(db_name + ".nsq"):
            subprocess.run([get_cmd("makeblastdb"), "-in", os.path.join(GENOME_DIR, g), "-dbtype", "nucl", "-out", db_name], stdout=subprocess.DEVNULL)

    print(f"Starting processing {len(ref_genes)} genes...")
    with mp.Pool(processes=THREADS) as pool:
        results = [r for r in pool.map(partial(process_gene, all_genomes=genomes), ref_genes) if r]
    
    print(f"\nFound {len(results)} core genes.")
    if results: concatenate_genes()