from pathlib import Path
import re
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# ============================================================
# USER SETTINGS
# ============================================================

FOLDER = Path(
    r"C:\Users\xenou\OneDrive\Desktop\Major REVISION AcinetobActer\CALCULATE"
)

# Bacterial genetic code
TRANSLATION_TABLE = 11

# Same protein-length filtering used in the previous analysis
MIN_PROTEIN_LENGTH = 30

# Number of bootstrap replicates for genome-level mean pI
BOOTSTRAP_REPS = 2000

# Reproducible results
RANDOM_SEED = 20260918


# ============================================================
# STANDARD AMINO ACIDS
# ============================================================

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ============================================================
# FIND FASTA FILES
# ============================================================

FASTA_EXTENSIONS = {
    ".fa",
    ".fasta",
    ".fna",
    ".ffn",
    ".fas",
    ".cds"
}

fasta_files = sorted(
    [
        f for f in FOLDER.iterdir()
        if f.is_file()
        and f.suffix.lower() in FASTA_EXTENSIONS
    ]
)


if len(fasta_files) == 0:
    raise FileNotFoundError(
        "\nNo FASTA files were found in:\n"
        + str(FOLDER)
    )


print("=" * 80)
print("ACINETOBACTER PROTEOME RECONSTRUCTION")
print("=" * 80)
print(f"Folder: {FOLDER}")
print(f"FASTA files found: {len(fasta_files)}")
print()


# ============================================================
# FUNCTION: GET NAME FROM FILE
# ============================================================

def get_genome_name(filename):

    name = Path(filename).stem

    # Replace underscores with spaces
    name = name.replace("_", " ")

    # Remove repeated spaces
    name = re.sub(r"\s+", " ", name).strip()

    return name


# ============================================================
# FUNCTION: TRANSLATE CDS
# ============================================================

def translate_cds(sequence):

    sequence = str(sequence).upper()

    # Remove gaps/spaces
    sequence = sequence.replace("-", "")
    sequence = sequence.replace(" ", "")

    if len(sequence) < 3:
        return None, "too_short_nt"

    # Remove incomplete terminal codon
    remainder = len(sequence) % 3

    if remainder != 0:
        sequence = sequence[:-remainder]

    if len(sequence) < 3:
        return None, "too_short_nt"

    try:

        protein = str(
            Seq(sequence).translate(
                table=TRANSLATION_TABLE,
                to_stop=False,
                cds=False
            )
        )

    except Exception:

        return None, "translation_error"

    # Remove terminal stop codon
    protein = protein.rstrip("*")

    # Reject internal stop codons
    if "*" in protein:
        return None, "internal_stop"

    # Reject proteins <=30 aa
    if len(protein) <= MIN_PROTEIN_LENGTH:
        return None, "short_protein"

    # Reject ambiguous/non-standard amino acids
    if not set(protein).issubset(STANDARD_AA):
        return None, "nonstandard_aa"

    return protein, "retained"


# ============================================================
# FUNCTION: BOOTSTRAP CI FOR MEAN pI
# ============================================================

def bootstrap_mean_ci(values, repetitions=2000, seed=1):

    values = np.asarray(values, dtype=float)

    if len(values) < 2:
        return np.nan, np.nan

    rng = np.random.default_rng(seed)

    bootstrap_means = np.empty(repetitions)

    for i in range(repetitions):

        sample = rng.choice(
            values,
            size=len(values),
            replace=True
        )

        bootstrap_means[i] = np.mean(sample)

    lower = np.percentile(
        bootstrap_means,
        2.5
    )

    upper = np.percentile(
        bootstrap_means,
        97.5
    )

    return lower, upper


# ============================================================
# STORAGE
# ============================================================

genome_results = []

protein_results = []

qc_results = []


# ============================================================
# PROCESS EACH GENOME
# ============================================================

for number, fasta_file in enumerate(fasta_files, start=1):

    genome_name = get_genome_name(fasta_file.name)

    print(
        f"[{number}/{len(fasta_files)}] "
        f"{fasta_file.name}"
    )

    pI_values = []
    mw_values = []

    total_cds = 0
    retained_proteins = 0

    qc = {
        "too_short_nt": 0,
        "translation_error": 0,
        "internal_stop": 0,
        "short_protein": 0,
        "nonstandard_aa": 0
    }

    # --------------------------------------------------------
    # Read CDS FASTA
    # --------------------------------------------------------

    for record in SeqIO.parse(
        fasta_file,
        "fasta"
    ):

        total_cds += 1

        protein, status = translate_cds(
            record.seq
        )

        # Protein rejected
        if protein is None:

            qc[status] += 1

            continue

        # ----------------------------------------------------
        # Calculate pI and molecular weight
        # ----------------------------------------------------

        try:

            analysis = ProteinAnalysis(
                protein
            )

            protein_pI = analysis.isoelectric_point()

            protein_mw = analysis.molecular_weight()

        except Exception:

            qc["nonstandard_aa"] += 1

            continue

        retained_proteins += 1

        pI_values.append(
            protein_pI
        )

        mw_values.append(
            protein_mw
        )

        # Save individual protein information
        protein_results.append(
            {
                "Genome_File": fasta_file.name,
                "Species": genome_name,
                "CDS_ID": record.id,
                "Protein_Length_aa": len(protein),
                "pI": protein_pI,
                "MolWt_Da": protein_mw
            }
        )

    # --------------------------------------------------------
    # No proteins retained
    # --------------------------------------------------------

    if retained_proteins == 0:

        print(
            "    WARNING: no proteins retained"
        )

        continue

    # --------------------------------------------------------
    # Genome-level statistics
    # --------------------------------------------------------

    mean_pI = np.mean(pI_values)

    median_pI = np.median(pI_values)

    sd_pI = np.std(
        pI_values,
        ddof=1
    )

    mean_mw = np.mean(mw_values)

    median_mw = np.median(mw_values)

    # --------------------------------------------------------
    # Bootstrap CI for genome-level mean pI
    # --------------------------------------------------------

    ci_lower, ci_upper = bootstrap_mean_ci(
        pI_values,
        repetitions=BOOTSTRAP_REPS,
        seed=RANDOM_SEED + number
    )

    # --------------------------------------------------------
    # Save genome result
    # --------------------------------------------------------

    genome_results.append(
        {
            "Genome_File": fasta_file.name,
            "Species": genome_name,

            "Avg_pI": mean_pI,
            "Median_Protein_pI": median_pI,
            "SD_Protein_pI": sd_pI,

            "pI_Bootstrap_95CI_Lower": ci_lower,
            "pI_Bootstrap_95CI_Upper": ci_upper,

            "Avg_MolWt_Da": mean_mw,
            "Median_Protein_MolWt_Da": median_mw,

            "Total_Proteins": retained_proteins,
            "Original_CDS_Count": total_cds
        }
    )

    # --------------------------------------------------------
    # QC information
    # --------------------------------------------------------

    qc_results.append(
        {
            "Genome_File": fasta_file.name,
            "Species": genome_name,
            "Original_CDS_Count": total_cds,
            "Retained_Proteins": retained_proteins,
            "Rejected_Total": (
                total_cds - retained_proteins
            ),

            "Too_Short_NT": qc["too_short_nt"],
            "Translation_Error": qc["translation_error"],
            "Internal_Stop": qc["internal_stop"],
            "Short_Protein": qc["short_protein"],
            "Nonstandard_AA": qc["nonstandard_aa"]
        }
    )

    print(
        f"    CDS: {total_cds:,} | "
        f"Proteins: {retained_proteins:,} | "
        f"Mean pI: {mean_pI:.4f} | "
        f"Mean MW: {mean_mw:,.2f}"
    )


# ============================================================
# CREATE DATAFRAMES
# ============================================================

genome_df = pd.DataFrame(
    genome_results
)

protein_df = pd.DataFrame(
    protein_results
)

qc_df = pd.DataFrame(
    qc_results
)


# ============================================================
# SAVE CSV FILES
# ============================================================

genome_csv = (
    FOLDER /
    "RECALCULATED_Genome_Summary.csv"
)

protein_csv = (
    FOLDER /
    "RECALCULATED_Individual_Proteins.csv"
)

qc_csv = (
    FOLDER /
    "RECALCULATED_QC.csv"
)

genome_df.to_csv(
    genome_csv,
    index=False
)

protein_df.to_csv(
    protein_csv,
    index=False
)

qc_df.to_csv(
    qc_csv,
    index=False
)


# ============================================================
# SAVE ONE EXCEL WORKBOOK
# ============================================================

excel_file = (
    FOLDER /
    "RECALCULATED_PROTEOME_ANALYSIS.xlsx"
)

with pd.ExcelWriter(
    excel_file,
    engine="openpyxl"
) as writer:

    genome_df.to_excel(
        writer,
        sheet_name="Genome_Summary",
        index=False
    )

    qc_df.to_excel(
        writer,
        sheet_name="QC",
        index=False
    )

    protein_df.to_excel(
        writer,
        sheet_name="Individual_Proteins",
        index=False
    )


# ============================================================
# PRINT FINAL SUMMARY
# ============================================================

print()
print("=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)

print(
    f"Genomes successfully analysed: "
    f"{len(genome_df)}"
)

print(
    f"Individual proteins analysed: "
    f"{len(protein_df):,}"
)

if len(genome_df) > 0:

    print()

    print(
        "Protein count range: "
        f"{genome_df['Total_Proteins'].min():,.0f} - "
        f"{genome_df['Total_Proteins'].max():,.0f}"
    )

    print(
        "Genome mean pI range: "
        f"{genome_df['Avg_pI'].min():.4f} - "
        f"{genome_df['Avg_pI'].max():.4f}"
    )

    print()

    # Look specifically for baumannii
    baumannii = genome_df[
        genome_df["Species"]
        .str.lower()
        .str.contains("baumannii", na=False)
    ]

    if len(baumannii) > 0:

        print("BAUMANNII CHECK:")
        print(
            baumannii[
                [
                    "Genome_File",
                    "Species",
                    "Avg_pI",
                    "Avg_MolWt_Da",
                    "Total_Proteins",
                    "pI_Bootstrap_95CI_Lower",
                    "pI_Bootstrap_95CI_Upper"
                ]
            ].to_string(index=False)
        )

print()
print("FILES CREATED:")
print(genome_csv)
print(protein_csv)
print(qc_csv)
print(excel_file)

print()
print("=" * 80)