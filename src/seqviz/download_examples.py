# scripts/download_examples.py
"""
One-time script to download cytochrome b sequences from NCBI.
Run this once to populate your examples/ folder.
Accessions from: Awad et al. (2015) Iran J Vet Res 16(2):218-222
"""

from Bio import Entrez, SeqIO
import os

Entrez.email = "your.email@example.com"  # NCBI requires this — they rate-limit by email

ACCESSIONS = {
    "gallus_gallus_cytb":      "KF964328",
    "columba_livia_cytb":      "KF964326",
    "coturnix_japonica_cytb":  "KF964327",
    "streptopelia_senegalensis_cytb": "KF964325",
}

os.makedirs("examples", exist_ok=True)

for filename, accession in ACCESSIONS.items():
    filepath = f"examples/{filename}.fasta"
    print(f"Downloading {accession} → {filepath}")

    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text"
    )
    record = SeqIO.read(handle, "fasta")
    handle.close()

    SeqIO.write(record, filepath, "fasta")
    print(f"  ✓ {len(record.seq)} bp written")

print("\nAll sequences downloaded.")
