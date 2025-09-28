from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd 
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import re


#%%

# Counts the number of sequences in a fasta file
def count_sequences(fasta):
    count = sum(1 for _ in SeqIO.parse(fasta, "fasta"))
    #print("Number of sequences:", count)
    return count

#%%

# Change the order of the sequences in a FASTA file
def reorder_sequences(input_fasta, output_fasta, reference_id, save_to):
    # Read all sequences
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    #print('Read sequences.')
    # Separate the reference genome
    ref_seq = None
    other_seqs = []

    for seq in sequences:
        if seq.id == reference_id:
            ref_seq = seq
        else:
            other_seqs.append(seq)
    #print('Separated reference genome and other sequences.')
    # Write the sequences back with the reference genome first
    if ref_seq:
        os.chdir(save_to)
        with open(output_fasta, "w") as output_handle:
            SeqIO.write([ref_seq] + other_seqs, output_handle, "fasta")
        #print(f"Reordered FASTA saved as {output_fasta}")
        os.chdir('..')
        return f"Reordered FASTA saved as {output_fasta}"
    else:
        return None
        #print("Reference genome not found in the FASTA file.")
        
#%%

# Remove duplicates in a FASTA file (try to have refseq as first entry)

def remove_dupes(input_fasta, output_fasta, save_to):
    unique_sequences = {}
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = str(record.seq)
        if sequence not in unique_sequences:
            unique_sequences[sequence] = record.id
    os.chdir(save_to)      
    with open(output_fasta, 'w') as file:
        for sequence, seq_id in unique_sequences.items():
            file.write(f'>{seq_id}\n{sequence}\n')
    os.chdir('..')
    return f"Removed duplicates and saved as {output_fasta}"

#%%

# For proteins: remove sequences with ambiguous (X) nucl, stop codons (*) or excessive gaps (-)

def remove_ambiguous(input_fasta, output_fasta, save_to):
    os.chdir(save_to)
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if "X" not in record.seq and "*" not in record.seq and record.seq.count("-") < len(record.seq) * 0.05:
                SeqIO.write(record, out_f, "fasta")

    print(f"Ambiguous sequences removed. Cleaned sequences saved to {output_fasta}")
    
#%%

def extract_sequence(fasta_file, output, sequence_ID, save_to):
    os.chdir(save_to)
    found = False
    with open(output, 'w') as extract:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id == sequence_ID:
                extract.write('>'+record.id+'\n')
                extract.write(str(record.seq.upper())+'\n')
                found = True
                break
    if not found:
        os.chdir('..')
        return None
    
    os.chdir('..')
    return f'Sequence extracted.'
    #print(f'Sequence extracted to {output}.')
    
#%%

# Export Primer IDs and sequences from a dictionary

def export_primers_to_fasta(primer_list, output_file):
    """
    Exports a list of primer dictionaries to a FASTA file using Primer ID as header and Sequence as sequence.

    Args:
        primer_list (list): List of primer dictionaries.
        output_file (str): Path to output FASTA file.
    """
    records = []
    for primer in primer_list:
        seq_id = primer.get("Primer ID", "unknown_id")
        seq_str = primer.get("Sequence", "")
        if seq_str:  # Avoid empty sequences
            record = SeqRecord(Seq(seq_str), id=seq_id, description="")
            records.append(record)
    
    SeqIO.write(records, output_file, "fasta")

#%%
def save_to_csv(data, path):
    df = pd.DataFrame(data)
    df.to_csv(path, index=False)

def primer_venn_diagram(file1, file2, label1="File1", label2="File2"):
    # Read sequences (use IDs or sequences for comparison)
    primers1 = {str(record.seq) for record in SeqIO.parse(file1, "fasta")}
    primers2 = {str(record.seq) for record in SeqIO.parse(file2, "fasta")}

    # Make the Venn diagram
    venn2([primers1, primers2], set_labels=(label1, label2))
    plt.title("Primer Overlap")
    plt.show()

def fetch_doi_dict(nbib_file):
    doi_dict = {}
    PMID_dict = {} # some articles don't have DOIs so we can try to save PMID instead (always has PMID bc it was fetched from PubMed db)
    pmid = None # so it exists outside of the loop for safety
    has_doi = False # initialize condition
        
    with open(nbib_file, 'r', encoding='utf-8') as file:
        for line in file:
            if line.startswith('PMID- '): 
                pmid = line.strip().split('- ')[1]
                has_doi = False # reset because new article
                    
            elif line.startswith('AID - ') and '[doi]' in line:
                match = re.search(r'10\.\d{4,9}/\S+', line)
                if match:
                    doi = match.group(0)
                    doi_dict.setdefault(f"https://doi.org/{doi}", [])
                    has_doi = True # mark article as having DOI
                        
            elif line.strip() == '': # spots the empty line inbetween entries
                if pmid and not has_doi:
                    PMID_dict.setdefault(f"https://pubmed.ncbi.nlm.nih.gov/{pmid}", [])
                    pmid = None
                    has_doi = False # reset back to default second time for safe measures


    return doi_dict, PMID_dict


#%%
# os.chdir(r"D:\MESTRADO\DISSERTACAO\Sequencias\EBOLA")
# ref_fasta = open('reference_genome_noalign.fasta', 'w')
# file = r"D:\MESTRADO\DISSERTACAO\Sequencias\EBOLA\reordered_sequences.fasta"
# for record in SeqIO.parse(file, 'fasta'):
#     if record.id == 'NC_002549.1':
#         ref_fasta.write('>'+record.id+'\n')
#         ref_fasta.write(str(record.seq.upper())+'\n')
# ref_fasta.close()

#%%
# sequence = next(SeqIO.parse(r"D:\MESTRADO\DISSERTACAO\Sequencias\EBOLA\reference_genome_noalign.fasta",'fasta'))

# region = sequence.seq[6038:8068]

# protein = region.translate()
# print(f"Translated Protein: {protein}")

#remove_ambiguous(r"\\wsl.localhost\Ubuntu\home\anasfplima\clustered_gp.fasta", 'test.fasta', r'D:\MESTRADO\DISSERTACAO\Sequencias\EBOLA\ZEBOVgp4\spike glycoprotein all sequences ncbi')

# count_sequences(r"D:\MESTRADO\DISSERTACAO\Sequencias\HIV\HIV1\3D MODELLING\gp160\all protein sequences ncbi\Alinhamentos\clustered_gp.fasta")

#extract_sequence(r"D:\MESTRADO\DISSERTACAO\Sequencias\EBOLA\Alinhamentos\ebolaMAFFTnodupes.fasta", 'refseqaligned.fasta', 'NC_002549.1', r'D:\MESTRADO\DISSERTACAO\Sequencias\EBOLA')