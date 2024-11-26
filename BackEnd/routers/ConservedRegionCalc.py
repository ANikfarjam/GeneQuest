import subprocess
from Bio import AlignIO

def run_mafft(input_fasta, output_fasta):
    """
    Run MAFFT to align sequences in a FASTA file, allowing unusual symbols.
    """
    command = f"mafft --auto --anysymbol {input_fasta} > {output_fasta}"
    subprocess.run(command, shell=True, check=True)

def define_cons_region(sequence_list):
    """
    Discover the conservative region in multiple aligned sequences.
    """
    if not sequence_list or not all(len(seq) == len(sequence_list[0]) for seq in sequence_list):
        raise ValueError("All sequences must have the same length.")
    
    # Transpose the sequences to compare column-wise
    transposed = zip(*sequence_list)
    
    # Collect positions where all characters match
    conservative_region = []
    for column in transposed:
        if all(char == column[0] for char in column):
            conservative_region.append(column[0])
        else:
            conservative_region.append('-')  # Mark non-conservative regions with a placeholder
    
    return ''.join(conservative_region)

def calculate_conserved_region(sequence_data):
    """
    Accepts raw sequence data, aligns sequences, and calculates the conservative region.
    """
    # Create temporary FASTA files
    input_fasta = "temp_input.fasta"
    output_fasta = "temp_output.fasta"

    with open(input_fasta, "w") as f:
        for i, record in enumerate(sequence_data):
            sequence = str(record.seq).replace('U', 'T')  # Convert RNA to DNA
            f.write(f">seq{i + 1}\n{sequence}\n")

    # Align sequences
    run_mafft(input_fasta, output_fasta)

    # Parse aligned sequences
    alignment = AlignIO.read(output_fasta, "fasta")
    aligned_sequences = [str(record.seq) for record in alignment]

    # Find the conservative region
    conservative_region = define_cons_region(aligned_sequences)

    return conservative_region
