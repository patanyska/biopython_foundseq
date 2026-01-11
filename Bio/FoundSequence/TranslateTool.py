# Copyright 2024 by Patricia Nogueira All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
"""Code to invoke the Expasy Translate Tool over the internet.

This module provides code to work with the Expasy Translate Tool
and return protein sequence formated in frames 5'3' and 3'5', and also detect the big ORF in all sequence
API Documentation: https://web.expasy.org/translate/programmatic_access.html
Variables:

    -file valid .Fasta file

"""

import os
import requests
from Bio import SeqIO #pip install Bio
import tempfile



EXPASY_URL="https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
NUCLEOTIDES = ["A", "C", "G", "T"]


""" Function to validate file format to guarantee that the uploaded file is .fasta extension
    Variables: 
    - file          File .fasta uploaded with the nucleotide sequence
"""
def validate_FileFormat(file):
    if not isinstance(file, str):
        return False
    return file.endswith(".fasta")

""" Function to validate if fasta file is empty
    parameters:
    - file          File .fasta uploaded with the nucleotide sequence
"""   
def validate_FileEmpty(file):
    try:
        SeqIO.read(file, "fasta")
        return True
    except Exception as e:
        return False
    
""" Function to validate file content to guarantee that all nucleotides are valid
    Variables:
    - file          File .fasta uploaded with the nucleotide sequence
    - web           true if is called for a web aplication; false if is called from desktop application
"""
def validate_Nucleotide_Sequence(file,web):
    if not web:   
        sequence = SeqIO.read(file, "fasta")
        for line in sequence:
            trim_line=str.strip(line)
            list_line=[]
            for l in trim_line:
                list_line.append(l)
                count=1
                for nuc in list_line:
                    if nuc.upper() not in NUCLEOTIDES:
                        return False
                    count=count+1
        return True
    else:
        file.seek(0)
        sequence_record=SeqIO.read(file,"fasta")
        sequence = sequence_record.seq
        if not sequence:
            return False # Sequence content itself is empty
        
        seq_upper = str(sequence).upper()
        for nucleotide in seq_upper:
            if nucleotide not in NUCLEOTIDES:
                return False # Invalid nucleotide found

            # If all checks pass
            return True

""" Function to invoke Expasy translate tool to protein
        Variables:
        - file  .FASTA   File .fasta uploaded with the nucleotide sequence
        - web   true if is called for a web aplication; false if is called from desktop application
"""
def expasy_Translate_Tool(file,web):
    
    if not web:   
        if not validate_FileFormat(file):
                raise ValueError(f"The file has a wrong format")
         
        if not validate_FileEmpty(file):
                 raise ValueError(f"The file cannot be empty (0 bytes).")

        if not validate_Nucleotide_Sequence(file,web):
                 raise ValueError(f"The file has invalid nucleotides.")
        
        
        sequence = SeqIO.read(file, "fasta")
        response = requests.post(EXPASY_URL,
                            data={
                                "dna_sequence": str(sequence.seq),
                                "output_format": "fasta"
                            })
        response.raise_for_status()
        output=response.content.decode("utf-8")
        return output
    else:
        with tempfile.NamedTemporaryFile(mode='w+', delete=True,encoding='utf-8',
             prefix='temp_fasta_', suffix='.fasta') as temp_fasta_file:
            
            for chunk in file.chunks():
                temp_fasta_file.write(chunk.decode('utf-8'))
                      
            temp_fasta_file.flush()
            if os.fstat(temp_fasta_file.fileno()).st_size == 0:
                raise ValueError("The file cannot be empty (0 bytes).")
            
            if not validate_Nucleotide_Sequence(temp_fasta_file,web):
                raise ValueError(f"The file has invalid nucleotides.")
            
            try:
                temp_fasta_file.seek(0)
                sequence_record=SeqIO.read(temp_fasta_file,"fasta")
                sequence = sequence_record.seq
                response = requests.post(
                                    EXPASY_URL,
                                    data={
                                        "dna_sequence": str(sequence),
                                        "output_format": "fasta"
                                    })
                response.raise_for_status()
                EXPASY_OUTPUT=response.content.decode("utf-8")
                return EXPASY_OUTPUT
            except requests.exceptions.RequestException as e:
                raise ValueError(f"Error communicating with the Expasy API: {e}")
            except Exception as e:
                raise ValueError(f"An unexpected error occurred: {e}")
                   
""" Function to read the protein file generated and save all open reading frames (ORFs) found in the sequence
      - protein  string  protein retrieved by Expasy the previous function
"""
def get_BigORF(protein):
    orfLength=0
    lastORFLength=0
    orf=""
    bigOrf=""
    start_orf=False
    clean_protein=protein.replace("\n","")
    if(len(str(clean_protein)))>0:
        for char in str(clean_protein):
            if((char=="M" or start_orf==True) and char!="-" and char!=">"):
                if(start_orf==False):
                    start_orf=True
                orfLength=orfLength+1
                orf=orf+char
            else:
                start_orf=False
                if(orfLength>lastORFLength):
                    lastORFLength=0
                    lastORFLength=orfLength
                    bigOrf=orf
                    orfLength=0
                    orf=""
                else:
                    orf=""
                    orfLength=0
        
    return bigOrf