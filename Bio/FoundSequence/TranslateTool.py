# Copyright 2024 by Patricia Nogueira All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
"""Code to invoke the Expasy Translate Tool over the internet.

This module provides code to work with the Expasy Translate Tool
and return protein sequence formated in frames 5'3' and 3'5', and also detect the big ORF in all sequence
Variables:

    -file valid .Fasta file

"""

import os
import requests
from Bio import SeqIO #pip install Bio
from datetime import datetime
from pathlib import PurePath
from Bio import SeqIO



EXPASY_URL="https://web.expasy.org/cgi-bin/translate/dna2aa.cgi"
NUCLEOTIDES = ["A", "C", "G", "T"]


""" Function to validate file format to guarantee that the uploaded file is .fasta extension
    Variables: 
    - file          File .fasta uploaded with the nucleotide sequence
"""
def validate_FileFormat(file):
   
    if (file.endswith(".fasta")):
        return True
    else:
        return False
""" Function to validate if fasta file is empty
    parameters:
    - file          File .fasta uploaded with the nucleotide sequence
"""   
def validate_FileEmpty(file):
    try:
        SeqIO.read(file, "fasta")
        return True
    except:
        return False
    
""" Function to validate file content to guarantee that all nucleotides are valid
    Variables:
    - file          File .fasta uploaded with the nucleotide sequence
"""
def validate_Nucleotide_Sequence(file):
    sequence = SeqIO.read(file, "fasta")
    for line in sequence:
        trim_line=str.strip(line)
        list_line=[]
        for l in trim_line:
            list_line.append(l)
            count=1
            for nuc in list_line:
                if (nuc.upper() not in NUCLEOTIDES):
                    return False
                count=count+1
    return True


""" Function to invoke Expasy translate tool to protein
        Variables:
        - file  .FASTA   File .fasta uploaded with the nucleotide sequence
"""
def expasy_Translate_Tool(file,web):
    
    if(web==False):   
        if(validate_FileFormat(file)==False):
                raise ValueError(
            f"The file has a wrong format"
        )
         
        if(validate_FileEmpty(file)==False):
                 raise ValueError(
            f"The file cannot be empty"
        )

        if(validate_Nucleotide_Sequence(file)==False):
                 raise ValueError(
            f"The file has invalid nucleotides."
        )
        
        
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
        mainpath=PurePath(__file__).parent
        folder_path=str(mainpath) + "\\temp_files\\"
        isExist = os.path.exists(folder_path)
        if not isExist:
             os.makedirs(folder_path)

        now=datetime.now()
        date_time = now.strftime("%d%m%Y_%H%M%S")
        FASTA_PATH=os.path.join(str(mainpath) + "\\temp_files\\","sequence_fasta"+date_time+".fasta")
        with open(FASTA_PATH, 'wb+') as destination:
            for chunk in file.chunks():
                destination.write(chunk)

        if(validate_FileEmpty(FASTA_PATH)==False):
                 raise ValueError(
            f"The file cannot be empty"
        )

        if(validate_Nucleotide_Sequence(FASTA_PATH)==False):
                 raise ValueError(
            f"The file has invalid nucleotides."
        )
        sequence = SeqIO.read(FASTA_PATH, "fasta")
        response = requests.post(
                            EXPASY_URL,
                            data={
                                "dna_sequence": str(sequence.seq),
                                "output_format": "fasta"
                            })
        response.raise_for_status()
        EXPASY_OUTPUT=response.content.decode("utf-8")
        return EXPASY_OUTPUT
                   
                    
        
 
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