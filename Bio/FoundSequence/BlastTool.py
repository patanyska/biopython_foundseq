# Copyright 2024 by Patricia Nogueira.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Code to invoke the Blast Tool over the internet.

This module provides code to work with the NCBI BLAST+ provided by EMBL-EBI
Details about API here https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?pageId=94147939#NCBIBLAST+HelpandDocumentation-RESTAPI
"""
from Bio import SeqIO
from Bio.Data import IUPACData
import requests
import json
import time
import re

EMBL_EBI_CREATE_JOB_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
EMBL_EBI_STATUS_JOB_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
EMBL_EBI_BLAST_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"

AMINO_ACIDS = set(IUPACData.extended_protein_letters)
NUCLEOTIDES = ["A", "C", "G", "T"]

"""Function to validate empty email
Variables: 
- email        email from user
"""
def validate_Empty_Email(email):
    return bool(email.strip())

"""Function to validate email format
Variables: 
- email        email from user
"""
def validate_Email_Format(email):
    regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,7}\b'
    return bool(re.fullmatch(regex, email))

"""Function to find variantes in sequence
Variables: 
- query         sequence in file
- subject       reference sequence to compare with query
"""
def find_Variants(query,subject):
    if not query or not subject:
        raise ValueError("Input sequences cannot be empty.")
    if len(query) != len(subject):
        raise ValueError("Input sequences must have the same length for simple mismatch detection.")

    variants = []
    for i, (q_base, s_base) in enumerate(zip(query, subject)):
        if q_base != s_base:
            variants.append({
                'original': s_base,
                'variation': q_base,
                'position': str(i + 1)
            })
    return variants

""" Validates if the first sequence in a FASTA file contains only nucleotides.
    Variables:
    - file          File .fasta uploaded with the nucleotide sequence
"""
def validate_Nucleotide_Sequence(sequence)-> bool:
    try:
        # Iterate through each character in the sequence
        for s in sequence:
            # Check if the uppercase version is not in our valid set
            if s.upper() not in NUCLEOTIDES:
                return False

        # If the loop completes without finding invalid characters, the sequence is valid
        return True

    except Exception as e:
        # Handle other potential parsing errors
        #print(f"An unexpected error occurred during parsing: {e}")
        return False

"""Validates if the first sequence in a FASTA file contains only valid amino acids.

    Uses the `VALID_AMINO_ACIDS` set defined globally in the script, which
    defaults to the 20 standard amino acids from Biopython's IUPACData.

    Args:
        file_handle: An open file handle to the FASTA file.

    Returns:
        True if the sequence contains only amino acids specified in VALID_AMINO_ACIDS,
        False otherwise. Also returns False if the file cannot be parsed or
        does not contain exactly one sequence. Returns True for an empty sequence.

    Raises:
        May re-raise exceptions from SeqIO during parsing if not handled.
"""
def validate_Protein_Sequence(sequence) -> bool:
    
    try:
        # Iterate through each character in the sequence
        for s in sequence:
           if s.upper() not in AMINO_ACIDS:
                return False

        # If the loop completes without finding invalid characters, the sequence is valid
        return True

    except Exception as e:
        # Handle other potential parsing errors
        return False

"""Function to call NCBI Blast+ API
Variables:
    -email        string    Valid email
    -program	  string	BLAST program to use to perform the search.
    -matrix	      string	Scoring matrix to be used in the search.
    -alignments	  int	    Maximum number of alignments displayed in the output.
    -scores	      int	    Maximum number of scores displayed in the output.
    -exp	      string	E-value threshold.
    -dropoff	  int	    Amount score must drop before extension of hits is halted.
    -match_scores string	Match/miss-match scores to generate a scoring matrix for for nucleotide searches.
    -gapopen	  int	    Penalty for the initiation of a gap.
    -gapext	      int	    Penalty for each base/residue in a gap.
    -filter	      string	Low complexity sequence filter to process the query sequence before performing the search.
    -seqrange	  string	Region of the query sequence to use for the search. Default: whole sequence.
    -gapalign	  boolean	Perform gapped alignments.
    -compstats	  string	Compositional adjustment or compositional statistics mode to use.
    -align	      int	    Alignment format to use in output.
    -stype	      string	Query sequence type. One of: dna, rna or protein.
    -sequence	  string	Query sequence. The use of fasta formatted sequence is recommended.
    -database	  list	    List of database names for search.

"""
def blast(email,
  program,
  matrix,
  alignments,
  scores,
  exp,
  dropoff,
  match_scores,
  gapopen,
  gapext,
  filter,
  seqrange,
  gapalign,
  compstats,
  align,
  stype,
  sequence,
  database):
    
    if not validate_Empty_Email(email):
        raise ValueError("The email is mandatory")
        
    if not validate_Email_Format(email):
          raise ValueError("The email format is not valid")
           
    
    programs = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
    matrixes = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250"]

    if program not in programs:
        program="blastp"
        
    if matrix not in matrixes:
        matrix='BLOSUM62'
        
    if len(sequence)==0:
        sequence=""
       
    if not sequence:
        raise ValueError("Sequence cannot be empty.")
    
    if(program=="blastp" and not validate_Protein_Sequence(sequence)):
         raise ValueError("Sequence with invalid amino-acids")
    
    if(program=="blastn" and not validate_Nucleotide_Sequence(sequence)):
         raise ValueError("Sequence with invalid nucleotides")

    files = {
            'email': email,
            'program': program,
            'matrix': matrix,
            'alignments': alignments if alignments is not None else '5',
            'scores': scores if scores is not None else '5',
            'exp': exp if exp is not None else '1e-3',
            'dropoff': dropoff if dropoff is not None else '0',
            'match-scores': match_scores if match_scores is not None else '50',
            'gapopen': gapopen if gapopen is not None else '-1',
            'gapext': gapext if gapext is not None else '-1',
            'filter': filter if filter is not None else 'F',
            'seqrange': seqrange if seqrange is not None else 'START-END',
            'gapalign': gapalign if gapalign is not None else 'true',
            'compstats': compstats if compstats is not None else 'F',
            'align': align if align is not None else '0',
            'stype': stype if stype is not None else 'protein',
            'sequence': sequence,
            'database': database if database is not None else 'uniprotkb_refprotswissprot',
        }

    try:
        # Request to create job
        job = requests.post(EMBL_EBI_CREATE_JOB_URL, files=files)
        job.raise_for_status()  # Raise an exception for HTTP errors

        job_id = job.text.strip()
        job_status = ""
        wait_time = 10  # Initial wait time in seconds
        max_attempts = 100 # Maximum number of status checks

        for attempt in range(max_attempts):
            status_url = EMBL_EBI_STATUS_JOB_URL + job_id
            job_status_request = requests.get(status_url)
            job_status_request.raise_for_status()
            job_status = job_status_request.text.strip()
           
            if job_status == "FINISHED":
                break
            elif job_status in ["NOT_FOUND", "FAILURE", "ERROR"]:
                raise ValueError(f"BLAST job failed with status: {job_status}")
            elif job_status == "RUNNING":
                time.sleep(wait_time)
                wait_time = min(wait_time * 2, 120) # Exponential backoff up to 2 minutes
            else:
                time.sleep(wait_time)

        else:
            raise TimeoutError(f"BLAST job did not finish within the maximum wait time.")

        # Get JSON results
        result_url = EMBL_EBI_BLAST_URL + job_id + '/json'
        response = requests.get(result_url)
        response.raise_for_status()
        json_data = response.json()  # Use response.json() directly

        return json_data

    except requests.exceptions.RequestException as e:
        raise ValueError(f"Error communicating with the BLAST API: {e}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Error decoding JSON response from the BLAST API: {e}")
    except Exception as e:
        raise ValueError(f"An unexpected error occurred: {e}")


    ''' #request to create job 
        job = requests.post(EMBL_EBI_CREATE_JOB_URL, files=files)

        #wait until job not finish and then get json file with data
        job_status=""
        while job_status!="FINISHED":
            job_status_request=requests.get(EMBL_EBI_STATUS_JOB_URL+job.text)
            job_status=job_status_request.text

            if(job_status=="NOT_FOUND" or job_status=="FAILURE" or job_status=="ERROR"):
                raise ValueError(f"Job had a problem and the status is: {job_status}")
                    
            time.sleep(60) #sleep for 60 seconds until Job is finished

        response=requests.get(EMBL_EBI_BLAST_URL+job.text+'/json')
        bytes_data = response.content
        json_data = json.loads(bytes_data.decode('utf-8'))
                    
        return json_data   
    except Exception as e:
        raise ValueError(f"An error occurred: {str(e)}")'''



