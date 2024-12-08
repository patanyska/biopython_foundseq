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

import requests
import json
import time
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


EMBL_EBI_CREATE_JOB_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
EMBL_EBI_STATUS_JOB_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
EMBL_EBI_BLAST_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"

"""Function to get match for sequence in Blast
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
    
    if(validate_Empty_Email(email)==False):
                 raise ValueError(
            f"The email is mandatory"
        )
    
    if(validate_Email_Format(email)==False):
                 raise ValueError(
            f"The email format is not valid"
        )
    
    try:
        programs = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
        matrixes=["BLOSUM45", "BLOSUM50","BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250"]
        if program not in programs:
            program="blastp"
            
        if matrix not in matrixes:
            matrix='BLOSUM62'
            
        if len(sequence)==0:
            sequence=""

        files = {
                    'email': email,
                    'program': 'blastp' if program is None else program,
                    'matrix': 'BLOSUM62' if matrix is None else matrix,
                    'alignments': '50' if alignments is None else alignments,
                    'scores': '5' if scores is None else scores,
                    'exp': '1e-3' if exp is None else exp,
                    'dropoff': '0' if dropoff is None else dropoff,
                    'match-scores': '50' if match_scores is None else match_scores,
                    'gapopen': '-1' if gapopen is None else gapopen,
                    'gapext': '-1' if gapext is None else gapext,
                    'filter': 'F' if filter is None else filter,
                    'seqrange': 'START-END' if seqrange is None else seqrange,
                    'gapalign': 'true' if gapalign is None else gapalign,
                    'compstats': 'F' if compstats is None else compstats,
                    'align': '0' if align is None else align,
                    'stype': 'protein' if stype is None else stype,
                    'sequence':sequence,
                    'database': 'uniprotkb_refprotswissprot' if database is None else database,
                }


        #request to create job 
        job = requests.post(EMBL_EBI_CREATE_JOB_URL, files=files)

        #wait until job not finish and then get json file with data
        job_status=""
        while job_status!="FINISHED":
            job_status_request=requests.get(EMBL_EBI_STATUS_JOB_URL+job.text)
            job_status=job_status_request.text

            if(job_status=="NOT_FOUND" or job_status=="FAILURE" or job_status=="ERROR"):
                raise ValueError(
                f"Job had a problem and the status is {job_status}"
            )
                    
            time.sleep(60) #sleep for 60 seconds until Job is finished

        response=requests.get(EMBL_EBI_BLAST_URL+job.text+'/json')
        bytes_data = response.content
        json_data = json.loads(bytes_data.decode('utf-8'))
                    
        return json_data   
    except:
        raise ValueError(
            f"A problem occurred getting match from Blast")

"""Function to find variantes in sequence
Variables: 
- query         sequence in file
- subject       reference sequence to compare with query
"""
def find_Variants(query,subject):
    try:
        variants=[]
        if(query!="" and subject!=""):
        
            qryRecord = SeqRecord(Seq(query), id="query")
            sbjRecord = SeqRecord(Seq(subject), id="refseq") 
            for i in range(len(qryRecord.seq)):
                if qryRecord.seq[i]!=sbjRecord.seq[i]:
                    m = {'original': sbjRecord.seq[i], 'variation': qryRecord.seq[i],'position':str(i+1)}
                    variants.append(m)
        return variants
    except Exception as e:
        raise ValueError(
            f"A problem occurred getting variants from sequence")

"""Function to validate empty email
Variables: 
- email        email from user
"""
def validate_Empty_Email(email):
    if(email!=""):
        return True
    else:
        return False

"""Function to validate email format
Variables: 
- email        email from user
"""
def validate_Email_Format(email):
    regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,7}\b'
    # pass the regular expression
    # and the string into the fullmatch() method
    if(re.fullmatch(regex, email)):
        return True
    else:
        return False