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
from Bio import SeqRecord, Seq


EMBL_EBI_CREATE_JOB_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
EMBL_EBI_STATUS_JOB_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
EMBL_EBI_BLAST_URL="https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"

"""Function to get match for sequence in Blast
Variables:

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

def blast(program,
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
        
            qry = SeqRecord(Seq(query), id="query")
            sbj = SeqRecord(Seq(subject), id="refseq") 
            for i in range(len(qry.seq)):
                if qry.seq[i]!=sbj.seq[i]:
                    m = {'original': sbj.seq[i], 'variation': qry.seq[i],'position':str(i+1)}
                    variants.append(m)
        return variants
    except:
        raise ValueError(
            f"A problem occurred getting variants from sequence")
