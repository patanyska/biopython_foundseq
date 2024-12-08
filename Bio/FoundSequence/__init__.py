# Copyright 2024 by Patricia Nogueira.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Code to invoke the FoundSequence Tool over the internet.

This module provides code to work with the several servics as:
 - Blast Tool
 - Translate Tool
 - FindProtein

 Invoking this module will allow get all sequence information in one call - function foundSequence
 The flow is:
 1 - Translate nucleotide sequence to protein sequence using Translate Tool - By ExPASy
 2 - Get biggest Open Reading Frame (ORF) from protein sequence
 3 - Send ORF to Blast Tool to discover witch protein is, possible varians and positions, and hit accession to UniProt- By EMBL_EBL Blast Tool
 4 - Get info about protein in Uniprot from hit accession - By UniProt

 If user wants, can get info only from separeted functions calling function from different 

 Variables:
   -file          bytes     valid .Fasta file
   -email         string    valid email
   -program	      string	BLAST program to use to perform the search.
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


import TranslateTool
import BlastTool
import FindProtein
import DrugBank


def foundSequence(file,
  email,
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
  database):
    
    dict={}
    raised = False

    try:
       expasy_result=TranslateTool.expasy_Translate_Tool(file)
       if(expasy_result==""):
            dict["expasy"] = {}
            dict["blast"]={}
            dict["uniprot"] = {}
       else:
        big_orf=TranslateTool.get_BigORF(expasy_result)
        if(big_orf==0):
            dict["expasy"] = {}
            dict["expasy"]['protein'] = expasy_result
            dict["expasy"]['bigORF']={}
        else:
            dict["expasy"] = {}
            dict["expasy"]['protein'] = expasy_result
            dict["expasy"]['bigORF']=big_orf
            
            blast_result=BlastTool.blast(email,program,matrix,alignments,scores,exp,dropoff,match_scores,gapopen,gapext,filter,seqrange,gapalign,compstats,align,stype,big_orf,database)
            struct_Blast=BlastTool.read_Json(blast_result)
            if(len(struct_Blast)==0):
                dict["blast"]={}
                dict["blast"]["variants"]={}
                dict["uniprot"] = {}
            else:
                dict["blast"]={}
                dict["blast"]["hit_id"]=struct_Blast[0]
                dict["blast"]["hit_def"]=struct_Blast[1]
                dict["blast"]["hit_acc"]=struct_Blast[2]
                dict["blast"]["hit_uni_de"]=struct_Blast[3]
                dict["blast"]["hit_uni_os"]=struct_Blast[4]
                dict["blast"]["hsp_gaps"]=struct_Blast[5]
                dict["blast"]["hsp_align_len"]=struct_Blast[6]
                dict["blast"]["hsp_qseq"]=struct_Blast[7]
                dict["blast"]["hsp_hseq"]=struct_Blast[8]
                variants=BlastTool.find_Variants(struct_Blast[7],struct_Blast[8])
                if(len(variants)==0):
                    dict["blast"]["variants"]={}
                    dict["uniprot"] = {}
                else:
                    dict["blast"]["variants"]=variants
                    uniprot_result=FindProtein.found_Uniprot_Protein(struct_Blast[2])
                    struct_Uniprot=FindProtein.read_Json(uniprot_result,variants)
                    if(len(struct_Uniprot)==0):
                        dict["uniprot"] = {}
                    else:
                        dict["uniprot"] = {}
                        dict["uniprot"]["entry_type"]=struct_Uniprot[0]["entry_type"]
                        dict["uniprot"]["scientific_name"]=struct_Uniprot[0]["scientific_name"]
                        dict["uniprot"]["common_name"]=struct_Uniprot[0]["common_name"]
                        dict["uniprot"]["taxon_id"]=struct_Uniprot[0]["taxon_id"]
                        dict["uniprot"]["lineage"]=struct_Uniprot[0]["lineage"]
                        dict["uniprot"]["full_name"]=struct_Uniprot[0]["full_name"]
                        dict["uniprot"]["short_name"]=struct_Uniprot[0]["short_name"]
                        dict["uniprot"]["protein_function"]=struct_Uniprot[0]["protein_function"]
                        dict["uniprot"]["catalytic_activity"]=struct_Uniprot[0]["catalytic_activity"]
                        dict["uniprot"]["disease"]=struct_Uniprot[0]["disease"]
                        dict["uniprot"]["acronym"]=struct_Uniprot[0]["acronym"]
                        dict["uniprot"]["description"]=struct_Uniprot[0]["description"]
                       
                        if(struct_Uniprot[0]["disease"]!="-"):
                            drugbank_result=DrugBank.found_Drug(struct_Uniprot[0]["disease"])
                            if(len(drugbank_result)==0):
                                dict["drugbank"] = {}
                            else:
                                dict["drugbank"] = {}
                                drugs={}
                                for d in drugbank_result:
                                    drug = {'name':d[0], 
                                            'description':d[1],
                                            'indication':d[2],
                                            'toxicity':d[3],
                                            'product_name':d[4],
                                            'labeller':d[5],
                                            'dosage_form':d[8],
                                            'strength':d[9],
                                            'route':d[10]}
                                    drugs.append(drug)
                                dict["drugbank"]=drugs     
        return dict
    except Exception as e:
        return str(e)

def main():
    my_file="C:\\Users\\PC\\Desktop\\Fastas\\mutseq10.fasta"
    program='blastp' 
    matrix='BLOSUM62'
    alignments='50'
    scores='5'
    exp='1e-3'
    dropoff='0'
    match_scores='50'
    gapopen='-1'
    gapext='-1'
    filter='F'
    seqrange= 'START-END'
    gapalign= 'true'
    compstats= 'F'
    align= '0' 
    stype= 'protein'
    database='uniprotkb_refprotswissprot'
    foundSequence(my_file,program,matrix,alignments,scores,exp,dropoff,match_scores,gapopen,gapext,filter,seqrange,gapalign,compstats,align,stype,database)

if __name__ == "__main__":
         main()