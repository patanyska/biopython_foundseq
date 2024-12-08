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
import DrugBankTool


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
        if(big_orf==""):
            dict["expasy"] = {}
            dict["expasy"]['protein'] = expasy_result
            dict["expasy"]['bigORF']={}
        else:
            dict["expasy"] = {}
            dict["expasy"]['protein'] = expasy_result
            dict["expasy"]['bigORF']=big_orf
            
            blast_result=BlastTool.blast(email,program,matrix,
                                         alignments,scores,exp,
                                         dropoff,match_scores,
                                         gapopen,gapext,filter,
                                         seqrange,gapalign,compstats,
                                         align,stype,big_orf,database)
            struct_Blast=read_Blast_Json(blast_result)
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
                    struct_Uniprot=read_Uniprot_Json(uniprot_result,variants)
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
                            drugbank_result=DrugBankTool.found_Drug(struct_Uniprot[0]["disease"])
                            if(len(drugbank_result)==0):
                                dict["drugbank"] = {}
                            else:
                                dict["drugbank"] = {}
                                drugs={}
                                for d in drugbank_result:
                                    drug = {'id':d[1],
                                            'name':d[0],
                                            'description':d[1],
                                            'state':d[5],
                                            'indication':d[6],
                                            'product_name':d[9],
                                            'labeller':d[10],
                                            'route':d[11],
                                            'country':d[12]}
                                    drugs.append(drug)
                                dict["drugbank"]=drugs
        return dict
    except Exception as e:
        return str(e)

def read_Blast_Json(json_file):
        try:        
            type=""
            result=[]
            if not 'hits' in json_file or len(json_file['hits']) == 0:
                return result
            else:
                
                struct_hit={}
                for hit in json_file['hits']:
                    #in case of protein
                    if 'hit_uni_os' in hit:
                        if hit["hit_uni_os"]== "Homo sapiens":
                            type="protein"
                            struct_hit=hit 
                            break 
                    else:
                        if hit["hit_db"]== "EM_HUM":
                            type="nucleotide"
                            struct_hit=hit 
                            break 
                
                hit_id=struct_hit["hit_id"]    
                hit_def=struct_hit["hit_def"] 
                hit_acc=struct_hit["hit_acc"].split("-")[0] #in case of protein is a uniprot id ex:P12345
                if type=="protein":
                    hit_uni_de=struct_hit["hit_uni_de"] 
                    hit_uni_os=struct_hit["hit_uni_os"] 
                else:
                    hit_uni_de="NA"
                    hit_uni_os="NA"
                hsp_gaps=struct_hit["hit_hsps"][0]["hsp_gaps"]
                hsp_align_len=struct_hit["hit_hsps"][0]["hsp_align_len"]
                hsp_qseq=struct_hit["hit_hsps"][0]["hsp_qseq"] #query
                hsp_hseq=struct_hit["hit_hsps"][0]["hsp_hseq"] #record
                
                result.append(hit_id)
                result.append(hit_def)
                result.append(hit_acc)
                result.append(hit_uni_de)
                result.append(hit_uni_os)
                result.append(hsp_gaps)
                result.append(hsp_align_len)
                result.append(hsp_qseq)
                result.append(hsp_hseq)
                
                return result
        except Exception as ex:
            return ex
        
def read_Uniprot_Json(json_file,variants):
    struct=[]
    try:
        entryType=""
        scientificName=""
        commonName=""
        taxonId=""
        lineage=""
        fullName=""
        shortName=""
        protein_function=""
        catalytic_activity=""
        struct_evidences_id=""
        disease=""
        acronym=""
        disease_description=""

    
        entryType=json_file["entryType"]

        json_file["features"]

        for f in json_file["features"]:
            if(f['type']=="Natural variant" or f['type']=="Mutagenesis"):
                for v in variants:
                    if(f["location"]["start"]["value"]==int(v['position']) and f['alternativeSequence']['originalSequence']==v["original"]):
                        for aseq in f['alternativeSequence']['alternativeSequences']:
                            if(aseq==v["variation"]):
                                struct_evidences_id=f["evidences"]["id"]

        scientificName=json_file["organism"]["scientificName"]

        if "commonName" in json_file["organism"].keys():
            commonName=json_file["organism"]["commonName"]
        else:
            commonName="-"

        taxonId=str(json_file["organism"]["taxonId"])

        for l in json_file["organism"]["lineage"]:
            lineage+=l+"; "

        fullName=json_file["proteinDescription"]["recommendedName"]["fullName"]["value"]

        if "shortNames" in json_file["proteinDescription"]["recommendedName"]:
            for sn in json_file["proteinDescription"]["recommendedName"]["shortNames"]:
                shortName+=sn["value"]+"; "
        else:
            shortName="-"

        for c in json_file["comments"]:
            if(c['commentType']=="FUNCTION"):
                for t in c["texts"]:
                    protein_function+=t["value"]+"\n"
            if(c['commentType']=="CATALYTIC ACTIVITY"):
                catalytic_activity=c["reaction"]["name"]


         
            
            if(c['commentType']=="DISEASE" and struct_evidences_id!=""):
                if ("disease" in c):
                    for ev in c['disease']['evidences']:
                        if(ev['id']==struct_evidences_id):
                            disease=c['disease']['diseaseId']
                            acronym=c['disease']['acronym']
                            disease_description=c['disease']['description']
            else:
                disease="-"
                acronym="-"
                disease_description="-"
                                            
        prot = {'entry_type':entryType,
                'scientific_name':scientificName,
                'common_name':commonName,
                'taxon_id':taxonId,
                'lineage':lineage,
                'full_name':fullName,
                'short_name':shortName,
                'protein_function':protein_function,
                'catalytic_activity':catalytic_activity,
                'disease':disease,
                'acronym':acronym,
                'disease_description':disease_description}
        struct.append(prot)
      

     
                                    
        return struct
    except:
         raise ValueError(f"A problem occurred during reading UniProt json"
        )
        