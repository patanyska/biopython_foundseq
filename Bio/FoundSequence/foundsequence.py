# Copyright 2024 by Patricia Nogueira.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Code to invoke the several services to get info about a sequence

- FASTA file

"""

import os
import requests
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ExPASy import TranslateTool
from Blast import EMBL_EBI
from UniProt import FindProtein


#region functions to call isolated services

def Expasy_Translate_Tool(file):
    return TranslateTool.expasy_Translate_Tool(file)



def Expasy_Big_ORFl(protein):
    return TranslateTool.get_BigORF(protein)



def Blast_EMBLEBI_Blast(program,
                        matrix,
                        alignments,
                        scores,
                        exp,
                        dropoff,
                        match_scores,
                        gapopne,gapext,
                        filter,
                        seqrange,
                        gapalign,
                        compstats,
                        align,
                        stype,
                        sequence,
                        database):
    
    return EMBL_EBI.blast(program,
                        matrix,
                        alignments,
                        scores,
                        exp,
                        dropoff,
                        match_scores,
                        gapopne,gapext,
                        filter,
                        seqrange,
                        gapalign,
                        compstats,
                        align,
                        stype,
                        sequence,
                        database)
#endregion

#region function to call all services and return one result will all info

def Found_All_About_Sequence(file):

    #definition of dictionary to add all info from different services
    dict={}
    dict["expasy"] = {}
    dict["blast"]={}
    dict["uniprot"] = {}

   
    expasy_Result=TranslateTool.expasy_Translate_Tool(file)
    if(len(expasy_Result)!=0):
        big_ORF=TranslateTool.get_BigORF(expasy_Result)
        if(len(big_ORF)!=0):
            dict["expasy"]['protein'] = expasy_Result
            dict["expasy"]['bigORF']=big_ORF
            blast_Result=EMBL_EBI.blast("blastp",
                                        "BLOSUM62",
                                        "50",
                                        "5",
                                        "1e-3",
                                        "0",
                                        "50",
                                        "-1",
                                        "-1",
                                        "F",
                                        "START-END",
                                        "true",
                                        "F",
                                        "0",
                                        "protein",
                                        ">EMBOSS_001\n"+big_ORF,
                                        "uniprotkb_refprotswissprot")
            
            struct_Blast=read_Blast_Json(blast_Result)
            if(len(struct_Blast)>0):
                dict["blast"]["hit_id"]=struct_Blast[0]
                dict["blast"]["hit_def"]=struct_Blast[1]
                dict["blast"]["hit_acc"]=struct_Blast[2]
                dict["blast"]["hit_uni_de"]=struct_Blast[3]
                dict["blast"]["hit_uni_os"]=struct_Blast[4]
                dict["blast"]["hsp_gaps"]=struct_Blast[5]
                dict["blast"]["hsp_align_len"]=struct_Blast[6]
                dict["blast"]["hsp_qseq"]=struct_Blast[7]
                dict["blast"]["hsp_hseq"]=struct_Blast[8]
                variants=find_Sequence_Variants(struct_Blast[7],struct_Blast[8])
                if(len(variants)==0):
                    dict["blast"]["variants"]={}
                    dict["uniprot"] = {}
                    dict["drugbank"] = {}
                                
                else:
                    dict["blast"]["variants"]=variants

                    uniprot_result=FindProtein.uniprot(struct_Blast[2])
                    struct_Uniprot=read_Uniprot_Json(uniprot_result,variants)
                    if(len(struct_Uniprot)>0):
                        dict["uniprot"]["protein_function"]=struct_Uniprot[0]["protein_function"]
                        dict["uniprot"]["catalytic_activity"]=struct_Uniprot[0]["catalytic_activity"]
                        dict["uniprot"]["disease"]=struct_Uniprot[0]["disease"]
                        dict["uniprot"]["acronym"]=struct_Uniprot[0]["acronym"]
                        dict["uniprot"]["description"]=struct_Uniprot[0]["description"]
          
            
        else:
            dict["expasy"]['protein'] = "TWas not possible to find ORF's in the protein."+"\n\n\n"+expasy_Result
            dict["expasy"]['bigORF']=""
     
    else:
        raise ValueError(
            f"The Expasy Translate tool was unable to convert the sequence to a protein. The translation result is empty. 
            Please check that the entered sequence is correct and try again.")
        
#endregion

#region Helper Functions
def read_Blast_Json(json_file):        
        f = open(json_file)
        data = json.load(f)
        result=[]
        if not 'hits' in data or len(data['hits']) == 0:
           return result
        else:
            struct_hit={}
            for hit in data['hits']:
                if hit["hit_uni_os"]== "Homo sapiens":
                    struct_hit=hit 
                    break 
            
            hit_id=struct_hit["hit_id"]    
            hit_def=struct_hit["hit_def"] 
            hit_acc=struct_hit["hit_acc"].split("-")[0] #Uniprot Accession Id
            hit_uni_de=struct_hit["hit_uni_de"] 
            hit_uni_os=struct_hit["hit_uni_os"] 
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
        
def find_Sequence_Variants(query,subject):
    variants=[]
    if(query!="" and subject!=""):
    
        qry = SeqRecord(Seq(query), id="query")
        sbj = SeqRecord(Seq(subject), id="refseq") 
           
        for i in range(len(qry.seq)):
            if qry.seq[i]!=sbj.seq[i]:
                m = {'original': sbj.seq[i], 
                    'variation': qry.seq[i],
                    'position':str(i+1)}
                variants.append(m)
    return variants

def read_Uniprot_Json(json_file,variants):
    try:
        with open(json_file) as f:
            data = json.load(f)

        struct_evidences_id=""
        protein_function=""
        catalytic_activity=""
        disease=""
        acronym=""
        disease_description=""
        struct_diaseases=[]

        for feat in data["features"]:
            if(feat['type']=="Natural variant" or feat['type']=="Mutagenesis"):
                for v in variants:
                    if(feat["location"]["start"]["value"]==int(v['position'])
                    and feat['alternativeSequence']['originalSequence']==v["original"]):
                        for aseq in feat['alternativeSequence']['alternativeSequences']:
                            if(aseq==v["variation"]):
                                struct_evidences_id=feat["evidences"][0]["id"]


        for c in data["comments"]:
            if(c['commentType']=="FUNCTION"):
                if ("texts" in c): 
                    protein_function=c["texts"][0]["value"]

            if(c['commentType']=="CATALYTIC ACTIVITY"):
                if ("reaction" in c): 
                    catalytic_activity=c["reaction"]["name"]

            if(c['commentType']=="DISEASE"):
                if ("disease" in c):
                    for ev in c['disease']['evidences']:
                        if(ev['id']==struct_evidences_id):
                            disease=c['disease']['diseaseId']
                            acronym=c['disease']['acronym']
                            disease_description=c['disease']['description']

        d = {'protein_function':protein_function,
             'catalytic_activity':catalytic_activity,
             'disease':disease if disease!="" else "No disease found",
             'acronym':acronym if acronym!="" else "N/A" ,
             'description':disease_description if disease_description!="" else "N/A" }
        struct_diaseases.append(d)

        return struct_diaseases
    except:
        raise ValueError(
            f"The UniProt was unabled to return info about yor protein.")
#endregion