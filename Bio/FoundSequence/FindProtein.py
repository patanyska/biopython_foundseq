# Copyright 2024 by Patricia Nogueira.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Code to invoke the Uniprot over the internet.

This module provides code to work with the NCBI BLAST+ provided by SWISS-PROT
Details about API here https://www.ebi.ac.uk/proteins/api/doc/
Variables:

    - hit accession - key to identify a protein page

"""

import requests
import json



UNIPROT_URL="https://rest.uniprot.org/uniprotkb/"

""" Function to invoke URL from UniProt and get info about a protein
    Variables:

    - accession    key to identify a protein page
"""
def found_Uniprot_Protein(hit_accession):
    try:

        requestURL = UNIPROT_URL+hit_accession+".json"
        browser = requests.Session()
        response = browser.post(
            requestURL,
            headers={'Accept': 'application/json'})
        response.raise_for_status()
        bytes_data = response.content
        json_data = json.loads(bytes_data.decode('utf-8'))
        return json_data
    except:
         raise ValueError(f"A problem occurred finding information about protein in UniProt"
        )


""" Function to read UniProt Json and find disease according to mutations
       Variables:
       
        - json_file          json file return from uniprot API
        - variants          variants found in sequence returned from Blastp        
"""
 
def read_Uniprot_Json(self,json_file,variants):
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
        