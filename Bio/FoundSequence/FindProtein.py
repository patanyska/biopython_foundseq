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
def read_Json(json_file,variants):
    try:
        breaker=False
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
                    if(feat["location"]["start"]["value"]==int(v['position']) and feat['alternativeSequence']['originalSequence']==v["original"]):
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
                 'disease': disease if disease!="" else "NF",
                 'acronym':acronym if acronym!="" else "NF" ,
                 'description':disease_description if disease_description!="" else "NF" }
            struct_diaseases.append(d)
                                    
            return struct_diaseases
    except:
           raise ValueError(
            f"A problem occurred during reading UniProt json"
        )