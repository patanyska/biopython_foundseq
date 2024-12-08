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
 
