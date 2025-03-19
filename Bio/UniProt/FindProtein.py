# Copyright 2024 by Patricia Nogueira.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Code to invoke the UNIPROT API over the internet.

This module provides code to work with the NCBI BLAST+ provided by EMBL-EBI
Variables:

    -hit_accession  string	BLAST unique code to identify the protein page.
"""


import os
import requests
import json


UNIPROT_URL="https://rest.uniprot.org/uniprotkb/"

def uniprot(hit_accession):
    if(len(hit_accession)>0):
        
        request = UNIPROT_URL+hit_accession+".json"
        browser = requests.Session()
        response = browser.post(request,headers={'Accept': 'application/json'})
        return response.content.decode('utf8')
    else:
        raise ValueError(
            f"Expected hit_accession. The value cannot be empty."
        )

        