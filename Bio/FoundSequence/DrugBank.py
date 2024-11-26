# Copyright 2024 by Patricia Nogueira.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

import DrugBankDataAccess

"""Code to invoke the DrugBank results.

This module provides code to work with drugs provided by DRUG BANK
Variables:

    - diseaase - disease identification or acronym

"""


def found_Drug(diseases):
    if(len(diseases)>0):
        connection=DrugBankDataAccess.openConnectionDrugBank()
        drugs=DrugBankDataAccess.get_drug_by_indication(connection,diseases)
        DrugBankDataAccess.closeConnectionDrugBank(connection)
        return drugs
    else:
        raise ValueError(
            f"The disease field is mandatpry to get drugs")
