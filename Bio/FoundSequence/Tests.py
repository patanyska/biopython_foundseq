# Copyright 2024 by Patricia Nogueira.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


'''Code to test the several tool developed in new Biopython library'''

import unittest
import TranslateTool
import BlastTool
import FindProtein
import DrugBankTool
import __init__

from pathlib import Path
root_folder = Path(__file__).parents[2]


file_path_empty = str(root_folder) +'\\Tests\FoundSequence\\empty.fasta'
file_path_error = str(root_folder) +'\\Tests\FoundSequence\\error.fasta'
file_path_malformated = str(root_folder) +'\\Tests\FoundSequence\\malformated.fasta'
file_path1 = str(root_folder) +'\\Tests\FoundSequence\\mutseq1.fasta'
file_path8 = str(root_folder) +'\\Tests\FoundSequence\\mutseq8.fasta'
file_path11 = str(root_folder) +'\\Tests\FoundSequence\\mustseq11.fasta'
file_path15 = str(root_folder) +'\\Tests\FoundSequence\\mutseq15.fasta'

class TestTranslateTool(unittest.TestResult):
    def testValidFiles(self):
   
       (TranslateTool.expasy_Translate_Tool(file_path15))


