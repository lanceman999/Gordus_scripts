#!/usr/bin/env python
import os

pileup_file = 'MaleProsoma_scaffold_5_.mpileup'

base_name = os.path.splitext(os.path.basename(pileup_file))[0]

scaffold_number = "_".join(base_name.split("_")[1:3])

print(str(scaffold_number))