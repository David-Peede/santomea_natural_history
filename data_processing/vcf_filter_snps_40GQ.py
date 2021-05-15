#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 14:04:25 2021

@author: davidpeede
"""

import gzip
import sys

def gq_filter(vcf):
    """
    vcf: VCF filtered for DP and MQ
    --------------------------------------------------------------------------
    output: A VCF filtered for SNPs with a GQ of 40 or higher
    """
    outfile = sys.stdout
    with gzip.open(vcf, "rt") as data:
        for line in data:
            if line.startswith("##") or line.startswith("#"):
                outfile.write(line)
            else:
                spline = line.split()
                FORMAT = spline[8]
                geno = spline[9]
                ref_len = len(spline[3])
                alt_len = len(spline[4])
                if (((ref_len + alt_len) > 2)):
                    continue
                else:
                    if (len(FORMAT) == 2):
                        outfile.write(line)
                    else:
                        GT, PL, GQ = geno.split(':')
                        if (int(GQ) < 40):
                            continue
                        else:
                            outfile.write(line)
    return

gq_filter(str(sys.argv[1]))