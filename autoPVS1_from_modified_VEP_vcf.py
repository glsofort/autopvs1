#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 19:13

import os
import re
import sys
import random
import string
from collections import namedtuple

from pvs1 import PVS1
from cnv import PVS1CNV, CNVRecord
from read_data import trans_gene, gene_trans, gene_alias, vep_cache
from read_data import transcripts_hg19, transcripts_hg38, genome_hg19, genome_hg38
from utils import vep2vcfv2, vep2vcf, get_transcript, vep_consequence_trans, VCFRecord


lof_type = ['frameshift', 'nonsense', 'splice-5', 'splice-3', 'init-loss']
vep_lof_list = ['frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 'start_lost']
VAR = namedtuple('VAR', ('varid', 'gene', 'trans', 'canonical', 'pick', 'record'))


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def main():
    genome_version = sys.argv[1]
    anno_fh = open(sys.argv[2])
    header = list()
    for line in anno_fh:
        if line.strip().startswith("#"):
            header = line.strip().split("\t")
            header[0] = header[0].replace("#", "")
            continue
        records = line.strip().split("\t")
        if len(records) == len(header):
            info = dict(zip(header, records))
        else:
            print(f"lenrecords: {len(records)}; lenheader: {len(header)}")
            print(f"records: {records}; header: {header}")
            raise Exception("Inconsistent length for line and header!")

        if genome_version == 'hg19':
            vcfrecord = vep2vcfv2(info['CHROM'], info['POS'], info['REF'], info['ALT'], genome_hg19)
            transcript = get_transcript(info['Feature'], transcripts_hg19)
        else:
            vcfrecord = vep2vcfv2(info['CHROM'], info['POS'], info['REF'], info['ALT'], genome_hg38)
            transcript = get_transcript(info['Feature'], transcripts_hg38)
        
        consequence = vep_consequence_trans(info['Consequence'])
        vcf_id = "-".join([vcfrecord.chrom, str(vcfrecord.pos), vcfrecord.ref, vcfrecord.alt])
        print(transcript)
        
        if consequence in lof_type and transcript:
            lof_pvs1 = PVS1(vcfrecord, consequence, info['HGVSc'], info['HGVSp'], transcript, genome_version)
            trans_name = lof_pvs1.transcript.full_name

            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  trans_name,
                  lof_pvs1.consequence,
                  lof_pvs1.strength_raw.name,
                  lof_pvs1.strength.name,
                  lof_pvs1.criterion,
                  sep="\t")
        elif consequence in lof_type:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_canonical',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")
        else:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_lof',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")

        '''
        if info['Function'] in ['splice-5', 'splice-3']:
            pvs1 = PVS1(vcfrecord, info['Function'], info['pHGVS1'], info['MIMInheritance'], transcript)
            splice = Splicing(vcfrecord, transcript)
            print(vcf_id,
                  pvs1.function,
                  pvs1.strength_raw.name,
                  pvs1.strength.name,
                  splice.is_undergo_NMD,
                  splice.has_cryptic_splice_site,
                  splice.is_exon_skipping,
                  splice.preserves_reading_frame,
                  splice.is_critical_to_protein_func,
                  splice.variant_removes_10_percent_of_protein,
                  splice.is_critical_to_protein_func_detail,
                  sep="\t")
        '''


if __name__ == '__main__':
    main()
