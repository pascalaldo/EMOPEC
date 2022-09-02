import pytest

import emopec

def test_library_generation():
    my_cds = 'ATGAAGTTTATCATTAAATTGTTCCCGGAAATCACCATCAAAAGCC'
    upstream = 'ATACAAGTCGCTTAAGGCTTGCCAAC'
    sd_seq = 'GAACCA'
    spacing = 'TTGCCGCC'

    sd_constraint = 'GWWCCANNNNNNNN'

    lib = emopec.make_library(upstream, sd_seq, spacing, my_cds, dg_seq=sd_constraint)
    lib_dict = {vals[0]: tuple(vals[1:]) for vals in lib}
    for seq in lib_dict:
        print(seq)
        assert (seq[1] in 'AT' and seq[2] in 'AT')

def test_library_generation_multi_constraints():
    my_cds = 'ATGAAGTTTATCATTAAATTGTTCCCGGAAATCACCATCAAAAGCC'
    upstream = 'ATACAAGTCGCTTAAGGCTTGCCAAC'
    sd_seq = 'GAACCA'
    spacing = 'TTGCCGCC'

    sd_constraint = ['GWWCCANNNNNNNN', 'GSSCCANNNNNNNN']

    lib = emopec.make_library(upstream, sd_seq, spacing, my_cds, dg_seq=sd_constraint)
    lib_dict = {vals[0]: tuple(vals[1:]) for vals in lib}
    for seq in lib_dict:
        print(seq)
        assert ((seq[1] in 'AT' and seq[2] in 'AT') or (seq[1] in 'CG' and seq[2] in 'CG'))