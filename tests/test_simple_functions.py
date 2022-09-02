import pytest

import emopec


def test_simple_expression_estimation():
    exp = emopec.get_expression('AGGAGA')
    assert pytest.approx(exp) == 9.41932973140501

def test_predict_spacing():
    my_leader = 'ATACAAGTCGCTTAAGGCTTGCCAACGAACCATTGCCGCC'
    upstream, sd_seq, spacing, expr = emopec.predict_spacing(my_leader)
    expr_2 = emopec.get_expression(sd_seq, sd_dist=len(spacing))
    assert upstream == 'ATACAAGTCGCTTAAGGCTTGCCAAC'
    assert sd_seq == 'GAACCA'
    assert spacing == 'TTGCCGCC'
    assert pytest.approx(expr) == 1.6219460336049762
    assert pytest.approx(expr_2) == 1.6219460336049762

def test_library_generation():
    my_cds = 'ATGAAGTTTATCATTAAATTGTTCCCGGAAATCACCATCAAAAGCC'
    upstream = 'ATACAAGTCGCTTAAGGCTTGCCAAC'
    sd_seq = 'GAACCA'
    spacing = 'TTGCCGCC'
    lib = emopec.make_library(upstream, sd_seq, spacing, my_cds)
    lib_dict = {vals[0]: tuple(vals[1:]) for vals in lib}
    control_dict = {
        "AGACAT": (2.6361, 0.30, 40),
        "TAGAGT": (3.6126, 0.70, 55),
        "CGAGTG": (4.6015, 2.50, 70),
        "AGAGTA": (5.6119, 0.00, 86),
        "AGGAGT": (6.3893, 1.20, 98),
        }
    for seq, ref_scores in control_dict.items():
        assert seq in lib_dict
        assert pytest.approx(lib_dict[seq][0], 0.0001) == ref_scores[0]
        assert pytest.approx(lib_dict[seq][1], 0.01) == ref_scores[1]
        assert pytest.approx(100*lib_dict[seq][2], 1) == ref_scores[2]
