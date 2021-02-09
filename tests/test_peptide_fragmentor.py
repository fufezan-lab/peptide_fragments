import pytest
import numpy as np
from pandas import DataFrame

from peptide_fragmentor import PeptideFragment0r


# def test_precursor_ion_calculated_correctly():
#     fragments = PeptideFragment0r('MK', charges=[1])
#     # no series, only precursor
#     fragments = fragger.fragment_peptide(ion_series=[])
#     assert isinstance(fragments, DataFrame)
#     assert len(fragments) == 1
#     row = fragments.iloc[0]



def test_fragment_two_aa_peptide_a_series():
    """Test a1 fragmentation"""
    fragments = PeptideFragment0r('MK', charges=[1], ions=['a']).df
    # fragments = fragger.fragment_peptide(ion_series=['a'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 2
    row = fragments.iloc[0]
    assert row['name'] == 'a1'
    assert row['hill'] == 'C(4)H(9)N(1)S(1)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 104.05284693456998


def test_fragment_three_aa_peptide_a_series():
    """Test a2 fragmentation"""
    fragments = PeptideFragment0r('MKK', charges=[1], ions=['a']).df
    # fragments = fragger.fragment_peptide(ion_series=['a'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 3
    row = fragments.iloc[2]
    assert row['name'] == 'a2'
    assert row['hill'] == 'C(10)H(21)N(3)O(1)S(1)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 232.14780994897004


def test_fragment_two_aa_peptide_b_series():
    """Test b2 fragmentation"""
    fragments = PeptideFragment0r('KK', charges=[1], ions=['b']).df
    # fragments = fragger.fragment_peptide(ion_series=['b'])
    assert isinstance(fragments, DataFrame)

    row = fragments.iloc[6]
    assert row['name'] == 'b2'
    assert row['hill'] == 'C(12)H(24)N(4)O(2)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 257.19720249556997


def test_fragment_three_aa_peptide_b_series():
    """Test b3 fragmentation"""
    fragments = PeptideFragment0r(
        'KKR',
        charges=[1],
        ions=['b'],
        neutral_losses={},
    ).df
    # fragments = fragger.fragment_peptide(ion_series=['b'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 3
    row = fragments.iloc[-1]
    assert row['name'] == 'b3'
    assert row['hill'] == 'C(18)H(36)N(8)O(3)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 413.29831351797  # mass of full peptide - H20 + 1 proton


def test_fragment_three_aa_peptide_c_series():
    """Test c1/2 fragmentation"""
    fragments = PeptideFragment0r('MKK', charges=[1],ions=['c']).df
    # fragments = fragger.fragment_peptide(ion_series=['c'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 3

    row = fragments.iloc[0]
    assert row['name'] == 'c1'
    assert row['hill'] == 'C(5)H(12)N(2)O(1)S(1)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 149.07436

    row = fragments.iloc[2]
    assert row['name'] == 'c2'
    assert row['hill'] == 'C(11)H(24)N(4)O(2)S(1)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 277.16932


def test_fragment_three_aa_peptide_x_series():
    """Test x2 fragmentation"""
    fragments = PeptideFragment0r('MKK', charges=[1], ions=['x']).df
    # fragments = fragger.fragment_peptide(ion_series=['x'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 3

    row = fragments.iloc[3]
    assert row['name'] == 'x2'
    assert row['hill'] == 'C(13)H(24)N(4)O(4)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 301.18708


def test_fragment_one_aa_peptide_y_series():
    """Test y1 fragmentation"""
    fragments = PeptideFragment0r('K', charges=[1], ions=['y']).df
    # fragments = fragger.fragment_peptide(ion_series=['y'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 1

    row = fragments.iloc[1]
    assert row['name'] == 'y1'
    assert row['hill'] == 'C(6)H(14)N(2)O(2)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 147.11280646


def test_fragment_two_aa_peptide_y_series():
    """Test y2 fragmentation"""
    fragments = PeptideFragment0r('KK', charges=[1], ions=['y']).df
    # fragments = fragger.fragment_peptide(ion_series=['y'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 2

    row = fragments.iloc[3]
    assert row['name'] == 'y2'
    assert row['hill'] == 'C(12)H(26)N(4)O(3)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        274.2004907132
    )


def test_fragment_three_aa_peptide_z_series():
    """Test z2 fragmentation"""
    fragments = PeptideFragment0r('MKK', charges=[1], ions=['z']).df
    # fragments = fragger.fragment_peptide(ion_series=['z'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 3

    row = fragments.iloc[3]
    assert row['name'] == 'z2'
    assert row['hill'] == 'C(12)H(24)N(3)O(3)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 259.189


def test_fragment_two_aa_peptide_b_series_with_mod():
    """Test b2 fragmentation with unimod"""
    fragments = PeptideFragment0r('MKK#Oxidation:1', charges=[1], ions=['b']).df
    # fragments = fragger.fragment_peptide(ion_series=['b'])
    assert isinstance(fragments, DataFrame)
    # assert len(fragments) == 3

    row = fragments.iloc[7]
    assert row['name'] == 'b2'
    assert row['hill'] == 'C(11)H(21)N(3)O(3)S(1)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 276.1376


def test_fragment_one_aa_peptide_internal():
    fragments = PeptideFragment0r('RKKR', charges=[1], ions=['I']).df
    # fragments = fragger.fragment_peptide(ion_series=['internal b-y'])
    assert isinstance(fragments, DataFrame)
    print(fragments)


# def test_fragment_two_aa_peptide_a_series_2():
#     """Test a1 fragmentation"""
#     fragger = PeptideFragment0r('KK', charges=[1])
#     fragments = fragger.fragment_peptide(ion_series=['a'])
#     assert isinstance(fragments, DataFrame)
#     # assert len(fragments) == 2
#     row = fragments.iloc[0]
#     assert row['name'] == 'a1'
#     assert row['cc'] == 'C(5)H(12)N(2)'
#     assert row['charge'] == 1
#     assert pytest.approx(
#         row['mz'],
#         5e-6
#     ) == 101.10732486116999


# def test_fragment_three_aa_peptide_a_series_2():
#     """Test a2 fragmentation"""
#     fragger = PeptideFragment0r('KKK', charges=[1])
#     fragments = fragger.fragment_peptide(ion_series=['a'])
#     assert isinstance(fragments, DataFrame)
#     # assert len(fragments) == 3
#     row = fragments.iloc[1]
#     assert row['name'] == 'a2'
#     assert row['cc'] == 'C(11)H(24)N(4)O(1)'
#     assert row['charge'] == 1
#     assert pytest.approx(
#         row['mz'],
#         5e-6
#     ) == 229.20228787557
