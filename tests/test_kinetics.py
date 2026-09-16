import math

import pytest

from radpoly.kinetics import AVOGADRO, RateConstants, compute_propensities, geometric_mean_cross_termination


def test_rate_constants_rejects_negative():
    with pytest.raises(ValueError):
        RateConstants(
            kd=-1, f=0.5, ki1=1, ki2=1, k11=1, k12=1, k21=1, k22=1,
            ktc11=1, ktc12=1, ktc22=1, ktd11=1, ktd12=1, ktd22=1,
        )


def test_rate_constants_rejects_bad_efficiency():
    with pytest.raises(ValueError):
        RateConstants(
            kd=1, f=1.5, ki1=1, ki2=1, k11=1, k12=1, k21=1, k22=1,
            ktc11=1, ktc12=1, ktc22=1, ktd11=1, ktd12=1, ktd22=1,
        )


def test_homopolymer_convenience_constructor():
    rc = RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=1e7, ktd=0.0)
    assert rc.k11 == 1e3
    assert rc.k12 == 0.0
    assert rc.k21 == 0.0
    assert rc.k22 == 0.0
    assert rc.ktc11 == 1e7
    assert rc.ktc22 == 0.0


def test_geometric_mean_cross_termination():
    assert geometric_mean_cross_termination(4.0, 9.0) == pytest.approx(6.0)


def test_compute_propensities_unimolecular_decomposition():
    rc = RateConstants.homopolymer(kd=2.0, f=1.0, ki=0.0, kp=0.0, ktc=1.0, ktd=0.0)
    prop = compute_propensities(
        rc, volume_L=1.0, n_initiator=100, n_radical=0, n_monomer1=0, n_monomer2=0, n_p1=0, n_p2=0
    )
    assert prop.decomposition == pytest.approx(200.0)
    assert prop.chain_initiation_1 == 0.0
    assert prop.prop_11 == 0.0


def test_compute_propensities_bimolecular_hand_check():
    # 手計算での検算: ki1 * nR * nM1 / (NA*V)
    rc = RateConstants.homopolymer(kd=0.0, f=1.0, ki=3.0, kp=0.0, ktc=1.0, ktd=0.0)
    volume_L = 1e-20  # NA*V を手頃な値にするため小さい体積を使う
    na_v = AVOGADRO * volume_L
    prop = compute_propensities(
        rc, volume_L, n_initiator=0, n_radical=10, n_monomer1=5, n_monomer2=0, n_p1=0, n_p2=0
    )
    expected = 3.0 * 10 * 5 / na_v
    assert prop.chain_initiation_1 == pytest.approx(expected)


def test_compute_propensities_same_species_termination_hand_check():
    rc = RateConstants.homopolymer(kd=0.0, f=1.0, ki=0.0, kp=0.0, ktc=5.0, ktd=0.0)
    volume_L = 1e-20
    na_v = AVOGADRO * volume_L
    prop = compute_propensities(
        rc, volume_L, n_initiator=0, n_radical=0, n_monomer1=0, n_monomer2=0, n_p1=7, n_p2=0
    )
    expected = 5.0 * 7 * 6 / na_v
    assert prop.term_comb_11 == pytest.approx(expected)
