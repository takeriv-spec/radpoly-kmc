import pytest

from radpoly.analysis import compute_stats, reactivity_ratios, theoretical_instantaneous_composition, theoretical_pdi_limit, theoretical_rp
from radpoly.kinetics import RateConstants


def test_compute_stats_uniform_distribution():
    stats = compute_stats([10, 10, 10, 10])
    assert stats.n_chains == 4
    assert stats.mn == pytest.approx(10.0)
    assert stats.mw == pytest.approx(10.0)
    assert stats.pdi == pytest.approx(1.0)


def test_compute_stats_manual_calc():
    lengths = [1, 2, 3]
    stats = compute_stats(lengths)
    assert stats.mn == pytest.approx(2.0)
    # Mw = sum(l^2)/sum(l) = (1+4+9)/6 = 14/6
    assert stats.mw == pytest.approx(14 / 6)
    assert stats.pdi == pytest.approx((14 / 6) / 2.0)


def test_compute_stats_empty():
    stats = compute_stats([])
    assert stats.n_chains == 0
    assert stats.pdi == 0.0


def test_theoretical_pdi_limit_combination_only():
    rc = RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=1e7, ktd=0.0)
    assert theoretical_pdi_limit(rc) == pytest.approx(1.5)


def test_theoretical_pdi_limit_disproportionation_only():
    rc = RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=0.0, ktd=1e7)
    assert theoretical_pdi_limit(rc) == pytest.approx(2.0)


def test_theoretical_pdi_limit_mixed_not_implemented():
    rc = RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=1e7, ktd=1e7)
    with pytest.raises(NotImplementedError):
        theoretical_pdi_limit(rc)


def test_theoretical_pdi_limit_copolymer_not_implemented():
    rc = RateConstants(
        kd=1e-5, f=0.6, ki1=1e3, ki2=1e3, k11=1e3, k12=5e2, k21=5e2, k22=1e3,
        ktc11=1e7, ktc12=1e7, ktc22=1e7, ktd11=0, ktd12=0, ktd22=0,
    )
    with pytest.raises(NotImplementedError):
        theoretical_pdi_limit(rc)


def test_theoretical_rp_scales_with_sqrt_kd():
    rc1 = RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=1e7, ktd=0.0)
    rc2 = RateConstants.homopolymer(kd=4e-5, f=0.6, ki=1e3, kp=1e3, ktc=1e7, ktd=0.0)
    rp1 = theoretical_rp(rc1, conc_I=0.01, conc_M=1.0)
    rp2 = theoretical_rp(rc2, conc_I=0.01, conc_M=1.0)
    # kd が4倍なら sqrt(kd) は2倍 -> Rpも2倍
    assert rp2 / rp1 == pytest.approx(2.0, rel=1e-6)


def test_theoretical_rp_copolymer_not_implemented():
    rc = RateConstants(
        kd=1e-5, f=0.6, ki1=1e3, ki2=1e3, k11=1e3, k12=5e2, k21=5e2, k22=1e3,
        ktc11=1e7, ktc12=1e7, ktc22=1e7, ktd11=0, ktd12=0, ktd22=0,
    )
    with pytest.raises(NotImplementedError):
        theoretical_rp(rc, conc_I=0.01, conc_M=1.0)


def test_reactivity_ratios():
    rc = RateConstants(
        kd=0, f=1.0, ki1=1, ki2=1, k11=10.0, k12=2.0, k21=3.0, k22=6.0,
        ktc11=1, ktc12=1, ktc22=1, ktd11=0, ktd12=0, ktd22=0,
    )
    r1, r2 = reactivity_ratios(rc)
    assert r1 == pytest.approx(5.0)
    assert r2 == pytest.approx(2.0)


def test_mayo_lewis_ideal_copolymerization():
    """r1*r2=1 の理想共重合の場合、F1 = r1*f1 / (r1*f1 + f2) という単純式に一致するはず。"""
    rc = RateConstants(
        kd=0, f=1.0, ki1=1, ki2=1, k11=4.0, k12=2.0, k21=2.0, k22=1.0,  # r1=2, r2=0.5 -> r1*r2=1
        ktc11=1, ktc12=1, ktc22=1, ktd11=0, ktd12=0, ktd22=0,
    )
    f1 = 0.3
    F1 = theoretical_instantaneous_composition(rc, f1)
    r1 = 2.0
    expected = r1 * f1 / (r1 * f1 + (1 - f1))
    assert F1 == pytest.approx(expected)


def test_mayo_lewis_azeotrope_when_r1_r2_equal_and_less_than_one():
    """r1=r2の場合、f1=0.5でF1=0.5(共沸点、フィードと組成が一致)になるはず。"""
    rc = RateConstants(
        kd=0, f=1.0, ki1=1, ki2=1, k11=1.0, k12=2.0, k21=2.0, k22=1.0,  # r1=r2=0.5
        ktc11=1, ktc12=1, ktc22=1, ktd11=0, ktd12=0, ktd22=0,
    )
    F1 = theoretical_instantaneous_composition(rc, 0.5)
    assert F1 == pytest.approx(0.5)
