"""シミュレータの結合テスト。質量収支と、既知の理論値との照合(検証用)。"""
import random

import pytest

from radpoly.analysis import (
    chain_length_distribution,
    compute_stats,
    copolymer_composition,
    theoretical_instantaneous_composition,
    theoretical_pdi_limit,
)
from radpoly.kinetics import AVOGADRO, RateConstants, compute_propensities
from radpoly.simulator import run_kmc
from radpoly.state import ReactorState


def _typical_homopolymer_rc(ktc=0.0, ktd=0.0):
    # AIBN/スチレン様の典型的なオーダーの速度定数(検証用、実測値そのものではない)
    return RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=ktc, ktd=ktd)


def test_homopolymer_mass_balance_holds_throughout():
    rc = _typical_homopolymer_rc(ktc=1e8)
    n_monomer0 = 20_000
    volume_L = n_monomer0 / (AVOGADRO * 5.0)
    result = run_kmc(rc, volume_L, n_initiator0=200, n_monomer1_0=n_monomer0, seed=42)
    result.state.check_mass_balance(n_monomer0, 0)


def test_homopolymer_conversion_is_monotonic_and_bounded():
    rc = _typical_homopolymer_rc(ktc=1e8)
    n_monomer0 = 20_000
    volume_L = n_monomer0 / (AVOGADRO * 5.0)
    result = run_kmc(rc, volume_L, n_initiator0=200, n_monomer1_0=n_monomer0, seed=1)
    conversions = [pt.conversion for pt in result.trajectory]
    assert all(0.0 <= c <= 1.0 for c in conversions)
    assert all(b >= a - 1e-9 for a, b in zip(conversions, conversions[1:]))


def test_homopolymer_reproducible_with_seed():
    rc = _typical_homopolymer_rc(ktc=1e8)
    n_monomer0 = 5_000
    volume_L = n_monomer0 / (AVOGADRO * 5.0)
    r1 = run_kmc(rc, volume_L, n_initiator0=50, n_monomer1_0=n_monomer0, seed=7)
    r2 = run_kmc(rc, volume_L, n_initiator0=50, n_monomer1_0=n_monomer0, seed=7)
    assert r1.trajectory[-1].conversion == r2.trajectory[-1].conversion
    assert r1.state.dead_chains == r2.state.dead_chains


def test_termination_only_matches_analytic_second_order_decay():
    """伝播・開始を止め、純粋な停止反応だけの減衰が d[P]/dt=-2kt[P]^2 の解析解に従うことを確認する。

    解析解: [P](t) = [P]0 / (1 + 2*kt*[P]0*t)
    (状態を直接組み立てて停止反応のみループさせる、シミュレータ全体は使わない)
    """
    ktc = 1e7  # L/(mol・s)
    volume_L = 1e-15
    na_v = AVOGADRO * volume_L
    n_p0 = 2000
    conc_p0 = n_p0 / na_v

    rc = RateConstants.homopolymer(kd=0.0, f=1.0, ki=0.0, kp=0.0, ktc=ktc, ktd=0.0)
    state = ReactorState(n_initiator=0, n_radical=0, n_monomer1=10**9, n_monomer2=0)
    state.live_1 = [[1, 0] for _ in range(n_p0)]

    rng = random.Random(123)
    target_fraction_remaining = 0.5
    # [P]/[P]0 = 1/(1+2*kt*[P]0*t) = target -> t = (1/target - 1)/(2*kt*[P]0)
    t_target = (1.0 / target_fraction_remaining - 1.0) / (2 * ktc * conc_p0)

    while state.time < t_target and state.n_p1 > 1:
        prop = compute_propensities(
            rc, volume_L, n_initiator=0, n_radical=0,
            n_monomer1=state.n_monomer1, n_monomer2=0,
            n_p1=state.n_p1, n_p2=0,
        )
        a_total = prop.term_comb_11  # このテストでは唯一有効な反応
        if a_total <= 0:
            break
        tau = rng.expovariate(a_total)
        state.time += tau
        state.do_term_comb(rng, state.live_1, state.live_1)

    analytic_remaining_fraction = 1.0 / (1.0 + 2 * ktc * conc_p0 * state.time)
    simulated_remaining_fraction = state.n_p1 / n_p0

    assert simulated_remaining_fraction == pytest.approx(analytic_remaining_fraction, rel=0.1)


def _run_constant_population_termination(kp, kt, mechanism, n_p_target, n_dead_target, seed):
    """生きているラジカル数 n_p1 を厳密に一定に保ちながら、伝播と単一機構の停止だけを
    競合させる直接検証。停止で消費された分だけ即座に長さ1の鎖を補充することで
    "疑似定常状態" を厳密に(近似ではなく)実現し、理論のPDI極限(1.5/2.0)がkMCの
    コア機構(propagation/termination の propensity 計算とランダム選択)から
    正しく再現されることを検証する。RateConstants.homopolymer や run_kmc の
    開始剤分解・連鎖開始は使わない(それらのタイミング効果を混入させないため)。

    mechanism: "comb" または "disp"
    """
    assert mechanism in ("comb", "disp")
    volume_L = 1.0 / AVOGADRO  # NA*V = 1 とし、分子数をそのまま無次元カウントとして扱う
    n_monomer0 = 10**8  # 十分多く、枯渇しない

    ktc = kt if mechanism == "comb" else 0.0
    ktd = kt if mechanism == "disp" else 0.0
    rc = RateConstants.homopolymer(kd=0.0, f=1.0, ki=0.0, kp=kp, ktc=ktc, ktd=ktd)

    state = ReactorState(n_initiator=0, n_radical=0, n_monomer1=n_monomer0, n_monomer2=0)
    state.live_1 = [[1, 0] for _ in range(n_p_target)]
    rng = random.Random(seed)

    while len(state.dead_chains) < n_dead_target:
        prop = compute_propensities(
            rc, volume_L, n_initiator=0, n_radical=0,
            n_monomer1=state.n_monomer1, n_monomer2=0,
            n_p1=state.n_p1, n_p2=0,
        )
        term_prop = prop.term_comb_11 if mechanism == "comb" else prop.term_disp_11
        a_total = prop.prop_11 + term_prop
        tau = rng.expovariate(a_total)
        state.time += tau
        r = rng.random() * a_total
        if r < prop.prop_11:
            state.do_prop_11(rng)
        else:
            if mechanism == "comb":
                state.do_term_comb(rng, state.live_1, state.live_1)
            else:
                state.do_term_disp(rng, state.live_1, state.live_1)
            # 個体数 n_p1 を一定に保つため、消費した2本分を長さ1の新しい鎖で補充する
            state.live_1.append([1, 0])
            state.live_1.append([1, 0])

    return state


def test_core_termination_mechanism_matches_theoretical_pdi_combination():
    state = _run_constant_population_termination(
        kp=1e-6, kt=0.0025, mechanism="comb", n_p_target=200, n_dead_target=6000, seed=2024,
    )
    lengths = [n1 + n2 for n1, n2 in state.dead_chains]
    stats = compute_stats(lengths)
    assert stats.pdi == pytest.approx(1.5, abs=0.05)


def test_core_termination_mechanism_matches_theoretical_pdi_disproportionation():
    state = _run_constant_population_termination(
        kp=1e-6, kt=0.0025, mechanism="disp", n_p_target=200, n_dead_target=6000, seed=2024,
    )
    lengths = [n1 + n2 for n1, n2 in state.dead_chains]
    stats = compute_stats(lengths)
    # 有限のXn(~100)による理論極限からの補正(PDI_exact=2-1/Xn≈1.99)と
    # 有限サンプル数による統計誤差を許容するトレランス
    assert stats.pdi == pytest.approx(2.0, abs=0.08)


def test_homopolymer_full_pipeline_pdi_is_within_plausible_range():
    """開始剤分解を含むフルパイプラインでは、序盤([M]大)と終盤([M]小)の鎖が混ざり
    バッチ全体の累積PDIは瞬間理論値(1.5/2.0)よりも大きくなる
    (「連鎖長ドリフト」、バッチ重合の既知の性質)。ここでは理論値への厳密な一致では
    なく、崩壊していないこと(妥当な範囲に収まること)だけを確認する
    厳密な機構検証は test_core_termination_mechanism_* 側で行っている。
    """
    rc = RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=1e8, ktd=0.0)
    n_monomer0 = 20_000
    volume_L = n_monomer0 / (AVOGADRO * 5.0)
    result = run_kmc(rc, volume_L, n_initiator0=200, n_monomer1_0=n_monomer0, seed=42)
    lengths = chain_length_distribution(result, include_live=False)
    stats = compute_stats(lengths)
    assert stats.n_chains > 10
    assert 1.0 <= stats.pdi <= 3.0


def test_copolymer_mass_balance_holds():
    rc = RateConstants(
        kd=1e-5, f=0.6, ki1=1e3, ki2=1e3,
        k11=1e3, k12=5e2, k21=8e2, k22=1e3,
        ktc11=1e8, ktc12=1e8, ktc22=1e8, ktd11=0, ktd12=0, ktd22=0,
    )
    n1_0, n2_0 = 10_000, 10_000
    volume_L = (n1_0 + n2_0) / (AVOGADRO * 5.0)
    result = run_kmc(rc, volume_L, n_initiator0=100, n_monomer1_0=n1_0, n_monomer2_0=n2_0, seed=3)
    result.state.check_mass_balance(n1_0, n2_0)


def test_copolymer_composition_matches_mayo_lewis_at_low_conversion():
    """低転化率(フィード組成がほぼ一定)では、生成共重合体の組成がMayo-Lewis式に近いはず。"""
    rc = RateConstants(
        kd=1e-4, f=0.6, ki1=1e3, ki2=1e3,
        k11=1e3, k12=2e2, k21=5e2, k22=1e3,  # r1=5, r2=2
        ktc11=1e9, ktc12=1e9, ktc22=1e9, ktd11=0, ktd12=0, ktd22=0,
    )
    # フィード比を大きく偏らせず、転化率を低く止めて「フィード組成一定」近似を満たす
    f1_feed = 0.5
    n1_0, n2_0 = 40_000, 40_000
    volume_L = (n1_0 + n2_0) / (AVOGADRO * 5.0)
    result = run_kmc(
        rc, volume_L, n_initiator0=20, n_monomer1_0=n1_0, n_monomer2_0=n2_0,
        max_events=2_000_000, seed=5,
    )
    conversion = result.trajectory[-1].conversion
    assert conversion < 0.15  # フィード組成一定近似が妥当な範囲に収まっていることを確認

    n1_made, n2_made = copolymer_composition(result, include_live=False)
    simulated_F1 = n1_made / (n1_made + n2_made)
    theory_F1 = theoretical_instantaneous_composition(rc, f1_feed)
    assert simulated_F1 == pytest.approx(theory_F1, abs=0.05)
