"""Gillespie Direct Method による二元共重合 kMC シミュレータ。

対象反応(開始・伝播・停止のみ、連鎖移動なし。端末モデル):
    I        --kd,f-->  2 R・                (開始剤分解、効率fでケージ効果を考慮)
    R・ + M1 --ki1-->    P1・                 (連鎖開始)
    R・ + M2 --ki2-->    P2・
    P1・+ M1 --k11-->    P1・                 (伝播、末端モノマーで速度定数が決まる)
    P1・+ M2 --k12-->    P2・
    P2・+ M1 --k21-->    P1・
    P2・+ M2 --k22-->    P2・
    P1・+ P1・ --ktc11--> D(結合)             (停止)
    P1・+ P1・ --ktd11--> D + D(不均化)
    P2・+ P2・ --ktc22--> D(結合)
    P2・+ P2・ --ktd22--> D + D(不均化)
    P1・+ P2・ --ktc12--> D(結合、交差)
    P1・+ P2・ --ktd12--> D + D(不均化、交差)

単一モノマー(ホモポリマー)は n_monomer2_0=0 とすれば、M2側の反応は
propensity=0 となり自然に発火しなくなる(RateConstants.homopolymer 参照)。
"""
from __future__ import annotations

import random
from dataclasses import dataclass

from .kinetics import RateConstants, compute_propensities
from .state import ReactorState


@dataclass
class TrajectoryPoint:
    time: float
    conversion: float
    n_p1: int
    n_p2: int


@dataclass
class SimulationResult:
    state: ReactorState
    trajectory: list
    n_monomer1_initial: int
    n_monomer2_initial: int
    stopped_reason: str


def run_kmc(
    rc: RateConstants,
    volume_L: float,
    n_initiator0: int,
    n_monomer1_0: int,
    n_monomer2_0: int = 0,
    max_time: float = float("inf"),
    max_events: int = 10_000_000,
    record_every: int = 1000,
    seed: int | None = None,
) -> SimulationResult:
    """kMC シミュレーションを実行する。

    Parameters
    ----------
    volume_L:
        シミュレーション体積 [L]。分子数と濃度(mol/L)を結びつける。
        初期モノマー濃度 [M]0(合計)にしたい場合は
        volume_L = (n_monomer1_0 + n_monomer2_0) / (AVOGADRO * [M]0) として求める。
    n_monomer2_0:
        モノマー2の初期分子数。0ならホモポリマー(単一モノマー)として振る舞う。
    max_events:
        無限ループ防止のための安全上限。
    record_every:
        何イベントごとに軌跡(転化率など)を記録するか。
    """
    if n_initiator0 <= 0:
        raise ValueError("n_initiator0 は正である必要があります")
    if n_monomer1_0 <= 0:
        raise ValueError("n_monomer1_0 は正である必要があります")
    if n_monomer2_0 < 0:
        raise ValueError("n_monomer2_0 は0以上である必要があります")

    rng = random.Random(seed)
    state = ReactorState(
        n_initiator=n_initiator0, n_radical=0,
        n_monomer1=n_monomer1_0, n_monomer2=n_monomer2_0,
    )
    n_total0 = n_monomer1_0 + n_monomer2_0
    trajectory: list[TrajectoryPoint] = [
        TrajectoryPoint(time=0.0, conversion=0.0, n_p1=0, n_p2=0)
    ]

    stopped_reason = "max_events"
    for event_idx in range(1, max_events + 1):
        prop = compute_propensities(
            rc, volume_L,
            n_initiator=state.n_initiator,
            n_radical=state.n_radical,
            n_monomer1=state.n_monomer1,
            n_monomer2=state.n_monomer2,
            n_p1=state.n_p1,
            n_p2=state.n_p2,
        )
        a_total = prop.total
        if a_total <= 0.0:
            stopped_reason = "no_reactions_left"
            break

        tau = rng.expovariate(a_total)
        if state.time + tau > max_time:
            stopped_reason = "max_time"
            break
        state.time += tau

        # 累積和によるルーレット選択(フラットなif/elifで、ループ内でのタプル/
        # クロージャ生成を避けて性能を確保する)
        r = rng.random() * a_total
        c = prop.decomposition
        if r < c:
            state.do_decomposition(rng, rc.f)
        elif r < (c := c + prop.chain_initiation_1):
            state.do_chain_initiation_1()
        elif r < (c := c + prop.chain_initiation_2):
            state.do_chain_initiation_2()
        elif r < (c := c + prop.prop_11):
            state.do_prop_11(rng)
        elif r < (c := c + prop.prop_12):
            state.do_prop_12(rng)
        elif r < (c := c + prop.prop_21):
            state.do_prop_21(rng)
        elif r < (c := c + prop.prop_22):
            state.do_prop_22(rng)
        elif r < (c := c + prop.term_comb_11):
            state.do_term_comb(rng, state.live_1, state.live_1)
        elif r < (c := c + prop.term_disp_11):
            state.do_term_disp(rng, state.live_1, state.live_1)
        elif r < (c := c + prop.term_comb_22):
            state.do_term_comb(rng, state.live_2, state.live_2)
        elif r < (c := c + prop.term_disp_22):
            state.do_term_disp(rng, state.live_2, state.live_2)
        elif r < (c := c + prop.term_comb_12):
            state.do_term_comb(rng, state.live_1, state.live_2)
        else:
            state.do_term_disp(rng, state.live_1, state.live_2)

        if event_idx % record_every == 0:
            n_remaining = state.n_monomer1 + state.n_monomer2
            conversion = 1.0 - n_remaining / n_total0
            trajectory.append(
                TrajectoryPoint(
                    time=state.time, conversion=conversion,
                    n_p1=state.n_p1, n_p2=state.n_p2,
                )
            )

        if state.n_monomer1 + state.n_monomer2 == 0:
            stopped_reason = "monomer_depleted"
            break

    n_remaining = state.n_monomer1 + state.n_monomer2
    trajectory.append(
        TrajectoryPoint(
            time=state.time,
            conversion=1.0 - n_remaining / n_total0,
            n_p1=state.n_p1,
            n_p2=state.n_p2,
        )
    )
    return SimulationResult(
        state=state,
        trajectory=trajectory,
        n_monomer1_initial=n_monomer1_0,
        n_monomer2_initial=n_monomer2_0,
        stopped_reason=stopped_reason,
    )
