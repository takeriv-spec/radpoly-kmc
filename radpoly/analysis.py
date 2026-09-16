"""シミュレーション結果からの物性量計算(Mn, Mw, PDI, 共重合体組成)と理論値との比較。

理論式の出典:
- Odian, G. "Principles of Polymerization" 4th ed., Wiley (2004), Ch.3, 6.
- Mayo, F.R.; Lewis, F.M. "Copolymerization. I. A Basis for Comparing the Behavior of
  Monomers in Copolymerization..." J. Am. Chem. Soc. 1944, 66, 1594.

方針: 検証できる理論式(単一機構のPDI極限、Mayo-Lewis組成式)のみ実装する。
結合・不均化混在時のPDIや、共重合系のRp(Walling式など交差項φを要する式)は
裏取りが不十分なため実装しない(analysis.theoretical_pdi_limit, theoretical_rp の
docstring参照)。
"""
from __future__ import annotations

from dataclasses import dataclass

from .kinetics import RateConstants
from .simulator import SimulationResult


@dataclass
class ChainLengthStats:
    n_chains: int
    mn: float  # 数平均重合度 (Xn)
    mw: float  # 重量平均重合度 (Xw)
    pdi: float  # 多分散度 Xw/Xn


def chain_length_distribution(result: SimulationResult, include_live: bool = True) -> list[int]:
    """全鎖の全長(M1+M2のモノマー単位数)のリストを返す。"""
    lengths = [n1 + n2 for n1, n2 in result.state.dead_chains]
    if include_live:
        lengths += [c[0] + c[1] for c in result.state.live_1]
        lengths += [c[0] + c[1] for c in result.state.live_2]
    return lengths


def compute_stats(lengths: list[int]) -> ChainLengthStats:
    n = len(lengths)
    if n == 0:
        return ChainLengthStats(n_chains=0, mn=0.0, mw=0.0, pdi=0.0)
    total = sum(lengths)
    mn = total / n
    mw = sum(l * l for l in lengths) / total if total > 0 else 0.0
    pdi = mw / mn if mn > 0 else 0.0
    return ChainLengthStats(n_chains=n, mn=mn, mw=mw, pdi=pdi)


def copolymer_composition(result: SimulationResult, include_live: bool = True) -> tuple[int, int]:
    """生成した全鎖に取り込まれたモノマー1・モノマー2の合計単位数 (n1, n2) を返す。"""
    n1 = sum(d[0] for d in result.state.dead_chains)
    n2 = sum(d[1] for d in result.state.dead_chains)
    if include_live:
        n1 += sum(c[0] for c in result.state.live_1) + sum(c[0] for c in result.state.live_2)
        n2 += sum(c[1] for c in result.state.live_1) + sum(c[1] for c in result.state.live_2)
    return n1, n2


def reactivity_ratios(rc: RateConstants) -> tuple[float, float]:
    """端末モデルの反応性比 r1=k11/k12, r2=k22/k21 を返す。"""
    if rc.k12 <= 0 or rc.k21 <= 0:
        raise ValueError("k12, k21 が0だと反応性比を定義できません")
    return rc.k11 / rc.k12, rc.k22 / rc.k21


def theoretical_instantaneous_composition(rc: RateConstants, f1: float) -> float:
    """Mayo-Lewis式による瞬間共重合体組成 F1(生成する共重合体中のモノマー1モル分率)。

    F1 = (r1*f1^2 + f1*f2) / (r1*f1^2 + 2*f1*f2 + r2*f2^2)
    (Mayo & Lewis, J. Am. Chem. Soc. 1944, 66, 1594)

    f1: フィード(未反応モノマー)中のモノマー1のモル分率(0-1)。転化率が低く
    フィード組成がほぼ一定とみなせる瞬間(重合初期)にのみ有効な近似であり、
    転化が進んだ系全体の累積組成とは一致しない。
    """
    if not (0.0 < f1 < 1.0):
        raise ValueError("f1 は 0 と 1 の間で指定してください")
    f2 = 1.0 - f1
    r1, r2 = reactivity_ratios(rc)
    numerator = r1 * f1 * f1 + f1 * f2
    denominator = r1 * f1 * f1 + 2 * f1 * f2 + r2 * f2 * f2
    return numerator / denominator


def theoretical_rp(rc: RateConstants, conc_I: float, conc_M: float) -> float:
    """定常状態近似による理論的な重合速度 Rp [mol/(L・s)](ホモポリマーのみ)。

    この kMC モデル自身の反応定義(kinetics.py の propensity 式)から導出する
    (外部文献の式番号を記憶だけで引用することは避け、自己完結した導出を示す):

    - ラジカル生成速度(開始): Ri = 2*f*kd*[I]
      (kinetics.py: 分解イベントは kd*[I] の速度で起こり、確率fで2本のラジカルを生む)
    - ラジカル消費速度(停止): Rt = 2*kt*[P・]^2  (kt = ktc11+ktd11)
      (kinetics.py の定義上、停止イベント速度は kt*n*(n-1)/(NA V) であり、
       1イベントにつき2本のラジカルが消費されるため Rt = 2*kt*[P・]^2 になる)
    - 疑似定常状態近似 Ri = Rt より:
        2*f*kd*[I] = 2*kt*[P・]_ss^2
        [P・]_ss = sqrt(f*kd*[I] / kt)
    - Rp = kp*[M]*[P・]_ss

    共重合系(k12>0 または k21>0)でのRpは交差項φを含むWalling式が必要だが、
    その式の裏取りができていないため未実装。ホモポリマー退化ケース
    (RateConstants.homopolymer で構築、またはk12=k21=0)でのみ使用可能。
    """
    if rc.k12 > 0.0 or rc.k21 > 0.0:
        raise NotImplementedError(
            "共重合系のRp理論式(Walling式)は裏取りが不十分なため未実装です。"
            "ホモポリマー退化ケース(k12=k21=0)でのみ使用してください。"
        )
    kt = rc.ktc11 + rc.ktd11
    if kt <= 0:
        raise ValueError("kt(=ktc11+ktd11) が 0 のため定常状態近似が使えません")
    p_ss = (rc.f * rc.kd * conc_I / kt) ** 0.5
    return rc.k11 * conc_M * p_ss


def theoretical_pdi_limit(rc: RateConstants) -> float:
    """停止機構のみで決まる瞬間PDIの理論極限(ホモポリマーのみ、連鎖移動なし、鎖長が十分大きい極限)。

    厳密解が確立しているのは純粋な機構のみ:
    - 結合停止のみ (ktd11=0): PDI -> 1.5
    - 不均化停止のみ (ktc11=0): PDI -> 2.0
    (Odian 4th ed., 3-4節)

    共重合系、および結合・不均化が混在する場合の閉じた式は未検証のため実装していない。
    """
    if rc.k12 > 0.0 or rc.k21 > 0.0:
        raise NotImplementedError(
            "共重合系のPDI理論式は未実装です。ホモポリマー退化ケースでのみ使用してください。"
        )
    if rc.ktd11 == 0.0 and rc.ktc11 > 0.0:
        return 1.5
    if rc.ktc11 == 0.0 and rc.ktd11 > 0.0:
        return 2.0
    raise NotImplementedError(
        "結合・不均化が混在する場合のPDI理論式は未検証のため実装していません。"
        "ktc11 または ktd11 のどちらかを0にして純粋な機構で検証してください。"
    )
