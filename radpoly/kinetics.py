"""反応速度定数と Gillespie propensity(遷移確率速度)の計算。

二元共重合(モノマー M1, M2)の「端末モデル」(terminal model / first-order Markov model)
を扱う。反応性は成長鎖の末端モノマーのみで決まると仮定する(前々末端の影響は無視)。
単一モノマー(ホモポリマー)は M2 を使わない(n_monomer2=0)退化ケースとして扱える。

理論的基盤:
- Gillespie, D.T. "Exact Stochastic Simulation of Coupled Chemical Reactions"
  J. Phys. Chem. 81, 2340 (1977). メソスコピック(確率的)速度定数への変換式を採用。
  異種二分子反応 A+B: propensity = k/(NA*V) * nA * nB
  同種二分子反応 A+A: propensity = k * nA*(nA-1) / (NA*V)
    (Gillespie の c=2k/(NA V) と組み合わせ数 n(n-1)/2 の積に一致)
  単分子反応: propensity = k * n
- Odian, G. "Principles of Polymerization" 4th ed., Wiley (2004), Ch.3, 6.
  古典的なフリーラジカル(共)重合の速度論の定義・慣例に準拠。
- Mayo, F.R.; Lewis, F.M. J. Am. Chem. Soc. 1944, 66, 1594.
  端末モデルの共重合体組成式(Mayo-Lewis式)。analysis.py で使用。

単位系:
- kd: 1/s (開始剤分解、単分子反応)
- f: 無次元 (0-1、開始剤効率)
- ki1, ki2, k11, k12, k21, k22, ktc*, ktd*: L/(mol・s) (二分子反応の"化学"速度定数)
  ktcXY, ktdXY は d[P・]/dt|term = -2*kt*[P・]^2 という Odian の慣例に従う定義。
"""
from __future__ import annotations

from dataclasses import dataclass

AVOGADRO = 6.02214076e23  # mol^-1


@dataclass(frozen=True)
class RateConstants:
    kd: float     # 開始剤分解速度定数 [1/s]
    f: float      # 開始剤効率 [-] (0-1)
    ki1: float    # 連鎖開始(R・+M1) [L/(mol・s)]
    ki2: float    # 連鎖開始(R・+M2) [L/(mol・s)]
    k11: float    # 伝播 P1・+M1 -> P1・ [L/(mol・s)]
    k12: float    # 伝播 P1・+M2 -> P2・ [L/(mol・s)]
    k21: float    # 伝播 P2・+M1 -> P1・ [L/(mol・s)]
    k22: float    # 伝播 P2・+M2 -> P2・ [L/(mol・s)]
    ktc11: float  # 停止(結合) P1-P1 [L/(mol・s)]
    ktc12: float  # 停止(結合) P1-P2(交差) [L/(mol・s)]
    ktc22: float  # 停止(結合) P2-P2 [L/(mol・s)]
    ktd11: float  # 停止(不均化) P1-P1 [L/(mol・s)]
    ktd12: float  # 停止(不均化) P1-P2(交差) [L/(mol・s)]
    ktd22: float  # 停止(不均化) P2-P2 [L/(mol・s)]

    def __post_init__(self) -> None:
        for name in (
            "kd", "ki1", "ki2", "k11", "k12", "k21", "k22",
            "ktc11", "ktc12", "ktc22", "ktd11", "ktd12", "ktd22",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} は非負である必要があります")
        if not (0.0 <= self.f <= 1.0):
            raise ValueError("f (開始剤効率) は 0〜1 の範囲で指定してください")

    @classmethod
    def homopolymer(cls, kd: float, f: float, ki: float, kp: float, ktc: float, ktd: float) -> "RateConstants":
        """単一モノマー系(共重合なし)のための便利コンストラクタ。

        モノマー2関連の速度定数はすべて0にする。呼び出し側で n_monomer2_0=0
        として run_kmc を呼べば、モノマー2側の反応は propensity=0 で発火せず、
        実質的に従来のホモポリマーkMCと同じ挙動になる。
        """
        return cls(
            kd=kd, f=f, ki1=ki, ki2=0.0,
            k11=kp, k12=0.0, k21=0.0, k22=0.0,
            ktc11=ktc, ktc12=0.0, ktc22=0.0,
            ktd11=ktd, ktd12=0.0, ktd22=0.0,
        )


def geometric_mean_cross_termination(kt11: float, kt22: float) -> float:
    """交差停止速度定数の幾何平均近似: kt12 = sqrt(kt11 * kt22)。

    化学制御モデル(chemically-controlled termination)でしばしば用いられる簡易則。
    実測の交差停止速度定数が無い場合の近似値として提供するが、拡散律速系
    (ゲル効果が強い系など)では成り立たない場合がある点に注意。
    """
    if kt11 < 0 or kt22 < 0:
        raise ValueError("kt11, kt22 は非負である必要があります")
    return (kt11 * kt22) ** 0.5


@dataclass
class Propensities:
    """各素反応の Gillespie propensity [events/s]。"""

    decomposition: float
    chain_initiation_1: float
    chain_initiation_2: float
    prop_11: float
    prop_12: float
    prop_21: float
    prop_22: float
    term_comb_11: float
    term_disp_11: float
    term_comb_22: float
    term_disp_22: float
    term_comb_12: float
    term_disp_12: float

    @property
    def total(self) -> float:
        return (
            self.decomposition
            + self.chain_initiation_1 + self.chain_initiation_2
            + self.prop_11 + self.prop_12 + self.prop_21 + self.prop_22
            + self.term_comb_11 + self.term_disp_11
            + self.term_comb_22 + self.term_disp_22
            + self.term_comb_12 + self.term_disp_12
        )


def compute_propensities(
    rc: RateConstants,
    volume_L: float,
    n_initiator: int,
    n_radical: int,
    n_monomer1: int,
    n_monomer2: int,
    n_p1: int,
    n_p2: int,
) -> Propensities:
    """現在の分子数から各反応の propensity を計算する。

    n_p1: 末端がM1の生きているラジカル鎖の本数。n_p2: 末端がM2の本数。
    """
    na_v = AVOGADRO * volume_L

    return Propensities(
        decomposition=rc.kd * n_initiator,
        chain_initiation_1=rc.ki1 * n_radical * n_monomer1 / na_v,
        chain_initiation_2=rc.ki2 * n_radical * n_monomer2 / na_v,
        prop_11=rc.k11 * n_p1 * n_monomer1 / na_v,
        prop_12=rc.k12 * n_p1 * n_monomer2 / na_v,
        prop_21=rc.k21 * n_p2 * n_monomer1 / na_v,
        prop_22=rc.k22 * n_p2 * n_monomer2 / na_v,
        term_comb_11=rc.ktc11 * n_p1 * (n_p1 - 1) / na_v,
        term_disp_11=rc.ktd11 * n_p1 * (n_p1 - 1) / na_v,
        term_comb_22=rc.ktc22 * n_p2 * (n_p2 - 1) / na_v,
        term_disp_22=rc.ktd22 * n_p2 * (n_p2 - 1) / na_v,
        term_comb_12=rc.ktc12 * n_p1 * n_p2 / na_v,
        term_disp_12=rc.ktd12 * n_p1 * n_p2 / na_v,
    )
