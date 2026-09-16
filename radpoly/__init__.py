"""radpoly: フリーラジカル(共)重合の kinetic Monte Carlo (kMC) シミュレータ。

開始・伝播・停止(結合/不均化)のみを扱う最小構成。二元共重合(端末モデル)対応。
単一モノマーは RateConstants.homopolymer() + n_monomer2_0=0 の退化ケースとして扱う。
詳細は CLAUDE.md / README.md 参照。
"""
from .kinetics import (
    AVOGADRO,
    Propensities,
    RateConstants,
    compute_propensities,
    geometric_mean_cross_termination,
)
from .state import ReactorState
from .simulator import SimulationResult, TrajectoryPoint, run_kmc
from .analysis import (
    ChainLengthStats,
    chain_length_distribution,
    compute_stats,
    copolymer_composition,
    reactivity_ratios,
    theoretical_instantaneous_composition,
    theoretical_pdi_limit,
    theoretical_rp,
)

__all__ = [
    "AVOGADRO",
    "RateConstants",
    "Propensities",
    "compute_propensities",
    "geometric_mean_cross_termination",
    "ReactorState",
    "SimulationResult",
    "TrajectoryPoint",
    "run_kmc",
    "ChainLengthStats",
    "chain_length_distribution",
    "compute_stats",
    "copolymer_composition",
    "reactivity_ratios",
    "theoretical_instantaneous_composition",
    "theoretical_pdi_limit",
    "theoretical_rp",
]
