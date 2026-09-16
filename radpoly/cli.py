"""コマンドラインからkMCシミュレーションを実行するエントリポイント。

2つのサブコマンド:
- homo: 単一モノマーのフリーラジカル重合
- copo: 二元共重合(端末モデル)
"""
from __future__ import annotations

import argparse
import csv
import sys

from .analysis import (
    chain_length_distribution,
    compute_stats,
    copolymer_composition,
    theoretical_pdi_limit,
    theoretical_rp,
)
from .kinetics import AVOGADRO, RateConstants, geometric_mean_cross_termination
from .simulator import run_kmc


def _add_common_output_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--n-monomer0", type=int, default=100_000, help="シミュレーション上のモノマー総分子数(既定10万)")
    p.add_argument("--max-time", type=float, default=float("inf"), help="最大シミュレーション時間 [s]")
    p.add_argument("--max-events", type=int, default=10_000_000, help="最大イベント数(安全上限)")
    p.add_argument("--seed", type=int, default=None, help="乱数シード")
    p.add_argument("--out-trajectory", type=str, default=None, help="転化率-時間データのCSV出力先")
    p.add_argument("--out-mwd", type=str, default=None, help="分子量分布(鎖長リスト)のCSV出力先")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ラジカル(共)重合の kinetic Monte Carlo シミュレーション")
    sub = p.add_subparsers(dest="mode", required=True)

    homo = sub.add_parser("homo", help="単一モノマーのフリーラジカル重合")
    homo.add_argument("--kd", type=float, required=True, help="開始剤分解速度定数 [1/s]")
    homo.add_argument("--f", type=float, default=0.6, help="開始剤効率 (既定値0.6)")
    homo.add_argument("--ki", type=float, required=True, help="連鎖開始速度定数 [L/(mol・s)]")
    homo.add_argument("--kp", type=float, required=True, help="伝播速度定数 [L/(mol・s)]")
    homo.add_argument("--ktc", type=float, default=0.0, help="停止(結合)速度定数 [L/(mol・s)]")
    homo.add_argument("--ktd", type=float, default=0.0, help="停止(不均化)速度定数 [L/(mol・s)]")
    homo.add_argument("--conc-I", type=float, required=True, help="初期開始剤濃度 [mol/L]")
    homo.add_argument("--conc-M", type=float, required=True, help="初期モノマー濃度 [mol/L]")
    _add_common_output_args(homo)

    copo = sub.add_parser("copo", help="二元共重合(端末モデル)")
    copo.add_argument("--kd", type=float, required=True, help="開始剤分解速度定数 [1/s]")
    copo.add_argument("--f", type=float, default=0.6, help="開始剤効率 (既定値0.6)")
    copo.add_argument("--ki1", type=float, required=True, help="連鎖開始 R・+M1 [L/(mol・s)]")
    copo.add_argument("--ki2", type=float, required=True, help="連鎖開始 R・+M2 [L/(mol・s)]")
    copo.add_argument("--k11", type=float, required=True, help="伝播 P1・+M1 [L/(mol・s)]")
    copo.add_argument("--k12", type=float, required=True, help="伝播 P1・+M2 [L/(mol・s)]")
    copo.add_argument("--k21", type=float, required=True, help="伝播 P2・+M1 [L/(mol・s)]")
    copo.add_argument("--k22", type=float, required=True, help="伝播 P2・+M2 [L/(mol・s)]")
    copo.add_argument("--ktc11", type=float, default=0.0, help="停止(結合) P1-P1 [L/(mol・s)]")
    copo.add_argument("--ktc22", type=float, default=0.0, help="停止(結合) P2-P2 [L/(mol・s)]")
    copo.add_argument("--ktc12", type=float, default=None, help="停止(結合) P1-P2交差 [L/(mol・s)] (省略時はktc11・ktc22の幾何平均)")
    copo.add_argument("--ktd11", type=float, default=0.0, help="停止(不均化) P1-P1 [L/(mol・s)]")
    copo.add_argument("--ktd22", type=float, default=0.0, help="停止(不均化) P2-P2 [L/(mol・s)]")
    copo.add_argument("--ktd12", type=float, default=None, help="停止(不均化) P1-P2交差 [L/(mol・s)] (省略時はktd11・ktd22の幾何平均)")
    copo.add_argument("--conc-I", type=float, required=True, help="初期開始剤濃度 [mol/L]")
    copo.add_argument("--conc-M1", type=float, required=True, help="初期モノマー1濃度 [mol/L]")
    copo.add_argument("--conc-M2", type=float, required=True, help="初期モノマー2濃度 [mol/L]")
    _add_common_output_args(copo)

    return p


def _write_trajectory_csv(path: str, trajectory) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "conversion", "n_p1", "n_p2"])
        for pt in trajectory:
            w.writerow([pt.time, pt.conversion, pt.n_p1, pt.n_p2])


def _write_mwd_csv(path: str, lengths: list[int]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["chain_length"])
        for l in lengths:
            w.writerow([l])


def _run_homo(args: argparse.Namespace) -> int:
    rc = RateConstants.homopolymer(kd=args.kd, f=args.f, ki=args.ki, kp=args.kp, ktc=args.ktc, ktd=args.ktd)

    n_initiator0 = round(args.conc_I / args.conc_M * args.n_monomer0)
    volume_L = args.n_monomer0 / (AVOGADRO * args.conc_M)

    result = run_kmc(
        rc, volume_L, n_initiator0, args.n_monomer0, n_monomer2_0=0,
        max_time=args.max_time, max_events=args.max_events, seed=args.seed,
    )

    final = result.trajectory[-1]
    print(f"停止理由: {result.stopped_reason}")
    print(f"到達時間: {final.time:.4g} s, 転化率: {final.conversion:.4f}")

    lengths = chain_length_distribution(result)
    stats = compute_stats(lengths)
    print(f"生成鎖数: {stats.n_chains}, Xn={stats.mn:.2f}, Xw={stats.mw:.2f}, PDI={stats.pdi:.3f}")

    try:
        rp = theoretical_rp(rc, args.conc_I, args.conc_M)
        print(f"[検証用] 定常状態理論によるRp(初期時点) = {rp:.4g} mol/(L・s)")
    except (ValueError, NotImplementedError) as e:
        print(f"[検証用] 理論Rp計算スキップ: {e}")

    try:
        pdi_theory = theoretical_pdi_limit(rc)
        print(f"[検証用] 停止機構による理論PDI極限 = {pdi_theory:.3f}")
    except NotImplementedError as e:
        print(f"[検証用] 理論PDI計算スキップ: {e}")

    if args.out_trajectory:
        _write_trajectory_csv(args.out_trajectory, result.trajectory)
    if args.out_mwd:
        _write_mwd_csv(args.out_mwd, lengths)
    return 0


def _run_copo(args: argparse.Namespace) -> int:
    ktc12 = args.ktc12 if args.ktc12 is not None else geometric_mean_cross_termination(args.ktc11, args.ktc22)
    ktd12 = args.ktd12 if args.ktd12 is not None else geometric_mean_cross_termination(args.ktd11, args.ktd22)
    if args.ktc12 is None:
        print(f"[注意] ktc12 未指定のため幾何平均近似を使用: ktc12={ktc12:.4g}")
    if args.ktd12 is None:
        print(f"[注意] ktd12 未指定のため幾何平均近似を使用: ktd12={ktd12:.4g}")

    rc = RateConstants(
        kd=args.kd, f=args.f, ki1=args.ki1, ki2=args.ki2,
        k11=args.k11, k12=args.k12, k21=args.k21, k22=args.k22,
        ktc11=args.ktc11, ktc12=ktc12, ktc22=args.ktc22,
        ktd11=args.ktd11, ktd12=ktd12, ktd22=args.ktd22,
    )

    conc_M_total = args.conc_M1 + args.conc_M2
    n_monomer1_0 = round(args.conc_M1 / conc_M_total * args.n_monomer0)
    n_monomer2_0 = args.n_monomer0 - n_monomer1_0
    n_initiator0 = round(args.conc_I / conc_M_total * args.n_monomer0)
    volume_L = args.n_monomer0 / (AVOGADRO * conc_M_total)

    result = run_kmc(
        rc, volume_L, n_initiator0, n_monomer1_0, n_monomer2_0=n_monomer2_0,
        max_time=args.max_time, max_events=args.max_events, seed=args.seed,
    )

    final = result.trajectory[-1]
    print(f"停止理由: {result.stopped_reason}")
    print(f"到達時間: {final.time:.4g} s, 転化率: {final.conversion:.4f}")

    lengths = chain_length_distribution(result)
    stats = compute_stats(lengths)
    print(f"生成鎖数: {stats.n_chains}, Xn={stats.mn:.2f}, Xw={stats.mw:.2f}, PDI={stats.pdi:.3f}")

    n1, n2 = copolymer_composition(result)
    if n1 + n2 > 0:
        print(f"共重合体組成(累積): M1モル分率={n1 / (n1 + n2):.4f} (M1={n1}, M2={n2})")

    if args.out_trajectory:
        _write_trajectory_csv(args.out_trajectory, result.trajectory)
    if args.out_mwd:
        _write_mwd_csv(args.out_mwd, lengths)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.mode == "homo":
        return _run_homo(args)
    return _run_copo(args)


if __name__ == "__main__":
    sys.exit(main())
