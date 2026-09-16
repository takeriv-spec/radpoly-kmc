# CLAUDE.md — mc-polymer-kmc

> 言語: このプロジェクトのドキュメント・コミット・Claude の応答はすべて **日本語**。

## このプロジェクトは何か

フリーラジカル(共)重合の反応挙動を **kinetic Monte Carlo (kMC)** でシミュレーションするライブラリ。
TU Clausthal の `mcPolymer`(ATRP/NMP/RAFT・分岐・共重合まで扱う本格的な研究用ソフトウェア、
<https://www.itc.tu-clausthal.de/en/research/mcpolymer>)に着想を得ているが、**そのポートやクローン
ではない**。開始・伝播・停止(結合/不均化)のみを扱う最小構成の独立実装。

対象:
- 単一モノマーのフリーラジカル重合(ホモポリマー)
- 二元共重合(端末モデル / Mayo-Lewis型)

対象外(意図的にスコープ外):
- ATRP・NMP・RAFTなどの制御ラジカル重合
- 連鎖移動(モノマー・溶媒・CTAへの移動)
- バックバイティング等の分岐反応
- 3種類以上のモノマーの共重合
- ゲル効果(拡散律速停止)、非等温系

## 理論的基盤

- Gillespie, D.T. "Exact Stochastic Simulation of Coupled Chemical Reactions"
  J. Phys. Chem. 81, 2340 (1977). kMCアルゴリズム(Direct Method)と、
  メソスコピック(確率的)速度定数への変換式の根拠。
- Odian, G. "Principles of Polymerization" 4th ed., Wiley (2004), Ch.3, 6.
  古典的なフリーラジカル(共)重合の速度論の定義・慣例。
- Mayo, F.R.; Lewis, F.M. J. Am. Chem. Soc. 1944, 66, 1594.
  端末モデルの瞬間共重合体組成式(Mayo-Lewis式)。

**方針: 裏取りできていない式はコードに載せない。** 単一機構(結合のみ/不均化のみ)の
PDI理論極限(1.5/2.0)とMayo-Lewis組成式は検証済みで実装している。一方、結合・不均化
混在時のPDI式や、共重合系のRp(Walling式など交差項φを要する式)は裏取りが不十分なため
未実装([radpoly/analysis.py](radpoly/analysis.py) の `theoretical_pdi_limit`, `theoretical_rp`
の docstring 参照、`NotImplementedError` を送出する)。

## 実装した反応(端末モデル、開始・伝播・停止のみ)

```
I        --kd,f-->  2 R・              (開始剤分解、効率fでケージ効果を考慮)
R・ + M1 --ki1-->    P1・
R・ + M2 --ki2-->    P2・
P1・+ M1 --k11-->    P1・              (伝播、末端モノマーで速度定数が決まる)
P1・+ M2 --k12-->    P2・
P2・+ M1 --k21-->    P1・
P2・+ M2 --k22-->    P2・
P1・+P1・ --ktc11/ktd11--> 停止(結合/不均化)
P2・+P2・ --ktc22/ktd22--> 停止(結合/不均化)
P1・+P2・ --ktc12/ktd12--> 停止(結合/不均化、交差)
```

単一モノマー(ホモポリマー)は `RateConstants.homopolymer()` + `n_monomer2_0=0` の
退化ケースとして扱う(M2側の反応はpropensity=0で自然に発火しなくなる)。

## 単位系(重要、変更時は要注意)

- `kd`: 1/s(単分子反応)。`f`: 無次元(0-1)。
- `ki1,ki2,k11,k12,k21,k22,ktcXY,ktdXY`: L/(mol・s)。
- `ktcXY, ktdXY` は `d[P・]/dt|term = -2*kt*[P・]^2` という慣例で定義している
  (kt = 対応するktc+ktdの和)。この慣例に基づき、[radpoly/kinetics.py](radpoly/kinetics.py) の
  `compute_propensities` で分子数ベースのpropensityに変換している。**この慣例を変えると
  停止関連のpropensity式・`theoretical_rp`の導出が両方とも影響を受ける**ので、
  変更時は両方を合わせて直すこと。

## 実装上の注意(kMC設計、変更時に踏まえること)

- **分子数は「シミュレーション体積」内の実際の分子数**を直接扱う(モル濃度に換算しない)。
  `volume_L` で分子数と濃度(mol/L)を結びつける。大きな `n_monomer0` ほど統計精度は
  上がるが計算コストが増える(純Pythonで概ね100万イベント/数秒〜数十秒)。
- **瞬間PDI理論値(1.5/2.0)とバッチ全体の累積PDIは別物。** バッチを転化率100%近くまで
  進めると、序盤([M]大、鎖が長く育つ)と終盤([M]小、鎖が短い)の鎖が混ざり、
  累積PDIは瞬間理論値よりずっと大きくなる(「連鎖長ドリフト」)。低転化率で打ち切れば
  瞬間理論値と比較できる([tests/test_simulator.py](tests/test_simulator.py) 参照)。
- **PDI理論値との照合には、生きているラジカル数(n_p1)が十分大きく安定した
  疑似定常状態が必要。** 実際の開始剤分解を使うテストではQSSに達するまでの
  過渡応答が混入しやすいため、`tests/test_simulator.py` の
  `test_core_termination_mechanism_matches_theoretical_pdi_*` では意図的に
  「消費された分だけ即座に長さ1の鎖を補充する」構成でn_p1を厳密に一定に保ち、
  伝播・停止のコア機構だけを直接検証している。

## ディレクトリ構成

```
radpoly/
  kinetics.py    速度定数(RateConstants)とGillespie propensity計算
  state.py       反応系の状態(分子数・生きている鎖/死鎖の組成)
  simulator.py   Gillespie Direct Methodのメインループ(run_kmc)
  analysis.py    Mn/Mw/PDI計算、理論値との比較(検証用)
  cli.py         CLIエントリポイント(homo/copoサブコマンド)
tests/           pytest。質量収支・理論値照合が必須ケース
examples/        典型的パラメータでの実行例
```

## テスト

```bash
python -m pytest tests -q
```

## labochro.tech との関係

このライブラリを `D:/CaudePJ/Research/labochro` に vendored コピーし、Web版ツールとして
公開する計画がある(進行中の場合は labochro 側の CLAUDE.md も参照)。**式の実装・修正は
必ずこのリポジトリ側で行い、labochro側は vendored コピーを更新するだけ**にすること
(viscosity-predictor・chain_fit と同じ運用)。
