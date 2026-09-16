# radpoly-kmc

フリーラジカル(共)重合の反応挙動を kinetic Monte Carlo (kMC) でシミュレーションするライブラリ。
開始・伝播・停止(結合/不均化)のみを扱う最小構成。単一モノマー(ホモポリマー)と
二元共重合(端末モデル)の両方に対応。

TU Clausthal の [mcPolymer](https://www.itc.tu-clausthal.de/en/research/mcpolymer)
(ATRP/NMP/RAFT・分岐・共重合まで扱う本格的な研究用ソフトウェア)に着想を得ているが、
**そのポートやクローンではない独立実装**。理論的基盤・スコープの詳細は [CLAUDE.md](CLAUDE.md) 参照。

## インストール

```bash
pip install -r requirements.txt
```

## 使い方(CLI)

### 単一モノマー(ホモポリマー)

```bash
python -m radpoly.cli homo \
  --kd 1e-5 --f 0.6 --ki 1e3 --kp 1e3 --ktc 1e8 --ktd 0 \
  --conc-I 0.01 --conc-M 5.0 \
  --n-monomer0 100000 --seed 42 \
  --out-trajectory traj.csv --out-mwd mwd.csv
```

### 二元共重合

```bash
python -m radpoly.cli copo \
  --kd 1e-5 --f 0.6 --ki1 1e3 --ki2 1e3 \
  --k11 1e3 --k12 5e2 --k21 8e2 --k22 1e3 \
  --ktc11 1e8 --ktc22 1e8 --ktd11 0 --ktd22 0 \
  --conc-I 0.01 --conc-M1 2.5 --conc-M2 2.5 \
  --n-monomer0 100000 --seed 42
```

`--ktc12` / `--ktd12`(交差停止速度定数)を省略すると、化学制御モデルの幾何平均近似
`sqrt(ktc11*ktc22)` が自動的に使われる(近似であることが警告表示される)。

### 実例: スチレン(AIBN開始、60℃)

裏取り済みの実測値に近い速度定数を使った例(出典・注意点は [CLAUDE.md](CLAUDE.md) の
「スチレンの実例」参照)。**AIBNは60℃では半減期が約20時間と遅いため、意味のある
転化率に達するには反応時間(シミュレーション上の化学時間)が数日相当になる
(計算自体は数秒で終わる。「反応時間が長い」= 遅い開始剤の物理として正しい)。**

```bash
python -m radpoly.cli homo \
  --kd 9.44e-6 --f 0.6 --ki 163 --kp 163 --ktc 1e8 --ktd 0 \
  --conc-I 0.01 --conc-M 8.7 \
  --n-monomer0 20000000 --max-events 30000000 --seed 42
```

もっと速く進む例が見たい場合は70℃相当の値を使う(CLAUDE.md参照)。

## 使い方(ライブラリ)

```python
from radpoly import RateConstants, run_kmc, chain_length_distribution, compute_stats, AVOGADRO

rc = RateConstants.homopolymer(kd=1e-5, f=0.6, ki=1e3, kp=1e3, ktc=1e8, ktd=0.0)
n_monomer0 = 100_000
volume_L = n_monomer0 / (AVOGADRO * 5.0)  # [M]0 = 5.0 mol/L

result = run_kmc(rc, volume_L, n_initiator0=200, n_monomer1_0=n_monomer0, seed=42)
lengths = chain_length_distribution(result)
stats = compute_stats(lengths)
print(f"転化率={result.trajectory[-1].conversion:.3f}, Xn={stats.mn:.1f}, PDI={stats.pdi:.2f}")
```

## 単位系

- `kd`: 1/s、`f`: 無次元(0-1)、その他の速度定数: L/(mol・s)
- 濃度: mol/L

詳細な反応スキーム・単位系の慣例・実装上の注意点は [CLAUDE.md](CLAUDE.md) を参照。

## 検証

```bash
python -m pytest tests -q
```

質量収支の保存、既知の解析解(純粋な二次減衰)、瞬間PDI理論極限(1.5/結合、2.0/不均化)、
Mayo-Lewis組成式との照合を含む。

## 既知の制限

- ATRP・NMP・RAFT・連鎖移動・分岐(バックバイティング等)・3種以上のモノマー共重合・
  ゲル効果は未対応(スコープ外)。
- バッチ重合を高転化率まで進めると、序盤と終盤で鎖長が変わる「連鎖長ドリフト」により
  累積PDIは瞬間理論値より大きくなる(既知の現象。CLAUDE.md参照)。
- 結合・不均化混在時のPDI理論値、共重合系のRp理論値(Walling式)は裏取り不足のため未実装。
- **実測に近い希薄なラジカル濃度では、計算可能な規模で疑似定常状態(QSS)に到達できない。**
  生きているラジカル数が1を大きく下回り、間欠的(バースト的)な挙動になるため、
  `theoretical_rp` はフルパイプラインでは大まかな目安以上の意味を持たない
  (詳細は CLAUDE.md の「重要な既知の限界」参照)。理論値との定量比較は、n_p1を
  厳密に一定に保つ検証テストの枠組みでのみ信頼できる。
