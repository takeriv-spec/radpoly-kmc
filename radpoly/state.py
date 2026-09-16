"""kMC シミュレーションの反応系の状態(分子数と鎖の組成)。

二元共重合を扱うため、各ラジカル鎖は末端モノマー種によって2つのプールに分ける:
- live_1: 末端がモノマー1(last unit = M1)の生きているラジカル鎖
- live_2: 末端がモノマー2(last unit = M2)の生きているラジカル鎖
各鎖は [n1, n2] (取り込んだM1・M2のモノマー単位数)のリストで表現する
(末端種はどちらのプールに入っているかで分かるので別途保持しない)。
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field


@dataclass
class ReactorState:
    n_initiator: int
    n_radical: int  # 開始剤分解直後、まだモノマーに付加していない一次ラジカルの数
    n_monomer1: int
    n_monomer2: int
    live_1: list = field(default_factory=list)  # list[list[int, int]] : [n1, n2]、末端=M1
    live_2: list = field(default_factory=list)  # list[list[int, int]] : [n1, n2]、末端=M2
    dead_chains: list = field(default_factory=list)  # list[tuple[int, int]] : (n1, n2)
    time: float = 0.0

    @property
    def n_p1(self) -> int:
        return len(self.live_1)

    @property
    def n_p2(self) -> int:
        return len(self.live_2)

    def monomer1_in_chains(self) -> int:
        return (
            sum(c[0] for c in self.live_1)
            + sum(c[0] for c in self.live_2)
            + sum(d[0] for d in self.dead_chains)
        )

    def monomer2_in_chains(self) -> int:
        return (
            sum(c[1] for c in self.live_1)
            + sum(c[1] for c in self.live_2)
            + sum(d[1] for d in self.dead_chains)
        )

    def check_mass_balance(self, n1_initial: int, n2_initial: int) -> None:
        """モノマー総量の保存を検証する(質量収支チェック、M1・M2それぞれ)。"""
        t1 = self.n_monomer1 + self.monomer1_in_chains()
        t2 = self.n_monomer2 + self.monomer2_in_chains()
        if t1 != n1_initial or t2 != n2_initial:
            raise AssertionError(
                f"モノマー質量収支が破れています: "
                f"M1 現在合計={t1}(初期{n1_initial}), M2 現在合計={t2}(初期{n2_initial})"
            )

    @staticmethod
    def _pop_random(chains: list, rng: random.Random) -> list:
        """指定リストからランダムに1本選んで取り除き、その鎖([n1,n2])を返す。

        O(1) の swap-remove(鎖の識別に意味は無いため順序保持は不要)。
        """
        idx = rng.randrange(len(chains))
        item = chains[idx]
        last = chains.pop()
        if idx < len(chains):
            chains[idx] = last
        return item

    def do_decomposition(self, rng: random.Random, f: float) -> None:
        self.n_initiator -= 1
        if rng.random() < f:
            self.n_radical += 2
        # (1-f) の確率でケージ内再結合により失活(何も生成しない)

    def do_chain_initiation_1(self) -> None:
        self.n_radical -= 1
        self.n_monomer1 -= 1
        self.live_1.append([1, 0])

    def do_chain_initiation_2(self) -> None:
        self.n_radical -= 1
        self.n_monomer2 -= 1
        self.live_2.append([0, 1])

    def do_prop_11(self, rng: random.Random) -> None:
        idx = rng.randrange(len(self.live_1))
        self.live_1[idx][0] += 1
        self.n_monomer1 -= 1

    def do_prop_22(self, rng: random.Random) -> None:
        idx = rng.randrange(len(self.live_2))
        self.live_2[idx][1] += 1
        self.n_monomer2 -= 1

    def do_prop_12(self, rng: random.Random) -> None:
        """P1・+ M2 -> P2・ (末端がM2に変わるため live_1 から live_2 へ移動)。"""
        chain = self._pop_random(self.live_1, rng)
        chain[1] += 1
        self.n_monomer2 -= 1
        self.live_2.append(chain)

    def do_prop_21(self, rng: random.Random) -> None:
        """P2・+ M1 -> P1・ (末端がM1に変わるため live_2 から live_1 へ移動)。"""
        chain = self._pop_random(self.live_2, rng)
        chain[0] += 1
        self.n_monomer1 -= 1
        self.live_1.append(chain)

    def do_term_comb(self, rng: random.Random, pool_a: list, pool_b: list) -> None:
        """結合停止。pool_a, pool_b は同じプール(同種停止)でも異なるプール(交差停止)でもよい。

        pool_a と pool_b が同一オブジェクトの場合、1本目を取り除いた後の残りから
        2本目を選ぶことになり正しく相異なる2本が選ばれる。
        """
        a = self._pop_random(pool_a, rng)
        b = self._pop_random(pool_b, rng)
        self.dead_chains.append((a[0] + b[0], a[1] + b[1]))

    def do_term_disp(self, rng: random.Random, pool_a: list, pool_b: list) -> None:
        a = self._pop_random(pool_a, rng)
        b = self._pop_random(pool_b, rng)
        self.dead_chains.append((a[0], a[1]))
        self.dead_chains.append((b[0], b[1]))
