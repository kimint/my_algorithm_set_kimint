"""
최적 볼록 분할: 새 점 없이 대각선만으로, 조각 수가 최소가 되도록 나눈다.

Keil 계열의 동적계획법(DP)
  - 현(chord) (i, j) 로 잘린 부분 다각형 P(i..j) = i, i+1, ..., j 마다
    최소 조각 수 OPT(i, j) 를 구한다.
  - (i, j) 를 변으로 갖는 조각 T 의 꼭짓점 i = t0 < t1 < ... < tk = j 를
    볼록성을 지키며 고르고, T 의 변 (t_m, t_m+1) 이 대각선이면 그 너머의
    부분 다각형 P(t_m..t_m+1) 은 재귀적으로 푼다.
        OPT(i, j) = 1 + min over T  sum OPT(t_m, t_m+1)
  - 최적해의 모든 대각선은 reflex 꼭짓점에 닿아 있다.
    (양 끝이 모두 convex 인 대각선은 지워도 두 조각을 합친 것이 볼록하므로)
    그래서 그런 대각선만 후보로 쓴다.

Keil–Snoeyink 의 O(n + r^2 min(r^2, n)) 구현과 속도는 다르지만,
최적해의 조각 수는 유일하므로 같은 값을 구한다.
"""

import sys

from convex_partition import ConvexPartition, cross, EPS

INF = float("inf")


def optimal_partition(points):
    """최소 조각 수가 되는 대각선 리스트를 돌려준다."""
    cp = ConvexPartition(points, verbose=False)
    n, P = cp.n, cp.P
    reflex = [cp.is_reflex(v) for v in range(n)]

    # out[a] : T 의 변 (a, b) 가 될 수 있는 b (> a)
    #          b == a + 1 이면 다각형의 변, 아니면 reflex 꼭짓점에 닿는 대각선
    out = [[] for _ in range(n)]
    for a in range(n):
        vis = cp.visible(a)
        for b in range(a + 1, n):
            if b == a + 1 or (b in vis and (reflex[a] or reflex[b])):
                out[a].append(b)

    def left(x, y, z):
        return cross(P[x], P[y], P[z]) >= -EPS  # 180도까지 허용

    memo = {}

    def solve(i, j):
        """(최소 조각 수, T 의 꼭짓점 사슬)"""
        if (i, j) in memo:
            return memo[(i, j)]

        # F[(a, b)] = (i 에서 시작해 ..., a, b 로 끝나는 볼록 사슬의 최소 비용, 직전 점)
        F = {}
        incoming = {}

        def add(a, b, cost, prev):
            if cost < F.get((a, b), (INF,))[0]:
                if (a, b) not in F:
                    incoming.setdefault(b, []).append(a)
                F[(a, b)] = (cost, prev)

        for b in out[i]:
            if b < j and left(j, i, b):
                c = 0 if b == i + 1 else solve(i, b)[0]
                add(i, b, c, None)

        for a in range(i + 1, j):
            if a not in incoming:
                continue
            for b in out[a]:
                if b > j:
                    break
                best, arg = INF, None
                for x in incoming[a]:
                    if left(x, a, b) and F[(x, a)][0] < best:
                        best, arg = F[(x, a)][0], x
                if arg is None:
                    continue
                c = 0 if b == a + 1 else solve(a, b)[0]
                add(a, b, best + c, arg)

        best, arg = INF, None
        for a in incoming.get(j, []):
            if left(a, j, i) and F[(a, j)][0] < best:
                best, arg = F[(a, j)][0], a

        if arg is None:
            memo[(i, j)] = (INF, None)
            return memo[(i, j)]

        chain = [j, arg]
        a, b = arg, j
        while True:
            prev = F[(a, b)][1]
            if prev is None:
                break
            a, b = prev, a
            chain.append(a)   # 마지막에 추가되는 점은 i
        chain.reverse()
        memo[(i, j)] = (1 + best, chain)
        return memo[(i, j)]

    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, 10 * n + 1000))
    try:
        total, _ = solve(0, n - 1)
        diagonals = []

        def collect(i, j):
            chain = memo[(i, j)][1]
            for a, b in zip(chain, chain[1:]):
                if b != a + 1:
                    diagonals.append((a, b))
                    collect(a, b)

        collect(0, n - 1)
    finally:
        sys.setrecursionlimit(old_limit)

    assert len(diagonals) + 1 == total
    return diagonals
