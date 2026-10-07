"""
스택 기반 볼록 분할(convex partition) 알고리즘

입력 : 반시계 방향(CCW) 순서로 나열된 단순 다각형의 꼭짓점 좌표 리스트
출력 : 다각형을 볼록 다각형들로 나누는 대각선 리스트 (꼭짓점 인덱스 쌍)

사용법:
    python convex_partition.py                  # 내장 예제 전부 실행
    python convex_partition.py polygon.txt      # 파일에서 다각형 읽기 (한 줄에 "x y")
    python convex_partition.py --save out.png   # 그림을 화면 대신 파일로 저장
    python convex_partition.py --quiet          # 단계별 로그 끄기
    python convex_partition.py --random 50      # 꼭짓점 50개짜리 무작위 다각형 100개로 테스트
    python convex_partition.py --compare 30     # 꼭짓점 30개짜리 무작위 다각형으로 최적해와 비교
        (옵션: --count 개수, --seed 시드)
"""

import math
import sys

EPS = 1e-9


# ---------------------------------------------------------------------------
# 기본 기하 함수
# ---------------------------------------------------------------------------

def cross(o, a, b):
    """(a - o) x (b - o). 양수면 o->a->b 가 왼쪽(반시계) 회전."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def signed_area(pts):
    s = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % len(pts)]
        s += x1 * y2 - x2 * y1
    return s / 2.0


def _on_segment(p, q, r):
    """p, q, r 이 한 직선 위에 있을 때 r 이 선분 pq 위에 있는지."""
    return (min(p[0], q[0]) - EPS <= r[0] <= max(p[0], q[0]) + EPS and
            min(p[1], q[1]) - EPS <= r[1] <= max(p[1], q[1]) + EPS)


def segments_intersect(p1, p2, q1, q2):
    """두 선분이 만나는지 (끝점이 닿는 경우도 포함)."""
    d1 = cross(q1, q2, p1)
    d2 = cross(q1, q2, p2)
    d3 = cross(p1, p2, q1)
    d4 = cross(p1, p2, q2)
    if (((d1 > EPS and d2 < -EPS) or (d1 < -EPS and d2 > EPS)) and
            ((d3 > EPS and d4 < -EPS) or (d3 < -EPS and d4 > EPS))):
        return True
    if abs(d1) <= EPS and _on_segment(q1, q2, p1):
        return True
    if abs(d2) <= EPS and _on_segment(q1, q2, p2):
        return True
    if abs(d3) <= EPS and _on_segment(p1, p2, q1):
        return True
    if abs(d4) <= EPS and _on_segment(p1, p2, q2):
        return True
    return False


def is_convex_polygon(pts):
    n = len(pts)
    return all(cross(pts[i - 1], pts[i], pts[(i + 1) % n]) >= -EPS for i in range(n))


def _in_front(a1, a2, b1, b2, p):
    """서로 교차하지 않는 선분 A=a1a2, B=b1b2 를 점 p 에서 볼 때 A 가 B 보다 앞(가까이)인지."""
    sp = cross(a1, a2, p)
    if abs(sp) > EPS:
        s = 1.0 if sp > 0 else -1.0
        # B 가 직선 A 를 기준으로 p 의 반대편에 있으면 A 가 앞
        if s * cross(a1, a2, b1) <= EPS and s * cross(a1, a2, b2) <= EPS:
            return True
    tp = cross(b1, b2, p)
    if abs(tp) > EPS:
        t = 1.0 if tp > 0 else -1.0
        # A 가 직선 B 를 기준으로 p 와 같은 편에 있으면 A 가 앞
        if t * cross(b1, b2, a1) >= -EPS and t * cross(b1, b2, a2) >= -EPS:
            return True
    return False


class _ActiveEdge:
    """회전 스윕에서 현재 광선과 만나는 변. 가까운 순서로 정렬된다."""
    __slots__ = ("a", "b", "p")

    def __init__(self, a, b, p):
        self.a, self.b, self.p = a, b, p

    def __lt__(self, other):
        return _in_front(self.a, self.b, other.a, other.b, self.p)


# ---------------------------------------------------------------------------
# 알고리즘
# ---------------------------------------------------------------------------

class ConvexPartition:
    def __init__(self, points, verbose=True, fast=True):
        if len(points) < 3:
            raise ValueError("꼭짓점이 3개 이상이어야 합니다.")
        self.P = [(float(x), float(y)) for x, y in points]
        if signed_area(self.P) <= 0:
            raise ValueError("꼭짓점은 반시계 방향(CCW) 순서로 입력해야 합니다.")
        self.n = len(self.P)
        if not self._is_simple():
            raise ValueError("변끼리 교차합니다. 단순 다각형이어야 합니다.")
        self.verbose = verbose
        self.diagonals = []                       # 그린 순서대로 (i, j)
        self.incident = {i: set() for i in range(self.n)}  # 꼭짓점별 연결된 대각선 끝점
        self.fallback_vertices = []               # 대각선 하나로 해결 못 한 꼭짓점

        # 현재 조각(face)들: 대각선을 그을 때마다 조각이 둘로 나뉜다
        self.faces = {0: list(range(self.n))}
        self.vfaces = [{0} for _ in range(self.n)]  # 꼭짓점별 속한 조각 번호
        self._next_face = 1

        # fast=True: 꼭짓점별 visibility 를 (처음 필요할 때 한 번) 계산해서 대각선 판정을 O(1) 로
        # fast=False: 매번 모든 변과 교차 검사 (예전 방식, 비교용)
        self.fast = fast
        self._vis = [None] * self.n

    def visible(self, i):
        """i 와 대각선으로 이을 수 있는 꼭짓점 집합 (캐시)."""
        if self._vis[i] is None:
            self._vis[i] = self._visible_from(i)
        return self._vis[i]

    # ----- 보조 함수 -------------------------------------------------------

    def _log(self, msg):
        if self.verbose:
            print(msg)

    def _is_simple(self):
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if j == i + 1 or (i == 0 and j == self.n - 1):
                    continue  # 이웃한 변
                if segments_intersect(self.P[i], self.P[(i + 1) % self.n],
                                      self.P[j], self.P[(j + 1) % self.n]):
                    return False
        return True

    def nxt(self, i):
        return (i + 1) % self.n

    def prv(self, i):
        return (i - 1) % self.n

    def adjacent(self, i, j):
        return j == self.nxt(i) or j == self.prv(i)

    def is_reflex(self, i):
        """원래 다각형 P 에서 꼭짓점 i 가 reflex(내각 > 180도) 인지."""
        return cross(self.P[self.prv(i)], self.P[i], self.P[self.nxt(i)]) < -EPS

    def _in_cone(self, i, j):
        """선분 i-j 가 꼭짓점 i 에서 다각형 내부 쪽으로 출발하는지."""
        a, b = self.P[i], self.P[j]
        a0, a1 = self.P[self.prv(i)], self.P[self.nxt(i)]
        if cross(a, a1, a0) >= 0:          # i 가 convex
            return cross(a, b, a0) > EPS and cross(b, a, a1) > EPS
        # i 가 reflex
        return not (cross(a, b, a1) >= -EPS and cross(b, a, a0) >= -EPS)

    def _visible_from(self, i):
        """
        꼭짓점 i 에서 다각형 내부를 통해 보이는 꼭짓점들 (= i 와 대각선으로 이을 수 있는 점).
        i 를 중심으로 광선을 반시계 방향으로 한 바퀴 돌리며(회전 스윕),
        광선과 만나는 변들을 가까운 순서로 유지한다.
        """
        import bisect

        n, P, p = self.n, self.P, self.P[i]
        base = math.atan2(P[self.nxt(i)][1] - p[1], P[self.nxt(i)][0] - p[0])
        rel = {}
        for k in range(n):
            if k != i:
                a = math.atan2(P[k][1] - p[1], P[k][0] - p[0])
                rel[k] = (a - base) % (2 * math.pi)
        rel[self.nxt(i)] = 0.0
        interior = rel[self.prv(i)]  # i 의 내각 범위: (0, interior)

        order = sorted(rel, key=lambda k: (rel[k], (P[k][0] - p[0]) ** 2 + (P[k][1] - p[1]) ** 2))
        pos = {k: idx for idx, k in enumerate(order)}

        def crosses_start(u, w):
            # 각도 범위가 시작 광선(각도 0)을 지나는 변
            return abs(rel[u] - rel[w]) > math.pi

        keys = {}
        active = []
        for k in range(n):
            u, w = k, self.nxt(k)
            if i in (u, w):
                continue
            keys[(u, w)] = _ActiveEdge(P[u], P[w], p)
            if crosses_start(u, w):
                bisect.insort(active, keys[(u, w)])

        visible = set()
        prev_rel = None
        for v in order:
            # 같은 방향에 더 가까운 꼭짓점이 있으면 선분이 그 꼭짓점을 지나므로 대각선이 아님
            on_line = prev_rel is not None and abs(rel[v] - prev_rel) <= 1e-12
            prev_rel = rel[v]
            inserts = []
            for (u, w) in ((self.prv(v), v), (v, self.nxt(v))):
                if i in (u, w):
                    continue
                other = u if w == v else w
                later = pos[other] > pos[v]
                if crosses_start(u, w):
                    later = not later
                if later:
                    inserts.append(keys[(u, w)])
                else:
                    active.remove(keys[(u, w)])

            if EPS < rel[v] < interior - EPS and not self.adjacent(i, v) and not on_line:
                if not active or not segments_intersect(p, P[v], active[0].a, active[0].b):
                    visible.add(v)

            for key in inserts:
                bisect.insort(active, key)
        return visible

    def is_diagonal(self, i, j):
        """선분 i-j 가 다각형 P 의 대각선인지 (내부를 지나고 변과 만나지 않음)."""
        if i == j or self.adjacent(i, j):
            return False
        if self.fast:
            return j in self.visible(i)
        return self._is_diagonal_naive(i, j)

    def _is_diagonal_naive(self, i, j):
        if i == j or self.adjacent(i, j):
            return False
        if not (self._in_cone(i, j) and self._in_cone(j, i)):
            return False
        pi, pj = self.P[i], self.P[j]
        for k in range(self.n):
            k2 = self.nxt(k)
            if k in (i, j) or k2 in (i, j):
                continue
            if segments_intersect(pi, pj, self.P[k], self.P[k2]):
                return False
        return True

    def _crosses_existing(self, i, j):
        """i-j 가 이미 그린 대각선과 교차하는지."""
        for u, v in self.diagonals:
            if len({i, j, u, v}) < 4:     # 끝점을 공유하면 교차 아님
                continue
            if segments_intersect(self.P[i], self.P[j], self.P[u], self.P[v]):
                return True
        return False

    def can_draw(self, i, j):
        if not self.is_diagonal(i, j) or j in self.incident[i]:
            return False
        if self.fast:
            # P 의 대각선이면서 두 끝점이 같은 조각에 있으면 기존 대각선과 교차하지 않는다
            return bool(self.vfaces[i] & self.vfaces[j])
        return not self._crosses_existing(i, j)

    def draw(self, i, j):
        self.diagonals.append((i, j))
        self.incident[i].add(j)
        self.incident[j].add(i)
        self._split_face(i, j)
        self._log(f"    >> 대각선 그림: {i} - {j}")

    def _split_face(self, i, j):
        for f in self.vfaces[i] & self.vfaces[j]:
            face = self.faces[f]
            pi, pj = sorted((face.index(i), face.index(j)))
            if pj - pi in (1, len(face) - 1):
                continue  # 이 조각에서는 변
            new = self._next_face
            self._next_face += 1
            self.faces[f] = face[pi:pj + 1]
            self.faces[new] = face[pj:] + face[:pi + 1]
            for v in face[pj + 1:] + face[:pi]:
                self.vfaces[v].discard(f)
                self.vfaces[v].add(new)
            self.vfaces[i].add(new)
            self.vfaces[j].add(new)
            return

    def _max_wedge(self, i, extra=None):
        """
        꼭짓점 i 의 내각이 (이미 그린 대각선 + extra 대각선) 으로 나뉘었을 때
        나뉜 각 중 가장 큰 각 (라디안).
        """
        def ang(k):
            return math.atan2(self.P[k][1] - self.P[i][1], self.P[k][0] - self.P[i][0])

        base = ang(self.nxt(i))  # 내각은 nxt 방향에서 반시계로 prv 방향까지
        others = [self.prv(i)] + list(self.incident[i])
        if extra is not None:
            others.append(extra)
        rel = sorted((ang(k) - base) % (2 * math.pi) for k in others)
        biggest, last = 0.0, 0.0
        for r in rel:
            biggest = max(biggest, r - last)
            last = r
        return biggest

    def is_convex_now(self, i, extra=None):
        """현재 그려진 대각선 (+ extra) 을 고려했을 때 i 가 convex vertex 인지."""
        return self._max_wedge(i, extra) <= math.pi + EPS

    # ----- 알고리즘 본체 ---------------------------------------------------

    def run(self):
        reflex = [i for i in range(self.n) if self.is_reflex(i)]
        self._log(f"꼭짓점 수: {self.n},  reflex 꼭짓점: {reflex}")
        if not reflex:
            self._log("reflex 꼭짓점이 없음 -> 이미 볼록 다각형이므로 대각선이 필요 없음")
            return self.diagonals

        stack = self._phase1(reflex[0])
        stack2 = self._phase2(stack)
        self._phase3(stack2)

        self._log(f"\n결과 대각선 ({len(self.diagonals)}개): {self.diagonals}")
        return self.diagonals

    def _phase1(self, start):
        """Step 1~4: 반시계 방향으로 돌면서 stack 을 만든다."""
        self._log("\n[Step 1~4] stack 만들기")
        vt = start
        stack = [vt]
        self._log(f"  첫 reflex 꼭짓점 {vt} 을 push -> stack = {stack}")

        w = self.nxt(start)
        while w != start:
            # vt-w 가 대각선이 아니면 w 의 이전 점을 vt 로
            if not self.adjacent(vt, w) and not self.is_diagonal(vt, w):
                vt = self.prv(w)
                stack.append(vt)
                self._log(f"  {stack[-2]}-{w} 는 대각선이 아님 -> 이전 점 {vt} push -> stack = {stack}")
                continue  # 새 vt 기준으로 w 를 다시 확인
            if self.is_reflex(w):
                vt = w
                stack.append(vt)
                self._log(f"  reflex 꼭짓점 {vt} 발견 -> push -> stack = {stack}")
            w = self.nxt(w)

        self._log(f"  한 바퀴 다 돌았음 -> 최종 stack = {stack}")
        return stack

    def _phase2(self, stack):
        """Step 5~8: stack 과 stack2 의 맨 위를 연결하며 대각선을 긋는다."""
        self._log("\n[Step 5~8] stack / stack2 처리")
        stack2 = [stack.pop()]
        self._log(f"  stack = {stack}, stack2 = {stack2}")

        while stack:
            a, b = stack[-1], stack2[-1]
            self._log(f"  stack 맨 위 = {a}, stack2 맨 위 = {b}")

            if self.is_convex_now(b):
                # 경우 3: b 가 이미 convex -> 대각선 없이 b 를 버림
                stack2.pop()
                self._log(f"    {b} 는 이미 convex -> stack2 에서 버림")
            elif self.can_draw(a, b):
                becomes_convex = self.is_convex_now(b, extra=a)
                self.draw(a, b)
                if becomes_convex:
                    # 경우 1: 대각선을 그으면 b 가 convex 가 됨
                    stack2.pop()
                    self._log(f"    {b} 가 convex 가 됨 -> stack2 에서 버림")
                else:
                    # 경우 2: 대각선을 그어도 b 가 여전히 reflex
                    self._log(f"    {b} 는 여전히 reflex -> stack2 에 남김")
            else:
                # a-b 가 다각형의 변인 경우 등 (대각선을 그을 수 없음)
                self._log(f"    {a}-{b} 는 대각선으로 그을 수 없음 -> {b} 는 reflex 로 stack2 에 남김")

            stack2.append(stack.pop())
            self._log(f"    stack = {stack}, stack2 = {stack2}")

        return stack2

    def _phase3(self, stack2):
        """Step 9~13: stack2 에 남은 reflex 꼭짓점을 대각선 하나씩으로 해결한다."""
        self._log("\n[Step 9~13] stack2 에 남은 꼭짓점 처리")
        while stack2:
            vt = stack2.pop()
            if self.is_convex_now(vt):
                self._log(f"  {vt}: 이미 convex -> 넘어감")
                continue

            self._log(f"  {vt}: reflex -> 대각선 찾는 중")
            vt2 = self.nxt(vt)
            while vt2 != vt:
                if vt2 in self.incident[vt]:
                    # vt-vt2 는 이미 그어져 있음 -> '대각선 하나로 안 되는 경우' 로 처리
                    self._log(f"    {vt}-{vt2} 는 이미 그어진 대각선 -> 하나로 안 되는 경우로 처리")
                    break
                if self.can_draw(vt, vt2) and self.is_convex_now(vt, extra=vt2):
                    self.draw(vt, vt2)
                    break
                vt2 = self.nxt(vt2)

            if not self.is_convex_now(vt):
                self._fix_dynamically(vt)

    def _fix_dynamically(self, vt):
        """
        대각선 하나로 vt 를 convex 로 만들 수 없을 때:
        vt 가 convex 가 될 때까지, 매번 vt 의 가장 큰 각을 가장 많이 줄이는 대각선을 추가한다.
        """
        self.fallback_vertices.append(vt)
        self._log(f"    대각선 하나로는 안 됨 -> {vt} 가 convex 가 될 때까지 대각선 추가")
        while not self.is_convex_now(vt):
            current = self._max_wedge(vt)
            best, best_wedge = None, current
            for w in range(self.n):
                if self.can_draw(vt, w):
                    wedge = self._max_wedge(vt, extra=w)
                    if wedge < best_wedge - EPS:
                        best, best_wedge = w, wedge
            if best is None:
                self._log(f"    경고: {vt} 의 각을 더 줄일 수 있는 대각선이 없음")
                return
            self.draw(vt, best)
            self._log(f"       {vt} 의 가장 큰 각: {math.degrees(current):.1f}도 -> "
                      f"{math.degrees(best_wedge):.1f}도")

    # ----- 결과 확인 -------------------------------------------------------

    def pieces(self):
        """대각선으로 나뉜 조각들 (꼭짓점 인덱스 리스트)."""
        return list(self.faces.values())

    def report(self):
        pcs = self.pieces()
        bad = [p for p in pcs if not is_convex_polygon([self.P[k] for k in p])]
        print(f"조각 수: {len(pcs)}")
        for p in pcs:
            ok = "볼록" if p not in bad else "볼록 아님!"
            print(f"  {p}  ({ok})")
        print("=> 모든 조각이 볼록함" if not bad else "=> 볼록하지 않은 조각이 있음")
        return not bad

    def plot(self, ax, title=""):
        import matplotlib.pyplot as plt

        cmap = plt.get_cmap("tab20")
        for idx, piece in enumerate(self.pieces()):
            xs = [self.P[k][0] for k in piece]
            ys = [self.P[k][1] for k in piece]
            ax.fill(xs, ys, color=cmap(idx % 20), alpha=0.45, linewidth=0)

        xs = [p[0] for p in self.P] + [self.P[0][0]]
        ys = [p[1] for p in self.P] + [self.P[0][1]]
        ax.plot(xs, ys, "k-", linewidth=1.8)

        for order, (i, j) in enumerate(self.diagonals, 1):
            (x1, y1), (x2, y2) = self.P[i], self.P[j]
            ax.plot([x1, x2], [y1, y2], "r--", linewidth=1.5)
            ax.text((x1 + x2) / 2, (y1 + y2) / 2, f"d{order}", color="red", fontsize=8,
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))

        for i, (x, y) in enumerate(self.P):
            color = "red" if self.is_reflex(i) else "black"
            ax.plot(x, y, "o", color=color, markersize=5)
            ax.annotate(str(i), (x, y), textcoords="offset points", xytext=(5, 5), fontsize=9)

        ax.set_title(f"{title}  ({len(self.diagonals)} diagonals)")
        ax.set_aspect("equal")
        ax.grid(alpha=0.3)


def convex_partition(points, verbose=True):
    """다각형을 볼록 분할하고 대각선 리스트를 돌려준다."""
    return ConvexPartition(points, verbose).run()


# ---------------------------------------------------------------------------
# 예제 & 실행
# ---------------------------------------------------------------------------

def _star(k=5, r_out=4.0, r_in=1.6):
    pts = []
    for t in range(2 * k):
        r = r_out if t % 2 == 0 else r_in
        a = math.pi / 2 + t * math.pi / k
        pts.append((round(r * math.cos(a), 4), round(r * math.sin(a), 4)))
    return pts


EXAMPLES = {
    "L-shape": [(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4)],
    "U-shape": [(0, 0), (6, 0), (6, 6), (4, 6), (4, 2), (2, 2), (2, 6), (0, 6)],
    "Comb": [(0, 0), (9, 0), (9, 5), (8, 5), (7, 2), (6, 5), (5, 5),
             (4, 2), (3, 5), (2, 5), (1, 2), (0, 5)],
    "Star": _star(),
    # 대각선 하나로는 꼭짓점 3 을 해결할 수 없는 예
    "Notch": [(-3, 0), (3, 0), (3, 4), (0, 1), (-3, 4)],
}


def random_polygon(n, rng):
    """원점 기준 각도 순으로 점을 이은 무작위 단순 다각형 (반시계 방향)."""
    while True:
        angles = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
        gaps = [angles[i + 1] - angles[i] for i in range(n - 1)]
        gaps.append(angles[0] + 2 * math.pi - angles[-1])
        if max(gaps) < math.pi:  # 원점이 다각형 안에 있어야 단순 다각형이 보장됨
            break
    pts = []
    for a in angles:
        r = rng.uniform(1, 10)
        pts.append((r * math.cos(a), r * math.sin(a)))
    return pts


def random_test(n, count, seed):
    """무작위 다각형 count 개로 알고리즘을 돌려 결과를 요약한다."""
    import random
    rng = random.Random(seed)
    results = []
    failed = 0
    total_reflex = total_diag = total_pieces = total_fallback = 0

    for k in range(count):
        pts = random_polygon(n, rng)
        cp = ConvexPartition(pts, verbose=False)
        cp.run()
        pcs = cp.pieces()
        ok = all(is_convex_polygon([cp.P[v] for v in p]) for p in pcs)
        if not ok:
            failed += 1
            print(f"  실패: {k}번째 다각형 (seed={seed})")
        total_reflex += sum(cp.is_reflex(i) for i in range(cp.n))
        total_diag += len(cp.diagonals)
        total_pieces += len(pcs)
        total_fallback += len(cp.fallback_vertices)
        results.append((f"random #{k}", cp))

    print("=" * 60)
    print(f"무작위 테스트: 꼭짓점 {n}개짜리 다각형 {count}개 (seed={seed})")
    print("=" * 60)
    print(f"  모든 조각이 볼록한 다각형 : {count - failed} / {count}")
    print(f"  평균 reflex 꼭짓점 수     : {total_reflex / count:.1f}")
    print(f"  평균 대각선 수            : {total_diag / count:.1f}")
    print(f"  평균 조각 수              : {total_pieces / count:.1f}")
    print(f"  대각선 하나로 안 된 꼭짓점: 다각형당 평균 {total_fallback / count:.1f}개")
    return results


def compare_with_optimal(n, count, seed):
    """우리 알고리즘과 최적해(optimal_partition.py)의 조각 수를 비교한다."""
    import random
    import time
    from optimal_partition import optimal_partition

    rng = random.Random(seed)
    ours_total = opt_total = same = 0
    worst = (0.0, None)
    t_ours = t_opt = 0.0
    results = []

    for k in range(count):
        pts = random_polygon(n, rng)
        t = time.perf_counter()
        cp = ConvexPartition(pts, verbose=False)
        cp.run()
        t_ours += time.perf_counter() - t

        t = time.perf_counter()
        opt_diags = optimal_partition(pts)
        t_opt += time.perf_counter() - t

        ours, opt = len(cp.diagonals) + 1, len(opt_diags) + 1
        ours_total += ours
        opt_total += opt
        same += ours == opt
        if ours / opt > worst[0]:
            worst = (ours / opt, k)

        if len(results) < 2:
            best = ConvexPartition(pts, verbose=False)
            for i, j in opt_diags:
                best.draw(i, j)
            results.append((f"#{k} ours", cp))
            results.append((f"#{k} optimal", best))

    print("=" * 60)
    print(f"최적해 비교: 꼭짓점 {n}개짜리 다각형 {count}개 (seed={seed})")
    print("=" * 60)
    print(f"  평균 조각 수   우리: {ours_total / count:.2f}   최적: {opt_total / count:.2f}"
          f"   (평균 {ours_total / opt_total:.3f}배)")
    print(f"  최적과 같은 다각형: {same} / {count}")
    print(f"  가장 나쁜 경우: 최적의 {worst[0]:.2f}배 ({worst[1]}번째 다각형)")
    print(f"  평균 실행 시간 우리: {t_ours / count * 1000:.1f} ms   최적: {t_opt / count * 1000:.1f} ms")
    return results


def read_polygon(path):
    pts = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip().replace(",", " ")
            if line and not line.startswith("#"):
                x, y = line.split()[:2]
                pts.append((float(x), float(y)))
    return pts


def main(argv):
    save_path = None
    verbose = True
    files = []
    random_n, compare_n, count, seed = None, None, 100, 0
    i = 0
    while i < len(argv):
        if argv[i] in ("--save", "--random", "--compare", "--count", "--seed"):
            value = argv[i + 1]
            if argv[i] == "--save":
                save_path = value
            elif argv[i] == "--random":
                random_n = int(value)
            elif argv[i] == "--compare":
                compare_n = int(value)
            elif argv[i] == "--count":
                count = int(value)
            else:
                seed = int(value)
            i += 2
            continue
        if argv[i] == "--quiet":
            verbose = False
        else:
            files.append(argv[i])
        i += 1

    if compare_n is not None:
        # 그림은 앞의 2개 다각형 (우리 / 최적)
        results = compare_with_optimal(compare_n, count, seed)
    elif random_n is not None:
        # 그림은 앞의 4개만
        results = random_test(random_n, count, seed)[:4]
    else:
        polygons = {f: read_polygon(f) for f in files} if files else EXAMPLES
        results = []
        for name, pts in polygons.items():
            print("=" * 60)
            print(f"다각형: {name}")
            print("=" * 60)
            cp = ConvexPartition(pts, verbose)
            cp.run()
            cp.report()
            print()
            results.append((name, cp))

    try:
        import matplotlib
        if save_path:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib 이 없어서 그림은 생략합니다.  (설치: pip install matplotlib)")
        return

    cols = min(len(results), 2)
    rows = (len(results) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows), squeeze=False)
    for ax, (name, cp) in zip(axes.flat, results):
        cp.plot(ax, name)
    for ax in list(axes.flat)[len(results):]:
        ax.axis("off")
    fig.tight_layout()

    # 창을 띄울 수 없는 환경(WSL 등)이면 파일로 저장
    if not save_path and matplotlib.get_backend().lower() == "agg":
        save_path = "convex_partition.png"

    if save_path:
        fig.savefig(save_path, dpi=120)
        print(f"그림 저장: {save_path}")
    else:
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
