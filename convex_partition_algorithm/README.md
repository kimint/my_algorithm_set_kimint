# Convex Partition Algorithm

한국어 | [English](README.en.md)

단순 다각형(simple polygon)을 **대각선만 그어서**(새 점 추가 없이) 볼록 다각형들로 나누는 스택 기반 알고리즘입니다.
최적해를 구하는 동적계획법도 함께 들어 있어서 결과를 비교할 수 있습니다.

## 파일

| 파일 | 내용 |
|---|---|
| `convex_partition.py` | 스택 기반 볼록 분할 알고리즘, 예제, 무작위 테스트, 그림 |
| `optimal_partition.py` | 최소 조각 수를 구하는 동적계획법 (비교용) |
| `requirements.txt` | 필요한 패키지 (`matplotlib`, 그림용) |

## 설치 및 실행

```bash
python3 -m venv myenv
source myenv/bin/activate
pip install -r requirements.txt
```

```bash
python3 convex_partition.py                  # 내장 예제 실행 (L-shape, U-shape, Comb, Star, Notch)
python3 convex_partition.py polygon.txt      # 파일에서 다각형 읽기 (한 줄에 "x y", 반시계 방향)
python3 convex_partition.py --random 50      # 꼭짓점 50개짜리 무작위 다각형 100개로 테스트
python3 convex_partition.py --compare 50     # 무작위 다각형으로 최적해와 조각 수 비교
```

옵션: `--quiet`(단계별 로그 끄기), `--save out.png`(그림 저장), `--count 개수`, `--seed 시드`.
그림 창을 띄울 수 없는 환경(WSL 등)에서는 `convex_partition.png` 로 자동 저장됩니다.

## 알고리즘

입력: 반시계 방향(CCW) 순서의 꼭짓점 리스트

**Step 1~4: stack 만들기**
1. 다각형을 반시계 방향으로 돌며 첫 reflex vertex 를 찾아 `vt` 로 두고 stack 에 push 한다.
2. `vt` 부터 반시계 방향으로 돌면서
   - `vt` 와 현재 점 `w` 를 잇는 선이 다각형의 대각선이 아니면, `w` 의 이전 점을 `vt` 로 바꾸고 push
   - `w` 가 reflex vertex 이면 `vt = w` 로 바꾸고 push
   - 한 바퀴를 다 돌면 다음 단계로

**Step 5~8: stack / stack2 처리**
- stack 에서 pop 한 값을 stack2 에 넣고, stack 맨 위 `a` 와 stack2 맨 위 `b` 를 비교한다.
  - `b` 가 이미 convex → `b` 를 버림
  - 대각선 `a-b` 를 그어서 `b` 가 convex 가 됨 → 긋고 `b` 를 버림
  - 대각선 `a-b` 를 그어도 `b` 가 reflex → 긋고 `b` 를 stack2 에 남김
  - 그 다음 `a` 를 stack2 로 옮긴다. stack 이 빌 때까지 반복.

**Step 9~13: 남은 reflex vertex 처리**
- stack2 에서 하나씩 pop 해서, 아직 reflex 이면
  1. 그 꼭짓점의 가장 큰 각을 반으로 나누는 **이등분선**을 쏘고, 경계와 만나는 곳에서
     오른쪽 / 왼쪽 꼭짓점을 번갈아 확인하며 꼭짓점을 convex 로 만드는 첫 대각선을 긋는다.
  2. 대각선 하나로 안 되면, convex 가 될 때까지 가장 큰 각을 가장 많이 줄이는 대각선을 하나씩 추가한다.

(예전 방식인 "반시계 방향으로 돌며 찾기"는 `_phase3` 안에 주석으로 남겨 두었습니다.)

### 구현 세부

- **Visibility**: 꼭짓점마다 보이는 꼭짓점 집합을 회전 스윕으로 계산(필요할 때 한 번)해서 대각선 판정을 O(1) 로 한다.
- **Face 관리**: 대각선을 그을 때마다 조각을 나눠 기록해서, "기존 대각선과 교차하는가"를 "두 점이 같은 조각에 있는가"로 판정한다.
- `ConvexPartition(points, fast=False)` 로 예전 방식(매번 모든 변과 교차 검사)을 쓸 수 있다.

## 결과

무작위 다각형(seed=0, 각 100개)에서 최적해와 비교:

| 꼭짓점 수 | 평균 조각 수 | 최적 | 비율 | 가장 나쁜 경우 |
|---|---|---|---|---|
| 20 | 9.80 | 8.00 | 1.23배 | 1.67배 |
| 50 | 28.17 | 22.82 | 1.23배 | 1.41배 |

- 모든 테스트에서 모든 조각이 볼록하게 나뉜다.
- 최적보다 평균 약 23% 많은 조각을 만든다. (Hertel–Mehlhorn 과 비슷한 수준)

## 시간복잡도

| | 시간 |
|---|---|
| 이 알고리즘 (visibility 사용) | 최악 약 O(n² log n) |
| 예전 방식 (`fast=False`) | 최악 O(n³) |
| Hertel–Mehlhorn | O(n log n) |
| 최적해 (Keil–Snoeyink) | O(n + r² min(r², n)) |

n: 꼭짓점 수, r: reflex vertex 수

## 최적해 (`optimal_partition.py`)

Keil 계열의 동적계획법입니다. 현 (i, j) 로 잘린 부분 다각형마다 최소 조각 수를 구하고,
reflex vertex 에 닿는 대각선만 후보로 사용합니다. 꼭짓점 4~10개짜리 다각형 400개에서
모든 조합을 확인하는 브루트포스와 결과가 일치하는 것을 확인했습니다.
