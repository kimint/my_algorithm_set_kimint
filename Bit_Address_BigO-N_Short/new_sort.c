#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef struct { uint64_t *b, *e; } Seg;

/* double -> 순서 보존 unsigned 키. IEEE-754: 양수는 비트패턴이 이미
   unsigned 비교 순서와 일치. 음수만 전부 반전(부호비트는 항상 반전). */
static uint64_t double_to_key(double d)
{
    uint64_t u;
    memcpy(&u, &d, sizeof(u));
    uint64_t mask = (uint64_t)(-(int64_t)(u >> 63)) | 0x8000000000000000ULL;
    return u ^ mask;
}

static double key_to_double(uint64_t u)
{
    uint64_t mask = ((u >> 63) - 1) | 0x8000000000000000ULL;
    u ^= mask;
    double d;
    memcpy(&d, &u, sizeof(d));
    return d;
}

/* int32 -> 순서 보존 unsigned 키. 2의 보수는 이미 양수/음수 구간이
   순서대로 배치돼 있어서, 부호비트 하나만 반전하면 그걸로 충분
   (double처럼 조건부 전체 반전이 필요 없음). */
static uint64_t int_to_key(int32_t v)
{
    return (uint32_t)v ^ 0x80000000u;
}

static int32_t key_to_int(uint64_t u)
{
    return (int32_t)((uint32_t)u ^ 0x80000000u);
}

static uint64_t *stable_part(uint64_t *b, uint64_t *e, int bit, uint64_t *buf)
{
    int n = (int)(e - b);
    int c1 = 0;
    uint64_t *p;

    for (p = b; p != e; ++p)
        c1 += (int)((*p >> bit) & 1);

    uint64_t *q1 = buf, *q0 = buf + c1;
    for (p = b; p != e; ++p) {
        if ((*p >> bit) & 1) *q1++ = *p;
        else                  *q0++ = *p;
    }

    memcpy(b, buf, (size_t)n * sizeof(uint64_t));
    return b + c1;
}

/* k: 살펴볼 비트 수 (double=64, int32=32로 호출). 데이터 값과 무관한
   상수이므로 bit_width() 같은 계산이 없고, 진짜 O(n)이 됨. */
static void bit_sort_desc_keys(uint64_t *arr, int n, int k)
{
    if (n <= 1) return;

    uint64_t *buf = malloc((size_t) n      * sizeof(uint64_t));
    Seg      *cur = malloc((size_t)(n + 1) * sizeof(Seg));
    Seg      *nxt = malloc((size_t)(n + 1) * sizeof(Seg));
    int cn = 0, nn = 0;

    cur[cn++] = (Seg){arr, arr + n};

    for (int bit = k - 1; bit >= 0; --bit) {
        nn = 0;
        for (int i = 0; i < cn; ++i) {
            uint64_t *b = cur[i].b, *e = cur[i].e;
            if (e - b <= 1) continue;

            uint64_t *sp = stable_part(b, e, bit, buf);

            if (sp != b && sp != e) {
                if (sp - b > 1) nxt[nn++] = (Seg){b,  sp};
                if (e  - sp > 1) nxt[nn++] = (Seg){sp, e};
            } else {
                if (e - b > 1) nxt[nn++] = (Seg){b, e};
            }
        }
        Seg *tmp = cur; cur = nxt; nxt = tmp;
        cn = nn;
    }

    free(buf); free(cur); free(nxt);
}

void bit_sort_desc_double(double *arr, int n)
{
    uint64_t *keys = malloc((size_t)n * sizeof(uint64_t));
    for (int i = 0; i < n; ++i) keys[i] = double_to_key(arr[i]);

    bit_sort_desc_keys(keys, n, 64);

    for (int i = 0; i < n; ++i) arr[i] = key_to_double(keys[i]);
    free(keys);
}

void bit_sort_desc_int(int32_t *arr, int n)
{
    uint64_t *keys = malloc((size_t)n * sizeof(uint64_t));
    for (int i = 0; i < n; ++i) keys[i] = int_to_key(arr[i]);

    bit_sort_desc_keys(keys, n, 32);

    for (int i = 0; i < n; ++i) arr[i] = key_to_int(keys[i]);
    free(keys);
}

int main(void)
{
    int mode;   /* 32 = int32_t 정렬, 64 = double 정렬 */
    scanf("%d", &mode);

    if (mode != 32 && mode != 64) {
        fprintf(stderr, "mode must be 32 or 64\n");
        return 1;
    }

    int n;
    scanf("%d", &n);

    clock_t t0, t1;

    if (mode == 32) {
        int32_t *arr = malloc((size_t)n * sizeof(int32_t));
        for (int i = 0; i < n; ++i)
            scanf("%d", &arr[i]);

        t0 = clock();
        bit_sort_desc_int(arr, n);
        t1 = clock();

        for (int i = 0; i < n; ++i)
            printf("%d ", arr[i]);
        printf("\n");

        free(arr);
    } else {
        double *arr = malloc((size_t)n * sizeof(double));
        for (int i = 0; i < n; ++i)
            scanf("%lf", &arr[i]);

        t0 = clock();
        bit_sort_desc_double(arr, n);
        t1 = clock();

        for (int i = 0; i < n; ++i)
            printf("%g ", arr[i]);
        printf("\n");

        free(arr);
    }

    printf("Time: %.6f s\n", (double)(t1 - t0) / CLOCKS_PER_SEC);
    return 0;
}
