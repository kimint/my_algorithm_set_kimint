#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct { int *b, *e; } Seg;

static int bit_width(const int *arr, int n)
{
    int mx = 0;
    for (const int *p = arr, *end = arr + n; p != end; ++p)
        if (*p > mx) mx = *p;
    int k = 0;
    while (mx >> k) ++k;
    return k;
}

static int *stable_part(int *b, int *e, int bit, int *buf)
{
    int  n  = (int)(e - b);   /* segment length: pointer subtraction  */
    int  c1 = 0;
    int *p;

    /* Pass 1: count set bits — pure pointer walk, no index */
    for (p = b; p != e; ++p)
        c1 += (*p >> bit) & 1;

    /* Pass 2: fill buf: 1-partition at front, 0-partition after */
    int *q1 = buf, *q0 = buf + c1;
    for (p = b; p != e; ++p) {
        if ((*p >> bit) & 1) *q1++ = *p;
        else                  *q0++ = *p;
    }

    memcpy(b, buf, (size_t)n * sizeof(int));
    return b + c1;   /* split point as pointer — no index conversion needed */
}

void bit_sort_desc(int *arr, int n)
{
    if (n <= 1) return;

    int   k   = bit_width(arr, n);
    int  *buf = malloc((size_t) n      * sizeof(int));  /* partition scratch */
    Seg  *cur = malloc((size_t)(n + 1) * sizeof(Seg));  /* active segments   */
    Seg  *nxt = malloc((size_t)(n + 1) * sizeof(Seg));  /* next segments     */
    int   cn  = 0, nn = 0;

    cur[cn++] = (Seg){arr, arr + n};   /* start: one segment = entire array  */

    for (int bit = k - 1; bit >= 0; --bit) {
        nn = 0;
        for (int i = 0; i < cn; ++i) {
            int *b = cur[i].b, *e = cur[i].e;
            if (e - b <= 1) continue;               /* singleton: done      */

            int *sp = stable_part(b, e, bit, buf);  /* partition in-place   */

            if (sp != b && sp != e) {               /* non-trivial split    */
                if (sp - b > 1) nxt[nn++] = (Seg){b,  sp};
                if (e  - sp > 1) nxt[nn++] = (Seg){sp, e};
            } else {                                /* no split this round  */
                if (e - b > 1)  nxt[nn++] = (Seg){b, e};
            }
        }

        /* O(1) swap: swap pointers instead of copying segment data */
        Seg *tmp = cur; cur = nxt; nxt = tmp;
        cn = nn;
    }

    free(buf); free(cur); free(nxt);
}

int main(void)
{
    int n;
    scanf("%d", &n);

    int *arr = malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; ++i)
        scanf("%d", &arr[i]);

    clock_t t0 = clock();
    bit_sort_desc(arr, n);
    clock_t t1 = clock();

    for (int i = 0; i < n; ++i)
        printf("%d ", arr[i]);
    printf("\n");
    printf("Time: %.6f s\n", (double)(t1 - t0) / CLOCKS_PER_SEC);

    free(arr);
    return 0;
}