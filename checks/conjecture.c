/* Exhaustive check of Conjecture 5.5 for all rankings sigma of n candidates:
   (t-1)|E_t| <= F(sigma) for every real t >= 1, where E_t is the set of pairs a < b with
   sigma(a) - sigma(b) >= t(b - a) and F is the footrule distance to the identity.
   |E_t| is a step function of t and (t-1)|E_t| increases between its steps, so it suffices to test
   t = h/g for the inverted pairs (h = sigma(a) - sigma(b), g = b - a). All arithmetic is in integers.
   Prints the largest value of (t-1)|E_t|/F as a fraction.
   Compile: cc -O2 -o conjecture conjecture.c      Run: ./conjecture 10            */
#include <stdio.h>
#include <stdlib.h>
static int nextperm(int *a, int m) {
  int i = m - 2; while (i >= 0 && a[i] >= a[i + 1]) i--; if (i < 0) return 0;
  int j = m - 1; while (a[j] <= a[i]) j--; int t = a[i]; a[i] = a[j]; a[j] = t;
  for (int l = i + 1, r = m - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
  return 1;
}
int main(int argc, char **argv) {
  int n = atoi(argv[1]), s[16];
  for (int i = 0; i < n; i++) s[i] = i;
  long long bestnum = 0, bestden = 1; long long count = 0;
  do {
    int F = 0; for (int c = 0; c < n; c++) F += abs(s[c] - c);
    if (F == 0) continue;
    count++;
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) {
      int h = s[a] - s[b], g = b - a;
      if (h <= g) continue;                         /* t = h/g must exceed 1 */
      int f = 0;
      for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++)
        if ((long long)(s[x] - s[y]) * g >= (long long)h * (y - x)) f++;
      long long num = (long long)(h - g) * f, den = (long long)g * F;   /* (t-1)|E_t| / F */
      if (num > den) { printf("COUNTEREXAMPLE n=%d t=%d/%d\n", n, h, g); return 1; }
      if (num * bestden > bestnum * den) { bestnum = num; bestden = den; }
    }
  } while (nextperm(s, n));
  long long a1 = bestnum, b1 = bestden; while (b1) { long long t = a1 % b1; a1 = b1; b1 = t; }
  if (a1) { bestnum /= a1; bestden /= a1; }
  printf("n=%d: %lld rankings checked, largest (t-1)|E_t|/F = %lld/%lld\n", n, count, bestnum, bestden);
  return 0;
}
