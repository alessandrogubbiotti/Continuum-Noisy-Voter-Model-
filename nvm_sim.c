/*
 * nvm_sim.c -- the noisy voter model on a ring, in the diffusive scaling.
 *
 * Every site has one of two colours, and time is discrete. At each step every
 * site, independently of the others,
 *
 *     with probability lambda/N^2    is marked: it takes a new colour, 0 or 1
 *                                    with probability 1/2;
 *     otherwise                      takes the colour that one of its two
 *                                    neighbours, the left or the right one with
 *                                    probability 1/2, had at the step before.
 *
 * So the colour of a site is found by following its genealogy backwards in
 * time, a random walk that at each step goes to the left or to the right, until
 * the first mark is met. The genealogies of two sites are independent until
 * they meet, and then they are the same: they are coalescing random walks,
 * marked with a geometric clock of parameter lambda/N^2. At time 0 all the
 * sites are marked, so the colours are independent.
 *
 * Sites are spaced by 1/N and a unit of time is N^2 steps. As N grows the
 * genealogies converge to coalescing Brownian motions, the Brownian web, and
 * the marks to a Poisson process of rate lambda along them: the limit is the
 * Continuum Noisy Voter Model of Fontes, Isopi, Newman and Ravishankar,
 * https://arxiv.org/abs/math/0404505. In the limit, at a large time, on the
 * line, the correlation of the colours of two sites at distance d is
 * exp(-sqrt(2 lambda) d), and the interfaces between the two colours have
 * density sqrt(lambda/2). On the ring of length L the density is
 * sqrt(lambda/2) tanh(L sqrt(lambda/2)): it is the one of the line when the
 * ring is long with respect to 1/sqrt(2 lambda).
 *
 * The interfaces are, in the limit, the ones of the Glauber dynamics with
 * intensity nu = sqrt(lambda/2), with the time of the Glauber dynamics equal
 * to half of the time of this program.
 *
 * A genealogy moves by one site at each step, so the genealogy of the site x
 * at the step k stays on the sites y, at the step j, with y + j of the parity
 * of x + k. The two sublattices never meet and carry two independent copies of
 * the model. Only one of them is simulated: at the step k the sites x with
 * x + k even. A row has N*L/2 sites, spaced by 2/N.
 *
 * Compile and run (the GSL is needed):
 *
 *     cc nvm_sim.c -o nvm_sim -lgsl -lgslcblas -lm
 *     ./nvm_sim              writes output.txt and output.pgm
 *     python3 plot.py        writes output.png
 *
 * The parameters are the #define below. They can also be given at compile
 * time, for instance
 *
 *     cc -DN=150 -DL=4 -DLAMBDA=2 nvm_sim.c -o nvm_sim -lgsl -lgslcblas -lm
 *
 * N, L and RESOLUTION are integers, T and LAMBDA are real numbers.
 *
 * output.txt has a first line with the parameters and then one line for each
 * saved time, RESOLUTION of them in a unit of time, from the time 0: the
 * colours of the sites of the row.
 * output.pgm is the same as a picture, with the time going upwards.
 */

/* A parameter given with -D is a macro also while the headers are read, and
   the header of the GSL has a variable called T: the one-letter parameters
   are put aside here and taken back after the headers. */
#pragma push_macro("N")
#pragma push_macro("L")
#pragma push_macro("T")
#undef N
#undef L
#undef T
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#pragma pop_macro("T")
#pragma pop_macro("L")
#pragma pop_macro("N")

/* Parameters */
#ifndef N
#define N 300           /* Sites in a unit of length */
#endif
#ifndef L
#define L 1             /* Length of the ring */
#endif
#ifndef T
#define T 2             /* Final time: N^2 steps for each unit */
#endif
#ifndef LAMBDA
#define LAMBDA 1.0      /* Rate of the marks along a genealogy */
#endif
#ifndef RESOLUTION
#define RESOLUTION 100  /* Saved rows in a unit of time */
#endif
#ifndef SEED
#define SEED 0          /* Random seed; with 0 it is taken from the clock */
#endif

/* The colours of the saved rows, one byte for each site. */
static void write_colors_txt(const char *filename, const unsigned char *rows, int n_rows, int n,
                             unsigned long seed) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    fprintf(f, "# N=%d L=%d T=%.15g lambda=%.15g resolution=%d seed=%lu\n",
            (int)(N), (int)(L), (double)(T), (double)(LAMBDA), (int)(RESOLUTION), seed);
    for (int t = 0; t < n_rows; ++t) {
        for (int x = 0; x < n; ++x) {
            fputc('0' + rows[(size_t)t * n + x], f);
            fputc(x < n - 1 ? ' ' : '\n', f);
        }
    }

    int failed = ferror(f);
    if (fclose(f) != 0) failed = 1;
    if (failed) {
        fprintf(stderr, "%s could not be written completely\n", filename);
        exit(EXIT_FAILURE);
    }
}

/* Output color field to .pgm grayscale image. The last time is the first line
   of the picture, so that the time goes upwards. */
static void write_colors_pgm(const char *filename, const unsigned char *rows, int n_rows, int n) {
    FILE *f = fopen(filename, "wb");
    if (!f) {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    fprintf(f, "P5\n%d %d\n255\n", n, n_rows);
    for (int t = n_rows - 1; t >= 0; --t) {
        for (int x = 0; x < n; ++x) {
            fputc(rows[(size_t)t * n + x] * 255, f);
        }
    }

    int failed = ferror(f);
    if (fclose(f) != 0) failed = 1;
    if (failed) {
        fprintf(stderr, "%s could not be written completely\n", filename);
        exit(EXIT_FAILURE);
    }
}

/* The step at which the row number j is saved: the first even step not before
   j N^2/RESOLUTION. The steps are even so that all the saved rows are on the
   same sites. */
static long long saving_step(long long j, long long steps_per_unit, int resolution) {
    return 2 * ((j * steps_per_unit + 2LL * resolution - 1) / (2LL * resolution));
}

/* The number of interfaces of a row: the pairs of neighbours of different colour. */
static int count_interfaces(const unsigned char *row, int n) {
    int count = 0;
    for (int x = 0; x < n; ++x) {
        count += (row[x] != row[x + 1 == n ? 0 : x + 1]);
    }
    return count;
}

int main(void) {
    const int big_n = (N), length = (L), resolution = (RESOLUTION);
    const double lambda = (double)(LAMBDA), final_time = (double)(T);

    /* The checks of the parameters. */
    if (big_n < 2 || length < 1 || resolution < 1 || !(final_time > 0.0) || !(lambda >= 0.0)) {
        fprintf(stderr, "N must be at least 2, L and RESOLUTION at least 1, T positive, LAMBDA not negative\n");
        return 1;
    }
    const long long steps_per_unit = (long long)big_n * big_n;     /* N^2 */
    if (((long long)big_n * length) % 2 != 0 || (long long)big_n * length < 4) {
        fprintf(stderr, "N*L must be even and at least 4: the walks live on a sublattice of the ring\n");
        return 1;
    }
    if (2LL * resolution > steps_per_unit) {
        fprintf(stderr, "RESOLUTION must be at most N^2/2\n");
        return 1;
    }
    const double mark = lambda / (double)steps_per_unit;           /* lambda/N^2 */
    if (mark > 1.0) {
        fprintf(stderr, "LAMBDA must be at most N^2\n");
        return 1;
    }

    const int n = (int)((long long)big_n * length / 2);            /* sites of a row */
    const long long steps = (long long)floor(final_time * (double)steps_per_unit + 0.5);
    int n_rows = 0;                                                 /* the time 0 is saved too */
    while (saving_step(n_rows, steps_per_unit, resolution) <= steps) n_rows += 1;

    /* Init GSL RNG */
    unsigned long seed = (unsigned long)(SEED);
    if (seed == 0) seed = (unsigned long)time(NULL);
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    unsigned char *old = malloc((size_t)n);
    unsigned char *new = malloc((size_t)n);
    unsigned char *rows = malloc((size_t)n_rows * n);
    if (!old || !new || !rows) {
        perror("malloc");
        return 1;
    }

    /* Time 0: all the sites are marked, the colours are independent. */
    for (int x = 0; x < n; ++x) {
        old[x] = (unsigned char)gsl_rng_uniform_int(rng, 2);
        rows[x] = old[x];
    }
    int saved = 1;
    long long next = saving_step(saved, steps_per_unit, resolution);
    double interfaces = 0.0;    /* summed over the saved rows of the second half of the run */
    int rows_counted = 0;

    /* At the step k the row has the sites 2x + (k mod 2), x = 0, ..., n-1, of
       the ring of N*L sites. The neighbours of the site x of an odd row are
       the sites x and x+1 of the row before; the ones of the site x of an
       even row are x-1 and x. */
    for (long long k = 1; k <= steps; ++k) {
        const int odd = (int)(k & 1);
        for (int x = 0; x < n; ++x) {
            if (gsl_rng_uniform(rng) < mark) {
                /* A mark: a new colour. */
                new[x] = (unsigned char)gsl_rng_uniform_int(rng, 2);
            } else {
                /* The colour of a neighbour: 0 = left, 1 = right. */
                const int left  = odd ? x : (x == 0 ? n - 1 : x - 1);
                const int right = odd ? (x + 1 == n ? 0 : x + 1) : x;
                new[x] = old[gsl_rng_uniform_int(rng, 2) == 0 ? left : right];
            }
        }
        unsigned char *swap = old;
        old = new;
        new = swap;

        if (k == next) {
            for (int x = 0; x < n; ++x) rows[(size_t)saved * n + x] = old[x];
            if (2 * k >= steps) {
                interfaces += count_interfaces(old, n);
                rows_counted += 1;
            }
            saved += 1;
            next = saving_step(saved, steps_per_unit, resolution);
        }
    }

    write_colors_txt("output.txt", rows, saved, n, seed);
    write_colors_pgm("output.pgm", rows, saved, n);

    printf("Simulation complete. Output written to output.txt and output.pgm\n");
    printf("N = %d, L = %d, T = %.15g, lambda = %.15g, seed = %lu: %d rows of %d sites\n",
           big_n, length, final_time, lambda, seed, saved, n);
    if (rows_counted > 0) {
        /* The limit at a large time: sqrt(lambda/2) on the line, and
           sqrt(lambda/2) tanh(L sqrt(lambda/2)) on the ring of length L. */
        const double on_the_line = sqrt(lambda / 2.0);
        printf("Interfaces in a unit of length, averaged over the second half of the run: %.4f"
               " (in the limit %.4f on this ring, sqrt(lambda/2) = %.4f on the line)\n",
               interfaces / rows_counted / length, on_the_line * tanh(length * on_the_line),
               on_the_line);
    }

    gsl_rng_free(rng);
    free(old);
    free(new);
    free(rows);
    return 0;
}
