#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define N 300    // space
#define T 180000    // time
#define IDX(x,t) ((t) * N + (x))  // flat indexing


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// The present program makes tries to simulate the Continuum Noisy Voter Model (CNVM))
// We simulate a system of coalescing random walks and we mark them with a 
// geometric (since we are discrete) clock of paramerer lamnda/N. 
// The results look like the first program I had made, even if here I am 
// rescaling the space diffusively, fully simularing the deiscrete dynamics 
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



//typedef struct { 
//	// The parameters needed to define
//	int N; // scaling parameter
//	double T;
//	int L; 
//	double annihilation;  // annihilation rates of the microscopic model
//	double creation;  // creation rates of the microscopical model
//	// The flipping rates for an interface spin are assumed to be 1; 
//	//  The parameters needed to discretize the dynamics
//
//	int N_simulations; // number of independent simulations
//	int resolution; // Number of subntervals in which to divide the unit of macroscopic time. Since we need N^2 time steps of the microscopic dynamics, we need to plot N*N/resolution micro evolutions
//	int Micro_n_steps; // Number of time steps in which to divide the microscopic dynamics. It should be made in such a way that annihilation/Micro_Time_Step << 1, and we must pay attention to the resolution for the Mersenne Twister, that is creation/Micro_Time_Step must be detectable by a uniformly generated random variable.  
//	
//	Initialization_function initialize; // Pointer to the function that initializes the spins at time 0.
//	void *initialization_param; // The void pointer will hold the third argument of the above function. It will be a double pointer if the initial distribution needs a parameter, a NULL pointer otherwise
//

typedef struct Node {
    struct Node *parents[2]; // backward edges
    struct Node *child;      // forward pointer
    int n_parents;
    int is_marked;
    int color;               // -1 = uninitialized, 0 or 1 = opinion
} Node;




//
//int set_N(Configuration *conf, const char *value);
//
//int set_T(Configuration *conf, const char *value);
//
//int set_L(Configuration *conf, const char *value); 
//
//int set_N_simulations(Configuration *conf, const char *value); 
//
//int set_resolution(Configuration *conf, const char *value); 
//
//int set_Micro_n_steps(Configuration *conf, const char *value); 
//
//typedef struct {
//    const char *key;
//    ConfigSetter setter;
//} ConfigEntry; // This structure will be used to associate to each character which field of conf to fill, that is, which  Config_setter to use 
//
//
//
//
//ConfigEntry config_table[] = {
//    {"N", set_N},
//    {"T", set_T},
//    {"L", set_L},
//    {"N_simulations", set_N_simulations}, 
//    {"resolution", set_resolution},
//    {"Micro_n_steps", set_Micro_n_steps},
//    {NULL, NULL} //A sentinel indicating the end of the structure. 
//}; // I need to validate the variables that I read!!!REMEMBER Ok, let's use the simple suggestion  and write a validating function 
//
Node grid[N * T];


void initialize_tree(double lambda, gsl_rng *rng) {
    // Initialize all nodes
    for (int t = 0; t < T; ++t) {
        for (int x = 0; x < N; ++x) {
            Node *n = &grid[IDX(x, t)];
            n->n_parents = 0;
            n->is_marked = 0;
            n->color = -1;
            n->child = NULL;
        }
    }

    // Each node at time t chooses a child at t+1
    for (int t = 0; t < T - 1; ++t) {
        for (int x = 0; x < N; ++x) {
            Node *parent = &grid[IDX(x, t)];

            int dir = gsl_rng_uniform_int(rng, 2); // 0 = left, 1 = right
            int child_x = wrap(x + (dir == 0 ? -1 : 1));
            Node *child = &grid[IDX(child_x, t + 1)];

            parent->child = child;
            child->parents[child->n_parents++] = parent;

            // Optional: mark nodes at random (except at t = T-1, handled separately)
            parent->is_marked = (gsl_rng_uniform(rng) < lambda / N);
            if (parent->is_marked)
                parent->color = gsl_rng_uniform_int(rng, 2);
        }
    }

    // Final row (leaves): force marking and random color
    for (int x = 0; x < N; ++x) {
        Node *leaf = &grid[IDX(x, T - 1)];
        leaf->is_marked = 1; // force mark
        leaf->color = gsl_rng_uniform_int(rng, 2); // initial opinion
        // no child
    }
}

// New propagate_color function (no recursion)
void propagate_color(Node *n, gsl_rng *rng) {
    if (n->color == -1) return;  // skip if not colored

    for (int i = 0; i < n->n_parents; ++i) {
        Node *p = n->parents[i];
        if (p->color == -1) {
            if (p->is_marked) {
                // Assign random color to marked parents
                p->color = gsl_rng_uniform_int(rng, 2);
            } else {
                // Inherit child's color if not marked
                p->color = n->color;
            }
            // No further propagation here
        }
    }
}

void assign_colors_backwards(gsl_rng *rng) {
    // Since no recursion, we do multiple passes backward through time to propagate colors
    for (int t = T - 1; t > 0; --t) {
        for (int x = 0; x < N; ++x) {
            Node *n = &grid[IDX(x, t)];
            propagate_color(n, rng);
        }
    }

    // Assign fallback color for uncolored nodes
    for (int i = 0; i < N * T; ++i) {
        if (grid[i].color == -1)
            grid[i].color = 0;
    }
}


// Output color field to .pgm grayscale image
void write_colors_pgm(const char *filename) {
    FILE *f = fopen(filename, "wb");
    if (!f) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "P5\n%d %d\n255\n", N, T);
    for (int t = 0; t < T; ++t) {
        for (int x = 0; x < N; ++x) {
            unsigned char pixel = (unsigned char)(grid[IDX(x, t)].color * 255);
            fwrite(&pixel, sizeof(unsigned char), 1, f);
        }
    }

    fclose(f);
}


void write_colors_txt(const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    for (int t = 0; t < T; ++t) {
        for (int x = 0; x < N; ++x) {
            fprintf(f, " %d ", grid[IDX(x, t)].color);
            if (x < N - 1)
                fputc(',', f);  // comma between values
        }
        fputc('\n', f);  // newline after each time slice
    }

    fclose(f);
}

void write_colors_bin(const char *filename) {
    FILE *f = fopen(filename, "wb");
    if (!f) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    for (int t = 0; t < T; ++t) {
        for (int x = 0; x < N; ++x) {
            unsigned char c = (unsigned char)grid[IDX(x, t)].color;
            fwrite(&c, sizeof(unsigned char), 1, f);
        }
    }

    fclose(f);
}

int main() {
    const double lambda = 1.0/N; 

    // Init GSL RNG
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    initialize_tree(lambda, rng);
	assign_colors_backwards(rng);
//	write_colors_pgm("output.pgm");
	write_colors_bin("output.bin");
	write_colors_txt("output.txt");
    gsl_rng_free(rng);
    printf("Simulation complete. Output written to output.pgm\n");
    return 0;
}

int set_N(Configuration *conf, const char *value) {
      printf("set_T called with value='%s'\n", value);
    int v = atoi(value);
    if (v <= 0) return -1;
    printf(" The value of N: %d", v); 
    conf->N = v;
    printf(" The value of N: %d", conf -> N); 
    return 0;
}

int set_T(Configuration *conf, const char *value) {
    float v = atof(value);
    if (v <= 0.0f) return -1;
    conf->T = v;
    return 0;
}

int set_L(Configuration *conf, const char *value) {
    float v = atof(value);
    if (v <= 0.0f) return -1;
    conf->L = v;
    return 0;
}

int set_N_simulations(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->N_simulations = v;
    return 0;
}

int set_resolution(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->resolution = v;
    return 0;
}

int set_Micro_n_steps(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->Micro_n_steps = v;
    return 0;
