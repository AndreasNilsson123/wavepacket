#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "wavepacket_bonus.h"

int main(int argc, char **argv)
{
    /* Initialize */
    double t;
    parameters params;
    bool user_output = false;
    char wf_output[50];
    if (argc < 2)
    {
        fprintf(stderr, "Not enough input arguments!\n");
        abort();
    }
    /* Read parameters from input file */
    read_parameters(&params, argv[1], &user_output, wf_output);
    /* Set up grid */
    make_grid(&params);
    /* Header for stdout */
    print_header(params);
    fprintf(stdout, "\nwf_output: %s\n", wf_output);
    
    /* Allocate memory for wavefunction */
    double complex *psi = (double complex *) malloc(params.nx*sizeof(double complex));
    double complex *psii = (double complex *) malloc(params.nx*sizeof(double complex));
    double *pot = (double *) malloc(params.nx*sizeof(double));
    /* Initialize potential, user supplied */
    initialize_potential(params, argc, argv);

    /* Initialize wavefunction */
    initialize_wf(params, argc, argv, psi);

    if(user_output)
        initialize_user_observe(params, argc, argv);

    /* Initialize differention matrix */
    gsl_matrix_complex *T = gsl_matrix_complex_calloc(params.nx,params.nx);
    initialize_differention_matrix(params, T);
    
    /* Initialize potential matrix */
    gsl_matrix_complex *V = gsl_matrix_complex_calloc(params.nx,params.nx);
    // We make use of time independent potentials so calculations are moved outside
    // to fasten up calculations.
    potential(params, 0, pot);
    initialize_potential_matrix(params, V, pot);
    /* Initialize Hamiltonian matrix */
    gsl_matrix_complex *H = gsl_matrix_complex_calloc(params.nx,params.nx);
    gsl_matrix_complex_add(H,T);
    gsl_matrix_complex_add(H,V);
    gsl_matrix_complex_free(T);
    gsl_matrix_complex_free(V);
    fprintf(stdout, "\nSTATUS: \n");
    /* Diagonolize hamiltonian */
    gsl_vector *eval = gsl_vector_calloc(params.nx);
    gsl_matrix_complex *evec = gsl_matrix_complex_calloc(params.nx, params.nx);
    gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(params.nx);
    gsl_eigen_hermv(H, eval, evec, w);
    gsl_eigen_hermv_free(w);
    //gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
    gsl_matrix_complex_free(H);
    fprintf(stdout, "Check 1: Diagonalization complete!\n");

    /* Allocate memory */
    gsl_vector_complex *exp_E = gsl_vector_complex_calloc(params.nx);

    // Retrive values for vector
    time_step_matrix(params, exp_E, eval);
    fprintf(stdout, "Check 2: Time step matrix!\n");
    /* Get start wf */
    get_start_wf(params, psi, evec);
    fprintf(stdout, "Check 3: Start wave function computed!\n");
    /* Time evolution */
    t = 0.; 
    fprintf(stdout, "Check 4: Begin time evolution!\n");
    for(size_t t_index = 0; t_index < params.nt; t_index++)
    {
        /* Computes psi_(t+1) = M*psi(t) */
        time_evo(params, psi, psii, exp_E, t);
        if(t_index%params.nprint==0){
            //fprintf(stdout, "%.2f\n", (double)t_index/(double)params.nt);
            //psii = psi -> we don't need to reconvert psi.
            recover_true_wf(params, psii, evec);
            user_observe(params, t, psii);
        }
        
        t += params.dt;
    }
    fprintf(stdout, "Check 5: Time evolution complete!\n");
    /* Get real wf */
    recover_true_wf(params, psi, evec);
    /* Write last wavefunction to file */
    FILE *fp_out_wf = fopen(wf_output, "wb");
    if(fp_out_wf==NULL)
    {
        fprintf(stderr, "Error: Failed to open %s\n",wf_output);
    }
    renormalize(params,psi);
    for(int i=0; i<params.nx; i++)
    {
        fprintf(fp_out_wf,"%f %13.6e %13.6e %13.6e\n",
        params.x[i],creal(psi[i]),cimag(psi[i]), cabs(psi[i])*cabs(psi[i])); 
    }

    free(psi);
    free(psii);
    free(pot);
    gsl_matrix_complex_free(evec);
    gsl_vector_free(eval);
    return 0;
}


void read_parameters(parameters *const params, char *filename, bool *user_output, char *wf_output)
{   
    extern char *out_file;
    FILE *fp_in = fopen(filename, "r");
    if(fp_in == NULL)
    {
        fprintf(stderr, "Error: Failed to open %s\n",filename);
        abort();
    }
    char key[30], value[80];
    while(fscanf(fp_in, "%s = %s", key, value) != EOF)
    {
        if(!strcmp(key, "mass"))
        {
            params->mass = atof(value);
        }
        else if(!strcmp(key,"nx"))
        {
            params->nx = atoi(value);
            params->nx_local = atoi(value);
        }
        else if(!strcmp(key,"x_min"))
        {
            params->x_min = atof(value);
        }
        else if(!strcmp(key,"x_max"))
        {
            params->x_max = atof(value);
        }
        else if(!strcmp(key,"nt"))
        {
            params->nt = atoi(value);
        }
        else if(!strcmp(key,"dt"))
        {
            params->dt = atof(value);
        }
        else if(!strcmp(key,"nprint"))
        {
            params->nprint = atoi(value);
        }
        else if(!strcmp(key,"user_output"))
        {
            *user_output = atoi(value);
        }
        else if(!strcmp(key,"wf_output"))
        {
            strcpy(wf_output, value);
        }
        else
        {
            fprintf(stderr,"Unknown parameter %s in file: %s\n", key, filename);
        }
    }
    /* Set other params */
    params->hbar = 1;
    if(params->nt < params->nprint)
    {
        fprintf(stderr,"nprint: %d > nt: %d\n", params->nprint, params->nt);
        abort();
    }
    return;
}

void print_header(const parameters params)
{
    fprintf(stdout, "*** PARAMETERS ***\n");
    fprintf(stdout, "x_min: %.2f\n", params.x_min);
    fprintf(stdout, "x_max: %.2f\n", params.x_max);
    fprintf(stdout, "nx: %ld\n", params.nx);
    fprintf(stdout, "dx: %.2f\n", params.dx);
    fprintf(stdout, "k_max: %.1f\n", 3.141592/params.dx);
    fprintf(stdout, "nt: %d\n", params.nt);
    fprintf(stdout, "dt: %.2e\n", params.dt);
}

/* -------------------------------------- */
/* ----------- Calculations ------------- */
/* -------------------------------------- */
void make_grid(parameters *params)
{
    params->x = (double *) malloc(params->nx * sizeof(double));
    params->x2 = (double *) malloc(params->nx * sizeof(double));
    params->dx = (params->x_max-params->x_min)/(double)(params->nx-1);

    for(size_t i=0; i<params->nx ; i++)
    {
        
        params->x[i] = params->x_min + (double)(i)*params->dx;
        params->x2[i] = (params->x[i])*(params->x[i]);
    }
    return;
}

void initialize_differention_matrix(const parameters params, gsl_matrix_complex *T)
{
    for (size_t i=0; i < params.nx-1; i++)
    {
        for(size_t j=0; j < params.nx; j++)
        {
            if(j==i){
                gsl_matrix_complex_set(T, i, j, gsl_complex_rect(M_PI*M_PI/3,0.)); // Set diagonal elements
            } else {
                int exponent = (i-j)%2; 
                gsl_matrix_complex_set(T, i, j, gsl_complex_rect(gsl_pow_int(-1,exponent)*2/gsl_pow_2(i-j),0.));
            }
        }   
    }
    // Scaling factor for T 
    double scale_T = params.hbar*params.hbar/(2*params.mass*params.dx*params.dx);
    // Apply scaling to T
    gsl_matrix_complex_scale(T, gsl_complex_rect(scale_T,0.));
    return;
}

void initialize_potential_matrix(const parameters params, gsl_matrix_complex *V, double *pot)
{
    for (size_t i=0; i < params.nx; i++)
    {
        gsl_matrix_complex_set(V, i, i, gsl_complex_rect(pot[i],0.));
    }
    return;
}

void time_step_matrix(parameters params, gsl_vector_complex *exp_E, gsl_vector *eval)
{
    // Set diagonal exp(...) matrix as a vector 
    for(int i=0; i<params.nx; i++)
    {
        gsl_vector_complex_set(exp_E, i, gsl_complex_polar(1., -params.dt*gsl_vector_get(eval,i)/params.hbar));
    }
}


void time_evo(parameters params, double complex *psi, double complex *psii, gsl_vector_complex *exp_E, double t)
{
    gsl_complex sum;
    for(int i = 0; i<params.nx; i++)
    {
        sum = GSL_COMPLEX_ZERO;

        sum = gsl_complex_add(sum, gsl_complex_mul_real(gsl_vector_complex_get(exp_E, i), creal(psi[i])));
        sum = gsl_complex_add(sum, gsl_complex_mul_imag(gsl_vector_complex_get(exp_E, i), cimag(psi[i])));
        psii[i] = GSL_REAL(sum) + I*GSL_IMAG(sum);
    }
    for(int i = 0; i<params.nx; i++)
    {
        psi[i] = psii[i];
    }

    //renormalize(params,psi);
    return;
}

void get_start_wf(parameters params, complex double *psi, gsl_matrix_complex *evec)
{
    gsl_complex sum;
    double complex *psi_temp = (double complex *) malloc(params.nx*sizeof(double complex));
    for(int i = 0; i<params.nx; i++)
    {
        sum = GSL_COMPLEX_ZERO;
        for(int j = 0; j<params.nx; j++)
        {
            sum = gsl_complex_add(sum, gsl_complex_mul_real(gsl_complex_conjugate(gsl_matrix_complex_get(evec, j, i)), creal(psi[j])));
            sum = gsl_complex_add(sum, gsl_complex_mul_imag(gsl_complex_conjugate(gsl_matrix_complex_get(evec, j, i)), cimag(psi[j])));
        }
        psi_temp[i] = GSL_REAL(sum) + I*GSL_IMAG(sum);
    }
    for(int i = 0; i<params.nx; i++)
    {
        psi[i] = psi_temp[i];
    }

    return;
}

void recover_true_wf(parameters params, complex double *psi, gsl_matrix_complex *evec)
{
    gsl_complex sum;
    double complex *psi_temp = (double complex *) malloc(params.nx*sizeof(double complex));
    for(int i = 0; i<params.nx; i++)
    {
        sum = GSL_COMPLEX_ZERO;
        for(int j = 0; j<params.nx; j++)
        {
            sum = gsl_complex_add(sum, gsl_complex_mul_real(gsl_matrix_complex_get(evec, i, j), creal(psi[j])));
            sum = gsl_complex_add(sum, gsl_complex_mul_imag(gsl_matrix_complex_get(evec, i, j), cimag(psi[j])));
        }
        psi_temp[i] = GSL_REAL(sum) + I*GSL_IMAG(sum);
    }
    for(int i = 0; i<params.nx; i++)
    {
        psi[i] = psi_temp[i];
    }

    return;
}

/* -------------------------------------- */
/* ----- Normalization calculations ----- */
/* -------------------------------------- */
double complex
integrate3D (const parameters params, const double complex * const f1,
	     const double complex * const f2)
{
  /** Integrates f1^* f2 dtau  **/
  
  double complex sum = 0., sum_total;
  
  for (size_t l = 0; l < params.nx; ++l)
    {
      sum += conj (f1[l]) * f2[l];
    }
  
  sum *= params.dx;

  return sum;
}

double
norm (const parameters params, const double complex * const psi)
{
  /** Calculates the norm of the wave function **/

  return sqrt (creal (integrate3D (params, psi, psi)));
}

void
renormalize (const parameters params, double complex *psi)
{
  /** Renormalize the wave function **/

  double fact = 1. / norm (params, psi);
  
  for (size_t l = 0; l < params.nx; ++l)
    psi[l] *= fact;

  return;
}

/* -------------------------------------- */
/* ---------- Utility routines ---------- */
/* -------------------------------------- */
void
abort ()
{
  /** Abort execution **/
  fprintf (stderr, "Execution aborted.\n");
  exit(-99);
}