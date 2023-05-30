#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavepacket_bonus.h"

// SIMPLE FIX SINCE M_PI WAS NOT RECOGINIZED
#define  M_PI  3.1415926535897932384626433
//Shared variables
double V0, a, sigma, nt, x0;
double k0;
double *V;
bool comp_analytical_wf, comp_RT;
FILE *user_fp;

void read_pot_parameters(const parameters parmas, FILE *fp_in)
{
    bool isRead = false;
    extern bool comp_analytical_wf, comp_RT;
    extern double V0,a,sigma,nt,x0;
    char key[30], value[80];
    
    /* Open file */
    //FILE *fp_in = fopen("val.in","r");
    if (!fp_in){
        fprintf(stderr,"Could not open file!\n");
        abort();
    }
    while(fscanf(fp_in, "%s = %s", key, value) != EOF)
    {
        if(!strcmp(key,"V0"))
        {
            V0 = atof(value);
        }
        else if(!strcmp(key,"a"))
        {
            a = atof(value);
        }
        else if(!strcmp(key,"sigma"))
        {
            sigma = atof(value);
        }
        else if(!strcmp(key,"nt"))
        {
            nt = atof(value);
        }
        else if(!strcmp(key,"x0"))
        {
            x0 = atof(value);
        }
        else if(!strcmp(key,"comp_analytical_wf"))
        {
            comp_analytical_wf = atoi(value);
        }
        else if(!strcmp(key,"comp_RT"))
        {
            comp_RT = atoi(value);
        }
        else
        {
            fprintf(stderr,"Could not find parameter %s\n", key);
        }         
    }
    fclose(fp_in);
}

void initialize_potential (const parameters params, const int argc, char ** const argv)
{
    /* Read file from in arguments */
    FILE *fp_param = fopen(argv[2], "r");
    if(fp_param == NULL)
    {
        fprintf(stderr,"Could not open file %s\n", argv[2]);
        abort();
    }
    /* Initialize parameters needed for potential */
    read_pot_parameters(params, fp_param);
    extern double *V;
    V = (double *) calloc(params.nx_local, sizeof(double));
    extern double a,V0;
    for(size_t i=0; i<params.nx_local; i++)
    {
        double x_val = params.x_min + i*params.dx;
        if(fabs(x_val) < a)
        {
            V[i] = V0;
        } 
        else 
        {
            V[i] = 0.;
        } 
    }
    return;
}

void potential (const parameters params, const double t, double * const pot)
{
    /* Calculate potential at all points */
    extern double *V;
    for(size_t i=0; i<params.nx_local; i++)
    {
        pot[i] = V[i];
    }
    return;
}


void initialize_wf (const parameters params, const int argc, char ** const argv, 
		    double complex *psi)
{
    extern double sigma, x0;

    size_t index;
    double psi_val, xi;   
    // Get k0 from argument
    double k0_wf = atof(argv[3]);
    for(size_t i=0; i<params.nx_local; i++)
    {
        xi = params.x_min + i*params.dx;
        psi[i] = sqrt(sqrt(1./(M_PI*pow(sigma,2))))
                *cexp(I*k0_wf*xi)
                *exp(-pow(xi-x0,2)/(2*pow(sigma,2)));
	
    }
    renormalize(params, psi);
}

void initialize_user_observe (const parameters params, const int argc, char ** const argv)
{
    extern double k0;
    k0 = atof(argv[3]);
    extern bool comp_analytical_wf, comp_RT;
    if(comp_analytical_wf==true)
    {
        user_fp = fopen("data/wf_exakt.dat","w");
        if(user_fp==NULL)
        {
            fprintf(stderr,"Could not open file wf_exakt.dat\n");
            abort();
        }
        printf("\noutput_file: %s\n", "data/wf_exakt.dat");
        fprintf(user_fp, "# x psi_exact_r psi_exact_i psi_exact_2\n");
    }
    else if(comp_RT==true) // Compute transmisson coefficients
    {
        static char fname[50];
        char *filename = fname;
        filename += sprintf(filename, "%s","data/xxx_"); //Change when switching to barrier
        filename += sprintf(filename,"%.0f", k0);
        filename += sprintf(filename,"%s", ".out");
        printf("\noutput_file: %s\n", fname);
        user_fp = fopen(fname,"w");
        if(user_fp==NULL)
        {
            fprintf(stderr,"Could not open file %s\n", fname);
            abort();
        }
        fprintf(user_fp, "# t_RT R T\n");
    
    }
    return;
}

void user_observe (const parameters params, const double t, 
		   const double complex * const psi)
{
    extern bool comp_analytical_wf, comp_RT;
    extern double nt; // External parameter for number of time steps
    extern double a;
    if(comp_analytical_wf==true){ // t >= params.dt*nt-params.dt && 
        /* Compute exakt psi */
        double complex psi_exakt[params.nx_local]; // Psi_exakt
        extern double sigma; // External parameter sigma0
        extern double k0; //External parameter k0
        extern double x0;
        double x_i = x0; // Value for x0

        /* Loop variables */
        double x_val;
        double complex psi_val;
        
        for(size_t i = 0; i < params.nx_local; i++)
        {
            x_val = params.x_min + i*params.dx; 
            double complex c = 2.*sigma*sigma + 2.*t*I; // Reaccuring value 
            double complex exponent_1 = I * (2. * sigma * sigma * (k0 * x_val - k0 * k0 * t / 2.)) / c;
            double complex exponent_2 = -((x_val - x_i) * (+x_val - x_i) + 2. * k0 * x_i * t) / c;
            double complex denominator = csqrt(1.0 / c);
            double complex numerator = cpow(4 * sigma * sigma / M_PI, 0.25);
            psi_exakt[i] = numerator * denominator * cexp(exponent_1) * cexp(exponent_2);
        }
        renormalize(params, psi_exakt);
        for(size_t i = 0; i < params.nx_local; i++)
        {
            x_val = params.x_min + i*params.dx;
            double psi2 = cabs(psi_exakt[i]);
            psi2 *= psi2;
            fprintf(user_fp, "%13.6e %13.6e %13.6e %13.6e\n", x_val, creal(psi_exakt[i]),cimag(psi_exakt[i]), psi2);     
        }
    } 
    else if(comp_RT==true)// Compute transmission coefficients 
    {
        double sum_R, sum_T;
        sum_R = 0.5*cabs(psi[0])*cabs(psi[0]);
        sum_T = 0.5*cabs(psi[params.nx_local-1])*cabs(psi[params.nx_local-1]);
        bool is_last = true;
        bool is_first = true;
        for(size_t i = 1; i < params.nx_local-1; i++)
        {
            // Compute reflection coefficient
            double xi = params.x_min + i*params.dx;
            if(xi < -a)
            {
                sum_R += cabs(psi[i])*cabs(psi[i]); 
            }else if(is_last==true){
                is_last = false;
                sum_R += 0.5*cabs(psi[i])*cabs(psi[i]);
            }
            // Compute reflection coefficient
            if(xi >= a)
            {
                sum_T += cabs(psi[i])*cabs(psi[i]);
            } else if(is_first==true && is_last == false){
                is_first = false;
                sum_T += 0.5*cabs(psi[i])*cabs(psi[i]);
            }
        }
        fprintf(user_fp,"%.6e %.6e %.6e\n", t, params.dx*sum_R, params.dx*sum_T);
    } 


    return;
}
