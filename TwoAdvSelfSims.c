/* TwoAdvSelfSims.c 

Simulation calculating fixation probability of second beneficial allele, given existing
ben allele at initial frequency p. For use in the study "Linkage and the limits to natural
selection in partially selfing populations".

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.

This program can be compiled with e.g. GCC using a command like:
gcc TwoAdvSelfSims -lm -lgsl -lgslcblas -I/usr/local/include -L/usr/local/lib TwoAdvSelfSims.c

Then run by executing:
./TwoAdvSelfSims N self rec ha sa hb sb p reps
Where:
- N is the population size
- self is the rate of self-fertilisation
- rec is recombination rate
- ha, hb is dominance at the original, introduced beneficial allele
- sa, sb is selection coefficient of the original, introduced beneficial allele
- p is the initial frequency of the first beneficial allele (when the second is introduced)
- reps is how many times the second allele should FIX before simulation stops 
(the number of actual runs is greater due to stochastic loss of second allele)

Note that haplotypes are defined as:
x1 = ab
x2 = Ab
x3 = aB
x4 = AB

Genotypes defined as:
g11 = g1 = ab/ab
g12 = g2 = Ab/ab
g13 = g3 = aB/ab
g14 = g4 = AB/ab
g22 = g5 = Ab/Ab
g23 = g6 = Ab/aB
g24 = g7 = Ab/AB
g33 = g8 = aB/aB
g34 = g9 = aB/AB
g44 = g10 = AB/AB

Output files are the parameters;
followed by number of times each haplotype fixed;
followed by average total generations elapsed in each case;
Then total number of simulations ran;
Then fixation prob of allele, both unscaled and scaled to unlinked case, 
along with 95% CI intervals for the latter case;
then number of allele fixations.

Note that 'fixation' DIFFERS depending on the inputs of sa, sb.
If sa >= sb (interference case) then 'fixation' counts as fixation of second allele on any genetic background.
If sa < sb (replacement case) then 'fixation' only considers fixation of second allele with neutral haplotype
(I.e. where the 'less fit' neutral allele at locus A fixes, instead of selected allele).

*/

/* Preprocessor statements */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/* Function prototypes */
void geninit(double *geninit, double FIS, const gsl_rng *r);
void selection(double *geninit);
void reproduction(double *geninit);
unsigned int hcheck(double *geninit, double *haps, unsigned int *hf, unsigned int stype);

/* Global variable declaration */
unsigned int N = 0;		/* Pop size */
double rec = 0;			/* Recombination rate */
double self = 0;		/* Rate of self-fertilisation */
double ha = 0;			/* Dominance of site A */
double sa = 0;			/* Fitness of site A */
double hb = 0;			/* Dominance of site B */
double sb = 0;			/* Fitness of site B */
double pee = 0;			/* Freq of initial sweep */

/* Main program */
int main(int argc, char *argv[]){
	unsigned int i; 				/* A counter */
	unsigned int reps;				/* Length of simulation (no. of introductions of neutral site) */
	unsigned int stype = 0;				/* What type of sim (replacement or hitch-hiking)? */
	unsigned int nsfix = 0;			/* sims where target type fixed */
	unsigned int nstot = 0;			/* total sims ran */	
	unsigned int gens = 0;			/* Number gens elapsed */
	unsigned int isfin = 0;			/* Is sim finished? */
	unsigned int hf = 0;			/* The hap that fixed */
	double pf = 0;					/* Overall fix prob */
	double FIS = 0;					/* Wright's FIS */
	double StdFix = 0;				/* Standard Fixation prob if unlinked */
	double citop = 0;
	double cibot = 0;
	double nsfix2 = 0;			
	double nstot2 = 0;			
	char selfchar[10];
	char recchar[15];
	char hchar[10];
	char pchar[10];
	char fname[64];
	FILE *ofp_tr;					/* Pointer for tree output */
	
	/* GSL random number definitions */
	const gsl_rng_type * T; 
	gsl_rng * r;
	
	/* This reads in data from command line. */
	if(argc != 10){
		fprintf(stderr,"Invalid number of input values.\n");
		exit(1);
	}
	N = strtod(argv[1],NULL);
	self = strtod(argv[2],NULL);
	rec = strtod(argv[3],NULL);
	ha = strtod(argv[4],NULL);
	sa = strtod(argv[5],NULL);
	hb = strtod(argv[6],NULL);
	sb = strtod(argv[7],NULL);
	pee = strtod(argv[8],NULL);
	reps = strtod(argv[9],NULL);
	
	if(sa >= sb){
		stype = 0;
	}else if(sa < sb){
		stype = 1;
	}
	FIS = self/(2.0-self);
	
	/* Arrays definition and memory assignment */
	double *genotype = calloc(10,sizeof(double));				/* Genotype frequencies */
	unsigned int *gensamp = calloc(10,sizeof(unsigned int));	/* New population samples */
	double *haps = calloc(4,sizeof(double));					/* Haplotypes */
	unsigned int *pfix = calloc(4,sizeof(unsigned int));		/* Haplotypes that fix */
	unsigned int *tfix = calloc(4,sizeof(unsigned int));		/* time that haps fix */
	  
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */
     
	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	nsfix = 0;
	nstot = 0;
	
    while(nsfix < reps){
    	nstot++;
	
		/* Initialising genotypes */
		geninit(genotype,FIS,r);
		gens = 0;
		isfin = 0;
		hf = 0;
		
		while(isfin == 0){
    		/* Selection routine */
    		selection(genotype);
    		
	    	/* Reproduction routine */
    		reproduction(genotype);
       		
    		/* Sampling based on new frequencies */
	       	gsl_ran_multinomial(r,10,N,genotype,gensamp);
       		for(i = 0; i < 10; i++){
    	   		*(genotype + i) = (*(gensamp + i))/(1.0*N);
	       	}
	       	
	       	gens++;
	       	isfin = hcheck(genotype, haps, &hf, stype);
       	}
       	
       	if(isfin == 1){
       		(*(pfix + hf))++;
    	  	(*(tfix + hf)) += gens;
	       	if(stype == 0){
       			nsfix += (*(haps + 2) + *(haps + 3));
       		}else if(stype == 1){
    	   		nsfix += (*(haps + 2));
	       	}
       	}
    
	}	/* End of simulation */
	
	pf = nsfix/(1.0*nstot);
	StdFix = 2*sb*((hb+FIS-hb*FIS)/(1+FIS));		/* Fix prob of new allele if unlinked */
	
	nstot2 = nstot + 3.84;
	nsfix2 = (1.0/(nstot2))*(nsfix + (3.84/2.0));
	citop = nsfix2 + 1.96*sqrt((1.0/(nstot2))*nsfix2*(1.0-nsfix2));
	cibot = nsfix2 - 1.96*sqrt((1.0/(nstot2))*nsfix2*(1.0-nsfix2));	
	citop = citop/(1.0*StdFix);
	cibot = cibot/(1.0*StdFix);	
	
	/* Printing solutions to file */
	/* First, converting values to strings */
    sprintf(selfchar, "%0.2lf",self);
    sprintf(recchar, "%0.5lf",rec);
    sprintf(hchar, "%0.3lf",ha);
    sprintf(pchar, "%0.5lf",pee);            

	strcpy(fname,"sim_self");
	strcat(fname,selfchar);
	strcat(fname,"rec");
	strcat(fname,recchar);
	strcat(fname,"h");
	strcat(fname,hchar);
	strcat(fname,"p");
	strcat(fname,pchar);
	strcat(fname,".sim");
	
	ofp_tr = fopen(fname,"a+");
	fprintf(ofp_tr,"%d %lf %lf %lf %lf %lf %lf %lf ",N,self,rec,ha,sa,hb,sb,pee);
	for(i = 0; i < 4; i++){
		fprintf(ofp_tr,"%d ",*(pfix + i));
	}
	for(i = 0; i < 4; i++){
		fprintf(ofp_tr,"%lf ",((*(tfix + i)))/(1.0*(*(pfix + i))));
	}
	fprintf(ofp_tr,"%d %lf %lf %lf %lf %d\n",nstot,pf,pf/(1.0*StdFix),citop,cibot,reps);
	fclose(ofp_tr);	
	
	/* Freeing memory and wrapping up */
 	gsl_rng_free(r);
 	free(tfix);
 	free(pfix);
 	free(haps);
 	free(gensamp);
	free(genotype);
	return 0;
}

/* Initialising genotypes */
void geninit(double *geninit, double FIS, const gsl_rng *r){

	unsigned int htype = 0;		/* Type of het assignment */
	
	double *ptype = calloc(3,sizeof(double));
	unsigned int *ctype = calloc(3,sizeof(double));	
		
	/* Routine to determine initial background of second mutant, given first is at frequency p. */
	*(ptype + 0) = (1-pee)*(1-pee) + FIS*pee*(1-pee);
	*(ptype + 1) = 2*pee*(1-pee)*(1-FIS);
	*(ptype + 2) = pee*pee + FIS*pee*(1-pee);
	
	/* First initialise baseline freqs */
	*(geninit + 0) = *(ptype + 0);
	*(geninit + 1) = *(ptype + 1);
	*(geninit + 2) = 0;
	*(geninit + 3) = 0;
	*(geninit + 4) = *(ptype + 2);
	*(geninit + 5) = 0;
	*(geninit + 6) = 0;
	*(geninit + 7) = 0;
	*(geninit + 8) = 0;
	*(geninit + 9) = 0;
	
	/* Then decide where to add new mutant */
	gsl_ran_multinomial(r,3,1,ptype,ctype);
	if(*(ctype + 0) == 1){
		*(geninit + 0) -= 1/(1.0*N);
		*(geninit + 2) += 1/(1.0*N);
	}else if(*(ctype + 1) == 1){
		*(geninit + 1) -= 1/(1.0*N);
	  	htype = gsl_ran_bernoulli(r,0.5);
	  	if(htype == 0){
	  		*(geninit + 5) += 1/(1.0*N);
	  	}else if(htype == 1){
	  		*(geninit + 3) += 1/(1.0*N);
	  	}
	}else if(*(ctype + 2) == 1){
		*(geninit + 4) -= 1/(1.0*N);
		*(geninit + 6) += 1/(1.0*N);
	}
	
	free(ctype);
	free(ptype);

}	/* End of gen initiation routine */

/* Selection routine */
void selection(double *geninit){
	/* Fitness of each genotype */
	double W11, W12, W13, W14, W22, W23, W24, W33, W34, W44;		
	double Wmean;				/* Mean fitness */
	
	W11 = 1;
	W12 = 1 + ha*sa;
	W13 = 1 + hb*sb;
	W14 = 1 + ha*sa + hb*sb;
	W22 = 1 + sa;
	W23 = 1 + ha*sa + hb*sb;
	W24 = 1 + sa + hb*sb;
	W33 = 1 + sb;
	W34 = 1 + ha*sa + sb;
	W44 = 1 + sa + sb;
	
	/* Mean fitness calculation */
	Wmean = ((*(geninit + 0))*W11) + ((*(geninit + 1))*W12) + ((*(geninit + 2))*W13) + ((*(geninit + 3))*W14) + ((*(geninit + 4))*W22) + ((*(geninit + 5))*W23) + ((*(geninit + 6))*W24) + ((*(geninit + 7))*W33) + ((*(geninit + 8))*W34) + ((*(geninit + 9))*W44);
	
	/* Changing frequencies by selection */
	*(geninit + 0) = ((*(geninit + 0))*W11)/Wmean;
	*(geninit + 1) = ((*(geninit + 1))*W12)/Wmean;
	*(geninit + 2) = ((*(geninit + 2))*W13)/Wmean;
	*(geninit + 3) = ((*(geninit + 3))*W14)/Wmean;
	*(geninit + 4) = ((*(geninit + 4))*W22)/Wmean;
	*(geninit + 5) = ((*(geninit + 5))*W23)/Wmean;
	*(geninit + 6) = ((*(geninit + 6))*W24)/Wmean;
	*(geninit + 7) = ((*(geninit + 7))*W33)/Wmean;
	*(geninit + 8) = ((*(geninit + 8))*W34)/Wmean;
	*(geninit + 9) = ((*(geninit + 9))*W44)/Wmean;
	
}	/* End of selection routine */

/* Reproduction routine */
void reproduction(double *geninit){
	/* Fed-in genotype frequencies (for ease of programming) */
	double g11s, g12s, g13s, g14s, g22s, g23s, g24s, g33s, g34s, g44s;
	/* Haplotypes */
	double x1, x2, x3, x4;
	
	/* Initial definition of genotypes */
	g11s = *(geninit + 0);
	g12s = *(geninit + 1);
	g13s = *(geninit + 2);
	g14s = *(geninit + 3);
	g22s = *(geninit + 4);
	g23s = *(geninit + 5);
	g24s = *(geninit + 6);
	g33s = *(geninit + 7);
	g34s = *(geninit + 8);
	g44s = *(geninit + 9);
	
	/* Baseline change in haplotype frequencies */
	x1 = g11s + (g12s + g13s + g14s)/2.0 - ((g14s - g23s)*rec)/2.0;
	x2 = g22s + (g12s + g23s + g24s)/2.0 + ((g14s - g23s)*rec)/2.0;
	x3 = g33s + (g13s + g23s + g34s)/2.0 + ((g14s - g23s)*rec)/2.0;
	x4 = g44s + (g14s + g24s + g34s)/2.0 - ((g14s - g23s)*rec)/2.0;
	
	/* Change in SEXUAL frequencies (both outcrossing and selfing) */
	*(geninit + 0) = (g11s + (g12s + g13s + g14s*pow((1 - rec),2) + g23s*pow(rec,2))/4.0)*self + (1 - self)*pow(x1,2);
	*(geninit + 4) = (g22s + (g12s + g24s + g23s*pow((1 - rec),2) + g14s*pow(rec,2))/4.0)*self + (1 - self)*pow(x2,2);
	*(geninit + 7) = (g33s + (g13s + g34s + g23s*pow((1 - rec),2) + g14s*pow(rec,2))/4.0)*self + (1 - self)*pow(x3,2);
	*(geninit + 9) = (g44s + (g24s + g34s + g14s*pow((1 - rec),2) + g23s*pow(rec,2))/4.0)*self + (1 - self)*pow(x4,2);
	*(geninit + 1) = ((g12s + (g14s + g23s)*(1 - rec)*rec)*self)/2.0 + 2.0*(1 - self)*x1*x2;
	*(geninit + 2) = ((g13s + (g14s + g23s)*(1 - rec)*rec)*self)/2.0 + 2.0*(1 - self)*x1*x3;
	*(geninit + 3) = ((g14s*pow((1 - rec),2) + g23s*pow(rec,2))*self)/2.0 + 2.0*(1 - self)*x1*x4;
	*(geninit + 5) = ((g23s*pow((1 - rec),2) + g14s*pow(rec,2))*self)/2.0 + 2.0*(1 - self)*x2*x3;
	*(geninit + 6) = ((g24s + (g14s + g23s)*(1 - rec)*rec)*self)/2.0 + 2.0*(1 - self)*x2*x4;
	*(geninit + 8) = ((g34s + (g14s + g23s)*(1 - rec)*rec)*self)/2.0 + 2.0*(1 - self)*x3*x4;
		
}	/* End of reproduction routine */

/* Has any allele fixed or not? */
unsigned int hcheck(double *geninit, double *haps, unsigned int *hf, unsigned int stype){
	/* Fed-in genotype frequencies (for ease of programming) */
	double g11s, g12s, g13s, g14s, g22s, g23s, g24s, g33s, g34s, g44s;
	unsigned int retval = 0;
	
	/* Initial definition of genotypes */
	g11s = *(geninit + 0);
	g12s = *(geninit + 1);
	g13s = *(geninit + 2);
	g14s = *(geninit + 3);
	g22s = *(geninit + 4);
	g23s = *(geninit + 5);
	g24s = *(geninit + 6);
	g33s = *(geninit + 7);
	g34s = *(geninit + 8);
	g44s = *(geninit + 9);
	
	/* Calculation of haplotypes */
	*(haps + 0) = g11s + (g12s + g13s + g14s)/2.0;
	*(haps + 1) = g22s + (g12s + g23s + g24s)/2.0;
	*(haps + 2) = g33s + (g13s + g23s + g34s)/2.0;
	*(haps + 3) = g44s + (g14s + g24s + g34s)/2.0;
	
/* 	printf("Haps are %lf %lf %lf %lf\n",*(haps + 0),*(haps + 1),*(haps + 2),*(haps + 3));*/
	
	if(*(haps + 0) == 1){
		retval = 1;
		*hf = 0;
	}
	else if(*(haps + 1) == 1){
		retval = 1;
		*hf = 1;
	}
	else if(*(haps + 2) == 1){
		retval = 1;
		*hf = 2;
	}
	else if(*(haps + 3) == 1){
		retval = 1;
		*hf = 3;
	}else if(stype == 1){
		if( (*(haps + 2) + *(haps + 3) ) == 0 ){
			retval = 2;
			*hf = 0;
		}
	}
	
	return retval;
		
}	/* End of hap check routine */

/* End of program */
