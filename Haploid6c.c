#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
//#include <sys/time.h>
//#include "nrutil.h"
/* compile using 
 * gcc -g hickeyHaploid6c.c  poidev.c ran1.c -lm -o hickeyHaploid6c
 *
 *  Modified July 24 2018 to also record the min as well as max sex fit
 *  Modified July 30 2018 to also record the stddev as well as max sex fit
 *
 */

#define WITHREPLACEMENT 0
// #define GENERATIONS 201
#define GENERATIONS 101
#define MAX_SIZE     400000
#define INITIAL_SIZE     100000
#define INITIAL_ASEXUALS 100000
#define NO_LOCI 100
#define MUT_RATE 1.0e-08     // rate per locus
#define RECOMB_RATE 0.01  // rate per individual
#define AVG_ALLELE_FREQ 0.01 // initial advantagous allele freq
#define FECUNDITY 2.0 // fecundity
#define SWAP(a,b) tmp1=(a);(a)=(b);(b)=tmp1;


extern double ran1(long *iseed);
extern double poidev(double xm, long *idum);
// extern double gammln(double xx);
// extern void quicksort(long *list, long i, long j);
long *iseed, idum;
time_t tp;

struct {
   long sex; 
   double fitness;
   double relativeFitness;
   long maternal_allele[NO_LOCI];
   long paternal_allele[NO_LOCI];
} n[MAX_SIZE], m[MAX_SIZE], o[MAX_SIZE], tmp, tmp0, tmp1, tmp2;
long  exceeded=1;
long  i, j, jj, k, l, iter, numberNew, gen, numberInd;
long  tmp_allele, save_k;
double arandom;
long above;
double fit_asexual, max_asexual, min_asexual, rel_asexual;
double fit_sexual, max_sexual, rel_sexual;
double max, max_fitness;
long no_asexual;
long no_sexual, no_recomb=0;
long count;
double meanNoGoodAllelesPerInd=0.0;
double sumSquares=0.0;

// Function Declarations
void usage(void);
void recombine(long i);
void calculateFitness(void);
long  find_a_mate(long i);
void join(long i, long j, long k, long flag);
void permute(void); // not currently used
double myselect(long i, long j); 
void countInGametes(void);  // for debugging
void countInZygotes(void);  // for debugging

int main(int argc, char *argv[]) {
    iseed = &idum;
    idum = (unsigned) time(&tp);
    idum = -idum;
    // idum = -12345;
// *****   initialize   *******************************************
    tmp0.sex = 0; tmp0.fitness=0.0; // set up a null value tmp0
    for(i=0; i<NO_LOCI; i++) { 
	tmp0.paternal_allele[i]=0; 
	tmp0.maternal_allele[i]=0; 
    }
    for(i=0; i<INITIAL_ASEXUALS; i++) { 
	n[i].sex=1;                         // asexual is 1
	n[i].fitness=1; 
	for(j=0; j<NO_LOCI; j++) {
	    count=0;
	    if(ran1(iseed) < AVG_ALLELE_FREQ) {
		n[i].paternal_allele[j]=1;  // advantagous is 1 else 0
		count++;
	    } else { n[i].paternal_allele[j]=0; }
	    if(ran1(iseed) < AVG_ALLELE_FREQ) {
		n[i].maternal_allele[j]=1;  // advantagous is 1 else 0
		count++;
	    } else { n[i].maternal_allele[j]=0; }
	    n[i].fitness=n[i].fitness*(1+count*0.05);  // multiplicative fitness
	}
    }
    for(i=INITIAL_ASEXUALS; i<INITIAL_SIZE; i++) { 
	n[i].sex=0;                         // sexual is 0 
	n[i].fitness=1; 
	for(j=0; j<NO_LOCI; j++) {
	    count=0;
	    if(ran1(iseed) < AVG_ALLELE_FREQ) {
		n[i].paternal_allele[j]=1;  // advantagous is 1 else 0
		count++;
	    } else { n[i].paternal_allele[j]=0; }
	    if(ran1(iseed) < AVG_ALLELE_FREQ) {
		n[i].maternal_allele[j]=1;  // advantagous is 1 else 0
		count++;
	    } else { n[i].maternal_allele[j]=0; }
	    n[i].fitness=n[i].fitness*(1+count*0.05);  // multiplicative fitness; dominance
	}
    }
    numberInd=INITIAL_SIZE;
    numberInd=100000;
    gen = 0; 
    // CALCULATE AVG ASEXUAL AND SEXUAL FITNESS
    calculateFitness();
    fit_asexual=(double) fit_asexual/(double) no_asexual;
    fit_sexual=(double) fit_sexual/(double) no_sexual;
    // standard deviation s^2 = \frac{1}{N-1}(ss-N*\mu^2)
    sumSquares = sqrt((sumSquares - no_asexual*fit_asexual*fit_asexual)/(no_asexual-1.0));

    //---------------------------------------------
    fprintf(stdout,"\nDiploid Model, single locus with selection based on twice median \n");
    arandom=MUT_RATE;
    arandom=0.0;
    fprintf(stdout,"relative fitness to determine who survives).  Mutation rate is %10.5g,\n",arandom);
    arandom=RECOMB_RATE;
    fprintf(stdout,"Recomb rate is NULL per sexual individual per generation\n");
    fprintf(stdout,"Starting values at generation %3d are; population size %6d,\n",gen, numberInd);
//  fprintf(stdout,"  with %5d asexuals with avg fitness of %10.5f versus ,\n",no_asexual,fit_asexual);
    fprintf(stdout,"%5d sexuals  with avg fitness of %10.5f.\n",no_sexual,fit_sexual);
    fprintf(stdout,"\n");
    fprintf(stdout,"                       avg         avg         min        max       max    ");
    fprintf(stdout," Alleles      rel         rel \n");
    fprintf(stdout,"gen popSize asexuals asexualFit  sexualFit  asexualFit  asexualFit sexualFit ");
    fprintf(stdout,"  perInd   asexualFit  sexualFit stddev\n");
    //---------------------------------------------
    
    for(gen=1; gen<GENERATIONS; gen++) {
	// MUTATE
	for(i=0; i<numberInd; i++) {
	    for(j=0; j<NO_LOCI; j++) {
		if(ran1(iseed) < MUT_RATE) n[i].paternal_allele[j]=0;
		if(ran1(iseed) < MUT_RATE) n[i].maternal_allele[j]=0;
	    }
	}
	//
	// CREATE RECOMBINED GAMETES PROPORTIONAL TO FECUNDITY*RelFitness
	for(i=0, j=0; i<numberInd; i++) {
	    if(n[i].sex==0) { // sexual; choose number of gametes k
		arandom = poidev(FECUNDITY*n[i].relativeFitness,iseed); 
		k = (long) arandom; 
		for(l=0; l<k; l++) {
		    recombine(i); // recombines n[i] and puts into tmp
		    m[j++]=tmp;
		}
	    } else { // asexual has twice fecundity but no recombination
		arandom = poidev(FECUNDITY*n[i].relativeFitness,iseed); 
		k = (long) arandom; 
		for(l=0; l<k; l++) m[j++]=n[i];
	    }
	}
	// MATE VIA RANDOM UNION (WITH REPLACEMENT) OF SEXUAL GAMETES
	numberNew=j;  // total number of gametes
	if(WITHREPLACEMENT) {
	    max=INITIAL_SIZE; if(numberNew < max) max=numberNew;
	    for(i=0; i<max; i++) { 
		j = (long) numberNew*ran1(iseed); // choose gamete at random
		if(m[j].sex==0) { // individual is sexual 
		    k = find_a_mate(j);
		    join(i,j,k,0);  // joins j and k into i
		} else { // is asexual
		    n[i]=m[j];
		}
	    }
	} else { // WITHOUT REPLACEMENT  // taking m[] and putting into next generation n[]
	    // premute first and then pair off the sexuals.
            // fprintf(stdout,"before permutation\n");
            // countInGametes();
	    for(above=0, i=0; i<numberNew; i++) {
		j = (long) numberNew * ran1(iseed);
		tmp = m[i];
		m[i] = m[j];
		m[j] = tmp;
	    }
            // fprintf(stdout,"after permutation\n");
            // countInGametes();
	    for(k=0, i=0; i<numberNew; i++) {
		if(m[i].sex==0) {
		    j=i;
		    while(j++ < numberNew) {
			if(m[j].sex==0) {
			    join(k,i,j,2); // join only copies alleles+sex
                            // fprintf(stdout,"Joined i=%5d, j=%5d to k=%5d\n",i,j,k); 
			    k++;
			    break;
			}
		    }
		    // if(j==numberNew) fprintf(stdout," No more sexuals \n");
		} else {
		    if(m[i].sex==1) n[k++] = m[i];
                    // fprintf(stdout,"Found asexual (sex=1) or an already mated gamete (sex=2): %5d\n",i);
		}
	    }
	    numberNew=k-1;
	    if(numberNew<INITIAL_SIZE) numberInd=numberNew;
	    else numberInd=INITIAL_SIZE;
            // countInZygotes();
	}
	// Now have new collection of m genotypes, reset to n
	// CALCULATE FITNESS
	// calculate avg asexual and sexual fitness
	calculateFitness();
	if(no_asexual>0) { fit_asexual=fit_asexual/no_asexual; } else fit_asexual=0;
	if(no_sexual>0) { fit_sexual=fit_sexual/no_sexual; } else fit_sexual=0;
        // standard deviation s^2 = \frac{1}{N-1}(ss-N*\mu^2)
        sumSquares = sqrt((sumSquares - no_asexual*fit_asexual*fit_asexual)/(no_asexual-1.0));
	//
	//  fprintf(stdout,"                       avg        avg        max        max   \n");
        //  fprintf(stdout,"gen popSize asexuals asexualFit sexualFit asexualFit sexualFit\n");
	fprintf(stdout,"%3d  %6d  %6d %10.5f %10.5f %10.5f %10.5f %10.5f",gen,numberInd,no_asexual,fit_asexual,fit_sexual,min_asexual, max_asexual,max_sexual);
        // if(max_sexual > max_asexual) max_fitness = max_sexual else max_fitness = max_asexual;
        max_fitness = max_sexual > max_asexual ? max_sexual : max_asexual;
        rel_asexual = fit_asexual/max_fitness;
        rel_sexual  = fit_sexual/max_fitness;
	fprintf(stdout,"%10.5f ",meanNoGoodAllelesPerInd);
	fprintf(stdout,"%10.5f ",rel_asexual);
	fprintf(stdout,"%10.5f ",rel_sexual);
	fprintf(stdout,"%10.5f \n",sumSquares);
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void usage(void) {
    fprintf(stdout,"\n");
    fprintf(stdout," This program does drift-over-time with \n");
    fprintf(stdout," with selection on offspring variation for Donal.\n");
    fprintf(stdout,"\n");
}

//----------------------------------------------------------------------------------------
void recombine(long i) { // RECOMBINE: 
//  takes an individual i and recombine it's chromosomes
//  return the recombined individual in struct tmp
    long j, k, tmp_allele;
    if(n[i].sex==0) { // double check it is sexual
	tmp=n[i];     // copy n[i] into tmp
	k=(long) (NO_LOCI-1)*ran1(iseed); // Always recombine; choose a spot for recomb
	for(j=k+1; j<NO_LOCI; j++) {
	    tmp_allele=tmp.paternal_allele[j];
	    tmp.paternal_allele[j]=tmp.maternal_allele[j];
	    tmp.maternal_allele[j]=tmp_allele;
	}
    } else { // not sexual; recombination doesn't make sense
	fprintf(stdout," In recombine: trying to recombine an asexual\n");
    }
    return;
}

//----------------------------------------------------------------------------------------
void calculateFitness(void) { // Calculate avg asexual and sexual fitness
    long i, j, midpoint;
    double median, countSum;
    fit_asexual=0.0;  max_asexual=0.0;  no_asexual=0;
    fit_sexual=0.0;   max_sexual=0.0;   no_sexual=0;
    sumSquares=0.0;
    min_asexual=100.0;
    max=0.0;  // Find max fitness of all n[]
    countSum=0.0;
    for(i=0; i<numberInd; i++) { 
	n[i].fitness=1.0;
	for(j=0; j<NO_LOCI; j++) { 
	    count=0;
	    if(n[i].paternal_allele[j]==1) count++;
	    if(n[i].maternal_allele[j]==1) count++;
	    n[i].fitness=n[i].fitness*(1+count*0.05);
            if(n[i].paternal_allele[j]==1) countSum++;
            if(n[i].maternal_allele[j]==1) countSum++;
	}
	if(n[i].fitness > max) max = n[i].fitness;
	if(n[i].sex) {   // only if asexual=1
	    if(n[i].fitness > max_asexual) max_asexual = n[i].fitness;
	    if(n[i].fitness < min_asexual) min_asexual = n[i].fitness;
            sumSquares += n[i].fitness * n[i].fitness;
	    fit_asexual+=n[i].fitness; 
	    no_asexual++;
	} else {
	    if(n[i].fitness > max_sexual) max_sexual = n[i].fitness;
	    fit_sexual+=n[i].fitness; 
	    no_sexual++;
	}
    }
    meanNoGoodAllelesPerInd = countSum/(numberInd);
    // Now make fitnesses relative
    // First need to find the median fitness -- use quickselect to do this
    midpoint = numberInd/2;
    // do myselect on a copy 'o' of 'n'
    for(i=0; i < numberInd; i++) o[i]=n[i];
    median = myselect(midpoint, numberInd-1);
    // Make fitnesses relative to the 2*(x/median)
    for(i=0; i < numberInd; i++) {  
	// n[i].relativeFitness = n[i].fitness/max;
	n[i].relativeFitness = 2.0*(n[i].fitness/median);
    }
    return;
}

long find_a_mate(long i) { // find another sexual other than i
    long j, k;
    j = (long) numberNew*ran1(iseed);
    k=j;
    while(j <= numberNew) {
	if(m[j].sex==0 && j!=i) { // individual is sexual
	    return(j);
	} else { 
	    j++;
	    if(j > numberNew) j=0;
	    if(j == k) { fprintf(stdout," No other sexual but %d\n",i); return(i); }
	}
    }
    fprintf(stdout," In find_a_mate -- I shouldn't be here \n"); return(i); 
}

void join(long i,long j,long k, long flag) {  // mates individuals j and k to produce i
    // i will have one chromosome from j and one from k
    // I use flag to enforce WITHOUT REPLACEMENT
    m[k].sex=flag;  // if flag=2 ... can't mate in another round
    n[i].sex=0; // must be sexual at this point
    // leave fitness and relativeFitness undefined
    arandom=ran1(iseed);
    if(arandom < 0.5) {
	for(l=0; l < NO_LOCI; l++) n[i].paternal_allele[l]=m[j].paternal_allele[l];
    } else {
	for(l=0; l < NO_LOCI; l++) n[i].paternal_allele[l]=m[j].maternal_allele[l];
    }
    arandom=ran1(iseed);
    if(arandom < 0.5) {
	for(l=0; l < NO_LOCI; l++) n[i].maternal_allele[l]=m[k].paternal_allele[l];
    } else {
	for(l=0; l < NO_LOCI; l++) n[i].maternal_allele[l]=m[k].maternal_allele[l];
    }
}

void permute(void) {
    //
    // SELECT A NEW POP BY SAMPLE WITHOUT REPLACMENT
    // Do this via a permutation of all samples
    // and pick the first n.
    //
    // Permute
    for(i=0; i<numberInd; i++) { 
        j = (long) numberInd * ran1(iseed);
        tmp = n[i];
        n[i] = n[j];
        n[j] = tmp;
    }
    if(numberInd > INITIAL_SIZE) { 
        for(i=INITIAL_SIZE+1; i<numberInd; i++) { n[i]=tmp0; } // any beyond are null
        numberInd=INITIAL_SIZE;
    }
}

   
// FROM NRC ...
// using partioning to find a quick moment (here the median) of an array
// Hate passing arrays and hence will avoid the big pass and make the
// array a global variable.
// double select(unsigned long k, unsigned long n, double arr[]) {
//    unsigned long i,ir,j,l,mid;
// 
double myselect(long k, long m) {
    long i,ir,j,l,mid;
    // structures defined in preamble; o, tmp1, tmp2;

    l=1;
    ir=m;
    for (;;) {
        if (ir <= l+1) {
            if (ir == l+1 && o[ir].fitness < o[l].fitness) { SWAP(o[l],o[ir]) }
            return o[k].fitness;
        } else {
            mid=(l+ir) >> 1;
            SWAP(o[mid],o[l+1])
            if (o[l].fitness > o[ir].fitness) { SWAP(o[l],o[ir]) }
            if (o[l+1].fitness > o[ir].fitness) { SWAP(o[l+1],o[ir]) }
            if (o[l].fitness > o[l+1].fitness) { SWAP(o[l],o[l+1]) }
            i=l+1;
            j=ir;
            tmp2=o[l+1];
            for (;;) {
                do i++; while (o[i].fitness < tmp2.fitness);
                do j--; while (o[j].fitness > tmp2.fitness);
                if (j < i) break;
                SWAP(o[i],o[j])
            }
            o[l+1]=o[j];
            o[j]=tmp2;
            if (j >= k) ir=j-1;
            if (j <= k) l=i;
        }
    }
}

void countInZygotes(void) {
    long i, sum=0; 
    double a;
    for(i=0; i< numberNew; i++) {
        sum += n[i].maternal_allele[0];
        sum += n[i].paternal_allele[0];
    }
    fprintf(stdout," Sum is %10d out of numberNew                %10d\n",sum,numberNew);
    sum=0;
    for(i=0; i< numberInd; i++) {
        sum += n[i].maternal_allele[0];
        sum += n[i].paternal_allele[0];
    }
    a=(double) sum/numberInd;
    fprintf(stdout," Sum is %10d out of numberInd                %10d; %10.4f %\n",sum,numberInd,a);
    sum=0;
    for(i=numberInd; i< numberNew; i++) {
        sum += n[i].maternal_allele[0];
        sum += n[i].paternal_allele[0];
    }
    a=(double) sum/ (double) (numberNew-numberInd);
    fprintf(stdout," Sum is %10d out of (numberNew - numberInd)  %10d; %10.4f %\n",sum,numberNew-numberInd,a);
}

void countInGametes(void) {
    long i, sum=0; 
    double a;
    for(i=0; i< numberNew; i++) {
        sum += m[i].maternal_allele[0];
        sum += m[i].paternal_allele[0];
    }
    fprintf(stdout," Sum is %10d out of numberNew                      %10d\n",sum,numberNew);
    sum=0;
    for(i=0; i< numberInd; i++) {
        sum += m[i].maternal_allele[0];
        sum += m[i].paternal_allele[0];
    }
    a=(double) sum/numberInd;
    fprintf(stdout," Sum in m is %10d out of numberInd                 %10d; %10.4f %\n",sum,numberInd,a);
    sum=0;
    for(i=numberInd; i< numberNew; i++) {
        sum += m[i].maternal_allele[0];
        sum += m[i].paternal_allele[0];
    }
    a=(double) sum/(numberNew-numberInd);
    fprintf(stdout," Sum in m is %10d out of (numberNew - numberInd)   %10d; %10.4f %\n",sum,numberNew-numberInd,a);
}
