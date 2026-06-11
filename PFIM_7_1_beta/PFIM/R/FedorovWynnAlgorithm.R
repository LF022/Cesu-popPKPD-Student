#' @description The class \code{FedorovWynnAlgorithm} implements the FedorovWynn algorithm.
#' @title FedorovWynnAlgorithm
#' @inheritParams Optimization
#' @param elementaryProtocols List of elementary protocols
#' @param numberOfSubjects Numeric vector specifying number of subjects
#' @param showProcess Logical indicating whether to show process
#' @param proportionsOfSubjects Numeric vector of subject proportions
#' @param FedorovWynnAlgorithmOutputs A list giving the output of the optimization algorithm.
#' @include Optimization.R
#' @export

FedorovWynnAlgorithm = new_class("FedorovWynnAlgorithm", package = "PFIM", parent = Optimization,

                                 properties = list(elementaryProtocols = new_property(class_list, default = list()),
                                                   numberOfSubjects = new_property(class_vector, default = 0.0),
                                                   proportionsOfSubjects = new_property(class_vector, default = 0.0),
                                                   showProcess = new_property(class_logical, default = FALSE),
                                                   FedorovWynnAlgorithmOutputs = new_property(class_list, default = list())
                                 ))

# ==============================================================================
#' Helper: Build an Evaluation from an Optimization and Design
#' @name .buildEvalFromOptimization
#' @param optimizationObject An Optimization object
#' @param design A Design object
#' @return An Evaluation object ready to be run
#' @keywords internal
# ==============================================================================

.buildEvalFromOptimization = function( optimizationObject, design ) {

  # Extract properties from optimization
  name                     = prop( optimizationObject, "name" )
  modelEquations           = prop( optimizationObject, "modelEquations" )
  modelParameters          = prop( optimizationObject, "modelParameters" )
  modelError               = prop( optimizationObject, "modelError" )
  outputs                  = prop( optimizationObject, "outputs" )
  fimType                  = prop( optimizationObject, "fimType" )
  odeSolverParameters      = prop( optimizationObject, "odeSolverParameters" )
  modelFromLibrary         = prop( optimizationObject, "modelFromLibrary" )
  modelCovariates          = prop( optimizationObject, "modelCovariates" )
  modelCovariatesEquation  = prop( optimizationObject, "modelCovariatesEquation" )

  # Create and return Evaluation
  Evaluation(
    name                    = paste0( name, "_", prop( design, "name" ) ),
    modelEquations          = modelEquations,
    modelParameters         = modelParameters,
    modelError              = modelError,
    outputs                 = outputs,
    designs                 = list( design ),
    fimType                 = fimType,
    odeSolverParameters     = odeSolverParameters,
    modelFromLibrary        = modelFromLibrary,
    modelCovariates         = modelCovariates,
    modelCovariatesEquation = modelCovariatesEquation
  )
}

plotFrequenciesFedorovWynnAlgorithm = new_generic( "plotFrequenciesFedorovWynnAlgorithm", c( "optimization", "optimizationAlgorithm" ) )

# ==============================================================================
#' Fedorov-Wynn algorithm in Rcpp.
#' @name FedorovWynnAlgorithm_Rcpp
#' @description Run the FedorovWynnAlgorithm in Rcpp
#' @param protocols_input parameter protocols_input
#' @param ndimen_input parameter ndimen_input
#' @param nbprot_input parameter nbprot_input
#' @param numprot_input parameter numprot_input
#' @param freq_input parameter freq_input
#' @param nbdata_input parameter nbdata_input
#' @param vectps_input parameter vectps_input
#' @param fisher_input parameter fisher_input
#' @param nok_input parameter nok_input
#' @param protdep_input parameter protdep_input
#' @param freqdep_input parameter freqdep_input
#' @return A list giving the results of the outputs of the FedorovWynn algorithm.
#' @export
# ==============================================================================

FedorovWynnAlgorithm_Rcpp = function( protocols_input,  ndimen_input, nbprot_input,
                                      numprot_input, freq_input, nbdata_input,
                                      vectps_input, fisher_input, nok_input,
                                      protdep_input, freqdep_input ){
  incltxtFedorovWynnAlgorithm = '

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h> 	/* Mathematical functions */
#include <time.h>	/* Function time used to initialise the random number generator */
#include <float.h>	/* Implementation related constants */
#include <signal.h>	/* Signal handling used to detect arithmetic errors */
#include <Rcpp.h> /* Rcpp */

using namespace Rcpp;

#define FTOL 1e-6 /* Minimal change in criterion during optimisation step*/
#define FREQMIN 1e-8 /* Minimal frequency for a protocol to be kept */
#define EPSILON DBL_EPSILON
#define CRI_MAX 1e30 /* maximum value for criterium, used to test correct comp*/
#define NRANSI
#define tol 0.0001
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20
#define SMALL 1.0e-7
#define IA 16807  	/*definition des constantes pour ran1 */
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define FREE_ARG char*

typedef struct PROTOC /* Structure defi
ning an individual or elementary protocol */
{
	int	ntps; /* Number of sampling times */
	double	*tps; /* Vector of the sampling times */
	int	ndim; /* Number of random effects (=> size of matrices)*/
	double	*fisher; /* Fisher information matrix of the individual protocol */
} PROTOC;

typedef struct POPROT /* Structure defining a population protocol */
{
	int 	maxnp; /* Max number of individual protocols, needed for memory allocation*/
	int	np; /* Number of individual protocols WARNING np<=maxnp !*/
	int	ndim; /* Number of random effects (=> size of matrices)*/
	int	*num; /* In this program we keep trace of the index of individual protocols*/
	PROTOC	*pind; /* In NPML we would only need the individual protocols themselves*/
	double	*freq; /* Frequency of each individual protocols */
	double 	*fisher; /* Fisher information matrix of population protocol */
	double 	*finv; /* Inverse of the Fisher information matrix */
	double	det; /* Determinant of the matrix */
} POPROT;

typedef struct matrix
{
	int 	nrow;
	int	ncol;
	double	**m;
} matrix;


/******************************************************************************
   Construction and destruction of structures
*******************************************************************************/

void PROTOC_alloc( PROTOC *p, int ntps, int ndim)
{
	p->ntps = ntps;
	p->tps = (double *)calloc( p->ntps, sizeof(double));
	assert(p->tps!=NULL);
	p->ndim=ndim;
	p->fisher = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double));
	assert(p->fisher!=NULL);
	return;
}

POPROT *POPROT_alloc(int np,int ndim,int maxnp)
{	POPROT *p;
	p = (POPROT *)malloc(sizeof(POPROT));
	assert(p!=NULL);

	p->maxnp=maxnp;
	p->np = np;
	p->ndim=ndim;
/* Keep trace of the index of the individual protocols */
	p->num = (int *)calloc( p->maxnp, sizeof(int));
	assert(p->num!=NULL);
	p->pind = (PROTOC *)calloc( p->maxnp, sizeof(PROTOC));
	assert(p->pind!=NULL);
	p->freq = (double *)calloc( p->maxnp, sizeof(double));
	assert(p->freq!=NULL);
	p->fisher = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double));
	assert(p->fisher!=NULL);
	p->finv = (double *)calloc( p->ndim*(p->ndim+1)/2, sizeof(double));
	assert(p->finv!=NULL);
	p->det=1.0;
	return p;
}

void PROTOC_copy(PROTOC *p,PROTOC *p1)
{
	int i,ndat;

	for(i=0;i<p->ntps;i++) {p1->tps[i]=p->tps[i];}
	ndat=(int)(p->ndim*(p->ndim+1))/2;
	for(i=0;i<ndat;i++) {p1->fisher[i]=p->fisher[i];

		}
	return;
}

matrix *matrix_create(int nrow, int ncol)
{	int i,j;
	matrix *mat;
	mat = (matrix *)malloc(sizeof(matrix));
	assert(mat!=NULL);
	mat->nrow=nrow;
	mat->ncol=ncol;
	mat->m = (double **)calloc( nrow, sizeof(double *));
	assert(mat->m!=NULL);
	for( i=0 ; i<nrow ; i++ )
	{
	mat->m[i] = (double *)calloc( ncol, sizeof(double));
	assert(mat->m[i]!=NULL);
	for( j=0 ; j<ncol ; j++ )
		mat->m[i][j]=0.0;
	}
	return mat;

}

void PROTOC_print(PROTOC *p)
{	int i,icas=0,j;
	fprintf(stderr,"Temps du protocole :");
	for (i=0;i<p->ntps;i++) {}
	for (i=0;i<p->ndim;i++) {
		for (j=0;j<=i;j++) {
			icas++;}
	}
	return;
}

void POPROT_print(POPROT *pop,int nofish)
{	int i,icas=0,j,k;

	for(i=0;i<pop->np;i++)
	{
		for (j=0;j<pop->pind[i].ntps;j++)
		{
		}

		if(nofish==1) {

		icas=0;
		for (j=0;j<pop->pind[i].ndim;j++)
		{
			for (k=0;k<=j;k++)
			{

				icas++;
			}

		}

		}
	}
	return;
}

void POPROT_printfisher(POPROT *pop)
{	int icas,j,k;


	icas=0;
	for (j=0;j<pop->ndim;j++)
	{
		for (k=0;k<=j;k++)
		{

		icas++;
		}

	}

	icas=0;
	for (j=0;j<pop->ndim;j++)
	{
		for (k=0;k<=j;k++)
		{

		icas++;
		}

	}
	return;
}

void POPROT_calculefisher(POPROT *pop)
{	int i,iprot;

	for(i=0;i<pop->ndim*(pop->ndim+1)/2;i++)
	{
		pop->fisher[i]=0.0;
		for(iprot=0;iprot<pop->np;iprot++)
		{
	pop->fisher[i]=pop->fisher[i]+pop->pind[iprot].fisher[i]*pop->freq[iprot];
		}
	}
}

void matrix_print(matrix *mat)
{	int i,j;

	for (i=0;i<mat->nrow;i++) {
		for (j=0;j<mat->ncol;j++) {
}
}
	return;
}

void PROTOC_destroy(PROTOC *p)
{
	free(p->fisher);
	free(p->tps);
}

void POPROT_destroy(POPROT *p)
{
	free(p->finv);
	free(p->fisher);
	free(p->freq);
	free(p->pind);
	free(p->num);
}

void matrix_destroy(matrix *mat)
{
	int i;
	for( i=0 ; i<mat->nrow ; i++ )
	{
	free(mat->m[i]);
	}
	free(mat->m);
}

/******************************************************************************
   Fonctions calculatoires
*******************************************************************************/

double matrace(int iprot,PROTOC *prot,POPROT *pop);
double lik(POPROT *pop);
void lubksb(double **a, int n, int *indx, double b[]);
int ludcmp(double **a, int n, int *indx, double *d);

/* Random number generator */
double randu(double x);
double gauss1(long *idum);
float ran1(long *idum);

/******************************************************************************
   Defining global variables (see variables_npml.txt for details)
*******************************************************************************/

volatile sig_atomic_t ierr=0;	/* Parameter signaling a SIGFPE error */
				/* to be used with a signal handler */
long seed=-10;     	/* graine du generateur aleatoire */

/******************************************************************************

/******************************************************************************
   Initialisation et ajout de protocole ?l?mentaire
*******************************************************************************/

POPROT *initprot(PROTOC *allprot, IntegerVector protdep, NumericVector freqdep)

{
	int i,ij,np,nmax;

	POPROT *mypop;

	np=protdep[0];

	nmax=allprot[0].ndim*(allprot[0].ndim+1)/2+1;
	mypop=POPROT_alloc(np,allprot[0].ndim,nmax);

	for(i=0;i<np;i++) {
		mypop->freq[i]=freqdep[i];
		mypop->num[i]=protdep[i+1]-1;
		ij=mypop->num[i];
		PROTOC_alloc(&mypop->pind[i],allprot[ij].ntps,allprot[ij].ndim);
		PROTOC_copy(&allprot[ij],&mypop->pind[i]);
	}

return mypop;

}

/******************************************************************************
   Optimisation de protocole
*******************************************************************************/

void tassement(POPROT *pop)
{	int i,idec=0,np;
	np=pop->np;
	for(i=0;i<pop->np;i++)
	{
		if(pop->freq[i]<FREQMIN)
		{
			idec++;
			np--;
		}
		else
		{
			if (idec>0)
			{
				pop->freq[i-idec]=pop->freq[i];
				pop->num[i-idec]=pop->num[i];

				PROTOC_destroy(&pop->pind[i-idec]);
	PROTOC_alloc(&pop->pind[i-idec],pop->pind[i].ntps,pop->pind[i].ndim);
				PROTOC_copy(&pop->pind[i],&pop->pind[i-idec]);

			}
		}
	}

	for(i=np;i<pop->np;i++) PROTOC_destroy(&pop->pind[i]);
	pop->np=np;
	return;
}




int ajout(PROTOC *allprot,POPROT *pop,int nprot)
{	double cri,xmul,xtest,ifin=0,xdim,rcalc;
	int iprot,j,deja,qajout=-1;

	cri=pop->det;

	xtest=(double)(pop->ndim);
	xdim=xtest;

	for(iprot=0;iprot<nprot;iprot++)
	{


		deja=0;
		for(j=0;j<pop->np;j++)
		{
			if (iprot==pop->num[j]) deja=1;
		}
		if(fabs(deja)<EPSILON)
		{
			xmul=matrace(iprot,allprot,pop);
			if(xmul>=xtest)
			{
				xtest=xmul;
				qajout=iprot;
			}
		}
	}
	if(fabs(xtest-xdim)<EPSILON)
	{
		ifin=1;
	}
	else
	{
		rcalc=(xtest-xdim)/(xdim*(xtest-1));

		PROTOC_alloc(&pop->pind[pop->np],allprot[qajout].ntps,allprot[qajout].ndim);
		PROTOC_copy(&allprot[qajout],&pop->pind[pop->np]);
		pop->num[pop->np]=qajout;
		pop->freq[pop->np]=rcalc;


		for(j=0;j<pop->np;j++) pop->freq[j]=pop->freq[j]*(1-rcalc);
		pop->np++;
	}

	cri=lik(pop);
	return ifin;
}

int takeout_k(POPROT *pop,PROTOC *allprot,int *nnul)
{
	int ir=-1,j;
	double xmul;
	xmul=lik(pop); /*update Mf-1 */
	for(j=0;j<pop->np;j++)
	{
		xmul=0;
		if(pop->freq[j]<=FREQMIN)
		{
			xmul=matrace(pop->num[j],allprot,pop);
			if(xmul>(double)pop->ndim)
			{
				ir=j;
				*nnul=*nnul-1;
				return ir;
			}
		}
	}
	return ir;
}

int project_grad(POPROT *pop, PROTOC *allprot, double *gal,int ir,double
		*vmgal,int nnul)
{
/* Computes the projection of the gradient of ln(V)=ln(det(H(pop))) on the
surface sum(frequencies)=1 and freq[j]=0 for all j such as pop->freq[j]<FREQMIN
gal[j] is the jth component of the projected gradient (the direction)
gal[j]=0 if pop->freq[j]<FREQMIN and j<>ir
ir is -1 or the index removed from the active set at the previous step
vmgal is the largest abs(gal[j])
returns in, the last index j such as gal[j]<0
*/
	int j,in=-1;
	double sg,sgal;

	sg=lik(pop); /* Check that Mf-1 has been computed */
	sg=0;
	for(j=0;j<pop->np;j++)
	{
		gal[j]=0;
		if((pop->freq[j]>=FREQMIN) | (j==ir))
			gal[j]=matrace(pop->num[j],allprot,pop);
		sg=sg+gal[j];
	}
	sgal=0;*vmgal=0;
	for(j=0;j<pop->np;j++)
	{
		if((pop->freq[j]>=FREQMIN) | (j==ir))
		{
			gal[j]=gal[j]-sg/(double)(pop->np-nnul);
			if(fabs(gal[j])>*vmgal) *vmgal=(double)fabs(gal[j]);
			if(gal[j]<0) in=j;
			sgal=sgal+gal[j];
		}
	}
	return in;
}


int optim_lambda(POPROT *pop, double *gal, double *alkeep, double *crit,
	double dmax, double crit2, double vmgal)
{
/* dmax is too large a move along the direction
=> move along the direction by dmax/2**k until
-the move is too small : no new frequencies have been found (return 1)
-a better likelihood is found => update pop with the new frequencies
	-if the criteria has changed by more than FTOL, return 0
	-else return 1; the frequencies correspond to the last (smallest) move
*/
	int j;
	double ro=1,crit1;
	do
	{
		ro=ro*0.5;
		for(j=0;j<pop->np;j++) pop->freq[j]=alkeep[j]+ro*dmax*gal[j];
		crit1=lik(pop); /* Updates Mf-1 for new frequencies */
		if(crit1>*crit)
		{
			for(j=0;j<pop->np;j++) alkeep[j]=pop->freq[j];
			*crit=crit1;
			if(fabs((crit2-*crit)/crit2)>FTOL) {return 0;}
			else {return 1;}
		}
	} while((dmax*ro*vmgal)>1e-4);
	for(j=0;j<pop->np;j++) pop->freq[j]=alkeep[j];
	return 1;
}



int doptimal(POPROT *pop,PROTOC *allprot)
{
	int nnul=0, ir=-1,in=-1,nit=0;
	int j,imax,iopt;
	double crit,crit1,crit2,vmgal,dmax;
	double *gal,*alkeep;
	double sal,sal1;

	gal=(double *)calloc(pop->np,sizeof(double));
	alkeep=(double *)calloc(pop->np,sizeof(double));
	crit=lik(pop); /*Calcul de Mf-1 */

	for(nit=0;nit<500;nit++)
	{
/* alkeep stores the initial frequencies of the elementary protocols
*/
		for(j=0;j<pop->np;j++)
			alkeep[j]=pop->freq[j];
		crit2=crit;
		if(ir==(-1)) /* First time through */
		{
			nnul=0; /* Computing the number of nul frequencies */
			for(j=0;j<pop->np;j++)
			{
				if(pop->freq[j]<FREQMIN) nnul++;
			}
		}
/* Compute the projected gradient (direction)*/
		in=project_grad(pop,allprot,gal,ir,&vmgal,nnul);
		if((in<0) | (gal[in]==0))
		{
			free(gal);
			free(alkeep);
			return nnul;
		}
/* Compute the optimal change along the gradient dmax (lambda)*/
		dmax=-pop->freq[in]/gal[in];
		imax=in;
		if((ir!=(-1)) & (gal[ir]<=0))
		{
			free(gal);
			free(alkeep);
			return nnul;
		}
		ir=-1;
		for(j=0;j<pop->np;j++)
		{
			if((gal[j]<0) & (fabs(pop->freq[j]/gal[j])<dmax))
			{
				dmax=-pop->freq[j]/gal[j];
				imax=j;
			}
		}
/* Updates the frequencies as beta*(t,j)=beta(t,j)+dmax*gal(j)
Compute the likelihood for the updated vector of frequencies
*/
		sal=0;
		for(j=0;j<pop->np;j++)
		{
			pop->freq[j]=pop->freq[j]+dmax*gal[j];
			sal=sal+pop->freq[j];
		}
		sal1=sal-pop->freq[imax];
		pop->freq[imax]=0;
		crit1=lik(pop);
/* If better likelihood, update alkeep (pop->freq has been updated already)
*/
		if(crit1>crit)
		{
			for(j=0;j<pop->np;j++) alkeep[j]=pop->freq[j];
			crit=crit1;
			nnul++;
		}
		else
		{
			for(j=0;j<pop->np;j++) pop->freq[j]=alkeep[j];
/* If not, try dmax*ro for ro=1/2**k until dmax*ro too small
*/
			iopt=optim_lambda(pop,gal,alkeep,&crit,dmax,crit2,vmgal);
			if(iopt==1)
			{
				ir=takeout_k(pop,allprot,&nnul);
/* ir=-1 if for all j, matrace(q[j]) <= ndim*(ndim+1)/2, ie optimal protocol
=> in that case exit the subroutine
*/
				if(ir==(-1))
				{
					free(gal);
					free(alkeep);
					return nnul;
				}
			}
		}
	}
	free(gal);
	free(alkeep);
	return -10;
}


void index_limit(POPROT *pop)
{
	return;
}

// ******************************************************************************************************************

int main_opt(PROTOC *allprot, POPROT *pop, IntegerVector nprot)
{
  int nit=0,ifin=0,nstop,ntest;
	double cri;

	cri=lik(pop);

	for(nit=0;nit<100;nit++)
	{

  nstop=0;
do
{
  ifin=ajout(allprot,pop,nprot[0]);
  nstop++;
} while ((cri<(-1e10))& (nstop<100));

if((nstop>99)& (cri<(-1e10)))

{
  return -10;
}
if(ifin==1) break; /* Si ifin=1, on sort de la boucle sur nit */
  ntest=doptimal(pop,allprot);
  if(ntest<0)
    return -100;
  tassement(pop);

  if(pop->np>(pop->ndim*(pop->ndim+1)/2))
  {
    index_limit(pop);
    return -1;
  }
  cri=lik(pop); /* Updates matrices */
    }
return nit;
}

// ******************************************************************************************************************

/******************************************************************************
  Fonctions calculatoires
	matrace:
	lik:determinant de hh et inverse de hh
*******************************************************************************/

double matrace(int iprot,PROTOC *prot,POPROT *pop)
{/* Computes tr(Mf(iprot)*Mf-1(poprot))
WARNING pop->finv must be computed beforehand (with a call to lik)*/
	int i,j,ij;
	double xcal=0;

	for(i=0;i<(prot->ndim);i++)
	{
		for(j=0;j<i;j++)
		{

			ij=i*(i+1)/2+j;
			xcal=xcal+2*(prot[iprot].fisher[ij])*(pop->finv[ij]);
			}
		ij=i*(i+3)/2;
		xcal=xcal+(prot[iprot].fisher[ij])*(pop->finv[ij]);
		}
	return xcal;
}

double lik(POPROT *pop)
{/* Function to compute the inverse of the Fisher information matrix of the
population protocol pop. Returns the determinant of the matrix, also stored in
pop->det. The Fisher  information matrix of the population protocol is stored
in pop->fisher and the inverse is stored in pop->finv
*/
	int ndim=pop->ndim;

	int ncase=(int)(ndim*(ndim+1)/2);

	int i,j,jj,ifail;
	int *indx;
	double *col;
	matrix *xa;
	double cri;

	indx=(int *)calloc(ndim,sizeof(int));
	col=(double *)calloc(ndim,sizeof(double));
	xa=matrix_create(ndim,ndim);


	for(i=0;i<ncase;i++)
	{
		pop->fisher[i]=0;
		for(j=0;j<pop->np;j++)
		{
pop->fisher[i]=pop->fisher[i]+pop->freq[j]*pop->pind[j].fisher[i];
		}
	}
	jj=0;
	for(i=0;i<ndim;i++)
	{
		for(j=0;j<=i;j++)
		{
			xa->m[i][j]=pop->fisher[jj];
			xa->m[j][i]=pop->fisher[jj];

			jj++;
		}
	}

	ifail=ludcmp(xa->m,ndim,indx,&cri);

	if(ifail==1) return CRI_MAX;
	for(i=0;i<ndim;i++) cri=cri*xa->m[i][i];
	pop->det=cri;
	for(i=0;i<ndim;i++)
	{
		for(j=0;j<ndim;j++) col[j]=0.0;
		col[i]=1.0;
		lubksb(xa->m,ndim,indx,col);
		for(j=0;j<ndim;j++)
		{
			jj=i*(i+1)/2+j;
			pop->finv[jj]=col[j]; /*=xb[j][i] using another matrix*/
			}
		}
	/* desallocation */
	matrix_destroy(xa);
	free(indx);
	return cri;
}

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii>=0)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

int ludcmp(double **a, int n,int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=(double *)calloc(n,sizeof(double));
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)

			if ((temp=fabs(a[i][j])) > big) big=temp;

		if (big == 0.0) {

			return 1;}
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != (n-1)) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	free(vv);
	return 0;
}

/* Random number generator */
double randu(double x)
{	double xrand;
	xrand=rand();
	return xrand;
}
double gauss1(long *idum)
{
/* 	float ran1(long *idum);
*/	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return (double)gset;
	}
}

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

/******************************************************************************
 FedorovInit_R
*******************************************************************************/
// [[Rcpp::export]]

List FedorovWynnAlgorithm_Rcpp( List protocols,
                                IntegerVector ndimen,
                                IntegerVector nbprot,
                                IntegerVector numprot,
                                NumericVector freq,
                                IntegerVector nbdata,
                                NumericVector vectps,
                                NumericVector fisher,
                                IntegerVector error,
                                IntegerVector protdep,
                                NumericVector freqdep)
{

IntegerVector nb_protocols = as<IntegerVector>(protocols[0]);
IntegerVector nb_times = as<IntegerVector>(protocols[1]);
IntegerVector nb_dimensions = as<IntegerVector>(protocols[2]);
IntegerVector total_cost = as<IntegerVector>(ndimen);
NumericMatrix samplingTimes = as<NumericMatrix >(protocols[4]);
NumericMatrix fisher_matrices = as<NumericMatrix>(protocols[5]);


// ***************************************************************************************

// ***************
// init parameters
// ***************

IntegerVector nprot = as<IntegerVector>(nb_protocols);
IntegerVector ndim = as<IntegerVector>(nb_dimensions);
int Ntot = total_cost[2];

// matrix for the optimal sampling times
NumericMatrix optimal_sampling_times(nprot[0],nb_times[0]);

int i,j,ij,ik,iprot,icas;

// *******************
//protocols
// *******************

// Vector of individual protocols
PROTOC *allprot;
allprot = (PROTOC *)calloc(nprot[0], sizeof(PROTOC));

// Population protocol
POPROT *pop;

// *********************************************
// sampling times
// allprot[i].tps = sampling time i
// allprot[i].tps[j] = element j of sampling time i
// *********************************************
// fisher matrices
// allprot[iprot].fisher[k]
// *********************************************

for(iprot=0;iprot<nprot[0];iprot++){

  PROTOC_alloc(&allprot[iprot],nb_times[0],ndim[0]);

  NumericVector row_matrix_sampling_times = samplingTimes.row(iprot);

  NumericVector fisher_matrix = fisher_matrices.row(iprot);

  for ( j=0; j< allprot[iprot].ntps;j++){
  allprot[iprot].tps[j] = row_matrix_sampling_times[j];
  }

  icas=0;

  for (ij=0;ij<ndim[0];ij++) {
    for (ik=0;ik<=ij;ik++) {
      allprot[iprot].fisher[icas]=fisher_matrix[icas]*Ntot/nb_times[0];
      icas++;
    }
  }
}

pop = initprot(allprot, protdep, freqdep);

// *********************************************
// optimisation des protocoles
// nok
// *********************************************

int nok = 0;
double sumf = 0;

nok = main_opt(allprot,pop,nprot);

  if(nok>=0) {

    sumf=0;

    for(i=0;i<pop->np;i++) {

      numprot[i]=pop->num[i]+1;
      nbdata[i]=pop->pind[i].ntps;

      freq[i]=pop->freq[i]*Ntot/(nb_times[0]);

      sumf=sumf+freq[i];

      freqdep[i]=pop->freq[i];
    }

    for(i=0;i<pop->np;i++) {
      freq[i]=freq[i]/sumf;

      for(j=0;j<pop->pind[i].ntps;j++) {
        vectps[j]=pop->pind[i].tps[j];

optimal_sampling_times(i,j) = pop->pind[i].tps[j];

      }
    }
    for(i=0;i<(pop->ndim*(pop->ndim+1)/2);i++) {
      fisher[i]=pop->fisher[i];
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("freq") = freq,
    Rcpp::Named("optimal_sampling_times") = optimal_sampling_times,
    Rcpp::Named("fisher") = fisher,
    Rcpp::Named("numprot") = numprot);


} // end FedorovInit_R
'
FedorovWynnAlgorithm_Rcpp = inline::cxxfunction( signature( protocols_input = "list",
                                                            ndimen_input = "integer",
                                                            nbprot_input = "integer",
                                                            numprot_input = "integer",
                                                            freq_input = "numeric",
                                                            nbdata_input = "integer",
                                                            vectps_input = "numeric",
                                                            fisher_input = "numeric",
                                                            nok_input = "integer",
                                                            protdep_input = "integer",
                                                            freqdep_input = "numeric"),
                                                 plugin = "Rcpp",
                                                 incl = incltxtFedorovWynnAlgorithm,
                                                 body = '
          List protocols = Rcpp::as<List>(protocols_input);
          IntegerVector ndimen  = Rcpp::as<IntegerVector>(ndimen_input);
          IntegerVector nbprot  = Rcpp::as<IntegerVector>(nbprot_input);
          IntegerVector numprot = Rcpp::as<IntegerVector>(numprot_input);
          NumericVector freq = Rcpp::as<NumericVector>(freq_input);
          IntegerVector nbdata = Rcpp::as<IntegerVector>(nbdata_input);
          NumericVector vectps = Rcpp::as<NumericVector>(vectps_input);
          NumericVector fisher = Rcpp::as<NumericVector>(fisher_input);
          IntegerVector nok = Rcpp::as<IntegerVector>(nok_input);
          IntegerVector protdep = Rcpp::as<IntegerVector>(protdep_input);
          NumericVector freqdep = Rcpp::as<NumericVector>(freqdep_input);

          return Rcpp::wrap( FedorovWynnAlgorithm_Rcpp(  protocols, ndimen, nbprot, numprot, freq, nbdata, vectps, fisher, nok, protdep, freqdep ) );')

output = FedorovWynnAlgorithm_Rcpp( protocols_input,  ndimen_input, nbprot_input,
                                    numprot_input, freq_input, nbdata_input,
                                    vectps_input, fisher_input, nok_input,
                                    protdep_input, freqdep_input )

return( output )

}


# ==============================================================================
#' FedorovWynnAlgorithm: S7 class
#' @description Implements the Fedorov-Wynn algorithm for optimal design.
#' @title FedorovWynnAlgorithm
#' @inheritParams Optimization
#' @param elementaryProtocols       List of elementary protocols.
#' @param numberOfSubjects          Total number of subjects.
#' @param proportionsOfSubjects     Initial proportion vector for each protocol.
#' @param showProcess               Logical; print iteration progress.
#' @param FedorovWynnAlgorithmOutputs Output list from the Rcpp routine.
#' @param modelCovariatesEquation modelCovariatesEquation
#' @param numberOfOccasions numberOfOccasions
#' @export

FedorovWynnAlgorithm = new_class( "FedorovWynnAlgorithm",
                                  package    = "PFIM",
                                  parent     = Optimization,
                                  properties = list(
                                    elementaryProtocols         = new_property( class_list,    default = list()  ),
                                    numberOfSubjects            = new_property( class_vector,  default = 0.0     ),
                                    proportionsOfSubjects       = new_property( class_vector,  default = 0.0     ),
                                    showProcess                 = new_property( class_logical, default = FALSE   ),
                                    FedorovWynnAlgorithmOutputs = new_property( class_list,    default = list()  )
                                  )
)

plotFrequenciesFedorovWynnAlgorithm = new_generic(
  "plotFrequenciesFedorovWynnAlgorithm",
  c( "optimization", "optimizationAlgorithm" )
)

# ==============================================================================
#' optimizeDesign: FedorovWynnAlgorithm
#' @name optimizeDesign
#' @param optimizationObject   An \code{Optimization} object.
#' @param optimizationAlgorithm A \code{FedorovWynnAlgorithm} object.
#' @return The \code{optimizationObject} with updated slots.
#' @export
# ==============================================================================

method( optimizeDesign, list( Optimization, FedorovWynnAlgorithm ) ) =
  function( optimizationObject, optimizationAlgorithm ) {

    p                          = prop( optimizationObject, "optimizerParameters" )
    initialSamplings           = p$elementaryProtocols
    totalNumberOfIndividuals   = p$numberOfSubjects
    proportionsOfSubjects      = p$proportionsOfSubjects
    # Flatten initial elementary protocols into a single vector for row matching
    initialElementaryProtocols = unlist( initialSamplings )
    totalCost                  = sum( lengths( initialSamplings ) * totalNumberOfIndividuals )

    design      = prop( optimizationObject, "designs" )[[ 1L ]]
    optimalDesign = design
    designName  = prop( design, "name" )

    fimsFromConstraints = generateFimsFromConstraints( optimizationObject )

    # Reduce list-of-row-vectors into matrices expected by the C routine
    samplingsForFedorovWynn = reduce( fimsFromConstraints$samplingsForFedorovWynnAlgo[[ designName ]], rbind )
    fisherMatrices          = reduce( fimsFromConstraints$listFimsAlgoFW[[            designName ]], rbind )
    listArms                = fimsFromConstraints$listArms[[ designName ]]

    # ── Assemble input structure ───────────────────────────────────────────────
    nProtocols = nrow( samplingsForFedorovWynn )
    nTimes     = ncol( samplingsForFedorovWynn )
    ndimFim    = fimsFromConstraints$dimFim

    elementaryProtocolsFW = list(
      numberOfprotocols = nProtocols,
      numberOfTimes     = nTimes,
      nbOfDimensions    = ndimFim,
      totalCost         = totalCost,
      samplingTimes     = samplingsForFedorovWynn,
      fisherMatrices    = fisherMatrices
    )

    ndimen  = c( nProtocols, ndimFim, totalCost )
    npInit  = nProtocols
    nBuf    = nProtocols * 2L
    numprot = rep( 0, nBuf )
    freq    = rep( 0, nBuf )
    nbdata  = rep( 0, nBuf )
    vectps  = rep( 0, nBuf )
    fisher  = rep( 0, nProtocols * ( nProtocols + 1L ) / 2L )
    nok     = 0L

    # Identify which rows of the sampling grid match the user-specified protocols
    indexElemProt = which(
      rowSums( samplingsForFedorovWynn ==
                 initialElementaryProtocols[ col( samplingsForFedorovWynn ) ] ) == nTimes
    )
    nElemProt = length( indexElemProt )

    zeprot = c( nElemProt, indexElemProt )
    zefreq = rep( 0, nBuf )
    zefreq[ seq_len( length( proportionsOfSubjects ) ) ] = proportionsOfSubjects

    # ── Run C routine ──────────────────────────────────────────────────────────
    output = FedorovWynnAlgorithm_Rcpp(
      elementaryProtocolsFW, ndimen, npInit,
      numprot, freq, nbdata, vectps, fisher, nok,
      zeprot, zefreq
    )

    # ── Parse outputs ──────────────────────────────────────────────────────────
    if ( !any( rowSums( output$optimal_sampling_times ) > 0 ) )
      stop( "FedorovWynn algorithm did not converge." )

    optimalFrequencies       = output$freq[   output$freq    > 0 ]
    indexOptimalArms         = output$numprot[ output$numprot > 0 ]
    listArms                 = listArms[ indexOptimalArms ]
    numberOfIndividuals      = as.double( p$numberOfSubjects * optimalFrequencies )  # FIX: ensure double type

    prop( optimizationAlgorithm, "FedorovWynnAlgorithmOutputs" ) = list(
      listArms            = listArms,
      optimalFrequencies  = optimalFrequencies,
      numberOfIndividuals = numberOfIndividuals
    )

    fim         = prop( optimizationObject, "fim" )
    optimalArms = setOptimalArms( fim, optimizationAlgorithm )
    prop( optimalDesign, "arms" ) = map( optimalArms, ~ .x$arm )

    prop( optimizationObject, "optimisationDesign" ) = list(
      evaluationInitialDesign = run( .buildEvalFromOptimization( optimizationObject, design        ) ),
      evaluationOptimalDesign = run( .buildEvalFromOptimization( optimizationObject, optimalDesign ) )
    )
    prop( optimizationObject, "optimisationAlgorithmOutputs" ) = list(
      optimizationAlgorithm = optimizationAlgorithm,
      optimalArms           = optimalArms,
      frequencies           = optimalFrequencies
    )

    optimizationObject
  }

# ==============================================================================
#' plotFrequenciesFedorovWynnAlgorithm
#' @name plotFrequenciesFedorovWynnAlgorithm
#' @param optimization          An \code{Optimization} object.
#' @param optimizationAlgorithm A \code{FedorovWynnAlgorithm} object.
#' @return A ggplot bar chart of optimal protocol frequencies.
#' @export
# ==============================================================================

method( plotFrequenciesFedorovWynnAlgorithm,
        list( Optimization, FedorovWynnAlgorithm ) ) =
  function( optimization, optimizationAlgorithm ) {

    out  = prop( optimization, "optimisationAlgorithmOutputs" )
    data = data.frame(
      arm       = map_chr( out$optimalArms, ~ prop( .x$arm, "name" ) ),
      frequency = out$frequencies
    )
    .algoBars( data, x = "arm", y = "frequency", xlab = "Arm", ylab = "Frequency" )
  }

# ==============================================================================
#' constraintsTableForReport: FedorovWynnAlgorithm
#' @name constraintsTableForReport
#' @param optimizationAlgorithm A \code{FedorovWynnAlgorithm} object.
#' @param arms List of arms.
#' @export
# ==============================================================================

method( constraintsTableForReport, FedorovWynnAlgorithm ) =
  function( optimizationAlgorithm, arms ) {
    .constraintsKbl(
      arms, optimizationAlgorithm,
      c( "Arms name", "Number of subjects", "Outcome",
         "Initial samplings", "Fixed times",
         "Number of samplings optimisable", "Dose constraints" )
    )
  }
