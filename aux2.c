#include "SW.h"
#include "math.h"
/**
   \section{AverageInteraction}
   \subsection{Description}
   calculates the average interaction between neighbors.
   \subsection{Input parameters}
   \begin{itemize}
   \item[J] the interaction array. J.p[i][j] is the interaction between point i
       and its jth neighbour.
   \end{itemize}
   \subsection{Return value}
   returns ``the Characteristic Distance''.
   \subsection{file}
   aux2.c
**/
float AverageInteraction( RaggedArray J )
{
int     i,k;
int     sum = 0;
float   J1 = 0;
   
   for(i = 0;  i < J.n; i++)
     for(k = 0; k<J.c[i]; k++) {
       J1 += J.p[i][k];
       sum++;
     }
   
   J1 /= ((float) sum);
   
   return(J1);
}

/**
   \section{GlobalCorrelation}
   \subsection{Description}
   Builds the correlations array. If two points belongs to the same cluster,
   one is added to the coresponding matrix element.
   \subsection{Input parameters}
   \begin{itemize}
   \item[CorrN] the correlations array. CorrN.p[i][j] is the number of times
       vertex i and vertex NK.p[i][j] were in the same cluster so far.
   \item[NK] nearest neighbours array. NK.p[i][j] is the j-th neighbour of
       vertex i.
   \item[Block] the cluster each vertex belongs to.
       vertex i belongs to cluster number Block[i].
   \end{itemize}
   \subsection{Output parameters}
   \begin{itemize}
   \item[CorrN] the updated correlation array.
   \end{itemize}
   \subsection{file}
   aux2.c
**/
void  GlobalCorrelation( UIRaggedArray CorrN, UIRaggedArray NK,
			 unsigned int *Block )
{
   int  i,k;

   for(i = 0; i < NK.n; i++)
      for(k = 0; k<NK.c[i]; k++)
	 if( Block[i] == Block[ NK.p[i][k] ] ) CorrN.p[i][k]++;
           
   return;
}


/*void  Cal_Corr(RaggedArray CorrN,RaggedArray MFTc,float width,float T)
{
   int  i,k;

   for(i = 0; i < CorrN.n; i++){
	   for(k = 0; k<CorrN.c[i]; k++){
		  CorrN.p[i][k] = 1/(exp((T - MFTc.p[i][k])/width) + 1);
	   }
   }
			           
   return;
}*/

void  cal_effectiveJ(RaggedArray MFtotJ,RaggedArray MFTc,RaggedArray MFwid_end,
					 UIRaggedArray MFmin,RaggedArray J,UIRaggedArray NK,UIRaggedArray KN,float T)
{
	int i,k,k1,Q,min_idx,max_idx;
	float wid;

	Q = IGetParam( "PottsSpins" );


	for (i=0 ; i < MFtotJ.n; i++){
		for( k = MFtotJ.c[i]-1; NK.p[i][k]>i && k>=0; k-- ) {
			if (MFmin.p[i][k] == i){
				min_idx = i;
				max_idx = NK.p[i][k];
			}
			else{
				min_idx = NK.p[i][k];
				max_idx = i;
			}
			if (MFTc.p[i][k]<T){
				MFtotJ.p[i][k] = J.p[i][k];
				for ( k1=0 ; k1 < MFtotJ.c[min_idx] ; k1++){
					if (/*MFwid_end.p[min_idx][k1]>T &&*/ NK.p[min_idx][k1]!= max_idx){
						wid = 0.5*(MFwid_end.p[min_idx][k1]-MFTc.p[min_idx][k1]);
						MFtotJ.p[i][k] += J.p[min_idx][k1]*
								(1/(exp((T-(MFwid_end.p[min_idx][k1]+wid))/(wid/4)) + 1));
						/*exp(MFwid_end.p[min_idx][k1]/T)/(exp(MFwid_end.p[min_idx][k1]/T) + (Q-1));*/
					}
				}
				MFtotJ.p[ NK.p[i][k] ][ KN.p[i][k] ] = MFtotJ.p[i][k];
			}
		}
	}
}
			

int  cal_order_parameter(float *s,float *s_old,RaggedArray J,UIRaggedArray NK,float T)
{
	int i,k,flag,count;
	float Q,Jeff;
	float *tmp;

	Q = (float)IGetParam( "PottsSpins" );
	count=0;
	flag=1;
	while (flag){
		flag=0;
		count++;
		for(i = 0; i < J.n; i++){
			if (s_old[i]>0.0001){
				Jeff = 0;
				for(k = 0; k<J.c[i]; k++)
					Jeff += s_old[ NK.p[i][k] ]*J.p[i][k];
				if (Jeff>0)
					s[i] = (1 - exp(-Jeff/T))/(1 + (Q-1)*exp(-Jeff/T));
				else
					s[i] = 0;

				if (s_old[i]-s[i]>0.01){
					flag=1;
				}
			}
			else
				s[i] = 0;
		}
		tmp = s_old;
		s_old=s;
		s=tmp;
	}

	return count;
}


float cal_free_energy(float *s_old,float *alpha,RaggedArray J,UIRaggedArray NK,float T)
{
	int i,k;
	float Q,fe;

	if (T>0){
		for(i = 0; i < NK.n; i++){
			alpha[i] = 0;
			for(k=0 ; k<NK.c[i] ; k++)
				alpha[i] += J.p[i][k]*s_old[ NK.p[i][k] ]/T;
		}


		Q = (float)IGetParam( "PottsSpins" );
		fe=0;
		for(i = 0; i < NK.n; i++){
			if (alpha[i]<50)
				fe += T*((alpha[i]*exp(alpha[i]))/(Q-1+exp(alpha[i])) - log(Q-1+exp(alpha[i])));
			for(k = NK.c[i]-1; NK.p[i][k]>i && k>=0 ; k--)
				fe -= J.p[i][k]*((Q-1)*exp(-(alpha[i]+alpha[ NK.p[i][k] ]))+1)/(((Q-1)*exp(-alpha[i])+1)*((Q-1)*exp(-alpha[ NK.p[i][k] ])+1));
		}
	}
	else
	{
		fe=0;
		for(i = 0; i < NK.n; i++){
			for(k = NK.c[i]-1; NK.p[i][k]>i && k>=0 ; k--)
				fe -= J.p[i][k];
		}
	}
	return fe;
}

int cal_local_sus(float *local_sus,float *div_alpha,float *div_alpha_old,
				  float *s_old,float *alpha,RaggedArray Coef,RaggedArray J,UIRaggedArray NK,float T)
{
	int i,k,flag,count;
	float Q,Jeff;
	float *tmp;

	if (T==0){
		for (i=0 ; i<NK.n ; i++)
			local_sus[i] = 0;
		return 0;
	}

	
	for(i = 0; i < J.n; i++){
		alpha[i] = 0;
		for(k = 0; k<J.c[i]; k++)
			alpha[i] += s_old[ NK.p[i][k] ]*J.p[i][k];
		alpha[i] /= T;
	}

	Q = (float)IGetParam( "PottsSpins" );

	
	for(i = 0; i < J.n; i++)
		for(k = 0; k<J.c[i]; k++)
			Coef.p[i][k] = J.p[i][k]*Q*exp(-alpha[ NK.p[i][k] ])/
					(((Q-1)*exp(-alpha[ NK.p[i][k] ])+1)*((Q-1)*exp(-alpha[ NK.p[i][k] ])+1));


	count=0;
	flag=1;

	while (flag){
		flag=0;
		count++;
		for(i = 0; i < NK.n; i++){
			if (div_alpha_old[i]>0.0001){
				div_alpha[i] = 1;
				for(k = 0; k<NK.c[i]; k++)
					div_alpha[i] += Coef.p[i][k]*div_alpha_old[ NK.p[i][k] ];

				div_alpha[i] /= T;
				if (div_alpha_old[i]-div_alpha[i]>0.01 || div_alpha_old[i]-div_alpha[i]<-0.01){
					flag=1;
				}
			}
		}
		tmp = div_alpha_old;
		div_alpha_old=div_alpha;
		div_alpha=tmp;
	}

	for(i = 0; i < NK.n; i++)
		local_sus[i] = T*div_alpha_old[i]*(Q-1)*exp(-alpha[i])/
					(((Q-1)*exp(-alpha[i])+1)*((Q-1)*exp(-alpha[i])+1));


	return count;
}


void update_max_local_sus(float *CritTemp,float *MaxSus,float *LocalSus,float T,int N)
{
	int i;

	for(i = 0; i < N; i++){
		if (LocalSus[i]>MaxSus[i] && T-CritTemp[i]<0.05){
				MaxSus[i] = LocalSus[i];
				CritTemp[i] = T;
			}
	}
	return;
}


void update_max_slope(float *crit_temp,float *max_slope,float *s_old,float *s_old_tmp,float T,int N)
{
	int i;

	for (i=0 ; i<N ; i++){
		if ((s_old_tmp[i]-s_old[i])>max_slope[i]){
			max_slope[i] = (s_old_tmp[i]-s_old[i]);
			crit_temp[i] = T;
		}
	}
	return;
}

void find_crit_temp(float *crit_temp,float *s_cur,float *s_next,float T,int N)
{
	int i;

	for (i=0 ; i<N ; i++){
		if (s_cur[i]>=0.5){/* && s_next[i]<0.5){*/
			crit_temp[i] = T;
		}
	}
	return;
}

void  Cal_Corr(RaggedArray CorrN,RaggedArray CritTemp,UIRaggedArray NK,UIRaggedArray NK_mst,float *s,RaggedArray J,float T,float threshold)
{
   int  i,k;
   float Q,pair_corr,mf_corr;

   Q = (float)IGetParam( "PottsSpins" );

   for(i = 0; i < CorrN.n; i++)
	   for(k = 0; k<CorrN.c[i]; k++){
		   pair_corr = (1-exp(-J.p[i][k]/T))/(1+(Q-1)*exp(-J.p[i][k]/T));
		   if (NK_mst.p[i][k]<CorrN.n){
			    mf_corr = s[i]*s[ NK.p[i][k] ];
				CorrN.p[i][k] = (mf_corr>pair_corr)? mf_corr : pair_corr;	
				/*CorrN.p[i][k] = mf_corr;	*/
		   }
		   else{
				CorrN.p[i][k] = pair_corr;
		   }
		   if (CorrN.p[i][k]>threshold)
				CritTemp.p[i][k] = T;
	   }

/*			CorrN.p[i][k] = 1/Q + s[i]*s[ NK.p[i][k] ]*(Q-1)/Q;*/
			           
   return;
}
/*
void  Cal_Corr(RaggedArray CorrN,RaggedArray MFTc,RaggedArray MFwid,RaggedArray MFtotJ,RaggedArray J,float T)
{
   int  i,k,Q;

   Q = IGetParam( "PottsSpins" );

   for(i = 0; i < CorrN.n; i++){
	   for(k = 0; k<CorrN.c[i]; k++){
			if (T<=MFTc.p[i][k])
					CorrN.p[i][k] = 1;
			else{
				if (T<MFTc.p[i][k]+2*MFwid.p[i][k] && MFtotJ.p[i][k]>MFTc.p[i][k])
					CorrN.p[i][k] = 1/(exp((T-MFtotJ.p[i][k])/((MFtotJ.p[i][k]-MFTc.p[i][k])/4)) + 1);					
					/*CorrN.p[i][k] = ((Q-1)*(1/(exp((T-MFtotJ.p[i][k])/((MFtotJ.p[i][k]-MFTc.p[i][k])/4)) + 1))+1)/Q;*/
					/*CorrN.p[i][k]=exp(MFtotJ.p[i][k]/T)/(exp(MFtotJ.p[i][k]/T) + (Q-1));
				/*else
					CorrN.p[i][k]=exp(J.p[i][k]/T)/(exp(J.p[i][k]/T) + (Q-1));*/
/*			}
			if (CorrN.p[i][k]< exp(J.p[i][k]/T)/(exp(J.p[i][k]/T) + (Q-1)))
					CorrN.p[i][k]=exp(J.p[i][k]/T)/(exp(J.p[i][k]/T) + (Q-1));

			
	   }
   }
			           
   return;
}*/


void  GlobalCorrelationNbyN( unsigned int *CorrNbyN, UIRaggedArray NK,
			 unsigned int *Block )
{
  int  i,j,cur;

   cur=0;
   for(i = 0; i < (NK.n-1); i++)
     for(j= i+1; j < NK.n; j++)
       if( Block[i] == Block[j] ) ++CorrNbyN[cur++];
       else cur++;
           
   return;
}

/**
   \section{FourPointCorrelation}
   \subsection{Description}
   Builds the four point correlations array. If four points belongs to
   the same cluster, one is added to the coresponding matrix element.
   \subsection{Input parameters}
   \begin{itemize}
   \item[FPCorr] the four point correlation array. CorrN.p[i][j].p[k][l]
   is the number of times vertices i, NK.p[i][j], k and NK.p[k][l]
   were in the same cluster so far.
   \item[NK] nearest neighbours array. NK.p[i][j] is the j-th neighbour of
       vertex i.
   \item[Block] the cluster each vertex belongs to.
       vertex i belongs to cluster number Block[i].
   \end{itemize}
   \subsection{Output parameters}
   \begin{itemize}
   \item[FPCorr] the updated four point correlation array.
   \end{itemize}
   \subsection{file}
   aux2.c
**/
void  FourPointCorrelation( RARaggedArray FPCorr, UIRaggedArray NK,
			    unsigned int *Block) {
  int i,j,k;
  int i1,j1,k1;
  
  for(i = 0; i < NK.n; i++) 
    for(k = 0; k<NK.c[i]; k++) 
      if(NK.p[i][k]>i) 
        for(i1 = 0; i1 < NK.n; i1++)
          for(k1 = 0;  k1<NK.c[i1]; k1++) 
            if (NK.p[i1][k1]>i1) 
              if(Block[i] == Block[NK.p[i][k]] 
		 && Block[i1]==Block[NK.p[i1][k1]])
                (FPCorr.p[i][k]).p[i1][k1]++;
}



/* an auxiliary function for magnetization() and OrderClusterSize(). */
/* returns the opposite result to this of uicomp in edge.c           */
int uicompare(const void *i, const void *j)
{ return (int)(*((unsigned int*)j)) - (int)(*((unsigned int*)i)); }

/**
   \section{Magnetization}
   \subsection{Description}
   Obtain how many points have each spin value. order the groups in
   discending order of size. Calculate the magnetization for each spin
   color as \[ m_q = {Q N_q - N \over N (Q-1)}, \] where $N_q$ in the
   number of points with spin color $q$.
   \subsection{Input parameters}
   \begin{itemize}
   \item[N] number of points.
   \item[Q] number of spin values (colors).
   \item[nc] number of clusters.
   \item[*ClusterSize] cluster sizes.
   \item[$N_q$] previously allocated workspace of 
       size $>=$ Q*sizeof(unsigned int).
   \end{itemize}
   \subsection{Output parameters}
   \begin{itemize}
   \item[mag] the magnetization vector.
   \end{itemize}
   \subsection{Auxiliary function}
   int uicompare(const void *i, const void *j)
   \subsection{file}
   aux2.c
**/
float Magnetization( int N, int Q, int nc, unsigned int *ClusterSize,
		     float* mag, unsigned int *N_q )
{
   int k, q;

   memset( N_q, 0, Q*sizeof(unsigned int) );
   for(k = 0; k < nc; k++)
      N_q[ IRAND(Q) ] += ClusterSize[k];
   qsort(N_q,Q,sizeof(unsigned int),uicompare);
   for(q = 0; q < Q; q++)
      mag[q] = (float)( (int)(Q * N_q[q] - N) ) / (N*(Q - 1.0));

   return mag[0];
}

/**
   \section{OrderClusterSize}
   \subsection{Description}
   Order cluster sizes list in discending order..
   \subsection{Input parameters}
   \begin{itemize}
   \item[nc] number of clusters.
   \item[*ClusterSize] cluster sizes.
   \end{itemize}
   \subsection{Output parameters}
   \begin{itemize}
   \item[*ClusterSize] ordered cluster sizes.
   \end{itemize}
   \subsection{Auxiliary function}
   int uicompare(const void *i, const void *j)
   \subsection{file}
   aux2.c
**/
void OrderClusterSize( int nc, unsigned int *ClusterSize )
{
   qsort(ClusterSize,nc,sizeof(unsigned int),uicompare);   
}   

/**
   \section{ClusterAverage}
   \subsection{Description}
   Calculate cluster size averages. Returns the mean value and
   variance of the larger, second larger, .., smaller cluster sizes.
   \subsection{Input parameters}
   \begin{itemize}
   \item[nc] number of clusters.
   \item[*Size1] cluster size conumlant. The i-th element is the comulant of
          the i-th larger cluster size.
   \item[*size2] cluster size$^2$ comulant.
   \end{itemize}
   \subsection{Output parameters}
   \begin{itemize}
   \item[*Size1] sizes mean value.
   \item[*Size2] sizes variance.
   \end{itemize}
   \subsection{file}
   aux2.c
**/
void ClusterAverage(int ncy, int N, float *Size1, float *Size2)
{
   int i;
   /* assume Size1[] is in descending order */
   for(i = 0; i<N && Size1[i]>0; i++) {
      Size1[i] /= (float)ncy;
      Size2[i] /= (float)ncy;
      Size2[i] -= (Size1[i] * Size1[i]);
   }
}

/**
   \section{Susceptibility}
   \subsection{Description}
   Calculate magnatizations mean value and susceptibilities.
   \subsection{Input parameters}
   \begin{itemize}
   \item[Q] number of spin values.
   \item[ncy] number of SW sweeps performed = numbers of samples taken.
   \item[*M1] magnetizations comulant.
   \item[*M2] magnetizations$^2$ comulant.
   \end{itemize}
   \subsection{Output parameters}
   \begin{itemize}
   \item[*M1] magnetizations mean value.
   \item[*xi] susceptibilities.
   \end{itemize}
   \subsection{file}
   aux2.c
**/
void Susceptibility( int Q, int ncy, float* M1, float* M2, float* xi )
{
   int q;
   for(q=0;q<Q;q++) {
      M1[q] /= (float)ncy;
      M2[q] /= (float)ncy;          
      xi[q] = (M2[q] -  M1[q] * M1[q]); 
   }
}
/*
int  Thresholding( int ncy, float T,   
		     float *crit_temp, UIRaggedArray NK,CRaggedArray Bond, 
			 unsigned int *Block, unsigned int *ClusterSize,
                  unsigned int *OldBlock, int *n_cols, unsigned int *ws)

{
   int i,k;
   int   nc;

   if (!GetParam("NoDendrogram")) {
     for(i=0; i< NK.n; i++)
		 for(k = 0; k<NK.c[i];k++){
		   if (crit_temp[i] >=T  && crit_temp[NK.p[i][k]] >=T){
			 if (OldBlock[i] == OldBlock[NK.p[i][k]]) 
				Bond.p[i][k] = 1;
			 else {
				Bond.p[i][k] = 0;
				(*n_cols)++;
			}
		   }
		   else
				Bond.p[i][k] = 0;
     }
   } else {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
		   if (crit_temp[i] >=T  && crit_temp[NK.p[i][k]] >=T){
	   Bond.p[i][k] = 1;
	 }
	 else Bond.p[i][k] = 0;

   }

   nc = Coarsening(Bond,Block,NK,ClusterSize,ws);
   OrderingClusters(NK.n,nc,Block,ClusterSize,ws);
   
   return(nc);
} */   

int  Thresholding(int ncy, float T,   
                  RaggedArray CritTemp, UIRaggedArray NK, CRaggedArray Bond,
                  unsigned int *Block, unsigned int *ClusterSize,
                  unsigned int *OldBlock, int *n_cols, unsigned int *ws)
{
   int i,k,nc;
   float  th;

   if (!GetParam("NoDendrogram")) {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];  k++)
         if(T <= CritTemp.p[i][k]) {
           if ( OldBlock[i] == OldBlock[ NK.p[i][k] ] ) {
             Bond.p[i][k] = 1;
           }
           else {
             Bond.p[i][k] = 0;
             (*n_cols)++;
           }
         }
         else
	   Bond.p[i][k] = 0;
   } else {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];  k++)
         if (T <= CritTemp.p[i][k]) 
	   Bond.p[i][k] = 1;
	 else 
	   Bond.p[i][k] = 0;
   }

   nc = Coarsening(Bond,Block,NK,ClusterSize,ws);
   OrderingClusters(NK.n,nc,Block,ClusterSize,ws);
   
   return(nc);
}



/* ---------------------------------------------------------------------- */
/*   We perform a directed growth in order to  determin  a partition      */
/*   Two points belongs to the same cluster if they have a common         */
/*   predesessor (see Fukunaga).                                          */
/*   The predessessor of each point is the point with maximal correlation */
/*   among its neighboors.                                                */
/*   The threshold is the number for which we consider that the           */
/*   correlation is "zero"                                                */
/* ---------------------------------------------------------------------- */
/*   actually, here each point is joined to the cluster of                */
/*   its predessessor   -- Guy */

int  DirectedGrowth( int ncy, float threshold,   
		     RaggedArray Corr, UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, unsigned int *thBlock,
		     unsigned int *ClusterSize,unsigned int *thClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws )
{
   int i,k,Q;
   float th;
   float dh;
   int   nc;
   float max;
   int kmax;

   /*Q = IGetParam( "PottsSpins" );*/

   th = threshold * ncy;
/*   dh = ((Q-1)*th/20+1)/Q;*/
   dh = th/20;
   
   if (!GetParam("NoDendrogram")) {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
         if ((th < Corr.p[i][k]) && 
	     (thOldBlock[i] == thOldBlock[NK.p[i][k]])) {
	   Bond.p[i][k] = 1;
	 }
	 else Bond.p[i][k] = 0;

/*	 nc = Coarsening(Bond,thBlock,NK,thClusterSize,ws);*/

     for(i=0; i< NK.n; i++){
       max = 0;
       for(k = 0; k<NK.c[i];  k++) 
	 if( max < Corr.p[i][k] ){
	   max = Corr.p[i][k];
	   kmax = k;
	 }
       if(max > dh){
         if ( dgOldBlock[i] == dgOldBlock[ NK.p[i][kmax] ] ) {
	   Bond.p[i][kmax] = 1;
	   Bond.p[ NK.p[i][kmax] ][ KN.p[i][kmax] ] = 1;
	 }
       }
     }
   } else {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
         if (th < Corr.p[i][k]) {
	   Bond.p[i][k] = 1;
	 }
	 else Bond.p[i][k] = 0;

/*     nc = Coarsening(Bond,thBlock,NK,thClusterSize,ws);*/
     
     for(i=0; i< NK.n; i++){
       max = 0;
       for(k = 0; k<NK.c[i];  k++) 
	 if( max < Corr.p[i][k] ){
	   max = Corr.p[i][k];
	   kmax = k;
	 }
       if(max > dh){
	 Bond.p[i][kmax] = 1;
	 Bond.p[ NK.p[i][kmax] ][ KN.p[i][kmax] ] = 1;
       }
     }
   }

/*   nc = DGCoarsening(Bond,Block,thBlock,NK,ClusterSize,thClusterSize,ws);*/
   nc = Coarsening(Bond,Block,NK,ClusterSize,ws);
   OrderingClusters(NK.n,nc,Block,ClusterSize,ws);
   
   return(nc);
}    


int  nDirectedGrowth(float T,float Tstep,int Q,RaggedArray MFTc, UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, 
		     unsigned int *ClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws )
{
   int i,k,k_max;
   float th,ln,mftc_max;
   float dh;
   int   nc;
   int   max1, kmax1, max2, kmax2;

	ln = log(Q-1);
    th = T*2*ln*(Q-1)/(Q-2);

   
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
         if ((T < MFTc.p[i][k]) && 
	     (thOldBlock[i] == thOldBlock[NK.p[i][k]])) {
	   Bond.p[i][k] = 1;
	 }
	 else Bond.p[i][k] = 0;

     nc = Coarsening(Bond,Block,NK,ClusterSize,ws);

     for(i=0; i< NK.n; i++){
		 if (ClusterSize[Block[i]]<5){
			max1 = 0;
			max2 = 0;
			for(k = 0; k<NK.c[i];  k++) {
				if( max1 < MFTc.p[i][k] && ClusterSize[ Block[NK.p[i][k]] ]>=5 ){
					max1 = MFTc.p[i][k];
					kmax1 = k;

				}
				if( max2 < MFTc.p[i][k] ){
					max2 = MFTc.p[i][k];
					kmax2 = k;
				}
			}
			if(T-Tstep>0 && max1 > T-Tstep){
				if ( dgOldBlock[i] == dgOldBlock[ NK.p[i][kmax1] ] ) {
					Bond.p[i][kmax1] = 1;
					Bond.p[ NK.p[i][kmax1] ][ KN.p[i][kmax1] ] = 1;
				}
			}
			else {
				if(T-Tstep>0 && max2 > T-Tstep){
					if ( dgOldBlock[i] == dgOldBlock[ NK.p[i][kmax2] ] ) {
						Bond.p[i][kmax2] = 1;
						Bond.p[ NK.p[i][kmax2] ][ KN.p[i][kmax2] ] = 1;
					}
				}
			}
		 }
     }

	nc = Coarsening(Bond,Block,NK,ClusterSize,ws);
	OrderingClusters(NK.n,nc,Block,ClusterSize,ws);
   
   return(nc);
}    


void FindTc(int ncy,float threshold,unsigned int *CorrNbyN,float *TcNbyN,float T,int N)
{
	float th;
	int i;

	th = 0.5 * ncy;

	for(i=0; i< N*(N-1)/2; i++){
         if ( th>CorrNbyN[i] && TcNbyN[i]==0 )
			 TcNbyN[i] = T;
		 CorrNbyN[i]=0;
	}
}
												

double SingleLinkage(int ncy,UIRaggedArray CorrN,UIRaggedArray NK)
{

	int i,j,k,edges_counter,Q;
	int node1,node2,node_counter;
	double *edges_corr;
	double G;
	unsigned int *edges_1;
	unsigned int *edges_2;
	double *dist_vec;
	unsigned int    *nodes;	/* auxiliar */
	unsigned int    *indx;	/* auxiliar */
	double mean;
	double var;


	Q = IGetParam( "PottsSpins" );
	edges_corr = InitDVector( IGetParam( "NumberOfEdges" ) );
	edges_1 = InitUIVector( IGetParam( "NumberOfEdges" ) );
	edges_2 = InitUIVector( IGetParam( "NumberOfEdges" ) );
	dist_vec = InitDVector( NK.n-1 );
	edges_counter=0;

	for(i=0; i< NK.n; i++){
		for(k = NK.c[i]-1; NK.p[i][k]>i && k>=0 ;  k--) {
			G = ((Q-1)*((double)CorrN.p[i][k]/(double)ncy)+1)/Q;
			if (G<1)
				edges_corr[edges_counter] = -(G*log(G)/log((double)2) + (1-G)*log((1-G)/(Q-1))/log((double)2));
			else
				edges_corr[edges_counter] = -(G*log(G)/log((double)2));
/*			edges_corr[edges_counter] = -log(G);*/
			edges_1[edges_counter] = i;
			edges_2[edges_counter] = NK.p[i][k];
			edges_counter++;
		}
	}

	indx = InitUIVector(edges_counter);
	DSortIndex(edges_counter,edges_corr,indx);
	nodes = InitUIVector(NK.n);
	for (i=0 ; i<NK.n ; i++)
		nodes[i] = i;

	for (node_counter=0 , i=0 ; i < edges_counter ; i++){
		node1 = nodes[edges_1[indx[i]]];
		node2 = nodes[edges_2[indx[i]]];
		if (node1!=node2){
			dist_vec[node_counter] = edges_corr[indx[i]];
			for (j=0 ; j<NK.n ; j++)
				if (nodes[j]==node1 || nodes[j]==node2)
					nodes[j] =  NK.n+node_counter;
			node_counter ++;
		}
	}

	for (mean=0 , i=0 ; i< NK.n-1 ; i++){
		if (dist_vec[i]<INF){
			mean += dist_vec[i];
			if (dist_vec[i]>0)
				mean=mean;
		}
			
	}
	mean = mean/(NK.n-1);
	for (var=0 ,  i=0 ; i< NK.n-1 ; i++)
		if (dist_vec[i]<INF)
			var += (dist_vec[i]-mean)*(dist_vec[i]-mean);
	var = var/(NK.n-2);

	return var;
}


int  nThresholding(float T,float Tstep,int Q,RaggedArray MFTc, UIRaggedArray NK, CRaggedArray Bond,
                  unsigned int *Block, unsigned int *ClusterSize,
                  unsigned int *OldBlock, int *n_cols, unsigned int *ws)
{
   int i,k,nc,k_max;
   float th,ln,mftc_max;

	ln = log(Q-1);
   th = T*2*ln*(Q-1)/(Q-2);

   if (!GetParam("NoDendrogram")) {
	for(i=0; i< NK.n; i++)
	{
		for(k = 0,mftc_max=0; k<NK.c[i];  k++)
		{
         if(T <= MFTc.p[i][k]) {
           if ( OldBlock[i] == OldBlock[ NK.p[i][k] ] ) {
             Bond.p[i][k] = 1;
           }
           else {
             Bond.p[i][k] = 0;
             (*n_cols)++;
           }
         }
         else
		 {
			Bond.p[i][k] = 0;
			if (MFTc.p[i][k]>mftc_max)
			{
				mftc_max = MFTc.p[i][k];
				k_max = k;
			}
		 }
		}
		if((T-Tstep)>0 && (T-Tstep)<mftc_max) 
			Bond.p[i][k_max] = 1;
	}
	} else {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];  k++)
         if(T <= MFTc.p[i][k])
	   Bond.p[i][k] = 1;
	 else 
	   Bond.p[i][k] = 0;
   }

   nc = Coarsening(Bond,Block,NK,ClusterSize,ws);
   OrderingClusters(NK.n,nc,Block,ClusterSize,ws);
   
   return(nc);
}


float SolveNumeric(int i, RaggedArray J, RaggedArray MFTc,RaggedArray MFwid)
{
	int k,Q;
	float x,y,step,cpl;

	Q = IGetParam( "PottsSpins" );
	cpl = (Q-2)/(2*(Q-1)*log(Q-1));

	x=0;
	for (step=1 ; step>0.01 ; step=step/2){
		x +=step;
		for (y=0, k=0 ; k<MFTc.c[i] ; k++){
			if (x<=MFTc.p[i][k])
				y+=J.p[i][k];
			else
				/*if (x<MFTc.p[i][k]+0.05)
					y+=J.p[i][k]*(1-(x-MFTc.p[i][k])/0.05);*/
				if (MFwid.p[i][k])
					y+=J.p[i][k]*(1/(exp((x-(MFTc.p[i][k]+MFwid.p[i][k]))/(MFwid.p[i][k]/4)) + 1));
					/*y+=J.p[i][k]*exp(MFTc.p[i][k]/(cpl*x))/(exp(MFTc.p[i][k]/(cpl*x)) + (Q-1));
					/*y+=J.p[i][k]*((Q-1)*(1/(exp((x-(MFTc.p[i][k]+MFwid.p[i][k]))/(MFwid.p[i][k]/4)) + 1))+1)/Q;*/
					/*y+=J.p[i][k]*exp((MFTc.p[i][k]-x)/MFwid.p[i][k]);*/
			
		}
		while (x<y){
			x +=step;
			for (y=0, k=0 ; k<MFTc.c[i] ; k++){
				if (x<=MFTc.p[i][k])
					y+=J.p[i][k];
				else
					/*if (x<MFTc.p[i][k]+0.05)
						y+=J.p[i][k]*(1-(x-MFTc.p[i][k])/0.05);*/
					if (MFwid.p[i][k])
						y+=J.p[i][k]*(1/(exp((x-(MFTc.p[i][k]+MFwid.p[i][k]))/(MFwid.p[i][k]/4)) + 1));
					/*	y+=J.p[i][k]*exp(MFTc.p[i][k]/(cpl*x))/(exp(MFTc.p[i][k]/(cpl*x)) + (Q-1));
						/*y+=J.p[i][k]*((Q-1)*(1/(exp((x-(MFTc.p[i][k]+MFwid.p[i][k]))/(MFwid.p[i][k]/4)) + 1))+1)/Q;*/
						/*y+=J.p[i][k]*exp((MFTc.p[i][k]-x)/MFwid.p[i][k]);*/
			}
		}
		x -=step;
	}

	return x;
}


/*
int  SteepestDescent( int ncy, float T,   
		     float *crit_temp, UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, unsigned int *thBlock,
		     unsigned int *ClusterSize,unsigned int *thClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws )
{
   int i,k,min_sz;
   float th;
   float dh;
   int   nc;
   float max;
   int kmax;

   min_sz = IGetParam("DirectedGrowthMinSize");

   if (!GetParam("NoDendrogram")) {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
         if (crit_temp[i] >=T  && crit_temp[NK.p[i][k]] >=T &&
	     (thOldBlock[i] == thOldBlock[NK.p[i][k]])) {
	   Bond.p[i][k] = 1;
	 }
	 else Bond.p[i][k] = 0;

	 nc = Coarsening(Bond,thBlock,NK,thClusterSize,ws);

     for(i=0; i< NK.n; i++){
		 if (thClusterSize[ thBlock[i] ]<min_sz){
       max = 0;
       for(k = 0; k<NK.c[i];  k++) 
	 if( max < crit_temp[ NK.p[i][k] ] && crit_temp[ NK.p[i][k] ]>=crit_temp[i]){
	   max = crit_temp[ NK.p[i][k] ];
	   kmax = k;
	 }
       if(max > 0){
	 Bond.p[i][kmax] = 1;
	 Bond.p[ NK.p[i][kmax] ][ KN.p[i][kmax] ] = 1;
       }
		 }
     }
   } else {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
		   if (crit_temp[i] >=T  && crit_temp[NK.p[i][k]] >=T){
	   Bond.p[i][k] = 1;
	 }
	 else Bond.p[i][k] = 0;

     nc = Coarsening(Bond,thBlock,NK,thClusterSize,ws);
     
     for(i=0; i< NK.n; i++){
		 if (thClusterSize[ thBlock[i] ]<min_sz){

       max = 0;
       for(k = 0; k<NK.c[i];  k++) 
	 if( max < crit_temp[ NK.p[i][k] ] && crit_temp[ NK.p[i][k] ]>=crit_temp[i]){
	   max = crit_temp[ NK.p[i][k] ];
	   kmax = k;
	 }
       if(max > 0){
	 Bond.p[i][kmax] = 1;
	 Bond.p[ NK.p[i][kmax] ][ KN.p[i][kmax] ] = 1;
       }
		 }
     }
   }

   nc = DGCoarsening(Bond,Block,thBlock,NK,ClusterSize,thClusterSize,ws);
   OrderingClusters(NK.n,nc,Block,ClusterSize,ws);
   
   return(nc);
}    */






int  SteepestDescent( int ncy, float T, RaggedArray CritTemp,
		     UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, unsigned int *thBlock,
		     unsigned int *ClusterSize,unsigned int *thClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws )
{
   int i,k,k1,min_sz;
   float th;
   float dh;
   int   nc;
   float max;
   int kmax,max_cnt;

   min_sz = IGetParam("DirectedGrowthMinSize");


   if (!GetParam("NoDendrogram")) {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
         if ((T <= CritTemp.p[i][k]) && 
	     (thOldBlock[i] == thOldBlock[NK.p[i][k]])) {
	   Bond.p[i][k] = 1;
	 }
	 else Bond.p[i][k] = 0;

	   
	 nc = Coarsening(Bond,thBlock,NK,thClusterSize,ws);

     for(i=0; i< NK.n; i++){
       max = 0;
	   max_cnt = 0;
       for(k = 0; k<NK.c[i];  k++) 
	 if( max <= CritTemp.p[i][k] ){
	   max = CritTemp.p[i][k];
	   kmax = k;
	   max_cnt++;
	 }
       if(max > 0 /*&& max_cnt<NK.c[i]*/){
         if ( dgOldBlock[i] == dgOldBlock[ NK.p[i][kmax] ] ) {
	   Bond.p[i][kmax] = 1;
	   Bond.p[ NK.p[i][kmax] ][ KN.p[i][kmax] ] = 1;
	 }
       }
	 }
/*
     for(i=0; i< NK.n; i++){
		 if (thClusterSize[ thBlock[i] ]<min_sz){
			 min=INF;
			 for (k=0, kmin=0; k<NK.c[i] ; k++){
				 if (crit_temp[ NK.p[i][k] ]>crit_temp[i] && (crit_temp[ NK.p[i][k] ]-crit_temp[i])/J.p[i][k]<min){
					 min = (crit_temp[ NK.p[i][k] ]-crit_temp[i])/J.p[i][k];
					 kmin = k;
				 }
			 }
			 if (min<INF){
				 Bond.p[i][kmin] = 1;
				 Bond.p[ NK.p[i][kmin] ][ KN.p[i][kmin] ] = 1;
			}
		 }
	 }*/
				 
   } else {
     for(i=0; i< NK.n; i++)
       for(k = 0; k<NK.c[i];k++)
         if (T <= CritTemp.p[i][k])
		   Bond.p[i][k] = 1;
		 else Bond.p[i][k] = 0;

     nc = Coarsening(Bond,thBlock,NK,thClusterSize,ws);
     
     for(i=0; i< NK.n; i++){
       max = 0;
	   max_cnt = 0;
       for(k = 0; k<NK.c[i];  k++) 
	 if( max <= CritTemp.p[i][k] ){
	   max = CritTemp.p[i][k];
	   kmax = k;
	   max_cnt++;
	 }
       if(max > 0 /*&& max_cnt<NK.c[i]*/){
	   Bond.p[i][kmax] = 1;
	   Bond.p[ NK.p[i][kmax] ][ KN.p[i][kmax] ] = 1;
	 }
	 }
   }

   nc = DGCoarsening(Bond,Block,thBlock,NK,ClusterSize,thClusterSize,ws);
        /*for(i=0; i< NK.n; i++){
		 if (ClusterSize[ Block[i] ]<min_sz){
			 for (k=0 ; k<NK.c[i] ; k++){
				 if (crit_temp[ NK.p[i][k] ]==crit_temp[i] && ClusterSize[ Block[ NK.p[i][k] ] ]>=min_sz){
					for (k1=0 ; k1<NK.c[ NK.p[i][k] ] ; k1++){
						if (crit_temp[ NK.p[i][k] ] <T  || crit_temp[ NK.p[ NK.p[i][k] ][k1] ] <T){
							Bond.p[ NK.p[i][k] ][ k1 ] = 0;
							Bond.p[ NK.p[ NK.p[i][k] ][k1] ][ KN.p[ NK.p[i][k] ][k1] ] = 0;
						}
					}
				 }
			 }
		 }
		}

   nc = DGCoarsening(Bond,Block,thBlock,NK,ClusterSize,thClusterSize,ws);*/
   OrderingClusters(NK.n,nc,Block,ClusterSize,ws);
   
   return(nc);
}    




