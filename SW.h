#ifndef SW_H
#define SW_H

#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include "RaggedArray.h"
#include "param.h"
#include "utilities.h"
#include "timer.h"

/* edge.c */
void ReadEdgeFile( int N, UIRaggedArray *nk );
UIRaggedArray InvertEdges(UIRaggedArray NK);
/*UIRaggedArray knn(int N, int D, double** X );*/
UIRaggedArray knn( int N, int D, double** X, char **mst_edges );
UIRaggedArray MstEdges(UIRaggedArray NK,char **mst_edges);
UIRaggedArray nknn( int N, int D, double** X );
void mstree(int N, int D,double** X, unsigned int** edg );
void OrderEdges( UIRaggedArray *nk );

/* distance.c */
double Distance(int D, double *X, double *Y);
double Distance_Linf(int D, double *X, double *Y);
RaggedArray EdgeDistance( int D, UIRaggedArray NK, double **X );
void DistanceToInteraction( RaggedArray M, UIRaggedArray NK, UIRaggedArray KN );

/* aux.c */
void InitialSpinConfig(int N, unsigned int *Spin, int Q);
void NewSpinConfig( int N, unsigned int  *Spin, unsigned int *Block, 
		    int NBlk, int Q, unsigned int *NewSpinValue);
void NewSpinConfigField( int N, unsigned int  *Spin, unsigned int *Block, 
		    int NBlk, int Q, unsigned int *NewSpinValue, 
                    unsigned int *StateN);
void DeletionProbabilities( float T, RaggedArray J, RaggedArray P );
int SetBond( RaggedArray P, unsigned int *Spin, CRaggedArray Bond,
             UIRaggedArray NK, UIRaggedArray KN);
int  Coarsening(CRaggedArray Bond, unsigned int *Block,
                UIRaggedArray NK, unsigned int *ClusterSize,
		unsigned int* Stack );
int  DGCoarsening(CRaggedArray Bond, unsigned int *Block,unsigned int *thBlock,
                UIRaggedArray NK, unsigned int *ClusterSize,unsigned int *thClusterSize,
		unsigned int* Stack );
void  OrderingClusters( int N, int nblock, unsigned int *Block,
			unsigned int *Size, unsigned int *Indx );
void CheckParam();
void DefaultParam();
double Energy( unsigned int* Block, RaggedArray J, UIRaggedArray NK  );

/* aux2.c */
float AverageInteraction( RaggedArray J );
void  GlobalCorrelation( UIRaggedArray CorrN, UIRaggedArray NK,
			 unsigned int *Block );
void  GlobalCorrelationNbyN( unsigned int *CorrNbyN, UIRaggedArray NK,
			 unsigned int *Block );
void  FourPointCorrelation( RARaggedArray FPCorr, UIRaggedArray NK,
			    unsigned int *Block);
float Magnetization( int N, int Q, int nc, unsigned int *ClusterSize,
		     float* mag, unsigned int *N_spin );
void OrderClusterSize( int nc, unsigned int *ClusterSize );
void ClusterAverage(int ncy, int N, float *Size1, float *Size2);
void Susceptibility( int Q, int ncy, float* M1, float* M2, float* xi );
int  Thresholding(int ncy, float T,   
                  RaggedArray CritTemp, UIRaggedArray NK, CRaggedArray Bond,
                  unsigned int *Block, unsigned int *ClusterSize,
                  unsigned int *OldBlock, int *n_cols, unsigned int *ws);
int  DirectedGrowth( int ncy, float threshold,   
		     RaggedArray Corr, UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, unsigned int *thBlock,
		     unsigned int *ClusterSize, unsigned int *thClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws );
int  nDirectedGrowth(float T,float Tstep,int Q,RaggedArray MFTc, UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, unsigned int *ClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws );

void FindTc(int ncy,float threshold,unsigned int *CorrNbyN,float *TcNbyN,float T,int N);
double SingleLinkage(int ncy,UIRaggedArray CorrN,UIRaggedArray NK);

int  nThresholding(float T,float Tstep,int Q,RaggedArray MFTc, UIRaggedArray NK, CRaggedArray Bond,
                  unsigned int *Block, unsigned int *ClusterSize,
                  unsigned int *OldBlock, int *n_cols, unsigned int *ws);
void  Cal_Corr(RaggedArray CorrN,RaggedArray CritTemp,UIRaggedArray NK,UIRaggedArray NK_mst,float *s,RaggedArray J,float T,float threshold);
void  cal_effectiveJ(RaggedArray MFtotJ,RaggedArray MFTc,RaggedArray MFwid_end,
					 UIRaggedArray MFmin,RaggedArray J,UIRaggedArray NK,UIRaggedArray KN,float T);

float SolveNumeric(int i, RaggedArray J, RaggedArray MFTc,RaggedArray MFwid);

int  cal_order_parameter(float *s,float *s_old,RaggedArray J,UIRaggedArray NK,float T);
void update_max_slope(float *crit_temp,float *max_slope,float *s_old,float *s_old_tmp,float T,int N);
float cal_free_energy(float *s_old,float *alpha,RaggedArray J,UIRaggedArray NK,float T);


/*int  SteepestDescent( int ncy, float T,   
		     float *crit_temp, UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, unsigned int *thBlock,
		     unsigned int *ClusterSize,unsigned int *thClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws );*/
int  SteepestDescent( int ncy, float T, RaggedArray CritTemp,
		     UIRaggedArray NK, UIRaggedArray KN,
		     CRaggedArray Bond, unsigned int *Block, unsigned int *thBlock,
		     unsigned int *ClusterSize,unsigned int *thClusterSize, unsigned int *dgOldBlock,
		     unsigned int *thOldBlock, unsigned int *ws );

void find_crit_temp(float *crit_temp,float *s_cur,float *s_next,float T,int N);


int cal_local_sus(float *local_sus,float *div_alpha,float *div_alpha_old,
				  float *s_old,float *alpha,RaggedArray Coef,RaggedArray J,UIRaggedArray NK,float T);
void update_max_local_sus(float *CritTemp,float *MaxSus,float *LocalSus,float T,int N);


/* io.c */
void ReadData(int N, int D, double **X);
void ReadPrevTempFiles(unsigned int *thOldBlock, 
		       unsigned int *dgOldBlock,int N);
void PrintAverages(int nT, float T, double e, double ee, float nc, float *Size); 
void PrintMagnet( int nT, float T, float* M, float *xi );
void PrintCurBlock( int nT, int cyc, unsigned int *Block, int N);
void PrintStateN( int nT, float T, unsigned int *StateN, int ncy, int N, int Q);
void PrintCorrNbyN( int nT, float T, unsigned int *CorrNbyN, int N, int ncy);
void PrintTcNbyN(float *TcNbyN,int N);
void NNPrintCorrN( int nT, float T, float dend_var, UIRaggedArray CorrN, int ncy,
		   UIRaggedArray NK );
void nNNPrintCorrN( int nT, float T, float dend_var, RaggedArray CorrN, int ncy,
		   UIRaggedArray NK );
void PrintFPointCorr( int nT, float T, RARaggedArray FPCorr,
		      UIRaggedArray NK, int ncy );
void PrintFPSum(int nT, RaggedArray J, UIRaggedArray Corr, RARaggedArray FPCorr,
		UIRaggedArray NK, int ncy);
void PrintSizes(char* pe,int nth,int nT,float T,int nc,
		unsigned int *Size,int n_cols);
void WriteLabels(char* pe,int nth,int nT,float T,int N,unsigned int *Block);
void ReadParam( int argc, char* argv[] );
void PrintParam();
void WriteEdges( UIRaggedArray nk );
void WriteJs( UIRaggedArray nk, RaggedArray J );
void PrintTime( float T );
void ReadCorrN(FILE* corr_file,RaggedArray CorrN, int ncy,UIRaggedArray NK,UIRaggedArray KN);

void NNPrintMFTc(int n,RaggedArray MFTc,UIRaggedArray NK );
void NNPrintOrderParam(int nT,float T,int iter_num,float fe,float *s,int N);
void NNPrintLocalSus(int nT,float T,int iter_num,float *sus,int N);

#endif

/**
\section{Environment Parameters}
To make the program simpler and clearer we use environment parameters, 
in which we keep data and flags to be use by different functions,
without the need to pass them as arguments. This also allows an easy way
to add more flags and options to the program wuthout littering the code
too much. The parameters used in the programm so far are:
\begin{description}
\item[AverageInteraction] $\:$\\ the average interaction between neighbours.
\item[ClustersReported] $\:$\\ how many cluster sizes are reported.
\item[ClusterMinSizeReported] $\:$\\ minimal size of cluster to be reported.
     all clusters above or equal this size will be reported, even if their
     number exceeds the \texttt{ClustersReported} value.
\item[CharDist] $\:$\\ characteristic distance between neighbours
        (used for calculating the interaction).
\item[DataFile] $\:$\\ the file containing the coordinate of the point, or
        the distances matrix.
\item[DataIsInteraction] $\:$\\ The distances obtained should be cosidered as
        the interaction.
\item[Dimensions] $\:$\\ dimention of the varctors describing the points.
        if 0, the data is expected to be a distances matrix.
\item[DirectedGrowth] $\:$\\ Report cluster obtained from the correlations
       by directed growth.
\item[DataIsMatrix] $\:$\\ the data file is organized as a distance matrix
     and not as a list.
\item[EdgeFile] $\:$\\ file containing edges to be added to the list of
       nearest neighbours.
\item[FourPointCorr] $\:$\\ Sample and report four point correlations.
\item[ForceNN] $\:$\\ Use this (false) number of nearest neighbours for
      calculating the interaction.
\item[ForceRandomSeed] $\:$\\ Random seed to start the program with, if we
     want this run to be identicle to a previous one. If not set, the seed
     is taken as the clock value.
\item[ForceChD] $\:$\\Use this (false) characteristic distance for calculating
      the interaction.
\item[InfMetric] $\:$\\ Use infinity metric to calculate distances.
\item[KNearestNeighbours] $\:$\\ maximal number of nearest neighbours
        (used in the knn algorithm).
\item[Lambda] $\:$\\ $\lambda$ value to be used if \texttt{UseZ} is set.
\item[MinTemp] $\:$\\ The lowest temperature to use..
\item[MaxTemp] $\:$\\ The maximal temp. to use.
\item[MSTree] $\:$\\ add the edges of the minimal spanning tree.
\item[NumberOfPoints] $\:$\\ number of points.
\item[NumberOfEdges] $\:$\\ the total number of edges.
\item[NearestNeighbrs] $\:$\\ average number of nearest neighbours
         (used for calculating the interaction).
\item[OutFile] $\:$\\ prefix for the output files.
\item[PottsSpins] $\:$\\ The number of different spin values.
\item[PrevTempFile] $\:$\\ prefix for former output files, used to obtain
     initial states for the threshold and directed growth calculations.
\item[RandomSeed] $\:$\\ The random seed with which the program has started.
\item[SusceptColors] $\:$\\ number of susciptibilities from the susceptibility
          vector to report.
\item[SaveAverages] $\:$\\ sample and report data on Swensen-Wang averages. 
\item[SaveSuscept] $\:$\\ sample and report magnetization and
        susceptibility data.
\item[SWCycles] $\:$\\ number Swensen-Wang sweeps. 
\item[SWFraction] $\:$\\ the fraction SW sweeps for which averages are
      calculated. The first (1-SWfract)*cyc sweeps are discarded.
\item[Threshold] $\:$\\ Report cluster obtained from the correlations
       by thresholding.
\item[ThresholdTheta] $\:$\\ The threshold on correlations between neighbours,
     above which they are assumed to belong to the same cluster.
\item[ThresholdMin]
\item[ThresholdMax]
\item[ThresholdStep] $\:$\\ If those variables are set, the programs report the
      results obtained for different thresholds, taken from min to max 
      values with the step indicated.
\item[TempStep] $\:$\\ the steps in temperature from min to max values in which
      the simulation should be ran.
\item[Timing] $\:$\\ report the time required by the programm for each step.
\item[UseZ] $\:$\\ Scale the interaction according to a given $\lambda$.
\item[WriteLabels] $\:$\\ Write the label of the cluster each spin belongs to.
\item[WriteCorFile] $\:$\\ Write the correlations between each
        pair of neighbours.

\end{description}
**/


