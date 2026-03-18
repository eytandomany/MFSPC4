#include "SW.h"
#include "math.h"

int main(int argc, char *argv[])
{
int  N;                  /* Number of points                          */
int      Q;	         /* Number of Potts Spins; Si = 0,...,Q-1     */
float    T;              /* Temperature                               */
float    Tmin,Tmax,dT;   /* Temp. min, max, and step                  */
UIRaggedArray NK;        /* Containns the neighboors of each point.   */
                         /* NK.n is the number of points and NK.c[i]  */
                         /* is the number of neighbours of point i.   */
                         /* NK[i][j], j=1..NK.c[i] are the labels of  */
                         /* neighboring points for each i=1..NK.n     */
UIRaggedArray KN;        /* KN[i][k] = m, means that point i is the   */
                         /* m-th neighbor of point j = N[i][k]        */
UIRaggedArray NK_mst;
char **mst_edges=NULL;
int      D;              /* Dimension of the point vectors, D=0 means */
                         /* that the distances are provided (instead  */
                         /* the coordinates of the vectors).          */
double **X;              /* if D=0 X[i][j] is the distance between    */
                         /* points i and j. if else, it is the j-th   */
                         /* coordinate of the  i-th point.            */
RaggedArray J;           /* J.p[i][j] is the interaction between the  */
                         /* i-th spin and its j-th neighbour          */

CRaggedArray Bond;       /* Bond[i][k] takes value 1 (0) if bond      */
	     	         /* between spins i and its k-th neighbor is  */
                         /* frozen (deleted).                         */
unsigned int *ClusterSize;/* ClusterSize[n] contains the number of    */
                         /* points belonging to cluster n in the      */
                         /* current configuration. The clusters are   */
                         /* ordered fron the biggest (n=1) to the     */
                         /* smallest                                  */
unsigned int   *Block;   /* Block[i] Is the number of the cluster to  */
		         /* which Spin[i] belongs                     */
unsigned int *thClusterSize;/* ClusterSize[n] contains the number of    */
                         /* points belonging to cluster n in the      */
                         /* current configuration. The clusters are   */
                         /* ordered fron the biggest (n=1) to the     */
                         /* smallest                                  */
unsigned int   *thBlock;   /* Block[i] Is the number of the cluster to  */
		         /* which Spin[i] belongs                     */
unsigned int *UIWorkSpc; /* auxiliary work space                      */
unsigned int *dgOldBlock;/* dgOldBlock[i] Is the number of the cluster*/
			 /* to which Spin[i] belonged in the previous */
                         /* temperature after directed growth         */
unsigned int *thOldBlock;/* OldBlock[i] Is the number of the cluster  */
			 /* to which Spin[i] belonged in the previous */
                         /* temperature after thresholding            */
RaggedArray CorrN;     /* Two Points Correlations comulant          */
RaggedArray CritTemp;     /* Two Points Correlations comulant          */
float	 corr_th;
int      nc;             /* present number of clusters                */
int      nc1;            /* cluster number comulant                   */
int      nb;             /* present number of frozen bonds            */
int      nT;             /* Temp. step counter                         */
int n_cols;              /* N times connecting points from different   */
                         /* labels in previous temp. at theta stage    */
int      i;              /* auxiliary loop index                       */
char cflag = 1;			/* Continue run flag */

char pn[PARAM_MAX_FIELD];
FILE *corr_file;

float *s;
float *s_old;
float **s_final;
float *s_old_tmp;
float *crit_temp;
float *max_slope;
float *alpha;
int *iter_num;
int flag,count;
float tmpT,tmpFE;
float *FreeEnergy;
int numTstep;


/*RaggedArray Coef;
float *local_sus;
float *div_alpha;
float *div_alpha_old;
float *max_sus;
int iter_num2;*/


  printf("\nSPC v2.1.1 (15/9/2003) -- Mean Field approximation");   
  printf("\n-------------------------------------------------------------------\n");
  if(argc < 2){
    printf("Usage: %s [Parameter file]\n",argv[0]);
    exit(1);
  }
  DefaultParam();
  ReadParam(argc,argv);
  /* spelling mistakes I made. better be removed in the future... */
  if(!GetParam("Dimensions")) SetParam("Dimensions",GetParam("Dimentions"));
  if(GetParam("WriteLables")) SetParam("WriteLabels",NULL);
  /* .............................................................*/ 
  CheckParam();

  if( GetParam( "Timing" ) ) start_timer();
 
  if( GetParam( "ForceRandomSeed" ) ) i = IGetParam( "ForceRandomSeed" );
  else i = time(NULL) % RAND_MAX % INT_MAX;
  srand( i );
  ISetParam( "RandomSeed", i ); 

  N = IGetParam( "NumberOfPoints" );
  Q = IGetParam( "PottsSpins" );
  Tmin = FGetParam( "MinTemp" );
  Tmax = FGetParam( "MaxTemp" );
  dT = FGetParam( "TempStep" );/*10;*/
  D = IGetParam( "Dimensions" );
  corr_th = FGetParam("ThresholdTheta" );    

  if( D==0 ) X = InitDMatrix(N,N);
  else X = InitDMatrix(N,D);
  ReadData(N,D,X);

  NK.n = 0;
  if( GetParam( "KNearestNeighbours" )  || GetParam( "MSTree" ) ){
	  if ( GetParam( "MSTree" ) ){
		mst_edges = InitCMatrix(N,N);
		ResetCMatrix(mst_edges,N,N);
	  }

     NK = knn( N, D, X, mst_edges );
  }
  if( GetParam( "EdgeFile") ) ReadEdgeFile( N, &NK );
  OrderEdges( &NK ); /* Edges *must* be ordered when calling SetBond() */
  KN = InvertEdges( NK );
  WriteEdges( NK );
  J = EdgeDistance( D, NK, X );
  assure( IGetParam( "NumberOfEdges" ) > 0, "no edges" );

  NK_mst = MstEdges(NK,mst_edges);
  if (mst_edges){
	FreeCMatrix(mst_edges,N);
  }

  FreeDMatrix(X,N);

  if ( !GetParam( "DataIsInteraction" ) )
     DistanceToInteraction( J, NK, KN );
  FSetParam( "AverageInteraction", AverageInteraction(J) );

  PrintParam();

  if( GetParam( "Timing" ) ) PrintTime(-1.0);

    /* Memory allocations: */
  CorrN = InitRaggedArray(NK);
  CritTemp = InitRaggedArray(NK);
  Bond = InitCRaggedArray(NK);
  ClusterSize = InitUIVector(N);
  Block = InitUIVector(N);
  thClusterSize = InitUIVector(N);
  thBlock = InitUIVector(N);
  UIWorkSpc = InitUIVector((2*N>Q)?2*N:Q);  /* bounds on UIWorkSpc size:
                                               >=Q for magnetization
					       >=2N for OrderingClusters  */

  dgOldBlock = InitUIVector(N);  
  thOldBlock = InitUIVector(N);
  memset( dgOldBlock, 0, N*sizeof(unsigned int) );
  memset( thOldBlock, 0, N*sizeof(unsigned int) );
  if( GetParam( "PrevTempFile" ) )
     ReadPrevTempFiles(thOldBlock,dgOldBlock,N);


	 if (GetParam ("ClBuilder")){
		sprintf(pn,"%s%s",GetParam("OutFile"),".cor");
		corr_file = fopen( pn, "r");
		assure( corr_file, "correlation file" );
	 }

	if (GetParam("TempStepMul"))
		numTstep = log(Tmax/Tmin)/log(dT)+2;
	else
		numTstep = (Tmax-Tmin)/dT + 2;

	iter_num = InitUIVector(numTstep);
	FreeEnergy = InitVector(numTstep);
	for (i=0 ; i<numTstep ; i++){
		iter_num[i] = 0;
	}
	s_final = InitMatrix(numTstep,N);
	s = InitVector(N);
	s_old = InitVector(N);
	
	crit_temp = InitVector(N);
	max_slope = InitVector(N);
	s_old_tmp = InitVector(N);

	alpha = InitVector(N);

	/*local_sus = InitVector(N);
	div_alpha = InitVector(N);
	div_alpha_old = InitVector(N);
	max_sus = InitVector(N);*/

	for (i=0 ; i<N ; i++){
		s_old[i] = 1;
		crit_temp[i] = 0;
		max_slope[i] = 0;
/*		max_sus[i] = 0;
		div_alpha_old[i] = 10000;*/
	}



/*********************** T LOOP **********************/
  for(T = Tmin, nT = 0; T <= Tmax && cflag; nT++, T=GetParam("TempStepMul")?T*dT:T+dT ){



	 iter_num[nT] += cal_order_parameter(s,s_old,J,NK,T);
/*	 update_max_slope(crit_temp,max_slope,s_old,s_old_tmp,T,N);*/
	 FreeEnergy[nT] = cal_free_energy(s_old,alpha,J,NK,T);
	 /*NNPrintOrderParam(nT,T,iter_num[nT],FreeEnergy[nT],s_old,N);*/

	for (i=0 ; i<N ; i++){
		s_final[nT][i] = s_old[i];
		s_old_tmp[i] = s_old[i];
	}
	 count = 0;
	 flag = 1;
	 tmpT = T;
	 while (flag && nT-count-1>=0)
	 {
		tmpT = GetParam("TempStepMul")?tmpT/dT:tmpT-dT;
		iter_num[nT] += cal_order_parameter(s,s_old_tmp,J,NK,tmpT);
		tmpFE = cal_free_energy(s_old_tmp,alpha,J,NK,tmpT);
/*		NNPrintOrderParam(count,tmpT,iter_num[nT],tmpFE,s_old_tmp,N);*/
		if (FreeEnergy[nT-count-1]-tmpFE>0.01){
			flag = 1;
			FreeEnergy[nT-count-1] = tmpFE;
			for (i=0 ; i<N ; i++){
				s_final[nT-count-1][i] = s_old_tmp[i];
				s_old[i] = s_old_tmp[i];
			}
		}
		else{
			flag = 0;
			nT -=count;
			T = GetParam("TempStepMul")?tmpT*dT:tmpT+dT;
		}
		count++;
	 }
  }

/*********************** T LOOP **********************/
/*  ResetRaggedArray(CritTemp);
  for(T = Tmin, nT = 0; T <= Tmax && cflag; nT++, T=GetParam("TempStepMul")?T*dT:T+dT ){

     ResetRaggedArray(CorrN);
	 Cal_Corr(CorrN,CritTemp,NK,NK_mst,s_final[nT],J,T,corr_th);
	if (nT>0)
	  find_crit_temp(crit_temp,s_final[nT],s_final[nT+1],T,N);

	
	NNPrintOrderParam(nT,T,iter_num[nT],FreeEnergy[nT],s_final[nT],N);
  }


    NNPrintOrderParam(nT,T,0,0,crit_temp,N);*/

  
  for(T = Tmin, nT = 0; T <= Tmax && cflag; nT++, T=GetParam("TempStepMul")?T*dT:T+dT ){
	 
/* new version */
     /* Memory reset: */
     ResetCRaggedArray(Bond);
     ResetRaggedArray(CorrN);
	 Cal_Corr(CorrN,CritTemp,NK,NK_mst,s_final[nT],J,T,corr_th);
	 if( GetParam( "WriteCorFile" ) ){
		nNNPrintCorrN(nT,T,0,CorrN,1,NK);
	 }

     if( GetParam( "Threshold") ) {
	   n_cols=0;
	   nc = Thresholding(1,T,CritTemp,NK,Bond,Block,ClusterSize,
			     thOldBlock,&n_cols,UIWorkSpc);
	   PrintSizes(".th_",1,nT,T,nc,ClusterSize,n_cols);
	   if ( GetParam("WriteLabels") )
	      WriteLabels(".th_",1,nT,T,N,Block);                
	memcpy( thOldBlock, Block, N*sizeof(unsigned int) );
     }


     /* threshod + directed growth */
     if( GetParam( "DirectedGrowth" ) || GetParam( "StopRunAtBreak" ) ) {
	   /*nc = SteepestDescent(1,T,CritTemp,NK,KN,Bond,Block,thBlock,
						ClusterSize,thClusterSize,dgOldBlock,thOldBlock,UIWorkSpc);*/

	   nc = DirectedGrowth(1,corr_th,CorrN,NK,KN,Bond,Block,thBlock,
						ClusterSize,thClusterSize,dgOldBlock,thOldBlock,UIWorkSpc);
	   /* Notice that above thOldBlock is the new thBlocks but */
	   /* it is OK to use them */  
	   if ( IGetParam( "StopRunAtBreak" )) {
	     int srcnt;
	     for (srcnt=0;ClusterSize[srcnt]>=IGetParam( "StopRunClusterSize") 
		    && srcnt<nc;srcnt++);
	     cflag = (srcnt<IGetParam( "StopRunAtBreak" ));
	   }
	   PrintSizes(".dg_",1, nT,T,nc,ClusterSize,0);
	   if ( GetParam("WriteLabels") )
	      WriteLabels(".dg_",1,nT,T,N,Block);
	if( !GetParam( "Threshold") )
	   memcpy( thOldBlock, Block, N*sizeof(unsigned int) );
	memcpy( dgOldBlock, Block, N*sizeof(unsigned int) );
     }


  }   /******************** END OF T LOOP ***********************/


  FreeMatrix(s_final,numTstep);
  free(FreeEnergy);
  free(iter_num);

  free(s);
  free(s_old);
  free(s_old_tmp);
  free(crit_temp);
  free(max_slope);

  free(alpha);

/*  free(local_sus);
  free(div_alpha);
  free(div_alpha_old);
  FreeRaggedArray(Coef);
  free(max_sus);*/

  FreeRaggedArray(CorrN);
  FreeRaggedArray(CritTemp);
  FreeCRaggedArray(Bond);
  free(ClusterSize);
  free(Block);
  free(thClusterSize);
  free(thBlock);
  free(UIWorkSpc);

  FreeUIRaggedArray(NK);
  FreeUIRaggedArray(NK_mst);
  FreeUIRaggedArray(KN);
  FreeRaggedArray(J);
  free(dgOldBlock);
  free(thOldBlock);

  return(0);  
}
