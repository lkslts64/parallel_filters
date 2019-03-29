#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "defs.h"

;int main (int argc, char* argv[])
{
  int rank, size;
 
  MPI_Init (&argc, &argv);      /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
	srand(time(NULL));
	int numprocs=sqrt(size);
	//find north,east,west,south procs
	int i;
	int west,north,south,east;
	if((rank%numprocs)==0){
		west=MPI_PROC_NULL;
		east=rank+1;
		compass(&north,&south,rank,numprocs);
	}
	else if((rank%numprocs)==(numprocs-1)){
		west=rank-1;
		east=MPI_PROC_NULL;
		compass(&north,&south,rank,numprocs);
	}
	else{
		west=rank-1;
		east=rank+1;
		compass(&north,&south,rank,numprocs);
	}
	
	int botleftdiag=MPI_PROC_NULL,botrightdiag=MPI_PROC_NULL,upleftdiag=MPI_PROC_NULL,uprightdiag=MPI_PROC_NULL;
	
	if(south>=0){
		if(west>=0)
			botleftdiag=south-1;
		if(east>=0)
			botrightdiag=south+1;
	}
	if(north>=0){
		if(west>=0)
			upleftdiag=north-1;
		if(east>=0)
			uprightdiag=north+1;
	}
	MPI_Status readst;
	MPI_File fh;
	MPI_Offset filesize;
	MPI_File_open(MPI_COMM_WORLD,"waterfall_grey_1920_2520.raw",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_get_size(fh,&filesize);
	int bufsize=filesize/size;
	
	
	
	
	//two arrays for storing image
	int arrsizeI=bufsize/100;		
	int arrsizeJ=100;
	
	unsigned char *arr1=malloc(arrsizeI*arrsizeJ*sizeof(unsigned char));
	MPI_File_seek(fh,rank*bufsize,MPI_SEEK_SET);
	MPI_File_read(fh,arr1,arrsizeI*arrsizeJ,MPI_UNSIGNED_CHAR,&readst);		//read image (gray)
	//for(i=0;i<bufsize/10000;i++)
		//printf("%d,",arr1[i]);
	MPI_File_close(&fh);
	unsigned char *arr2=malloc(arrsizeI*arrsizeJ*sizeof(unsigned char));		//allocate 2nd array
	
	//making column type for communication.
	MPI_Datatype column;
	MPI_Type_vector(arrsizeI, 1, arrsizeJ, MPI_UNSIGNED_CHAR , &column);
	MPI_Type_commit(&column);
	//making a row type .
	MPI_Datatype row;
	MPI_Type_contiguous(arrsizeJ,MPI_UNSIGNED_CHAR,&row);
	MPI_Type_commit(&row);
	
	MPI_Status status[16];
	MPI_Request req[16];
	unsigned char *Nrow=malloc(arrsizeJ*sizeof(unsigned char));
	unsigned char *Srow=malloc(arrsizeJ*sizeof(unsigned char));
	unsigned char *Ecol=malloc(arrsizeI*sizeof(unsigned char));
	unsigned char *Wcol=malloc(arrsizeI*sizeof(unsigned char));
	
	memset(Nrow,'\0',arrsizeJ);
	memset(Srow,'\0',arrsizeJ);
	memset(Ecol,'\0',arrsizeI);
	memset(Wcol,'\0',arrsizeI);
	
	unsigned char upleft='\0',upright='\0',botleft='\0',botright='\0';
	
	
	int j,m=0;
	int sum=0;
	int count=0;
	double filter[9]={1,2,1,2,4,2,1,2,1};	//idio filtro me tin ekfwnisi
	sum=16;
	for(i=0;i<9;i++){
		filter[i]=filter[i]/(double)sum;
	}
	int flag;
	int count1=0;
	//helping struct for MPI_Test
	char checked[8]={'\0'};
	int lflag=1;
	int loopcount=0;
	unsigned char *tmp;
	//loop:
	int result;
		//prepare send and recv 
		MPI_Send_init(arr1,1,row,north,0,MPI_COMM_WORLD,&req[0]);		//sending upper-row to northern process.
		MPI_Send_init(arr1+(arrsizeI-1)*arrsizeJ,1,row,south,0,MPI_COMM_WORLD,&req[1]);	//sending bottom-row to southern process.
		MPI_Send_init(arr1+arrsizeJ-1, 1, column,east,0, MPI_COMM_WORLD,&req[2]); //sending right-column to east
		MPI_Send_init(arr1, 1, column, west, 0, MPI_COMM_WORLD,&req[3]);	//sending left column to west
		MPI_Send_init(&arr1[0],1,MPI_UNSIGNED_CHAR,upleftdiag,0,MPI_COMM_WORLD,&req[4]); 	//sending upper-left-diagonal element (one)
		MPI_Send_init(&arr1[arrsizeJ],1,MPI_UNSIGNED_CHAR,uprightdiag,0,MPI_COMM_WORLD,&req[5]);			//single element diagonal
		MPI_Send_init(&arr1[(arrsizeI-1)*arrsizeJ],1,MPI_UNSIGNED_CHAR,botleftdiag,0,MPI_COMM_WORLD,&req[6]);	// 	 ^^		^^		^^
		MPI_Send_init(&arr1[arrsizeI*arrsizeJ-1],1,MPI_UNSIGNED_CHAR,botrightdiag,0,MPI_COMM_WORLD,&req[7]);	// 	 ^^		^^		^^
		//receive
		MPI_Recv_init(Nrow,arrsizeJ, MPI_UNSIGNED_CHAR,north, 0, MPI_COMM_WORLD,&req[8]);
		MPI_Recv_init(Srow,arrsizeJ, MPI_UNSIGNED_CHAR,south, 0, MPI_COMM_WORLD,&req[9]);
		MPI_Recv_init(Ecol,arrsizeI, MPI_UNSIGNED_CHAR,east, 0, MPI_COMM_WORLD,&req[10]);
		MPI_Recv_init(Wcol,arrsizeI, MPI_UNSIGNED_CHAR,west, 0, MPI_COMM_WORLD,&req[11]);
		MPI_Recv_init(&upleft,1,MPI_UNSIGNED_CHAR,upleftdiag,0,MPI_COMM_WORLD,&req[12]);
		MPI_Recv_init(&upright,1,MPI_UNSIGNED_CHAR,uprightdiag,0,MPI_COMM_WORLD,&req[13]);
		MPI_Recv_init(&botleft,1,MPI_UNSIGNED_CHAR,botleftdiag,0,MPI_COMM_WORLD,&req[14]);
		MPI_Recv_init(&botright,1,MPI_UNSIGNED_CHAR,botrightdiag,0,MPI_COMM_WORLD,&req[15]);
	
	double startime,endtime;
	
	MPI_Barrier(MPI_COMM_WORLD);	
	startime=MPI_Wtime();
	while(1){
		//start all send-recv
		MPI_Start(&req[0]);
		MPI_Start(&req[1]);
		MPI_Start(&req[2]);
		MPI_Start(&req[3]);
		MPI_Start(&req[4]);
		MPI_Start(&req[5]);
		MPI_Start(&req[6]);
		MPI_Start(&req[7]);
		MPI_Start(&req[8]);
		MPI_Start(&req[9]);
		MPI_Start(&req[10]);
		MPI_Start(&req[11]);
		MPI_Start(&req[12]);
		MPI_Start(&req[13]);
		MPI_Start(&req[14]);
		MPI_Start(&req[15]);
		
		
		//apply filter at the inner array
		for(i=arrsizeJ+1;i<(arrsizeJ*(arrsizeI-1));i++){
			if(count==arrsizeJ-2){
				count++;
				continue;
			}
			else if(count==arrsizeJ-1){
				count=0;
				continue;
			}
			arr2[i]=(double)arr1[i-arrsizeJ-1]*filter[8]+(double)arr1[i-arrsizeJ]*filter[7]+(double)arr1[i-arrsizeJ+1]*filter[6]+(double)arr1[i-1]*filter[5]+(double)arr1[i]*filter[4]+(double)arr1[i+1]*filter[3]+(double)arr1[i+arrsizeJ-1]*filter[2]+(double)arr1[i+arrsizeJ]*filter[1]+(double)arr1[i+arrsizeJ+1]*filter[0];
			count++;
		}
		i=0;
		//test and if ready make some processesing with the received data
		while(i<8){
			MPI_Test(&req[8],&flag,&status[8]);
			if(flag && checked[0]=='\0'){
				i++;
				checked[0]=1;
				if(Nrow[0]!='\0'){
					for(j=1;j<arrsizeJ-1;j++){
						arr2[j]=(double)Nrow[j-1]*filter[8]+(double)Nrow[j]*filter[7]+(double)Nrow[j+1]*filter[6]+(double)arr1[j-1]*filter[5]+(double)arr1[j]*filter[4]+(double)arr1[j+1]*filter[3]+(double)arr1[j+arrsizeJ-1]*filter[2]+(double)arr1[j+arrsizeJ]*filter[1]+(double)arr1[j+arrsizeJ+1]*filter[0];
					}
				}
				else{
					for(j=1;j<arrsizeJ-1;j++)
						arr2[j]=arr1[j];
				}
			}
			MPI_Test(&req[9],&flag,&status[9]);
			if(flag && checked[1]=='\0'){
				i++;
				checked[1]=1;
				if(Srow[0]!='\0'){
					for(j=((arrsizeI-1)*arrsizeJ+1);j<(arrsizeJ*arrsizeI-1);j++){
						arr2[j]=(double)Srow[m]*filter[2]+(double)Srow[m+1]*filter[1]+(double)Srow[m+2]*filter[0]+(double)arr1[j-1]*filter[5]+(double)arr1[j]*filter[4]+(double)arr1[j+1]*filter[3]+(double)arr1[j-arrsizeJ+1]*filter[6]+(double)arr1[j-arrsizeJ]*filter[7]+(double)arr1[j-arrsizeJ-1]*filter[8];
						m++;
					}
				}
				else{
					for(j=((arrsizeI-1)*arrsizeJ+1);j<(arrsizeJ*arrsizeI-1);j++)
						arr2[j]=arr1[j];
				}
			}
			MPI_Test(&req[10],&flag,&status[10]);
			if(flag && checked[2]=='\0'){
				i++;
				checked[2]=1;
				if(Wcol[0]!='\0'){
					m=0;
					for(j=arrsizeJ;j<((arrsizeI-1)*arrsizeJ);j+=arrsizeJ){
						arr2[j]=(double)Wcol[m]*filter[8]+(double)Wcol[m+1]*filter[5]+(double)Wcol[m+2]*filter[2]+(double)arr1[j-arrsizeJ]*filter[7]+(double)arr1[j-arrsizeJ+1]*filter[6]+(double)arr1[j]*filter[4]+(double)arr1[j+1]*filter[3]+(double)arr1[j+arrsizeJ]*filter[1]+(double)arr1[j+arrsizeJ+1]*filter[0];
						m++;
					}
				}
				else{
					for(j=arrsizeJ;j<((arrsizeI-1)*arrsizeJ);j+=arrsizeJ)
						arr2[j]=arr1[j];
				}
			}
			
			MPI_Test(&req[11],&flag,&status[11]);
			if(flag && checked[3]=='\0'){
				i++;
				checked[3]=1;
				if(Ecol[0]!='\0'){
					m=0;
					for(j=2*arrsizeJ-1;j<=(arrsizeJ*(arrsizeI-1)-1);j+=arrsizeJ){
						arr2[j]=(double)Ecol[m]*filter[6]+(double)Ecol[m+1]*filter[3]+(double)Ecol[m+2]*filter[0]+(double)arr1[j-arrsizeJ-1]*filter[8]+(double)arr1[j-arrsizeJ]*filter[7]+(double)arr1[j-1]*filter[5]+(double)arr1[j]*filter[4]+(double)arr1[j+arrsizeJ-1]*filter[2]+(double)arr1[j+arrsizeJ]*filter[1];
						m++;
					}
				}
				else{
					for(j=2*arrsizeJ-1;j<=(arrsizeJ*(arrsizeI-1)-1);j+=arrsizeJ)
						arr2[j]=arr1[j];
				}
			}
			MPI_Test(&req[12],&flag,&status[12]);
			if(flag && checked[4]=='\0'){
				i++;
				checked[4]=1;
				if(upleft!='\0'){
					arr2[0]=(double)upleft*filter[8]+(double)Nrow[0]*filter[7]+(double)Nrow[1]*filter[6]+(double)Wcol[0]*filter[5]+(double)Wcol[1]*filter[2]+(double)arr1[0]*filter[4]+(double)arr1[1]*filter[3]+(double)arr1[arrsizeJ]*filter[1]+(double)arr1[arrsizeJ+1]*filter[0];
				}
				else{
					arr2[0]=arr1[0];
				}
			}
			
			MPI_Test(&req[13],&flag,&status[13]);
			if(flag && checked[5]=='\0'){
				i++;
				checked[5]=1;
				if(upright!='\0'){
					arr2[arrsizeJ-1]=(double)upright*filter[6]+(double)Nrow[arrsizeJ-2]*filter[8]+(double)Nrow[arrsizeJ-1]*filter[7]+(double)Ecol[0]*filter[3]+(double)Ecol[1]*filter[0]+(double)arr1[arrsizeJ-2]*filter[5]+(double)arr1[arrsizeJ-1]*filter[4]+(double)arr1[2*arrsizeJ-2]*filter[2]+(double)arr1[2*arrsizeJ-1]*filter[1];
				}
				else{
					arr2[arrsizeJ-1]=arr1[arrsizeJ-1];
				}
			}
			MPI_Test(&req[14],&flag,&status[14]);
			if(flag && checked[6]=='\0'){
				i++;
				checked[6]=1;
				if(botleft!='\0'){
					arr2[(arrsizeI-1)*arrsizeJ]=(double)Wcol[arrsizeI-2]*filter[8]+(double)Wcol[arrsizeI-1]*filter[5]+(double)Srow[0]*filter[1]+(double)Srow[1]*filter[0]+(double)botleft*filter[2]+(double)arr1[(arrsizeI-2)*arrsizeJ]*filter[7]+(double)arr1[(arrsizeI-2)*arrsizeJ+1]*filter[6]+(double)arr1[(arrsizeI-1)*arrsizeJ]*filter[4]+(double)arr1[(arrsizeI-2)*arrsizeJ+1]*filter[3];
				}
				else{
					arr2[(arrsizeI-1)*arrsizeJ]=arr1[(arrsizeI-1)*arrsizeJ];
				}
			}
			MPI_Test(&req[15],&flag,&status[15]);
			if(flag && checked[7]=='\0'){
				i++;
				checked[7]=1;
				if(botright!='\0'){
					arr2[arrsizeI*arrsizeJ-1]=(double)Srow[arrsizeJ-2]*filter[2]+(double)Srow[arrsizeJ-1]*filter[1]+(double)botright*filter[0]+(double)Ecol[arrsizeI-2]*filter[6]+(double)Ecol[arrsizeI-1]*filter[3]+(double)arr1[(arrsizeI-1)*arrsizeJ-2]*filter[8]+(double)arr1[(arrsizeI-1)*arrsizeJ-1]*filter[7]+(double)arr1[arrsizeI*arrsizeJ-2]*filter[5]+(double)arr1[arrsizeI*arrsizeJ-1]*filter[4];
				}
				else{
					arr2[arrsizeI*arrsizeJ-1]=arr1[arrsizeI*arrsizeJ-1];
				}
			}
		}
		i=0;
		while(i<8){
			MPI_Wait(&req[i],&status[i]);
			i++;
		}
		tmp=arr1;
		arr1=arr2;
		arr2=tmp;
		for(i=0;i<8;i++)
			checked[i]='\0';
		
		loopcount++;
		
	}
	MPI_Barrier(MPI_COMM_WORLD);	
	endtime=MPI_Wtime();
		if(rank==0)
			printf("LOOPCOUNT IS %d \n",loopcount);
		printf("time of proc %d is %lf\n",rank,endtime-startime);
	MPI_Finalize();
	return 0;
}

//finding north and south pole 
void compass(int *north,int *south,int rank,int numprocs)
{
	if((rank/numprocs)==0){
		*north=MPI_PROC_NULL;
		*south=rank+numprocs;
	}
	else if((rank/numprocs)==(numprocs-1)){
		*north=rank-numprocs;
		*south=MPI_PROC_NULL;
	}
	else{
		*north=rank-numprocs;
		*south=rank+numprocs;
	}		
}
