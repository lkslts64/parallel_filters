#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"

;int main (int argc, char* argv[])
{
  int rank, size;
 
  MPI_Init (&argc, &argv);      /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */

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
	
	//two arrays for storing image
	long bufsize=4838400/size;
	int arrsizeI=bufsize/100;
	int arrsizeJ=100;
	struct RGB *arr1=malloc(arrsizeI*arrsizeJ*sizeof(struct RGB));		
	struct RGB *arr2=malloc(arrsizeI*arrsizeJ*sizeof(struct RGB));
	
	for(i=0;i<arrsizeI*arrsizeJ;i++){
		arr1[i].rgb[0]=rand()%256;
		arr1[i].rgb[1]=rand()%256;
		arr1[i].rgb[2]=rand()%256;
	}

	MPI_Datatype rgb;
	MPI_Type_contiguous(3,MPI_UNSIGNED_CHAR,&rgb);
	MPI_Type_commit(&rgb);
	
	
	//making a row type .
	MPI_Datatype row;
	MPI_Type_contiguous(arrsizeJ,rgb,&row);
	MPI_Type_commit(&row);
	//making column type for communication.
	MPI_Datatype column;
	MPI_Type_vector(arrsizeI, 1, arrsizeJ,rgb, &column);
	MPI_Type_commit(&column);	
	
	
	
	
	MPI_Status status[16];
	MPI_Request req[16];
	struct RGB *Nrow=malloc(arrsizeJ*sizeof(struct RGB));
	struct RGB *Srow=malloc(arrsizeJ*sizeof(struct RGB));
	struct RGB *Ecol=malloc(arrsizeI*sizeof(struct RGB));
	struct RGB *Wcol=malloc(arrsizeI*sizeof(struct RGB));
	
	Nrow[0].rgb[0]='\0';
	Srow[0].rgb[0]='\0';
	Ecol[0].rgb[0]='\0';
	Wcol[0].rgb[0]='\0';
	
	struct RGB upleft,upright,botleft,botright;
	upleft.rgb[0]='\0';
	upright.rgb[0]='\0';
	botleft.rgb[0]='\0';
	botright.rgb[0]='\0';
	
	int j,m=0;
	int sum=0;
	int count=0;
	
	double filter[9]={0,0,2,0,1,0,0,0,1};
	sum=4;
	
	
	for(i=0;i<9;i++){
		filter[i]=filter[i]/(double)sum;
	}
	
	int flag;
	int count1=0;
	//helping struct for MPI_Test
	char checked[8]={'\0'};
	int lflag=1;
	int loopcount=0;
	struct RGB *tmp;
	int result;
	int k;
		
		
		
		MPI_Send_init(arr1,1,row,north,0,MPI_COMM_WORLD,&req[0]);		//sending upper-row to northern process.
		MPI_Send_init(arr1+90,1,row,south,0,MPI_COMM_WORLD,&req[1]);	//sending bottom-row to southern process.

		MPI_Send_init(arr1+9, 1, column,east,0, MPI_COMM_WORLD,&req[2]); //sending right-column to east
		MPI_Send_init(arr1, 1, column, west, 0, MPI_COMM_WORLD,&req[3]);	//sending left column to west
		MPI_Send_init(&arr1[0],3,MPI_UNSIGNED_CHAR,upleftdiag,0,MPI_COMM_WORLD,&req[4]); 	//sending upper-left-diagonal element (one)
		MPI_Send_init(&arr1[arrsizeJ],3,MPI_UNSIGNED_CHAR,uprightdiag,0,MPI_COMM_WORLD,&req[5]);			//single element diagonal
		MPI_Send_init(&arr1[arrsizeI*arrsizeJ-10],3,MPI_UNSIGNED_CHAR,botleftdiag,0,MPI_COMM_WORLD,&req[6]);	// 	 ^^		^^		^^
		MPI_Send_init(&arr1[arrsizeI*arrsizeJ-1],3,MPI_UNSIGNED_CHAR,botrightdiag,0,MPI_COMM_WORLD,&req[7]);	// 	 ^^		^^		^^
		//receive
		MPI_Recv_init(Nrow,3*arrsizeJ, MPI_UNSIGNED_CHAR,north, 0, MPI_COMM_WORLD,&req[8]);
		//MPI_Irecv(Nrow,1, row,north, 0, MPI_COMM_WORLD,&req[8]);
		MPI_Recv_init(Srow,3*arrsizeJ, MPI_UNSIGNED_CHAR,south, 0, MPI_COMM_WORLD,&req[9]);
		MPI_Recv_init(Ecol,3*arrsizeI, MPI_UNSIGNED_CHAR,east, 0, MPI_COMM_WORLD,&req[10]);
		MPI_Recv_init(Wcol,3*arrsizeI, MPI_UNSIGNED_CHAR,west, 0, MPI_COMM_WORLD,&req[11]);
		MPI_Recv_init(&upleft,3,MPI_UNSIGNED_CHAR,upleftdiag,0,MPI_COMM_WORLD,&req[12]);
		MPI_Recv_init(&upright,3,MPI_UNSIGNED_CHAR,uprightdiag,0,MPI_COMM_WORLD,&req[13]);
		MPI_Recv_init(&botleft,3,MPI_UNSIGNED_CHAR,botleftdiag,0,MPI_COMM_WORLD,&req[14]);
		MPI_Recv_init(&botright,3,MPI_UNSIGNED_CHAR,botrightdiag,0,MPI_COMM_WORLD,&req[15]);
		
	double startime,endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	startime=MPI_Wtime();
	while(1){
		
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
		
		
		
		for(i=arrsizeJ+1;i<(arrsizeJ*(arrsizeI-1));i++){
			if(count==arrsizeJ-2){
				count++;
				continue;
			}
			else if(count==arrsizeJ-1){
				count=0;
				continue;
			}
			for(k=0;k<3;k++)
				arr2[i].rgb[k]=(double)arr1[i-arrsizeJ-1].rgb[k]*filter[8]+(double)arr1[i-arrsizeJ].rgb[k]*filter[7]+(double)arr1[i-arrsizeJ+1].rgb[k]*filter[6]+(double)arr1[i-1].rgb[k]*filter[5]+(double)arr1[i].rgb[k]*filter[4]+(double)arr1[i+1].rgb[k]*filter[3]+(double)arr1[i+arrsizeJ-1].rgb[k]*filter[2]+(double)arr1[i+arrsizeJ].rgb[k]*filter[1]+(double)arr1[i+arrsizeJ+1].rgb[k]*filter[0];
			count++;
		}
		i=0;
		while(i<8){
			MPI_Test(&req[8],&flag,&status[8]);
			if(flag && checked[0]=='\0'){
				i++;
				checked[0]=1;
				if(Nrow[0].rgb[0]!='\0'){
					for(j=1;j<arrsizeJ-1;j++){
						for(k=0;k<3;k++)
							arr2[j].rgb[k]=(double)Nrow[j-1].rgb[k]*filter[8]+(double)Nrow[j].rgb[k]*filter[7]+(double)Nrow[j+1].rgb[k]*filter[6]+(double)arr1[j-1].rgb[k]*filter[5]+(double)arr1[j].rgb[k]*filter[4]+(double)arr1[j+1].rgb[k]*filter[3]+(double)arr1[j+arrsizeJ-1].rgb[k]*filter[2]+(double)arr1[j+arrsizeJ].rgb[k]*filter[1]+(double)arr1[j+arrsizeJ+1].rgb[k]*filter[0];
					}
				}
				else{
					for(j=1;j<arrsizeJ-1;j++)
						for(k=0;k<3;k++)
							arr2[j].rgb[k]=arr1[j].rgb[k];
				}
			}
			MPI_Test(&req[9],&flag,&status[9]);
			if(flag && checked[1]=='\0'){
				i++;
				checked[1]=1;
				if(Srow[0].rgb[0]!='\0'){
					for(j=((arrsizeI-1)*arrsizeJ+1);j<(arrsizeJ*arrsizeI-1);j++){
						for(k=0;k<3;k++){
							arr2[j].rgb[k]=(double)Srow[m].rgb[k]*filter[2]+(double)Srow[m+1].rgb[k]*filter[1]+(double)Srow[m+2].rgb[k]*filter[0]+(double)arr1[j-1].rgb[k]*filter[5]+(double)arr1[j].rgb[k]*filter[4]+(double)arr1[j+1].rgb[k]*filter[3]+(double)arr1[j-arrsizeJ+1].rgb[k]*filter[6]+(double)arr1[j-arrsizeJ].rgb[k]*filter[7]+(double)arr1[j-arrsizeJ-1].rgb[k]*filter[8];
							m++;
						}
					}
				}
				else{
					for(j=((arrsizeI-1)*arrsizeJ+1);j<(arrsizeJ*arrsizeI-1);j++)
						for(k=0;k<3;k++)
							arr2[j].rgb[k]=arr1[j].rgb[k];
				}
			}
			MPI_Test(&req[10],&flag,&status[10]);
			if(flag && checked[2]=='\0'){
				i++;
				checked[2]=1;
				if(Wcol[0].rgb[0]!='\0'){
					m=0;
					for(j=arrsizeJ;j<((arrsizeI-1)*arrsizeJ);j+=arrsizeJ){
						for(k=0;k<3;k++){
							arr2[j].rgb[k]=(double)Wcol[m].rgb[k]*filter[8]+(double)Wcol[m+1].rgb[k]*filter[5]+(double)Wcol[m+2].rgb[k]*filter[2]+(double)arr1[j-arrsizeJ].rgb[k]*filter[7]+(double)arr1[j-arrsizeJ+1].rgb[k]*filter[6]+(double)arr1[j].rgb[k]*filter[4]+(double)arr1[j+1].rgb[k]*filter[3]+(double)arr1[j+arrsizeJ].rgb[k]*filter[1]+(double)arr1[j+arrsizeJ+1].rgb[k]*filter[0];
							m++;
						}
					}
				}
				else{
					for(j=arrsizeJ;j<((arrsizeI-1)*arrsizeJ);j+=arrsizeJ)
						for(k=0;k<3;k++)
							arr2[j].rgb[k]=arr1[j].rgb[k];
				}
			}
			
			MPI_Test(&req[11],&flag,&status[11]);
			if(flag && checked[3]=='\0'){
				i++;
				checked[3]=1;
				if(Ecol[0].rgb[0]!='\0'){
					m=0;
					for(j=2*arrsizeJ-1;j<=(arrsizeJ*(arrsizeI-1)-1);j+=arrsizeJ){
						for(k=0;k<3;k++){
							arr2[j].rgb[k]=(double)Ecol[m].rgb[k]*filter[6]+(double)Ecol[m+1].rgb[k]*filter[3]+(double)Ecol[m+2].rgb[k]*filter[0]+(double)arr1[j-arrsizeJ-1].rgb[k]*filter[8]+(double)arr1[j-arrsizeJ].rgb[k]*filter[7]+(double)arr1[j-1].rgb[k]*filter[5]+(double)arr1[j].rgb[k]*filter[4]+(double)arr1[j+arrsizeJ-1].rgb[k]*filter[2]+(double)arr1[j+arrsizeJ].rgb[k]*filter[1];
							m++;
						}
					}
				}
				else{
					for(j=2*arrsizeJ-1;j<=(arrsizeJ*(arrsizeI-1)-1);j+=arrsizeJ)
						for(k=0;k<3;k++)
							arr2[j].rgb[k]=arr1[j].rgb[k];
				}
			}
			MPI_Test(&req[12],&flag,&status[12]);
			if(flag && checked[4]=='\0'){
				i++;
				checked[4]=1;
				if(upleft.rgb[0]!='\0'){
					for(k=0;k<3;k++)
						arr2[0].rgb[k]=(double)upleft.rgb[k]*filter[8]+(double)Nrow[0].rgb[k]*filter[7]+(double)Nrow[1].rgb[k]*filter[6]+(double)Wcol[0].rgb[k]*filter[5]+(double)Wcol[1].rgb[k]*filter[2]+(double)arr1[0].rgb[k]*filter[4]+(double)arr1[1].rgb[k]*filter[3]+(double)arr1[arrsizeJ].rgb[k]*filter[1]+(double)arr1[arrsizeJ+1].rgb[k]*filter[0];
				}
				else{
					for(k=0;k<3;k++)
						arr2[0].rgb[k]=arr1[0].rgb[k];
				}
			}
			
			MPI_Test(&req[13],&flag,&status[13]);
			if(flag && checked[5]=='\0'){
				i++;
				checked[5]=1;
				if(upright.rgb[0]!='\0'){
					for(k=0;k<3;k++)
						arr2[arrsizeJ-1].rgb[k]=(double)upright.rgb[k]*filter[6]+(double)Nrow[arrsizeJ-2].rgb[k]*filter[8]+(double)Nrow[arrsizeJ-1].rgb[k]*filter[7]+(double)Ecol[0].rgb[k]*filter[3]+(double)Ecol[1].rgb[k]*filter[0]+(double)arr1[arrsizeJ-2].rgb[k]*filter[5]+(double)arr1[arrsizeJ-1].rgb[k]*filter[4]+(double)arr1[2*arrsizeJ-2].rgb[k]*filter[2]+(double)arr1[2*arrsizeJ-1].rgb[k]*filter[1];
				}
				else{
					for(k=0;k<3;k++)
						arr2[arrsizeJ-1].rgb[k]=arr1[arrsizeJ-1].rgb[k];
				}
			}
			MPI_Test(&req[14],&flag,&status[14]);
			if(flag && checked[6]=='\0'){
				i++;
				checked[6]=1;
				if(botleft.rgb[0]!='\0'){
					for(k=0;k<3;k++)
						arr2[(arrsizeI-1)*arrsizeJ].rgb[k]=(double)Wcol[arrsizeI-2].rgb[k]*filter[8]+(double)Wcol[arrsizeI-1].rgb[k]*filter[5]+(double)Srow[0].rgb[k]*filter[1]+(double)Srow[1].rgb[k]*filter[0]+(double)botleft.rgb[k]*filter[2]+(double)arr1[(arrsizeI-2)*arrsizeJ].rgb[k]*filter[7]+(double)arr1[(arrsizeI-2)*arrsizeJ+1].rgb[k]*filter[6]+(double)arr1[(arrsizeI-1)*arrsizeJ].rgb[k]*filter[4]+(double)arr1[(arrsizeI-2)*arrsizeJ+1].rgb[k]*filter[3];
				}
				else{
					for(k=0;k<3;k++)
						arr2[(arrsizeI-1)*arrsizeJ].rgb[k]=arr1[(arrsizeI-1)*arrsizeJ].rgb[k];
				}
			}
			MPI_Test(&req[15],&flag,&status[15]);
			if(flag && checked[7]=='\0'){
				i++;
				checked[7]=1;
				if(botright.rgb[0]!='\0'){
					for(k=0;k<3;k++)
						arr2[arrsizeI*arrsizeJ-1].rgb[k]=(double)Srow[arrsizeJ-2].rgb[k]*filter[2]+(double)Srow[arrsizeJ-1].rgb[k]*filter[1]+(double)botright.rgb[k]*filter[0]+(double)Ecol[arrsizeI-2].rgb[k]*filter[6]+(double)Ecol[arrsizeI-1].rgb[k]*filter[3]+(double)arr1[(arrsizeI-1)*arrsizeJ-2].rgb[k]*filter[8]+(double)arr1[(arrsizeI-1)*arrsizeJ-1].rgb[k]*filter[7]+(double)arr1[arrsizeI*arrsizeJ-2].rgb[k]*filter[5]+(double)arr1[arrsizeI*arrsizeJ-1].rgb[k]*filter[4];
				}
				else{
					for(k=0;k<3;k++)
						arr2[arrsizeI*arrsizeJ-1].rgb[k]=arr1[arrsizeI*arrsizeJ-1].rgb[k];
				}
			}
		}
		i=0;
		while(i<8){
			MPI_Wait(&req[i],&status[i]);
			i++;
		}
		i=0;
		for(i=0;i<arrsizeJ*arrsizeI;i++)
			for(k=0;k<3;k++)
				if(arr1[i].rgb[k]==arr2[i].rgb[k])
					count1++;
		if(count1==3*arrsizeJ*arrsizeI)
			lflag=0;
		
		
		MPI_Allreduce(&lflag,&result,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		if(!result)
			break;
		count1=0;
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
			printf("LOOPCOUNT IS %d\n",loopcount);
		printf("TIME IS %lf\n",endtime-startime);
		
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
