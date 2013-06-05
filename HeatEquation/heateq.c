#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"


#include "heateq.h"


#define a 0.1

REAL f(REAL x,REAL t){
	return 0.0;
}

REAL u_x0(REAL x){
	return 0.0;
}

REAL u_0t(REAL t){
	return 1.0-cos(t);
}

REAL u_1t(REAL t){
	return 0.0;
}

REAL u_next(REAL u_k_n_1,REAL u_k_n,REAL u_k_n1,REAL f_k_n,REAL dt,REAL h){
	return ((a*a*dt)/(h*h))*(u_k_n_1-2*u_k_n+u_k_n1)+dt*f_k_n+u_k_n;
}


void comp_n_n0(int rank,int N,int numprocs,int*n,int*n0){
	*n=(N-1)/numprocs;
	if(rank<(N-1)%numprocs){
		*n0 = *n*rank+rank+1;
	}else{
		*n0 = *n*rank+(N-1)%numprocs+1;
	}
	if(rank<(N-1)%numprocs)*n=*n+1;
}

int main(int argc, char *argv[]){

	REAL t,x;
	int i,j,k;

	int numprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);

	//printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

	if(argc<3){
		printf("Please enter two params (N, T)\n");
		MPI_Finalize();
		return;
	};

	int N; // ���������� ������-���� ���������� ���� ��
	int n; // ������-���� ���������� ���� x ������ �������������� ����������������
	int n0; // ���������� ������������ ��������
	REAL h; // ������ ����������

	int n_in; //  ������-���� ���������� ���� x �������������������� �������������� ���������������������� ���������������������� 0
	int n0_in;// ���������� ������������ �������� �������������������� �������������� ���������������������� ���������������������� 0

	REAL T; // ������������������ ����������
	REAL dt; // ������ ���� ��������������
	int nT; // ������-���� ���������� ���� ��������������

	REAL*u; // ������������ �� ���������������� ���� �������������� ��������
	REAL*w; // ������������ �� ���������������� ���� �������� ��������
	REAL*s; // ������������ ������ ���������� �������� ������������
	REAL bufoutL,bufoutR;
	REAL bufinL,bufinR;
	int RequestCount;
	MPI_Request requests[4];
	MPI_Status statuses[4];
	MPI_Status status;

	double Time0,Time1;

	Time0 = MPI_Wtime();


	N = atoi(argv[1]);
	T = atof(argv[2]);

	//printf("#%d:\tN = %d\tT = %d\n",rank,N,T);

	comp_n_n0(rank,N,numprocs,&n,&n0);

	h = 1.0/N;

	//printf("%d:\t %d %d %d\n",rank,n,n0,n+n0-1);

	// n ���������� ���������� + 2 �� ����������
	u = (REAL*)malloc((n+2)*sizeof(REAL));
	w = (REAL*)malloc((n+2)*sizeof(REAL));


	dt = h*h/(3*a*a); // ������ �������������� ���� �������������� ������������������������
	nT = (int)(T/dt)+1;
	dt = T/nT;

	// ������������������ ������������������ ������������ ������ ��=0
	for(i=0;i<=n+1;i++)u[i]=u_x0((n0-1+i)*h);

	for(k=0;k<nT;k++){
		//if(rank==0)printf("%d:\t time=%f (%d/%d)\n",rank,k*dt,k,nT);

		// ���������������������� �������������������� ����������
		for(i=1;i<=n;i++)w[i]= u_next(u[i-1],u[i],u[i+1],f((n0+i)*h,k*dt),dt,h);

		RequestCount = 0;

		// �������� �������� ���� ����������, ���������������� ������������ ������������ �������� �������������� ����������
		if(rank!=0){
			bufoutL=w[1];
			MPI_Isend(&bufoutL,1,MPIREAL,rank-1,0,MPI_COMM_WORLD,&requests[RequestCount]);
			RequestCount++;
		}
		// �������� �������� ���� ������������, ���������������� �������������� ������������ �������� ������������ �������������� ����������
		if(rank!=numprocs-1){
			bufoutR=w[n];
			MPI_Isend(&bufoutR,1,MPIREAL,rank+1,0,MPI_COMM_WORLD,&requests[RequestCount]);
			RequestCount++;
		}
		// �������� �������� ���� ����������, ������������������ ���������� ���� ������������ ������������
		if(rank!=0){
			MPI_Irecv(&bufinL,1,MPIREAL,rank-1,0,MPI_COMM_WORLD,&requests[RequestCount]);
			RequestCount++;
		}
		// �������� �������� ���� ������������, ������������������ ���������� ���� �������������� ������������
		if(rank!=numprocs-1){
			MPI_Irecv(&bufinR,1,MPIREAL,rank+1,0,MPI_COMM_WORLD,&requests[RequestCount]);
			RequestCount++;
		}

		// �������������������� �������������������� ���������������� �� ����������������
		MPI_Waitall(RequestCount,requests,statuses);

		if(rank!=0){
			w[0]=bufinL;
		}else{
			w[0]=u_0t(dt*k);
		}
		if(rank!=numprocs-1){
			w[n+1]=bufinR;
		}else{
			w[n+1]=u_1t(dt*k);
		}

		// ������������������ ������������ �� ������������ �������� ���� ������������;
		memcpy(u,w,(n+2)*sizeof(REAL));
		for(i=0;i<=n+1;i++)u[i]=w[i];

		//printf("%d:\t sendrecv complited\n",rank);

	}

	// ���������� �������������������� ������ ������������ ���� ������������������ 0

	free(w);
	if(rank!=0){
		MPI_Send(u+1,n,MPIREAL,0,0,MPI_COMM_WORLD);
		//printf("send from %d to %d (p=%d c=%d)\n",n0,n0+n-1,rank,n);
	}else{
		s=(REAL*)malloc((N+1)*sizeof(REAL));
		for(i=0;i<=N;i++)s[i]=-1.0;

		s[0]=u_0t(T);
		s[N]=u_1t(T);


		//printf("copy from %d to %d (p=%d c=%d)\n",n0,n0+n-1,i,n);
		for(i=1;i<=n;i++)s[i]=u[i];
		j=n;

		for(i=1;i<numprocs;i++){
			comp_n_n0(i,N,numprocs,&n,&n0);
			MPI_Recv(s+n0,n,MPIREAL,i,0,MPI_COMM_WORLD,&status);
			//printf("recv from %d to %d (p=%d c=%d)\n",n0,n0+n-1,i,n);
		}

		//printf("m: [\n");
		//for(i=0;i<=N;i++)printf("\t%f,\n",s[i]);
		//printf("]\n");

		Time1 = MPI_Wtime();

		FILE*fd;
		fd=fopen("report.txt","w");
		for(i=0;i<=N;i++)fprintf(fd,"%lf\n",s[i]);
		fclose(fd);
		if(rank==0)printf("Time of computing %f\n",Time1-Time0);

	}



	//for(i=1;i<n-1;i++)


	MPI_Finalize();


	return 0;
}






