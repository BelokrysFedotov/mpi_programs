/**
 * Вычисление интеграла sin(1/x) на отрезке [A;B] с точностью delta
 * Параметры A,B,delta передаются через аргументы программы
 *
 * Интеграл вычисляется на неравномерной сетке, которая строится в два этапа.
 * На первом этапе выбираются нули функции: sin(1/x)=0; получается N отрезков
 * На втором этапе, каждый отрезок разбивается на равномерную подсетку так,
 * что бы погршность полученая на этом отрезке была меньше delta/N
 * Полученая сетка разделяется равномерно между процессорами геометрическим методом.
 * Интеграл вычисляется методом трапеций.
 * После вычисления значений интегралов на своем учаске, процессоры передают значения
 * на нулевой процессор, который их складывает и получает полное значение интеграла.
 */
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"


#include "integral.h"

// Интегрируемая функция
REAL f(REAL x){
	return sin(1/x);
}

// Квадрат числа
REAL pow2(REAL X){return X*X;}
// Куб числа
REAL pow3(REAL X){return pow2(X)*X;}
// Четвертая степень числа
REAL pow4(REAL X){REAL Y;Y=pow2(X);return pow2(Y);}

/**
 * вычисление кол-ва узлов сетки, которое необходимо для интегрирования
 * функции на отрезке [A;B] с точностью D
 */
INT comp_N(REAL A,REAL B,REAL D){
	return (INT)ceil(sqrt(((B+B+1)*pow3(B-A))/(12*pow4(A)*D)))+1;
}

// Метод трапеций для отрезка [A;B]
REAL comp_S(REAL A,REAL B){
	return 	0.5*(f(A)+f(B))*(B-A);
}

int main(int argc, char *argv[]){

	INT i,j,k;
	REAL A,B,delta;
	REAL x1,x2,s,S;
	INT na,nb,n,N,N1,N2,n1,n2,m;

	int numprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double Time0,Time1;

	// Инициализация MPI и получение начальных параметров
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);

	// сообщение пользователю о процессе
	printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

	// проверяем кол-во переданных параметров, если их меньше 3 то заканчиваем работу
	if(argc<3){
		if(rank==0)printf("Please enter 3 params (A, B, delta)\n");
		MPI_Finalize();
		return;
	};

	// время начала расчетов
	Time0 = MPI_Wtime();

	// получаем начальные параметры
	A = atof(argv[1]);
	B = atof(argv[2]);
	delta = atof(argv[3]);

	// выводим полученые параметры
	if(rank==0)printf("Integral_%.2f^%.2f sin(1/x) dx. Delta = %.2e\n",A,B,delta);

	// находя нули функции, мы получаем, что x = 1/(pi*n); где n - целое
	// находим крайние значения n лежащие в отрезке [A;B]
	na = ((INT)ceil(1_pi/A))-1;
	nb = ((INT)floor(1_pi/B))+1;

	// выводим эти значения
	if(rank==0)printf("n = %ld..%ld\n",nb,na);
	// узлы первой сетки лежат в точках 1/(pi*n) где n = nb..na

	// N - общее кол-во узлов итоговой сетки
	N=0;
	// проходим по всем отрезкам первой сетки
	for(n=nb-1;n<=na;n++){
		// с учетом крайних отрезков
		if(n==nb-1){x2=B;}else{x2=1_pi/n;}
		if(n==na){x1=A;}else{x1=1_pi/(n+1);}
		// суммируем кол-во подузлов в этом отрезке
		N+=comp_N(x1,x2,delta/(na-nb+2));
	}
	if(rank==0)printf("N = %ld\n",N);

	// k - кол-во узлов итоговой сетки, которое получит каждый процессор, кроме последнего, он получает все что осталось
	k = N/numprocs;
	if(N%numprocs)k++;
	// вычисляем границы куска сетки для данного процессора
	N1 = k*rank;
	N2 = N1+k;if(rank+1==numprocs)N2=N;
	//printf("%d:\t[%ld,%ld]\t%ld\n",rank,N1,N2,N2-N1);

	//ищем интеграл на текущем процессоре
	s=0;
	// проходим все отрезки первой сетки
	for(n=nb-1,m=0;n<=na;n++,m=j){
		if(n==nb-1){x2=B;}else{x2=1_pi/n;}
		if(n==na){x1=A;}else{x1=1_pi/(n+1);}

		// вычисляем кол-во узлов в текущем отрезке
		k=comp_N(x1,x2,delta/(na-nb+2));
		// m - порядковый номер левого узла отрезка в общей сетке
		// j - правого узла
		j=m+k;
		//printf("see [%ld;%ld]\n",m,j);

		/**
		 * Проверяем стоит ли интегрировать данный отрезок или его часть
		 * n1 - локальный номер левого узла в этом отрезке
		 * n2 - правого узла
		 * если n2>n1 то этот отрезок не нужно интегрировать
		 */
		n1 = m;
		n2 = j;
		if(n1<N1)n1 = N1;
		n1 = n1 - m;
		if(n2>N2)n2 = N2;
		n2 = n2 - m;
		//printf("work[%ld;%ld]\n",n1,n2);

		//интегрирование отрезка
		for(i=n1;i<n2;i++)s+=comp_S(x1*(k-i)/k+x2*i/k,x1*(k-i-1)/k+x2*(i+1)/k);
	}
	//printf("S = %f\n",s);

	// отсылка суммы на первый процессор с паралельным сумированием
	MPI_Reduce(&s,&S,1,MPIREAL,MPI_SUM,0,MPI_COMM_WORLD);

	// вывод значения интеграла
	if(rank==0)printf("Integral = %f\n",S);

	// вычисление времени работы
	Time1 = MPI_Wtime();
	if(rank==0)printf("Time of computing %f\n",Time1-Time0);

	MPI_Finalize();
	return 0;
}






