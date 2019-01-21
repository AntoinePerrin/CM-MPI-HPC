#include <mpi.h>
#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include "matrix.h"

using namespace std;


int npes, myrank;

const double dx = 0.5;
const double u = 250;
const double xTot = 400;
const double PI = 3.14;


/*!
 *  \brief Send a vector to another node
 *
 *  Send a vector to another node with the function MPI_Send
 *
 *  \param vector the vector to send
 *  \param dest the destination of the vector
 *  \param comm the value of MPI_COMM_WORLD
 */

void send(vector<double> const& vector, int dest, MPI_Comm comm){								//declaration of the function send
	unsigned len = vector.size();																//declare a variable len of the vector
	MPI_Send(&len, 1, MPI_UNSIGNED, dest, 1, comm);												//send the lenght of the vector to the desire node
	if(len !=0){																				//if the lenght is not null =>
		MPI_Send(vector.data(), len, MPI_DOUBLE, dest, 1, comm);								//send the vector with MPI_Send
	}
}

/*!
 *  \brief Receive a vector from another node
 *
 *  Receive a vector from another node with the function MPI_Recv
 *
 *  \param vector the vector to receive
 *  \param dest the source of the vector
 *  \param comm the value of MPI_COMM_WORLD
 */

vector<double> recv(vector<double>& res, int src, MPI_Comm comm){								//declaration of the function send
	MPI_Status status;																			//get the status of MPI
	unsigned len;																				//create a variable len
	MPI_Recv(&len, 1, MPI_UNSIGNED, src, 1, comm, &status);										//get the lenght of the vector
	if(len !=0){																				//if the lenght is superior than 0 =>
		res.resize(len);																		//resize the lenght of the vector to len
		MPI_Recv(res.data(), len, MPI_DOUBLE, src, 1, comm, &status);							//receive the vector
	}
	return res;																					//return the vector res
}

/*!
 *  \brief Print the vector given
 *
 *  Print every composant of the vector
 *
 *  \param vect A vector to print
 */

void showVector(vector <double> vect) {															//decalaration of the function which can show all the value of a vector of double
	for (unsigned int i = 0; i < vect.size(); i++) {											//create a loop with i which will be used as an index
		printf("%0.1f \n", vect[i]);															//print the value of the vector at the index i
	}
	printf("\n");																				//print a enter for more clarity
}

/*!
 *  \brief function f(x)
 *
 *  Return the value of f(x)
 *
 *  \param x Value of x
 *  \return Return a double, the value of f(x)
 */
double finitx(double x) {																		//declaration of the function f(x)
	if (x <= 50 || x >= 110)																	//if x is lower than 50 or higher than 110 =>
		return 0.0;																				//return 0
	return 100 * (sin(PI*(x - 50) / 60));														//if x is between 50 and 110 return the wave
}

/*!
 *  \brief Create the vector f(x)
 *
 *  Get the value of f(x) for each x and put it in the return vector
 *
 *  \return Return a vector of every f(x)
 */

vector<double> finit() {																		//declaration of the function which create the vector with all the value of f(x)
	vector <double> res;																		//decalartion of the vector of double
	double x = 0;																				//initiate the variable x at 0
	while (x <= xTot) {																			//while x is lower than xtot (400) =>
		res.push_back(finitx(x));																//get the value of f(x) and add it to the vector res
		x += dx;																				//add dx to x
	}
	return res;																					//return the vector with all the value of f(x)
}

/*!
 *  \brief Handle the Thomas algorithm
 *
 *  Function called to calculate the Thomas algorithm in serial
 *
 *  \param A The matrix A
 *  \param f The vector f
 *  \return Return a vector of the solution at time n + 1
 */

vector <double> ThomasAlgorithmSerial(Matrix A, vector <double> f) {							//declare the function which handle the thomas algorithm in serial
	int sizeOfA = A.getNcols();																	//declare an int to the number of columunm of the matrix A
	vector <double> x(sizeOfA);																	//declare a vector x
	vector <double> m(sizeOfA);																	//declare a vector m
	vector <double> bPrime(sizeOfA);															//declare a vector bPrime
	vector <double> dPrime(sizeOfA);															//declare a vector dPrime

	// We create vectors from the non null coefficients of the matrix
	vector <double> a(sizeOfA);																	//declare a vector a
	vector <double> b(sizeOfA);																	//declare a vector b
	vector <double> c(sizeOfA);																	//declare a vector c
																								//all of this vector have the same size which is sizeofA
	for (int i = 0; i < sizeOfA; i++) {															//create a loop with an index i until i reach the sizeofA
		b[i] = A[i][i];																			//set the vector b at the index i the value of the matrix A at [i][i]
		if (i > 0)																				//if i higher than 0 =>
			a[i] = A[i][i - 1];																	//set the vector a at the index i the value of the matrix A at [i][i-1]
		if (i < sizeOfA-1)																		//if i is lower than sizeOfA - 1 =>
			c[i] = A[i][i + 1];																	//set the vector c at the index i the value of the matrix A at [i][i+1]
	}
	a[0] = 0;																					//set the first value of a to 0
	c[sizeOfA - 1] = 0;																			//set last value of c to 0


	// First values of the vectors
	m[0] = 0;																					//set the first value of m to 0
	m[1] = a[1] / b[0];																			//set the second value of m to the second value of a divided by the first value of b
	bPrime[0] = b[0];																			//set the first value of bPrime to the fisrt value of b
	dPrime[0] = f[0];																			//set the first value of dPrime to the fisrt value of f

	// We calculate all values of the vectors
	for (int i = 1; i < sizeOfA; i++) {															//create a loop with an index i until i reach sizeOfA
		bPrime[i] = b[i] - m[i] * c[i - 1];														//set bPrime at index i
		dPrime[i] = f[i] - m[i] * dPrime[i - 1];												//set dPrime at index i
		if (i < sizeOfA - 1)																	//if i is lower than sizeOfA - 1 =>
			m[i + 1] = a[i + 1] / bPrime[i];													//set m at index i+1
	}

	x[sizeOfA - 1] = dPrime[sizeOfA - 1] / bPrime[sizeOfA - 1];									//set x at index sizeOfA - 1
	for (int i = sizeOfA - 2; i >= 0; i--)														//create a loop in reverse, the index i start at sizeOfA and decrease until reach 0
		x[i] = (dPrime[i] - c[i] * x[i + 1]) / bPrime[i];										//set x at index i

	return x;																					//return vector x
}

/*!
 *  \brief Handle the Thomas algorithm
 *
 *  Function called to calculate the Thomas algorithm in parallel
 *
 *  \param A The matrix A
 *  \param f The vector f
 *  \return Return a vector of the solution at time n + 1
 */

vector <double> ThomasAlgorithmPara(Matrix A, vector <double> f) {								//declare the function which handle the thomas algorithm in parallel
	int sizeOfA = A.getNcols();																	//declare an int x to the number of columunm of the matrix A
	vector <double> x(sizeOfA);																	//declare a vector x
	vector <double> m(sizeOfA);																	//declare a vector m
	vector <double> bPrime(sizeOfA);															//declare a vector bPrime
	vector <double> dPrime(sizeOfA);															//declare a vector dPrime

	// We create vectors from the non null coefficients of the matrix
	vector <double> a(sizeOfA);																	//declare a vector a
	vector <double> b(sizeOfA);																	//declare a vector b
	vector <double> c(sizeOfA);																	//declare a vector c
																								//all of this vector have the same size which is sizeofA
	if (myrank == 0){																			//if the function run on the first node
		for (int i = 0; i < sizeOfA; i++) {														//create a loop with an index i until i reach the sizeofA
			b[i] = A[i][i];																		//set the vector b at the index i the value of the matrix A at [i][i]
			if (i > 0){																			//if i higher than 0 =>
				a[i] = A[i][i - 1];																//set the vector a at the index i the value of the matrix A at [i][i-1]
			}
			if (i < sizeOfA-1){																	//if i is lower than sizeOfA - 1 =>
				c[i] = A[i][i + 1];																//set the vector c at the index i the value of the matrix A at [i][i+1]
			}
		}
		a[0] = 0;																				//set the first value of a to 0
		c[sizeOfA - 1] = 0;																		//set last value of c to 0
		send(a,1,MPI_COMM_WORLD);																//send the vector a to the node 1
		send(b,1,MPI_COMM_WORLD);																//send the vector b to the node 1
		send(c,1,MPI_COMM_WORLD);																//send the vector c to the node 1
		x = recv(x, npes - 1, MPI_COMM_WORLD);													//stuck the node 0 until the other node complete
	} else if (myrank == 1){																	//if the function run on the second node
		a = recv(a,0,MPI_COMM_WORLD);															//receive the vector a from the node 0
		b = recv(b,0,MPI_COMM_WORLD);															//receive the vector b from the node 0
		c = recv(c,0,MPI_COMM_WORLD);															//receive the vector c from the node 0
		// First values of the vectors
		m[0] = 0;																				//set the first value of m to 0
		m[1] = a[1] / b[0];																		//set the second value of m to the second value of a divided by the first value of b
		bPrime[0] = b[0];																		//set the first value of bPrime to the fisrt value of b
		dPrime[0] = f[0];																		//set the first value of dPrime to the fisrt value of f
	
		// We calculate all values of the vectors
		for (int i = 1; i < sizeOfA; i++) {														//create a loop with an index i until i reach sizeOfA
			bPrime[i] = b[i] - m[i] * c[i - 1];													//set bPrime at index i
			dPrime[i] = f[i] - m[i] * dPrime[i - 1];											//set dPrime at index i
			if (i < sizeOfA - 1)																//if i is lower than sizeOfA - 1 =>
				m[i + 1] = a[i + 1] / bPrime[i];												//set m at index i+1
		}
		
		send(dPrime,2,MPI_COMM_WORLD);															//send the vector dPrime to the node 2
		send(bPrime,2,MPI_COMM_WORLD);															//send the vector bPrime to the node 2
		send(c,2,MPI_COMM_WORLD);																//send the vector c to the node 2
		x = recv(x, npes - 1, MPI_COMM_WORLD);													//stuck the node 1 until the other node complete
	} else if (myrank == 2){																	//if the function run on the third node
		dPrime = recv(dPrime,1,MPI_COMM_WORLD);													//receive the vector dPrime from the node 1
		bPrime = recv(bPrime,1,MPI_COMM_WORLD);													//receive the vector bPrime from the node 1
		c = recv(c,1,MPI_COMM_WORLD);															//receive the vector c from the node 1
		
		x[sizeOfA - 1] = dPrime[sizeOfA - 1] / bPrime[sizeOfA - 1];								//set x at index sizeOfA - 1
		for (int i = sizeOfA - 2; i >= 0; i--)													//create a loop in reverse, the index i start at sizeOfA and decrease until reach 0
			x[i] = (dPrime[i] - c[i] * x[i + 1]) / bPrime[i];									//set x at index i
		
		send(x,npes - 1,MPI_COMM_WORLD);														//send the vector x to the last node
		x = recv(x, npes - 1, MPI_COMM_WORLD);													//stuck the node 2 until the other node complete
	} else if (myrank == npes - 1){															//if the function run on the last node
		x = recv(x,2,MPI_COMM_WORLD);															//receive the vector x from the node 2
		send(x,0,MPI_COMM_WORLD);																//send the vector x to the first node
		send(x,1,MPI_COMM_WORLD);																//send the vector x to the second node
		send(x,2,MPI_COMM_WORLD);																//send the vector x to the third node
		return x;																				//return vector x
	}																									
}

/*!
 *  \brief Function f(x,t)
 *
 *  Return the value of f(x,t)
 *
 *  \param t the value of t
 *  \param x the value of x
 *  \return Return the value of f(x,t)
 */
 
double analyticalSolutionx(double t, double x) {												//declaration of the function f(x,t)
	if (x <= 50 + 250 * t || x >= 110+250*t)													//if x is lower than 50 + 250 * t or higher than 110 + 250 * t =>
		return 0;																				//return 0
	return 100 * (sin(PI*(x - 50 - 250 * t) / 60));												//if x is between 50 and 110 return the wave
}

/*!
 *  \brief Create the vector f(x,t) for the given time
 *
 *  Get the value of f(x) for each x at time t and put it in the return vector
 *
 *  \param t the value of t
 *  \return Return a vector of every f(x,t)
 */
 
vector<double> analyticalSolution(double t) {													//declaration of the function which create the vector with all the value of f(x) analytical
	vector <double> res;																		//decalartion of the vector of double
	double x = 0;																				//initiate the variable x at 0
	while (x <= xTot) {																			//while x is lower than xtot (400) =>
		res.push_back(analyticalSolutionx(t, x));												//get the value of f(x) analytical and add it to the vector res
		x += dx;																				//add dx (5) to x
	}
	return res;																					//return the vector with all the value of f(x) analytical
}

//Scheme definitions

/*!
 *  \brief Explicit Upwind FTBS function
 *
 *  Function called to calculate the function f at time n + 1 in serial
 *
 *  \param previousSolution A vector of f at time n
 *  \param Dt Delta t
 *  \return Return a vector of the solution at time n + 1
 */

vector<double> ExplicitUpwindFTBSSerial(vector <double> previousSolution, double Dt) {			//declaration of the only function of the class which return a vector of the value of f at n+1
	vector <double> res;																		//initiate a vector of double to res
	const double c = (double)u*Dt / dx;															//define the value of c
	res.push_back((1 - c)*previousSolution[0]);													//for the first one, the use of left boundary is require
	double x = dx;																				//define a double x as dx
	int index = 1;																				//define a int index to 1
	while (x <= xTot) {																			//create loop while x is lower than Xtot (400)
		res.push_back((1 - c)*previousSolution[index] + c*previousSolution[index - 1]);			//add the value of the scheme to the vector res
		index++;																				//add one to the index
		x += dx;																				//add the value of delta x to x
	}
	return res;																					//return the vector of double res
}

/*!
 *  \brief Explicit Upwind FTBS function
 *
 *  Function called to calculate the function f at time n + 1 in parallel
 *
 *  \param previousSolution A vector of f at time n
 *  \param Dt Delta t
 *  \return Return a vector of the solution at time n + 1
 */

vector<double> ExplicitUpwindFTBSPara(vector <double> previousSolution, double Dt) {

	vector <double> res;																		//initiate a vector of double to res
	const double c = (double)u*Dt / dx;															//define the value of c
	const int nbPoint = (xTot/dx) / npes;														//define the number of point each node will have to calculate
	
	if(myrank == 0){																			//if the function run on the first node =>
		res.push_back((1 - c)*previousSolution[0]);												//calculate the first value
		double x = dx;																			//define a double x as dx
		for(unsigned int index = 1; index < nbPoint; index++) {									//create loop
			res.push_back((1 - c)*previousSolution[index] + c*previousSolution[index - 1]);		//add the value of the scheme to the vector res
			x += dx;																			//add the value of delta x to x
		}
		send(res,1,MPI_COMM_WORLD);																//send the first part of the vector to the seconde node
	} else {																					//if the function run on the other node =>
		res = recv(res,myrank-1,MPI_COMM_WORLD);												//receive the previous part from the previous node
		double x = dx*nbPoint* myrank;															//define a double x as dx
		for(unsigned int index = nbPoint *(myrank); index < nbPoint*(1+myrank); index++) {		//create loop
			res.push_back((1 - c)*previousSolution[index] + c*previousSolution[index - 1]);		//add the value of the scheme to the vector res
			x += dx;																			//add the value of delta x to x
		}	
		if (myrank == npes-1){																	//if the function run on the last node =>
			return res;																			//return the vector res
		} else {																				//else =>
			send(res,myrank+1,MPI_COMM_WORLD);													//send the vector res to next node
		}
	}
}


/*!
 *  \brief Impplicit Upwind FTBS function
 *
 *  Function called to calculate the function f at time n + 1 in serial
 *
 *  \param previousSolution A vector of f at time n
 *  \param Dt Delta t
 *  \return Return a vector of the solution at time n + 1
 */

vector<double> ImplicitUpwindFTBSSerial(vector <double> previousSolution, double Dt) {			//declaration of the function of the Implicit FTCS
	//declaration of all the variable
	vector <double> resX;																		//create a vector of double resX
	Matrix A(previousSolution.size(), previousSolution.size());									//create a matrix which as the same size of the last solution
	double c = u * Dt / dx;																		//define the c variable

	//create the matrix A

	for (unsigned int i = 0; i < previousSolution.size(); i++) {								//create a loop to go trought all the column
		for (unsigned int j = 0;j < previousSolution.size();j++) {								//create a loop to go trought all the row
			if (i == j)																			//if we are on the diagonal =>
				A[i][j] = 1 + c;																//set the value to 1+c
			else if (j == i - 1)																//if we are on the lower diagonal =>
				A[i][j] = -c;																	//set the value to -c
			else																				//in every other case =>
				A[i][j] = 0; 																	//set the value to 0
		}
	}
	return ThomasAlgorithmSerial(A, previousSolution);											//call the Thomas algorithm with the A Matrix and the previous solution
}


/*!
 *  \brief Implicit Upwind FTBS function
 *
 *  Function called to calculate the function f at time n + 1 in parallel
 *
 *  \param previousSolution A vector of f at time n
 *  \param Dt Delta t
 *  \return Return a vector of the solution at time n + 1
 */

vector<double> ImplicitUpwindFTBSPara(vector <double> previousSolution, double Dt) {			//declaration of the function of the Implicit Upwind FTBS
	vector<double> res(previousSolution.size());												//initiate a vector of double to res
	const double c = (double)u*Dt / dx;															//define the value of c
	const int nbPoint = (xTot/dx) / npes;														//define the number of point each node will have to calculate
	const double k = ((-c)/(1 + c));															//define a constant k
	
	if(myrank == 0){																			//if the function run on the fist node =>
		res[0] = previousSolution[0];															//set the first element of the vector 
		for(unsigned int index = 1; index < nbPoint; index++) {									//create loop
			res[index] = (previousSolution[index] - k * res[index - 1]);						//set the element of the vector res
		}
	} else {																					//else =>
		res = recv(res,myrank-1,MPI_COMM_WORLD);												//receive the previous part the vector res
		double x = dx*nbPoint* myrank;															//define a double x as dx
		for(unsigned int index = nbPoint *(myrank); index < nbPoint*(1+myrank); index++) {		//create loop
			res[index] = (previousSolution[index] - k * res[index - 1]);						//set dPrime at index i
		}
	}	
	if (myrank == npes-1){																		//if the function is run the last node =>
		for (unsigned int i = 0; i < res.size(); i++){											//create a loop for all the element of the vector res
			res [i] = res[i] / (1+c);															//divide all the element of res by 1+c
		}
		return res;																				//return the vector res
	} else {																					//else =>
		send(res,myrank+1,MPI_COMM_WORLD);														//send the vector res to the next node
		for (unsigned int i = 0; i < res.size(); i++){											//create a loop for all the element of the vector res
			res [i] = res[i] / (1+c);															//divide all the element of res by 1+c
		}
		return res;																				//return the vector res
	}
}


/*!
 *  \brief Implicit FTCS function
 *
 *  Function called to calculate the function f at time n + 1 in serial
 *
 *  \param previousSolution A vector of f at time n
 *  \param Dt Delta t
 *  \return Return a vector of the solution at time n + 1
 */

vector<double> ImplicitFTCSSerial(vector <double> previousSolution, double Dt) {				//declaration of the function of the Implicit FTCS
	//declaration of all the variable
	vector <double> resX;																		//create a vector of double resX
	Matrix A(previousSolution.size(), previousSolution.size());									//create a matrix which as the same size of the last solution
	double c = u * Dt / dx;																		//define the c variable

	//create the matrix A

	for (unsigned int i = 0; i < previousSolution.size(); i++) {								//create a loop to go trought all the column
		for (unsigned int j = 0;j < previousSolution.size();j++) {								//create a loop to go trought all the row
			if (i == j)																			//if we are on the diagonal =>
				A[i][j] = 1;																	//set the value to 1
			else if (j == i - 1)																//if we are on the lower diagonal =>
				A[i][j] = -c / 2;																//set the value to -c/2
			else if (j == i + 1)																//if we are on the higher diagonal =>
				A[i][j] = c / 2;																//set the value to c/2
			else																				//in every other case =>
				A[i][j] = 0;																	//set the value to 0
		}
	}
	return ThomasAlgorithmSerial(A, previousSolution);											//call the Thomas algorithm with the A Matrix and the previous solution
}

/*!
 *  \brief Implicit FTCS function
 *
 *  Function called to calculate the function f at time n + 1 in parallel
 *
 *  \param previousSolution A vector of f at time n
 *  \param Dt Delta t
 *  \return Return a vector of the solution at time n + 1
 */

vector<double> ImplicitFTCSPara(vector <double> previousSolution, double Dt) {					//declaration of the function of the Implicit FTCS
	//declaration of all the variable
	vector <double> resX;																		//create a vector of double resX
	Matrix A(previousSolution.size(), previousSolution.size());									//create a matrix which as the same size of the last solution
	double c = u * Dt / dx;																		//define the c variable

	//create the matrix A

	for (unsigned int i = 0; i < previousSolution.size(); i++) {								//create a loop to go trought all the column
		for (unsigned int j = 0;j < previousSolution.size();j++) {								//create a loop to go trought all the row
			if (i == j)																			//if we are on the diagonal =>
				A[i][j] = 1;																	//set the value to 1
			else if (j == i - 1)																//if we are on the lower diagonal =>
				A[i][j] = -c / 2;																//set the value to -c/2
			else if (j == i + 1)																//if we are on the higher diagonal =>
				A[i][j] = c / 2;																//set the value to c/2
			else																				//in every other case =>
				A[i][j] = 0;																	//set the value to 0
		}
	}
	return ThomasAlgorithmPara(A, previousSolution);											//call the Thomas algorithm with the A Matrix and the previous solution																//return the value of resX
}


/*!
 *  \brief Main function
 *
 *  Two for loop to go trough each case
 *	Delta t is affected in this function.
 *	call each function to calculate the correct vector
 *
 *  \return If everythings works correctly, return 0
 */

int main(int argc, char *argv []) {

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	setprecision(10);
	
	double time1 , time2;
	for(unsigned int parallele = 0; parallele <= 1; parallele++){
		if(parallele == 0 ) {printf("The programm run in serial");} else {printf("The programm run in parallel");}
		//1 for EXPLICIT UPWIND FTBS  2 for IMPLICIT UPWIND FTBS 3 for IMPLICIT FTCS
		for (unsigned int type = 1; type <= 3 ; type++){											//declare a int type
		//1 for dt = 0.002 2 for dt = 0.001 3 for dt = 0.0005
			for (unsigned int dtIndex = 1; dtIndex <= 3; dtIndex++){								//declare an int dtIndex
				double Dt;																			//declare a double dt
				if (dtIndex == 1)																	//if dtIndex equal 1 =>
					Dt = 0.002;																		//dt equal 0.002
				else if (dtIndex == 2)																//if dtIndex equal 2 =>
					Dt = 0.001;																		//dt equal 0.001
				else if (dtIndex == 3)																//if dtIndex equal 3 =>
					Dt = 0.0005;																	//dt equal 0.0005
			
				if (type == 1) {																	//if the type is 1 =>
					if(myrank == npes - 1){														//if the function run the last node =>
						printf("----------- EXPLICIT UPWIND FTBS -------------------\n");			//Print a the scheme choosen
						printf("Delta t: %f \n", Dt);												//print the value of Delta t
					}
					time1 = MPI_Wtime();															//get the start time
					// cout << "\nResult for dt = " << Dt << ": \n";								//print the value of the delta t
					vector<double> solution;														//intiate a vector of double solution
					// cout << "DELTA X: " << dx << "---- X TOT: " << xTot << "\n";					//print the value of delta x
					solution = finit();																//affect the result of the first solution to the vector solution
					vector <double> error(solution.size());											//create a vector of double error which has the same size of the solution
					double n = 0;																	//initate a double n at 0
					for(unsigned int nbLoop = 0 ; nbLoop < 5 ; nbLoop ++){							//creat a loop of 5 iteration
						double nLoop = 0;															//define the n of the loop to 0
						while (nLoop <= 0.1) {														//while nloop is lower than 0.1
							if(Parallele == 1){														//if we want parallel calculate =>
								solution = ExplicitUpwindFTBSPara(solution, Dt);					//call the function which calculate the solution at n+1
							} else {																//else =>
								solution = ExplicitUpwindFTBSSerial(solution, Dt);					//call the function which calculate the solution at n+1
							}
							n += Dt;																//add delta t to n
							nLoop += Dt;															//add delta t to nLoop
						}
						if(myrank == npes - 1){													//if the function run the last node =>
							// cout << "NUMERICAL solution for n = " << n << ":\n";					//print a title
							//printf("n = %f \n", n);												//print the value of n and an enter
							//showVector(solution);													//print the vector of the solution
							/*cout << "\n";															//print an enter
							cout << "ANALYTICAL solution for n = " << n << ":\n";					//print a title before print a solution
							showVector(analyticalSolution(n));										//print the vector of analytical solution
							cout << "\n";	*/														//print an enter
							/*
							for (unsigned int i = 0; i < solution.size();i++){						//create a loop with i which will be used as an index
								error[i] = solution[i] - analyticalSolution(n)[i];					//set to the vector at index i the diffrence between the solution and the analytical solution
							}
						
							cout << "ERROR for n = " << n << ":\n";									//print a title before print vector
							showVector(error);														//show the vector error
							*/
						}
					}
					if(myrank == npes - 1){														//if the function run the last node =>
						time2 = MPI_Wtime();														//get the finish time
						printf("The time taken is %f\n", time2-time1);								//print the duration
					}
				} else if (type == 2) {																//if the type is 2 =>
					// cout << "\nResult for dt = " << Dt << ": \n";								//print the value of delta t
					vector<double> solution;														//create a vector of double solution
					// cout << "DELTA X: " << fx.dx << "---- X TOT: " << fx.xTot << "\n";			//print the value of delta x
					if (myrank == npes - 1) {														//if the function run the last node =>
						printf("----------- IMPLICIT UPWIND FTBS -------------------\n");			//print a title
						printf("Delta t: %f \n", Dt);												//print the value of Delta t
					}
					time1 = MPI_Wtime();															//get the start time
					solution = finit();																//affect the result of the first solution to the vector solution
					//if(myrank == npes - 1 )showVector(solution);									//affect the vector from the function finit in commons to solution
					vector <double> error(solution.size());											//create a vector of double error which has the same size as solution 
					double n = 0;																	//create a double n equal to 0
					for(unsigned int nbLoop = 0 ; nbLoop < 5 ; nbLoop ++){							//creat a loop of 5 iteration
						double nLoop = 0;															//define the n of the loop to 0
						while (nLoop <= 0.1) {														//create loop until nLoop <0.1
							if(Parallele == 1){														//if we want parallel calculate =>
								solution = ImplicitUpwindFTBSPara(solution, Dt);					//call the function which calculate the solution at n+1
							}else{																	//else =>
								solution = ImplicitUpwindFTBSSerial(solution, Dt);					//call the function which calculate the solution at n+1
							}
							n += Dt;																//add delta t to n
							nLoop += Dt;															//add delta t to nLoop
						}
						
						if(myrank == npes - 1){													//if the function run the last node =>
							// cout << "NUMERICAL solution for n = " << n << ":\n";					//print a title
							//printf("n = %f \n", n);												//print the value of n and an enter
							//showVector(solution);													//print the vector of the solution
							/*cout << "\n";															//print an enter
							cout << "ANALYTICAL solution for n = " << n << ":\n";					//print a title before print a solution
							showVector(analyticalSolution(n));										//print the vector of analytical solution
							cout << "\n";	*/														//print an enter
							/*
							for (unsigned int i = 0; i < solution.size();i++){						//create a loop with i which will be used as an index
								error[i] = solution[i] - analyticalSolution(n)[i];					//set to the vector at index i the diffrence between the solution and the analytical solution
							}
						
							cout << "ERROR for n = " << n << ":\n";									//print a title before print vector
							showVector(error);														//show the vector error
							*/
						}
					}
					if(myrank == npes - 1){														//if the function run the last node =>
						time2 = MPI_Wtime();														//get the finish time
						printf("The time taken is %f\n", time2-time1);								//print the duration
					}																				//call the function resultDt with delta t
				} else if (type == 3) {																//if the type is 3 =>
					// cout << "\nResult for dt = " << Dt << ": \n";								//print the value of delta t
					vector<double> solution;														//create a vector of double solution
					// cout << "DELTA X: " << fx.dx << "---- X TOT: " << fx.xTot << "\n";			//print the value of delta x
					if (myrank == npes - 1) {														//if the function run the last node =>
						printf("----------- IMPLICIT FTCS -------------------\n");					//print a title
						printf("Delta t: %f \n", Dt);												//print the value of Delta t
					}
					time1 = MPI_Wtime();															//get the start time
					solution = finit();																//affect the vector from the function finit in commons to solution
					vector <double> error(solution.size());											//create a vector of double error which has the same size as solution 
					double n = 0;																	//create a double n equal to 0
					for(unsigned int nbLoop = 0 ; nbLoop < 5 ; nbLoop ++){							//creat a loop of 5 iteration
						double nLoop = 0;															//define the n of the loop to 0
						while (nLoop <= 0.1) {														//create loop until nLoop <0.1
							if(Parallele == 1){														//if we want parallel calculate =>
								solution = ImplicitFTCSPara(solution, Dt);							//call the function which calculate the solution at n+1
							} else {																//else =>
								solution = ImplicitFTCSSerial(solution, Dt);						//call the function which calculate the solution at n+1
							}
							n += Dt;																//add delta t to n
							nLoop += Dt;															//add delta t to nLoop
						}
						if(myrank == npes - 1){													//if the function run the last node =>
							// cout << "NUMERICAL solution for n = " << n << ":\n";					//print a title
							//printf("n = %f \n", n);												//print the value of n and an enter
							//showVector(solution);													//print the vector of the solution
							/*cout << "\n";															//print an enter
							cout << "ANALYTICAL solution for n = " << n << ":\n";					//print a title before print a solution
							showVector(analyticalSolution(n));										//print the vector of analytical solution
							cout << "\n";	*/														//print an enter
							
							/*for (unsigned int i = 0; i < solution.size();i++){					//create a loop with i which will be used as an index
								error[i] = solution[i] - analyticalSolution(n)[i];					//set to the vector at index i the diffrence between the solution and the analytical solution
							}

							cout << "ERROR for n = " << n << ":\n";									//print a title before print vector
							showVector(error);														//show the vector error
							*/
						}
					}
					if(myrank == npes - 1){														//if the function run the last node =>
						time2 = MPI_Wtime();														//get the finish time
						printf("The time taken is %f\n", time2-time1);								//print the duration
					}
				}
			}
		}
	}
	return 0;																					//if everything works correctly return 0
}
