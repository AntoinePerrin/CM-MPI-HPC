#include <mpi.h>
#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include "matrix.h"

using namespace std;


int npes, myrank, nNode = 4, Parallele = 0;

const double dx = 0.5;
const double u = 250;
const double xTot = 400;
const double PI = 3.14;



void send(vector<double> const& vector, int dest, MPI_Comm comm){
	unsigned len = vector.size();
	MPI_Send(&len, 1, MPI_UNSIGNED, dest, 1, comm);
	if(len !=0){
		MPI_Send(vector.data(), len, MPI_DOUBLE, dest, 1, comm);
	}
}

vector<double> recv(vector<double>& res, int src, MPI_Comm comm){
	MPI_Status status;
	unsigned len;
	MPI_Recv(&len, 1, MPI_UNSIGNED, src, 1, comm, &status);
	if(len !=0){
		res.resize(len);
		MPI_Recv(res.data(), len, MPI_DOUBLE, src, 1, comm, &status);
	}
	return res;
}

void showVector(vector <double> vect) {									//decalaration of the function which can show all the value of a vector of double
	for (unsigned int i = 0; i < vect.size(); i++) {						//create a loop with i which will be used as an index
		printf("%0.1f \n", vect[i]);								//print the value of the vector at the index i
	}
	printf("\n");											//print a enter for more clarity
}

double finitx(double x) {										//declaration of the function f(x)
	if (x <= 50 || x >= 110)									//if x is lower than 50 or higher than 110 =>
		return 0.0;										//return 0
	return 100 * (sin(PI*(x - 50) / 60));								//if x is between 50 and 110 return the wave
}

vector<double> finit() {										//declaration of the function which create the vector with all the value of f(x)
	vector <double> res;										//decalartion of the vector of double
	double x = 0;											//initiate the variable x at 0
	while (x <= xTot) {										//while x is lower than xtot (400) =>
		res.push_back(finitx(x));								//get the value of f(x) and add it to the vector res
		x += dx;										//add dx (5) to x
	}
	return res;											//return the vector with all the value of f(x)
}

vector <double> ThomasAlgorithmSerial(Matrix A, vector <double> f) {							//declare the function which handle the thomas algorithm
	int sizeOfA = A.getNcols();																			//declare an int x to the number of columunm of the matrix A
	vector <double> x(sizeOfA);																			//declare a vector x
	vector <double> m(sizeOfA);																			//declare a vector m
	vector <double> bPrime(sizeOfA);																	//declare a vector bPrime
	vector <double> dPrime(sizeOfA);																	//declare a vector dPrime

	// We create vectors from the non null coefficients of the matrix
	vector <double> a(sizeOfA);																			//declare a vector a
	vector <double> b(sizeOfA);																			//declare a vector b
	vector <double> c(sizeOfA);																			//declare a vector c
																										//all of this vector have the same size which is sizeofA
	for (int i = 0; i < sizeOfA; i++) {																	//create a loop with an index i until i reach the sizeofA
		b[i] = A[i][i];																					//set the vector b at the index i the value of the matrix A at [i][i]
		if (i > 0)																						//if i higher than 0 =>
			a[i] = A[i][i - 1];																			//set the vector a at the index i the value of the matrix A at [i][i-1]
		if (i < sizeOfA-1)																				//if i is lower than sizeOfA - 1 =>
			c[i] = A[i][i + 1];																			//set the vector c at the index i the value of the matrix A at [i][i+1]
	}
	a[0] = 0;																							//set the first value of a to 0
	c[sizeOfA - 1] = 0;																					//set last value of c to 0


	// First values of the vectors
	m[0] = 0;																							//set the first value of m to 0
	m[1] = a[1] / b[0];																					//set the second value of m to the second value of a divided by the first value of b
	bPrime[0] = b[0];																					//set the first value of bPrime to the fisrt value of b
	dPrime[0] = f[0];																					//set the first value of dPrime to the fisrt value of f

	// We calculate all values of the vectors
	for (int i = 1; i < sizeOfA; i++) {																	//create a loop with an index i until i reach sizeOfA
		bPrime[i] = b[i] - m[i] * c[i - 1];																//set bPrime at index i
		dPrime[i] = f[i] - m[i] * dPrime[i - 1];														//set dPrime at index i
		if (i < sizeOfA - 1)																			//if i is lower than sizeOfA - 1 =>
			m[i + 1] = a[i + 1] / bPrime[i];															//set m at index i+1
	}

	x[sizeOfA - 1] = dPrime[sizeOfA - 1] / bPrime[sizeOfA - 1];											//set x at index sizeOfA - 1
	for (int i = sizeOfA - 2; i >= 0; i--)																//create a loop in reverse, the index i start at sizeOfA and decrease until reach 0
		x[i] = (dPrime[i] - c[i] * x[i + 1]) / bPrime[i];												//set x at index i

	return x;																							//return vector x
}

vector <double> ThomasAlgorithmPara(Matrix A, vector <double> f) {							//declare the function which handle the thomas algorithm
	int sizeOfA = A.getNcols();																			//declare an int x to the number of columunm of the matrix A
	vector <double> x(sizeOfA);																			//declare a vector x
	vector <double> m(sizeOfA);																			//declare a vector m
	vector <double> bPrime(sizeOfA);																	//declare a vector bPrime
	vector <double> dPrime(sizeOfA);																	//declare a vector dPrime

	// We create vectors from the non null coefficients of the matrix
	vector <double> a(sizeOfA);																			//declare a vector a
	vector <double> b(sizeOfA);																			//declare a vector b
	vector <double> c(sizeOfA);																			//declare a vector c
	if (myrank == 0){																									//all of this vector have the same size which is sizeofA
		for (int i = 0; i < sizeOfA; i++) {																	//create a loop with an index i until i reach the sizeofA
			b[i] = A[i][i];																					//set the vector b at the index i the value of the matrix A at [i][i]
			if (i > 0){																		//if i higher than 0 =>
				a[i] = A[i][i - 1];																			//set the vector a at the index i the value of the matrix A at [i][i-1]
			}
			if (i < sizeOfA-1){																			//if i is lower than sizeOfA - 1 =>
				c[i] = A[i][i + 1];																		//set the vector c at the index i the value of the matrix A at [i][i+1]
			}
		}
		a[0] = 0;																							//set the first value of a to 0
		c[sizeOfA - 1] = 0;
		send(a,1,MPI_COMM_WORLD);
		send(b,1,MPI_COMM_WORLD);
		send(c,1,MPI_COMM_WORLD);
		x = recv(x, nNode - 1, MPI_COMM_WORLD);																				//set last value of c to 0
	} else if (myrank == 1){
		a = recv(a,0,MPI_COMM_WORLD);
		b = recv(b,0,MPI_COMM_WORLD);
		c = recv(c,0,MPI_COMM_WORLD);
		// First values of the vectors
		m[0] = 0;																							//set the first value of m to 0
		m[1] = a[1] / b[0];																					//set the second value of m to the second value of a divided by the first value of b
		bPrime[0] = b[0];																					//set the first value of bPrime to the fisrt value of b
		dPrime[0] = f[0];																					//set the first value of dPrime to the fisrt value of f
	
		// We calculate all values of the vectors
		for (int i = 1; i < sizeOfA; i++) {																	//create a loop with an index i until i reach sizeOfA
			bPrime[i] = b[i] - m[i] * c[i - 1];																//set bPrime at index i
			dPrime[i] = f[i] - m[i] * dPrime[i - 1];														//set dPrime at index i
			if (i < sizeOfA - 1)																			//if i is lower than sizeOfA - 1 =>
				m[i + 1] = a[i + 1] / bPrime[i];															//set m at index i+1
		}
		
		send(dPrime,2,MPI_COMM_WORLD);
		send(bPrime,2,MPI_COMM_WORLD);
		send(c,2,MPI_COMM_WORLD);
		x = recv(x, nNode - 1, MPI_COMM_WORLD);
	} else if (myrank == 2){
		dPrime = recv(dPrime,1,MPI_COMM_WORLD);
		bPrime = recv(bPrime,1,MPI_COMM_WORLD);
		c = recv(c,1,MPI_COMM_WORLD);

		
		x[sizeOfA - 1] = dPrime[sizeOfA - 1] / bPrime[sizeOfA - 1];											//set x at index sizeOfA - 1
		for (int i = sizeOfA - 2; i >= 0; i--)																//create a loop in reverse, the index i start at sizeOfA and decrease until reach 0
			x[i] = (dPrime[i] - c[i] * x[i + 1]) / bPrime[i];												//set x at index i
		
		send(x,nNode - 1,MPI_COMM_WORLD);
		x = recv(x, nNode - 1, MPI_COMM_WORLD);
	} else if (myrank == nNode - 1){
		x = recv(x,2,MPI_COMM_WORLD);
		send(x,0,MPI_COMM_WORLD);
		send(x,1,MPI_COMM_WORLD);
		send(x,2,MPI_COMM_WORLD);
		return x;
	}																						//return vector x
}

double analyticalSolutionx(double t, double x) {							//declaration of the function f(x,t)
	if (x <= 50 + 250 * t || x >= 110+250*t)							//if x is lower than 50 + 250 * t or higher than 110 + 250 * t =>
		return 0;														//return 0
	return 100 * (sin(PI*(x - 50 - 250 * t) / 60));							//if x is between 50 and 110 return the wave
}

vector<double> analyticalSolution(double t) {								//declaration of the function which create the vector with all the value of f(x) analytical
	vector <double> res;										//decalartion of the vector of double
	double x = 0;											//initiate the variable x at 0
	while (x <= xTot) {										//while x is lower than xtot (400) =>
		res.push_back(analyticalSolutionx(t, x));						//get the value of f(x) analytical and add it to the vector res
		x += dx;										//add dx (5) to x
	}
	return res;											//return the vector with all the value of f(x) analytical
}

vector<double> ExplicitUpwindFTBSSerial(vector <double> previousSolution, double Dt) {			//declaration of the only function of the class which return a vector of the value of f at n+1
	vector <double> res;										//initiate a vector of double to res
	const double c = (double)u*Dt / dx;								//define the value of c
	res.push_back((1 - c)*previousSolution[0]);							//for the first one, the use of left boundary is require
	double x = dx;											//define a double x as dx
	int index = 1;											//define a int index to 1
	while (x <= xTot) {										//create loop while x is lower than Xtot (400)
		res.push_back((1 - c)*previousSolution[index] + c*previousSolution[index - 1]);		//add the value of the scheme to the vector res
		index++;										//add one to the index
		x += dx;										//add the value of delta x to x
	}
	return res;											//return the vector of double res
}

vector<double> ExplicitUpwindFTBSPara(vector <double> previousSolution, double Dt) {

	vector <double> res;										//initiate a vector of double to res
	const double c = (double)u*Dt / dx;								//define the value of c
	const int nbPoint = (xTot/dx) / nNode;
	
	if(myrank == 0){
		res.push_back((1 - c)*previousSolution[0]);
		double x = dx;											//define a double x as dx
		for(unsigned int index = 1; index < nbPoint; index++) {										//create loop while x is lower than Xtot (400)
			res.push_back((1 - c)*previousSolution[index] + c*previousSolution[index - 1]);		//add the value of the scheme to the vector res										//add one to the index
			x += dx;										//add the value of delta x to x
		}
		send(res,1,MPI_COMM_WORLD);
	} else {
		res = recv(res,myrank-1,MPI_COMM_WORLD);
		double x = dx*nbPoint* myrank;											//define a double x as dx
		for(unsigned int index = nbPoint *(myrank); index < nbPoint*(1+myrank); index++) {										//create loop while x is lower than Xtot (400)
			res.push_back((1 - c)*previousSolution[index] + c*previousSolution[index - 1]);		//add the value of the scheme to the vector res										//add one to the index
			x += dx;										//add the value of delta x to x
		}	
		if (myrank == nNode-1){
			return res;
		} else {
			send(res,myrank+1,MPI_COMM_WORLD);
		}
	}
}

vector<double> ImplicitUpwindFTBSSerial(vector <double> previousSolution, double Dt) {	//declaration of the function of the Implicit FTCS
	//declaration of all the variable
	vector <double> resX;																			//create a vector of double resX
	Matrix A(previousSolution.size(), previousSolution.size());										//create a matrix which as the same size of the last solution
	double c = u * Dt / dx;																	//define the c variable

	//create the matrix A

	for (unsigned int i = 0; i < previousSolution.size(); i++) {									//create a loop to go trought all the column
		for (unsigned int j = 0;j < previousSolution.size();j++) {									//create a loop to go trought all the row
			if (i == j)																						//if we are on the diagonal =>
				A[i][j] = 1 + c;																			//set the value to 1+c
			else if (j == i - 1)																			//if we are on the lower diagonal =>
				A[i][j] = -c;																				//set the value to -c
			else																							//in every other case =>
				A[i][j] = 0; 																		//set the value to 0
		}
	}
	return ThomasAlgorithmSerial(A, previousSolution);																		//call the Thomas algorithm in the Implicit Scheme class with the A Matrix and the previous solution																//return the value of resX
}

vector<double> ImplicitUpwindFTBSPara(vector <double> previousSolution, double Dt) {		//declaration of the function of the Implicit Upwind FTBS																								//create a variable of common to use all the common variable
	
	vector<double> res(previousSolution.size());										//initiate a vector of double to res
	const double c = (double)u*Dt / dx;								//define the value of c
	const int nbPoint = (xTot/dx) / nNode;
	const double k = ((-c)/(1 + c));
	
	if(myrank == 0){
		res[0] = previousSolution[0];
		for(unsigned int index = 1; index < nbPoint; index++) {										//create loop while x is lower than Xtot (400)
			res[index] = (previousSolution[index] - k * res[index - 1]);														//set dPrime at index i
		}
	} else {
		res = recv(res,myrank-1,MPI_COMM_WORLD);
		double x = dx*nbPoint* myrank;											//define a double x as dx
		for(unsigned int index = nbPoint *(myrank); index < nbPoint*(1+myrank); index++) {										//create loop while x is lower than Xtot (400)
			res[index] = (previousSolution[index] - k * res[index - 1]);														//set dPrime at index i

		}
	}	
	if (myrank == nNode-1){
		for (unsigned int i = 0; i < res.size(); i++){
			res [i] = res[i] / (1+c);
		}
		return res;
	} else {
		send(res,myrank+1,MPI_COMM_WORLD);
		for (unsigned int i = 0; i < res.size(); i++){
			res [i] = res[i] / (1+c);
		}
		return res;
	}
}



vector<double> ImplicitFTCSSerial(vector <double> previousSolution, double Dt) {	//declaration of the function of the Implicit FTCS
	//declaration of all the variable
	vector <double> resX;																			//create a vector of double resX
	Matrix A(previousSolution.size(), previousSolution.size());										//create a matrix which as the same size of the last solution
	double c = u * Dt / dx;																	//define the c variable

	//create the matrix A

	for (unsigned int i = 0; i < previousSolution.size(); i++) {									//create a loop to go trought all the column
		for (unsigned int j = 0;j < previousSolution.size();j++) {									//create a loop to go trought all the row
			if (i == j)																				//if we are on the diagonal =>
				A[i][j] = 1;																		//set the value to 1
			else if (j == i - 1)																	//if we are on the lower diagonal =>
				A[i][j] = -c / 2;																	//set the value to -c/2
			else if (j == i + 1)																	//if we are on the higher diagonal =>
				A[i][j] = c / 2;																	//set the value to c/2
			else																					//in every other case =>
				A[i][j] = 0;																		//set the value to 0
		}
	}
	return ThomasAlgorithmSerial(A, previousSolution);																		//call the Thomas algorithm in the Implicit Scheme class with the A Matrix and the previous solution																//return the value of resX
}


vector<double> ImplicitFTCSPara(vector <double> previousSolution, double Dt) {	//declaration of the function of the Implicit FTCS
	//declaration of all the variable
	vector <double> resX;																			//create a vector of double resX
	Matrix A(previousSolution.size(), previousSolution.size());										//create a matrix which as the same size of the last solution
	double c = u * Dt / dx;																	//define the c variable

	//create the matrix A

	for (unsigned int i = 0; i < previousSolution.size(); i++) {									//create a loop to go trought all the column
		for (unsigned int j = 0;j < previousSolution.size();j++) {									//create a loop to go trought all the row
			if (i == j)																				//if we are on the diagonal =>
				A[i][j] = 1;																		//set the value to 1
			else if (j == i - 1)																	//if we are on the lower diagonal =>
				A[i][j] = -c / 2;																	//set the value to -c/2
			else if (j == i + 1)																	//if we are on the higher diagonal =>
				A[i][j] = c / 2;																	//set the value to c/2
			else																					//in every other case =>
				A[i][j] = 0;																		//set the value to 0
		}
	}
	return ThomasAlgorithmPara(A, previousSolution);																		//call the Thomas algorithm in the Implicit Scheme class with the A Matrix and the previous solution																//return the value of resX
}


/*!
 *  \brief Main function
 *
 *  Let the user choosen the scheme and the delta t.
 *	Delta t is affected in this function.
 *	call showResultsInConsole(const int type, const double dt)
 *
 *  \return If everythings works correctly, return 1
 */

int main(int argc, char *argv []) {

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	setprecision(10);
	
	double time1 , time2;

	//1 for EXPLICIT UPWIND FTBS  2 for IMPLICIT UPWIND FTBS 3 for IMPLICIT FTCS
	for (unsigned int type = 1; type <= 3 ; type++){											//declare a int type
	//1 for dt = 0.002 2 for dt = 0.001 3 for dt = 0.0005
		for (unsigned int dtIndex = 1; dtIndex <= 3; dtIndex++){											//declare an int dtIndex
			double Dt;											//declare a double dt
			if (dtIndex == 1)										//if dtIndex equal 1 =>
				Dt = 0.002;										//dt equal 0.002
			else if (dtIndex == 2)										//if dtIndex equal 2 =>
				Dt = 0.001;										//dt equal 0.001
			else if (dtIndex == 3)										//if dtIndex equal 3 =>
				Dt = 0.0005;										//dt equal 0.0005
		
			if (type == 1) {													//if the type is 1 =>
				if(myrank == nNode - 1){
					printf("----------- EXPLICIT UPWIND FTBS -------------------\n");		//Print a the scheme choosen
					printf("Delta t: %f \n", Dt);							//print the value of Delta t
				}
				time1 = MPI_Wtime();
				// cout << "\nResult for dt = " << Dt << ": \n";					//print the value of the delta t
				vector<double> solution;								//intiate a vector of double solution
				// cout << "DELTA X: " << dx << "---- X TOT: " << xTot << "\n";				//print the value of delta x
				solution = finit();									//affect the result of the first solution to the vector solution
				vector <double> error(solution.size());							//create a vector of double error which has the same size of the solution
				double n = 0;										//initate a double n at 0
				for(unsigned int nbLoop = 0 ; nbLoop < 5 ; nbLoop ++){
					double nLoop = 0;
					while (nLoop <= 0.1) {
						if(Parallele == 1){									//create loop until n <0.5
							solution = ExplicitUpwindFTBSPara(solution, Dt);					//call the function which calculate the solution at n+1
						}else {
							solution = ExplicitUpwindFTBSSerial(solution, Dt);					//call the function which calculate the solution at n+1
						}
						n += Dt;									//add delta t to n
						nLoop += Dt;
					}
					if(myrank == nNode - 1){
						// cout << "NUMERICAL solution for n = " << n << ":\n";				//print a title
						//printf("n = %f \n", n);								//print the value of n and an enter
						//showVector(solution);								//print the vector of the solution
						/*cout << "\n";									//print an enter
						cout << "ANALYTICAL solution for n = " << n << ":\n";				//print a title before print a solution
						showVector(analyticalSolution(n));						//print the vector of analytical solution
				 		cout << "\n";	*/								//print an enter
						/*
						for (unsigned int i = 0; i < solution.size();i++){				//create a loop with i which will be used as an index
							error[i] = solution[i] - analyticalSolution(n)[i];			//set to the vector at index i the diffrence between the solution and the analytical solution
						}
					
						cout << "ERROR for n = " << n << ":\n";						//print a title before print vector
						showVector(error);								//show the vector error
						*/
					}
				}
				if(myrank == nNode - 1){
					time2 = MPI_Wtime();
					printf("The time taken is %f\n", time2-time1);
				}												//call the function resultDt with delta t
			} else if (type == 2) {													//if the type is 2 =>

				// cout << "\nResult for dt = " << Dt << ": \n";													//print the value of delta t
				vector<double> solution;																		//create a vector of double solution
				// cout << "DELTA X: " << fx.dx << "---- X TOT: " << fx.xTot << "\n";							//print the value of delta x
				if (myrank == nNode - 1) {
					printf("----------- IMPLICIT UPWIND FTBS -------------------\n");
					printf("Delta t: %f \n", Dt);
				}								//print the value of Delta t
				time1 = MPI_Wtime();
				solution = finit();
				//if(myrank == nNode - 1 )showVector(solution);																			//affect the vector from the function finit in commons to solution
				vector <double> error(solution.size());															//create a vector of double error which has the same size as solution 
				double n = 0;																					//create a double n equal to 0
				for(unsigned int nbLoop = 0 ; nbLoop < 5 ; nbLoop ++){
					double nLoop = 0;
					while (nLoop <= 0.1) {									//create loop until n <0.5
						if(Parallele == 1){
							solution = ImplicitUpwindFTBSPara(solution, Dt);					//call the function which calculate the solution at n+1
						}else{
							solution = ImplicitUpwindFTBSSerial(solution, Dt);					//call the function which calculate the solution at n+1
						}
						n += Dt;									//add delta t to n
						nLoop += Dt;
					}
					
					if(myrank == nNode - 1){
						// cout << "NUMERICAL solution for n = " << n << ":\n";				//print a title
						//printf("n = %f \n", n);								//print the value of n and an enter
						//showVector(solution);								//print the vector of the solution
						/*cout << "\n";									//print an enter
						cout << "ANALYTICAL solution for n = " << n << ":\n";				//print a title before print a solution
						showVector(analyticalSolution(n));						//print the vector of analytical solution
						cout << "\n";	*/								//print an enter
						/*
						for (unsigned int i = 0; i < solution.size();i++){				//create a loop with i which will be used as an index
							error[i] = solution[i] - analyticalSolution(n)[i];			//set to the vector at index i the diffrence between the solution and the analytical solution
						}
					
						cout << "ERROR for n = " << n << ":\n";						//print a title before print vector
						showVector(error);								//show the vector error
						*/
					}
				}
				if(myrank == nNode - 1){
					time2 = MPI_Wtime();
					printf("The time taken is %f\n", time2-time1);
				}												//call the function resultDt with delta t
			} else if (type == 3) {													//if the type is 4 =>
				// cout << "\nResult for dt = " << Dt << ": \n";													//print the value of delta t
				vector<double> solution;																		//create a vector of double solution
				// cout << "DELTA X: " << fx.dx << "---- X TOT: " << fx.xTot << "\n";							//print the value of delta x
				if (myrank == nNode - 1) {
					printf("----------- IMPLICIT FTCS -------------------\n");
					printf("Delta t: %f \n", Dt);								//print the value of Delta t
				}
				time1 = MPI_Wtime();
				solution = finit();																			//affect the vector from the function finit in commons to solution
				vector <double> error(solution.size());															//create a vector of double error which has the same size as solution 
				double n = 0;																					//create a double n equal to 0
				for(unsigned int nbLoop = 0 ; nbLoop < 5 ; nbLoop ++){
					double nLoop = 0;
					while (nLoop <= 0.1) {									//create loop until n <0.5
						if(Parallele == 1){
							solution = ImplicitFTCSPara(solution, Dt);					//call the function which calculate the solution at n+1
						} else {
							solution = ImplicitFTCSSerial(solution, Dt);
						}
						n += Dt;									//add delta t to n
						nLoop += Dt;
					}
					if(myrank == nNode - 1){
						// cout << "NUMERICAL solution for n = " << n << ":\n";				//print a title
						//printf("n = %f \n", n);								//print the value of n and an enter
						//showVector(solution);								//print the vector of the solution
						/*cout << "\n";									//print an enter
						cout << "ANALYTICAL solution for n = " << n << ":\n";				//print a title before print a solution
						showVector(analyticalSolution(n));						//print the vector of analytical solution
						cout << "\n";	*/								//print an enter
						
						/*for (unsigned int i = 0; i < solution.size();i++){				//create a loop with i which will be used as an index
							error[i] = solution[i] - analyticalSolution(n)[i];			//set to the vector at index i the diffrence between the solution and the analytical solution
						}
					
						cout << "ERROR for n = " << n << ":\n";						//print a title before print vector
						showVector(error);								//show the vector error
						*/
					}
				}
				if(myrank == nNode - 1){
					time2 = MPI_Wtime();
					printf("The time taken is %f\n", time2-time1);
				}												//call the function resultDt with delta t
			}
		}
	}
	return 0;											//if everything works correctly return 0
}
