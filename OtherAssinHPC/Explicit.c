#include <mpi.h>
#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;


int npes, myrank;

MPI_Status status;
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &npes);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

const double dx = 0.5;
const double u = 250;
const double xTot = 400;
const double PI = 3.14;


void showVector(vector <double> vect) {									//decalaration of the function which can show all the value of a vector of double
	for (unsigned int i = 0; i < vect.size(); i++) {						//create a loop with i which will be used as an index
		cout << vect[i] << "; ";								//print the value of the vector at the index i
	}
	cout << "\n";											//print a enter for more clarity
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


vector<double> ExplicitScheme_nplus1(vector <double> previousSolution, double Dt) {			//declaration of the only function of the class which return a vector of the value of f at n+1
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


void resultDtExp(const double Dt) {									//declare the only function of this class which handle the calculation and the print of the scheme
	// cout << "\nResult for dt = " << Dt << ": \n";						//print the value of the delta t
	vector<double> solution;								//intiate a vector of double solution
	// cout << "DELTA X: " << dx << "---- X TOT: " << xTot << "\n";				//print the value of delta x
	cout << "Delta t: " << Dt << "\n";							//print the value of Delta t
	solution = finit();									//affect the result of the first solution to the vector solution
	vector <double> error(solution.size());							//create a vector of double error which has the same size of the solution
	double n = 0;										//initate a double n at 0
	while (n <= 0.52) {									//create loop until n <0.5
		// cout << "NUMERICAL solution for n = " << n << ":\n";				//print a title
		cout << "n = " << n << "\n";							//print the value of n and an enter
		showVector(solution);								//print the vector of the solution
		/*
		cout << "\n";									//print an enter
		cout << "ANALYTICAL solution for n = " << n << ":\n";				//print a title before print a solution
		showVector(analyticalSolution(n));						//print the vector of analytical solution
		cout << "\n";*/									//print an enter
		for (unsigned int i = 0; i < solution.size();i++){				//create a loop with i which will be used as an index
			error[i] = solution[i] - analyticalSolution(n)[i];			//set to the vector at index i the diffrence between the solution and the analytical solution
		}
		/*
		cout << "ERROR for n = " << n << ":\n";						//print a title before print vector
		showVector(error);								//show the vector error
		*/
		solution = ExplicitScheme_nplus1(solution, Dt);					//call the function which calculate the solution at n+1
		n += Dt;									//add delta t to n
	}
}





/*!
 *  \brief showResultsInConsole function
 *
 *  Handle the choice of scheme and call the right one
 *
 *  \param type An integer which represent the scheme chosen
 *  \param dt Delta t choosen
 */


void showResultsInConsole(const int type, const double dt) {									//declare the function which select the correct scheme
	if (type == 1) {													//if the type is 1 =>
		cout << "----------- EXPLICIT UPWIND FTBS -------------------\n";						//Print a the scheme choosen
		resultDtExp(dt);												//call the function resultDt with delta t
	}	
	if (type == 2) {													//if the type is 2 =>
		cout << "----------- IMPLICIT UPWIND FTBS -------------------\n";						//Print a the scheme choosen
		//resultDtImpUpFTBS(dt);												//call the function resultDt with delta t
	}
	if (type == 3) {													//if the type is 4 =>
		cout << "----------- IMPLICIT FTCS -------------------\n";							//Print a the scheme choosen
		//resultDtImpFTCS(dt);												//call the function resultDt with delta t
	}
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
	setprecision(10);

	cout << "1 for EXPLICIT UPWIND FTBS \n2 for IMPLICIT UPWIND FTBS \n3 for IMPLICIT FTCS \n";	//print the text to explain how to acces to each scheme
	int type;											//declare a int type
	cin >> type;											//get the value of given by the user and store it in the int type
	cout << "1 for dt = 0.002\n2 for dt = 0.001\n3 for dt = 0.0005\n";				//print to choose between each dt
	int dtIndex;											//declare an int dtIndex
	cin >> dtIndex;											//get the value of given by the user and store it in the int dtIndex
	double dt;											//declare a double dt
	if (dtIndex == 1)										//if dtIndex equal 1 =>
		dt = 0.002;										//dt equal 0.002
	else if (dtIndex == 2)										//if dtIndex equal 2 =>
		dt = 0.001;										//dt equal 0.001
	else if (dtIndex == 3)										//if dtIndex equal 3 =>
		dt = 0.0005;										//dt equal 0.0005

	showResultsInConsole(type, dt);									//call the function showResultInConsole with the type of scheme and the delta t choosen
	return 0;											//if everything works correctly return 0
}

/*int main (int argc, char *argv []){
	int npes, myrank;
	int a[10], b[10];

	a[0] = 1;
	b[0] = 2;

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &npes);

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	//first one
	
	double start = MPI_Wtime();
	for(int i = 0; i < 10000; i++){
		for(int p = 0; p<npes; p++){
			if(p==myrank){
				MPI_Send(a, 10, MPI_INT, (p+1)%npes, 1, MPI_COMM_WORLD);
			} else if((p+1)%npes == myrank) {
				MPI_Recv(b, 10, MPI_INT, p, 1, MPI_COMM_WORLD, &status);
			}
		//printf("A value %d and b value %d \n", a[0], b[0]);
		}
	}
	double stop = MPI_Wtime();
	double elapse = (stop - start);
	printf("Time taken for an explicit scheme %f\n",elapse);
	

	MPI_Finalize();
	return 0;
}*/
