//test program
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "point.h"

using namespace std;

int main(int argc, char* argv[])
{	
	//verbose: true for more output; false for less output
	const bool verbose = true;
	
	//check for name of data file
	if(argc == 1)
	{
		cout << "USAGE: run <filename>\n";
		return 1;
	}
	
	//integer dimension of data
	int d;	
	
	//create vector for points
	vector<Point> points;
		
	//read points from file
	if(verbose) { cout << "READING FILE:\n"; }
	string line;
	ifstream myfile(argv[1]);
	if(myfile.is_open())
	{
		//get dimension of the points from the first line of the file
		getline(myfile,line);
		stringstream(line) >> d;
		if(verbose) { cout << "  dimension: " << d << "\n"; }
		
		while( getline(myfile,line) )
		{
			//cout << "point: " << line << '\n'; //TESTING
			
			//parse current point from string
			istringstream iss(line);
			double* n = new double[d];
			//cout << "address of n[]: " << n << " " << &n << "\n"; //TESTING
			for(int i=0; i<d; i++)
			{
				iss >> n[i];	//extract the next double from the string
			}
			double t;	//time of birth for this point
			iss >> t; 
			
			Point p (n, t);
			if(verbose) 
			{
				double *m = p.get_coords();
				cout << "  point: (";
				for(int i=0; i<d; i++)
				{
					cout << m[i];
					if(i<d-1) { cout << ", "; }
				}
				cout << ") born at time " << p.get_birth() << "\n";
			}
			
			//add current point to the vector
			points.push_back(p);
		}
		myfile.close();
		
		cout << "Input finished.\n";
	}
	else
	{
		cout << "Error: Unable to open file " << argv[1] << ".\n";
		return 1;
	}
	
	//test points vector
	cout << "TESTING VECTOR:\n";
	for(int i=0; i<points.size(); i++)
	{
		Point p = points.at(i);
		double *m = p.get_coords();
		cout << " point: (";
		for(int i=0; i<d; i++)
		{
			cout << m[i];
			if(i<d-1) { cout << ", "; }
		}
		cout << ") born at time " << p.get_birth() << "\n";		
	}
	
	
	
	
	
	
	
	cout << "\n";
	
	
}
