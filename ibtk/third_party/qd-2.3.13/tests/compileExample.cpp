

// simple exaple of QD usage to illustrate linking process
// Alex Kaiser, LBNL, 6/3/2010



#include <iostream>
#include <qd/qd_real.h>
#include <qd/fpu.h>

using namespace std; 

int main() {
	
	// ensure that 80-bit arithmetic is not in place
	// this call forces 64-bit arithmetic
	unsigned int old_cw;
	fpu_fix_start(&old_cw);

	cout.precision(60); 
	
	// simple read example
	/*
	qd_real readTest ; 
	cin >> readTest ; 	
	cout << "readTest = " << readTest << endl ;
	 */ 
	
	// simple demo
	qd_real x = "1.0" ;  
	x /= 3.0 ;
	qd_real y ; 
	y = pow( qd_real(2.0) , 3 ) ; 
	cout << "y = " << y << endl;
	cout << "x = " << x << endl; 
	
	
	qd_real a ; 
	qd_real b = qd_real("0.1");
	
	a = sqrt(b);
	cout << " sqrt(0.1) = " << a << endl;
	cout << " sqrt(0.1) * sqrt(0.1) = " << a * a << endl; 
	
	fpu_fix_end(&old_cw); 
	return 0;
}


