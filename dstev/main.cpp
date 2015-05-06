#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <complex>
#include <cstdlib>
#include <vector>
#include "matrix.h"
#include "diag.h"

using namespace std;

int main(){

//Matrix a(1,2); //( x, x)
//Matrix b(2,1); //( x )
               //( x )

Matrix a(3,3);
vector<double > eigenvalues;

for(int i=0; i<3; i++){
   for(int j=0; j<3; j++){
      if (i==j) a.set(i,j,1);
      else if(i==j+1 || i==j-1) a.set(i,j,1);
      else a.set(i,j,0);
   }
}

cout << "matrix is:" << endl;
a.print();
cout << "eigenvalues from full diagonalization are:" << endl;
eigenvalues=a.diag();
for(int i=0; i<eigenvalues.size(); i++) cout << eigenvalues[i] << "\t";
cout << endl;

cout << "eigenvalues from tridiagonalization are:" << endl;

vector<double> diag(3,1);
vector<double> offdiag(2,1);

int N = diag.size();
//vector<double> eigenvalues(N);
double eigenvectors[N*N];
char v = 'V';
int info, lwork = max(1, 3*N-1);
double *work = NULL;

work = new double[lwork];
utils::dstev_(&v, &N, &diag.at(0), &offdiag.at(0), eigenvectors , &N ,work, &info);

if (info != 0)
  cout << "diagonalization routine dsyev_ error: info = " << info << endl;Matrix c;

for(int i=0; i<diag.size(); i++) cout << diag[i] << "\t";
cout << endl;

cout << "The eigenvectors are:" << endl;
for(int i=0; i<N*N;i++) cout << eigenvectors[i] <<"\t";
cout<< endl;

/*Matrix d;

a.set(0,0,1);
a.set(0,1,2);

b.set(0,0,2);
b.set(1,0,1);


//cout<<"yo"<<endl;
c=a*b;
d=b*a;

cout<<"c=a*b"<<endl;
c.print();
cout<<"d=b*a"<<endl;
d.print();*/

return 0;
}

