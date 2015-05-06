/* Exact Diagonalization */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include "matrix.h"
#include <cstdlib>

using namespace std;

typedef enum
{
	firstSite = 0, hubbard = 1
} Model;

typedef struct Parameters
{
    int N;// length of the chain
    int omegaPoints;// gridpoints for the frequency of A(omega) and G(omega)
    double u;// Coulomb interaction
    double eps;// on-site energy
    double t;// hopping strength
    double broadening;// broadening parameter for the delta-peaks
    double bandWidth;// range for the frequency
    Model model;// defines the model (either U_i = U*delta_{i,1} or U_i = U)
} Parameters;

typedef map<int, map<int, int> > QSzCount;// total number, Sz, number of state
typedef map<int, map<int, vector<double> > > Energies;// total number, Sz, energies
typedef map<int, map<int, Matrix> > States;// total number, Sz, eigenstates
typedef map<int, map<int, vector<vector<int> > > > Basis;
// total number, Sz, basis |-1,0,2,1,1.....> as occupation -1 spin dn,1 spin up, 0 empty, 2 double occupied

bool newConfiguration(vector<int> &s, int lower, int upper);// generate new configuration from lower(-1) spin dn to upper(2) double occupied
void broadeningPoles(const vector<double> &poles, const vector<double> &weights, vector<double> &newGrid, vector<double> &smoothFunction, Parameters &p);
void kramersKronig(const vector<double> &x, const vector<double> &fin, vector<double> &fout, Parameters &p);
double gaussian(double omega, double b);

int main(int argc, const char* argv[])
{
    int N = 6, charge, spin;
    map<int, int> countSubspaces;// the number of matrix with the same size
    Parameters parameters = {N, 5000, 4., -2., 0.5, 0.5, 20, hubbard/*firstSite*/};
    time_t start, end;
    time(&start);
    ofstream info("info.dat");
    Energies energies;
    States states;
    // s[i] labels a state at site i: s[i] = 0 means empty state,
    // -1 means spin down state, 1 means spin up state and 2 means doubly occupied state
    vector<int> s(parameters.N, -1);// initialize as all spin down
    Matrix hamiltonian;
    // qszcount counts the number of subspaces per quantum numbers Q, Sz
    QSzCount qszcount;
    Basis basis;
    
    // reads in command line arguments, it is not essential
    switch (argc)
    {
        case 7:
            parameters.model = (Model) atoi(argv[6]);
        case 6:
            parameters.broadening = atof(argv[5]);
        case 5:
            parameters.t = atof(argv[4]);
        case 4:
            parameters.eps = atof(argv[3]);
        case 3:
            parameters.u = atof(argv[2]);
        case 2:
            parameters.N = atoi(argv[1]);
    }

    int site = 0;//parameters.N/2-1;
    
    for (int q = -parameters.N; q <= parameters.N; q++)// the real charge should plus N 
        for (int sz = -parameters.N; sz <= parameters.N; sz++)
            qszcount[q][sz] = 0;
    
    info << parameters.N << "-site chain\n" << endl;
    if (parameters.model == firstSite)
        info << "tight-binding with interaction on first site only\n" << endl;
    else
        info << "Hubbard model\n" << endl;
    info << "U = " << parameters.u << "\neps = " << parameters.eps << "\nt = ";
    info << parameters.t << "\n" << endl;
    
    // Set up the basis
    //  setupBasis(parameters.N, parameters.N, basis, s, qszcount);
    bool hasNext = true;
    for (;;)
    {
        charge = 0;
        spin = 0;
        // loop over all sites for looking the Q and Sz, start from state
        for (unsigned int i = 0; i < s.size(); i++)
        {
            // Calculate spin and charge of each configuration...
            charge += abs(s[i])-1;
            spin += s[i]*(2-abs(s[i]));// the term *(2-abs(s[i])) is set for double occupied state, i.e. spin=0.
        }
        // and store the configuration in the corresponding Q, Sz subspace
        basis[charge][spin].push_back(s);
        qszcount[charge][spin]++; 
        // leave the loop if there exists no further basis state
        if (!hasNext)
            break;
        // determine next configuration and if there exists another one after that
        hasNext = newConfiguration(s, -1, 2);
    }
    
    // find and print largest subspace
    int largestSubspace = 0;
    for (int q = -parameters.N; q <= parameters.N; q++)
        for (int sz = -parameters.N; sz <= parameters.N; sz++)
        {
            if (qszcount[q][sz] > largestSubspace)
                largestSubspace = qszcount[q][sz];
            countSubspaces[qszcount[q][sz]]++;
        }
    
    // count number and sizes of subspaces
    for (map<int, int>::iterator it = countSubspaces.begin(); it != countSubspaces.end(); it++)
        info << (*it).second/*number of countSubspaces*/ << " subspaces of size " << (*it).first/*subsoace size (key qszcount)*/ << endl;
    info << "\nLargest subspace = " << largestSubspace << endl << endl;
    
    ofstream energiesOfStream("energies.dat");
    // loop through subspaces:
    for (int q = -parameters.N; q <= parameters.N; q++)
    {
        for (int sz = -parameters.N; sz <= parameters.N; sz++)
        {
            if (qszcount[q][sz] == 0)
                continue;
            
            // set up Hamiltonian in this subspace:
            hamiltonian.resize(qszcount[q][sz], qszcount[q][sz]);
            hamiltonian.zero();
            
            for (int r = 0; r < qszcount[q][sz]; r++) // ket state
            {
                if (parameters.model == firstSite)
                {
                    // Problem 2 a):
                    // =============
                    // Fill in the diagonal matrix elements of the Hamiltonian
                    // for the model with only first-site interaction
                    /*if ((basis[q][sz][r][r] == -1) )
                            hamiltonian.set(r, r, parameters.t);
                    if ((basis[q][sz][r][r] == -1) )
                            hamiltonian.set(r, r, parameters.t);
                    if ((basis[q][sz][r][r] == -1) )
                            hamiltonian.set(r, r, parameters.t);
                    if ((basis[q][sz][r][r] == -1) )
                            hamiltonian.set(r, r, parameters.t);*/
                }
                else
                {
                    // Problem 2 b):
                    // =============
                    // Fill in the diagonal matrix elements of the Hamiltonian
                    // for the Hubbard model
                    double temp=0; //temperory memory for diagonal element
                    for(int m=0; m<parameters.N; m++)
                    { 
                      //cout << basis[q][sz][r][m] << endl;
                      if ((basis[q][sz][r][m] == -1) )
                            temp += parameters.eps;
                      if ((basis[q][sz][r][m] ==  1) )
                            temp += parameters.eps;
                      if ((basis[q][sz][r][m] ==  0) )
                            temp =  temp;
                      if ((basis[q][sz][r][m] == 2) )
                            temp += 2*parameters.eps+parameters.u;

                      hamiltonian.set(r,r, temp);
                    }
                    //cout << "diag element=" << hamiltonian.get(r,r) <<endl;
                                        
                }
                
                // hopping between sites:
                for (int rp = 0; rp < qszcount[q][sz]; rp++) // bra state
                {
                    for (int m = 0; m < parameters.N-1; m++)// searching hoping term from the basis
                    {
                        bool p = false;
                        for (int mp = 0; mp < parameters.N; mp++)
                        {
                            // if anything but two neighbouring sites...
                            if ((mp == m) || (mp == m+1))
                                continue;
                            // ... are different from each other... 
                            if (basis[q][sz][r][mp] != basis[q][sz][rp][mp])
                                p = true;
                        }//mp
                        // ... then there couldn't be a non-vanishing matrix element in the Hamiltonian
                        if (p)
                            continue;
                        
                        // Problem 2 c):
                        // ==========
                        // In the following, fill in all the missing matrix elements
                        
                        if ((basis[q][sz][r][m] == 0) && (basis[q][sz][r][m+1] == 1) && (basis[q][sz][rp][m] == 1) && (basis[q][sz][rp][m+1] == 0))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == -1) && (basis[q][sz][r][m+1] == 1) && (basis[q][sz][rp][m] == 2) && (basis[q][sz][rp][m+1] == 0))
                            hamiltonian.set(r, rp, -parameters.t);
                        if ((basis[q][sz][r][m] == 0) && (basis[q][sz][r][m+1] == -1) && (basis[q][sz][rp][m] == -1) && (basis[q][sz][rp][m+1] == 0))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == 1) && (basis[q][sz][r][m+1] == -1) && (basis[q][sz][rp][m] == 2) && (basis[q][sz][rp][m+1] == 0))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == 0) && (basis[q][sz][r][m+1] == 2) && (basis[q][sz][rp][m] == 1) && (basis[q][sz][rp][m+1] == -1))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == -1) && (basis[q][sz][r][m+1] == 2) && (basis[q][sz][rp][m] == 2) && (basis[q][sz][rp][m+1] == -1))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == 0) && (basis[q][sz][r][m+1] == 2) &&(basis[q][sz][rp][m] == -1) && (basis[q][sz][rp][m+1] == 1))
                            hamiltonian.set(r, rp, -parameters.t);
                        if ((basis[q][sz][r][m] == 1) && (basis[q][sz][r][m+1] == 2) && (basis[q][sz][rp][m] == 2) && (basis[q][sz][rp][m+1] == 1))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == 1) && (basis[q][sz][r][m+1] == 0) &&(basis[q][sz][rp][m] == 0) && (basis[q][sz][rp][m+1] == 1))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == 1) && (basis[q][sz][r][m+1] == -1) &&(basis[q][sz][rp][m] == 0) && (basis[q][sz][rp][m+1] == 2))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == -1) && (basis[q][sz][r][m+1] == 0) &&(basis[q][sz][rp][m] == 0) && (basis[q][sz][rp][m+1] == -1))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == -1) && (basis[q][sz][r][m+1] == 1) &&(basis[q][sz][rp][m] == 0) && (basis[q][sz][rp][m+1] == 2))
                            hamiltonian.set(r, rp, -parameters.t);
                        if ((basis[q][sz][r][m] == 2) && (basis[q][sz][r][m+1] == 0) &&(basis[q][sz][rp][m] == -1) && (basis[q][sz][rp][m+1] == 1))
                            hamiltonian.set(r, rp, -parameters.t);
                        if ((basis[q][sz][r][m] == 2) && (basis[q][sz][r][m+1] == -1) && (basis[q][sz][rp][m] == -1) && (basis[q][sz][rp][m+1] == 2))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == 2) && (basis[q][sz][r][m+1] == 0) && (basis[q][sz][rp][m] == 1) && (basis[q][sz][rp][m+1] == -1))
                            hamiltonian.set(r, rp, parameters.t);
                        if ((basis[q][sz][r][m] == 2) && (basis[q][sz][r][m+1] == 1) && (basis[q][sz][rp][m] == 1) && (basis[q][sz][rp][m+1] == 2))
                            hamiltonian.set(r, rp, parameters.t);
                    }//m
                }//rp
            }//r
            
            // diagonalize the Hamiltonian in this subspace:
            cout << "Diagonalizing in subspace: Q = " << q << ", Sz = " << sz << ", subspace size: " << qszcount[q][sz] << endl;
//            hamiltonian.print();
            energies[q][sz] = hamiltonian.diag();
            states[q][sz] = hamiltonian;
            
            // write energies to file:
            for (int r = 0; r < qszcount[q][sz]; r++)
                energiesOfStream << q << "\t" << sz << "\t" << r << "\t" << energies[q][sz][r] << endl;

            //hamiltonian.print();            
            // deallocate hamiltonian
            hamiltonian.erase();
        }//sz
    }//q
    energiesOfStream.close();
    
    
    // Problem 3:
    // ==========
    // Determine the ground state
    bool flag = true;
    charge = 0, spin = 0;
    int rr = 0;
    double lowestEnergy = 0;

    // loop through subspaces:
    for (int q = -parameters.N; q <= parameters.N; q++)
    {
        for (int sz = -parameters.N; sz <= parameters.N; sz++)
        {
            for ( int r=0; r<energies[q][sz].size(); r++)
            {
                if(energies[q][sz][r] < lowestEnergy)
                {
                   lowestEnergy=energies[q][sz][r];
                   charge=q;
                   spin=sz;
                   rr=r;
                }
            }
        }
    }
    
    info << "Ground state quantum numbers and energy:" << endl;
    info << "Q = " << charge << ", Sz = " << spin  << ", r = " << rr << ", energy = " << lowestEnergy << endl << endl;
    info << "basis config= (" << "\t";
    for(int i=0; i<basis[charge][spin][rr].size(); i++)
        info << basis[charge][spin][rr][i] << "\t";
    info << ")" << endl; 
   
    // Problem 5:
    // ==========
    /* Calculation of the spectral weights:
     * Using the Lehmann representation one has to calculate the matrix elements for the
     * creation and annihilation operator (here we have chosen spin up) in the old basis
     * and then transform them to the new basis with the orthogonal transformation matrix
     * states. The weights are positioned at the corresonding eigenenergies.
     */
    // print the frequencies and the spectral weights in the file spectralWeights.dat
    ofstream spectralWeightOut("spectralWeights.dat");
    vector<double> weights, frequencies, exp_e;
    vector<int> temp_conf; //temp for configuration compared
    
    // Calculate weights and frequencies for the creation operator
    for (int i=0; i < energies[charge+1][spin+1].size() ; i++) // i -> loop through the right eigenstates
    {
        double temp_weight=0; // temp memory for weight
        for ( int j=0; j< energies[charge][spin].size() ; j++ ) // j -> loop through the left eigenstates's basis
        {   // compute the basis <0,0,0|c1_up to temp_conf
            if( basis[charge][spin][j][site] == 1 || basis[charge][spin][j][site] ==2 )
              continue;
            if( basis[charge][spin][j][site] == -1)
            {
              temp_conf = basis[charge][spin][j];
              temp_conf[site] = 2;
              //for(int k=0; k<temp_conf.size() ;k++)
              //   cout<< temp_conf[k];
              //cout<<endl;
              for ( int k=0; k< energies[charge+1][spin+1].size() ; k++) // k -> loop through the right eigenstates's basis
              {
                  if( temp_conf == basis[charge+1][spin+1][k] )// find complementary |1,1,k>
                      temp_weight +=  states[charge][spin].get(rr,j)*states[charge+1][spin+1].get(i,k);
              }
            }
            if( basis[charge][spin][j][site] == 0 )
            {
              temp_conf = basis[charge][spin][j];
              temp_conf[site] = 1;
              //for(int k=0; k<temp_conf.size() ;k++)
              //   cout<< temp_conf[k];
              //cout<<endl;
              for ( int k=0; k< energies[charge+1][spin+1].size() ; k++) // k -> loop through the right eigenstates's basis
              {
                  if( temp_conf == basis[charge+1][spin+1][k] )// fin complementary |1,1,r>
                      temp_weight +=  states[charge][spin].get(rr,j)*states[charge+1][spin+1].get(i,k);
              }
            }
        }
        weights.push_back( pow(temp_weight, 2 ) );
        frequencies.push_back(  energies[charge+1][spin+1][i] -energies[charge][spin][j] );
        exp_e.push_back( exp(-beta*energies[charge+1][spin+1][i]) + exp(-beta*energies[charge][spin][j]) )
    }    
    // Calculate weights and frequencies for the annihilation operator


    for (int i=0; i < energies[charge-1][spin-1].size() ; i++) // i -> loop through the right eigenstates
    {
        double temp_weight=0; // temp memory for weight
        for ( int j=0; j< energies[charge][spin].size() ; j++ ) // j -> loop through the left eigenstates's basis
        {   // compute the basis <0,0,0|c1_up to temp_conf
            if( basis[charge][spin][j][site] == 0 || basis[charge][spin][j][site] ==-1 )
              continue;
            if( basis[charge][spin][j][site] == 1)
            {
              temp_conf = basis[charge][spin][j];
              temp_conf[site] = 0;
              //for(int k=0; k<temp_conf.size() ;k++)
              //   cout<< temp_conf[k];
              //cout<<endl;
              for ( int k=0; k< energies[charge-1][spin-1].size() ; k++) // k -> loop through the right eigenstates's basis
              {
                  if( temp_conf == basis[charge-1][spin-1][k] )// find complementary |1,1,k>
                      temp_weight +=  states[charge][spin].get(rr,j)*states[charge-1][spin-1].get(i,k);
              }
            }
            if( basis[charge][spin][j][site] == 2 )
            {
              temp_conf = basis[charge][spin][j];
              temp_conf[site] = -1;
              //for(int k=0; k<temp_conf.size() ;k++)
              //   cout<< temp_conf[k];
              //cout<<endl;
              for ( int k=0; k< energies[charge-1][spin-1].size() ; k++) // k -> loop through the right eigenstates's basis
              {
                  if( temp_conf == basis[charge-1][spin-1][k] )// fin complementary |1,1,r>
                      temp_weight +=  states[charge][spin].get(rr,j)*states[charge-1][charge-1].get(i,k);
              }
            }
        }
        weights.push_back( pow(temp_weight, 2 ) );
        frequencies.push_back(  energies[charge][spin][rr] -energies[charge-1][spin-1][i] );
    }

    for (unsigned int i = 0; i < frequencies.size(); i++)
        spectralWeightOut << frequencies[i] << "\t" << weights[i] << endl;
    spectralWeightOut.close();
    

	
    
    // Check if the spectral weight is really 1
    double weight = 0;
    for (unsigned int i = 0; i < frequencies.size(); i++)
        weight += weights[i];
    info << "Spectral weight = " << weight << endl;
    
    // broadening of the spectral function weights
    vector<double> omega, specFunc;
    broadeningPoles(frequencies, weights, omega, specFunc, parameters);
    ofstream spectralFunctionOut("spectralFunction.dat");
    for (unsigned int i = 0; i < omega.size(); i++)
        spectralFunctionOut << omega[i] << "\t" << specFunc[i] << endl;
    spectralFunctionOut.close();
    
    // calculation of the Greensfunction and its real part via Kramers Kronig
    vector<double> imGreensFct(parameters.omegaPoints), reGreensFct(parameters.omegaPoints);
    for (int i = 0; i < parameters.omegaPoints; i++)
        imGreensFct[i] = -specFunc[i]*M_PI;
    kramersKronig(omega, imGreensFct, reGreensFct, parameters);
    ofstream printGF("greenFunction.dat");
    for (int i = 0; i < parameters.omegaPoints; i++)
        printGF << omega[i] << "\t" << imGreensFct[i] << "\t" << reGreensFct[i] << endl;
    printGF.close();
    
    int numberOfStates = 0;
    for (int q = -parameters.N; q <= parameters.N; q++)
        for (int sz = -parameters.N; sz <= parameters.N; sz++)
            numberOfStates += qszcount[q][sz];
    info << numberOfStates << "\t" << "states" << endl;
    time(&end);
    info << "program execution needed " << difftime(end, start) << " seconds" << endl;
    
    return 0;
}

bool newConfiguration(vector<int> &s, int lower, int upper)
{
    for (unsigned int i = 0; i < s.size(); i++)
    {
        if (s[i] < upper)
        {
            // increase one state, then leave the loop
            s[i]++;
            break;
        } else
            s[i] = lower;
    }
    // if there is any state not doubly occupied, we have some more
    // states to build and return true, ... 
    for (unsigned int i = 0; i < s.size(); i++)
        if (s[i] != upper)
            return true;
    // ... else we return false
    return false;
}

void broadeningPoles(const vector<double> &poles, const vector<double> &weights, vector<double> &newGrid, vector<double> &smoothFunction, Parameters &p)
{
    double stepWidth = 2*p.bandWidth*p.t/p.omegaPoints;
    // set up the grid for the frequency values
    newGrid.resize(p.omegaPoints);
    smoothFunction.resize(p.omegaPoints);
    for (int i = 0; i < p.omegaPoints; i++)
        newGrid[i] = -p.bandWidth*p.t + i*stepWidth;
    
    // summation of contributions of all gaussians times their weight at each frequency
    for (int i = 0; i < p.omegaPoints; i++)
    {
        smoothFunction[i] = 0;
        for (unsigned int j = 0; j < poles.size(); j++)
        {
            smoothFunction[i] += weights[j] * gaussian(newGrid[i]-poles[j], p.broadening);
        }
    }
}

double gaussian(double omega, double b)
{
    b = 1/b;
    return b*exp(-omega*omega*b*b)/sqrt(M_PI);
}


void kramersKronig(const vector<double> &x, const vector<double> &fin, vector<double> &fout, Parameters &p)
{
    /* This function calculates the Kramers-Kronig-Transform of the function fin
     * using the trapezian integration method and returns the result in function
     * fout.
     * Kramers Kronig int_{-/infty}^/infty fin(y)/(x-y)dy
     * In order not to divide by zero -> i != j.
     */
    
    double deltax = x[1] - x[0];
    for (int i = 0; i < p.omegaPoints; i++)
    {
        if (i != 0 && i != p.omegaPoints-1)
            fout[i] = 0.5*(fin[0]/(x[0]-x[i]) + fin[p.omegaPoints-1]/(x[p.omegaPoints-1]-x[i]));
        for (int j = 1; j < p.omegaPoints-1; j++)
        {
            if (i != j)
                fout[i] += fin[j]/(x[j]-x[i]);
        }
        fout[i] *= -deltax/M_PI;
    }
}




