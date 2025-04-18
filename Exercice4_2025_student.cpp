#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ConfigFile.tpp"


using namespace std;


const double epsilon_0 = 8.85418782e-12 ; 

const double PI=3.1415926535897932384626433832795028841971e0;
// Résolution d'un système d'équations linéaires par élimination de
// Gauss-Jordan:
template<class T>
vector<T>
solve(const vector<T>& diag,
      const vector<T>& lower,
      const vector<T>& upper,
      const vector<T>& rhs)
{
    vector<T> solution(diag.size());
    vector<T> new_diag(diag);
    vector<T> new_rhs(rhs);

    for (int i = 1; i < diag.size(); ++i) {
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution[diag.size() - 1] =
      new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i(diag.size() - 2); i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    return solution;
}

//TODO build the epsilon function
double epsilon(double r, double r1, double eps_a , double eps_b ) {
   
    if ( r < r1 )
    { return eps_a ; }
	else 
	{ return eps_b ; }
}


//TODO build the rho_epsilon function (rho_lib / epsilon_0)
double rho_epsilon(double r , double r1 , bool unif_rho , double rho_lib ,  double eps_0 = epsilon_0 )
{
	if ( unif_rho ) 
	{ return rho_lib / eps_0 ; }
	else // cas rho non-uniforme ; attention rho-lib = 10e4 dans ce cas ( changer dans le configfile ) 
	{ 
		if ( r < r1 )
		{ return rho_lib * sin(PI*r/r1) ; }
		else 
		{ return 0.0 ; }
		
	}
}

int
main(int argc, char* argv[])
{    

    // USAGE: Exercise4 [configuration-file] [<settings-to-overwrite> ...]

    // Read the default input
    string inputPath = "/Users/a-x-3/Desktop/Exercice4_2025_student/configuration.in.example";
    // Optionally override configuration file.
    if (argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    // Override settings
    for (int i = 2; i < argc; i++)
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Read geometrical inputs
    const double R  = configFile.get<double>("R");
    const double r1 = configFile.get<double>("r1"); 

    const double rho0 = configFile.get<double>("rho0");

    // For the analytical comparison
    const bool uniform_rho_case = configFile.get<bool>("uniform_rho_case");
    
    // Dielectric relative permittivity
    const double epsilon_a = configFile.get<double>("epsilon_a");
    const double epsilon_b = configFile.get<double>("epsilon_b");
    
    // Boundary conditions
    const double VR = configFile.get<double>("VR");
    
    // Discretization
    //const int alpha = configFile.get<int>("alpha") ; // coeff de proportionnalité 
    const int N1 = configFile.get<int>("N1");
    const int N2 = configFile.get<int>("N2"); //alpha * N1 ; // 
    cout << N1 << endl ; 
    cout << N2 << endl ; 
    
    // Fichiers de sortie:
    string fichier = configFile.get<string>("output");
    string fichier_phi = fichier + "_phi.out";  //"/Users/a-x-3/Desktop/Exercice4_2025_student/" + fichier + "_phi.out"; // /Users/a-x-3/Desktop/Exercice4_2025_student/
    string fichier_E   = fichier + "_E.out";   //"/Users/a-x-3/Desktop/Exercice4_2025_student/" + fichier + "_E.out";
    string fichier_D   = fichier + "_D.out";  //"/Users/a-x-3/Desktop/Exercice4_2025_student/" + fichier + "_D.out";

    // Create our finite elements
    const int pointCount = N1 + N2 + 1; // Number of finite elements
    const double h1 = (r1 - 0) / N1;
    const double h2 = (R  - r1) / N2;

    // Position of elements
    vector<double> r(pointCount);

    //TODO build the nodes vector r
    
    // Arrays initialization
    vector<double> h(pointCount-1); 	// Distance between grid points
    vector<double> midPoint(pointCount-1);  // Midpoint of each grid element
        
    // cout << epsilon_a << endl ; 	    
        
    // TODO build the h vector and midpoint vector
            
	for ( size_t i = 0 ; i < N1 ; ++i  ) // construction de h (première moitié)
	{ h[i] = h1 ; }
		
	for ( size_t i = N1 ; i < pointCount - 1 ; ++i  ) // construction de h (deuxième moitié)
	{ h[i] = h2 ; }
	
	for (size_t i = 0 ; i < pointCount ; ++i) // construction du vecteur position
	{
		if ( i == 0 )
		{ r[i] = 0 ; }
		else
		{ r[i] = r[i-1] + h[i-1] ; }
	}	
	
	for ( size_t i = 0 ; i < pointCount - 1 ; ++i  ) // construction du mid point 
	{ 
		if ( i == 0 )
		{ midPoint[i] = h1/2 ; }
		else 
		{
			if (i < N1)
			{ midPoint[i] = midPoint[i-1] + h1 ; }
			else 
			{
				if ( i == N1 )
				{ midPoint[i] = ( midPoint[i-1] + h1/2 + h2 + r[i] ) / 2 ; } // on fait le milieu entre les deux points des différents maillages
				else 
				{ 
					if ( i == N1 + 1 )
					{ midPoint[i] = r[i] + h2/2 ; } // on se replace au bon endroit en prenant le mileu entre les deux premiers points du nouveau maillage 
					else
					{ midPoint[i] = midPoint[i-1] + h2 ; }
				}
			}
		}
	 }
	
	//for ( auto const & ri : r  )
	//{ cout << ri << endl ; }
	
	//cout << 'h' << endl ; 
	
	//for ( auto const & hi : h )
	//{ cout << hi << endl ; }
	
	//for ( size_t i = 0 ; i < pointCount - 1  ; ++i  )
	//{ cout << i << ' ' << midPoint[i] << endl ; }
	
	//cout << "rho_eps : " << rho_epsilon(2) << endl ; 
    
    
    // N1 intervalles équidistants entre r = 0 et r = r1
    // N2 intervalles équidistants entre r = r1 et r = R 
    
    // Construct the matrix and right-hand side
    vector<double> diagonal(pointCount, 1.0);  // Diagonal
    vector<double> lower(pointCount - 1, 0.0); // Lower diagonal
    vector<double> upper(pointCount - 1, 0.0); // Upper diagonal
    vector<double> rhs(pointCount, 0.0);       // Right-hand-side
    
    // Loop over the intervals: add the contributions to matrix and rhs   
	for (int k = 0; k < pointCount-1; ++k) {
        // TODO build the vectors diagonal, lower, upper, rhs
        
        if ( k == 0 )
        {
			diagonal[k] = midPoint[k] * epsilon(midPoint[k],r1,epsilon_a,epsilon_b) / ( h[k]) ; // pas de k-1 => intégrale gauche nulle 
			//lower[k]    = 0 ; // pas de k-1 => intégrale nulle 
			upper[k]    = - midPoint[k] * epsilon(midPoint[k],r1,epsilon_a,epsilon_b) / ( h[k]) ; 
			rhs[k] 		= h[k] * midPoint[k] * rho_epsilon (midPoint[k],r1,uniform_rho_case,rho0) / 2 ; // pas de k-1 => intégrale gauc
		}
		else 
		{
			//cout << "k-1 : " << epsilon(midPoint[k-1],r1,epsilon_a,epsilon_b) << endl ; 
			//cout << "k : " << epsilon(midPoint[k],r1,epsilon_a,epsilon_b) << endl ; 
			diagonal[k] = midPoint[k-1] * epsilon(midPoint[k-1],r1,epsilon_a,epsilon_b) / ( h[k-1]) + midPoint[k] * epsilon(midPoint[k],r1,epsilon_a,epsilon_b) / ( h[k]) ; 
			lower[k-1]    = - midPoint[k-1] * epsilon(midPoint[k-1],r1,epsilon_a,epsilon_b) / ( h[k-1]) ; 
			upper[k]    = - midPoint[k] * epsilon(midPoint[k],r1,epsilon_a,epsilon_b) / ( h[k]) ; 
			rhs[k] 		= h[k-1] * midPoint[k-1] * rho_epsilon (midPoint[k-1],r1,uniform_rho_case,rho0) / 2 + h[k] * midPoint[k] * rho_epsilon (midPoint[k],r1,uniform_rho_case,rho0) / 2   ; 
		}

    }
        
    // TODO boundary condition at r=R (modify the lines below)
    lower[lower.size() - 1]       = 0.0;
    diagonal[diagonal.size() - 1] = 1.0;
    rhs[rhs.size() - 1] = VR;

	//cout << "diagonal" << endl ; 
	//for ( auto const & dia : diagonal )
	//{ cout << dia << endl ; } 
	
	//cout << "lower" << endl ; 
	//for ( auto const & low : lower )
	//{ cout << low << endl ; } 

	//cout << "upper" << endl ; 
	//for ( auto const & up : upper )
	//{ cout << up << endl ; } 
	
    // Solve the system of equations
    vector<double> phi = solve(diagonal, lower, upper, rhs);

    // Calculate electric field E and displacement vector D
    vector<double> E(pointCount - 1, 0);
    vector<double> D(pointCount - 1, 0);
    vector<double> alphadr(pointCount -2 , 0); // dérivée de alpha par rapport à r
    for (int i = 0; i < E.size(); ++i) {
        // TODO calculate E and D
        E[i] = (phi[i] - phi[i+1]) / h[i] ; 
        // cout << E[i] << endl ; 
        D[i] = epsilon_0 * epsilon(midPoint[i], r1 , epsilon_a , epsilon_b ) * E[i] ; 
        //cout << midPoint[i] << " : " << epsilon(midPoint[i], r1 , epsilon_a , epsilon_b ) << endl ; 
    }
    
    //for (int i = 0; i < alpha.size(); ++i) 
    //{
		//alphadr[i] = E[i] * epsilon; 
	//}

    // Export data
    {
        // Electric potential phi
        ofstream ofs(fichier_phi);
        ofs.precision(15);

        if (r.size() != phi.size())
            throw std::runtime_error("error when writing potential: r and "
                                     "phi does not have size");

        for (int i = 0; i < phi.size(); ++i) {
            ofs << r[i] << " " << phi[i] << " " << rhs[i] << endl;
            cout << rhs[i] << endl ; 
        }
    }

    {
        // Electric field E
        ofstream ofs(fichier_E);
        ofs.precision(15);

        if (r.size() != (E.size() + 1))
            throw std::runtime_error("error when writing electric field: size of "
                                     "E should be 1 less than r");

        for (int i = 0; i < E.size(); ++i) {
            ofs << midPoint[i] << " " << E[i] << endl;
        }
    }
    {
        // Displacement field D
        ofstream ofs(fichier_D);
        ofs.precision(15);

        if (E.size() != D.size())
            throw std::runtime_error("error when writing displacement field: size of "
                                     "D should be equal to E");

        for (int i = 0; i < D.size(); ++i) {
            ofs << midPoint[i] << " " << D[i] << endl;
        }
    }

    return 0;
}

