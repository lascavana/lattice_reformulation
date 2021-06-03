#include <stdio.h>

// basic file operations
#include <iostream>
#include <fstream>
// getrusage
#include <sys/time.h>
#include <sys/resource.h>
// NTL library
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/LLL.h>

using namespace std;
using namespace NTL;

// global variables
int n, m;
mat_ZZ Aext;
long c[1000], u[1000];
ZZ N1=to_ZZ(10000000);
ZZ N2=to_ZZ(1000000000);


void ReadConstraints(const char *finput)
{
  int i, j;

  // open file for reading
  ifstream input_file(finput);

  // read A(ext), c, u
  if (input_file.is_open())
  {
    input_file >> n >> m;

    // read extended matrix
    Aext.SetDims(m,n+1);
    for (i = 0; i < m; i++)
  	{
  	    for (j = 0; j < (n); j++)
          input_file >> Aext[i][j];
        input_file >> Aext[i][n];
    }
    // read objective function
    for (j=0; j<n; j++) input_file >> c[j];
    // read upper bounds
    for (j=0; j<n; j++) input_file >> u[j];

    input_file.close();
  }
  else cout << "Unable to open file";
}


void check_basis(mat_ZZ Q)
{
  int i, j;
  mat_ZZ Atrans, X;

  // create transposed and non-extended version of matrix
  Atrans.SetDims(n,m);
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++) Atrans[j][i] = Aext[i][j];
    }

  mul(X, Q, Atrans);
  if (IsZero(X)) cout << "  Basis was correctly generated." << endl;
  else
  {
    cout << "Error in reduction." << endl;
    cout << "Q*A^T is not zero." << endl;
    throw  std::bad_function_call();
  }
}


void MakeAHLLattice(mat_ZZ &Q, vec_ZZ &x0)
{
  int i,j;
  mat_ZZ L, U, X1, Atrans;
  ZZ determ;

  // create L matrix
  L.SetDims(n+1,n+m+1);
  for (i=0;i<n+1;i++)
  {
     for(j=0;j<n+m+1;j++)
         L[i][j]=0;
     L[i][i]=to_ZZ(1);
  }
  L[n][n]= N1;
  for (j=0;j<n;j++){
     for (i=0; i<m; i++)
        L[j][n+i+1] = to_ZZ(Aext[i][j])*N2;
  }
  for (i=0; i<m; i++)
    L[n][n+i+1] = to_ZZ(Aext[i][n])*(N2);

  // lattice reduction
  LLL(determ, L, U, 99, 100, 0);

  // check for successful reduction
  if (L[n-m][n] != N1 && L[n-m][n] != -N1)
  {
    cout << "N1 does not appear in position [n-m, n]" << endl;
    cout << "L[n-m][n] = " << L[n-m][n] << endl;
    cout << "N1 = " << N1 << endl;
    throw  std::bad_function_call();
  }

  // get AHL basis
  Q.SetDims(n-m,n);
  for (i=0;i<n-m; i++){
    for (j=0;j<n;j++)
      Q[i][j]=L[i][j];
   }
   check_basis(Q);

   x0.SetLength(n);
   for (i=0; i<n; i++)
      x0[i] = L[n-m][i];
}


void MakeCEALattice(mat_ZZ &Q, vec_ZZ &x0)
{
  int i,j;
  mat_ZZ U, L, X1;
  ZZ determ;

  // create L matrix
  L.SetDims(n+1,n+m+1);
  for (i=0;i<n;i++)
    L[i][i]=to_ZZ(2);
  L[n][n]= to_ZZ(N1);
  for (j=0;j<n+1;j++){
    for (i=0; i<m; i++)
      L[j][n+i+1] = to_ZZ(2*Aext[i][j])*to_ZZ(N2);
  }
  for (j=0; j<n; j++)
    L[n][j] = to_ZZ(-1);

  // lattice reduction
  LLL(determ, L, U, 99, 100, 0);

  // get CEA basis
  Q.SetDims(n-m,n);
  for (i=0; i<n-m; i++){
    for (j=0; j<n; j++)
      Q[i][j]=L[i][j]/to_ZZ(2);
  }
  check_basis(Q);

  x0.SetLength(n);
  for (j=0; j<n; j++){
      x0[j] = to_ZZ(0.5*to_RR(L[n-m][j]) + 0.5);
  }
}


vec_ZZ make_objective(mat_ZZ basis, vec_ZZ x0)
{
  Vec<ZZ> obj;
  int i,j;

  obj.SetLength(n-m+1);

  for (i=0; i<(n-m); i++)
  {
      obj[i] = 0;
      for (j=0; j<n; j++)
        obj[i] += to_ZZ(c[j])*basis[i][j];
  }  

  obj[n-m]=0;
  for (j=0; j<n; j++)
    obj[n-m]+=to_ZZ(c[j])*x0[j];

  return obj;
}


void print_original(const char *filename)
{
  int i, j, k = 0;
  ofstream output_file(filename);

  // objective function
  output_file << "min ";
  for (j=0; j<(n-1); j++)
  {
    if (k>9){k=0; output_file << "\n";}
    output_file << c[j] << "x" << j+1 << " +";
    k++;
  }
  output_file << c[n-1] << "x" << n << "\n";

  // constraints
  output_file << "subject to" << "\n";
  for (j=0;j<m;j++)
  {
      k=0;
      output_file << Aext[j][0] << "x1 ";
      for (i=1;i<n;i++)
      {
        if (Aext[j][i] != 0)
        {
          if (k>9){k=0; output_file << "\n";}
          if (Aext[j][i] >= 0) output_file << "+";
          output_file << Aext[j][i];
          output_file << "x" << i+1 << " ";
          k++;
        }
       }
      output_file << "= ";
      output_file << (Aext[j][n])*to_ZZ(-1) << "\n";
  }

  // bounds
  output_file << "Bounds " << "\n";
  for (i=0; i<n; i++)
      output_file << "x" << i+1 << " <= " << u[i] << "\n";

  // variable definitions
  output_file << "Generals " << "\n";
  k = 0;
  for (i=0; i<n; i++)
  {
      if (k>15){k=0; output_file << "\n";}
      output_file << "x" << i+1 << " ";
      k++ ;
  }

  output_file << "\nEnd \n";
  output_file << "\n\n";
  output_file.close();
}


void print_reformulation(const char *filename, mat_ZZ basis, vec_ZZ x0)
{
  int i, j, k = 0, poscounter = 0;
  ofstream output_file(filename);

  // objective function
  vec_ZZ obj = make_objective(basis, x0);
  output_file << "min ";
  for (j=0; j<(n-m); j++)
  {
    if (obj[j] != 0)
    {
      if (k>9){k=0; output_file << "\n";}
      if (poscounter > 0 && obj[j]>0) output_file << " + ";
      output_file << obj[j] << "k" << j+1;
      poscounter++;
      k++;
    }
  }
  if (obj[n-m] > 0) output_file << "+ " << obj[n-m] << "\n";
  if (obj[n-m] < 0) output_file << obj[n-m] << "\n";


  // constraints
  output_file << "subject to" << "\n";
  for (j=0;j<n;j++)
  {
    k=0; poscounter = 0;
    for (i=0; i<n-m; i++)
    {
      if (basis[i][j] != 0)
      {
        if (k>9){k=0; output_file << "\n";}
        if (basis[i][j] >= 0 && poscounter > 0) output_file << "+";
        output_file << basis[i][j] << "k" << i+1 << " ";
        poscounter++;
        k++;
      }
    }
    output_file << " >= ";
    output_file << (x0[j])*(to_ZZ(-1)) << "\n";
  }
  for (j=0;j<n;j++)
  {
    k=0; poscounter = 0;
    for (i=0; i<n-m; i++)
    {
      if (basis[i][j] != 0)
      {
        if (k>9){k=0; output_file << "\n";}
        if (basis[i][j] >= 0 && poscounter > 0) output_file << "+";
        output_file << basis[i][j] << "k" << i+1 << " ";
        poscounter++;
        k++;
      }
    }
    output_file << " <= ";
    output_file << to_ZZ(u[j]) - x0[j] << "\n";
  }


  // bounds
  output_file << "Bounds " << "\n";
  for (i=0; i<(n-m); i++)
    output_file << "-inf <= k" << i+1 << " <= +inf \n";

  // variable definitions
  output_file << "Generals " << "\n";
  k=0;
  for (i=0; i<(n-m); i++)
  {
    if (k>15){k=0; output_file << "\n";}
    output_file << "k" << i+1 << " ";
    k++ ;
  }


  output_file << "\nEnd \n";
  output_file << "\n\n";
  output_file.close();
}



void print_translation(mat_ZZ basis, vec_ZZ x0)

{
  int i,j,k;
  ZZ determ;
  mat_ZZ L, L1, basisT, test_matrix; 
  mat_ZZ translBasis;


  translBasis.SetDims(n-m,n);

  transpose(basisT, basis);

  
  L.SetDims(n+1,2*n-m+1);
  for (i=0;i<n+1;i++){
     for(j=0;j<2*n-m+1;j++)
         L[i][j]=0;
  }
  for (i=0;i<n;i++)
     L[i][i]=to_ZZ(1);

  L[n][n]= N1;
  for (j=0;j<n;j++){
     for (i=0; i<n-m; i++){
        L[j][n+i+1] = (basis[i][j])*N2;
     }
  }

  L1.SetDims(n+1,2*n-m+1);
  cout << "  Getting translation matrix W" << endl;
  for (k=0; k<(n-m); k++)
  {
    cout << "  Step " << k << endl;
      for (i=0; i<n+1; i++){
          for (j=0; j<2*n-m+1; j++)
        L1[i][j] = L[i][j];
      }
      for (j=n+1; j<2*n-m+1; j++)
          L1[n][j]=to_ZZ(0);
      L1[n][n+k+1] = to_ZZ(-1)*(N2);
      
      LLL(determ, L1, 99, 100, 0);

      if (L1[m][n] != N1 && L1[m][n] != -N1)
      {
         cout << "N1 does not appear in position [m, n]" << endl;
         cout << L1[m][n] << endl;
         cout << "N1 = " << N1 << endl;
         throw  std::bad_function_call();
      }


      for (j=0; j<n; j++)
         translBasis[k][j] = L1[m][j];
     
  }

 
   mul(test_matrix,translBasis,basisT);
   cout << "  The translated matrix * transpose of the basis is the following matrix" << endl;
   cout << test_matrix << endl;
   long flag = IsIdent(test_matrix,n-m);
   if (flag != 1){
      cout << "  Translation is not correct" << endl;
      throw  std::bad_function_call();
   }

  ofstream output_file1("translation_W.txt");
  for (i=0;i<n-m;i++)
  {
    for (j=0;j<n;j++)
      output_file1 << translBasis[i][j] << " ";
    output_file1 << "\n";
  }
  output_file1.close();

  ofstream output_file2("translation_x0.txt");
  for (i=0;i<n;i++)
  {
    output_file2 << x0[i] << "\n";
  }
  output_file2.close();

}


int main(int argc, char *argv[])
{
  long ttime;
  string opt;

  // parse command line parameters
  if (argc < 3)
  {
    cout << "Input as follows: reduce <inputfile> <outputfile>\n";
    return 0;
  }
  if (argc == 4)
  {
    opt = argv[3];
    if (opt != "--translate")
    {
      cout << "Unidentified command line argument " << argv[3] << endl;
      cout << "Valid arguments: " << endl;
      cout << "--translate :: prints translation matrix " << endl;
      throw std::invalid_argument("Invalid argument");
    }
  }

  // read input file
  ReadConstraints(argv[1]);


  // track resource usage
  struct rusage tusage;
  getrusage(RUSAGE_SELF, &tusage);
  ttime = tusage.ru_utime.tv_usec + tusage.ru_stime.tv_usec;


  // get both reformulations
  mat_ZZ Q_ahl, Q_cea;
  vec_ZZ x0_ahl, x0_cea;
  cout << "~ Reducing instance " << argv[1] << " with AHL method" << endl;
  MakeAHLLattice(Q_ahl, x0_ahl);
  cout << "~ Reducing instance " << argv[1] << " with CEA method" << endl;
  MakeCEALattice(Q_cea, x0_cea);

  getrusage (RUSAGE_SELF, &tusage);
  ttime = (tusage.ru_utime.tv_usec + tusage.ru_stime.tv_usec) - ttime;
  cout << "~ Total time: " << ttime << " microseconds" << endl;


  // generate output file names
  string outputfile_orig = ((string)argv[2]).append("_orig.lp");
  string outputfile_ahl  = ((string)argv[2]).append("_ahl.lp");
  string outputfile_cea  = ((string)argv[2]).append("_cea.lp");


  // print original lp
  print_original((char*)outputfile_orig.c_str());


  // print both reformulations
  print_reformulation((char*)outputfile_ahl.c_str(), Q_ahl, x0_ahl);
  print_reformulation((char*)outputfile_cea.c_str(), Q_cea, x0_cea);

  // get translation of AHL if requested
  if (opt == "--translate")
  {
    cout << "~ Computing AHL translation matrix" << endl;
    print_translation(Q_ahl, x0_ahl);
  }


}
