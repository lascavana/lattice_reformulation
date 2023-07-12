#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <unordered_map>

// NTL library
#include <NTL/ZZ.h>
#include <NTL/LLL.h>
#include <NTL/mat_ZZ.h>

// SCIP
#include <scip/scip.h>

#include "Eventhdlr.hpp"
#include "utils.hpp"

using namespace std;
using namespace NTL;
using matrix = vector<vector<int>>;


void check_basis(mat_ZZ Aext, mat_ZZ Q)
{
  int i, j;
  mat_ZZ A, X;

  int m = Aext.NumRows();
  int n = Aext.NumCols() - 1;

  // create non-extended version of matrix
  A.SetDims(m,n);
  for (i = 0; i < m; i++)
  {
      for (j = 0; j < n; j++) A[i][j] = Aext[i][j];
  }

  mul(X, A, Q);
  if (IsZero(X)) cout << "    basis was correctly generated." << endl;
  else
  {
    cout << "Error in reduction." << endl;
    cout << "A*Q is not zero." << endl;
    throw  std::bad_function_call();
  }
}


/* get constraint matrix */
SCIP_RETCODE GetInstanceData(
  SCIP* scip,
  mat_ZZ &consmat,
  vector<double> &upper,
  vector<double> &lower,
  vector<double> &objfun,
  bool &maximization
)
{
  SCIPinfoMessage(scip, NULL, "    retreiving constraint matrix\n");

  /* check stage */
	assert(scip != nullptr);
	if (SCIPgetStage(scip) != SCIP_STAGE_SOLVING) {
    SCIPerrorMessage("cannot branch when not solving\n");
		return SCIP_INVALIDCALL;
	}

  /* get constraints */
  SCIP_CONS** conss = SCIPgetConss(scip);

  /* set dimensions */
  int n = SCIPgetNVars(scip);
  int m = SCIPgetNConss(scip);
  mat_ZZ Aext;
  Aext.SetDims(m,n+2); // Aext = [lhs|A|rhs]
  SCIPinfoMessage(scip, NULL, "    Aext dim: (%d, %d)\n", m, n+2);

  /* assign var index to a matrix colum */
  unordered_map<int, int> idx2col;
  SCIP_VAR** allvars = SCIPgetVars(scip);
  for (int i = 0; i < n; ++i)
  {
    	int index = SCIPvarGetIndex(allvars[i]);
      idx2col.insert({index, i});
  }

  /* get contraint data */
  int n2 = n, m2 = m;
  unsigned int success = TRUE;
  for (int i = 0; i < m; ++i)
  {
    SCIP_CONS* cons = conss[i];
    /* get number of variables */
    int nvars;
    SCIPgetConsNVars(scip, cons, &nvars, &success);
    assert(success);

    /* get values and variables */
    SCIP_VAR** vars;
    SCIP_Real* vals;
    SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
    SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
    SCIP_CALL( SCIPgetConsVals(scip, cons, vals, nvars, &success) );
    assert(success);
    SCIP_CALL( SCIPgetConsVars(scip, cons, vars, nvars, &success) );
    assert(success);

    /* get lhs and rhs */
    SCIP_Real lhs = SCIPconsGetLhs(scip, cons, &success);
    SCIP_Real rhs = SCIPconsGetRhs(scip, cons, &success);

    /* handle constraint */
    if (rhs-lhs > 0) // different
    {
      n2++; // will need a slack variable
      if (lhs>-1e15 and rhs<1e15) // both finite, need extra row and var
      {
        n2++; m2++;
      }
    }

    /* fill in matrix entries */
    for (int j = 0; j < nvars; ++j)
    {
      // check that variables are not continuous (not implemented yet)
      // SCIP_VARTYPE vartype = SCIPvarGetType(vars[j]);
      //assert( vartype !=  SCIP_VARTYPE_CONTINUOUS );
      if (vals[j] == 0) continue;
      int col = idx2col[ SCIPvarGetIndex(vars[j]) ];
      try
      {
        Aext[i][col+1] = double2zz(vals[j]);
      }
      catch (const std::invalid_argument& e)
      {
        cout << e.what() << endl;
        return SCIP_ERROR;
      }
    }

    try
    {
      Aext[i][0] = double2zz(lhs);
      Aext[i][n+1] = double2zz(rhs);
    }
    catch (const std::invalid_argument& e)
    {
      cout << e.what() << endl;
      return SCIP_ERROR;
    }

    /* release arrays */
    SCIPfreeBufferArray(scip, &vals);
    SCIPfreeBufferArray(scip, &vars);
  }
  SCIPinfoMessage(scip, NULL, "    Added %d slack variables ", n2-n);
  SCIPinfoMessage(scip, NULL, "and %d constraints. \n", m2-m);

  /* post-process data and save into consmat */
  int idx = 0;
  int sign = 1;
  int slack_counter= 0;
  consmat.SetDims(m2,n2+1); // consmat = [A|b]
  for (int i = 0; i < min(m,m2); ++i)
  {
    ZZ b;
    ZZ lhs = Aext[i][0];
    ZZ rhs = Aext[i][n+1];

    // change sign if necessary
    if (rhs > 1e15) // lhs =< a*x
    {
      sign = -1;
      b = -lhs;
    }
    else
    {
      sign = 1;
      b = rhs;
    }

    for (int j = 0; j < n; ++j) consmat[idx][j] = sign*Aext[i][j+1];
    for (int j = n; j < n2; ++j) consmat[idx][j] = 0; // slack variables
    consmat[idx][n2] = b;

    // inequality contraints
    if (rhs-lhs > 0)
    {
      consmat[idx][n+slack_counter] = 1;
      slack_counter++;

      if (lhs > -1e15 && rhs < 1e15) // both finite, need two constraints
      {
        idx++;
        for (int j = 0; j < n; ++j) consmat[idx][j] = -Aext[i][j+1];
        for (int j = n; j < n2; ++j) consmat[idx][j] = 0; // slack variables
        consmat[idx][n+slack_counter] = 1;
        consmat[idx][n2] = -lhs;
        slack_counter++;

      }
    }

    idx++;
  }
  assert(idx==m2);


  // get objective function and bounds//
  /*
    note: if orig problem is maximization, the obj fun coefficients
    returned by SCIPvarGetObj will have opposite sign. The variable
    s corrects for this.
  */
  SCIP_Real objscale = SCIPgetTransObjscale(scip);
  SCIP_OBJSENSE objsense = SCIPgetObjsense(scip);
  if ( objsense == SCIP_OBJSENSE_MINIMIZE )
  {
    maximization = FALSE;
  }
  else
  {
    maximization = TRUE;
    objscale *= -1.0;
  }
  objfun.resize(n2+1); upper.resize(n2); lower.resize(n2);
  for (int i = 0; i < n; ++i)
  {
    	int col = idx2col[ SCIPvarGetIndex(allvars[i]) ] ;
      objfun[col] = objscale*SCIPvarGetObj(allvars[i]);
      upper[col] = SCIPvarGetUbLocal(allvars[i]);
      lower[col] = SCIPvarGetLbLocal(allvars[i]);
  }
  objfun[n2] = objscale*SCIPgetTransObjoffset(scip);
  for (int i = n; i < n2; ++i) { lower[i] = 0; }


  return SCIP_OKAY;
}


/* lattice reduction */
void reduce_ahl(
  const mat_ZZ  &Aext,
  mat_ZZ        &Q,
  vec_ZZ        &x0,
  ZZ            &determ,
  bool          highquality
)
{
  mat_ZZ L, U, X1, Atrans;
  ZZ N1=to_ZZ(10000000);
  ZZ N2=to_ZZ(1000000000);

  int m = Aext.NumRows();
  int n = Aext.NumCols() - 1;

  /* create L matrix */
  L.SetDims(n+1,n+m+1);
  for (int i=0;i<n+1;i++)
  {
     for(int j=0;j<n+m+1;j++)
         L[i][j]=0;
     L[i][i]=to_ZZ(1);
  }
  L[n][n]= N1;
  for (int j=0;j<n;j++){
     for (int i=0; i<m; i++)
        L[j][n+i+1] = to_ZZ(Aext[i][j])*N2;
  }
  for (int i=0; i<m; i++)
    L[n][n+i+1] = to_ZZ(-Aext[i][n])*(N2);

  /* lattice reduction */
  if (highquality)
  {
    LLL(determ, L, U, 99, 100, 0);
  }
  else
  {
    LLL(determ, L, U, 40, 100, 0);
  }
  

  /* find N1 in L matrix */
    /* Note: it should appear in position [n-p][n], where
       p is the number of linearly independent rows (therefore
       smaller or equal to m). */
  int p;
  for (int i=0;i<n+1;i++)
  {
    if (L[i][n] == 0) continue;
    else
    {
      if (L[i][n] != N1 && L[i][n] != -N1)
        throw std::runtime_error("ERROR: N1 not found. Instance infeasible");
      p = n-i;
      break;
    }
  }
  assert(p<=m);

  /* get kernel lattice basis */
  Q.SetDims(n,n-p);
  for (int i=0;i<n-p; i++){
    for (int j=0;j<n;j++)
      Q[j][i]=L[i][j];
   }
  check_basis(Aext, Q);

  x0.SetLength(n);
  for (int i=0; i<n; i++) x0[i] = L[n-p][i];

}

mat_ZZ regularize_Q(
  mat_ZZ a,
  mat_ZZ Q
)
{
  int m = a.NumRows();
  int n = a.NumCols()-1;
  assert(n==Q.NumRows());

  /* get (D*Q)^trans */
  mat_ZZ DQT;
  DQT.SetDims(n-1,n);
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n-1; j++)
    {
      DQT[j][i] = a[0][i]*Q[i][j];
    }
  }

  /* reduce */
  ZZ determ; mat_ZZ U;
  LLL(determ, DQT, U, 99, 100, 0);

  /* create new mat */
  mat_ZZ Qnew;
  Qnew.SetDims(n,n-1);
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n-1; j++)
    {
      Qnew[i][j] = DQT[j][i] / a[0][i];
    }
  }
  check_basis(a, Qnew);

  return Qnew;
}

void get_volume_bound(
  mat_ZZ a
)
{
  long m = a.NumRows();
  long n = a.NumCols()-1;
  assert(m==1);

  double c = 1/0.74;

  double f1 = (n-1)*(n-2)*log(3*sqrt(c)/2)/2 + (n-1)*log(2) - log(n)/2;
  double f2 = (n-1)*log(conv<double>(a[0][n]));
  for (int j=0; j<n; j++) f2 -= log(conv<double>(a[0][j]));

  cout << "First term: " << f1 << endl;
  cout << "Second term: " << f2 << endl;
}

/* get filename of reformulated instance */
string get_new_filename(
  string path,
  bool highquality,
  bool regularize
)
{
  size_t dirpos = path.find_last_of("/");
  // get directory
  string dir = path.substr(0, dirpos);
  // get file
  string file = path.substr(dirpos+1);
  // remove extension
  size_t extpos = file.find_last_of(".");
  file = file.substr(0, extpos);
  string filename;
  if (regularize)
  {
    filename = dir + "/ahl_diag_" + file + ".lp";
  }
  else if (!highquality)
  {
    filename = dir + "/ahl_poor_" + file + ".lp";
  }
  else
  {
    filename = dir + "/ahl_" + file + ".lp";
  }
   
  return filename;
}


/** destructor of event handler to free user data (called when SCIP is exiting) */
SCIP_DECL_EVENTFREE(Eventhdlr_AHL::scip_free)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */

SCIP_DECL_EVENTINIT(Eventhdlr_AHL::scip_init)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */

SCIP_DECL_EVENTEXIT(Eventhdlr_AHL::scip_exit)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The event handler may use this call to initialize its branch and bound specific data.
 *
 */

SCIP_DECL_EVENTINITSOL(Eventhdlr_AHL::scip_initsol)
{
   SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The event handler should use this call to clean up its branch and bound data.
 */

SCIP_DECL_EVENTEXITSOL(Eventhdlr_AHL::scip_exitsol)
{
   return SCIP_OKAY;
}


/** frees specific constraint data */

SCIP_DECL_EVENTDELETE(Eventhdlr_AHL::scip_delete)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** execution method of event handler
 *
 *  Processes the event. The method is called every time an event occurs, for which the event handler
 *  is responsible. Event handlers may declare themselves resposible for events by calling the
 *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
 *  given event handler and event data.
 */
SCIP_DECL_EVENTEXEC(Eventhdlr_AHL::scip_exec)
{  /*lint --e{715}*/

  /* read problem data */
  mat_ZZ Aext;
  bool maximization;
  vector<double> objfun(10000, 0.0);
  vector<double> upper(10000, SCIPinfinity(scip));
  vector<double> lower(10000, -SCIPinfinity(scip));
  SCIP_CALL( GetInstanceData(scip, Aext, upper, lower, objfun, maximization) );
  int m = Aext.NumRows();
  int n = Aext.NumCols() - 1;
  SCIPinfoMessage(scip, NULL, "Shape of Aext: (%d,%d) \n", m,n+1);


  /* get kernel basis */
  mat_ZZ Q; vec_ZZ x0; ZZ determ;
  SCIPinfoMessage(scip, NULL, "    reducing kernel basis...\n");
  auto start = chrono::high_resolution_clock::now();
  try
  {
    reduce_ahl(Aext, Q, x0, determ, highquality);
  }
  catch (const std::runtime_error& e)
  {
    cout << e.what() << endl;
    return SCIP_ERROR;
  }
  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
  SCIPinfoMessage(scip, NULL, "    took %lld ms to reduce \n", duration.count());

  /* update m. This is now the number of l.i. rows! */
  m = n - Q.NumCols(); 

  /* regularize */
  if (diag)
  {
    assert(m==1);
    Q = regularize_Q(Aext, Q);
    get_volume_bound(Aext);
  }

  /* print Q */
  bool print_Q = false; 
  if (print_Q)
  {
    ofstream output_file2("Q.txt");
    for (int i=0;i<n;i++)
    {
      for (int j=0;j<n-m;j++)
        output_file2 << Q[i][j] << " ";
      output_file2 << "\n";
    }
  }

  /* convert to int type */
  vector<int> x(n, 0);
  matrix basis(n, vector<int>(n-m, 0));
  for (int i=0;i<n;i++)
  {
    x[i] = conv<int>(x0[i]);
    for (int j=0;j<n-m;j++)
    {
      basis[i][j] = conv<int>(Q[i][j]);
    }
  }

  /* print the reformulated problem */
  string filename = get_new_filename(instancepath, highquality, diag);
  print_reformulation(scip, filename.c_str() , basis, x, upper, lower, objfun, maximization);


  SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1) );
  SCIP_CALL( SCIPinterruptSolve(scip)	);
  return SCIP_OKAY;
}
