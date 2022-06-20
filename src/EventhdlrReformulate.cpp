#include <cmath>
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
#include <scip/scipdefplugins.h>

#include "EventhdlrReformulate.hpp"

using namespace std;
using namespace NTL;
using matrix = vector<vector<int>>;

ZZ double2zz(double value)
{
  double frac = value - std::floor(value);
  if (frac>1e-8)
    throw std::invalid_argument("ERROR: coefficient not integer");
  return to_ZZ(value);
}

int double2int(double value)
{
  double frac = value - std::floor(value);
  if (frac>1e-8)
    throw std::invalid_argument("ERROR: coefficient not integer");
  return static_cast<int>(value);
}

matrix matrix_multiply(matrix A, matrix B)
{
    int m = A.size();
    int n = B[0].size();
    int p = B.size();
    assert(A[0].size() == p);

    matrix C(m, vector<int>(n, 0)); // mxn matrix of zeros
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < p; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
    return C;
}

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
  if (IsZero(X)) cout << "  Basis was correctly generated." << endl;
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
    /* get number of variables */
    int nvars;
    SCIPgetConsNVars(scip, conss[i], &nvars, &success);
    assert(success);

    /* get values and variables */
    SCIP_VAR** vars;
    SCIP_Real* vals;
    SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
    SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
    SCIP_CALL( SCIPgetConsVals(scip, conss[i], vals, nvars, &success) );
    assert(success);
    SCIP_CALL( SCIPgetConsVars(scip, conss[i], vars, nvars, &success) );
    assert(success);

    /* get lhs and rhs */
    SCIP_Real lhs = SCIPconsGetLhs(scip, conss[i], &success);
    SCIP_Real rhs = SCIPconsGetRhs(scip, conss[i], &success);

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
void Reduce(
  mat_ZZ Aext,
  mat_ZZ &Q,
  vec_ZZ &x0,
  ZZ &determ
)
{
  int i,j;
  mat_ZZ L, U, X1, Atrans;
  ZZ N1=to_ZZ(10000000);
  ZZ N2=to_ZZ(1000000000);

  int m = Aext.NumRows();
  int n = Aext.NumCols() - 1;

  /* create L matrix */
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
    L[n][n+i+1] = to_ZZ(-Aext[i][n])*(N2);

  /* lattice reduction */
  LLL(determ, L, U, 99, 100, 0);

  /* find N1 in L matrix */
    /* Note: it should appear in position [n-p][n], where
       p is the number of linearly independent rows (therefore
       smaller or equal to m). */
  int p;
  for (i=0;i<n+1;i++)
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
  for (i=0;i<n-p; i++){
    for (j=0;j<n;j++)
      Q[j][i]=L[i][j];
   }
  check_basis(Aext, Q);

  x0.SetLength(n);
  for (i=0; i<n; i++) x0[i] = L[n-p][i];

}

/* get new variable bounds */
SCIP_RETCODE get_new_varbounds(
  matrix basis,
  vector<double> lhs,
  vector<double> rhs,
  vector<int> &upper_bounds,
  vector<int> &lower_bounds
)
{
  int m = basis.size();
  int n = basis[0].size();

  SCIP* scip;
  SCIPcreate(&scip);
  SCIP_CALL( SCIPincludeDefaultPlugins(scip) ); /* include default plugins */
  SCIP_CALL( SCIPcreateProbBasic(scip, "bounds") ); /* creating the SCIP Problem. */
  SCIPsetIntParam(scip, "display/verblevel", 0); /* turn off terminal display */
  SCIPenableReoptimization(scip, TRUE);


  /* create n variables */
  vector<SCIP_VAR*> xvars(n);
  SCIP_Real inf = SCIPinfinity(scip);
  for( int k = 0; k < n; ++k )
  {
    SCIP_VAR* var = nullptr;
    string name = "x"+to_string(k+1);
    SCIP_CALL( SCIPcreateVarBasic(
                                   scip,                        /* SCIP environment */
                                   &var,                        /* reference to the variable */
                                   name.c_str(),                /* name of the variable */
                                   -inf,                        /* lower bound of the variable */
                                   inf,                         /* upper bound of the variable */
                                   0.0,                         /* obj. coefficient. */
                                   SCIP_VARTYPE_CONTINUOUS      /* variable is binary */
                                   ) );
    SCIP_CALL( SCIPaddVar(scip, var) );
    xvars[k] = var;
  }

  /* create m constraints */
  vector<SCIP_CONS*> constraints(m);
  for( int k = 0; k < m; ++k )
  {
     SCIP_CONS* cons = nullptr;
     string name = "C"+to_string(k+1);

     int nvars = 0;
     vector<double> vals;
     vector<SCIP_VAR*> vars;
     for( int j = 0; j < n; ++j )
     {
       if (basis[k][j] == 0) { continue; }
       vals.push_back( conv<double>(basis[k][j]) );
       vars.push_back( xvars[j] );
       nvars++;
     }

     SCIP_CALL( SCIPcreateConsBasicLinear(
                                           scip,
                                           &cons,                 /* pointer to hold the created constraint */
                                           name.c_str(),          /* name of constraint */
                                           nvars,                 /* number of nonzeros in the constraint */
                                           vars.data(),           /* array with variables of constraint entries */
                                           vals.data(),           /* array with coefficients of constraint entries */
                                           lhs[k],                /* left hand side of constraint */
                                           rhs[k]) );             /* right hand side of constraint */

     SCIP_CALL( SCIPaddCons(scip, cons) );
     constraints[k] = cons;
  }


  /* get upper bounds */
  for( int k = 0; k < n; ++k )
  {
    //SCIPinfoMessage(scip, NULL, "getting upper bound for var %d\n", k);
    vector<SCIP_Real> objfun(n, 0.0);
    objfun[k] = 1.0;
    SCIP_CALL( SCIPchgReoptObjective(scip,
                                     SCIP_OBJSENSE_MAXIMIZE,
                                     xvars.data(),
                                     objfun.data(),
                                     n) );
    // solve //
    SCIP_RETCODE retcode = SCIPsolve(scip);
    if ( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL )
    {
      assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
      SCIP_SOL* sol = SCIPgetBestSol(scip);
      SCIP_Real b = SCIPgetSolVal(scip, sol, xvars[k]);
      upper_bounds.push_back( ceil(b) );
    }
    else
    {
      upper_bounds.push_back( SCIPinfinity(scip) );
    }

    // free //
    SCIPfreeReoptSolve(scip);
  }

  /* get lower bounds */
  for( int k = 0; k < n; ++k )
  {
    vector<SCIP_Real> objfun(n, 0.0);
    objfun[k] = 1.0;
    SCIP_CALL( SCIPchgReoptObjective(scip,
                                     SCIP_OBJSENSE_MINIMIZE,
                                     xvars.data(),
                                     objfun.data(),
                                     n) );
    // solve //
    SCIP_RETCODE retcode = SCIPsolve(scip);
    if ( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL )
    {
      assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
      SCIP_SOL* sol = SCIPgetBestSol(scip);
      SCIP_Real b = SCIPgetSolVal(scip, sol, xvars[k]);
      lower_bounds.push_back( floor(b) );
    }
    else
    {
      lower_bounds.push_back( -SCIPinfinity(scip) );
    }

    // free //
    SCIPfreeReoptSolve(scip);
  }

  /*freeing the variables */
  for( int k = 0; k < n; ++k )
  {
     SCIP_CALL( SCIPreleaseVar(scip, &xvars[k]) );
  }
  xvars.clear();

  /* freeing the constraints */
  for( auto &constr : constraints )
  {
     SCIP_CALL( SCIPreleaseCons(scip, &constr) );
  }
  constraints.clear();

  SCIPfree(&scip);
  return SCIP_OKAY;
}

/* get new objective function */
vector<int> get_new_objfun(
  matrix basis,
  vector<int> x0,
  vector<double> c)
{
  int n = basis.size();
  int m = n - basis[0].size();

  vector<int> obj(n-m+1, 0);

  for (int i=0; i<(n-m); i++)
  {
      for (int j=0; j<n; j++)
        obj[i] += double2int(c[j])*basis[j][i];
  }

  obj[n-m]=0;
  for (int j=0; j<n; j++)
    obj[n-m]+=double2int(c[j])*x0[j];
  obj[n-m]+=double2int(c[n]);

  return obj;
}


/* get filename of reformulated instance */
string get_new_filename(
  string path
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
  string filename = dir + "/ahl_" + file + ".lp";
  return filename;
}

void nnegative_quad_transform(
  SCIP* scip,
  matrix &basis,
  vector<double> &lhs,
  vector<double> &rhs,
  vector<int> &upper_bounds,
  vector<int> &lower_bounds,
  vector<int> &obj
)
{
  int n = basis.size();
  int p = lower_bounds.size();

  /* mu' =  A*mu + b */
  vector<int> b(p, 0);
  matrix A(p, vector<int>(p, 0)); // pxp matrix of all zeros

  for (int i=0; i<p; i++)
  {
    if (lower_bounds[i] > -1e9)
    {
      A[i][i] = 1;
      b[i] = -lower_bounds[i];
      upper_bounds[i] = upper_bounds[i] - lower_bounds[i];
      lower_bounds[i] = 0;
    }
    else
    {
      assert(upper_bounds[i]<1e9);
      A[i][i] = -1;
      b[i] = upper_bounds[i];
      upper_bounds[i] = SCIPinfinity(scip);
      lower_bounds[i] = 0;
      //SCIPinfoMessage(scip, NULL, "infinite lower bound. New upper bound %d\n", upper_bounds[i]);
    }
  }

  /* update basis */
  basis = matrix_multiply(basis, A);

  /* update rhs and lhs */
  for (int i=0; i<n; i++)
  {
    int shift = 0;
    for (int j=0; j<p; j++)
    {
      shift += basis[i][j]*b[j];
    }
    lhs[i] += shift;
    rhs[i] += shift;
  }

  /* update objective function */
  int shift = 0;
  for (int i=0; i<p; i++)
  {
    obj[i] = obj[i]*A[i][i];
    shift += obj[i]*b[i];
  }
  obj[p] = obj[p] - shift;

}

/* print reformulation */
void print_original(
  SCIP* scip,
  const char *filename,
  mat_ZZ Aext,
  vector<double> upper,
  vector<double> lower,
  vector<double> objfun,
  bool maximization
)
{
  int i, j, k = 0;
  ofstream output_file(filename);

  int n = Aext.NumCols() - 1;
  int m = Aext.NumRows();


  /* write objective function */
  if (maximization) { output_file << "maximize "; }
  else { output_file << "minimize "; }
  for (j=0; j<n; j++)
  {
    if (objfun[j] != 0)
    {
      if (k>20){k=0; output_file << "\n";}
      if (objfun[j]>0)
      {
        output_file << "+"<<objfun[j]<<" k"<<j+1<< " ";
      }
      else
      {
        output_file << objfun[j]<<" k"<<j+1<< " ";
      }
      k++;
    }
  }
  if (objfun[n] > 0) output_file << "+ " << objfun[n];
  if (objfun[n] < 0) output_file << objfun[n];
  output_file << "\n";


  /* write constraints */
  int conss_counter = 0;
  output_file << "subject to" << "\n";
  for (i=0;i<m;i++)
  {
    if (Aext[i][n]>1e9) continue;
    k=0;
    conss_counter++;
    output_file << "C" << conss_counter << ": ";
    for (j=0; j<n; j++)
    {
      if (Aext[i][j] != 0)
      {
        if (k>20){k=0; output_file << "\n";}
        if (Aext[i][j] > 0)
        {
          output_file << "+"<<Aext[i][j]<<" k"<<j+1<<" ";
        }
        else
        {
          output_file << Aext[i][j]<<" k"<<j+1<<" ";
        }
        k++;
      }
    }
    output_file << " = ";
    output_file << Aext[i][n] << "\n";
  }

  /* write variable bounds */
  output_file << "Bounds\n";
  for (i=0; i<n; i++)
  {
    if (lower[i]<-1e9)
    {
      output_file << "-inf";
    }
    else
    {
      output_file << lower[i];
    }
    output_file << "<= k" << i+1 << "<=";
    if (upper[i]>1e9)
    {
      output_file << "inf" << "\n";
    }
    else
    {
      output_file << upper[i] << "\n";
    }
  }

  /* write variables */
  output_file << "Generals " << "\n";
  k=0;
  for (i=0; i<n; i++)
  {
    if (k>20){k=0; output_file << "\n";}
    output_file << "k" << i+1 << " ";
    k++ ;
  }


  output_file << "\nEnd \n";
  output_file << "\n\n";
  output_file.close();
}



/* print reformulation */
void print_reformulation(
  SCIP* scip,
  const char *filename,
  matrix basis,
  vector<int> x0,
  vector<double> upper,
  vector<double> lower,
  vector<double> objfun,
  bool maximization
)
{
  int i, j, k = 0;
  ofstream output_file(filename);

  int n = basis.size();
  int m = n - basis[0].size();

  /* get new objective function */
  vector<int> obj = get_new_objfun(basis, x0, objfun);

  /* get rhs and lhs of constraints */
  vector<double> lhs, rhs;
  for (j = 0; j < n; j++)
  {
    lhs.push_back( lower[j] - conv<double>(x0[j]) );
    rhs.push_back( upper[j] - conv<double>(x0[j]) );
  }

  /* get bounds of new variables */
  vector<int> newupper, newlower;
  SCIPinfoMessage(scip, NULL, "   calculating bounds for new variables.\n");
  SCIP_RETCODE retcode = get_new_varbounds(basis, lhs, rhs, newupper, newlower);

  /* transform variables to nonnegative quadrant */
  //nnegative_quad_transform(scip, basis, lhs, rhs, newupper, newlower, obj);

  /* write objective function */
  if (maximization) { output_file << "maximize "; }
  else { output_file << "minimize "; }
  for (j=0; j<(n-m); j++)
  {
    if (obj[j] != 0)
    {
      if (k>20){k=0; output_file << "\n";}
      if (obj[j]>0)
      {
        output_file << "+"<<obj[j]<<" k"<<j+1<< " ";
      }
      else
      {
        output_file << obj[j]<<" k"<<j+1<< " ";
      }
      k++;
    }
  }
  if (obj[n-m] > 0) output_file << "+ " << obj[n-m];
  if (obj[n-m] < 0) output_file << obj[n-m];
  output_file << "\n";


  /* write constraints */
  int conss_counter = 0;
  output_file << "subject to" << "\n";
  for (j=0;j<n;j++) //  Q\mu >= l-x0
  {
    if (lhs[j]<-1e9) continue;
    k=0;
    conss_counter++;
    output_file << "C" << conss_counter << ": ";
    for (i=0; i<n-m; i++)
    {
      if (basis[j][i] != 0)
      {
        if (k>20){k=0; output_file << "\n";}
        if (basis[j][i] > 0)
        {
          output_file << "+"<<basis[j][i]<<" k"<<i+1<<" ";
        }
        else
        {
          output_file << basis[j][i]<<" k"<<i+1<<" ";
        }
        k++;
      }
    }
    output_file << " >= ";
    output_file << lhs[j] << "\n";
  }
  for (j=0;j<n;j++) //  Q\mu <= u-x0
  {
    if (rhs[j]>1e9) continue;
    k=0;
    conss_counter++;
    output_file << "C" << conss_counter << ": ";
    for (i=0; i<n-m; i++)
    {
      if (basis[j][i] != 0)
      {
        if (k>20){k=0; output_file << "\n";}
        if (basis[j][i] > 0)
        {
          output_file << "+"<<basis[j][i]<<" k"<<i+1<<" ";
        }
        else
        {
          output_file << basis[j][i]<<" k"<<i+1<<" ";
        }
        k++;
      }
    }
    output_file << " <= ";
    output_file << rhs[j] << "\n";
  }

  /* write variable bounds */
  output_file << "Bounds\n";
  for (i=0; i<(n-m); i++)
  {
    if (newlower[i]<-1e9)
    {
      output_file << "-inf";
    }
    else
    {
      output_file << newlower[i];
    }
    output_file << "<= k" << i+1 << "<=";
    if (newupper[i]>1e9)
    {
      output_file << "inf" << "\n";
    }
    else
    {
      output_file << newupper[i] << "\n";
    }
  }

  /* write variables */
  output_file << "Generals " << "\n";
  k=0;
  for (i=0; i<(n-m); i++)
  {
    if (k>20){k=0; output_file << "\n";}
    output_file << "k" << i+1 << " ";
    k++ ;
  }


  output_file << "\nEnd \n";
  output_file << "\n\n";
  output_file.close();
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
SCIP_DECL_EVENTFREE(EventhdlrReformulate::scip_free)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */

SCIP_DECL_EVENTINIT(EventhdlrReformulate::scip_init)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */

SCIP_DECL_EVENTEXIT(EventhdlrReformulate::scip_exit)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The event handler may use this call to initialize its branch and bound specific data.
 *
 */

SCIP_DECL_EVENTINITSOL(EventhdlrReformulate::scip_initsol)
{
   SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The event handler should use this call to clean up its branch and bound data.
 */

SCIP_DECL_EVENTEXITSOL(EventhdlrReformulate::scip_exitsol)
{
   return SCIP_OKAY;
}


/** frees specific constraint data */

SCIP_DECL_EVENTDELETE(EventhdlrReformulate::scip_delete)
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
SCIP_DECL_EVENTEXEC(EventhdlrReformulate::scip_exec)
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


  /* save matrix */
  bool save_mat = FALSE;
  if (save_mat)
  {
    ofstream output_file("contraint_matrix.txt");
    output_file << n << " " << m << "\n";
    for (int i=0;i<m;i++)
    {
      for (int j=0;j<n;j++)
        output_file << Aext[i][j] << " ";
      output_file << -Aext[i][n];
      output_file << "\n";
    }
    for (int j=0;j<n;j++)
      output_file << lower[j] << " ";
    output_file << "\n";
    for (int j=0;j<n;j++)
      output_file << upper[j] << " ";
    output_file << "\n";
    output_file.close();
    string filename2 = "test.lp";
    print_original(scip, filename2.c_str() , Aext, upper, lower, objfun, maximization);
  }


  /* get kernel basis */
  mat_ZZ Q; vec_ZZ x0; ZZ determ;
  SCIPinfoMessage(scip, NULL, "    reducing kernel basis\n");
  try
  {
    Reduce(Aext, Q, x0, determ);
  }
  catch (const std::runtime_error& e)
  {
    cout << e.what() << endl;
    return SCIP_ERROR;
  }
  m = n - Q.NumCols(); // update m. This is now the number of l.i. rows!
  bool print_Q = FALSE; bool write_determ = FALSE;
  if (print_Q)
  {
    ofstream output_file2("Q.txt");
    for (int i=0;i<n;i++)
    {
      for (int j=0;j<n-m;j++)
        output_file2 << -Q[i][j] << " ";
      output_file2 << "\n";
    }
    for (int j=0;j<n;j++)
      output_file2 << x0[j] << " ";
    output_file2 << "\n";
  }
  if (write_determ)
  {
    fstream log;
    log.open("determinants.txt", fstream::app);
    log << instancepath << " " << determ << "\n";
    log.close();
  }


  // get tranlation of Q //
  // long ret;
  // vec_ZZ x, e_i;
  // mat_ZZ branchmat, test_matrix;
  // branchmat.SetDims(n-m,n);
  // x.SetLength(n);
  // e_i.SetLength(n-m);
  //
  // for (int i=0;i<n-m;i++)
  // {
  //   for (int j=0;j<n-m;j++) e_i[j] = 0;
  //   e_i[i] = to_ZZ(1);
  //   ret = LatticeSolve(x, Q, e_i, 1);
  //   assert(ret);
  //   for (int j=0;j<n;j++)
  //   {
  //     branchmat[i][j] = x[j];
  //   }
  // }
  //
  // mul(test_matrix, branchmat, Q);
  // long flag = IsIdent(test_matrix, n-m);
  // assert(flag);
  //
  // ofstream output_file3("W_matrix.txt");
  // for (int i=0;i<n-m;i++)
  // {
  //   for (int j=0;j<n;j++)
  //     output_file3 << branchmat[i][j] << " ";
  //   output_file3 << "\n";
  // }
  // output_file3.close();

  /* print the reformulated problem */
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
  string filename = get_new_filename(instancepath);
  print_reformulation(scip, filename.c_str() , basis, x, upper, lower, objfun, maximization);


  SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1) );
  SCIP_CALL( SCIPinterruptSolve(scip)	);
  return SCIP_OKAY;
}
