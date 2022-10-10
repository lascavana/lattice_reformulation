#include <vector>
#include <fstream>
#include <iostream>

#include "utils.hpp"

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

void nnegative_quad_transform(
  SCIP* scip,
  matrix &basis,
  vector<double> &lhs,
  vector<double> &rhs,
  vector<int> &upper_bounds,
  vector<int> &lower_bounds,
  vector<double> &obj
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

/* get new variable bounds */
SCIP_RETCODE get_new_varbounds(
  matrix basis,
  vector<int> lhs,
  vector<int> rhs,
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
vector<double> get_new_objfun(
  matrix basis,
  vector<int> x0,
  vector<double> c)
{
  int n = basis.size();
  int m = n - basis[0].size();

  vector<double> obj(n-m+1, 0);

  for (int i=0; i<(n-m); i++)
  {
      for (int j=0; j<n; j++)
        obj[i] += c[j]*conv<double>(basis[j][i]);
  }

  obj[n-m]=0;
  for (int j=0; j<n; j++)
    obj[n-m]+=c[j]*conv<double>(x0[j]);
  obj[n-m]+=c[n];

  return obj;
}

void print_for_ls(
  string filename,
  matrix basis,
  vector<int> lhs,
  vector<int> rhs,
  vector<int> upper,
  vector<int> lower
)
{
  int n = lhs.size();
  int nvars = basis[0].size();

  int nconss = 0;
  for (int i=0;i<n;i++)
  {
    if (rhs[i]<=1e9) nconss++;
    if (lhs[i]>=-1e9) nconss++;
  }

  ofstream output_file2(filename);
  output_file2 << nvars << " " << nconss << endl;

  /* write constraint matrix */
  for (int i=0;i<n;i++)
  {
    if (rhs[i]>1e9) continue;
    for (int j=nvars-1;j>-1;j--)
      output_file2 << basis[i][j] << " ";
    output_file2 << "\n";
  }
  for (int i=0;i<n;i++)
  {
    if (lhs[i]<-1e9) continue;
    for (int j=nvars-1;j>-1;j--)
      output_file2 << -basis[i][j] << " ";
    output_file2 << "\n";
  }

  /* write b */
  for (int i=0;i<n;i++) 
  {
    if (rhs[i]>1e9) continue;
    output_file2 << rhs[i] << " ";
  }
  for (int i=0;i<n;i++) 
  {
    if (lhs[i]<-1e9) continue;
    output_file2 << -lhs[i] << " ";
  }
  output_file2 << "\n";

  /* write upper and lower bounds */
  for (int j=nvars-1;j>-1;j--) output_file2 << lower[j] << " ";
  output_file2 << "\n";
  for (int j=nvars-1;j>-1;j--) output_file2 << upper[j] << " ";
  output_file2 << "\n";

  output_file2.close();

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
  vector<double> obj = get_new_objfun(basis, x0, objfun);

  /* get rhs and lhs of constraints */
  vector<int> lhs, rhs;
  for (j = 0; j < n; j++)
  {
    if (SCIPisInfinity(scip, -lower[j]))
    {
      lhs.push_back(-SCIPinfinity(scip));
    }
    else
    {
      lhs.push_back( conv<int>(lower[j]) - x0[j] );
    }

    if (SCIPisInfinity(scip, upper[j]))
    {
      rhs.push_back(SCIPinfinity(scip));
    }
    else
    {
      rhs.push_back( conv<int>(upper[j]) - x0[j] );
    }
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

  print_for_ls("LS.txt", basis, lhs, rhs, newupper, newlower);
}


