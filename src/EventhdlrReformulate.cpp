#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>

// NTL library
#include <NTL/ZZ.h>
#include <NTL/LLL.h>
#include <NTL/mat_ZZ.h>

// SCIP
#include <scip/scip.h>
#include "scip/type_cons.h"
#include <scip/scipdefplugins.h>

#include "EventhdlrReformulate.hpp"

using namespace std;
using namespace NTL;

/** constraint data for linear constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of row (for ranged rows) */
   SCIP_Real             rhs;                /**< right hand side of row */
   SCIP_Real             maxabsval;          /**< maximum absolute value of all coefficients */
   SCIP_Real             minabsval;          /**< minimal absolute value of all coefficients */
   SCIP_Real             minactivity;        /**< minimal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             maxactivity;        /**< maximal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             lastminactivity;    /**< last minimal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             lastmaxactivity;    /**< last maximal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             glbminactivity;     /**< minimal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             glbmaxactivity;     /**< maximal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             lastglbminactivity; /**< last global minimal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             lastglbmaxactivity; /**< last global maximal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             maxactdelta;        /**< maximal activity contribution of a single variable, or SCIP_INVALID if invalid */
   SCIP_VAR*             maxactdeltavar;     /**< variable with maximal activity contribution, or NULL if invalid */
   uint64_t              possignature;       /**< bit signature of coefficients that may take a positive value */
   uint64_t              negsignature;       /**< bit signature of coefficients that may take a negative value */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   SCIP_NLROW*           nlrow;              /**< NLP row, if constraint has been added to NLP relaxation */
   SCIP_VAR**            vars;               /**< variables of constraint entries */
   SCIP_Real*            vals;               /**< coefficients of constraint entries */
   SCIP_EVENTDATA**      eventdata;          /**< event data for bound change events of the variables */
   int                   minactivityneginf;  /**< number of coefficients contributing with neg. infinite value to minactivity */
   int                   minactivityposinf;  /**< number of coefficients contributing with pos. infinite value to minactivity */
   int                   maxactivityneginf;  /**< number of coefficients contributing with neg. infinite value to maxactivity */
   int                   maxactivityposinf;  /**< number of coefficients contributing with pos. infinite value to maxactivity */
   int                   minactivityneghuge; /**< number of coefficients contributing with huge neg. value to minactivity */
   int                   minactivityposhuge; /**< number of coefficients contributing with huge pos. value to minactivity */
   int                   maxactivityneghuge; /**< number of coefficients contributing with huge neg. value to maxactivity */
   int                   maxactivityposhuge; /**< number of coefficients contributing with huge pos. value to maxactivity */
   int                   glbminactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbminactivity */
   int                   glbminactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbminactivity */
   int                   glbmaxactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbmaxactivity */
   int                   glbmaxactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbmaxactivity */
   int                   glbminactivityneghuge;/**< number of coefficients contrib. with huge neg. value to glbminactivity */
   int                   glbminactivityposhuge;/**< number of coefficients contrib. with huge pos. value to glbminactivity */
   int                   glbmaxactivityneghuge;/**< number of coefficients contrib. with huge neg. value to glbmaxactivity */
   int                   glbmaxactivityposhuge;/**< number of coefficients contrib. with huge pos. value to glbmaxactivity */
   int                   varssize;           /**< size of the vars- and vals-arrays */
   int                   nvars;              /**< number of nonzeros in constraint */
   int                   nbinvars;           /**< the number of binary variables in the constraint, only valid after
                                              *   sorting in stage >= SCIP_STAGE_INITSOLVE
                                              */
   unsigned int          boundstightened:2;  /**< is constraint already propagated with bound tightening? */
   unsigned int          rangedrowpropagated:2; /**< did we perform ranged row propagation on this constraint?
                                                 *   (0: no, 1: yes, 2: with potentially adding artificial constraint */
   unsigned int          validmaxabsval:1;   /**< is the maximum absolute value valid? */
   unsigned int          validminabsval:1;   /**< is the minimum absolute value valid? */
   unsigned int          validactivities:1;  /**< are the activity bounds (local and global) valid? */
   unsigned int          validminact:1;      /**< is the local minactivity valid? */
   unsigned int          validmaxact:1;      /**< is the local maxactivity valid? */
   unsigned int          validglbminact:1;   /**< is the global minactivity valid? */
   unsigned int          validglbmaxact:1;   /**< is the global maxactivity valid? */
   unsigned int          presolved:1;        /**< is constraint already presolved? */
   unsigned int          removedfixings:1;   /**< are all fixed variables removed from the constraint? */
   unsigned int          validsignature:1;   /**< is the bit signature valid? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          normalized:1;       /**< is the constraint in normalized form? */
   unsigned int          upgradetried:1;     /**< was the constraint already tried to be upgraded? */
   unsigned int          upgraded:1;         /**< is the constraint upgraded and will it be removed after preprocessing? */
   unsigned int          indexsorted:1;      /**< are the constraint's variables sorted by type and index? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the constraint already extracted? */
   unsigned int          implsadded:1;       /**< were the implications of the constraint already extracted? */
   unsigned int          coefsorted:1;       /**< are variables sorted by type and their absolute activity delta? */
   unsigned int          varsdeleted:1;      /**< were variables deleted after last cleanup? */
   unsigned int          hascontvar:1;       /**< does the constraint contain at least one continuous variable? */
   unsigned int          hasnonbinvar:1;     /**< does the constraint contain at least one non-binary variable? */
   unsigned int          hasnonbinvalid:1;   /**< is the information stored in hasnonbinvar and hascontvar valid? */
   unsigned int          checkabsolute:1;    /**< should the constraint be checked w.r.t. an absolute feasibilty tolerance? */
};

int to_int(float value)
{
  float frac = value - std::floor(value);
  assert(frac<1e-8);
  return static_cast<int>(value);
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

  // check stage //
	assert(scip != nullptr);
	if (SCIPgetStage(scip) != SCIP_STAGE_SOLVING) {
    SCIPerrorMessage("cannot branch when not solving\n");
		return SCIP_INVALIDCALL;
	}

  // get linear constraints //
  SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, "linear");
  SCIP_CONS** conss = SCIPconshdlrGetConss(conshdlr);

  // set dimensions //
  int m = SCIPconshdlrGetNConss(conshdlr);
  int n = SCIPgetNVars(scip);
  mat_ZZ Aext;
  Aext.SetDims(m,n+2); // Aext = [lhs|A|rhs]

  // shift variable indexing //
  int shift = 99999999;
  SCIP_VAR** probvars = SCIPgetVars(scip);
  for (int i = 0; i < n; ++i)
  {
    	int index = SCIPvarGetIndex(probvars[i]);
      if (index < shift) shift = index;
  }

  /* get contraint data */
  int n2 = n, m2 = m;
  for (int i = 0; i < m; ++i)
  {
    SCIP_CONSDATA* data = SCIPconsGetData(conss[i]);
    int nvars = data->varssize;
    SCIP_Real* vals = data->vals;
    SCIP_VAR** vars = data->vars;
    SCIP_Real rhs = data->rhs;
    SCIP_Real lhs = data->lhs;

    // handle constraint
    if (rhs-lhs > 0) // different
    {
      n2++; // will need a slack variable
      if (lhs>-1e15 and rhs<1e15) // both finite, need extra row and var
      {
        n2++; m2++;
      }
    }

    // fill in matrix entries
    for (int j = 0; j < nvars; ++j)
    {
      // check that variables are not continuous (not implemented yet)
      // SCIP_VARTYPE vartype = SCIPvarGetType(vars[j]);
      // assert( vartype !=  SCIP_VARTYPE_CONTINUOUS );

      int val = to_int(vals[j]);
      int idx = SCIPvarGetIndex(vars[j]) - shift;

      Aext[i][idx+1] = val;
    }
    Aext[i][0] = lhs;
    Aext[i][n+1] = rhs;
    // SCIPinfoMessage(scip, NULL, "lhs=%ld , rhs=%ld\n", conv<long>(Aext[i][0]), conv<long>(Aext[i][n+1]));

  }


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
  objfun.resize(n2); upper.resize(n2); lower.resize(n2);
  for (int i = 0; i < n; ++i)
  {
    	int index = SCIPvarGetIndex(probvars[i]) - shift;
      objfun[index] = SCIPvarGetObj(probvars[i]);
      upper[index] = SCIPvarGetUbLocal(probvars[i]);
      lower[index] = SCIPvarGetLbLocal(probvars[i]);
  }
  for (int i = n; i < n2; ++i) { lower[i] = 0; }
  SCIP_OBJSENSE objsense = SCIPgetObjsense(scip);
  if ( objsense == SCIP_OBJSENSE_MINIMIZE )
  {
    maximization = FALSE;
  }
  else
  {
    maximization = TRUE;
  }


  return SCIP_OKAY;
}


/* lattice reduction */
void Reduce(
  mat_ZZ Aext,
  mat_ZZ &Q,
  vec_ZZ &x0
)
{
  int i,j;
  ZZ determ;
  mat_ZZ L, U, X1, Atrans;
  ZZ N1=to_ZZ(10000000);
  ZZ N2=to_ZZ(1000000000);

  int m = Aext.NumRows();
  int n = Aext.NumCols() - 1;

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
    L[n][n+i+1] = to_ZZ(-Aext[i][n])*(N2);

  // lattice reduction
  LLL(determ, L, U, 99, 100, 0);

  // check for successful reduction
  if (L[n-m][n] != N1 && L[n-m][n] != -N1)
  {
    cout << "N1 does not appear in position [n-m, n]" << endl;
    cout << "L[n-m][n] = " << L[n-m][n] << endl;
    cout << "N1 = " << N1 << endl;
    throw  bad_function_call();
  }

  // get AHL basis
  Q.SetDims(n,n-m);
  for (i=0;i<n-m; i++){
    for (j=0;j<n;j++)
      Q[j][i]=L[i][j];
   }
   check_basis(Aext, Q);

  x0.SetLength(n);
  for (i=0; i<n; i++)
    x0[i] = L[n-m][i];


}

/* get new variable bounds */
SCIP_RETCODE get_new_varbounds(
  mat_ZZ basis,
  vector<double> lhs,
  vector<double> rhs,
  vector<int> &upper_bounds,
  vector<int> &lower_bounds
)
{
  int m = basis.NumRows();
  int n = basis.NumCols();

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
    SCIPinfoMessage(scip, NULL, "getting upper bound for var %d\n", k);
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
    if ( SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED )
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
vec_ZZ get_new_objfun(mat_ZZ basis, vec_ZZ x0, vector<double> c)
{
  int n = basis.NumRows();
  int m = n - basis.NumCols();

  vec_ZZ obj;
  obj.SetLength(n-m+1);

  for (int i=0; i<(n-m); i++)
  {
      obj[i] = 0;
      for (int j=0; j<n; j++)
        obj[i] += to_ZZ(c[j])*basis[j][i];
  }

  obj[n-m]=0;
  for (int j=0; j<n; j++)
    obj[n-m]+=to_ZZ(c[j])*x0[j];

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
  string filename = dir + "/ahl_" + file;
  return filename;
}

/* print reformulation */
void print_reformulation(
  SCIP* scip,
  const char *filename,
  mat_ZZ basis,
  vec_ZZ x0,
  vector<double> upper,
  vector<double> lower,
  vector<double> objfun,
  bool maximization
)
{
  int i, j, k = 0;
  ofstream output_file(filename);

  int n = basis.NumRows();
  int m = n - basis.NumCols();

  /* get new objective function */
  vec_ZZ obj = get_new_objfun(basis, x0, objfun);

  /* get bounds of new variables */
  vector<double> lhs, rhs;
  for (j = 0; j < n; j++)
  {
    if (lower[j] < -1e10) { lhs.push_back( -SCIPinfinity(scip) ); }
    else { lhs.push_back( lower[j] - conv<double>(x0[j]) ); }
    if (upper[j] > 1e10) { rhs.push_back( SCIPinfinity(scip) ); }
    else { rhs.push_back( upper[j] - conv<double>(x0[j]) ); }
  }
  vector<int> newupper, newlower;
  SCIPinfoMessage(scip, NULL, "   calculating bounds for new variables.\n");
  SCIP_RETCODE retcode = get_new_varbounds(basis, lhs, rhs, newupper, newlower);


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
    if (lower[j] < -1e10) continue;
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
    if (upper[j] > 1e10) continue;
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
    output_file << newlower[i] << "<= k" << i+1 << "<=" << newupper[i] << "\n";
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

  mat_ZZ Aext;
  bool maximization;
  vector<double> objfun(10000, 0.0);
  vector<double> upper(10000, SCIPinfinity(scip));
  vector<double> lower(10000, -SCIPinfinity(scip));
  SCIP_CALL( GetInstanceData(scip, Aext, upper, lower, objfun, maximization) );

  int m = Aext.NumRows();
  int n = Aext.NumCols() - 1;
  SCIPinfoMessage(scip, NULL, "Shape of Aext: (%d,%d) \n", m,n+1);

  // Save matrix //
  ofstream output_file("contraint_matrix.txt");
  for (int i=0;i<m;i++)
  {
    for (int j=0;j<n+1;j++)
      output_file << Aext[i][j] << " ";
    output_file << "\n";
  }
  for (int j=0;j<n;j++)
    output_file << lower[j] << " ";
  output_file << "\n";
  for (int j=0;j<n;j++)
    output_file << upper[j] << " ";
  output_file << "\n";
  output_file.close();


  // get kernel basis //
  mat_ZZ Q; vec_ZZ x0;
  SCIPinfoMessage(scip, NULL, "    reducing kernel basis\n");
  Reduce(Aext, Q, x0);

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

  // print the reformulated problem //
  string filename = get_new_filename(instancepath);
  print_reformulation(scip, filename.c_str() , Q, x0, upper, lower, objfun, maximization);


  SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1) );
  SCIP_CALL( SCIPinterruptSolve(scip)	);
  return SCIP_OKAY;
}
