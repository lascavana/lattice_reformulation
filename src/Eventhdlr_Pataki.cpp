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

#include "Eventhdlr.hpp"
#include "utils.hpp"

using namespace std;
using namespace NTL;
using matrix = vector<vector<int>>;



/* get constraint matrix */
SCIP_RETCODE GetInstanceData(
  SCIP* scip,
  mat_ZZ &consmat,
  vector<int> &lhs,
  vector<int> &rhs,
  vector<int> &upper,
  vector<int> &lower,
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
  consmat.SetDims(m,n); // Aext = [lhs|A|rhs]
  lhs.resize(m); rhs.resize(m);
  SCIPinfoMessage(scip, NULL, "    Aext dim: (%d, %d)\n", m, n);

  /* assign var index to a matrix colum */
  unordered_map<int, int> idx2col;
  SCIP_VAR** allvars = SCIPgetVars(scip);
  for (int i = 0; i < n; ++i)
  {
    	int index = SCIPvarGetIndex(allvars[i]);
      idx2col.insert({index, i});
  }

  /* get contraint data */
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
    lhs[i] = double2int(SCIPconsGetLhs(scip, conss[i], &success));
    rhs[i] = double2int(SCIPconsGetRhs(scip, conss[i], &success));

    /* fill in matrix entries */
    for (int j = 0; j < nvars; ++j)
    {
      if (vals[j] == 0) continue;
      int col = idx2col[ SCIPvarGetIndex(vars[j]) ];
      try
      {
        consmat[i][col] = double2zz(vals[j]);
      }
      catch (const std::invalid_argument& e)
      {
        cout << e.what() << endl;
        return SCIP_ERROR;
      }
    }

    /* release arrays */
    SCIPfreeBufferArray(scip, &vals);
    SCIPfreeBufferArray(scip, &vars);
  }


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
  objfun.resize(n+1); upper.resize(n); lower.resize(n);
  for (int i = 0; i < n; ++i)
  {
    	int col = idx2col[ SCIPvarGetIndex(allvars[i]) ] ;
      objfun[col] = objscale*SCIPvarGetObj(allvars[i]);
      upper[col] = double2int(SCIPvarGetUbLocal(allvars[i]));
      lower[col] = double2int(SCIPvarGetLbLocal(allvars[i]));
  }
  objfun[n] = objscale*SCIPgetTransObjoffset(scip);
  for (int i = n; i < n; ++i) { lower[i] = 0; }


  return SCIP_OKAY;
}

/* extend */
void extend_mat(
  mat_ZZ &A,
  vector<int> &lhs,
  vector<int> &rhs,
  vector<int> upper,
  vector<int> lower
)
{
  int n = A.NumCols();
  int m = A.NumRows();

  int ncons = m;
  for (int j=0; j<n; j++)
  {
    if (lower[j] > -1e9 || upper[j] < 1e9) ncons++;
  }

  /* initialize matrix */
  mat_ZZ Aext; Aext.SetDims(ncons, n);
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<m; i++)
      Aext[i][j] = A[i][j];
    for (int i=m; i<ncons; i++)
      Aext[i][j] = to_ZZ(0);
  }

  /* fill in identity columns */
  int idx = m;
  for (int j=0; j<n; j++)
  {
    if (lower[j] > -1e20 || upper[j] < 1e20)
    {
      Aext[idx][j] = to_ZZ(1);
      lhs.push_back(lower[j]);
      rhs.push_back(upper[j]);
      idx++;
    }
  }
  assert(idx==ncons);

  A = Aext;
}


/* reduction */
mat_ZZ reduce_pataki(
  const mat_ZZ &A,
  mat_ZZ &U
)
{
  mat_ZZ Atrans;
  transpose(Atrans, A);

  /* reduce */
  ZZ determ;
  LLL(determ, Atrans, U, 99, 100, 0);

  /* transpose matrix */
  mat_ZZ Ared;
  transpose(Ared, Atrans);

  return Ared;
}

void transform_obj(
  vector<double> &objfun,
  const mat_ZZ &U
)
{
  int n = objfun.size() - 1;
  vector<double> old_objfun = objfun;

  for (int j=0; j<n; j++)
  {
    objfun[j] = 0.0;
    for (int i=0; i<n; i++)
    {
      objfun[j] += old_objfun[i]*conv<double>(U[j][i]);
    }   
  }
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
  string filename = dir + "/pat_" + file + ".lp";
   
  return filename;
}

void print_pataki(
  SCIP* scip,
  const char *filename,
  matrix A,
  vector<int> lhs,
  vector<int> rhs,
  vector<double> objfun,
  bool maximization
)
{
  ofstream output_file(filename);

  int n = A[0].size();
  int m = A.size();

  /* get new variable bounds */
  vector<int> newupper, newlower;
  SCIP_RETCODE retcode = get_new_varbounds(A, lhs, rhs, newupper, newlower);

  /* write objective function */
  if (maximization) { output_file << "maximize "; }
  else { output_file << "minimize "; }
  for (int j=0; j<n; j++)
  {
    if (objfun[j] != 0)
    {
      if (objfun[j]>0)
      {
        output_file << "+"<<objfun[j]<<" k"<<j+1<< " ";
      }
      else
      {
        output_file << objfun[j]<<" k"<<j+1<< " ";
      }
    }
  }
  if (objfun[n] > 0) output_file << "+ " << objfun[n];
  if (objfun[n] < 0) output_file << objfun[n];
  output_file << "\n";

  int conss_counter = 0;
  output_file << "subject to" << "\n";

  /* write constraints */
  for (int i=0; i<m; i++)
  {
    if (lhs[i] < -1e9) continue;

    conss_counter++;
    output_file << "C" << conss_counter << ": ";
    for (int j=0; j<n; j++)
    {
      if (A[i][j] != 0)
      {
        if (A[i][j] > 0)
        {
          output_file << "+"<<A[i][j]<<" k"<<j+1<<" ";
        }
        else
        {
          output_file << A[i][j]<<" k"<<j+1<<" ";
        }
      }
    }
    output_file << " >= ";
    output_file << lhs[i] << "\n";

  }
  for (int i=0; i<m; i++)
  {
    if (rhs[i] > 1e9) continue;

    conss_counter++;
    output_file << "C" << conss_counter << ": ";
    for (int j=0; j<n; j++)
    {
      if (A[i][j] != 0)
      {
        if (A[i][j] > 0)
        {
          output_file << "+"<<A[i][j]<<" k"<<j+1<<" ";
        }
        else
        {
          output_file << A[i][j]<<" k"<<j+1<<" ";
        }
      }
    }
    output_file << " <= ";
    output_file << rhs[i] << "\n";

  }

  /* write variable bounds */
  output_file << "Bounds\n";
  for (int i=0; i<n; i++)
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
  int k = 0;
  output_file << "Generals " << "\n";
  for (int i=0; i<n; i++)
  {
    if (k>20){k=0; output_file << "\n";}
    output_file << "k" << i+1 << " ";
    k++ ;
  }

  print_for_ls("LS.txt", A, lhs, rhs, newupper, newlower);
}


/** destructor of event handler to free user data (called when SCIP is exiting) */
SCIP_DECL_EVENTFREE(Eventhdlr_Pataki::scip_free)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */

SCIP_DECL_EVENTINIT(Eventhdlr_Pataki::scip_init)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */

SCIP_DECL_EVENTEXIT(Eventhdlr_Pataki::scip_exit)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The event handler may use this call to initialize its branch and bound specific data.
 *
 */

SCIP_DECL_EVENTINITSOL(Eventhdlr_Pataki::scip_initsol)
{
   SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The event handler should use this call to clean up its branch and bound data.
 */

SCIP_DECL_EVENTEXITSOL(Eventhdlr_Pataki::scip_exitsol)
{
   return SCIP_OKAY;
}


/** frees specific constraint data */

SCIP_DECL_EVENTDELETE(Eventhdlr_Pataki::scip_delete)
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
SCIP_DECL_EVENTEXEC(Eventhdlr_Pataki::scip_exec)
{  /*lint --e{715}*/

  /* read problem data */
  mat_ZZ A;
  bool maximization;
  vector<double> objfun(10000, 0.0);
  vector<int> lhs(10000, SCIPinfinity(scip));
  vector<int> rhs(10000, SCIPinfinity(scip));
  vector<int> upper(10000, SCIPinfinity(scip));
  vector<int> lower(10000, -SCIPinfinity(scip));
  SCIP_CALL( GetInstanceData(scip, A, lhs, rhs, upper, lower, objfun, maximization) );
  int m = A.NumRows();
  int n = A.NumCols();

  /* extend constraint mat to include variable bounds */
  extend_mat(A, lhs, rhs, upper, lower);
  m = A.NumRows();

  /* get reduced A */
  mat_ZZ U;
  auto start = chrono::high_resolution_clock::now();
  mat_ZZ Apat = reduce_pataki(A, U);
  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
  SCIPinfoMessage(scip, NULL, "    took %lld ms to reduce \n", duration.count());

  /* get new objective function */
  transform_obj(objfun, U);

  /* convert to int type */
  matrix basis(m, vector<int>(n, 0));
  for (int i=0;i<m;i++)
  {
    for (int j=0;j<n;j++)
      basis[i][j] = conv<int>(Apat[i][j]);

  }

  /* print the reformulated problem */
  string filename = get_new_filename(instancepath);
  print_pataki(scip, filename.c_str(), basis, lhs, rhs, objfun, maximization);


  SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1) );
  SCIP_CALL( SCIPinterruptSolve(scip)	);
  return SCIP_OKAY;
}
