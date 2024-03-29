#include <iostream>
#include <fstream>
#include <string>

// SCIP
#include "objscip/objscip.h"
#include "objscip/objconshdlr.h"
#include "objscip/objscipdefplugins.h"

#include "Eventhdlr.hpp"

using namespace scip;



/** reads parameters */
static
SCIP_RETCODE readParams(
   SCIP*                      scip,               /**< SCIP data structure */
   const char*                filename            /**< parameter file name, or NULL */
   )
{
   if( filename != NULL )
   {
      if( SCIPfileExists(filename) )
      {
         std::cout << "reading parameter file <" << filename << ">" << std::endl;
         SCIP_CALL( SCIPreadParams(scip, filename) );
      }
      else
         std::cout << "parameter file <" << filename << "> not found - using default parameters" << std::endl;
   }
   else if( SCIPfileExists("scipmip.set") )
   {
      std::cout << "reading parameter file <scipmip.set>" << std::endl;
      SCIP_CALL( SCIPreadParams(scip, "scipmip.set") );
   }

   return SCIP_OKAY;
}

/** starts SCIP */
static
SCIP_RETCODE fromCommandLine(
   SCIP*                      scip,               /**< SCIP data structure */
   const char*                filename            /**< input file name */
   )
{
   /********************
    * Problem Creation *
    ********************/

   std::cout << std::endl << "read problem <" << filename << ">" << std::endl;
   std::cout << "============" << std::endl << std::endl;
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );

   /*******************
    *    Pre-solve    *
    *******************/
   SCIP_CALL( SCIPpresolve(scip) );
   int m = SCIPgetNConss(scip);
   int n = SCIPgetNVars(scip);
   std::cout << "Number of constraints: " << m << "\n";
   std::cout << "Number of variables: " << n << "\n\n";


   /*******************
    * Problem Solving *
    *******************/

   SCIP_CALL( SCIPsolve(scip) );
   std::cout << std::endl << "FINISHED" << std::endl;
   std::cout << "================" << std::endl << std::endl;
   // SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );


   /**************
    * Statistics *
    **************/

   // std::cout << std::endl << "Statistics" << std::endl;
   // std::cout << "==========" << std::endl << std::endl;
   //
   // SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}


/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runSCIP(
   int                        argc,               /**< number of shell parameters */
   char**                     argv                /**< array with shell parameters */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(scip, NULL);
   std::cout << std::endl;

   /* include the right event handler */
   std::string strategy = argv[2];
   Eventhdlr_AHL EventhdlrAHL(scip, argv[1]);
   Eventhdlr_Pataki EventhdlrPataki(scip, argv[1]);
   if (strategy == "ahl")
   {
      SCIP_CALL( SCIPincludeObjEventhdlr(scip, &EventhdlrAHL, FALSE) );
   }
   else if (strategy == "ahl_diag")
   {
      EventhdlrAHL.diag = true;
      SCIP_CALL( SCIPincludeObjEventhdlr(scip, &EventhdlrAHL, FALSE) );
   }
   else if (strategy == "ahl_poor")
   {
      EventhdlrAHL.highquality = false;
      SCIP_CALL( SCIPincludeObjEventhdlr(scip, &EventhdlrAHL, FALSE) );
   }
   else if (strategy == "kp")
   {
      SCIP_CALL( SCIPincludeObjEventhdlr(scip, &EventhdlrPataki, FALSE) );
   }
   else
   {
      std::cout << "Invalid strategy. Choose among {ahl, ahl_diag, ahl_poor, kp}" << std::endl;
      return SCIP_LPERROR;
   }

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /**************
    * Parameters *
    **************/

   if( argc >= 4 )
   {
      SCIP_CALL( readParams(scip, argv[3]) );
   }
   else
   {
      SCIP_CALL( readParams(scip, NULL) );
   }
   /*CHECK_OKAY( SCIPwriteParams(scip, "scipmip.set", TRUE) );*/


   /**************
    * Start SCIP *
    **************/

   SCIP_CALL( fromCommandLine(scip, argv[1]) );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}


/** main method starting SCIP */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   /* print usage */
   std::cout << " usage: project9 filename strategy settingsfile" << std::endl;

   /* check usage */
   if( argc < 2 )
   {
     std::cout << " No file provided. Please enter a filename. " << std::endl;
     return 1;
   }
   if( argc < 3 )
   {
     std::cout << " No strategy. Please choose a reformulation strategy. " << std::endl;
     return 1;
   }

   /* run SCIP */
   SCIP_RETCODE retcode;
   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
