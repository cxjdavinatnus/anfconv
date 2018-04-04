/*****************************************************************************
ANFConv -- Copyright (c) 2011 Mate Soos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "anf.h"
#include "cnf.h"
#include "xlsimplifier.h"
#include <fstream>
#include <sys/wait.h>
#include "time_mem.h"
#include <fstream>
#include "linearizer.h"
#include "satsolve.h"
#include "simplifybysat.h"
#include "replacer.h"
#include <boost/scoped_ptr.hpp>

//#define DEBUG_EXTRACTION
using std::cout;
using std::cerr;
using std::endl;

std::string fileToRead;
std::string fileToWriteANF;
std::string fileToWriteCNF;
std::string programName;

//Writing options
bool writeANF;
bool writeCNF;

//Solving options
bool doSolveSAT; //Solve using CryptoMiniSat
bool doSolveGauss; //Solve using Gaussian elimination
bool doLinearize;

//Dependency
bool doCalcDep;
int doMinimise;

//How to simplify
bool doXLSimplify;
bool doANFSimplify = true; //NOTE: Things WILL break if this is not enabled
bool doSATSimplify;

//Parameters
ConfigData config;
uint64_t numConfl;
vector<string> extractString;

void parseOptions(int argc, char *argv[])
{
    // Declare the supported options.
    po::options_description generalOptions("Allowed options");
    generalOptions.add_options()
    ("help,h", "produce help message")
    ("read,r", po::value<std::string>(&fileToRead)
        , "Read ANF from this file")
    ("anfwrite,a", po::value<std::string>(&fileToWriteANF)
        , "Write ANF output to file")
    ("cnfwrite,c", po::value<std::string>(&fileToWriteCNF)
        , "Write CNF output to file")
    ("solvesat,s", po::bool_switch(&doSolveSAT)
        , "Solve with SAT solver")
    ("solvegauss,g", po::bool_switch(&doSolveGauss)
        , "Solve with Gaussian elimination")
    ("depend", po::bool_switch(&doCalcDep)
        , "Calculate and print dependecy of variables")
    ("karn", po::value<int>(&config.useKarn)->default_value(config.useKarn)
        , "Use Karnaugh table minimisation")
    ("program,p", po::value<std::string>(&programName)->default_value("/usr/local/bin/cryptominisat")
        , "SAT solver to use with full path")
    ("verbosity,v", po::value<int>(&config.verbosity)->default_value(1)
        , "Verbosity setting (0 = silent)")
    ("lin,l", po::bool_switch(&doLinearize)
        , "Solve by linearization")
    ("dump,d", po::bool_switch(&config.writePNG)
         , "Dump XL's and linearization's matrixes as PNG files")
    ("satsimp", po::bool_switch(&doSATSimplify)
        , "Simplify using SAT")
    ("xlsimp", po::bool_switch(&doXLSimplify)
        , "Simplify using XL")
    ("confl", po::value<uint64_t>(&numConfl)->default_value(20000)
        , "Conflict limit for built-in SAT solver")
    ("cutnum", po::value<uint32_t>(&config.cutNum)->default_value(config.cutNum)
        , "Cutting number when not using XOR clauses")
    ("extract,e", po::value<vector<string> >(&extractString)->multitoken()
        , "Extract the values of these variables as binary string from final ANF. \
             For example, '10-20' extracts x10...x20")
    ("minim", po::value<int>(&doMinimise)->default_value(0)
        , "Minimise (renumber) ANF before attempting dependency calculation")

    ;

    po::positional_options_description p;
    p.add("read", 1);

    po::variables_map vm;
    po::options_description cmdline_options;
    cmdline_options
    .add(generalOptions)
    ;

    try {
        po::store(
            po::command_line_parser(argc, argv)
            .options(cmdline_options)
            .positional(p).run()
            , vm
        );

        if (vm.count("help"))
        {
            cout << generalOptions << endl;
            exit(0);
        }

        po::notify(vm);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::unknown_option> >& c
    ) {
        cout << "Some option you gave was wrong. Please give '--help' to get help" << endl;
        cout << "Unkown option: " << c.what() << endl;
        exit(-1);
    } catch (boost::bad_any_cast &e) {
        cerr << e.what() << endl;
        exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_option_value> > what
    ) {
        cerr
        << "Invalid value '" << what.what() << "'"
        << " given to option '" << what.get_option_name() << "'"
        << endl;

        exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> > what
    ) {
        cerr
        << "Error: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> > what
    ) {
        cerr
        << "You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        exit(-1);
    }

    //Writing methods
    if (vm.count("anfwrite")) {
        writeANF = true;
    }
    if (vm.count("cnfwrite")) {
        writeCNF = true;
    }

    if (config.cutNum < 3 || config.cutNum > 10) {
        cout
        << "ERROR! For sanity, cutting number must be between 3 and 10"
        << endl;

        exit(-1);
    }

    if (!vm.count("read")) {
        cout
        << "You must give an ANF file to read in"
        << endl;

        exit(-1);
    }
}

int main(int argc, char *argv[])
{
    parseOptions(argc, argv);
    double myTime = cpuTime();

    //test();

    //ANF generation
    BoolePolyRing* ring = new BoolePolyRing(1);
    ANF* anf = new ANF(ring, config);
    const size_t ringSize = anf->readFile(fileToRead, false);
    cout << "--> Needed ring size is " << ringSize+1 << endl;
    delete anf;
    delete ring;

    ring = new BoolePolyRing(ringSize+1);
    anf = new ANF(ring, config);
    const size_t ringSize2 = anf->readFile(fileToRead, true);
    assert(ringSize == ringSize2);

    if (config.verbosity >= 1) {
        anf->printStats(config.verbosity);
    }

    //Simplification(s), recursively
    ANF origANF(*anf);
    if (doANFSimplify || doSATSimplify|| doXLSimplify) {
        if (config.verbosity>=1)
            cout << "Simple simplifying of ANF..." << std::flush;

        bool changed = true;
        uint32_t numIters = 0;
        while (changed) {
            changed = false;
            if (doANFSimplify) {
                anf->simplify(config.anfReplaceVars, true);
                //anf->simplify(config.anfReplaceVars, true);

                if (config.verbosity>=1) {
                    cout << "c Done. New stats:" << endl;
                    if (config.verbosity >= 1) {
                        cout << "c Time to read&simplify ANF: "
                        << std::fixed << std::setprecision(2) << (cpuTime() - myTime)
                        << endl;
                        anf->printStats(config.verbosity);
                    }
                }
            }

            if (doSATSimplify) {
                cout << "c Simplifying using SAT solver" << endl;
                SimplifyBySat simpsat(*anf, origANF, config);
                changed |= simpsat.simplify(numConfl);
                cout << "c Done simplifying with SAT solver." << endl;
            }

            //XL simplification
            if (doXLSimplify) {
                XLSimplifier xl(*anf, 2);
                xl.simplify(*anf, numIters, config);
                cout << "Simplifying ANF with XL..." << std::flush;
                changed |= xl.simplify(*anf, numIters, config);
                cout << "Done." << endl;
            }
            numIters++;
        }
    }

    //Writing simplified ANF
    if (writeANF) {
        std::ofstream ofs;
        ofs.open(fileToWriteANF.c_str());
        if (!ofs) {
            std::cerr
            << "Error opening file \"" << fileToWriteANF << "\" for writing"
            << endl;
            exit(-1);
        }

        ANF* newanf = NULL;
        if (doMinimise) {
            anf->preferLowVars();
            vector<uint32_t> oldToNew;
            vector<uint32_t> newToOld;
             newanf = anf->minimise(oldToNew, newToOld);
        } else {
            newanf = anf;
        }

        //Write to file
        ofs << *newanf << endl;
        ofs.close();

        if (doMinimise) {
            delete newanf;
        }
    }

    //Linearization-based solving
    if (doLinearize) {
        //Set up fullmatrix
        Linearizer linMatrix(*anf, config);

        //Create matrix
        linMatrix.calcFullEqMap();

        //Dump matrix in unechelonized form
        linMatrix.dumpMatrix("matrix.data");
        if (config.writePNG) {
            string filename = "matrixorig.png";
            cout
            << "Writing linearized matrix before echelonization to file '"
            << filename
            << endl;

            linMatrix.writePNG("matrixorig.png");
        }

        //Echelonize matrix
        linMatrix.echelonizeMatrix();

        //Dump echelonized matrix
        if (config.writePNG) {
            string filename = "matrixfinal.png";
            cout
            << "Writing linearized matrix after echelonization to file '"
            << filename
            << endl;
            linMatrix.writePNG(filename);
        }
        linMatrix.getSolutionFromMatrix();
        linMatrix.writeMatrixInfosToFile();

        if (config.verbosity >= 1) {
            cout << "Number of full equations for internal variables calculated: "
            << linMatrix.getFullVarEqsCalculated() << endl;
        }
    }

    //We need to extract data
    if (!extractString.empty()) {

        if (!anf->getOK()) {
            //If UNSAT, there is no solution to extract
            cout << "UNSAT, so no solution to extract" << endl;
        } else {
            //Go through each piece of data that needs extraction
            for(vector<string>::const_iterator
                it = extractString.begin(), end = extractString.end()
                ; it != end
                ; it++
            ) {
                anf->extractVariables(*it);
            }
        }
    }

    //Solve with gaussian elimination
    if (doSolveGauss) {
        Linearizer full2(*anf, config);
        full2.calcFullEqMap();
    }

    if (doCalcDep) {
        assert(config.anfReplaceVars);
        ANF* newanf = NULL;
        if (doMinimise) {
            anf->preferLowVars();
            vector<uint32_t> oldToNew;
            vector<uint32_t> newToOld;
             newanf = anf->minimise(oldToNew, newToOld);
        } else {
            newanf = anf;
        }
        DependencyCalc calcDep(*newanf, config);
        vector<uint32_t> var_order = calcDep.vars_to_compute_for_all_outputs();
        calcDep.calc_polys_in_order(var_order);
        calcDep.get_output_poly_depend();
        newanf->printIndependentVars(config.verbosity);
        //cout << "map of 430: " << oldToNew[430] << endl;

        if (doMinimise) {
            delete &(newanf->getRing());
            delete newanf;
        }
    }

    //CNF conversion and solving
    if (doSolveSAT || writeCNF) {
        myTime = cpuTime();
        CNF cnf(*anf, config);
        cnf.init();
        cnf.addAllEquations();
        cout << "c ---- CNF stats -----" << endl;
        cout << "c CNF conversion time  : " << (cpuTime() - myTime) << endl;
        cout << "c Clauses              : " << cnf.getNumClauses() << endl;
        cout << "c XorClauses           : " << cnf.getNumXors() << endl;
        cout << "c Added as CNF         : " << cnf.getAddedAsCNF() << endl;
        cout << "c Added as simple ANF  : " << cnf.getAddedAsSimpleANF() << endl;
        cout << "c Added as complex  ANF: " << cnf.getAddedAsComplexANF() << endl;

        if (writeCNF) {
            std::ofstream ofs;
            ofs.open(fileToWriteCNF.c_str());
            if (!ofs) {
                cout << "Error opening file \""
                << fileToWriteCNF
                << "\" for writing" << endl;
                exit(-1);
            }
            ofs << cnf << endl;
            ofs.close();
        }

        if (doSolveSAT) {
            SATSolve solver(config.verbosity, programName);
            vector<lbool> sol = solver.solveCNF(origANF, *anf, cnf);

            if (config.verbosity >= 1) {
                cout << "c CPU time unknown " << endl;
            }

            if (!sol.empty()) {
                for(vector<string>::const_iterator
                    it = extractString.begin(), end = extractString.end()
                    ; it != end
                    ; it++
                ) {
                    anf->extractVariables(*it, &sol);
                }
            }
        }
    }

    return 0;
}
