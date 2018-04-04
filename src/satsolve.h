/*****************************************************************************
anfconv
Copyright (C) 2016  Security Research Labs

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
******************************************************************************/

#ifndef SATSOLVE__H
#define SATSOLVE__H

#include "anf.h"
#include "cnf.h"

#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <string>
using std::string;

inline void testSolution(const ANF& anf, const vector<lbool>& solution)
{
    bool goodSol = anf.evaluate(solution);
    if (!goodSol) {
        std::cout << "ERROR! Solution found is incorrect!" << std::endl;
        exit(-1);
    } else {
        //if (verbosity >= 1)
            std::cout << "c Solution checked, found correct." << std::endl;
    }
}

inline void printSolution(const vector<lbool>& solution)
{
    size_t num = 0;
    std::stringstream toWrite;
    toWrite << "v ";
    for(vector<lbool>::const_iterator
        it = solution.begin(), end = solution.end()
        ; it != end
        ; it++, num++
    ) {
        if (*it != l_Undef) {
            toWrite << ((*it == l_True) ? "" : "-") << num << " ";
        }
    }
    std::cout << toWrite.str() << std::endl;
}

class SATSolve
{
    public:
        SATSolve(const int verbosity = 1, std::string solverExecutable = "./cmsat");

        vector<lbool> solveCNF(const ANF& orig_anf, const ANF& anf, const CNF& cnf);
        const vector<lbool>& getSolution() const
        {
            return solution;
        }

        const lbool getSatisfiable() const
        {
            return satisfiable;
        }

    private:
        void createChildProcess();
        string readchild();
        void error(const char *s);

        //Managing processes
        int in[2];
        int out[2];
        int pid; //process ID

        //Solution
        vector<lbool> solution;
        lbool satisfiable;

        //Parameters
        string solverExecutable;
        int verbosity;
};

#endif //SATSOLVE__H
