/*****************************************************************************
anfconv
Copyright (C) 2016  Security Research Labs

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

#ifndef _ITERATIVE_SOLVER_H_
#define _ITERATIVE_SOLVER_H_

#include "anf.h"
#include "configdata.h"
#include "cnf.h"
#include "cryptominisat5/cryptominisat.h"
#include <string>
using std::string;

class IterativeSolver
{
public:
    IterativeSolver(const ANF& _anf, const ConfigData& _config);
    ~IterativeSolver();

    //Iteratively solve
    void solve();

private:

    //Main constant data
    const ANF& anf;
    const ConfigData config;

    //Claculated data
    CMSat::SATSolver* solver;
    CNF cnf;

    const ClausesAdded all_clauses();
    void addClauses(const ClausesAdded& toAdd);
    void addNewVarsIfNeedBe(const vector<Lit>& lits);
    void addNewVarsIfNeedBe(const vector<uint32_t>& vars);

};

#endif //_ITERATIVE_SOLVER_H_
