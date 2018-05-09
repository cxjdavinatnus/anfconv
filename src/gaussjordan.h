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

#ifndef __GAUSSJORDAN_H__
#define __GAUSSJORDAN_H__

#include "anf.h"
#include "m4ri.h"
#include "configdata.h"
#include "time_mem.h"

#include <map>
using std::map;

class GaussJordan
{
    public:
        // upToDeg = -1 means use all polynomials
        GaussJordan(ANF& anf, const ConfigData& config, const int upToDeg = -1) :
                ring(anf.getRing()), config(config), nextVar(0) {
            printf("RUNNING GAUSS\n");
            buildMaps(anf, upToDeg);
            mat = mzd_init(anf.getEqs().size(), nextVar+1);
            assert(mzd_is_zero(mat));

            size_t num = 0;
            for (const BoolePolynomial& poly : anf.getEqs()) {
                if (upToDeg == -1 || poly.deg() <= upToDeg) {
                    for(const BooleMonomial& mono : poly) {
                        // Process non-empty monomials (aka non-constant)
                        if (mono.deg() != 0) {
                            map<BooleMonomial, uint32_t>::const_iterator it3 = monomMap.find(mono);
                            assert(it3 != monomMap.end());
                            uint32_t intVar = it3->second;
                            mzd_write_bit(mat, num, intVar, true);
                        }
                    }
                    // Constant goes here
                    mzd_write_bit(mat, num, nextVar, poly.hasConstantPart());
                    num++;
                }
            }
        }

        ~GaussJordan() {
            mzd_free(mat);
        }

        const mzd_t* getMatrix() const {
            return mat;
        }

        void run(ANF& addNewTruthsHere) {
            // Perform Gauss Jordan elimination
            double myTime = cpuTime();

            vector<BoolePolynomial> newTruths;
            findTruths(newTruths);

            //Add new truths, and simplify ANF
            if (newTruths.empty()) {
                cout << "No new truths found!" << endl;
            }
            cout << "Found " << newTruths.size() << " new truth(s)" << endl;

            //Add and print new truths found
            for(const BoolePolynomial& poly : newTruths) {
                if (config.verbosity >= 2) {
                    cout << "New truth: " << poly << endl;
                }
                addNewTruthsHere.addBoolePolynomial(poly);
            }

            cout << "c Time to perform Gauss Jordan: " << (cpuTime() - myTime) << endl;
            addNewTruthsHere.printStats(config.verbosity);
        }

    private:
        void buildMaps(const ANF& anf, const int upToDeg) {
            for(const BoolePolynomial& poly : anf.getEqs()) {
                if (upToDeg == -1 || poly.deg() <= upToDeg) {
                    for(const BooleMonomial& mono : poly) {
                        addMonom(mono);
                    }
                }
            }
        }

        void addMonom(const BooleMonomial& mono) {
            map<BooleMonomial, uint32_t>::const_iterator it = monomMap.find(mono);
            if (it == monomMap.end()) {
                monomMap[mono] = nextVar;
                revMonomMap.insert(make_pair(nextVar, mono));
                nextVar++;
            }
        }

        void findTruths(vector<BoolePolynomial>& truths) {
            //matrix is extended with assignement column, so at least 1 col large
            assert(mat->ncols > 0);

            mzd_echelonize(mat, 1);

            vector<int> ones;
            #ifdef DUMP_SIMP_EQ
            std::ofstream simpEQFile("simpEq.txt");
            #endif

            for (int row = mat->nrows-1; row >= 0; row--) {
                ones.clear();
                for (int thisCol = 0; thisCol < mat->ncols-1; thisCol++) {
                    if (mzd_read_bit(mat, row, thisCol)) {
                        ones.push_back(thisCol);
                        if (ones.size() > 2) {
                            break;
                        }
                    }
                }

                // Not trivial truth. Ignore
                if (ones.size() > 2) {
                    continue;
                }

                if (ones.size() == 0) {
                    if (mzd_read_bit(mat, row, mat->ncols-1) == 1) {
                        // Row is "1 = 0", UNSAT
                        truths.push_back(BoolePolynomial(true, ring));
                        return;
                    } else {
                        // Row is "0 = 0", skip row
                        continue;
                    }
                }

                #ifdef DUMP_SIMP_EQ
                {
                    BoolePolynomial eq(mzd_read_bit(mat, row, mat->ncols-1), ring);
                    for (int thisCol = 0; thisCol < mat->ncols-1; thisCol++) {
                        if (mzd_read_bit(mat, row, thisCol)) {
                            eq += revMonomMap.find(thisCol)->second;
                        }
                    }
                    simpEQFile << eq;
                }
                #endif

                BoolePolynomial eq(mzd_read_bit(mat, row, mat->ncols-1), ring);
                uint32_t deg0 = 0;
                uint32_t deg1 = 0;
                assert(1 <= ones.size() && ones.size() <= 2);

                map<uint32_t, BooleMonomial>::const_iterator it = revMonomMap.find(ones[0]);
                assert(it != revMonomMap.end());
                deg0 = it->second.deg();
                eq += BoolePolynomial(it->second);

                if (ones.size() == 2) {
                    map<uint32_t, BooleMonomial>::const_iterator it = revMonomMap.find(ones[1]);
                    assert(it != revMonomMap.end());
                    deg1 = it->second.deg();
                    eq += BoolePolynomial(it->second);
                }

                //
                // IGNORE non-plain truths
                //
                // E.g. stuff like "a*b + c = 1"
                if (ones.size() == 2 && (deg0 > 1 || deg1 > 1)) {
                    continue;
                }

                // E.g. stuff like "a*b = 0"
                // Note: "a*b + 1 = 0" is useful because simplification will yield a = b = 1
                if (ones.size() == 1 && deg0 > 1 && !eq.hasConstantPart()) {
                    continue;
                }

                truths.push_back(eq);
            }
        }

        //Parity matrix
        mzd_t *mat;
        BoolePolyRing& ring;
        const ConfigData& config;

        //map: outside monom -> inside var
        std::map<BooleMonomial, uint32_t> monomMap;
        //map: inside var -> outside monom
        std::map<uint32_t, BooleMonomial> revMonomMap;
        uint32_t nextVar;
};

#endif //__GAUSSJORDAN_H__
