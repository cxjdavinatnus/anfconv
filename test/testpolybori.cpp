#include <vector>
#include <list>
#include <map>
#include <set>
#include "cryptominisat5/solvertypesmini.h"
using CMSat::lbool;
using CMSat::Lit;
using std::map;
using std::vector;

int main()
{
    BoolePolyRing ring(1000);

    BoolePolynomial poly(false, ring);
    poly += BooleMonomial(ring);
    poly += BooleVariable(5, ring);
    poly += BooleVariable(7, ring);
    poly += BooleVariable(9, ring);
    BooleMonomial m(ring);
    m *= BooleVariable(5, ring);
    m *= BooleVariable(9, ring);
    poly += m;

    m = BooleMonomial(ring);
    m *= BooleVariable(7, ring);
    m *= BooleVariable(11, ring);
    m *= BooleVariable(9, ring);
    std::cout << "Monomial: " << m << std::endl;
    std::cout << "first var: " << m.firstVariable() << std::endl;
    std::cout << "first var's index: " << m.firstVariable().index() << std::endl;
    std::cout << "first index: " << m.firstIndex() << std::endl;
    poly += m;
    std::cout << "poly: " << poly << std::endl;
    std::cout << "poly deg: " << poly.deg() << std::endl;
    std::cout << "term size: " << poly.terms().size() << std::endl;
    std::cout << "first term: " << poly.terms()[0] << std::endl;
    std::cout << "second term: " << poly.terms()[1] << std::endl;

    std::cout << "Testing evaluation:" << std::endl;
    std::cout << "Poly: " << poly << std::endl;
    std::cout << "diagram: " << poly.diagram() << std::endl;

    BoolePolynomial poly2 = poly;

    vector<char> val(12, 0);
    val[5] = 1;
    while (!poly.isConstant()) {
        BoolePolynomial::navigator nav = poly.navigation();
        size_t var = *nav;
        std::cout << "var: " << var << " is equivalent to :" << val[var] << std::endl;
        if (val[var] == 0) {
            poly = BoolePolynomial(nav.thenBranch(), ring);
        } else {
            poly = BoolePolynomial(nav.thenBranch(),ring) + BoolePolynomial(nav.elseBranch(), ring);
        }
    }
    std::cout << "final poly:" << poly << std::endl;


    std::cout << "before change: " << poly2 << std::endl;
    size_t var_to_set = 5;
    bool toSet = true;
    poly2.set().change(5);
    std::cout << "after change: " << poly2 << std::endl;

    exit(-1);
    //anf.readFile("myanf.anf");

    /*BooleMonomial m1(5);
    BooleMonomial m2(2);
    BooleMonomial x = m1*m2*m2;

    std::cout << "Monom: " << x << std::endl;

    BoolePolynomial eq;
    eq += m1;
    eq += m2*(m1*m2);
    eq += x;
    std::cout << "Eq: " << eq << std::endl;

    eq += m1;
    std::cout << "Eq: " << eq << std::endl;*/

    /*Replacer repl;
    repl.newVar(0);
    repl.newVar(1);
    repl.newVar(2);
    repl.newVar(3);

    repl.setReplace(0, Lit(1, true));
    repl.setReplace(2, Lit(3, true));
    repl.setReplace(1, Lit(2, true));
    std::cout << repl << std::endl;
    repl.setValue(0, true);
    std::cout << repl << std::endl;*/

    /*ANF anf;
    BoolePolynomial eq;
    eq += BoolePolynomial(BooleMonomial(3));
    eq += BoolePolynomial(BooleMonomial(4));
    std::cout << "eq:" << eq << std::endl;
    anf.addBoolePolynomial(eq);

    BoolePolynomial eq3;
    eq3 += BoolePolynomial(BooleMonomial(2));
    eq3 += BoolePolynomial(BooleMonomial(4));
    eq3 += BoolePolynomial(true);
    anf.addBoolePolynomial(eq3);

    BoolePolynomial eq4;
    eq4 += BoolePolynomial(BooleMonomial(1));
    eq4 += BoolePolynomial(BooleMonomial(3));
    anf.addBoolePolynomial(eq4);

    std::cout << anf << std::endl;*/
}
