/*****************************************************************************

    test/test4.cpp

    This file is a part of the Arageli library.

    Copyright (C) 2005--2007 Sergey S. Lyalin
    University of Nizhni Novgorod, Russia

    The Arageli Library is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License version 2
    as published by the Free Software Foundation.

    The Arageli Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

    We are also open for dual licensing for the whole library or
    for its particular part. If you are interested to get the library
    in this way, i.e. not under the GNU General Public License,
    please contact Arageli Support Service support.arageli@gmail.com.

*****************************************************************************/

#include "stdafx.hpp"


//#include <arageli/simplex_method.hpp>
//#include <arageli/rand.hpp>
//#include <arageli/ctrl_slog.hpp>
//#include <arageli/ctrl_latexlog.hpp>
//#include <arageli/residue.hpp>
//#include <arageli/algebraic.hpp>
//#include <arageli/bareiss.hpp>
//
//#include <arageli/big_float.hpp>
//#include <arageli/intcount_barvinok.hpp>
//
//#include <arageli/cone.hpp>
//#include <arageli/polyhedron.hpp>
//#include <arageli/sideset.hpp>

#include "test1.hpp"

#pragma warning (disable : 4503)


using namespace Arageli;
//using namespace std;
using namespace Arageli::simplex_method;
using namespace Arageli::ctrl;
using namespace Arageli::ctrl::simplex_method;
using Arageli::vector;

using std::ostream;
using std::cout;
using std::cin;
using std::ofstream;
using std::ifstream;
using std::cerr;
using std::endl;
using std::string;


//string output_file_name;    // see test4_40
//ofstream output_file;

//#define cout output_file


//struct itflt : public intconvex_triangulation_idler
//{
//    /// WARNING! TEMPORARY MEMBER FOR DEBUGGING.
//    std::ostream& stream () const { return output_file; }
//};


const char* slr_name (solve_linsys_result status)
{
    static const char* names[] =
    {
        "SLR_EMPTY",
        "SLR_UNIQUE",
        "SLR_MULTIPLE"
    };

    return names[status];
}


const char* neir_name (nei_result status)
{
    static const char* names[] =
    {
        "NEIR_EMPTY",
        "NEIR_ALL_VERTS",
        "NEIR_AMONG_VERTS",
        "NEIR_NEW"
    };

    return names[status];
}


struct the_second_test_Shevchenko_ctrl_idler
{
    template <typename A, typename B, typename C>
    void preamble (const A& a, const B& b, const C& c) const {}

    basis_artificial_idler ctrl_for_basis_artificial () const
    { return basis_artificial_idler(); }

    primal_row_iters_idler ctrl_for_primal_row_iters () const
    { return primal_row_iters_idler(); }

    template <typename T, typename X>
    void primal_opt (const T& res, const X& x_opt) const {}

    gomory1_iters_idler ctrl_for_primal_task_gomory1_iters () const
    { return gomory1_iters_idler(); }

    gomory1_iters_idler ctrl_for_dual_task_gomory1_iters () const
    { return gomory1_iters_idler(); }

    dual_col_iters_idler ctrl_for_dual_col_iters () const
    { return dual_col_iters_idler(); }

    template <typename T, typename X>
    void primal_opt_int (const T& res, const X& x_opt_int) const {}

    template
    <
        typename Da, typename Db, typename Dc,
        typename T, typename Bsuba, typename Bsubc
    >
    void dual
    (
        const Da& da, const Db& db, const Dc& dc,
        const T& res_offset,
        const Bsuba& bsuba, const Bsubc& bsubc
    ) const {}

    template <typename T, typename Y>
    void dual_y_opt (const T& res, const Y& y_opt) const {}

    template
    <
        typename Invbsuba, typename Bsubc,
        typename Y, typename U
    >
    void y_to_u_opt
    (
        const Invbsuba& invbsuba, const Bsubc& bsubc,
        const Y& y_opt, const U& u_opt
    ) const {}

    template <typename T, typename Y>
    void dual_y_opt_int (const T& res, const Y& y_opt) const {}

    template
    <
        typename Invbsuba, typename Bsubc,
        typename Y, typename U
    >
    void y_to_u_opt_int
    (
        const Invbsuba& invbsuba, const Bsubc& bsubc,
        const Y& y_opt_int, const U& u_opt_int
    ) const {}

};


// Решается задача целочисленного линейного программирования
// и двойственная к ней симплек методом с отсечением Гомори.
// Ограничения на исходную задачу: в векторе b не должно быть
// отрицательных компонент; и прямая, и двойственная задачи
// должные иметь целочисленное решение.
template
<
    typename T1, bool REFCNT1,
    typename T2, bool REFCNT2,
    typename T3, bool REFCNT3,
    typename Ctrler
>
void the_second_test_Shevchenko
(
    const matrix<T1, REFCNT1>& a,
    const vector<T2, REFCNT2>& b,
    const vector<T3, REFCNT3>& c,
    Ctrler ctrler
)
{
    typedef matrix<T1, false> Matrix;
    typedef vector<T1, false> Vector;
    typedef vector<size_t, false> Basis;

    ctrler.preamble(a, b, c);
    Matrix q;
    row_table_create(a, b, q);
    Basis basis;
    result_kind rk = basis_artificial
        (q, basis, ctrler.ctrl_for_basis_artificial());
    ARAGELI_ASSERT_0(rk == rk_found);

    //Matrix bsuba;
    //fill_submatrix_col(a, basis-1, bsuba);
    //fill_subvector(c, basis-1, bsubc);
    Basis first_basis = basis;    // сохраняем на будущее

    // Решаем прямую задачу прямым строчечным симплекс методом.

    row_table_place_c(c, q);
    row_table_pivot_basis_c(q, basis);
    rk = primal_row_iters(q, basis, ctrler.ctrl_for_primal_row_iters());
    ARAGELI_ASSERT_0(rk == rk_found);

    Vector x_opt;
    row_table_extract_solution(q, basis, x_opt);
    ctrler.primal_opt(q(0, 0), x_opt);

    // Находим решение прямой задачи в целых числах.

    Matrix t;
    Basis nonbasis;
    row_to_col_table(q, basis, t, nonbasis);
    gomory1_iters
        (t, nonbasis, ctrler.ctrl_for_primal_task_gomory1_iters());
    Vector x_opt_int;
    col_table_extract_solution(t, a.ncols(), x_opt_int);
    ctrler.primal_opt_int(t(0, 0), x_opt_int);

    // Строим двойственную задачу.

    Matrix bsuba, da;
    Vector bsubc, db, dc;
    T1 primal_res_offset;
    primal_to_dual_standard_discr
        (a, b, c, first_basis, da, db, dc, primal_res_offset, bsuba, bsubc);
    cout << "***** " << det(bsuba) << " *******";
    ctrler.dual(da, db, dc, primal_res_offset, bsuba, bsubc);

    // Решаем двойственную задачу столбцовым двойственным симплекс методом.

    col_table_create_by_standard(da, db, dc, t, nonbasis);
    t(0, 0) = primal_res_offset;
    rk = dual_col_iters
        (t, nonbasis, ctrler.ctrl_for_dual_col_iters());
    ARAGELI_ASSERT_0(rk == rk_found);
    Vector y_opt;
    col_table_extract_solution(t, da.ncols(), y_opt);
    ctrler.dual_y_opt(-t(0, 0), y_opt);

    Matrix invbsuba = inverse(bsuba);
    Vector u_opt = (y_opt + bsubc)*inverse(bsuba);
    ctrler.y_to_u_opt(invbsuba, bsubc, y_opt, u_opt);

    // Решаем двойственную задачу в целых числах.

    gomory1_iters
        (t, nonbasis, ctrler.ctrl_for_dual_task_gomory1_iters());
    Vector y_opt_int;
    col_table_extract_solution(t, da.ncols(), y_opt_int);
    ctrler.dual_y_opt_int(-t(0, 0), y_opt_int);

    Vector u_opt_int = (y_opt_int + bsubc)*inverse(bsuba);
    ctrler.y_to_u_opt_int(invbsuba, bsubc, y_opt_int, u_opt_int);

    cout << "\n*****************************************\n";
    cout
        << "\nOptimal point: " << x_opt
        << "\nOptimal integer point: " << x_opt_int
        << "\nOptimal point for dual: " << u_opt
        << "\nOptimal integer point for dual: " << u_opt_int << '\n';
    cout << "\n*****************************************\n";
}


struct the_second_test_Shevchenko_ctrl_latexlog
{
    ostream& stream;

    the_second_test_Shevchenko_ctrl_latexlog(ostream& stream_a) : stream(stream_a) {}

    template <typename A, typename B, typename C>
    void preamble (const A& a, const B& b, const C& c) const
    {
        stream << "Дана задача линейного программирования:\n$$\\max cx$$\n"
            "$$\\left\\{ \\begin{tabular}{l}$Ax = b,$\\\\ $x \\ge 0,$ \\end{tabular}\\right.$$"
            "\nгде\n";
        stream << "$$A = ";
        output_latex(stream, a, true);
        stream << ",$$\n$$b = ";
        output_latex(stream, b, true);
        stream << "^T,$$\n$$c = ";
        output_latex(stream, c, true);
        stream << ".$$\nТребуется решить эту задачу при $x \\in R^{" << a.ncols()
            << "}$ и для $x \\in Z^{" << a.ncols() << "}$."
            " Построить двойственную и также "
            "решить для этих двух случаев.\\par";
        stream << endl;
    }

    basis_artificial_latexlog<ostream> ctrl_for_basis_artificial () const
    { return basis_artificial_latexlog<ostream>(stream, false, false); }

    primal_row_iters_latexlog<ostream> ctrl_for_primal_row_iters () const
    { return primal_row_iters_latexlog<ostream>(stream, false, false); }

    template <typename T, typename X>
    void primal_opt (const T& res, const X& x_opt) const
    {
        stream
            << "Итак, оптимальный вектор для прямой задачи $x = " << x_opt
            << "^T$, значение функции $" << res << "$.\n";
    }

    gomory1_iters_latexlog<ostream> ctrl_for_primal_task_gomory1_iters () const
    { return gomory1_iters_latexlog<ostream>(stream, false, false); }

    gomory1_iters_latexlog<ostream> ctrl_for_dual_task_gomory1_iters () const
    { return gomory1_iters_latexlog<ostream>(stream, false, false); }

    dual_col_iters_latexlog<ostream> ctrl_for_dual_col_iters () const
    { return dual_col_iters_latexlog<ostream>(stream, false, false); }

    template <typename T, typename X>
    void primal_opt_int (const T& res, const X& x_opt_int) const
    {
        stream
            << "Оптимальный вектор для прямой задачи в целых числах $x = " << x_opt_int
            << "^T$, значение функции $" << res << "$.\n";
    }

    template
    <
        typename Da, typename Db, typename Dc,
        typename T, typename Bsuba, typename Bsubc
    >
    void dual
    (
        const Da& da, const Db& db, const Dc& dc,
        const T& res_offset,
        const Bsuba& bsuba, const Bsubc& bsubc
    ) const
    {
        stream
            << "Двойственная задача:\n$$\\min ub$$\n"
            "$$uA \\ge c,$$\n";
        stream << "где $u = (u_1";
        for(size_t i = 2; i <= db.size(); ++i)
            stream << ", u_{" << i << "}";
        stream << ")$. Требуется решить эту задачу при $u \\in R^{" << da.nrows()
            << "}$ и для $x \\in Z^{" << da.nrows() << "}$."
            " Введём переменные незязки $y = uB - c_B \\ge 0$, где\n$$B = ";
        output_latex(stream, bsuba, true);
        stream << ",$$\n$$c_B = ";
        output_latex(stream, bsubc, true);
        stream << " $$ \nТогда исходную двойственную задачу можно переписать в новых"
            " переменных так:\n$$\\min c'x$$\n"
            "$$\\left\\{ \\begin{tabular}{l}$A'y^T \\le b',$\\\\ $y \\ge 0,$ \\end{tabular}\\right.$$"
            "\nгде\n";
        stream << "$$A' = ";
        output_latex(stream, da, true);
        stream << ",$$\n$$b' = ";
        output_latex(stream, db, true);
        stream << "^T,$$\n$$c' = ";
        output_latex(stream, dc, true);
        stream << "$$Довесок к целевой функции $" << res_offset << "$. Введением"
            " слабых переменных приведём задачу к каноническому виду.";
    }

    template <typename T, typename Y>
    void dual_y_opt (const T& res, const Y& y_opt) const
    {
        stream
            << "Таким образом, найден оптимальный вектор"
            " для двойственной задачи $y = " << y_opt
            << "^T$, значение функции $" << res << "$.\n";
    }

    template
    <
        typename Invbsuba, typename Bsubc,
        typename Y, typename U
    >
    void y_to_u_opt
    (
        const Invbsuba& invbsuba, const Bsubc& bsubc,
        const Y& y_opt, const U& u_opt
    ) const
    {
        stream <<
            "Восстановим значения оригинальных переменных $u$:\n"
            "$$u = (y + c_B)B^{-1} = \\left(y + ";
        output_latex(stream, bsubc, true);
        stream << "^T\\right)";
        output_latex(stream, invbsuba, true);
        stream << " = ";
        output_latex(stream, u_opt, true, eep_alone, false);
        stream << "^T.$$\n";
    }

    template <typename T, typename Y>
    void dual_y_opt_int (const T& res, const Y& y_opt) const
    {
        stream
            << "Таким образом, найден оптимальный целочисленный вектор"
            " для двойственной задачи $y = " << y_opt
            << "^T$, значение функции $" << res << "$.\n";
    }

    template
    <
        typename Invbsuba, typename Bsubc,
        typename Y, typename U
    >
    void y_to_u_opt_int
    (
        const Invbsuba& invbsuba, const Bsubc& bsubc,
        const Y& y_opt_int, const U& u_opt_int
    ) const
    {
        stream <<
            "Восстановим значения оригинальных переменных $u$:\n"
            "$$u = (y + c_B)B^{-1} = \\left(y + ";
        output_latex(stream, bsubc, true);
        stream << "^T\\right)";
        output_latex(stream, invbsuba, true);
        stream << " = ";
        output_latex(stream, u_opt_int, true, eep_alone, false);
        stream << "^T.$$\n";
    }

};


template <typename A, typename B, typename C>
void task_from_file (const char* filename, A& a, B& b, C& c)
{
    ifstream task(filename);

    if(!task)
    {
        cout << "\nERROR. Can't open file '" << filename << "'\n";
        exit(1);
    }

    task >> a >> b >> c;
}

template <typename A, typename B, typename C>
void task_from_file (const string& filename, A& a, B& b, C& c)
{ task_from_file(filename.c_str(), a, b, c); }


template <typename A, typename B, typename C, typename Basis>
void task_basis_from_file (const char* filename, A& a, B& b, C& c, Basis& basis)
{
    ifstream task(filename);

    if(!task)
    {
        cout << "\nERROR. Can't open file '" << filename << "'\n";
        exit(1);
    }

    task >> a >> b >> c >> basis;
}


template <typename A, typename B, typename C>
void print_task (const A& a, const B& b, const C& c)
{
    output_aligned(cout << "\nA =\n", a);
    cout << "b = " << b << "\nc = " << c << '\n';
    cout << "A*x = b, x >= 0.\n";
}


template <typename A, typename B>
void print_task (const A& a, const B& b)
{
    output_aligned(cout << "\nA =\n", a);
    cout << "b = " << b << '\n';
    cout << "A*x = b, x >= 0.\n";
}


void test4_1 (const std::string& file_with_problem)
{
    ofstream out("the_second_test_Shevchenko2.output.tex");

    typedef rational<> T;
    matrix<T, false> a;
    vector<T, false> b, c;

    task_from_file(file_with_problem.c_str(), a, b, c);

    the_second_test_Shevchenko
        (a, b, c, the_second_test_Shevchenko_ctrl_latexlog(out));
}


template <typename A, typename AA, typename B, typename Det, typename Basis>
void partial_rref_int_order
(
    const A& a,
    AA& aa, B& b, Det& det,
    const Basis& basis
)
{
    AA q; Basis basis_out;
    rref_order(a, aa, q, basis, basis_out, det);
    aa *= std::abs(det);
    // B = *b|N
    b = aa;
    output_aligned(cout << "\n|det(B)|*B^(-1)*(b|A) = \n", b);
    b.erase_cols(basis);
}


template <typename A1, typename A2>
void create_vector_order (const A1& a1, A2& a2, size_t n)
{
    a2 = a1;
    for(size_t i = 0; i < n; ++i)
        if(std::find(a1.begin(), a1.end(), i) == a1.end())
            a2.push_back(i);
}


void test4_2 ()
{
    // Применение алгоритма Моцкина-Бургера для области нашей задачи ЗЛП.

    typedef big_int T;
    typedef matrix<T, false> Matrix;
    typedef vector<T, false> Vector;

    Matrix a;
    Vector b, c;

    task_from_file("../../../../samples/old-samples/Shevchenko_course-3_v-896.task.txt", a, b, c);
    print_task(a, b, c);

    Vector basis = "(1, 3, 4)";

    //a.swap_cols(2, 3);
        rational<> det;

    {
        Matrix ab = a;
        ab.insert_col(0, b);

        matrix<rational<> > newa, bn;
        partial_rref_int_order(matrix<rational<> >(ab), newa, bn, det, basis);
        a = bn;
        output_aligned(cout << "PARTIAL_RREFINT(b|A) = \n", newa);
        output_aligned(cout << "\nbn = \n", bn);
    }

    //Vector nonbasis;
    //basis_to_nonbasis(basis, nonbasis, 0, a.ncols());
    //cout << "\nbasis = " << basis << '\n';
    //Matrix bmat = a.take_cols(basis);
    //Matrix Bm1N = a;

    a.opposite();
    //a.insert_col(0, b+1);
    a.mult_col(0, -1);
    Matrix ba_orig = a;
    for(int i = 0; i < a.ncols(); ++i)
    {
        Vector v(a.ncols());
        v[i] = 1;
        a.insert_row(i/*a.nrows()*/, v);
    }


    output_aligned(cout << "\nA =\n", a);
    //output_aligned(cout << "\nB =\n", bmat);
    //cout << "\nnonbasis = " << nonbasis;

    Matrix f, q, e;
    skeleton(a, f, q, e, ctrl::make_skeleton_slog(cout, false));

    matrix<rational<> > ff = f;
    for(int i = 0; i < ff.nrows(); ++i)
        if(!is_null(ff(i, 0)))ff.div_row(i, safe_reference(ff(i, 0)));

    output_aligned(cout << "\nff = \n", ff);

    //ff.erase_col(0);
    matrix<rational<> > xbdelta = transpose(matrix<rational<> >(ba_orig)*transpose(ff));
    output_aligned(cout << "\nxbdelta = \n", xbdelta);
    matrix<rational<> > xb = xbdelta/std::abs(det);
    output_aligned(cout << "\nxb = \n", xb);

    matrix<rational<> > pre_x = xb;

    for(size_t i = 0; i < ff.ncols(); ++i)
        pre_x.insert_col(pre_x.ncols(), ff.copy_col(i));

    output_aligned(cout << "\nx preorder = \n", pre_x);
    matrix<rational<> > x;

    Vector order;
    create_vector_order(basis, order, pre_x.ncols());
    Vector invorder;
    vec_inverse_permutation(order, invorder);

    mat_col_copy_order(pre_x, x, invorder);
    x.erase_col(0);

    output_aligned(cout << "\nx =\n", x);

    cout << "\n******* TRIANGULATION *********\n";
    output_aligned(cout << "\nQ from Skeleton = \n", q);

    matrix<int> tr;
    //triangulate_simple_1(q, tr);

    //output_aligned(cout << "\nsimplexes: \n", tr);
    cout << "This code cannot triangulate because a modification in base algorithm. Sorry.";

    //{
    //    typedef matrix<rational<> > MR;

    //    MR xN, v, q;

    //    xN = ff;
    //    xN.erase_col(0);

    //    MR xB = - MR(Bm1N)*transpose(xN);
    //    for(size_t i = 0; i < xB.nrows(); ++i)
    //        for(size_t j = 0; j < xB.ncols(); ++j)
    //            xB(i, j) += b[i];

    //    xB = transpose(xB);

    //    output_aligned(cout << "\nxN =\n", xN);
    //    output_aligned(cout << "\nxB =\n", xB);

    //}

    //output_aligned(cout << "\nff = \n", ff);
    //output_aligned(cout << "\nfff =\n", matrix<double>(ff));

    //output_aligned(cout << "\nB^(-1) = \n", inverse(matrix<rational<T> >(bmat)));
    //output_aligned(cout << "\n(b'|-A') =\n", ba_orig);
    //output_aligned(cout << "\nB^(-1)*(b'|-A') =\n", inverse(matrix<rational<T> >(bmat))*(ba_orig));
    //matrix<rational<T> > xout = transpose(inverse(matrix<rational<T> >(bmat))*(ba_orig)*transpose(f));
    //output_aligned(cout << "\n(B^(-1)*(b'|-A')*X')^T =\n", xout);

    //xout.insert_col(0, f.copy_col(0));
    //for(int i = 1; i < f.ncols(); ++i)
    //    xout.insert_col(xout.ncols(), f.copy_col(i));

    //output_aligned(cout << "\nX = \n", xout);
    //for(int i = 0; i < xout.nrows(); ++i)
    //    if(!is_null(xout(i, 0)))xout.div_row(i, safe_reference(xout(i, 0)));

    //output_aligned(cout << "\nX normalized =\n", xout);

}



void find_in_bounding_box
(
    const matrix<big_int> simplex, const big_int& det,
    Arageli::vector<big_int>& row_coefs, Arageli::vector<big_int>& res
)
{

    matrix<rational<> > rsimplex = simplex;
    vector<rational<> > up(rsimplex.ncols() - 1, -10000000), bot(rsimplex.ncols() - 1, 10000000);
    vector<rational<> > recup(rsimplex.ncols() - 1, 0), recbot(rsimplex.ncols() - 1, 0);
    for(size_t i = 0; i < rsimplex.nrows(); ++i)
    {
        if(rsimplex(i, 0) >= 1)
        {
            if(rsimplex(i, 0) > 1)rsimplex.div_row(i, safe_reference(rsimplex(i, 0)));

            for(size_t j = 0; j < up.size(); ++j)
            {
                rational<> curup = rsimplex(i, j+1), curbot = rsimplex(i, j+1);
                //std::cerr << "\nup = " << up << ", bot = " << bot;
                if(up[j] < curup)up[j] = curup;
                if(bot[j] > curbot)bot[j] = curbot;
            }
        }
        else if(rsimplex(i, 0) == 0)
        {
            // рецессивное направление:

            for(size_t j = 0; j < up.size(); ++j)
            {
                rational<> curup = rsimplex(i, j+1), curbot = rsimplex(i, j+1);
                if(recup[j] < curup)recup[j] = curup;
                if(recbot[j] > curbot)recbot[j] = curbot;
            }

        }
    }



    std::cerr << "\nup = " << up + recup << ", bot = " << bot + recbot << "\n";

    vector<big_int> rup = floor(up + recup), rbot = ceil(bot + recbot);

    std::cerr << "\nrup = " << rup << ", rbot = " << rbot << "\n";
    std::cerr << "Number of potential points: " << product(rup - bot + 1) << "\n";

    matrix<big_int> ainv = det*inverse(matrix<rational<> >(simplex));

    //output_aligned(cout << "\nSIMPLEX\n", ainv);

    vector<big_int> point = rbot;
    point.push_front(1);    // always

    for(;;)
    {
        vector<big_int> z = point*ainv;
        if(z == 0){ cerr << "ERROR1"; return; }
        if(Arageli::allcmp_vec_by_val(z, 0, std::greater_equal<big_int>()))
        {
            bool good = false;
            for(size_t i = 0; i < simplex.nrows(); ++i)
                if(simplex(i, 0) > 1 && z[i] != 0){ good = true; break; }

            if(good)
            {
                row_coefs = z;
                res = point;
                return;
            }
        }

        for(size_t i = 1; i < point.size(); ++i)
            if(++point[i] > rup[i-1])
                point[i] = rbot[i-1];
            else break;

        std::cerr << point;

        vector<big_int> tpoint = point;
        tpoint.erase(0);
        if(tpoint == rbot)break;
    }

    row_coefs.resize(0);
}


void picking_coefs
(
    const matrix<big_int> simplex, const big_int& det,
    Arageli::vector<big_int>& row_coefs, Arageli::vector<big_int>& res
)
{
    Arageli::vector<big_int> col_0 = simplex.copy_col(0), col_0_orig = col_0;
    std::replace(col_0.begin(), col_0.end(), big_int(), det+1);
    Arageli::vector<big_int> ups = det/col_0;
    if(is_null(col_0_orig))
    {
        cout << "ОШИБКА: одни рецессивные направления.#########################################";
        row_coefs.resize(0);
        return;
    }

    row_coefs.assign(col_0.size(), 0);
    row_coefs[0] = 1;

    Arageli::vector<int> nonintegers;
    for(size_t i = 0; i < simplex.nrows(); ++i)
        if(simplex(i, 0) > 1)nonintegers.push_back(i);

    cout << "\ndet = " << det << '\n';

    for(;;)
    {
        //cout << row_coefs << "\n";
        big_int sum = dotprod(row_coefs, col_0_orig);

        for(size_t i = 0; i < row_coefs.size(); ++i)
        {
            if(is_null(col_0_orig[i]))continue;
            ARAGELI_ASSERT_0(sum <= det);
            big_int delta = (det - sum)/col_0_orig[i];
            if(is_null(delta))
            {
                    sum -= row_coefs[i]*col_0_orig[i];
                    ARAGELI_ASSERT_0(!is_negative(sum));
                    row_coefs[i] = 0;
            }
            else
            {
                if(i == 0)
                    row_coefs[i] += delta;
                else
                    ++row_coefs[i];

                break;
            }
        }

        cerr << "\n" << row_coefs;

        if(is_null(row_coefs))break;

        if
        (
            dotprod(col_0_orig, row_coefs) == det &&
            !is_null(row_coefs.copy_subvector(nonintegers))
        )
        {
            // Проверяем систему сравнений.
            res.assign(simplex.ncols(), 0);
            for(size_t i = 0; i < simplex.nrows(); ++i)
                res += simplex.copy_row(i)*row_coefs[i];

            bool is_div = true;
            for(size_t i = 0; i < res.size(); ++i)
                if(!is_divisible(res[i], det)){ is_div = false; break; }

            if(is_div)
            {
                res /= det;
                return;
            }
        }

    }

    row_coefs.resize(0);
}




void solve_mod_2 (vector<big_int> cureq, big_int curb, big_int curmodule, matrix<big_int>& curres)
{
    ARAGELI_ASSERT_0(cureq.size() == 2);
    big_int first = cureq[0], alpha = cureq[1], gamma = curb;
    cout << "\nДвухмерная задача: ";
    for(;;)
    {
        cout << "\n(" << first << ", " << alpha << ") = " << gamma << " (mod " << curmodule << ")" << std::flush;

        if(gamma == 0)break;

        //alpha -= curmodule/2;
        if(first != 1)
        {
            // нормировка по первой переменной:
            big_int invfirst = inverse_mod(first, curmodule);
            cout
                << "\nНормируем по первой переменной: (" << first << "^(-1) = "
                << invfirst << "(mod " << curmodule << "))" << std::flush;
            first = first*invfirst%curmodule;
            alpha = alpha*invfirst%curmodule;
            gamma = gamma*invfirst%curmodule;

            cout << "\n(" << first << ", " << alpha << ") = " << gamma << " (mod " << curmodule << ")" << std::flush;
        }

        if(alpha == 1 || alpha == 0)break;
            // The Algorithm :)

            if(alpha > curmodule/2)
            {
                alpha -= curmodule;
                cout << "\n(" << first << ", " << alpha << ") = " << gamma << " (mod " << curmodule << ")" << std::flush;

                big_int newalpha = prrem(curmodule, -alpha);
                big_int newgamma = prrem(gamma - curmodule, -alpha);
                curmodule = -alpha;
                alpha = newalpha;
                gamma = newgamma;
            }
            else
            {
                big_int newalpha = prrem(-curmodule, alpha);
                gamma = prrem(gamma, alpha);
                curmodule = alpha;
                alpha = newalpha;
            }
    }
}



void test4_3 ()
{

    Arageli::vector<big_int> vect(2);
    vect[0] = 0;

}


void test4_4 ()
{
    typedef /*sparse_polynom<rational<> >*/rational<> P;
    typedef matrix<P> M;

    M a = "((-1, 1, 1), (2, -1, -1), (1, 0, 0), (0, 1, 0), (0, 0, 1))";
    M f, q, e;
    skeleton(a, f, q, e);

    output_aligned(cout << "\nF =\n", f);
    output_aligned(cout << "\nQ =\n", q);
    output_aligned(cout << "\nE =\n", e);

    a = "((1, 1, 0), (1, 2, 0), (1, 0, 1), (1, 0, 2))";
    skeleton(a, f, q, e);

    output_aligned(cout << "\n*********************\nF =\n", f);
    output_aligned(cout << "\nQ =\n", q);
    output_aligned(cout << "\nE =\n", e);
}


void test4_5 ()
{
    typedef sparse_polynom<rational<> > P;
    typedef matrix<P> M;

    M a = "((1, 1), (1, x), (1, 0))";
    M f, q, e;
    skeleton(a, f, q, e, ctrl::make_skeleton_slog(cout));

    output_aligned(cout << "\nF =\n", f);
    output_aligned(cout << "\nQ =\n", q);
    output_aligned(cout << "\nE =\n", e);

    a = "((-1, 1), (x, -1))";
    skeleton(a, f, q, e, ctrl::make_skeleton_slog(cout));

    output_aligned(cout << "\n*********************\nF =\n", f);
    output_aligned(cout << "\nQ =\n", q);
    output_aligned(cout << "\nE =\n", e);
}


void test4_6 ()
{
    typedef matrix<rational<int> > M;
    //M a = "((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, -1), (-1/2, 1, 1))";
    M a = "((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, -1), (-1/2, 1, 1))";
    M f, q, e;
    matrix<size_t> tr;
    skeleton(a, f, q, e);
    triangulate_simple_1(q, 3, tr);

    output_aligned(cout << "\nF =\n", f);
    output_aligned(cout << "\nQ =\n", q);
    output_aligned(cout << "\nE =\n", e);
    output_aligned(cout << "\nTr =\n", tr);

}


void test4_7 ()
{
    typedef residue<sparse_polynom<residue<int> > > T;
    T a = "( ((1(mod 2))*x) (mod ((1(mod 2))*x^3)) )";


    for(int i = 0; i < 30; ++i)
    {
        cout << a << '\n';
        a += T("1(mod 2) (mod ((1(mod 2))*x^3))");
    }

}


void test4_8 ()
{
    typedef std::map<std::string, Arageli::matrix<rational<> > > Vars;
    std::map<std::string, Arageli::matrix<rational<> > > vars;

    {
        std::ifstream varsfile("easymath.vars.txt");

        for(;;)
        {
            std::string line;
            std::getline(varsfile, line);
            if(line.empty())break;
            std::istringstream buffer(line);
            buffer >> line;
            buffer >> vars[line];
        }
    }

    try
    {
    for(;;)
    {
        {
            std::ofstream varsfile("easymath.vars.txt");

            for(Vars::iterator i = vars.begin(); i != vars.end(); ++i)
                varsfile << i->first << " " << i->second << "\n";
        }

        std::string line;
        cout << "\n>";
        std::getline(std::cin, line);
        std::istringstream buffer(line);
        std::string command;
        buffer >> command;

        if(command == "exit" || command == "quit")break;
        else if(command == "new")
        {
            buffer >> command;
            if(command == "matrix")
            {
                buffer >> command;
                buffer >> vars[command];
                if(!buffer && !buffer.eof())
                    cout << "Warning. Value is unset.";
            }
            else
                cout << "Error. Unknown type \"" << command << "\".";
        }
        else if(command == "delete")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
                vars.erase(i);
        }
        else if(command == "list")
        {
            for(Vars::iterator i = vars.begin(); i != vars.end(); ++i)
                Arageli::output_aligned(cout << i->first << " =\n", i->second);
        }
        else if(command == "print")
        {
            buffer >> command;

            if(command == "add" || command == "mul" || command == "sub")
            {
                std::string arg1, arg2;
                buffer >> arg1 >> arg2;
                if(arg2.empty())arg2 = arg1;
                Vars::mapped_type temp;

                Vars::iterator iarg1 = vars.find(arg1), iarg2 = vars.find(arg2);
                if(iarg1 == vars.end() || iarg2 == vars.end())
                {
                    cout
                        << "Error. Variable \"" << arg1 << "\" or \""
                        << arg2 << "\" are undefined.";
                }
                else
                    if(command == "add")
                        temp = iarg1->second + iarg2->second;
                    else if(command == "mul")
                        temp = iarg1->second * iarg2->second;
                    else if(command == "sub")
                        temp = iarg1->second - iarg2->second;
                    else if(command == "div")
                        temp = iarg1->second / iarg2->second;

                output_aligned(cout, temp);
            }
            else
            {
                Vars::iterator i = vars.find(command);
                if(i == vars.end())
                    cout << "Error. Variable \"" << command << "\" is undefined.";
                else
                    Arageli::output_aligned(cout << i->first << " =\n", i->second);
            }
        }
        else if(command == "set")
        {
            buffer >> command;

            if(command == "add" || command == "mul" || command == "sub" || command == "mov")
            {
                std::string arg1, arg2;
                buffer >> arg1 >> arg2;
                if(arg2.empty())arg2 = arg1;

                Vars::iterator iarg1 = vars.find(arg1), iarg2 = vars.find(arg2);
                if(iarg1 == vars.end() || iarg2 == vars.end())
                {
                    cout
                        << "Error. Variable \"" << arg1 << "\" or \""
                        << arg2 << "\" are undefined.";
                }
                else
                    if(command == "add")
                        iarg1->second += iarg2->second;
                    else if(command == "mul")
                        iarg1->second *= iarg2->second;
                    else if(command == "sub")
                        iarg1->second -= iarg2->second;
                    else if(command == "div")
                        iarg1->second /= iarg2->second;
                    else if(command == "mov")
                        iarg1->second = iarg2->second;

                cout << iarg1->first << " =\n";
                output_aligned(cout, iarg1->second);
            }
            else
            {
                buffer >> command;
                Vars::iterator i = vars.find(command);
                if(i == vars.end())
                    cout << "Error. Variable \"" << command << "\" is undefined.";
                else
                    buffer >> i->second;
            }
        }
        else if(command == "addmult")
        {
            size_t i1, i2;
            rational<> alpha, beta;
            std::string var;
            buffer >> command >> var >> i1 >> alpha >> i2 >> beta;

            Vars::iterator i = vars.find(var);
            if(i == vars.end())
                cout << "Error. Variable \"" << var << "\" is undefined.";
            else
            {
                if(command == "row")
                {
                    i->second.addmult_rows(i1, alpha, i2, beta);
                }
                else if(command == "col")
                {
                    i->second.addmult_cols(i1, alpha, i2, beta);
                }

                output_aligned(cout << i->first << " =\n", i->second);
            }

        }
        else if(command == "save")
        {
            std::ofstream varsfile("easymath.vars.txt");

            for(Vars::iterator i = vars.begin(); i != vars.end(); ++i)
                varsfile << i->first << " " << i->second << "\n";
        }
        else if(command == "det")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
                cout << det(i->second);
        }
        else if(command == "smith")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
            {
                matrix<big_int> P, Q, res;
                size_t rank;
                big_int det;
                smith(matrix<big_int>(i->second), res, P, Q, rank, det);
                output_aligned(cout << "Smith\n", res);
                output_aligned(cout << "P = \n", P);
                output_aligned(cout << "Q = \n", Q);
            }
        }
        else if(command == "solve_linsys_integer")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
            {
                buffer >> command;
                Vars::iterator j = vars.find(command);
                if(j == vars.end())
                    cout << "Error. Variable \"" << command << "\" is undefined.";
                else
                {
                    matrix<big_int> gx;
                    vector<big_int> offset;
                    solve_linsys_result status = solve_linsys_integer
                    (
                        matrix<big_int>(i->second),
                        vector<big_int>(j->second.size(), j->second.begin(), fromseq),
                        offset, gx
                    );
                    output_aligned(cout << "Solve integer linear system\ngx =\n", gx);
                    cout << "offset = " << offset << "\n";
                    cout << "status = " << slr_name(status) << '\n';
                }
            }
        }
        else if(command == "solve_linsys_integer_mod")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
            {
                matrix<big_int> gx;
                vector<big_int> order;
                solve_linsys_result status = solve_linsys_integer_mod
                (
                    matrix<big_int>(i->second),
                    abs(det(matrix<big_int>(i->second))),
                    gx, order
                );
                output_aligned(cout << "Solve integer linear system modulo\ngx =\n", gx);
                cout << "order = " << order << "\n";
                cout << "status = " << slr_name(status) << '\n';
            }
        }
        else if(command == "near_extreme_intpoint")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
            {
                vector<big_int> v;
                vector<std::size_t> rows;
                nei_result status = near_extreme_intpoint
                    (matrix<big_int>(i->second), v, rows);
                cout << "Find integer point in simplex\nv = " << v << "\nrows = " << rows;
                cout << "\nstatus = " << neir_name(status) << '\n';
            }
        }
        else if(command == "intconvex_triangulate_simple_1_par")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
            {
                matrix<big_int> gens = i->second;
                cone<> p(gens, fromgen);
                intconvex_triangulate_simple_1(p.generatrix(), p.relation(), p.dim(), gens);
                cone<> p1(gens, fromineq);
                output_aligned(cout << "\nintgens =\n", p1.generatrix()) << '\n';
            }
        }
        else if(command == "zgm2")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
            {
                matrix<big_int> agd = i->second;
                if(agd.nrows() != 1 || agd.ncols() != 3)
                    cout << "Error. Need 1x3 vector-row.";
                else
                {
                    matrix<big_int, false> v;
                    near_extreme_2d_mod(agd(0, 0), agd(0, 1), agd(0, 2), v, ctrl::intconvex_triangulation_idler());
                    output_aligned(cout << "\npoints =\n", v) << '\n';
                }
            }
        }
        else if(command == "ne23mod")
        {
            buffer >> command;
            Vars::iterator i = vars.find(command);
            if(i == vars.end())
                cout << "Error. Variable \"" << command << "\" is undefined.";
            else
            {
                buffer >> command;
                Vars::iterator j = vars.find(command);
                if(j == vars.end())
                    cout << "Error. Variable \"" << command << "\" is undefined.";
                else
                {
                    buffer >> command;
                    Vars::iterator k = vars.find(command);
                    if(k == vars.end())
                        cout << "Error. Variable \"" << command << "\" is undefined.";
                    else
                    {
                        matrix<big_int> v;
                        near_extreme_2_3_mod
                        (
                            matrix<big_int>(i->second),
                            vector<big_int>(j->second.size(), j->second.begin(), fromseq),
                            vector<big_int>(k->second.size(), k->second.begin(), fromseq),
                            product(vector<big_int>(k->second.size(), k->second.begin(), fromseq)),
                            v
                            , ctrl::intconvex_triangulation_idler()
                        );
                        output_aligned(cout << "\nv =\n", v);
                    }
                }
            }
        }
        else if(command == "clear")
        {
            vars.clear();
        }
        else if(command.length())
        {
            cout << "Error. Unknown command.";
        }
    }

    }
    catch(...)
    {
        std::ofstream varsfile("easymath.vars.txt");

        for(Vars::iterator i = vars.begin(); i != vars.end(); ++i)
            varsfile << i->first << " " << i->second << "\n";
    }

}


template <typename A>
void norm_first_col (A& a)
{
    for(size_t i = 0; i < a.nrows(); ++i)
        a.div_row(i, safe_reference(a(i, 0)));
}


void test4_9 ()
{
    typedef big_int T;
    matrix<T> a = "((18, 0, 0), (0, 18, 0), (0, 0, 18), (37, 7, 6), (40, 1, -21), (41, -1, 12))",
        f, q, e;
    skeleton(a, f, q, e, ctrl::make_skeleton_slog(cout, false));
    matrix<rational<> > qq;

    qq.insert_col(0, q.copy_col(3));
    qq.insert_col(1, q.copy_col(1));
    qq.insert_col(2, q.copy_col(4));
    qq.insert_col(3, q.copy_col(5));
    qq.insert_col(4, q.copy_col(2));

    qq.insert_col(0, q.copy_col(0));

    norm_first_col(qq);
    output_aligned(cout, qq);
}


void test4_10 ()
{
    cout << "\n\n ******* Integer Ring ******** \n\n";

    {
        typedef int T;

        matrix<T> a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;
        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det);

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

    cout << "\n\n ******* Rational Field ******** \n\n";

    {
        typedef rational<int> T;

        matrix<T> a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;
        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det);

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

    cout << "\n\n ******* Finite Field ******** \n\n";

    {
        typedef residue<int> T;

        matrix<T> a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;

        for(matrix<T>::iterator i = a.begin(); i < a.end(); ++i)
        {
            i->module() = (1 << next_mersen_prime_degree(10)) - 1;//11;
            i->normalize();
        }

        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det);

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

    cout << "\n\n *********** double *********** \n\n";

    {
        typedef double T;

        matrix<T> a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;

        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det);

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

    cout << "\n\n *********** big_float *********** \n\n";

    {
        //big_float::set_global_precision(100);
        //big_float::set_round_mode(TO_NEAREST);

        cout.setf(std::ios::scientific, std::ios::floatfield);

        typedef big_float T;

        matrix<T> a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;

        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det);

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nis it valid RREF?: " << (b == q * a);
    }
}


void test4_11 ()
{
    vector<rational<> > a = "(1, 1/2, 1/3, 1/4,    1/5)";
    vector<int> index = "(1, 2, 3, 3, 0, 1)";
    output_list(cout << '\n', a);
    output_list(cout << '\n', index);
    output_list(cout << '\n', a[index]);

    index = "(0, 2)";
    a[index] = a[index + 1];

    output_list(cout << "\n\n", a);

    //output_list(cout << '\n', a[index][index]);
    //output_list(cout << '\n', a[index][index][index]);

    //a[index][2] = 100;

    //output_list(cout << "\n\n", a);
    //output_list(cout << "\n\n", a[index]);
    //output_list(cout << '\n', a[index][index]);
    //output_list(cout << '\n', a[index][index][index]);

}


void test4_12 ()
{
    typedef residue<int> T;
    typedef matrix<T, false> M;

    M qi =
        "((1,  0,  0,  0,  0,  0,  0,  0),"
        " (2,  1,  7, 11, 10, 12,  5, 11),"
        " (3,  6,  4,  3,  0,  4,  7,  2),"
        " (4,  3,  6,  5,  1,  6,  2,  3),"
        " (2, 11,  8,  8,  3,  1,  3, 11),"
        " (6, 11,  8,  6,  2,  7, 10,  9),"
        " (5, 11,  7, 10,  0, 11,  7, 12),"
        " (3,  3, 12,  5,  0, 11,  9, 12))";

    for(M::iterator i = qi.begin(); i < qi.end(); ++i)
    {
        i->module() = 13;
        i->normalize();
    }

    qi -= M(qi.nrows(), eye);

    output_aligned(cout << "Q - I=\n", qi);

    M rreftqi, invqi;
    vector<int> basis;

    rref(transpose(qi), rreftqi, invqi, basis);
    rreftqi.erase_rows(vector<int>("(5, 6, 7)"));

    output_aligned(cout << "\nRREF((Q - I)^T) =\n", rreftqi);

    M bm = rreftqi.take_cols(basis);
    rreftqi.opposite();

    output_aligned(cout << "\nbasis matrix =\n", bm);
    output_aligned(cout << "\nnonbasis matrix =\n", rreftqi);


    for(size_t i = 0; i < rreftqi.ncols(); ++i)
    {
        vector<T> x(rreftqi.ncols());
        x[i] = 1;
        cout << "\nxbasis = " << x;
        cout << "\nxnonbasis = " << rreftqi*x;
    }
}


void test4_13 ()
{
    cout
        << "\n" << vector<int>("()")
        << "\n" << vector<int>("(1)")
        << "\n" << vector<int>("(1, 2, 3)")
        << "\n" << vector<int>("(1:5, 100)")
        << "\n" << vector<int>("(200, 5:-2)")
        << "\n" << vector<int>("(100, 0:5:25, 200)")
        << "\n" << vector<sparse_polynom<int> >("(0, x : x^2 : 3*x^2+x, 1)")
        << "\n" << vector<sparse_polynom<int> >("(x:x^2:3*x^2+x)")
        << "\n" << vector<sparse_polynom<int> >("(x-2:x+3)")
        << "\n" << vector<sparse_polynom<int> >("(x : x)");
    cout << "\n";

    vector<rational<> > v;
    cout << "\n" << v;
    v.insert(0, "(3, 2, 1)");
    cout << "\n" << v;
    v.insert(1, vector<rational<> >("(4/5:-1/5:1/5)") + 2);
    cout << "\n" << v;
    v.insert(v.size(), 3, rational<big_int>(-1));
    cout << "\n" << v;
}


void test4_14 ()
{
    //big_float::set_global_precision(100);
    //big_float::set_round_mode(TO_NEAREST);

    big_float
        a = "1.1234567890987654321",
        b = "82736418265418217823000000.0000000000000000000000000000000001111111";

    cout.setf ( std::ios::scientific, std::ios::floatfield );
    cout << "\na = " << a << "\nb = " << b << "\na + b = " << a + b;

    //cout << big_int(std::numeric_limits<__int64>::min());
}


void test4_15 ()
{
    matrix<int> m(32, 3200), b, p, q;
    std::size_t rank;
    int det;

    m(0, 0) = 345;
    m(0, 1) = -23;
    m(1, 0) = 77;
    m(1, 199) = 31;
    m(29, 310) = 11;

    smith(m, b, p, q, rank, det);

    output_aligned(cout, b.copy_rows(vector<int>("(1, 2, 3, 4, 5)")-1).copy_cols(vector<int>("(1, 2, 3, 4, 5)")-1));
}


void test4_16 ()
{
    typedef sparse_polynom<big_int, big_int> P1; typedef sparse_polynom<big_int> P2;
    output_aligned(cout, sylvester(P1("3*x^2-1"), P2("x^4+10*x^2-7")));
    cout << "\nres = " << resultant(P1("3*x^2-1"), P2("x^4+10*x^2-7"));

    typedef sparse_polynom<rational<> > P3;
    vector<P3> f;
    cout << "\n" << sqrfree_factorize_poly_rational(pow(P3("x^10-12*x^5+4*x^4-3*x^3-190"), 1), f);
    cout << "\n" << sqrfree_factorize_poly_rational(pow(P3("x^10-12*x^5+4*x^4-3*x^3-190"), 2), f);
    cout << "\n" << sqrfree_factorize_poly_rational(pow(P3("x^10-12*x^5+4*x^4-3*x^3-190"), 3), f);

    typedef sparse_polynom<P3> P4;

    P3 a = "2*x^3-x+17", b = "x^2+5";

    cout << "\n" << a.subs(P4("((x)-x)"));
    output_aligned(cout << "\n", sylvester(a.subs(P4("((x)-x)")), b));
    cout << resultant(a.subs(P4("((x)-x)")), b);


}


//template <typename P, typename Seg>
//double algebraic_to_double (const P& p, const Seg& seg)
//{
//    interval<double> fseg = seg;
//    interval_root_precise(p, fseg, 0.00000000001);
//    return fseg.first;
//}


void test4_17 ()
{
    typedef rational<big_int> T;
    typedef sparse_polynom<T> P;
    //typedef sparse_polynom<P> Pxy;
    typedef interval<T> Seg;
    //typedef vector<P, false> Polyvec;
    //typedef vector<Seg, false> Segvec;

    //Polyvec polyvec = "(x^2-2, x^2-2)";
    //Segvec segvec = "((-1, 2), (-2, 0))";
    //P p; Seg seg;
    //
    //P rslt = resultant(polyvec[0].subs(Pxy("((x)-x)")), polyvec[1]);

    //algebraic_func_rslt(polyvec, segvec, rslt, p, seg, my_func<Segvec>());

    //P p; Seg seg;

    //P ap = "x^2-2"; Seg as = "(1, 2)";
    //P bp = "x^3-3"; Seg bs = "(-2, -1)";

    //cout << std::setprecision(15);

    //algebraic_plus(ap, as, bp, bs, p, seg);
    //cout << "\n" << algebraic_to_double(p, seg);
    //algebraic_multiplies(p, seg, P("11*x-10"), Seg("10/11", "10/11"), p, seg);
    //cout << "\n" << algebraic_to_double(p, seg);
    //algebraic_divides(p, seg, P("x^5-12"), Seg(0, 12), p, seg);
    //cout << "\n" << algebraic_to_double(p, seg);

    std::time_t tm = time(0);

    algebraic<T, T>
        a("x^2-2", "(1, 2)"),
        b("x^3+3", "(-2, -1)"),
        c("11*x-10", "(10/11, 10/11)"),
        d("x^5-12", "(1, 12)");

    cout << "\na + b = " << a << " + " << b;
    a += b;
    cout << "\na + b = " << a << " = " << double(a);

    cout << "\na * c = " << a << " * " << c;
    a *= c;
    cout << "\na * c = " << a << " = " << double(a);

    cout << "\na / d = " << a << " / " << d;
    a /= d;
    cout << "\na / d = " << a << " = " << double(a);

    cout << "\na / a = " << a << " / " << a;
    a /= a;
    cout << "\na / a = " << a << " = " << double(a);
    cout << "\nis_null(a) = " << a.is_null();

    cout << "\n\nTime: " << time(0) - tm;
}


void test4_18 ()
{
    cout << "\n7->" << next_mersen_prime_degree(7);
    cout << "\n6->" << next_mersen_prime_degree(6) << "\n";

    typedef residue<int, Arageli::_Internal::module_2pm1<unsigned int, int> > T;

    {

        matrix<T> a = "((1, 2, 3), (4, 5, 6), (7, 8, 9))", b, q;

        for(matrix<T>::iterator i = a.begin(); i < a.end(); ++i)
        {
            i->module() = next_mersen_prime_degree(10);
            i->normalize();
        }

        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det);

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

}

void test4_19 ()
{
    cout << "\n\n ******* Intervals on double ******** \n\n";

    {
        typedef interval<double> T;

        matrix<T> a = "(((0.999, 1.001), (1.999, 2.001), (2.999, 3.001)), ((3.999, 4.001), (4.999, 5.001), (5.999, 6.001)), ((6.999, 7.001), (7.999, 8.001), (8.999, 9.001)))", b, q;
        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det, ctrl::make_rref_slog(cout));

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        output_aligned(cout << "\nq * a =\n", q * a);
        output_aligned(cout << "\nb - q * a =\n", b - q * a);
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

    cout << "\n\n ******* Intervals on rational ******** \n\n";

    {
        typedef interval<rational<> > T;

        matrix<rational<> > aa = "((1, 2, 3), (4, 5, 6), (7, 8, 9))";

        matrix<T> a = aa, b, q;

        for(matrix<T>::iterator i = a.begin(); i != a.end(); ++i)
        {
            i->first() -= rational<>(1, 1000);
            i->second() += rational<>(1, 1000);
        }

        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det, ctrl::make_rref_slog(cout));

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        output_aligned(cout << "\nq * a =\n", q * a);
        output_aligned(cout << "\nb - q * a =\n", b - q * a);
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

    cout << "\n\n ******* Intervals on rational functions ******** \n\n";

    {
        typedef interval<rational<sparse_polynom<rational<> > > > T;

        matrix<double> aa = "((1, 2, 3), (4, 5, 6), (7, 8, 9))";

        matrix<T> a = aa, b, q;

        for(matrix<T>::iterator i = a.begin(); i != a.end(); ++i)
        {
            i->first() -= rational<sparse_polynom<rational<> > >("x");
            i->second() += rational<sparse_polynom<rational<> > >("x");
        }

        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det, ctrl::make_rref_slog(cout));

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        output_aligned(cout << "\nq * a =\n", q * a);
        output_aligned(cout << "\nb - q * a =\n", b - q * a);
        cout << "\nis it valid RREF?: " << (b == q * a);
    }

    cout << "\n\n *********** big_float *********** \n\n";

    {
        //big_float::set_global_precision(100);
        //big_float::set_round_mode(TO_NEAREST);

        cout.setf(std::ios::scientific, std::ios::floatfield);

        typedef big_float T;

        //matrix<T> a = "(((0.999, 1.001), (1.999, 2.001), (2.999, 3.001)), ((3.999, 4.001), (4.999, 5.001), (5.999, 6.001)), ((6.999, 7.001), (7.999, 8.001), (8.999, 9.001)))", b, q;
        matrix<T> a = "(((1, 1), (2, 2), (3, 3)), ((4, 4), (5, 5), (6, 6)), ((7, 7), (8, 8), (9, 9)))", b, q;
        vector<size_t> basis;
        T det;

        rref(a, b, q, basis, det, ctrl::make_rref_slog(cout));

        output_aligned(cout << "a =\n", a);
        output_aligned(cout << "\nb =\n", b);
        output_aligned(cout << "\nq =\n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        output_aligned(cout << "\nq * a =\n", q * a);
        output_aligned(cout << "\nb - q * a =\n", b - q * a);
        cout << "\nis it valid RREF?: " << (b == q * a);
    }
}


void test4_20 ()
{
    typedef algebraic<> T;

    vector<interval<rational<> > > segs;
    sparse_polynom<rational<> > p = "x^3+6*x^2+5*x-5";
    //sparse_polynom<rational<> > p = "x^2-2";
    sturm<rational<> >(p, segs);

    cout << "\np = " << p;
    cout << "\nthe number of roots = " << segs.size();

    vector<T> roots(segs.size());

    for(size_t i = 0; i < segs.size(); ++i)
        roots[i] = T(p, segs[i]);

    output_aligned(cout << "\nroots =\n", roots);

    vector<double> froots = roots;

    output_aligned(cout << "\ndouble(roots) =\n", froots);

    vector<double> frootvals = p.subs(froots);

    output_aligned(cout << "\np(double(roots)) =\n", frootvals);

    vector<T> rootvals = p.subs(roots);

    output_aligned(cout << "\np(roots) =\n", rootvals);

    vector<double> rootfvals = rootvals;

    output_aligned(cout << "\ndouble(p(roots)) =\n", rootfvals);

    cout << "\nis_null(p(roots)) = " << is_null(rootvals);


}


void test4_21 ()
{
    matrix<big_int> a = "((-2, 1, 1), (2, -1, -1))", f, q, e;
    skeleton(a, f, q, e);

    output_aligned(cout << "\na = \n", a);
    output_aligned(cout << "\nf = \n", f);
    output_aligned(cout << "\nq = \n", q);
    output_aligned(cout << "\ne = \n", e);
}


void test4_22 ()
{
    typedef big_int T;
    typedef matrix<T> MT;
    MT a = "((-1, 1, 1), (2, -1, -1), (0, 1, 0), (0, 0, 1))"; T res;
    intcount_barvinok(a, res);

    cout << "barvinok: " << res;
}


void test4_23 ()
{
    typedef sparse_polynom<big_int> P;
    P p = 1; P x("x");
    for(int a = 1; a <= 50; ++a)
        p *= x - a;

    cout << p;
    vector<interval<rational<> > > lims;
    sturm<rational<> >(p, lims);

    output_aligned(cout << "\nlims =\n", lims);
}


void test4_24 ()
{
    typedef sparse_polynom<big_int> P;
    typedef interval<rational<> > SR;
    vector<P> ss;
    P p = "2*x-1";
    sturm_diff_sys(p, ss);
    SR seg;

    seg = "(1, -1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, -1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(0, 1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, 0)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1, 2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1/2, 2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, 1/2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1/2, 1/2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    p *= P("x-10");

    seg = "(1, -1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, -1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(0, 1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, 0)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1, 2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1/2, 2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, 1/2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1/2, 1/2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    p *= P("x+10");

    seg = "(1, -1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, -1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(0, 1)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, 0)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1, 2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1/2, 2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(-1, 1/2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);

    seg = "(1/2, 1/2)";
    cout << "\nseg = " << seg << ",  nr = " << sturm_number_roots(ss, seg);
}


void test4_25 ()
{
    typedef big_int T;
    typedef matrix<T, false> M;
    typedef vector<M::size_type, false> Basis;

    M a;

    a.resize(20, 20);
    for(M::iterator i = a.begin(); i != a.end(); ++i)
        *i = std::rand()/5000;

    output_aligned(cout << "A =\n", a);

    for(int i = 0; i < 2; ++i)
    {
        M b, q;
        Basis basis;
        M::size_type rank;
        T det;
        cout << "\n+++++++++ rref_gauss_bareiss +++++++++++";

        std::time_t tm = std::clock();

        rref_gauss_bareiss(a, b, q, basis, det);
        rank = basis.size();

        cout << "\ntime = " << std::clock() - tm;

        output_aligned(cout << "\nB = \n", b);
        output_aligned(cout << "\nQ = \n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nrank = " << rank;
    }

    for(int i = 0; i < 2; ++i)
    {
        M b, q;
        Basis basis;
        M::size_type rank;
        T det;
        cout << "\n+++++++++ bareiss +++++++++++";

        std::time_t tm = std::clock();

        bareiss(a, b, q, basis, det);
        rank = basis.size();

        cout << "\ntime = " << std::clock() - tm;

        output_aligned(cout << "\nB = \n", b);
        output_aligned(cout << "\nQ = \n", q);
        cout << "\nbasis = " << basis;
        cout << "\ndet = " << det;
        cout << "\nrank = " << rank;
    }

}


void test4_26 ()
{
    matrix<int> a = "((1, 1), (-1, -1))", f, q, e;
    //matrix<int> a, f, q, e;
    //a.resize(0, 2);

    skeleton(a, f, q, e);

    output_aligned(cout << "\na = \n", a);
    output_aligned(cout << "\nf = \n", f);
    output_aligned(cout << "\nq = \n", q);
    output_aligned(cout << "\ne = \n", e);
}


double rand_in_range (double x, double y)
{ return (double(rand())/RAND_MAX)*(y-x) + x; }

template <typename V1, typename V2>
double angle (const V1& v1, const V2& v2)
{
    return dotprod(v1, v2)/(std::sqrt(dotprod(v1, v1))*sqrt(dotprod(v2, v2)));
}


std::string num_to_mystr (int i)
{
    std::ostringstream buf; buf << i;
    std::string s = buf.str();
    for(int j = 0; j < s.length(); ++j)
        s[j] = s[j] - '0' + 'a';
    return s;
}


void test4_27 ()
{
    typedef big_int T;
    typedef polyhedron<T> Polyhedra;
    typedef matrix<rational<T> > Vertices;

    Vertices vertices(0, 3);

    const double prec = 1000000;

    #if 1

    for(Vertices::size_type i = 0; i < 200; ++i)
    {
        double phi = rand_in_range(0, 3.14159*2);
        double theta = rand_in_range(-3.14159/2, 3.14159/2);

        vector<rational<T> > vertex(3);

        double x = std::cos(phi)*std::cos(theta);
        double y = std::sin(phi)*std::cos(theta);
        double z = std::sin(theta);

        vertex[0] = rational<T>(prec*x, prec);
        vertex[1] = rational<T>(prec*y, prec);
        vertex[2] = rational<T>(prec*z, prec);

        vertices.insert_row(vertices.nrows(), vertex);
    }

    {
        Polyhedra p1(vertices, fromvert);
        p1.normalize_all();

        ofstream vrmlfile("sphere.wrl");

        output_vrml(vrmlfile, p1);
    }

    vertices.resize(0, 3);


    #endif

    #if 0

    matrix<rational<T> > rot = "((1/2, 1/2, 0), (-1/2, 1/2, 0), (0, 0, 1))";
    rot *= matrix<rational<T> >("((1/3, 0, 2/3), (0, 1, 0), (-2/3, 1, 1/3))");

    for(int i = 0; i < 3; ++i)
        rot.div_row(i, rational<T>(prec*std::sqrt(double(dotprod(rot.copy_row(i), rot.copy_row(i)))), prec));

    for(Vertices::size_type i = 0; i < 200; ++i)
    {
        vector<rational<T> > vertex(3);

        vertex[0] = rational<T>(prec*rand_in_range(-1, 1), prec);
        vertex[1] = rational<T>(prec*rand_in_range(-0.5, 0.5), prec);
        vertex[2] = rational<T>(prec*rand_in_range(-1, 1), prec);

        vertices.insert_row(vertices.nrows(), rot*vertex);
    }

    #endif


    Polyhedra p(vertices, fromvert);
    p.normalize_all();

    {
        ofstream vrmlfile("polytope1.wrl");

        output_vrml(vrmlfile, p);
    }

#if 0
    {
        matrix<double> vertices = p.vertices();

        output_aligned(cout << "vertices:\n", vertices, "ISMVec3f(", "),", ", ");

        vertices *= 5;
        vertices += 5;

        //output_aligned(cout << "vertices:\n", vertices, "ISMVec3f(", "),", ", ");

        Polyhedra::sideset sides = p.facets();    // all sides with dimention 2

        //matrix<double> ps(0, 2);
        std::ofstream file("polytope.tex");
        std::ofstream ssm("sphere-200.txt");

        output_aligned(ssm, vertices, "ISMVec3f(", "),", ", ");
        ssm << "\n\n";

        file << "\\documentclass[a4paper, 12pt]{article}\n\\usepackage[cp1251]{inputenc}"
            "\n\\usepackage[russian]{babel}\n\\usepackage{pstricks}\n\\begin{document}"
            "\n\\title{Arageli:\\\\Случайная триангуляция сферы\\\\300 вершин}"
            "\n\\maketitle\n";

        // making colors

        vector<double> lightdir = "(2, 1, -1)";

        double begin_col = 0.4, end_col = 0.95;
        int ncols = 200;

        for(int i = 0; i <= ncols; ++i)
        {
            double c = begin_col + i*(end_col-begin_col)/ncols;
            file << "\\definecolor{mygray" << num_to_mystr(i) << "}{rgb}{" << c << "," << c << "," << c << "}";
        }

        file << "\n\\begin{figure}[h]\n\\centering\n\\begin{pspicture}(10,10)\n";

        for(Polyhedra::sideset::const_iterator i = sides.begin(); i != sides.end(); ++i)
        {
            Polyhedra::sideset::vertex_indices_set vi = *i;
            for
            (
                Polyhedra::sideset::vertex_indices_set::const_iterator j = vi.begin();
                j != vi.end(); ++j
            )
            {
                ssm << *j << ", ";
            }
            ssm << "\n";

            std::size_t idx = i - sides.begin();
            if(is_positive(p.inequation(idx, 1)))
            {
                vector<double> n = p.facet_normal(idx);
                int ci = ncols*(angle(n, lightdir) + 1)/2;
                cout << ci << " ";
                if(ci < 0)ci = 0;
                if(ci > ncols)ci = ncols;

                file << "\\pspolygon[fillstyle=solid,linewidth=0,fillcolor=mygray" << num_to_mystr(ci) << "]";
                for
                (
                    Polyhedra::sideset::vertex_indices_set::const_iterator j = vi.begin();
                    j != vi.end(); ++j
                )
                {
                    file << "(" << vertices(*j, 1) << "," << vertices(*j, 2) << ")";
                }
                file << "\n";


                file << "\\psline[arrows=c-c](" << vertices(vi.back(), 1) << "," << vertices(vi.back(), 2) << ")";
                for
                (
                    Polyhedra::sideset::vertex_indices_set::const_iterator j = vi.begin();
                    j != vi.end()-1; ++j
                )
                {
                    file << "(" << vertices(*j, 1) << "," << vertices(*j, 2) << ")";
                    file << "\n\\psline[arrows=c-c](" << vertices(*j, 1) << "," << vertices(*j, 2) << ")";
                }
                file << "(" << vertices(vi.back(), 1) << "," << vertices(vi.back(), 2) << ")";
                file << "\n";
            }
        }

        file << "\n\\end{pspicture}\n\\end{figure}\n\\end{document}\n";

    }

#endif

}


void test4_28 ()
{
    typedef matrix<rational<> > Exactmat;
    typedef matrix<double> Roundmat;
    typedef polyhedron<> P;

    const size_t n = 20, m = 10;

    Roundmat rm(0, 3);
    //rm.reserve(3*n);

    for(size_t i = 0; i < n; ++i)
    {
        vector<double> x(3);
        x[0] = rand_in_range(-1, 1);
        x[1] = rand_in_range(-1, 1);
        x[2] = rand_in_range(-1, 1);

        rm.insert_row(rm.nrows(), x);
    }


    for(int i = 0;;)
    {
        // norm

        for(size_t j = 0; j < n; ++j)
            rm.div_row(j, std::sqrt(dotprod(rm.copy_row(j), rm.copy_row(j))));

        // view

        std::ostringstream buf;
        buf << std::setw(2) << std::setfill('0') << i;
        std::ofstream file(("mut." + buf.str() + ".wrl").c_str());
        P p(Exactmat(rm), fromvert);
        p.normalize_all();
        output_vrml(file, p);

        if(++i >= m)break;

        // shufle

        for(int k = 0; k < 10; ++k)
        {
            vector<double> sum(3);
            for(size_t j = 0; j < n; ++j)
                sum += rm.copy_row(j);

            vector<double> x(3);

            for(size_t j = 0; j < n; ++j)
                if(rand_in_range(0, 10) < 2)
                {
                    x[0] = rand_in_range(-0.1, 0.1);
                    x[1] = rand_in_range(-0.1, 0.1);
                    x[2] = rand_in_range(-0.1, 0.1);
                    rm.assign_row(j, (2*rm.copy_row(j) - sum + x));
                }
        }
    }
}


void test4_29 ()
{
    typedef big_int T;
    typedef cone<T> Cone;

    for(;;)
    {
        Cone::inequation_type ineq;

        cout << "\ninequations>";
        cin >> ineq;
        if(ineq.size() == 1 && ineq(0, 0) == -17)break;

        Cone c1(ineq, fromineq);

        output_list(cout << "\ncone =\n", c1);

        c1.normalize_all();

        output_list(cout << "\ncone after normalization =\n", c1);
        output_list(cout << "\ninequation =\n", c1.inequation());
        output_list(cout << "\ngeneratrix =\n", c1.generatrix());
    }
}


void test4_30 ()
{
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "((1))";
        matrix<int> tr;
        triangulate_simple_1(q, 1, tr);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "()";
        matrix<int> tr;
        triangulate_simple_1(q, 0, tr);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "((0 , 1), (1 , 0))";
        matrix<int> tr;
        triangulate_simple_1(q, 2, tr);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "((1 , 0), (0 , 1))";
        matrix<int> tr;
        triangulate_simple_1(q, 2, tr);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "()";
        matrix<int> tr;
        triangulate_simple_1(q, 1, tr, 1);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "()";
        matrix<int> tr;
        triangulate_simple_1(q, 2, tr, 2);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "()";
        matrix<int> tr;
        triangulate_simple_1(q, 3, tr, 3);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "((1))";
        matrix<int> tr;
        triangulate_simple_1(q, 1, tr, 1);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q = "((1))";
        matrix<int> tr;
        triangulate_simple_1(q, 2, tr, 2);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr);
    }
    cout << "\n***********\n";
    {
        matrix<bool, false> q =
            "((0 , 1 , 1 , 0),"
            " (0 , 0 , 1 , 1),"
            " (1 , 0 , 0 , 1),"
            " (1 , 1 , 0 , 0))";
        matrix<int> tr;
        triangulate_simple_1(q, 4, tr);
        output_aligned(cout << "q = (" << q.nrows() << ", " << q.ncols() << ") =\n", q);
        output_aligned(cout << "tr = (" << tr.nrows() << ", " << tr.ncols() << ") =\n", tr + 1);
    }
}


void test4_31 ()
{
    //int a, b;
    //char c;

    //std::cin >> a >> c >> b;
    //cout << a << c << b;

    cout << std::cin.getloc().name();
}


void test4_32 ()
{
    vector<big_int> a = "(1, 2, 3)";

    output_aligned(cout, a);
    output_aligned(cout, a, "[");
    output_aligned(cout, a, "[", "]");
}


void test4_33 ()
{
    {
        matrix<big_int> a = "((2, -3, -3), (2, 1, 0), (2, 0, 1), (1, 0, 0))", f, q, e;

        skeleton(a, f, q, e);

        a = f;
        if(!e.is_empty())
        {
            a.insert_matrix_bottom(e);
            vector<big_int> row(0, e.ncols());
            for(std::size_t i = 0; i < e.nrows(); ++i)
                row += e.copy_row(i);
            a.insert_row(a.nrows(), -row);
        }

        output_aligned(cout << "\nvertices =\n", a);
        intconvex_simple(a, f);
        output_aligned(cout << "\nipoints =\n", f);
        skeleton(f, a, q, e);
        output_aligned(cout << "\niineqs = \n", a);
        skeleton(a, f, q, e);
        output_aligned(cout << "\niverts = \n", f);
        output_aligned(cout << "\nibasis = \n", e);
    }
    {
        matrix<big_int> a = "((9, -2, -2), (1, 2, 0), (7, 0, -2), (1, 0, 0))", f, q, e;

        skeleton(a, f, q, e);

        a = f;
        if(!e.is_empty())
        {
            a.insert_matrix_bottom(e);
            vector<big_int> row(0, e.ncols());
            for(std::size_t i = 0; i < e.nrows(); ++i)
                row += e.copy_row(i);
            a.insert_row(a.nrows(), -row);
        }

        output_aligned(cout << "\nvertices =\n", a);
        intconvex_simple(a, f);
        output_aligned(cout << "\nipoints =\n", f);
        skeleton(f, a, q, e);
        output_aligned(cout << "\niineqs = \n", a);
        skeleton(a, f, q, e);
        output_aligned(cout << "\niverts = \n", f);
        output_aligned(cout << "\nibasis = \n", e);
    }
}


void test4_34 ()
{
    matrix<big_int, false> reces = "((1, 1, 1), (1, -1, 1), (1, 1, -1), (1, -1, -1), (1, 0, 2))";
    vector<bool, false> reces_mask(reces.nrows());
    reces_mask[0] = true;
    matrix<big_int, false> vert= "((0, 0, 0))";
    std::size_t n = power(2, reces_mask.size());

    for(std::size_t i = 1; i < n; ++i)
    {
        cout << "\nreces_mask = " << reces_mask;
        vector<big_int> v = "(0, 0, 0)";
        for(std::size_t j = 0; j < reces_mask.size(); ++j)
            if(reces_mask[j])
                v += reces.copy_row(j);
        vert.insert_row(vert.nrows(), v);

        // move to the next combination
        for(std::size_t j = 0; j < reces_mask.size(); ++j)
            if(reces_mask[j])reces_mask[j] = false;
            else
            {
                reces_mask[j] = true;
                break;
            }
    }

    polyhedron<> p(vert, fromivert);
    std::ofstream file("parallelepiped.wrl");
    p.normalize_all();
    output_vrml(file, p);
}


void test4_35 ()
{
    matrix<big_int> a = "((1, 0, 0), (0, 0, 1), (-3, 1, 4))", b, q;
    output_aligned(cout << "a =\n", a);
    cout << "\ndet(a) = " << det(a);
    cout << "\nRREF";
    vector<size_t, false> basis;
    big_int det;
    rref_gauss_bareiss(a, b, q, basis, det, ctrl::make_rref_gauss_bareiss_slog(cout));
}


void test4_36 ()
{
    matrix<sparse_polynom<rational<big_int> > > a = "((1, 0), (0, 13*x-7))", b, p, q;
    size_t rank;
    sparse_polynom<rational<big_int> > det;
    smith(a, b, p, q, rank, det);
    cout << "b = " << b;
    cout << "It is valid: " << (b == p*a*q);
}


void test4_37 ()
{
    typedef sparse_polynom<rational<> > P;

    {
        P a = "17*x-23", b = "6*x-23", u, v;
        cout << "\neuclid_bezout(a, b, u, v) = " << euclid_bezout(a, b, u, v);
        cout << "\na*u + b*v = " << a*u + b*v;
    }
    {
        P a = "17*x", b = "6*x", u, v;
        cout << "\neuclid_bezout(a, b, u, v) = " << euclid_bezout(a, b, u, v);
        cout << "\na*u + b*v = " << a*u + b*v;
    }
    {
        P a = "16*x^2-4", b = "4*x+2", u, v;
        cout << "\neuclid_bezout(a, b, u, v) = " << euclid_bezout(a, b, u, v);
        cout << "\na*u + b*v = " << a*u + b*v;
    }
}


void test4_38 ()
{
    // Working in ${\bf Q}/(x^2 + 1){\bf Q}$

    typedef residue<sparse_polynom<rational<> > > T;

    T a = "((1+x) (mod (x^2+1)))";
    T b = "((1-x) (mod (x^2+1)))";

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "a + b = " << a + b << endl;
    cout << "a - b = " << a - b << endl;
    cout << "a * b = " << a * b << endl;
    cout << "a / b = " << a / b << endl << endl;
}


void test4_39 ()
{
    typedef big_int T;
    matrix<T> a = "((1, 2), (3, 4))";
    vector<T> b = "(5, 6)";

    {
        typedef int TT;
        matrix<TT> aa = a, gx;
        vector<TT> bb = b, offset;
        solve_linsys_result status = solve_linsys_integer(aa, bb, offset, gx);
        output_aligned(cout << "aa =\n", aa);
        cout << "bb = " << b;
        cout << "\noffset = " << offset;
        output_aligned(cout << "\ngx =\n", gx);
        cout << "status = " << slr_name(status) << "\n";
    }

    {
        typedef big_int TT;
        matrix<TT> aa = a, gx;
        vector<TT> bb = b, offset;
        solve_linsys_result status = solve_linsys_integer(aa, bb, offset, gx);
        output_aligned(cout << "aa =\n", aa);
        cout << "bb = " << b;
        cout << "\noffset = " << offset;
        output_aligned(cout << "\ngx =\n", gx);
        cout << "status = " << slr_name(status) << "\n";
    }

    {
        typedef rational<int> TT;
        matrix<TT> aa = a, gx;
        vector<TT> bb = b, offset;
        solve_linsys_result status = solve_linsys_rational(aa, bb, offset, gx);
        output_aligned(cout << "aa =\n", aa);
        cout << "bb = " << b;
        cout << "\noffset = " << offset;
        output_aligned(cout << "\ngx =\n", gx);
        cout << "status = " << slr_name(status) << "\n";
    }

    {
        typedef rational<big_int> TT;
        matrix<TT> aa = a, gx;
        vector<TT> bb = b, offset;
        solve_linsys_result status = solve_linsys_rational(aa, bb, offset, gx);
        output_aligned(cout << "aa =\n", aa);
        cout << "bb = " << b;
        cout << "\noffset = " << offset;
        output_aligned(cout << "\ngx =\n", gx);
        cout << "status = " << slr_name(status) << "\n";
    }

}


void test4_40 ()
{
    cout
        << "double(0) = " << double(0)
        << "\nbig_float(0) = " << big_float(0) << std::endl;
}


namespace Arageli
{

namespace ctrl
{
    // Класс контролёра, который ничего не делает.
    // Нужен для того, что бы пользователь мог вызывать функцию так,
    // как будтно она неконтролируемая.
    struct sqrt_approx_idler
    {
        class abort : public ctrl::abort {};

        // Вызывается в начале работы функции; принимает аргументы.
        template <typename T1, typename T2>
        void preamble (const T1& x, const T2& tolerance) const {}

        // Вызывается на каждой итерации извлечения корня;
        // принимает текущее и следующее приближения.
        template <typename T1>
        void iter
        (
            const T1& cur,    // текущее приближение a[n]
            const T1& next,   // следующее приближение a[n+1]
            const T1& curtol  // |a[n+1] - a[n]|
        ) const {}

        // Вызывается в конце работы функции; принимает результат.
        template <typename T1>
        void conclusion (const T1& res) const {}
    };
}

// Контролируемая функция извлечения квадратного корня.
template <typename T1, typename T2, typename Ctrler>
T1 sqrt_approx (const T1& x, const T2& tolerance, Ctrler ctrler)
{
    ctrler.preamble(x, tolerance);

    if(is_null(x))
    {
        ctrler.conclusion(x);
        return x;
    }

    T1 next = x, cur, curtol;

    do
    {
        cur = next;
        next = (cur + x/cur)/2;
        curtol = abs(next - cur);
        ctrler.iter(cur, next, curtol);
    }while(curtol > tolerance);

    ctrler.conclusion(next);
    return next;
}

// Вызов функции без контроля.
template <typename T1, typename T2>
inline T1 sqrt_approx (const T1& x, const T2& tolerance)
{ return sqrt_approx(x, tolerance, ctrl::sqrt_approx_idler()); }

namespace ctrl
{
    // Контролёр для sqrt_approx_slog выводящий промежуточные
    // результаты в поток в простом текстовом формате
    template <typename Out>
    class sqrt_approx_slog : public sqrt_approx_idler
    {
        Out& out;
        std::size_t curiter;

    public:

        sqrt_approx_slog (Out& out_a) : out(out_a), curiter(1) {}

        template <typename T1, typename T2>
        void preamble (const T1& x, const T2& tolerance) const
        {
            out << "Находим приближение к квадратному корню из "
                << x << '\n'
                << "константа, контролирующая точность: "
                << tolerance << '\n';
        }

        template <typename T1>
        void iter
        (
            const T1& cur,
            const T1& next,
            const T1& curtol
        )
        {
            out << "a[" << curiter++ << "] = " << cur
                << ", отклонение: " << curtol << '\n';
        }

        template <typename T1>
        void conclusion (const T1& res) const
        {
            out << "Найдено приближение a[" << curiter << "] = "
                << res << '\n'
                << "всего итераций: " << curiter-1 << '\n';
        }
    };

    // Для удобства: make-функция.
    template <typename Out>
    inline sqrt_approx_slog<Out> make_sqrt_approx_slog (Out& out)
    { return sqrt_approx_slog<Out>(out); }
}
} // namespace Arageli


void test4_41 ()
{
    std::cout
        << std::setprecision(8) << sqrt_approx(5.0, 1.0/100000) << '\n'
        << sqrt_approx(rational<>(5), rational<>(1, 100000));

    std::cout << "\n****************\n";

    sqrt_approx
    (
        5.0, 1.0/100000,
        ctrl::make_sqrt_approx_slog(std::cout)
    );

    std::cout << '\n';

    sqrt_approx
    (
        rational<>(5), rational<>(1, 100000),
        ctrl::make_sqrt_approx_slog(std::cout)
    );
}


void test4_42 ()
{
    {
        typedef big_int T;

        for(int i = 1; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::ipositive<T> S;
            S s(i);
            rnd::equiprob<S> rnd(s);
            std::map<int, int> hist;
            for(int j = 0; j < 100*i; ++j)
                ++hist[rnd()];

            for(std::map<int, int>::iterator j = hist.begin(); j != hist.end(); ++j)
                cout << "[" << j->first << "] = " << j->second << "\n";

            cout << "\n";
        }
    }
    {
        typedef big_int T;

        for(int i = 0; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::inonnegative<T> S;
            S s(i);
            rnd::equiprob<S> rnd(s);
            std::map<int, int> hist;
            for(int j = 0; j < 100*(i+1); ++j)
                ++hist[rnd()];

            for(std::map<int, int>::iterator j = hist.begin(); j != hist.end(); ++j)
                cout << "[" << j->first << "] = " << j->second << "\n";

            cout << "\n";
        }
    }
    {
        typedef unsigned T;

        for(int i = 1; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::ipositive<T> S;
            S s(i);
            rnd::equiprob<S> rnd(s);
            std::map<int, int> hist;
            for(int j = 0; j < 100*i; ++j)
                ++hist[rnd()];

            for(std::map<int, int>::iterator j = hist.begin(); j != hist.end(); ++j)
                cout << "[" << j->first << "] = " << j->second << "\n";

            cout << "\n";
        }
    }
    {
        typedef unsigned T;

        for(int i = 0; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::inonnegative<T> S;
            S s(i);
            rnd::equiprob<S> rnd(s);
            std::map<int, int> hist;
            for(int j = 0; j < 100*(i+1); ++j)
                ++hist[rnd()];

            for(std::map<int, int>::iterator j = hist.begin(); j != hist.end(); ++j)
                cout << "[" << j->first << "] = " << j->second << "\n";

            cout << "\n";
        }
    }
    {
        typedef rational<int> T;

        for(int i = 1; i <= 10; ++i)
        {
            cout << "i = " << i;
            typedef set::grid1<T> S;
            S s(-1, 1, T(3, i));
            cout << "\nstep = " << s.step() << std::endl;
            rnd::equiprob<S> rnd(s);
            std::map<T, int> hist;
            for(int j = 0; j < 100*i; ++j)
                ++hist[rnd()];

            for(std::map<T, int>::iterator j = hist.begin(); j != hist.end(); ++j)
                cout << "[" << j->first << "] = " << j->second << "\n";

            cout << "\n";
        }
    }
}


void test4_43 ()
{
    {
        typedef big_int T;

        for(int i = 0; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::ipositive<T> S;
            S s(i);
            typedef enumer::ordered<S> E;
            for(E j(enumer::begin, s); !j.is_end(); ++j)
                cout << *j << ", ";
            cout << "\n\n";
        }
    }
    {
        typedef big_int T;

        for(int i = 0; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::inonnegative<T> S;
            S s(i);
            typedef enumer::ordered<S> E;
            for(E j(enumer::begin, s); !j.is_end(); ++j)
                cout << *j << ", ";
            cout << "\n\n";
        }
    }
    {
        typedef unsigned T;

        for(int i = 0; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::ipositive<T> S;
            S s(i);
            typedef enumer::ordered<S> E;
            for(E j(enumer::begin, s); !j.is_end(); ++j)
                cout << *j << ", ";
            cout << "\n\n";
        }
    }
    {
        typedef unsigned T;

        for(int i = 0; i <= 10; ++i)
        {
            cout << "i = " << i << "\n" << std::endl;
            typedef set::inonnegative<T> S;
            S s(i);
            typedef enumer::ordered<S> E;
            for(E j(enumer::begin, s); !j.is_end(); ++j)
                cout << *j << ", ";
            cout << "\n\n";
        }
    }
    {
        typedef rational<int> T;

        for(int i = 1; i <= 10; ++i)
        {
            cout << "i = " << i;
            typedef set::grid1<T> S;
            S s(-1, 1, T(3, i));
            cout << "\nstep = " << s.step() << std::endl;
            typedef enumer::ordered<S> E;
            for(E j(enumer::begin, s); !j.is_end(); ++j)
                cout << *j << ", ";
            cout << "\n\n";
        }
    }
    {
        typedef rational<big_int> T;

        for(int i = 1; i <= 10; ++i)
        {
            cout << "i = " << i;
            typedef set::grid1<T> S;
            S s(-1, 1, T(3, i));
            cout << "\nstep = " << s.step() << std::endl;
            typedef enumer::ordered<S> E;
            for(E j(enumer::begin, s); !j.is_end(); ++j)
                cout << *j << ", ";
            cout << "\n\n";
        }
    }
    {
        typedef rational<big_int> T;

        for(int i = 1; i <= 10; ++i)
        {
            cout << "i = " << i;
            typedef set::grid1<T> S;
            S s(-1, 1, T(3, i));
            cout << "\nstep = " << s.step() << std::endl;
            typedef enumer::ordered<S> E;
            for(E j(enumer::begin, s); j != E(enumer::end, s); ++j)
                cout << *j << ", ";
            cout << "\n\n";
        }
    }
}


void test4_44 ()
{

    typedef rational<big_int> T;
    typedef set::grid1<T> BS;
    typedef vector<T> V;
    typedef set::vector_upto_size<V, BS> S;
    typedef enumer::fast<S> E;

    for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
    {
        S s(i, BS(-1, 1, T(1, (j+1))));
        cout << "i = " << i << "\n";

        for(E e(enumer::begin, s); !e.is_end(); ++e)
            cout << *e << "\n";

        cout << "\n";

        rnd::equiprob<S> r(s);
        std::map<V, int> hist;
        int n = 1000;
        for(int k = 0; k < n; ++k)
            ++hist[r()];

        for(std::map<V, int>::iterator e = hist.begin(); e != hist.end(); ++e)
            cout << "[" << e->first << "] = " << e->second << "\n";

        cout << "\n";
    }
}


void test4_45 ()
{
    typedef rational<big_int> T;
    typedef set::grid1<T> BS;
    typedef vector<T> V;
    typedef set::vector_fixed_size<V, BS> VS;
    typedef set::byvector<sparse_polynom<T>, VS> S;
    typedef enumer::fast<S> E;

    for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
    {
        S s(VS(i, BS(-1, 1, T(1, (j+1)))));
        cout << "i = " << i << "\n";

        for(E e(enumer::begin, s); !e.is_end(); ++e)
            cout << *e << "\n";

        cout << "\n";

        rnd::equiprob<S> r(s);
        std::map<sparse_polynom<T>, int> hist;
        int n = 1000;
        for(int k = 0; k < n; ++k)
            ++hist[r()];

        for(std::map<sparse_polynom<T>, int>::iterator e = hist.begin(); e != hist.end(); ++e)
            cout << "[" << e->first << "] = " << e->second << "\n";

        cout << "\n";
    }
}


void test4_46 ()
{
    sparse_polynom<residue <big_int> > P = "((1 (mod 3))*x)";
    std::cout << P;
}


void test4_47 ()
{
    typedef int T;
    rnd::equiprob_alter<set::inonnegative<T> > myrnd;

    for(int i = 0; i < 5; ++i)
    {
        std::map<int, int> hist;
        for(int j = 0; j < 100; ++j)
            ++hist[myrnd(i)];

        for(std::map<int, int>::iterator e = hist.begin(); e != hist.end(); ++e)
            cout << "[" << e->first << "] = " << e->second << "\n";

        cout << "\n";
    }
}


void test4_48 ()
{
    typedef int T;
    typedef set::inonnegative<T> Set;
    typedef set::subsets_fixnum<Set> Subsets;
    typedef rnd::equiprob<Subsets> Rnd_subsets;

    for(int i = 0; i <= 6; ++i)
    {
        Subsets ss(Set(6), i);
        Rnd_subsets myrnd(ss);

        typedef std::map<Subsets::value_type, int> Hist;
        Hist hist;
        for(int j = 0; j < 10000; ++j)
            ++hist[myrnd()];

        for(Hist::iterator e = hist.begin(); e != hist.end(); ++e)
            cout << "[" << e->first << "] = " << e->second << "\n";

        cout << "\n";
    }
}


void test4_49 ()
{
    {
        typedef int T;
        typedef set::inonnegative<T> Set;
        typedef set::subsets_fixnum<Set> Subsets;
        typedef rnd::equiprob<Subsets>::base_type Basernd;
        typedef rnd::equiprob<Subsets, Basernd> Rnd_subsets;
        typedef rnd::equiprob<Set, Basernd> Rnd_other;

        Set os(100);
        Subsets ss(Set(15), 4);

        Basernd basernd;
        Rnd_subsets myrnd1(basernd, ss);
        Rnd_other myrnd2(basernd, os);

        for(int i = 0; i < 5; ++i)
            cout << "myrnd1() = " << myrnd1() << "\nmyrnd2() = " << myrnd2() << "\n";

        cout << "\n";
    }

    {
        typedef int T;
        typedef set::inonnegative<T> Set;
        typedef set::subsets_fixnum<Set> Subsets;
        typedef rnd::equiprob<Subsets>::base_type Basernd;
        typedef rnd::equiprob<Subsets, Basernd> Rnd_subsets;
        typedef rnd::equiprob<Set, Basernd> Rnd_other;

        Set os(100);
        Subsets ss(Set(15), 4);

        Basernd basernd(2);
        Rnd_subsets myrnd1(basernd, ss);
        Rnd_other myrnd2(basernd, os);

        for(int i = 0; i < 5; ++i)
            cout << "myrnd1() = " << myrnd1() << "\nmyrnd2() = " << myrnd2() << "\n";

        cout << "\n";
    }

    {
        typedef int T;
        typedef set::inonnegative<T> Set;
        typedef set::subsets_fixnum<Set> Subsets;
        typedef rnd::equiprob<Subsets>::base_type Basernd;
        typedef rnd::equiprob<Subsets, rnd::ref<Basernd> > Rnd_subsets;
        typedef rnd::equiprob<Set, rnd::ref<Basernd> > Rnd_other;

        Set os(100);
        Subsets ss(Set(15), 4);

        Basernd basernd;
        Rnd_subsets myrnd1(rnd::ref<Basernd>(basernd), ss);
        Rnd_other myrnd2(rnd::ref<Basernd>(basernd), os);

        for(int i = 0; i < 5; ++i)
            cout << "myrnd1() = " << myrnd1() << "\nmyrnd2() = " << myrnd2() << "\n";

        cout << "\n";
    }

    {
        typedef int T;
        typedef set::inonnegative<T> Set;
        typedef set::subsets_fixnum<Set> Subsets;
        typedef rnd::equiprob<Subsets>::base_type Basernd;
        typedef rnd::equiprob<Subsets, rnd::ref<Basernd> > Rnd_subsets;
        typedef rnd::equiprob<Set, rnd::ref<Basernd> > Rnd_other;

        Set os(100);
        Subsets ss(Set(15), 4);

        Basernd basernd(2);
        Rnd_subsets myrnd1(rnd::ref<Basernd>(basernd), ss);
        Rnd_other myrnd2(rnd::ref<Basernd>(basernd), os);

        for(int i = 0; i < 5; ++i)
            cout << "myrnd1() = " << myrnd1() << "\nmyrnd2() = " << myrnd2() << "\n";

        cout << "\n";
    }
}


void test4_50 ()
{
    typedef algebraic<> T;
    T
        a("x^4+10*x^2-10", "(-1, 0)"),
        b("x^4+11*x^2-10", "(0, 1)");

    std::cout << "a = " << a << " = " << double(a) << "\n";
    std::cout << "b = " << b << " = " << double(b) << "\n";
    std::cout << "-a = " << -a << " = " << double(-a) << "\n";
    std::cout << "-b = " << -b << " = " << double(-b) << "\n";
    std::cout << "a + b = " << a + b << " = " << double(a + b) << "\n";
    std::cout << "a - b = " << a - b << " = " << double(a - b) << "\n";
    std::cout << "a * b = " << a * b << " = " << double(a*b) << "\n";
    std::cout << "a / b = " << a / b << " = " << double(a/b) << "\n";
}


void test4_51 ()
{
    generic_type x;
    std::cout << x << "\n";
    x = new generic_big_int(5);
    generic_type y = new generic_big_int(6);

    std::cout << -x << "\n";
    std::cout << -y << "\n";
    std::cout << x + y << "\n";
    std::cout << x - y << "\n";
    std::cout << x * y << "\n";
    std::cout << x / y << "\n";
    std::cout << x % y << "\n";

    Arageli::vector<generic_type> vx, vy(15, y);
    std::cout << vx << "\n" << vy << "\n";
    vx.resize(15, x);
    std::cout << vx + vy << "\n";

    big_int xx = 2, yy = 2;
    std::cout << Arageli::power(xx, yy);
}


void test4_52 ()
{
    typedef std::complex<rational<> > C;
    typedef sparse_polynom<rational<> > P;
    typedef sparse_polynom<sparse_polynom<C> > PPC;

    PPC ppc;
    ppc += PPC::monom(PPC::coef_type::monom(C(1)), 1);
    ppc += PPC::monom(PPC::coef_type::monom(C(0, 1), 1), 0);
    P p = "2*x^3-17*x+5";
    std::cout
        << "ppc = " << ppc << "\n"
        << "p = " << p << "\n";
    PPC subres = p.subs(ppc);
    std::cout << "p.subs(ppc) = " << subres;
}

void test4_53 ()
{
    timer tm;
    for(big_int i = 0; i < big_int("10000"); ++i);
    std::cout << tm.time() << '\n' << tm.precision();
    tm.stop();
    std::cout << "\n" << tm;
    std::cin >> tm;
    std::cout << tm.time() << '\n' << tm;

}

void test4_54 ()
{
    vector<int> b, c = "(1, 2)";
    matrix<int> matr = "((1, 2), (3, 4))";
    b = solve_linsys(matr, c);
    output_aligned(std::cout, b);
    std::cout << "(matr*b == c) = " << (matr*b == c);
}

void test4_55 ()
{
    matrix<big_int> m = "((1, 0), (0, 1))";
    cone<> c1(m, fromgen);
    std::cout << f_vector_cone(c1);
}

void do_some_time (int i)
{
    static volatile int x = 0;
    for(int j = 0; j < i; ++j)
        x += j;
}

void test4_56 ()
{
    double precision;
    std::cout << "precision: ";
    std::cin >> precision;

    const int SIZE = 10;
    typedef Arageli::vector<double> V;
    V times(SIZE), x(SIZE), prec(SIZE);

    for(int i = 0; i < SIZE; ++i)
    {
        int j = 0;
        timer tm;
        do { do_some_time(1000000*i); ++j; } while(j < precision/*tm.precision() > precision*/);
        tm.stop();
        times[i] = tm.time()/j;
        prec[i] = tm.precision();
        x[i] = i;
    }

    typedef cartesian2d_chart_line<V, V> Line;
    typedef Arageli::vector<Line> Lines;
    Lines lines;
    lines.push_back(Line("time", x, times, "[linecolor=black]"));
    lines.push_back(Line("precision", x, times*prec, "[linecolor=green]"));
    std::ofstream file("timer.times.chart.tex");

    cartesian2d_chart chart
    (
        16, 20, 1, 1.5,
        "промежуток", "Время, сек."
    );

    output_pstricks_cartesian2d(file, 0, 0, chart, lines);
}

void test4_57 ()
{
    multimonom<rational<int> >
        mm1(1, vector<rational<> >("(1, 2, 3)")),
        mm2(1, vector<rational<> >("(3, 2, 1)")),
        mm_gcd = gcd(mm1, mm2),
        mm_lcm = lcm(mm1, mm2);

    std::cout
        << "mm1.multidegree() = " << mm1.multidegree() << '\n'
        << "mm2.multidegree() = " << mm2.multidegree() << '\n'
        << "mm_gcd.multidegree() = " << mm_gcd.multidegree() << '\n'
        << "mm_lcm.multidegree() = " << mm_lcm.multidegree();
}


void test4_58 ()
{
    typedef vector<int> MD;
    typedef std::complex<rational<> > F;
    typedef sparse_multipolynom<F> P;

    Arageli::vector<P> vp;
    //// P(F(), MD("()"))
    //vp.push_back(P(F(1), MD("(2)")) + P(F(2), MD("(1)")) + P(F(1), MD("(0)")));
    //vp.push_back(P(F(1), MD("(1)")) + P(F(1), MD("(0)")));

    //// P(F(1), MD("(, , , )"))
    //vp.push_back(P(F(1), MD("(4, 0, 0, 0)")) - P(F(1), MD("(0, 1, 0, 0)")));
    //vp.push_back(P(F(1), MD("(3, 0, 0, 0)")) - P(F(1), MD("(0, 0, 1, 0)")));
    //vp.push_back(P(F(1), MD("(2, 0, 0, 0)")) - P(F(1), MD("(0, 0, 0, 1)")));

    //// P(F(1), MD("(, )"))
    //vp.push_back(P(F(1), MD("(7, 0)")) - P(F(1), MD("(1, 0)")) - P(F(1), MD("(0, 0)")));
    //vp.push_back(P(F(1), MD("(3, 1)")) - P(F(4), MD("(1, 0)")) + P(F(1), MD("(0, 0)")));

    // P(F(1), MD("(, , )"))
    //vp.push_back(P(F(12), MD("(1, 0, 0)")) - P(F(1), MD("(1, 1, 0)")));
    //vp.push_back(P(F(0, 7), MD("(0, 0, 3)")) - P(F(1), MD("(3, 1, 1)")));
    //vp.push_back(P(F(1), MD("(4, 0, 1)")) - P(F(1), MD("(0, 1, 0)")));

    vp.push_back(P("-x^(1, 1, 0)+(12,0)*x^(1, 0, 0)"));
    vp.push_back(P("-x^(3, 1, 1)+(0,7)*x^(0, 0, 3)"));
    vp.push_back(P("x^(4, 0, 1)-x^(0, 1, 0)"));

    output_aligned(std::cout, vp);
    //std::cout << "\n";

    //typedef Arageli::_Internal::gbasis<P> GB;
    //Arageli::_Internal::gbasis<P> gb(vp);

    std::cout << "\n\n";

    output_aligned(std::cout, groebner_basis(vp));
    //for(GB::iterator i = gb.begin(); i != gb.end(); ++i)
    //    std::cout << i->first << " : " << i->second/P(i->second.leading_coef(), MD("(0, 0, 0)")) << "\n";
}


void test4_59 ()
{
    typedef sparse_polynom<double, double> PD;
    //PD p1 = "x^14-2.8+3*x^-0.3", p2 = "5*x^7.2+2.4*x^-3.1", p3 = "x^550";
    PD p1 = "x^2-1", p2 = "x-2", p3 = "x^50";
    p1 *= p3;
    std::cout
        << "p1 = " << p1 << "\n"
        << "p2 = " << p2 << "\n"
        << "p1 + p2 = " << p1 + p2 << "\n"
        << "p1 * p2 = " << p1 * p2 << "\n"
        << "p1 / p2 = " << p1 / p2 << "\n";
    PD p4 = p1/p2;
    std::cout << "\n\n";
    std::copy(p4.degrees_begin(), p4.degrees_end(), std::ostream_iterator<double>(std::cout, "\n"));
    std::cout << "\n\n";
    vector<double> degrees(p4.size(), p4.degrees_begin(), fromseq);
    std::cout << degrees;
    std::cout << "\n\n";
    std::adjacent_difference(degrees.begin(), degrees.end(), degrees.begin());
    std::cout << "\n\n";
    std::cout << degrees;
}


void test4_60 ()
{
    typedef Arageli::rational<> R;
    typedef Arageli::rational<sparse_multipolynom<R> > F;
    typedef sparse_multipolynom<F> P;

    Arageli::vector<P> vp;

    //Arageli::vector<P> vp =
    //"("
    //    "-x^(1, 1, 0)+(12,0)*x^(1, 0, 0),"
    //    "-x^(3, 1, 1)+(0,7)*x^(0, 0, 3)),"
    //    "x^(4, 0, 1)-x^(0, 1, 0)"
    //")";

    std::cin >> vp;

    output_aligned(std::cout, vp);
    std::cout << "\n\n";
    output_aligned(std::cout, groebner_basis(vp));
}


void test4_61 ()
{
    typedef Arageli::rational<> R;
    typedef sparse_multipolynom<R> P;

    for(;;)
    {
        Arageli::vector<P> vp;
        cin.clear();
        cin >> vp;
        if(vp.is_empty())break;
        output_aligned(std::cout, vp);
        cout << "BG = \n";
        output_aligned(std::cout, groebner_basis(vp));
    }
}


void test4_62 ()
{
    big_int a1 = "0", a2 = "-123", a3 = "292345375937495234088340238", a4 = "-9238402384028420394";
    std::cout << a1 << '\n';
    std::cout << a2 << '\n';
    std::cout << a3 << '\n';
    std::cout << a4 << '\n';

    std::ostringstream buf;
    output_binary(buf, a1);
    output_binary(buf, a2);
    output_binary(buf, a3);
    output_binary(buf, a4);

    std::cout << "\n" << buf.str() << "\n\n";

    big_int b1, b2, b3, b4;

    std::istringstream ibuf(buf.str());
    input_binary(ibuf, b1);
    input_binary(ibuf, b2);
    input_binary(ibuf, b3);
    input_binary(ibuf, b4);

    std::cout << b1 << '\n';
    std::cout << b2 << '\n';
    std::cout << b3 << '\n';
    std::cout << b4 << '\n';
}

void test4_63 ()
{
    sparse_polynom<rational<> > p = "x^2-2*x+1";
    std::cout << "p(x) = " << p;
    std::cout << "\np(x^-1+1) = " << p.subs(sparse_polynom<rational<> >("x^-1+1"));
}

void test4_64 ()
{
    const int m = 1000, n = 10;

    std::ofstream test("skeleton-test-1.ine");
    test << m << ' ' << n << '\n';
    for(int i = 0; i < m; ++i)
    {
        test << std::rand()/10 << ' ';
        for(int j = 1; j < n; ++j)
            test << std::rand() - RAND_MAX/2 << ' ';
        test << '\n';
    }
}




void test4_65 ()
{
    // Magic squares.

    const int size = 5;
    const int nvars = size*size;

    // Build the system.

    typedef big_int T;

    matrix<T> ineq(nvars + 1 + 2*(2*size + 2), nvars + 1 + 1, fromsize);

    // basic inequalities x_i >= 1
    ineq(0, 0) = 1;
    for(int i = 1; i <= nvars; ++i)
    {
        ineq(i, 0) = -1;
        ineq(i, i) =  1;
    }

    // horizontal and vertical magic constraints
    for(int i = 0; i < size; ++i)
    {
        // horizontal:
        for(int j = 0; j < size; ++j)
        {
            ineq(nvars+1 + 4*i + 0, 1 + i*size + j) =  1;
            ineq(nvars+1 + 4*i + 1, 1 + i*size + j) = -1;
        }
        ineq(nvars+1 + 4*i + 0, nvars + 1) = -1;
        ineq(nvars+1 + 4*i + 1, nvars + 1) =  1;

        // vertical:
        for(int j = 0; j < size; ++j)
        {
            ineq(nvars+1 + 4*i + 2, 1 + j*size + i) =  1;
            ineq(nvars+1 + 4*i + 3, 1 + j*size + i) = -1;
        }
        ineq(nvars+1 + 4*i + 2, nvars + 1) = -1;
        ineq(nvars+1 + 4*i + 3, nvars + 1) =  1;
    }

    // diagonal constraints
    for(int i = 0; i < size; ++i)
    {
        ineq(nvars+1 + 4*size + 0, 1 + i*size + i) =  1;
        ineq(nvars+1 + 4*size + 1, 1 + i*size + i) = -1;

        ineq(nvars+1 + 4*size + 2, 1 + i*size - i + size - 1) =  1;
        ineq(nvars+1 + 4*size + 3, 1 + i*size - i + size - 1) = -1;
    }
    ineq(nvars+1 + 4*size + 0, nvars + 1) = -1;
    ineq(nvars+1 + 4*size + 1, nvars + 1) =  1;
    ineq(nvars+1 + 4*size + 2, nvars + 1) = -1;
    ineq(nvars+1 + 4*size + 3, nvars + 1) =  1;

    std::ofstream file("magicNxN.ine");
    file << ineq.nrows() << " " << ineq.ncols() << "\n";
    output_aligned(file, ineq, "", "", " ");

    //Arageli::cone<> ascone(ineq, fromineq);
    //ascone.normalize_all();
    //std::cout << "\nspace = " << ascone.space_dim();
    //std::cout << "\ndim = " << ascone.dim();
    //std::cout << "\nvertex count = " << ascone.min_generatrix().nrows();
    //output_aligned(std::cout << "\n", ascone.min_generatrix());
}



int test4 (int argc, const char** argv)
{
    //test4_1(task);
    //test4_2();
    //test4_3();
    //test4_4();
    //test4_6();
    //test4_7();
    //test4_8();
    //test4_9();
    //test4_10 ();
    //test4_11 ();
    //test4_12 ();
    //test4_13 ();
    //test4_14();
    //test4_15();
    //test4_16();
    //test4_17();
    //test4_18();
    //test4_19();
    //test4_20();
    //test4_21();
    //test4_22();
    //test4_23();
    //test4_24();
    //test4_25();
    //test4_26();
    //test4_27();
    //test4_28();
    //test4_29();
    //test4_30();
    //test4_32();
    //test4_33();
    //test4_34();
    //test4_35();
    //test4_36();
    //test4_37();
    //test4_38();
    //test4_39();
    //test4_40();
    //test4_41();
    //test4_42();
    //test4_43();
    //test4_44();
    //test4_45();
    //test4_46();
    //test4_47();
    //test4_48();
    //test4_49();
    //test4_50();
    //test4_51();
    //test4_52();
    //test4_53();
    //test4_54();
    //test4_55();
    //test4_56();
    //test4_57();
    //test4_58();
    //test4_59();
    //test4_60();
    //test4_61();
    //test4_62();
    //test4_64();
    test4_65();

    return 0;
}
