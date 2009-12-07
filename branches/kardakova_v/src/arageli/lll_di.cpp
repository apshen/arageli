#include "config.hpp"

#if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE) ||    \
    defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_LLL_DI)


#include <cmath>

#include "vector.hpp"
#include "matrix.hpp"
#include "lll_di.hpp"

namespace Arageli
{

// The following 3 functions are used only internally in lll_di_reduction

/*
    lll_di_size_reduction routine reduces Gram-Shmidt coeffisient Mufl(k,j)
    to make it grater than -1/2 and less than 1/2

    This routine chenges vectors of the lattice B correspondingly
*/
template <typename B_type, typename H_type>
inline bool lll_di_size_reduction (
    int k,
    int j,
    B_type& B,
    H_type& H,
    // See review note from Sergey S. Lyalin below.
    matrix<double, false>& Bfl,
    // See review note from Sergey S. Lyalin below.
    matrix<double, false>& Mufl,
    int tau
)
{
    typedef typename B_type::difference_type index;

    int mu;
    bool Fc = false;

    mu = floor(Mufl(k, j) + 0.5);
    if (mu > pow(2.0, (double)tau/2.0)) { Fc = true; }
    for (index i = 0; i < j; ++i)
    {
        Mufl(k,i) -= mu*Mufl(j,i);
    }
    Mufl(k, j) -= mu;
    B.assign_col(k, B.copy_col(k) - (mu*B.copy_col(j)));
    H.assign_col(k, H.copy_col(k) - (mu*H.copy_col(j)));
    Bfl.assign_col(k, B.copy_col(k));

    return Fc;
}

template <typename B_type>
inline void lll_di_rotate(int i, int k, B_type& B)
{
    B.insert_col(i, B.copy_col(k));
    B.erase_cols(k+1, 1);
}

/*
    lll_di_deep_insertion routine performs "deep insertion" step
    of the algorithm
*/
template <typename B_type, typename H_type, typename c_type>
inline bool lll_di_deep_insertion (
    int& k,
    B_type& B,
    H_type& H,
    c_type delta,
    // See review note from Sergey S. Lyalin below.
    matrix<double, false>& Bfl,
    // See review note from Sergey S. Lyalin below.
    matrix<double, false>& Mufl
)
{
    typedef typename B_type::difference_type index;
    double c = dotprod(Bfl.copy_col(k), Bfl.copy_col(k));
    index i = 0;
    while (i < k)
    {
        if (delta*Mufl(i, i) < c)
        {
            c -= Mufl(k, i)*Mufl(k, i)*Mufl(i, i);
            i++;
        }
        else
        {
            // place vector k on i-th position.
            lll_di_rotate(i, k, B);
            lll_di_rotate(i, k, H);
            lll_di_rotate(i, k, Bfl);
            if (i > 0)
            {
                k = i-1;
                c = dotprod(Bfl.copy_col(k), Bfl.copy_col(k));
                i = 0;
                continue;
            }
            else
            {
                k = 1;
                return true;
            }
        }
    }
    k++;
    return false;
}


/*
    lll_di_reduction routine implements LLL reduction algorithm with
    "deep insertions".

    Parameters:
    B - matrix containing basis that will be reduced. Basis vectors are
        stored in columns. This matrix will be modified in the course of
        algorithm.
    H - transformation matrix that transforms the old basis into the new one.
    c - parameter of the reduction. Must be grater or equal to 0.5.
        Otherwise the algorithm could fall into infinite loop.
*/
template <typename B_type, typename H_type, typename c_type>
bool lll_di_reduction (B_type& B, H_type& H, c_type c)
{
    typedef typename B_type::difference_type index;
    typedef typename B_type::value_type T;

    index n = B.ncols();
    index m = B.nrows();

    // --------------------------------------------------------
    // Review note from Sergey S. Lyalin
    //
    // The following two definition introduces `double' as
    // a type for intermediate calculations. Please consider
    // to use such type which you can pass through template
    // arguments for lll_di_reduction function. It should be
    // passed through additional argument that is called
    // configurator. Please isolate all separate cases of
    // that usage to one or several typedefs inside the file
    // and then work with me to implement it through
    // Arageli::configurator.
    // Desirability of this change was confirmed with Zolotykh.
    // --------------------------------------------------------

    // Bfl is a floating point approximation to matrix B
    matrix<double, false> Bfl(m, n, fromsize);
    // Mufl represents floating point approximation to the coefficients of
    // Gram-Shmidt procedure for matrix B
    matrix<double, false> Mufl(n);

    c_type delta = (4 + c)/(4*c);
    int tau = 52; //prescision bits in double
    int k = 1;
    bool Fc = false;

label:
    while (k < n)
    {
        //      Initialization:
        // Make floating point approximation of vectors B[0],...,B[k]
        for (index i = 0; i <= k; ++i)
        {
            Bfl.assign_col(i, B.copy_col(i));
        }
        // Compute values that will be used on k-th iteration
        Mufl(k, k) = dotprod(Bfl.copy_col(k), Bfl.copy_col(k));
        if (k == 1) { Mufl(0, 0) = dotprod(Bfl.copy_col(0), Bfl.copy_col(0)); }
        for (index j = 0; j < k; ++j)
        {
            double dot_Bfl_k_j = dotprod(Bfl.copy_col(k), Bfl.copy_col(j));
            double len_Bfl_k_j = sqrt(dotprod(Bfl.copy_col(k), Bfl.copy_col(k)) *
                dotprod(Bfl.copy_col(j), Bfl.copy_col(j)));
            double s;
            if (fabs(dot_Bfl_k_j) < pow(2.0, -((double)tau/2.0))*len_Bfl_k_j )
            {
                // If computational errors are too large then compute
                // coefficients of Gram-Shmidt ortogonalization more
                // precisely
                s = dotprod(B.copy_col(k), B.copy_col(j));
            }
            else
            {
                s = dot_Bfl_k_j;
            }
            Mufl(k, j) = s;
            for (index i = 0; i < j; ++i)
            {
                Mufl(k, j) -= Mufl(j, i)*Mufl(k, i)*Mufl(i, i);
            }
            Mufl(k, j) /= Mufl(j, j);
            Mufl(k, k) -= Mufl(k, j)*Mufl(k, j)*Mufl(j, j);
            if (Mufl(k, k) == 0.0)
            {
                // Basic vectors are linear dependent
                return false;
            }
        }

        // Size-reduction:
        for (index j = k-1; j >= 0; --j)
        {
            if (fabs(Mufl(k,j)) > 0.5)
            {
                Fc = lll_di_size_reduction(k, j, B, H, Bfl, Mufl, tau);
            }
        }

        // If computational errors are too large perform additional step to
        // reduce them
        if (Fc)
        {
            Fc = false;
            k = (k-1 < 1) ? 1 : k-1;
            goto label;
        }

        // Perform "deep insertion" to find shorter vector of the lattice
        if (lll_di_deep_insertion(k, B, H, delta, Bfl, Mufl))
        {
            goto label;
        }
    } // while (k < n)

    return true;
}

} // namespace Arageli

#endif // #ifndef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE

