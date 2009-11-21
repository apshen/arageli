#include "config.hpp"

#if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE) ||    \
    defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_SC_REDUCTION)


#include "vector.hpp"
#include "matrix.hpp"
#include "sc_reduction.hpp"

namespace Arageli
{

template <typename B_type, typename H_type, typename c_type>
bool sc_reduction
(
    B_type &B,
    H_type &H,
    c_type c
)
{
    typedef typename B_type::difference_type index;
    typedef typename B_type::value_type T;

    index m = B.nrows();
    index n = B.ncols();

    bool res_ok = false;
    matrix<T, false> Bst(m, n, fromsize);
    matrix<T, false> Mu(n);
    matrix<T, false> P(n, eye);

    Bst.assign_col(0, B.copy_col(0));
    Mu(0, 0) = dotprod(B.copy_col(0), B.copy_col(0));
    H = H_type(n, eye);

    for(index i = 1; i < n; ++i)
    {
        Bst.assign_col(i, B.copy_col(i));
        for(index j = 0; j < i; ++j)
        {
            Mu(i, j) = dotprod(B.copy_col(i), Bst.copy_col(j))/
                Mu(j, j);
            Bst.assign_col(i, Bst.copy_col(i) -
                Mu(i, j)*Bst.copy_col(j));
        }
        Mu(i, i) = dotprod(Bst.copy_col(i), Bst.copy_col(i));
        if (Mu(i, i) == 0)
        {
            // Basic vectors are linear dependent
            return res_ok;
        }
    }

    T v_1_2(1, 2);
    index k_max = 0;
    index k;
    while(res_ok == false)
    {
        // Size reduction
        for(index s = 1; s < n; ++s)
        {
            for(index i = 0; i < n-s; ++i)
            {
                P(i, i+s) = Mu(i+s, i);
                for(index k = i+1; k < i+s; ++k)
                {
                    P(i, i+s) = P(i, i+s) +
                        Mu(k, i)*P(k, i+s);
                }
                P(i, i+s) = -floor(P(i, i+s) + v_1_2);
            }
        }

        B = B * P;
        H = H * P;

        // Re-counting Mu matrix after size-reduction
        matrix<T, false> Mu_new(n);

        for(index i = 0; i < n; ++i)
        {
            Mu_new(i, i) = Mu(i, i);
            for(index j = i+1; j < n; ++j)
            {
                T tmp = Mu(j, i);
                for(index k = i+1; k < j; ++k)
                {
                    tmp = tmp + Mu(k, i)*P(k, j);
                }
                Mu_new(j, i) = P(i, j) + tmp;
            }
        }
        Mu = Mu_new;

        k = 1;
        while((k < n) && (Mu(k-1, k-1) <= c*Mu(k, k)))
        {
            k = k+1;
        }

        if(k < n)
        {
            // Interchanging k-1 and k columns in B
            B.swap_cols(k-1, k);
            H.swap_cols(k-1, k);
            if(k > 1)
            {
                for(index j = 0; j < k-1; ++j)
                {
                    T tmp = Mu(k, j);
                    Mu(k, j) = Mu(k-1, j);
                    Mu(k-1, j) = tmp;
                }
            }
            T mu = Mu(k, k-1);
            T b2 =  Mu(k, k) + mu*mu*Mu(k-1, k-1);
            Mu(k,k-1) = mu*Mu(k-1, k-1)/b2;
            vector<T, false> b = Bst.copy_col(k-1);
            Bst.assign_col(k-1, Bst.copy_col(k) + mu*b);
            Bst.assign_col(k, (Mu(k, k)/b2)*b -
                Mu(k,k-1)*Bst.copy_col(k));
            Mu(k, k) = Mu(k, k)*Mu(k-1, k-1)/b2;
            Mu(k-1, k-1) = b2;

            for(index i = k+1; i < n; ++i)
            {
                T t = Mu(i,k);
                Mu(i,k) = Mu(i,k-1) - mu*t;
                Mu(i,k-1) = t + Mu(k,k-1)*Mu(i,k);
            }
        }
        else
        {
            res_ok = true;
        }
    }

    return res_ok;
}

} // namespace Arageli

#endif // #ifndef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE


