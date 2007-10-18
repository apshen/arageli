/*****************************************************************************

    hgcd.cpp

    This file is a part of the Arageli library.

    Copyright (C) 2007 Sergey V. Lobanov
    REFERENCE ADDITIONAL COPYRIGHTS HERE
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

/**
    \file hgcd.cpp
    \brief The hgcd.hpp file stuff implementation.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#include "config.hpp"

#if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE) ||    \
    defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE__hgcd)

// REFERENCE ADDITIONAL HEADERS HERE

#include "hgcd.hpp"
#include "factory.hpp"


namespace Arageli
{


//x and y can be changed
template<typename X,typename Y,typename M>
void emgcd_0(X& x, Y& y,M& u1,M& u2,M& v1,M& v2,M& w1,M& w2)
{
    typedef typename X::degree_type X_deg;
    typedef typename Y::degree_type Y_deg;
    ////if(prn) std::cout<<"x="<<x<<"\ny="<<y<<"\n";

    X_deg y_deg=y.degree();
    X_deg x_deg=x.degree();
    //* how to proof correctness with this modification?
    if(is_null(y_deg)&&(x_deg<=2))//with "&&x_deg<=2" algorithm is correct
    {
        //std::cout<<" | x="<<x<<" | y="<<y<<"\n";
        u1=unit<M>();//1*x^0
        u2=null<M>();//0*x^0
        w1=null<M>();
        v1=inverse(y);//(1/y)*x^0
        v2.swap(-x);//v2=-x;//is it optimal construction?
        w2.swap(y);
        return;
    }
    /**/
    ////if(prn) std::cout<<"Start emgcd.x="<<x<<" y="<<y<<"\n";
    ARAGELI_ASSERT_0(x_deg>y_deg);//y=0 is correct
    X_deg m=(x_deg/2)+(x_deg%2);

    if(y_deg<m/*||is_opposite_unit(y_deg)*/)
    {
        //exit
        u1.swap(x);//u1=x;
        u2.swap(y);
        w1=unit<M>();
        w2=null<M>();
        v1=null<M>();
        v2=unit<M>();
        //return unit matrix
    }
    else
    {
        //if(prn) std::cout<<"(";
        X::monom x_m(unit<X::monom::coef_type>(),m);//x^m
        //X b0=x/x_m;//TODO: optimize it!
        X c0=x%x_m;//TODO: optimize it!
        x/=x_m;//TODO: optimize it!
        //Y b1=y/x_m;
        Y c1=y%x_m;//TODO: optimize it!
        y/=x_m;//TODO: optimize it!
        M u11,u12,v11,v12,w11,w12;
        //if(prn) std::cout<<")\n";
        emgcd_0(x,y,u11,u12,v11,v12,w11,w12);//(x,y) was (b0,b1)
        //if(prn) std::cout<<"&";
        u11*=x_m;//TODO: optimize it!
        u11+=w11*c0;
        u11+=v11*c1;
        u12*=x_m;//TODO: optimize it!
        c0*=w12;//TODO: maybe change to u12+=w12*c0,...
        c1*=v12;
        u12+=c0;
        u12+=c1;
        //M e=u12;
        M::degree_type e_deg=u12.degree();
        //if(prn) std::cout<<"@\n";
        if(e_deg<m)
        {
            u1.swap(u11);//u1=u11;
            u2.swap(u12);//u2=u12;
            v1.swap(v11);//v1=v11;
            v2.swap(v12);//v2=v12;
            w1.swap(w11);//w1=w11;
            w2.swap(w12);//w2=w12;
        }
        else
        {
            //if(prn) std::cout<<"$";
            M q=u11/u12;//quotient
            u11-=q*u12;
            //M f=u11-q*u12;//remainder
            X_deg k=m*2-e_deg;
            X::monom x_k(unit<X::monom::coef_type>(),k);
            M g0=u12/x_k;//TODO: optimize it!
            M g1=u11/x_k;//TODO: optimize it!
            //M h0=u12%x_k;//TODO: optimize it!
            //M h1=u11%x_k;//TODO: optimize it!
            u12%=x_k;//TODO: optimize it!
            u11%=x_k;//TODO: optimize it!
            M u21,u22,v21,v22,w21,w22;
            //if(prn) std::cout<<"?\n";
            emgcd_0(g0,g1,u21,u22,v21,v22,w21,w22);
            //writing results: matrix [w21 w22; v21 v22]*[0 1;1 -q]*[w11 w12;v21 v22]
            //if(prn) std::cout<<"!";
            u21*=x_k;//TODO: optimize it!
            u22*=x_k;//TODO: optimize it!
            u21+=w21*u12;
            u21+=v21*u11;
            u1.swap(u21);//u1=u21;
            u12*=w22;
            u11*=v22;
            u22+=u12;
            u22+=u11;
            u2.swap(u22);//u2=u22;
            //M w11mqw12=w11-q*w12;
            //M v11mqv12=v11-q*v12;
            w11-=q*w12;
            v11-=q*v12;
            //
            w1=w21*w12;
            w1+=v21*w11;
            //
            w12*=w22;
            w11*=v22;
            w11+=w12;
            w2.swap(w11);//w2=w11;
            //
            v1=w21*v12;
            v1+=v21*v11;
            //
            v12*=w22;
            v11*=v22;
            v11+=v12;
            v2=v11;
            v2.swap(v11);//v2=v11;
            //if(prn) std::cout<<"%\n";
            //debug in release mode
            //__asm int 3
            //exit
        }
    }
    ////if(prn) std::cout<<"Finish emgcd.x="<<x<<" y="<<y<<" "<<" u1="<<u1<<" u2="<<u2<<"\n";
    ARAGELI_ASSERT_1((u1.degree()>=m)&&(m>u2.degree()));//additional verification for algorithm correctness
}

//x and y can be changed
template<typename X,typename Y,typename M>
void emgcd_0(X& x, Y& y,M& u1,M& u2)//truncated version
{
    typedef typename X::degree_type X_deg;
    typedef typename Y::degree_type Y_deg;
    ////if(prn) std::cout<<"x="<<x<<"\ny="<<y<<"\n";

    X_deg y_deg=y.degree();
    X_deg x_deg=x.degree();
    //* how to proof correctness with this modification?
    if(is_null(y_deg)&&(x_deg<=2))//with "&&x_deg<=2" algorithm is correct
    {
        //std::cout<<" | x="<<x<<" | y="<<y<<"\n";
        u1=unit<M>();//1*x^0
        u2=null<M>();//0*x^0
        return;
    }
    /**/
    ////if(prn) std::cout<<"Start emgcd.x="<<x<<" y="<<y<<"\n";
    ARAGELI_ASSERT_0(x_deg>y_deg);//y=0 is correct
    X_deg m=(x_deg/2)+(x_deg%2);

    if(y_deg<m/*||is_opposite_unit(y_deg)*/)
    {

        u1.swap(x);//u1=x;
        u2.swap(y);
        //exit
    }
    else
    {
        //if(prn) std::cout<<"(";
        X::monom x_m(unit<X::monom::coef_type>(),m);//x^m
        //X b0=x/x_m;//TODO: optimize it!
        X c0=x%x_m;//TODO: optimize it!
        x/=x_m;//TODO: optimize it!
        //Y b1=y/x_m;
        Y c1=y%x_m;//TODO: optimize it!
        y/=x_m;//TODO: optimize it!
        M u11,u12,v11,v12,w11,w12;
        //if(prn) std::cout<<")\n";
        emgcd_0(x,y,u11,u12,v11,v12,w11,w12);//(x,y) was (b0,b1)
        //if(prn) std::cout<<"&";
        u11*=x_m;//TODO: optimize it!
        u11+=w11*c0;
        u11+=v11*c1;
        //M d=u11;//copy?
        u12*=x_m;//TODO: optimize it!
        c0*=w12;//TODO: maybe change to u12+=w12*c0,...
        c1*=v12;
        u12+=c0;
        u12+=c1;
        //M e=u12;
        M::degree_type e_deg=u12.degree();
        //if(prn) std::cout<<"@\n";
        if(e_deg<m)
        {
            u1.swap(u11);//u1=u11;
            u2.swap(u12);//u2=u12;
        }
        else
        {
            //if(prn) std::cout<<"$";
            M q=u11/u12;//quotient
            u11-=q*u12;
            //M f=u11-q*u12;//remainder
            X_deg k=m*2-e_deg;
            X::monom x_k(unit<X::monom::coef_type>(),k);
            M g0=u12/x_k;//TODO: optimize it!
            M g1=u11/x_k;//TODO: optimize it!
            u12%=x_k;//TODO: optimize it!
            u11%=x_k;//TODO: optimize it!
            M u21,u22,v21,v22,w21,w22;
            ////if(prn) std::cout<<"?\n";
            emgcd_0(g0,g1,u21,u22,v21,v22,w21,w22);
            //writing results: matrix [w21 w22; v21 v22]*[0 1;1 -q]*[w11 w12;v21 v22]
            ////if(prn) std::cout<<"!";
            u21*=x_k;//TODO: optimize it!
            u22*=x_k;//TODO: optimize it!
            u21+=w21*u12;
            u21+=v21*u11;
            u1.swap(u21);//u1=u21;
            u12*=w22;
            u11*=v22;
            u22+=u12;
            u22+=u11;
            u2.swap(u22);//u2=u22;
            ////if(prn) std::cout<<"%\n";
            //exit
        }
    }
    ////if(prn) std::cout<<"Finish emgcd.x="<<x<<" y="<<y<<" "<<" u1="<<u1<<" u2="<<u2<<"\n";
    ARAGELI_ASSERT_1((u1.degree()>=m)&&(m>u2.degree()));//additional verification for algorithm correctness
}


//p0 and p1 can be changed
template<typename P0,typename P1,typename M>
void gcd_emgcd_0(P0& p0,P1& p1,M& p2)//,M& p3,M& v1,M& v2,M& w1,M& w2)
{

    P0::degree_type p0_deg=p0.degree();
    P1::degree_type p1_deg=p1.degree();

    //std::cout<<"Start gcd_emgcd.p0="<<p0<<" p1="<<p1<<"\n";
    ARAGELI_ASSERT_0(p0_deg>p1_deg);//y=0 is correct
    //step 1: call emgcd
    M p3;
    emgcd_0(p0,p1,p2,p3);//main action
    //step 2: reserved
    //step 3: check if p3==0 then exit (emgcd_0 finished calculating gcd), else perform one step of the Euclidean algorithm
    if(!is_null(p3))//else exit
    {
        //std::cout<<"gcd_emgcd.Perform a step of the Euclidian Algorithm. p2="<<p2<<" p3="<<p3<<"\n";
        M q=p2/p3;//it's one step of the Euclidian algorithm
        p2-=q*p3;
        //step 4: if p2==0 then finish, else recursively call gcd_emgcd_0
        if(is_null(p2))
            p3.swap(p2);//not copy, just link
        else
        {
            gcd_emgcd_0(p3,p2,p1);//recursion.
            p2.swap(p1);
        }
    }
}

//it's wrapper. How to return value by reference?
//
template<typename P>
P gcd_emgcd(P p0,P p1)
{
    //std::cout<<" * ";
    P p2;//,p3,v1,v2,w1,w2;
    P::degree_type d0=p0.degree(),d1=p1.degree();
    if(d1<d0)
        gcd_emgcd_0(p0,p1,p2);//,p3,v1,v2,w1,w2);
    else if(d1>d0)
        gcd_emgcd_0(p1,p0,p2);//,p3,v1,v2,w1,w2);
    else//d1==d0
    {
        P tmp=p1%p0;//it is not cheap operation
        if(is_null(tmp))
        {
            return p1/(p1.leading_coef_cpy());
            //std::cout<<" # ";
        }
        else
            gcd_emgcd_0(p1,tmp,p2);//,p3,v1,v2,w1,w2);
    }

    ARAGELI_ASSERT_1(is_null(p3));
    //std::cout<<"!!!"<<p2<<"\n";
    p2/=p2.leading_coef_cpy();
    //std::cout<<"!!!!"<<p2<<"\n";
    //std::cout<<" ^ ";
    return p2;

}



// PLACE ALL TEMPLATE NOT INLINE IMPLEMENTATIONS HERE

}


#else    // #if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE) || ...


namespace Arageli
{

// PLACE ALL NOT TEMPLATE NOT INLINE IMPLEMENTATIONS HERE

}


#endif    // #if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE) || ...
