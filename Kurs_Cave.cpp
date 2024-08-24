#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include<limits>
#include <fstream>
using namespace std;

template<typename T>
T lenght(vector<T> a, vector<T> b)
{
    return sqrt(((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2])));
}

template<typename T>
T lenghtV(vector<T> a, vector<T> b)
{
    return sqrt(((b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2])));
}

template<typename T>
T lenghtS(vector<T> a, vector<T> b)
{
    return sqrt(((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1])));
}


template<typename T>
vector < vector <vector<T>>> Dots_interpenter(vector < vector <vector<T>>>& picket)
{
    vector < vector < vector <T> > > Dot;
    vector < vector <T> > Tdot;
    vector<T> dot{ 0, 0, 0 };
    vector<T> dot1{ 0, 0, 0 };
    vector<T> dot2{ 0, 0, 0 };
    T in_gr = M_PI / 180;

 

    for (int i = 0; i < size(picket); i++)
    {
        dot2[0] +=  picket[i][0][0] * cos(picket[i][0][1] * in_gr) * cos(picket[i][0][2] * in_gr);
        dot2[1] +=  picket[i][0][0] * sin(picket[i][0][1] * in_gr) * cos(picket[i][0][2] * in_gr);
        dot2[2] +=  picket[i][0][0] * sin((picket[i][0][2] * in_gr));


        for (int j = 1; j < size(picket[i]); j++)
        {
            dot[0] = dot2[0] + picket[i][j][0] * cos(picket[i][j][1] * in_gr) * cos(picket[i][j][2] * in_gr);
            dot[1] = dot2[1] + picket[i][j][0] * sin(picket[i][j][1] * in_gr) * cos(picket[i][j][2] * in_gr);
            dot[2] = dot2[2] + picket[i][j][0] * sin((picket[i][j][2] * in_gr));
            Tdot.push_back(dot);
        }
        Dot.push_back(Tdot);
        Tdot.clear();
        dot1 = dot2;
    }



    return Dot;
}


template<typename T>
T area_of_triangle_Ger(vector<T> dot1, vector<T> dot2, vector<T> dot3)
{
    T a, b, c, p;
    a = sqrt((dot1[0] - dot2[0]) * (dot1[0] - dot2[0]) + (dot1[1] - dot2[1]) * (dot1[1] - dot2[1]));

    b = sqrt((dot2[0] - dot3[0]) * (dot2[0] - dot3[0]) + (dot2[1] - dot3[1]) * (dot2[1] - dot3[1]));

    c = sqrt((dot1[0] - dot3[0]) * (dot1[0] - dot3[0]) + (dot1[1] - dot3[1]) * (dot1[1] - dot3[1])
);

    p = (a + b + c) * 0.5;
    double d = sqrt(p * (p - a) * (p - b) * (p - c));
    return sqrt(p * (p - a) * (p - b) * (p - c));
}

template<typename T>
vector<T> max_elm(vector <vector<T>> dots, int j)
{
    int numb = 0;
    T max = dots[0][j];
    for (int i = 1; i < dots.size(); i++)
    {
        if (dots[i][j] > max)
        {
            numb = i;
            max = dots[i][j];
        }
    }
    return dots[numb];
}

template<typename T>
vector<T> max_elm_exp(vector <vector<T>> dots,  vector<T> a, int j)
{
    int numb = 0;
    T max;
    if(dots[0] != a)
        max = dots[0][j];
    else
        max = dots[1][j];

    for (int i = 1; i < dots.size(); i++)
    {
        if (dots[i][j] > max && dots[i] != a)
        {
            numb = i;
            max = dots[i][j];
        }
    }
    return dots[numb];
}

template<typename T>
vector<T> min_elm(vector <vector<T>> dots, int j)
{
    int numb = 0;
    T min = dots[0][j];
    for (int i = 1; i < dots.size(); i++)
    {
        if (dots[i][j] <  min)
        {
            numb = i;
            min = dots[i][j];
        }
    }
    return dots[numb];
}

template<typename T>
vector<T> min_elm_exp( vector <vector<T>> dots, vector<T> a, int j)
{
    int numb = 0;
    T min;

    if (dots[0] != a)
        min = dots[0][j];
    else
        min = dots[1][j];

    for (int i = 1; i < dots.size(); i++)
    {
        if (dots[i][j] < min && dots[i] != a)
        {
            numb = i;
            min = dots[i][j];
        }
    }
    return dots[numb];
}




template<typename T>
T Cave_S(vector < vector <vector<T>>> Dot)
{
    T S = 0;
    vector<T> a, b,c,d,e,f;

    for (int i = 0; i < Dot.size() - 2; i++)
    {
        a = max_elm(Dot[i], 1);
        b = min_elm(Dot[i], 1);
        c = max_elm(Dot[i+1], 1);
        d = min_elm(Dot[i+1], 1);
        e = min_elm(Dot[i], 0);
        f = min_elm(Dot[i + 1], 0);
        a[0] = b[0] = e[0];
        c[0] = d[0] = f[0];
        S += area_of_triangle_Ger(b, a, d);
        S += area_of_triangle_Ger(d,c, a);

    }
    a = max_elm(Dot[Dot.size() - 2], 1);
    b = min_elm(Dot[Dot.size() - 2], 1);
    c = max_elm(Dot[Dot.size() - 1], 1);
    d = min_elm(Dot[Dot.size() - 1], 1);
    e = min_elm(Dot[Dot.size() - 2], 0);
    f = max_elm(Dot[Dot.size() - 1], 0);
    a[0] = b[0] = e[0];
    c[0] = d[0] = f[0];
    S += area_of_triangle_Ger(b, a, d);
    S += area_of_triangle_Ger(d, c, a);

    return S;
}

template<typename T>
 vector <T>  V_Zero(vector<T> Zero, vector <T>  dots)
{

    for (int k = 0; k < 3; k++)
    {
        dots[k] =   dots[k] - Zero[k];
    }

    return dots;
}

template<typename T>
T Determinant(vector<T> a, vector<T> b, vector<T> c)
{
    T l = abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - b[0] * (a[1] * c[2] - a[2] * c[1]) +
        c[0] * (a[1] * b[2] - a[2] * b[1]));
    return abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - b[0] * (a[1] * c[2] - a[2] * c[1]) +
        c[0] * (a[1] * b[2] - a[2] * b[1]));
}



template<typename T>
bool AreSame(T a, T b)
{
    return fabs(a - b) < 0.0000001;
}

template<typename T>
T Sonor(vector<T> A, vector<T> B)
{
    T temp1, temp2, temp3;
    temp1 = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
    temp2 = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
    temp3 = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
    temp1 = temp1 / (temp2 * temp3);
    if(AreSame(temp1,1.0))
        temp1 = 1.0;
    else if(AreSame(temp1, -1.0))
        temp1 = -1.0;
    return acos(temp1);
    
}

template<typename T>
vector<T> Vector_P(vector<T> A, vector<T> B)
{
    vector<T> Ans{ 0, 0, 0 };
    Ans[0] = A[1] * B[2] - A[2] * B[1];
    Ans[1] = -(A[0] * B[2] - A[2] * B[0]);
    Ans[2] = A[0] * B[1] - A[1] * B[0];
    return Ans;
}


template<typename T>
T  Sorted_dots(vector<T> A, vector<T> B, vector<T> C, vector<T> Cen1, vector<T> Cen2)
{
    vector<T> Zero{ 0, 0, 0 };
    vector<T> T1{ 0, 0, 0 };
    vector<T> T2{ 0, 0, 0 };
    vector<T> T3{ 0, 0, 0 };
    Zero[0] = (Cen1[0] + Cen2[0]) * 0.5;
    Zero[1] = (Cen1[1] + Cen2[1]) * 0.5;
    Zero[2] = (Cen1[2] + Cen2[2]) * 0.5;

    T1[0] = A[0] - Zero[0]; T1[1] = A[1] - Zero[1]; T1[2] = A[2] - Zero[2];
    T2[0] = B[0] - Zero[0]; T2[1] = B[1] - Zero[1]; T2[2] = B[2] - Zero[2];
    T3[0] = C[0] - Zero[0]; T3[1] = C[1] - Zero[1]; T3[2] = C[2] - Zero[2];
    //T1[0] = Zero[0] - A[0]; T1[1] = Zero[1] - A[1]; T1[2] = Zero[2] - A[2];
    //T2[0] = B[0] - A[0]; T2[1] = B[1] - A[1]; T2[2] = B[2] - A[2];
    //T3[0] = C[0] - A[0]; T3[1] = C[1] - A[1]; T3[2] = C[2] - A[2];

    return Determinant(T1, T2, T3);


}

template<typename T>
vector < vector < vector <T> > > Sort_lenghts(vector < vector < vector <T> > > dots , T(*lenght)(vector<T>, vector<T>))
{
    vector <T>  a;
    int n = 0;
    for (int i = 0; i < dots.size(); i++)
    {
        a = dots[i][0];
        for (int j = 1; j < dots[i].size(); j++)
        {

            for (int k = 1; k < dots[i].size()-1; k++)
            {
                if (lenght(dots[i][k], a) > lenght(dots[i][k+1], a))
                {
                    swap(dots[i][k], dots[i][k+1]);
                }
            }
        }

    }
    return dots;
}

template<typename T>
int Sort_one(vector<T> a, vector<T> b, vector<T> c, vector<T> d , vector<T> dot , T(*lenght)(vector<T>, vector<T>))
{
    T S_lenght = lenght(a, dot);
    int number = 0;
    if (lenght(c, dot) < S_lenght)
    {
        S_lenght = lenght(c, dot);
        number = 1;
    }
    if (lenght(b, dot) < S_lenght)
    {
        S_lenght = lenght(b, dot);
        number = 2;
    }
    if (lenght(d, dot) < S_lenght)
    {
        S_lenght = lenght(d, dot);
        number = 3;
    }
    return number;
}

template<typename T>
vector < vector < vector <T> > > dots_original_S(vector < vector < vector <T> > > dots, vector < vector <T> >  Dot)
{
    if (dots[0][0] == dots[3][0])
        dots[3][0] = min_elm_exp(Dot, dots[0][0], 1);

    if (dots[0][0] == dots[1][0])
        dots[1][0] = max_elm_exp(Dot, dots[0][0], 1);

    if (dots[2][0] == dots[3][0])
        dots[3][0] = min_elm_exp(Dot, dots[2][0], 1);

    if (dots[2][0] == dots[1][0])
        dots[1][0] = max_elm_exp(Dot, dots[2][0], 1);

    return dots;
}

template<typename T>
vector < vector < vector <T> > > dots_original_V(vector < vector < vector <T> > > dots, vector < vector <T> >  Dot)
{
    if (dots[0][0] == dots[3][0])
        dots[3][0] = min_elm_exp(Dot, dots[0][0], 1);

    if (dots[0][0] == dots[1][0])
        dots[1][0] = max_elm_exp(Dot, dots[0][0], 1);

    if (dots[2][0] == dots[3][0])
        dots[3][0] = min_elm_exp(Dot, dots[2][0], 1);

    if (dots[2][0] == dots[1][0])
        dots[1][0] = max_elm_exp(Dot, dots[2][0], 1);

    return dots;
}

template<typename T>
T Cave_V(vector < vector <vector<T>>> Dot)
{
    T V = 0;

    vector<T> a, b, c, d;
    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector<T> V1{ 0, 0, 0 }, V2{ 0, 0, 0 };
    vector<T> A{ 0, 0, 0 };
    vector<T> B{ 0, 0, 0 };
    vector<T> C{ 0, 0, 0 };
    vector < vector < vector <T> >> Numbers;
    for (int i = 0; i < Dot.size() - 1; i++)
    {

        if (Dot[i].size() == 4 && Dot[i+1].size() == 4)
        {
            dots1 = { {max_elm(Dot[i], 2)}, {max_elm(Dot[i], 1)}, {min_elm(Dot[i], 2)}, {min_elm(Dot[i], 1)} };
            dots1 = dots_original_V(dots1, Dot[i]);
            dots2 = { {max_elm(Dot[i + 1], 2)}, {max_elm(Dot[i + 1], 1)} , {min_elm(Dot[i + 1], 2)}, {min_elm(Dot[i + 1], 1)} };
            dots2 = dots_original_V(dots2, Dot[i+1]);

            for (int j = 0; j < 4; j++)
            {

                V1[0] = (dots1[0][0][0] + dots1[1][0][0] + dots1[3][0][0] + dots1[2][0][0]) * 0.25; V1[1] = (dots1[0][0][1] + dots1[1][0][1] + dots1[3][0][1] + dots1[2][0][1]) * 0.25; V1[2] = (dots1[0][0][2] + dots1[1][0][2] + dots1[3][0][2] + dots1[2][0][2]) * 0.25;
                V2[0] = (dots2[0][0][0] + dots2[1][0][0] + dots2[3][0][0] + dots2[2][0][0]) * 0.25; V2[1] = (dots2[0][0][1] + dots2[1][0][1] + dots2[3][0][1] + dots2[2][0][1]) * 0.25; V2[2] = (dots2[0][0][2] + dots2[1][0][2] + dots2[3][0][2] + dots2[2][0][2]) * 0.25;

                A = dots1[j][0];
                B = dots2[j][0];
                C = dots1[(j + 1) % 4][0];
                V += Sorted_dots(A, B, C, V1, V2);



                A = dots2[j][0];
                B = dots1[(j + 1) % 4][0];
                C = dots2[(j + 1) % 4][0];
                V += Sorted_dots(A, C, B, V1, V2);
            }

            V += Sorted_dots(dots1[0][0], dots1[3][0], dots1[1][0], V1, V2);
            V += Sorted_dots(dots1[2][0], dots1[3][0], dots1[1][0], V1, V2);

            V += Sorted_dots(dots2[0][0], dots2[3][0], dots2[1][0], V1, V2);
            V += Sorted_dots(dots2[2][0], dots2[3][0], dots2[1][0], V1, V2);

        }
        else
        {
            if (Dot[i].size() > 4)
            {
                Numbers = { {},{},{},{} };
                a = max_elm(Dot[i], 2);
                b = min_elm(Dot[i], 2);
                c = max_elm(Dot[i], 1);
                d = min_elm(Dot[i], 1);
                dots1 = { {a}, {c}, {b}, {d} };
                dots1 = dots_original_V(dots1, Dot[i]);
                a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
                for (int j = 0; j < Dot[i ].size(); j++) 
                {
                    if (Dot[i][j] != a && Dot[i ][j] != b && Dot[i][j] != c && Dot[i][j] != d)
                        Numbers[Sort_one(a, b, c, d, Dot[i][j], lenghtV)].push_back(Dot[i ][j]);

                }

                for (int j = 0; j < Numbers[0].size(); j++)
                {
                    if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
                    {
                        dots1[0].push_back(Numbers[0][j]);
                    }
                    else
                    {
                        dots1[3].push_back(Numbers[0][j]);
                    }
                }

                for (int j = 0; j < Numbers[1].size(); j++)
                {
                    if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
                    {
                        dots1[1].push_back(Numbers[1][j]);
                    }
                    else
                    {
                        dots1[0].push_back(Numbers[1][j]);
                    }
                }

                for (int j = 0; j < Numbers[2].size(); j++)
                {
                    if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
                    {
                        dots1[2].push_back(Numbers[2][j]);
                    }
                    else
                    {
                        dots1[1].push_back(Numbers[2][j]);
                    }
                }

                for (int j = 0; j < Numbers[3].size(); j++)
                {
                    if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
                    {
                        dots1[3].push_back(Numbers[3][j]);
                    }
                    else
                    {
                        dots1[2].push_back(Numbers[3][j]);
                    }
                }

            }
            else
            {
                dots1 = { {max_elm(Dot[i], 2)}, {max_elm(Dot[i], 1)}, {min_elm(Dot[i], 2)}, {min_elm(Dot[i], 1)} };
                dots1 = dots_original_V(dots1, Dot[i]);
            }




            if (Dot[i+1].size() > 4)
            {
                Numbers = { {},{},{},{} };
                a = max_elm(Dot[i + 1], 2);
                b = min_elm(Dot[i + 1], 2);
                c = max_elm(Dot[i + 1], 1);
                d = min_elm(Dot[i + 1], 1);
                dots2 = { {a}, {c}, {b}, {d} };
                dots2 = dots_original_V(dots2, Dot[i + 1]);
                a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
                for (int j = 0; j < Dot[i + 1].size(); j++) 
                {
                    if(Dot[i + 1][j] !=a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                        Numbers[Sort_one(a,b,c,d, Dot[i + 1][j], lenghtV)].push_back(Dot[i + 1][j]);
                    
                }

                for (int j = 0; j < Numbers[0].size(); j++)
                {
                    if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
                    {
                        dots2[0].push_back(Numbers[0][j]);
                    }
                    else
                    {
                        dots2[3].push_back(Numbers[0][j]);
                    }
                }

                for (int j = 0; j < Numbers[1].size(); j++)
                {
                    if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
                    {
                        dots2[1].push_back(Numbers[1][j]);
                    }
                    else
                    {
                        dots2[0].push_back(Numbers[1][j]);
                    }
                }

                for (int j = 0; j < Numbers[2].size(); j++)
                {
                    if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
                    {
                        dots2[2].push_back(Numbers[2][j]);
                    }
                    else
                    {
                        dots2[1].push_back(Numbers[2][j]);
                    }
                }

                for (int j = 0; j < Numbers[3].size(); j++)
                {
                    if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
                    {
                        dots2[3].push_back(Numbers[3][j]);
                    }
                    else
                    {
                        dots2[2].push_back(Numbers[3][j]);
                    }
                }
                Numbers.clear();
               
            }
            else
            {
                dots2 = { {max_elm(Dot[i + 1], 2)}, {max_elm(Dot[i + 1], 1)} , {min_elm(Dot[i + 1], 2)}, {min_elm(Dot[i + 1], 1)} };
                dots2 = dots_original_V(dots2, Dot[i + 1]);
            }

            dots1 = Sort_lenghts(dots1, lenghtV);
            dots2 = Sort_lenghts(dots2, lenghtV);
            
            
            V1[0] = (dots1[0][0][0] + dots1[1][0][0] + dots1[3][0][0] + dots1[2][0][0]) * 0.25; V1[1] = (dots1[0][0][1] + dots1[1][0][1] + dots1[3][0][1] + dots1[2][0][1]) * 0.25; V1[2] = (dots1[0][0][2] + dots1[1][0][2] + dots1[3][0][2] + dots1[2][0][2]) * 0.25;
            V2[0] = (dots2[0][0][0] + dots2[1][0][0] + dots2[3][0][0] + dots2[2][0][0]) * 0.25; V2[1] = (dots2[0][0][1] + dots2[1][0][1] + dots2[3][0][1] + dots2[2][0][1]) * 0.25; V2[2] = (dots2[0][0][2] + dots2[1][0][2] + dots2[3][0][2] + dots2[2][0][2]) * 0.25;
            for (int j = 0; j < 4; j++)
            {
               
                if (dots1[j].size() > 1)
                {
                    
                    for (int k = 0; k < dots1[j].size() - 1; k++)
                    {
                        A = dots1[j][k];
                        B = dots2[j][0];
                        C = dots1[j][k+1];

                        
                        


                        V += Sorted_dots(A,B,C, V1, V2);
                        V += Sorted_dots(A, C, dots1[(j + 1) % 4][0], V1, V2);
                        
                    }
                    A = dots1[j ][dots1[j ].size()-1];
                    B = dots2[j ][0];
                    C = dots1[(j + 1) % 4][0];

                    
                   
                    V += Sorted_dots(A, B, C, V1, V2);
                }
                else
                {
                    A = dots1[j ][0];
                    B = dots2[j ][0];
                    C = dots1[(j + 1) % 4][0];

                   
                    

                    V += Sorted_dots(A, B, C, V1, V2);;
                }

                if (dots2[j].size() > 1)
                {
                    for (int k = 0; k < dots2[j].size() - 1; k++)
                    {
                        A = dots2[j][k];
                        B = dots1[(j + 1) % 4][0];
                        C = dots2[j][k + 1];

                        
                       

                        V += Sorted_dots(A, B, C, V1, V2);
                        V += Sorted_dots(A, C, dots2[(j + 1) % 4][0], V1, V2);

                    }

                }

                A = dots2[j][dots2[j].size()-1];
                B = dots1[(j + 1) % 4][0];
                C = dots2[(j + 1) % 4][0];

                
                

                V += Sorted_dots(A, B, C, V1, V2);

            }


            V += Sorted_dots(dots1[0][0], dots1[3][0], dots1[1][0], V1, V2);
            V += Sorted_dots(dots1[2][0], dots1[3][0], dots1[1][0], V1, V2);

            V += Sorted_dots(dots2[0][0], dots2[3][0], dots2[1][0], V1, V2);
            V += Sorted_dots(dots2[2][0], dots2[3][0], dots2[1][0], V1, V2);
        }
        dots1.clear();
        dots2.clear();

    }
    return V/6;
}

template<typename T>
T Cave_V_BN1(vector < vector <vector<T>>> Dot)
{
    T V = 0;

    vector<T> a, b, c, d;
    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector<T> V1{ 0, 0, 0 }, V2{ 0, 0, 0 };
    vector<T> A{ 0, 0, 0 };
    vector<T> B{ 0, 0, 0 };
    vector<T> C{ 0, 0, 0 };
    vector < vector < vector <T> >> Numbers;
    for (int i = 0; i < Dot.size() - 1; i++)
    {
            dots1 = { {max_elm(Dot[i], 2)}, {max_elm(Dot[i], 1)}, {min_elm(Dot[i], 2)}, {min_elm(Dot[i], 1)} };
            dots1 = dots_original_V(dots1, Dot[i]);
            dots2 = { {max_elm(Dot[i + 1], 2)}, {max_elm(Dot[i + 1], 1)} , {min_elm(Dot[i + 1], 2)}, {min_elm(Dot[i + 1], 1)} };
            dots2 = dots_original_V(dots2, Dot[i + 1]);

            for (int j = 0; j < 4; j++)
            {

                V1[0] = (dots1[0][0][0] + dots1[2][0][0]) * 0.5; V1[1] = (dots1[0][0][1] + dots1[2][0][1]) * 0.5; V1[2] = (dots1[0][0][2] + dots1[2][0][2]) * 0.5;
                V2[0] = (dots2[0][0][0] + dots2[2][0][0]) * 0.5; V2[1] = (dots2[0][0][1] + dots2[2][0][1]) * 0.5; V2[2] = (dots2[0][0][2] + dots2[2][0][2]) * 0.5;

                A = dots1[j][0];
                B = dots2[j][0];
                C = dots1[(j + 1) % 4][0];
                V += Sorted_dots(A, B, C, V1, V2);



                A = dots2[j][0];
                B = dots1[(j + 1) % 4][0];
                C = dots2[(j + 1) % 4][0];
                V += Sorted_dots(A, C, B, V1, V2);
            }

            V += Sorted_dots(dots1[0][0], dots1[3][0], dots1[1][0], V1, V2);
            V += Sorted_dots(dots1[2][0], dots1[3][0], dots1[1][0], V1, V2);

            V += Sorted_dots(dots2[0][0], dots2[3][0], dots2[1][0], V1, V2);
            V += Sorted_dots(dots2[2][0], dots2[3][0], dots2[1][0], V1, V2);

            
        
    }

    return V / 6;
}

template<typename T>
T Cave_V_BN2(vector < vector <vector<T>>> Dot)
{
    T V = 0;

    vector<T> a, b, c, d, e1,e2;
    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector<T> V1{ 0, 0, 0 }, V2{ 0, 0, 0 };
    vector<T> A{ 0, 0, 0 };
    vector<T> B{ 0, 0, 0 };
    vector<T> C{ 0, 0, 0 };
    vector < vector < vector <T> >> Numbers;
    for (int i = 0; i < Dot.size() - 2; i++)
    {
        e1 = min_elm(Dot[i], 0);
        e2 = min_elm(Dot[i+1], 0);
        dots1 = { {max_elm(Dot[i], 2)}, {max_elm(Dot[i], 1)}, {min_elm(Dot[i], 2)}, {min_elm(Dot[i], 1)} };
        dots1 = dots_original_V(dots1, Dot[i]);
        dots2 = { {max_elm(Dot[i + 1], 2)}, {max_elm(Dot[i + 1], 1)} , {min_elm(Dot[i + 1], 2)}, {min_elm(Dot[i + 1], 1)} };
        dots2 = dots_original_V(dots2, Dot[i + 1]);
        for (int j = 0; j < 4; j++)
        {
            dots1[j][0][0] = e1[0];
            dots2[j][0][0] = e2[0];
        }

        for (int j = 0; j < 4; j++)
        {

            V1[0] = (dots1[0][0][0] + dots1[2][0][0]) * 0.5; V1[1] = (dots1[0][0][1] + dots1[2][0][1]) * 0.5; V1[2] = (dots1[0][0][2] + dots1[2][0][2]) * 0.5;
            V2[0] = (dots2[0][0][0] + dots2[2][0][0]) * 0.5; V2[1] = (dots2[0][0][1] + dots2[2][0][1]) * 0.5; V2[2] = (dots2[0][0][2] + dots2[2][0][2]) * 0.5;

            A = dots1[j][0];
            B = dots2[j][0];
            C = dots1[(j + 1) % 4][0];
            V += Sorted_dots(A, B, C, V1, V2);



            A = dots2[j][0];
            B = dots1[(j + 1) % 4][0];
            C = dots2[(j + 1) % 4][0];
            V += Sorted_dots(A, C, B, V1, V2);
        }

        V += Sorted_dots(dots1[0][0], dots1[3][0], dots1[1][0], V1, V2);
        V += Sorted_dots(dots1[2][0], dots1[3][0], dots1[1][0], V1, V2);

        V += Sorted_dots(dots2[0][0], dots2[3][0], dots2[1][0], V1, V2);
        V += Sorted_dots(dots2[2][0], dots2[3][0], dots2[1][0], V1, V2);



    }


    e1 = min_elm(Dot[Dot.size() - 2], 0);
    e2 = max_elm(Dot[Dot.size() - 2 + 1], 0);
    dots1 = { {max_elm(Dot[Dot.size() - 2], 2)}, {max_elm(Dot[Dot.size() - 2], 1)}, {min_elm(Dot[Dot.size() - 2], 2)}, {min_elm(Dot[Dot.size() - 2], 1)} };
    dots1 = dots_original_V(dots1, Dot[Dot.size() - 2]);
    dots2 = { {max_elm(Dot[Dot.size() - 2 + 1], 2)}, {max_elm(Dot[Dot.size() - 2 + 1], 1)} , {min_elm(Dot[Dot.size() - 2 + 1], 2)}, {min_elm(Dot[Dot.size() - 2 + 1], 1)} };
    dots2 = dots_original_V(dots2, Dot[Dot.size() - 2 + 1]);
    
    for (int j = 0; j < 4; j++)
    {
        dots1[j][0][0] = e1[0];
        dots2[j][0][0] = e2[0];
    }

    for (int j = 0; j < 4; j++)
    {

        V1[0] = (dots1[0][0][0] + dots1[2][0][0]) * 0.5; V1[1] = (dots1[0][0][1] + dots1[2][0][1]) * 0.5; V1[2] = (dots1[0][0][2] + dots1[2][0][2]) * 0.5;
        V2[0] = (dots2[0][0][0] + dots2[2][0][0]) * 0.5; V2[1] = (dots2[0][0][1] + dots2[2][0][1]) * 0.5; V2[2] = (dots2[0][0][2] + dots2[2][0][2]) * 0.5;

        A = dots1[j][0];
        B = dots2[j][0];
        C = dots1[(j + 1) % 4][0];
        V += Sorted_dots(A, B, C, V1, V2);



        A = dots2[j][0];
        B = dots1[(j + 1) % 4][0];
        C = dots2[(j + 1) % 4][0];
        V += Sorted_dots(A, C, B, V1, V2);
    }

    V += Sorted_dots(dots1[0][0], dots1[3][0], dots1[1][0], V1, V2);
    V += Sorted_dots(dots1[2][0], dots1[3][0], dots1[1][0], V1, V2);

    V += Sorted_dots(dots2[0][0], dots2[3][0], dots2[1][0], V1, V2);
    V += Sorted_dots(dots2[2][0], dots2[3][0], dots2[1][0], V1, V2);

    return V / 6;
}

template<typename T>
T Cave_S2(vector < vector <vector<T>>> Dot)
{
    T S = 0;
    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector<T> a, b,c,d;
    vector < vector < vector <T> > >  Numbers;
    for (int i = 0; i < Dot.size() - 2; i++)
    {
        a = max_elm(Dot[i], 0);
        b = min_elm(Dot[i], 0);
        c = max_elm(Dot[i], 1);
        d = min_elm(Dot[i], 1);
        dots1 = { {a},{c}, {b}, {d} };

        if (Dot[i].size() > 4)
        {
            Numbers = { {},{},{},{} };

            dots1 = { {a}, {c}, {b}, {d} };
            //dots1 = dots_original_V(dots1, Dot[i]);
            a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
            for (int j = 0; j < Dot[i].size(); j++) 
            {
                if (Dot[i][j] != a && Dot[i][j] != b && Dot[i][j] != c && Dot[i][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i][j], lenghtS)].push_back(Dot[i][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) > lenghtS(Numbers[0][j], d))
                {
                    dots1[3].push_back(Numbers[0][j]);
                }

            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots1[1].push_back(Numbers[1][j]);
                }

            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
                {
                    dots1[2].push_back(Numbers[2][j]);
                }
                else
                {
                    dots1[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots1[3].push_back(Numbers[3][j]);
                }
                else
                {
                    dots1[2].push_back(Numbers[3][j]);
                }
            }
            Numbers.clear();

            dots1 = Sort_lenghts(dots1, lenghtS);
            for (int i = 0; i < dots1[1].size()-1; i++)
            {
                S += area_of_triangle_Ger(c, dots1[1][i+1], dots1[1][i]);
            }
            S += area_of_triangle_Ger(c, b, d);
            for (int i = 0; i < dots1[2].size()-1; i++)
            {
                S += area_of_triangle_Ger(d, dots1[2][i + 1], dots1[2][i]);
            }
        }
        else
        {
            dots1 = { {max_elm(Dot[i], 0)}, {max_elm(Dot[i], 1)} , {min_elm(Dot[i], 0)}, {min_elm(Dot[i ], 1)} };

            S += area_of_triangle_Ger(dots1[1][0], dots1[2][0], dots1[3][0]);
        }
        if (Dot[i + 1].size() > 4)
        {
            a = max_elm(Dot[i + 1], 0);
            b = min_elm(Dot[i + 1], 0);
            c = max_elm(Dot[i + 1], 1);
            d = min_elm(Dot[i + 1], 1);
            dots2 = { {a}, {c}, {b}, {d} };
            dots2 = dots_original_S(dots2, Dot[i+1]);
            a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
            Numbers = { {},{},{},{} };


            for (int j = 0; j < Dot[i + 1].size(); j++) //Переписать
            {
                if (Dot[i + 1][j] != a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i + 1][j], lenghtS)].push_back(Dot[i + 1][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) > lenghtS(Numbers[0][j], d))
                {
                    dots2[3].push_back(Numbers[0][j]);
                }

            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots2[1].push_back(Numbers[1][j]);
                }

            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
                {
                    dots2[2].push_back(Numbers[2][j]);
                }
                else
                {
                    dots2[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots2[3].push_back(Numbers[3][j]);
                }
                else
                {
                    dots2[2].push_back(Numbers[3][j]);
                }
            }

            Numbers.clear();

            dots2 = Sort_lenghts(dots2, lenghtS);
            for (int i = 0; i < dots2[1].size() - 1; i++)
            {
                S -= area_of_triangle_Ger(c, dots2[1][i + 1], dots2[1][i]);
            }
            S -= area_of_triangle_Ger(c, b, d);
            for (int i = 0; i < dots2[2].size() - 1; i++)
            {
                S -= area_of_triangle_Ger(d, dots2[2][i + 1], dots2[2][i]);
            }
          
        }
        else
        {
            dots2 = { {max_elm(Dot[i+1], 1)}, {max_elm(Dot[i + 1], 0)} , {min_elm(Dot[i + 1], 1)}, {min_elm(Dot[i + 1], 0)} };

            S -= area_of_triangle_Ger(dots2[2][0], dots2[0][0], dots2[3][0]);
        }
        S += area_of_triangle_Ger(dots1[1][0], dots1[3][0], dots2[3][0]);
        S += area_of_triangle_Ger(dots2[1][0], dots2[3][0], dots1[1][0]);
    }

    a = max_elm(Dot[Dot.size() - 2], 0);
    b = min_elm(Dot[Dot.size() - 2], 0);
    c = max_elm(Dot[Dot.size() - 2], 1);
    d = min_elm(Dot[Dot.size() - 2], 1);
    dots1 = { {a},{c}, {b}, {d} };
    if (Dot[Dot.size() - 2].size() > 4)
    {
        Numbers = { {},{},{},{} };

        dots1 = { {a}, {c}, {b}, {d} };
        dots1 = dots_original_S(dots1, Dot[Dot.size() - 2]);
        a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
        for (int j = 0; j < Dot[Dot.size() - 2].size(); j++) //Переписать
        {
            if (Dot[Dot.size() - 2][j] != a && Dot[Dot.size() - 2][j] != b && Dot[Dot.size() - 2][j] != c && Dot[Dot.size() - 2][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 2][j], lenghtS)].push_back(Dot[Dot.size() - 2][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtS(Numbers[0][j], c) > lenghtS(Numbers[0][j], d))
            {
                dots1[3].push_back(Numbers[0][j]);
            }

        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
            {
                dots1[1].push_back(Numbers[1][j]);
            }

        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
            {
                dots1[2].push_back(Numbers[2][j]);
            }
            else
            {
                dots1[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
            {
                dots1[3].push_back(Numbers[3][j]);
            }
            else
            {
                dots1[2].push_back(Numbers[3][j]);
            }
        }
        Numbers.clear();

        dots1 = Sort_lenghts(dots1, lenghtS);
        for (int i = 0; i < dots1[1].size() - 1; i++)
        {
            S += area_of_triangle_Ger(c, dots1[1][i + 1], dots1[1][i]);
        }
        S += area_of_triangle_Ger(c, b, d);
        for (int i = 0; i < dots1[2].size() - 1; i++)
        {
            S += area_of_triangle_Ger(d, dots1[2][i + 1], dots1[2][i]);
        }
    }
    else
    {
        dots1 = { {max_elm(Dot[Dot.size() - 2], 0)}, {max_elm(Dot[Dot.size() - 2], 1)} , {min_elm(Dot[Dot.size() - 2], 0)}, {min_elm(Dot[Dot.size() - 2], 1)} };

        S += area_of_triangle_Ger(dots1[1][0], dots1[2][0], dots1[3][0]);
    }

    if (Dot[Dot.size() - 1].size() > 4)
    {
        a = max_elm(Dot[Dot.size() - 1], 0);
        b = min_elm(Dot[Dot.size() - 1], 0);
        c = max_elm(Dot[Dot.size() - 1], 1);
        d = min_elm(Dot[Dot.size() - 1], 1);
        dots2 = { {a}, {c}, {b}, {d} };
        dots2 = dots_original_S(dots2, Dot[Dot.size() - 1]);
        a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];

        Numbers = { {},{},{},{} };


        for (int j = 0; j < Dot[Dot.size() - 1].size(); j++) //Переписать
        {
            if (Dot[Dot.size() - 1][j] != a && Dot[Dot.size() - 1][j] != b && Dot[Dot.size() - 1][j] != c && Dot[Dot.size() - 1][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 1][j], lenghtS)].push_back(Dot[Dot.size() - 1][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtS(Numbers[0][j], c) < lenghtS(Numbers[0][j], d))
            {
                dots2[0].push_back(Numbers[0][j]);
            }
            else
            {
                dots2[3].push_back(Numbers[0][j]);
            }
        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
            {
                dots2[1].push_back(Numbers[1][j]);
            }
            else
            {
                dots2[0].push_back(Numbers[1][j]);
            }
        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtS(Numbers[2][j], d) > lenghtS(Numbers[2][j], c))
            {
                dots2[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
            {
                dots2[3].push_back(Numbers[3][j]);
            }
        }

        Numbers.clear();

        dots2 = Sort_lenghts(dots2, lenghtS);
        for (int i = 0; i < dots2[3].size() - 1; i++)
        {
            S += area_of_triangle_Ger(a, dots2[3][i+ 1], dots2[3][i]);
        }
        S += area_of_triangle_Ger(d, a, b);
        for (int i = 0; i < dots2[0].size() - 1; i++)
        {
            S +=  area_of_triangle_Ger(b, dots2[0][i + 1], dots2[0][i]);
        }
    }
    else
    {
        dots2 = { {max_elm(Dot[Dot.size() - 1], 1)}, {max_elm(Dot[Dot.size() - 1], 0)} , {min_elm(Dot[Dot.size() - 1], 1)}, {min_elm(Dot[Dot.size() - 1], 0)} };
        
        S +=   area_of_triangle_Ger(dots2[1][0], dots2[0][0], dots2[3][0]);
    }
    S += area_of_triangle_Ger(dots1[1][0], dots1[3][0], dots2[3][0]);
    S += area_of_triangle_Ger(dots2[1][0], dots2[3][0], dots1[1][0]);

    return S;
}

template<typename T>
T Vector_P_Nap(vector<T> A, vector<T> B)
{
    T Ans = 0;
    Ans += A[0]*B[1] - B[0] * A[1];
    return Ans;
}

template<typename T>
T Vector_P_D(vector<T> A, vector<T> B)
{
    T Ans = 0;
    Ans += abs(A[0] * B[1] - B[0] * A[1]);
    return Ans;
}

template<typename T>
T Vector_P_Nap_V(vector<T> A, vector<T> B)
{
    T Ans = 0;
    Ans += A[2] * B[1] - B[2] * A[1];
    return Ans;
}

template<typename T>
T Vector_P_D_V(vector<T> A, vector<T> B)
{
    T Ans = 0;
    Ans += abs(A[2] * B[1] - B[2] * A[1]);
    return Ans;
}

template<typename T>
T S_nap(vector < vector <vector<T>>>& dots)
{
    T S = 0;
  

    for (int i = 0; i < dots.size() - 1; i++)
    {
        for (int j = 0; j < dots[i].size() - 1; j++)
        {
            S += Vector_P_Nap(dots[i][j], dots[i][j + 1]);
        }
        S += Vector_P_Nap(dots[i][dots[i].size() - 1], dots[i +1 ][0]);
    }

    for (int j = 0; j < dots[dots.size() - 1].size() - 1; j++)
    {
        S += Vector_P_Nap(dots[dots.size() - 1][j], dots[dots.size() - 1][j + 1]);
    }
    S += Vector_P_Nap(dots[dots.size() - 1][dots[dots.size() - 1].size() - 1], dots[0][0]);
    
    return 0.5*S;
}

template<typename T>
T S_nap2(vector < vector <vector<T>>>& dots, vector<T> V)
{
    T S = 0;

    vector <T> T1{ 0,0 };
    vector <T> T2{ 0,0 };

    for (int i = 0; i < dots.size() - 1; i++)
    {
        for (int j = 0; j < dots[i].size() - 1; j++)
        {
            T1[0] = dots[i][j][0] - V[0]; T1[1] = dots[i][j][1] - V[1];
            T2[0] = dots[i][j + 1][0] - V[0]; T2[1] = dots[i][j + 1][1] - V[1];
            S += Vector_P_D(T1, T2);
        }
        T1[0] = dots[i][dots[i].size() - 1][0] - V[0]; T1[1] = dots[i][dots[i].size() - 1][1] - V[1];
        T2[0] = dots[i + 1][0][0] - V[0]; T2[1] = dots[i + 1][0][1] - V[1];
        S += Vector_P_D(T1, T2);
    }

    for (int j = 0; j < dots[dots.size() - 1].size() - 1; j++)
    {
        T1[0] = dots[dots.size() - 1][j][0] - V[0]; T1[1] = dots[dots.size() - 1][j][1] - V[1];
        T2[0] = dots[dots.size() - 1][j + 1][0] - V[0]; T2[1] = dots[dots.size() - 1][j + 1][1] - V[1];
        S += Vector_P_D(T1, T2);
    }
    T1[0] = dots[dots.size() - 1][dots[dots.size() - 1].size() - 1][0] - V[0]; T1[1] = dots[dots.size() - 1][dots[dots.size() - 1].size() - 1][1] - V[1];
    T2[0] = dots[0][0][0] - V[0]; T2[1] = dots[0][0][1] - V[1];
    S += Vector_P_D(T1, T2);

    return 0.5 * S;
}

template<typename T>
T S_nap_V(vector < vector <vector<T>>>& dots)
{
    T S = 0;

    for (int i = 0; i < dots.size() - 1; i++)
    {
        for (int j = 0; j < dots[i].size() - 1; j++)
        {
            S += Vector_P_Nap_V(dots[i][j], dots[i][j + 1]);
        }
        S += Vector_P_Nap_V(dots[i][dots[i].size() - 1], dots[i + 1][0]);
    }

    for (int j = 0; j < dots[dots.size() - 1].size() - 1; j++)
    {
        S += Vector_P_Nap_V(dots[dots.size() - 1][j], dots[dots.size() - 1][j + 1]);
    }
    S += Vector_P_Nap_V(dots[dots.size() - 1][dots[dots.size() - 1].size() - 1], dots[0][0]);

    return 0.5 * S;
}

template<typename T>
T S_nap2_V(vector < vector <vector<T>>>& dots, vector<T> V)
{
    T S = 0;

    vector <T> T1{ 0,0, 0 };
    vector <T> T2{ 0,0 ,0};

    for (int i = 0; i < dots.size() - 1; i++)
    {
        for (int j = 0; j < dots[i].size() - 1; j++)
        {
            T1[0] = dots[i][j][1] - V[1]; T1[1] = dots[i][j][2] - V[2];
            T2[0] = dots[i][j + 1][1] - V[1]; T2[1] = dots[i][j + 1][2] - V[2];
            S += Vector_P_D_V(T1, T2);
        }
        T1[1] = dots[i][dots[i].size() - 1][1] - V[1]; T1[2] = dots[i][dots[i].size() - 1][2] - V[2];
        T2[1] = dots[i + 1][0][1] - V[1]; T2[2] = dots[i + 1][0][2] - V[2];
        S += Vector_P_D_V(T1, T2);
    }

    for (int j = 0; j < dots[dots.size() - 1].size() - 1; j++)
    {
        T1[1] = dots[dots.size() - 1][j][1] - V[1]; T1[2] = dots[dots.size() - 1][j][2] - V[2];
        T2[1] = dots[dots.size() - 1][j + 1][1] - V[1]; T2[2] = dots[dots.size() - 1][j + 1][2] - V[2];
        S += Vector_P_D_V(T1, T2);
    }
    T1[1] = dots[dots.size() - 1][dots[dots.size() - 1].size() - 1][1] - V[1]; T1[2] = dots[dots.size() - 1][dots[dots.size() - 1].size() - 1][2] - V[2];
    T2[1] = dots[0][0][1] - V[1]; T2[2] = dots[0][0][2] - V[2];
    S += Vector_P_D_V(T1, T2);

    return 0.5 * S;
}

template<typename T>
T Cave_S2_2(vector < vector <vector<T>>> Dot)
{
    T S = 0;
    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector < vector < vector <T> > >  Numbers;
    vector<T> a, b, c, d;
    for (int i = 0; i < Dot.size() - 2; i++)
    {
        a = max_elm(Dot[i], 0);
        b = min_elm(Dot[i], 0);
        c = max_elm(Dot[i], 1);
        d = min_elm(Dot[i], 1);
        dots1 = { {a},{c}, {b}, {d} };
        if (Dot[i].size() > 4)
        {
            Numbers = { {},{},{},{} };

            dots1 = { {a}, {c}, {b}, {d} };
            dots1 = dots_original_S(dots1, Dot[i]);
            a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
            for (int j = 0; j < Dot[i].size(); j++) //Переписать
            {
                if (Dot[i][j] != a && Dot[i][j] != b && Dot[i][j] != c && Dot[i][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i][j], lenghtS)].push_back(Dot[i][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) > lenghtS(Numbers[0][j], d))
                {
                    dots1[3].push_back(Numbers[0][j]);
                }

            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots1[1].push_back(Numbers[1][j]);
                }

            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
                {
                    dots1[2].push_back(Numbers[2][j]);
                }
                else
                {
                    dots1[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots1[3].push_back(Numbers[3][j]);
                }
                else
                {
                    dots1[2].push_back(Numbers[3][j]);
                }
            }
            Numbers.clear();

            dots1 = Sort_lenghts(dots1, lenghtS);
            dots1.erase(dots1.begin());
            dots1[2] = { d };
            S += S_nap(dots1);

        }
        else
        {
            dots1 = { {max_elm(Dot[i], 0)}, {max_elm(Dot[i], 1)} , {min_elm(Dot[i], 0)}, {min_elm(Dot[i], 1)} };

            S += area_of_triangle_Ger(dots1[1][0], dots1[2][0], dots1[3][0]);
            dots1.erase(dots1.begin());
        }
        if (Dot[i + 1].size() > 4)
        {
            a = max_elm(Dot[i + 1], 0);
            b = min_elm(Dot[i + 1], 0);
            c = max_elm(Dot[i + 1], 1);
            d = min_elm(Dot[i + 1], 1);
            dots2 = { {a}, {c}, {b}, {d} };
            dots2 = dots_original_S(dots2, Dot[i + 1]);
            a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];

            Numbers = { {},{},{},{} };


            for (int j = 0; j < Dot[i + 1].size(); j++) //Переписать
            {
                if (Dot[i + 1][j] != a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i + 1][j], lenghtS)].push_back(Dot[i + 1][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) < lenghtS(Numbers[0][j], d))
                {
                    dots2[0].push_back(Numbers[0][j]);
                }
                else
                {
                    dots2[3].push_back(Numbers[0][j]);
                }
            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots2[1].push_back(Numbers[1][j]);
                }
                else
                {
                    dots2[0].push_back(Numbers[1][j]);
                }
            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) > lenghtS(Numbers[2][j], c))
                {
                    dots2[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots2[3].push_back(Numbers[3][j]);
                }
            }

            Numbers.clear();

            dots2 = Sort_lenghts(dots2, lenghtS);
            dots2.erase(dots2.begin());
            dots2[2] = { d };
            S -= S_nap(dots2);
        }
        else
        {
            dots2 = { {max_elm(Dot[i + 1], 1)}, {max_elm(Dot[i + 1], 0)} , {min_elm(Dot[i + 1], 1)}, {min_elm(Dot[i + 1], 0)} };

            S -= area_of_triangle_Ger(dots2[1][0], dots2[2][0], dots2[3][0]);
            dots2.erase(dots2.begin());
        }

        S += area_of_triangle_Ger(dots1[0][0], dots1[2][0], dots2[2][0]);
        S += area_of_triangle_Ger(dots2[0][0], dots2[2][0], dots1[0][0]);

        
    }

    a = max_elm(Dot[Dot.size() - 2], 0);
    b = min_elm(Dot[Dot.size() - 2], 0);
    c = max_elm(Dot[Dot.size() - 2], 1);
    d = min_elm(Dot[Dot.size() - 2], 1);
    dots1 = { {a},{c}, {b}, {d} };
    if (Dot[Dot.size() - 2].size() > 4)
    {
        Numbers = { {},{},{},{} };

        dots1 = { {a}, {c}, {b}, {d} };
        dots1 = dots_original_S(dots1, Dot[Dot.size() - 2]);
        a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
        for (int j = 0; j < Dot[Dot.size() - 2].size(); j++) //Переписать
        {
            if (Dot[Dot.size() - 2][j] != a && Dot[Dot.size() - 2][j] != b && Dot[Dot.size() - 2][j] != c && Dot[Dot.size() - 2][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 2][j], lenghtS)].push_back(Dot[Dot.size() - 2][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtS(Numbers[0][j], c) > lenghtS(Numbers[0][j], d))
            {
                dots1[3].push_back(Numbers[0][j]);
            }

        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
            {
                dots1[1].push_back(Numbers[1][j]);
            }

        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
            {
                dots1[2].push_back(Numbers[2][j]);
            }
            else
            {
                dots1[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
            {
                dots1[3].push_back(Numbers[3][j]);
            }
            else
            {
                dots1[2].push_back(Numbers[3][j]);
            }
        }
        Numbers.clear();

        dots1 = Sort_lenghts(dots1, lenghtS);
        dots1.erase(dots1.begin());
        dots1[2] = { d };
        S += S_nap(dots1);

    }
    else
    {
        dots1 = { {max_elm(Dot[Dot.size() - 2], 0)}, {max_elm(Dot[Dot.size() - 2], 1)} , {min_elm(Dot[Dot.size() - 2], 0)}, {min_elm(Dot[Dot.size() - 2], 1)} };

        S += area_of_triangle_Ger(dots1[1][0], dots1[2][0], dots1[3][0]);
        dots1.erase(dots1.begin());
    }
    if (Dot[Dot.size() - 1].size() > 4)
    {
        a = max_elm(Dot[Dot.size() - 1], 0);
        b = min_elm(Dot[Dot.size() - 1], 0);
        c = max_elm(Dot[Dot.size() - 1], 1);
        d = min_elm(Dot[Dot.size() - 1], 1);
        dots2 = { {a}, {c}, {b}, {d} };
        dots2 = dots_original_S(dots2, Dot[Dot.size() - 1]);
        a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];

        Numbers = { {},{},{},{} };


        for (int j = 0; j < Dot[Dot.size() - 1].size(); j++) //Переписать
        {
            if (Dot[Dot.size() - 1][j] != a && Dot[Dot.size() - 1][j] != b && Dot[Dot.size() - 1][j] != c && Dot[Dot.size() - 1][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 1][j], lenghtS)].push_back(Dot[Dot.size() - 1][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtS(Numbers[0][j], c) < lenghtS(Numbers[0][j], d))
            {
                dots2[0].push_back(Numbers[0][j]);
            }
            else
            {
                dots2[3].push_back(Numbers[0][j]);
            }
        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
            {
                dots2[1].push_back(Numbers[1][j]);
            }
            else
            {
                dots2[0].push_back(Numbers[1][j]);
            }
        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtS(Numbers[2][j], d) > lenghtS(Numbers[2][j], c))
            {
                dots2[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
            {
                dots2[3].push_back(Numbers[3][j]);
            }
        }

        Numbers.clear();

        dots2 = Sort_lenghts(dots2, lenghtS);
        dots2.erase(dots2.begin());
        dots2[2] = { d };
        S += S_nap(dots2);
    }
    else
    {
        dots2 = { {max_elm(Dot[Dot.size() - 1], 1)}, {max_elm(Dot[Dot.size() - 1], 0)} , {min_elm(Dot[Dot.size() - 1], 1)}, {min_elm(Dot[Dot.size() - 1], 0)} };

        S += area_of_triangle_Ger(dots2[1][0], dots2[0][0], dots2[3][0]);
    }
    S += area_of_triangle_Ger(dots1[0][0], dots1[2][0], dots2[2][0]);
    S += area_of_triangle_Ger(dots2[0][0], dots2[2][0], dots1[0][0]);

    return S;
}

template<typename T>
T Cave_S2_3(vector < vector <vector<T>>> Dot)
{
    T S = 0;
    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector < vector < vector <T> > >  Numbers;
    vector<T> a, b, c, d;
    vector<T> V = { 0,0 };
    for (int i = 0; i < Dot.size() - 2; i++)
    {
        a = max_elm(Dot[i], 0);
        b = min_elm(Dot[i], 0);
        c = max_elm(Dot[i], 1);
        d = min_elm(Dot[i], 1);
        dots1 = { {a},{c}, {b}, {d} };
        if (Dot[i].size() > 4)
        {
            Numbers = { {},{},{},{} };

            dots1 = { {a}, {c}, {b}, {d} };
            for (int j = 0; j < Dot[i].size(); j++) //Переписать
            {
                if (Dot[i][j] != a && Dot[i][j] != b && Dot[i][j] != c && Dot[i][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i][j], lenghtS)].push_back(Dot[i][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) > lenghtS(Numbers[0][j], d))
                {
                    dots1[3].push_back(Numbers[0][j]);
                }

            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots1[1].push_back(Numbers[1][j]);
                }

            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
                {
                    dots1[2].push_back(Numbers[2][j]);
                }
                else
                {
                    dots1[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots1[3].push_back(Numbers[3][j]);
                }
                else
                {
                    dots1[2].push_back(Numbers[3][j]);
                }
            }
            Numbers.clear();

            dots1 = Sort_lenghts(dots1, lenghtS);
            dots1.erase(dots1.begin());
            dots1[2] = { d };
            V[0] = (dots1[0][0][0] + dots1[1][0][0]) * 0.5; V[1] = (dots1[0][0][1] + dots1[1][0][1]) * 0.5;
            S += S_nap2(dots1, V);

        }
        else
        {
            dots1 = { {max_elm(Dot[i], 0)}, {max_elm(Dot[i], 1)} , {min_elm(Dot[i], 0)}, {min_elm(Dot[i], 1)} };
            S += area_of_triangle_Ger(dots1[1][0], dots1[2][0], dots1[3][0]);
            dots1.erase(dots1.begin());
        }
        if (Dot[i + 1].size() > 4)
        {
            a = max_elm(Dot[i + 1], 0);
            b = min_elm(Dot[i + 1], 0);
            c = max_elm(Dot[i + 1], 1);
            d = min_elm(Dot[i + 1], 1);
            dots2 = { {a}, {c}, {b}, {d} };

            Numbers = { {},{},{},{} };


            for (int j = 0; j < Dot[i + 1].size(); j++) //Переписать
            {
                if (Dot[i + 1][j] != a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i + 1][j], lenghtS)].push_back(Dot[i + 1][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) < lenghtS(Numbers[0][j], d))
                {
                    dots2[0].push_back(Numbers[0][j]);
                }
                else
                {
                    dots2[3].push_back(Numbers[0][j]);
                }
            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots2[1].push_back(Numbers[1][j]);
                }
                else
                {
                    dots2[0].push_back(Numbers[1][j]);
                }
            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) > lenghtS(Numbers[2][j], c))
                {
                    dots2[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots2[3].push_back(Numbers[3][j]);
                }
            }

            Numbers.clear();

            dots2 = Sort_lenghts(dots2, lenghtS);
            dots2.erase(dots2.begin());
            dots2[2] = { d };
            V[0] = (dots2[0][0][0] + dots2[1][0][0] ) * 0.5; V[1] = (dots2[0][0][1] + dots2[1][0][1]) * 0.5;
            S -= S_nap2(dots2,V);
        }
        else
        {
            dots2 = { {max_elm(Dot[i + 1], 1)}, {max_elm(Dot[i + 1], 0)} , {min_elm(Dot[i + 1], 1)}, {min_elm(Dot[i + 1], 0)} };
            S -= area_of_triangle_Ger(dots2[1][0], dots2[2][0], dots2[3][0]);
            dots2.erase(dots2.begin());
        }

        S += area_of_triangle_Ger(dots1[0][0], dots1[2][0], dots2[2][0]);
        S += area_of_triangle_Ger(dots2[0][0], dots2[2][0], dots1[0][0]);


    }

    a = max_elm(Dot[Dot.size() - 2], 0);
    b = min_elm(Dot[Dot.size() - 2], 0);
    c = max_elm(Dot[Dot.size() - 2], 1);
    d = min_elm(Dot[Dot.size() - 2], 1);
    dots1 = { {a},{c}, {b}, {d} };
    if (Dot[Dot.size() - 2].size() > 4)
    {
        Numbers = { {},{},{},{} };

        dots1 = { {a}, {c}, {b}, {d} };
        for (int j = 0; j < Dot[Dot.size() - 2].size(); j++) //Переписать
        {
            if (Dot[Dot.size() - 2][j] != a && Dot[Dot.size() - 2][j] != b && Dot[Dot.size() - 2][j] != c && Dot[Dot.size() - 2][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 2][j], lenghtS)].push_back(Dot[Dot.size() - 2][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtS(Numbers[0][j], c) > lenghtS(Numbers[0][j], d))
            {
                dots1[3].push_back(Numbers[0][j]);
            }

        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
            {
                dots1[1].push_back(Numbers[1][j]);
            }

        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
            {
                dots1[2].push_back(Numbers[2][j]);
            }
            else
            {
                dots1[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
            {
                dots1[3].push_back(Numbers[3][j]);
            }
            else
            {
                dots1[2].push_back(Numbers[3][j]);
            }
        }
        Numbers.clear();

        dots1 = Sort_lenghts(dots1, lenghtS);
        dots1.erase(dots1.begin());
        dots1[2] = { d };
        V[0] = (dots1[0][0][0] + dots1[1][0][0] ) * 0.5; V[1] = (dots1[0][0][1] + dots1[1][0][1] ) * 0.5;
        S += S_nap2(dots1,V);

    }
    else
    {
        dots1 = { {max_elm(Dot[Dot.size() - 2], 0)}, {max_elm(Dot[Dot.size() - 2], 1)} , {min_elm(Dot[Dot.size() - 2], 0)}, {min_elm(Dot[Dot.size() - 2], 1)} };
        S += area_of_triangle_Ger(dots1[1][0], dots1[2][0], dots1[3][0]);
        dots1.erase(dots1.begin());
    }
    if (Dot[Dot.size() - 1].size() > 4)
    {
        a = max_elm(Dot[Dot.size() - 1], 0);
        b = min_elm(Dot[Dot.size() - 1], 0);
        c = max_elm(Dot[Dot.size() - 1], 1);
        d = min_elm(Dot[Dot.size() - 1], 1);
        dots2 = { {a}, {c}, {b}, {d} };

        Numbers = { {},{},{},{} };


        for (int j = 0; j < Dot[Dot.size() - 1].size(); j++) //Переписать
        {
            if (Dot[Dot.size() - 1][j] != a && Dot[Dot.size() - 1][j] != b && Dot[Dot.size() - 1][j] != c && Dot[Dot.size() - 1][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 1][j], lenghtS)].push_back(Dot[Dot.size() - 1][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtS(Numbers[0][j], c) < lenghtS(Numbers[0][j], d))
            {
                dots2[0].push_back(Numbers[0][j]);
            }
            else
            {
                dots2[3].push_back(Numbers[0][j]);
            }
        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
            {
                dots2[1].push_back(Numbers[1][j]);
            }
            else
            {
                dots2[0].push_back(Numbers[1][j]);
            }
        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtS(Numbers[2][j], d) > lenghtS(Numbers[2][j], c))
            {
                dots2[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
            {
                dots2[3].push_back(Numbers[3][j]);
            }
        }

        Numbers.clear();

        dots2 = Sort_lenghts(dots2, lenghtS);
        dots2.erase(dots2.begin());
        dots2[2] = { d };
        V[0] = (dots2[0][0][0] + dots2[1][0][0] ) * 0.5; V[1] = (dots2[0][0][1] + dots2[1][0][1]) * 0.5;
        S += S_nap2(dots2,V);
    }
    else
    {
        dots2 = { {max_elm(Dot[Dot.size() - 1], 1)}, {max_elm(Dot[Dot.size() - 1], 0)} , {min_elm(Dot[Dot.size() - 1], 1)}, {min_elm(Dot[Dot.size() - 1], 0)} };
        S += area_of_triangle_Ger(dots2[1][0], dots2[0][0], dots2[3][0]);
        dots2.erase(dots2.begin());
    }
    S += area_of_triangle_Ger(dots1[0][0], dots1[2][0], dots2[2][0]);
    S += area_of_triangle_Ger(dots2[0][0], dots2[2][0], dots1[0][0]);

    return S;
}




template<typename T>
T Cave_V2(vector < vector <vector<T>>> Dot)
{
    T V = 0;
    T S1 = 0;
    T S2 = 0;
    T H = 0;
    vector<T> a, b, c, d;

    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector < vector < vector <T> > > Numbers;
    vector<T> V1{ 0, 0, 0 }, V2{ 0, 0, 0 };
 

    for (int i = 0; i < Dot.size() - 2; i++)
    {

            if (Dot[i].size() > 4)
            {
                Numbers = { {},{},{},{} };
                a = max_elm(Dot[i], 2);
                b = min_elm(Dot[i], 2);
                c = max_elm(Dot[i], 1);
                d = min_elm(Dot[i], 1);
                dots1 = { {a}, {c}, {b}, {d} };
                dots1 = dots_original_V(dots1, Dot[i]);
                a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
                for (int j = 0; j < Dot[i].size(); j++) 
                {
                    if (Dot[i][j] != a && Dot[i][j] != b && Dot[i][j] != c && Dot[i][j] != d)
                        Numbers[Sort_one(a, b, c, d, Dot[i][j], lenghtV)].push_back(Dot[i][j]);

                }

                for (int j = 0; j < Numbers[0].size(); j++)
                {
                    if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
                    {
                        dots1[0].push_back(Numbers[0][j]);
                    }
                    else
                    {
                        dots1[3].push_back(Numbers[0][j]);
                    }
                }

                for (int j = 0; j < Numbers[1].size(); j++)
                {
                    if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
                    {
                        dots1[1].push_back(Numbers[1][j]);
                    }
                    else
                    {
                        dots1[0].push_back(Numbers[1][j]);
                    }
                }

                for (int j = 0; j < Numbers[2].size(); j++)
                {
                    if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
                    {
                        dots1[2].push_back(Numbers[2][j]);
                    }
                    else
                    {
                        dots1[1].push_back(Numbers[2][j]);
                    }
                }

                for (int j = 0; j < Numbers[3].size(); j++)
                {
                    if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
                    {
                        dots1[3].push_back(Numbers[3][j]);
                    }
                    else
                    {
                        dots1[2].push_back(Numbers[3][j]);
                    }
                }

            }
            else
            {
                dots1 = { {max_elm(Dot[i], 2)}, {max_elm(Dot[i], 1)}, {min_elm(Dot[i], 2)}, {min_elm(Dot[i], 1)} };
                dots1 = dots_original_V(dots1, Dot[i]);
            }


            if (Dot[i + 1].size() > 4)
            {
                Numbers = { {},{},{},{} };
                a = max_elm(Dot[i + 1], 2);
                b = min_elm(Dot[i + 1], 2);
                c = max_elm(Dot[i + 1], 1);
                d = min_elm(Dot[i + 1], 1);
                dots2 = { {a}, {c}, {b}, {d} };
                dots2 = dots_original_V(dots2, Dot[i + 1]);
                a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
                for (int j = 0; j < Dot[i + 1].size(); j++)
                {
                    if (Dot[i + 1][j] != a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                        Numbers[Sort_one(a, b, c, d, Dot[i + 1][j], lenghtV)].push_back(Dot[i + 1][j]);

                }

                for (int j = 0; j < Numbers[0].size(); j++)
                {
                    if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
                    {
                        dots2[0].push_back(Numbers[0][j]);
                    }
                    else
                    {
                        dots2[3].push_back(Numbers[0][j]);
                    }
                }

                for (int j = 0; j < Numbers[1].size(); j++)
                {
                    if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
                    {
                        dots2[1].push_back(Numbers[1][j]);
                    }
                    else
                    {
                        dots2[0].push_back(Numbers[1][j]);
                    }
                }

                for (int j = 0; j < Numbers[2].size(); j++)
                {
                    if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
                    {
                        dots2[2].push_back(Numbers[2][j]);
                    }
                    else
                    {
                        dots2[1].push_back(Numbers[2][j]);
                    }
                }

                for (int j = 0; j < Numbers[3].size(); j++)
                {
                    if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
                    {
                        dots2[3].push_back(Numbers[3][j]);
                    }
                    else
                    {
                        dots2[2].push_back(Numbers[3][j]);
                    }
                }
                Numbers.clear();

            }
            else
            {
                dots2 = { {max_elm(Dot[i + 1], 2)}, {max_elm(Dot[i + 1], 1)} , {min_elm(Dot[i + 1], 2)}, {min_elm(Dot[i + 1], 1)} };
                dots2 = dots_original_V(dots2, Dot[i + 1]);
            }

            dots1 = Sort_lenghts(dots1, lenghtV);
            dots2 = Sort_lenghts(dots2, lenghtV);
            S1 = S_nap_V(dots1);
            S2 = S_nap_V(dots2);
            H = lenghtS(min_elm(Dot[i], 0), min_elm(Dot[i+1], 0));

            V +=  H*  (S1+ sqrt(S1*S2)+S2)/3;
        

            dots1.clear();
            dots2.clear();

    }

    if (Dot[Dot.size() - 2].size() > 4)
    {
        Numbers = { {},{},{},{} };
        a = max_elm(Dot[Dot.size() - 2], 2);
        b = min_elm(Dot[Dot.size() - 2], 2);
        c = max_elm(Dot[Dot.size() - 2], 1);
        d = min_elm(Dot[Dot.size() - 2], 1);
        dots1 = { {a}, {c}, {b}, {d} };
        dots1 = dots_original_S(dots1, Dot[Dot.size() - 2]);
        a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
        for (int j = 0; j < Dot[Dot.size() - 2].size(); j++) 
        {
            if (Dot[Dot.size() - 2][j] != a && Dot[Dot.size() - 2][j] != b && Dot[Dot.size() - 2][j] != c && Dot[Dot.size() - 2][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 2][j], lenghtV)].push_back(Dot[Dot.size() - 2][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
            {
                dots1[0].push_back(Numbers[0][j]);
            }
            else
            {
                dots1[3].push_back(Numbers[0][j]);
            }
        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
            {
                dots1[1].push_back(Numbers[1][j]);
            }
            else
            {
                dots1[0].push_back(Numbers[1][j]);
            }
        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
            {
                dots1[2].push_back(Numbers[2][j]);
            }
            else
            {
                dots1[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
            {
                dots1[3].push_back(Numbers[3][j]);
            }
            else
            {
                dots1[2].push_back(Numbers[3][j]);
            }
        }

    }
    else
    {
        dots1 = { {max_elm(Dot[Dot.size() - 2], 2)}, {max_elm(Dot[Dot.size() - 2], 1)}, {min_elm(Dot[Dot.size() - 2], 2)}, {min_elm(Dot[Dot.size() - 2], 1)} };
        dots1 = dots_original_V(dots1, Dot[Dot.size() - 2]);
    }

    if (Dot[Dot.size() - 1].size() > 4)
    {
        Numbers = { {},{},{},{} };
        a = max_elm(Dot[Dot.size() - 1], 2);
        b = min_elm(Dot[Dot.size() - 1], 2);
        c = max_elm(Dot[Dot.size() - 1], 1);
        d = min_elm(Dot[Dot.size() - 1], 1);
        dots2 = { {a}, {c}, {b}, {d} };
        dots2 = dots_original_V(dots2, Dot[Dot.size() - 1]);
        a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
        for (int j = 0; j < Dot[Dot.size() - 1].size(); j++) 
        {
            if (Dot[Dot.size() - 1][j] != a && Dot[Dot.size() - 1][j] != b && Dot[Dot.size() - 1][j] != c && Dot[Dot.size() - 1][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 1][j], lenghtV)].push_back(Dot[Dot.size() - 1][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
            {
                dots2[0].push_back(Numbers[0][j]);
            }
            else
            {
                dots2[3].push_back(Numbers[0][j]);
            }
        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
            {
                dots2[1].push_back(Numbers[1][j]);
            }
            else
            {
                dots2[0].push_back(Numbers[1][j]);
            }
        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
            {
                dots2[2].push_back(Numbers[2][j]);
            }
            else
            {
                dots2[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
            {
                dots2[3].push_back(Numbers[3][j]);
            }
            else
            {
                dots2[2].push_back(Numbers[3][j]);
            }
        }
        Numbers.clear();

    }
    else
    {
        dots2 = { {max_elm(Dot[Dot.size() - 1], 2)}, {max_elm(Dot[Dot.size() - 1], 1)} , {min_elm(Dot[Dot.size() - 1], 2)}, {min_elm(Dot[Dot.size() - 1], 1)} };
        dots2 = dots_original_S(dots2, Dot[Dot.size() - 1]);
    }
    dots1 = Sort_lenghts(dots1, lenghtV);
    dots2 = Sort_lenghts(dots2, lenghtV);
    S1 = S_nap_V(dots1);
    S2 = S_nap_V(dots2);
    H = lenghtS(min_elm(Dot[Dot.size() - 2], 0), max_elm(Dot[Dot.size() - 1], 0));

    V += H * (S1 + sqrt(S1 * S2) + S2) / 3.0;
    return V;
}


template<typename T>
T Cave_V2_2(vector < vector <vector<T>>> Dot)
{
    T V = 0;
    T S1 = 0;
    T S2 = 0;
    T H = 0;
    vector<T> a, b, c, d;

    vector < vector < vector <T> > > dots1;
    vector < vector < vector <T> > > dots2;
    vector < vector < vector <T> > > Numbers;
    vector<T> V1{ 0,0, 0 }; vector<T> V2{ 0,0, 0 };


    for (int i = 0; i < Dot.size() - 2; i++)
    {

        if (Dot[i].size() > 4)
        {
            Numbers = { {},{},{},{} };
            a = max_elm(Dot[i], 2);
            b = min_elm(Dot[i], 2);
            c = max_elm(Dot[i], 1);
            d = min_elm(Dot[i], 1);
            dots1 = { {a}, {c}, {b}, {d} };
            dots1 = dots_original_V(dots1, Dot[i]);
            a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
            for (int j = 0; j < Dot[i].size(); j++)
            {
                if (Dot[i][j] != a && Dot[i][j] != b && Dot[i][j] != c && Dot[i][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i][j], lenghtV)].push_back(Dot[i][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
                {
                    dots1[0].push_back(Numbers[0][j]);
                }
                else
                {
                    dots1[3].push_back(Numbers[0][j]);
                }
            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
                {
                    dots1[1].push_back(Numbers[1][j]);
                }
                else
                {
                    dots1[0].push_back(Numbers[1][j]);
                }
            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
                {
                    dots1[2].push_back(Numbers[2][j]);
                }
                else
                {
                    dots1[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
                {
                    dots1[3].push_back(Numbers[3][j]);
                }
                else
                {
                    dots1[2].push_back(Numbers[3][j]);
                }
            }

        }
        else
        {
            dots1 = { {max_elm(Dot[i], 2)}, {max_elm(Dot[i], 1)}, {min_elm(Dot[i], 2)}, {min_elm(Dot[i], 1)} };
            dots1 = dots_original_V(dots1, Dot[i]);
        }


        if (Dot[i + 1].size() > 4)
        {
            Numbers = { {},{},{},{} };
            a = max_elm(Dot[i + 1], 2);
            b = min_elm(Dot[i + 1], 2);
            c = max_elm(Dot[i + 1], 1);
            d = min_elm(Dot[i + 1], 1);
            dots2 = { {a}, {c}, {b}, {d} };
            dots2 = dots_original_V(dots2, Dot[i + 1]);
            a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
            for (int j = 0; j < Dot[i + 1].size(); j++)
            {
                if (Dot[i + 1][j] != a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i + 1][j], lenghtV)].push_back(Dot[i + 1][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
                {
                    dots2[0].push_back(Numbers[0][j]);
                }
                else
                {
                    dots2[3].push_back(Numbers[0][j]);
                }
            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
                {
                    dots2[1].push_back(Numbers[1][j]);
                }
                else
                {
                    dots2[0].push_back(Numbers[1][j]);
                }
            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
                {
                    dots2[2].push_back(Numbers[2][j]);
                }
                else
                {
                    dots2[1].push_back(Numbers[2][j]);
                }
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
                {
                    dots2[3].push_back(Numbers[3][j]);
                }
                else
                {
                    dots2[2].push_back(Numbers[3][j]);
                }
            }
            Numbers.clear();

        }
        else
        {
            dots2 = { {max_elm(Dot[i + 1], 2)}, {max_elm(Dot[i + 1], 1)} , {min_elm(Dot[i + 1], 2)}, {min_elm(Dot[i + 1], 1)} };
            dots2 = dots_original_V(dots2, Dot[i + 1]);
        }


        dots1 = Sort_lenghts(dots1, lenghtV);
        dots2 = Sort_lenghts(dots2, lenghtV);
        V1[2] = (dots1[1][0][2] + dots1[0][0][2] + dots1[2][0][2] + dots1[3][0][2]) * 0.25; V1[1] = (dots1[1][0][1] + dots1[0][0][1] + dots1[2][0][1] + dots1[3][0][1]) * 0.25;
        V2[2] = (dots1[1][0][2] + dots1[0][0][2] + dots1[2][0][2] + dots1[3][0][2]) * 0.25; V2[1] = (dots1[1][0][1] + dots1[0][0][1] + dots1[2][0][1] + dots1[3][0][1]) * 0.25;
        S1 = S_nap2_V(dots1, V1);
        S2 = S_nap2_V(dots2, V2);
        H = lenghtS(min_elm(Dot[i], 0), min_elm(Dot[i + 1], 0));

        V += H * (S1 + sqrt(S1 * S2) + S2) / 3;


        dots1.clear();
        dots2.clear();

    }

    if (Dot[Dot.size() - 2].size() > 4)
    {
        Numbers = { {},{},{},{} };
        a = max_elm(Dot[Dot.size() - 2], 2);
        b = min_elm(Dot[Dot.size() - 2], 2);
        c = max_elm(Dot[Dot.size() - 2], 1);
        d = min_elm(Dot[Dot.size() - 2], 1);
        dots1 = { {a}, {c}, {b}, {d} };
        dots1 = dots_original_V(dots1, Dot[Dot.size() - 2]);
        a = dots1[0][0]; c = dots1[1][0]; b = dots1[2][0]; d = dots1[3][0];
        for (int j = 0; j < Dot[Dot.size() - 2].size(); j++)
        {
            if (Dot[Dot.size() - 2][j] != a && Dot[Dot.size() - 2][j] != b && Dot[Dot.size() - 2][j] != c && Dot[Dot.size() - 2][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 2][j], lenghtV)].push_back(Dot[Dot.size() - 2][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
            {
                dots1[0].push_back(Numbers[0][j]);
            }
            else
            {
                dots1[3].push_back(Numbers[0][j]);
            }
        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
            {
                dots1[1].push_back(Numbers[1][j]);
            }
            else
            {
                dots1[0].push_back(Numbers[1][j]);
            }
        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
            {
                dots1[2].push_back(Numbers[2][j]);
            }
            else
            {
                dots1[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
            {
                dots1[3].push_back(Numbers[3][j]);
            }
            else
            {
                dots1[2].push_back(Numbers[3][j]);
            }
        }

    }
    else
    {
        dots1 = { {max_elm(Dot[Dot.size() - 2], 2)}, {max_elm(Dot[Dot.size() - 2], 1)}, {min_elm(Dot[Dot.size() - 2], 2)}, {min_elm(Dot[Dot.size() - 2], 1)} };
        dots1 = dots_original_V(dots1, Dot[Dot.size() - 2]);
    }

    if (Dot[Dot.size() - 1].size() > 4)
    {
        Numbers = { {},{},{},{} };
        a = max_elm(Dot[Dot.size() - 1], 2);
        b = min_elm(Dot[Dot.size() - 1], 2);
        c = max_elm(Dot[Dot.size() - 1], 1);
        d = min_elm(Dot[Dot.size() - 1], 1);
        dots2 = { {a}, {c}, {b}, {d} };
        dots2 = dots_original_V(dots2, Dot[Dot.size() - 1]);
        a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
        for (int j = 0; j < Dot[Dot.size() - 1].size(); j++)
        {
            if (Dot[Dot.size() - 1][j] != a && Dot[Dot.size() - 1][j] != b && Dot[Dot.size() - 1][j] != c && Dot[Dot.size() - 1][j] != d)
                Numbers[Sort_one(a, b, c, d, Dot[Dot.size() - 1][j], lenghtV)].push_back(Dot[Dot.size() - 1][j]);

        }

        for (int j = 0; j < Numbers[0].size(); j++)
        {
            if (lenghtV(Numbers[0][j], c) < lenghtV(Numbers[0][j], d))
            {
                dots2[0].push_back(Numbers[0][j]);
            }
            else
            {
                dots2[3].push_back(Numbers[0][j]);
            }
        }

        for (int j = 0; j < Numbers[1].size(); j++)
        {
            if (lenghtV(Numbers[1][j], b) < lenghtV(Numbers[1][j], a))
            {
                dots2[1].push_back(Numbers[1][j]);
            }
            else
            {
                dots2[0].push_back(Numbers[1][j]);
            }
        }

        for (int j = 0; j < Numbers[2].size(); j++)
        {
            if (lenghtV(Numbers[2][j], d) < lenghtV(Numbers[2][j], c))
            {
                dots2[2].push_back(Numbers[2][j]);
            }
            else
            {
                dots2[1].push_back(Numbers[2][j]);
            }
        }

        for (int j = 0; j < Numbers[3].size(); j++)
        {
            if (lenghtV(Numbers[3][j], a) < lenghtV(Numbers[3][j], b))
            {
                dots2[3].push_back(Numbers[3][j]);
            }
            else
            {
                dots2[2].push_back(Numbers[3][j]);
            }
        }
        Numbers.clear();

    }
    else
    {
        dots2 = { {max_elm(Dot[Dot.size() - 1], 2)}, {max_elm(Dot[Dot.size() - 1], 1)} , {min_elm(Dot[Dot.size() - 1], 2)}, {min_elm(Dot[Dot.size() - 1], 1)} };
        dots2 = dots_original_V(dots2, Dot[Dot.size() - 1]);
    }
    dots1 = Sort_lenghts(dots1, lenghtV);
    dots2 = Sort_lenghts(dots2, lenghtV);
    V1[2] = (dots1[1][0][2] + dots1[0][0][2] + dots1[2][0][2] + dots1[3][0][2]) * 0.25; V1[1] = (dots1[1][0][1] + dots1[0][0][1] + dots1[2][0][1] + dots1[3][0][1]) * 0.25;
    V2[2] = (dots1[1][0][2] + dots1[0][0][2] + dots1[2][0][2] + dots1[3][0][2]) * 0.25; V2[1] = (dots1[1][0][1] + dots1[0][0][1] + dots1[2][0][1] + dots1[3][0][1]) * 0.25;
    S1 = S_nap2_V(dots1, V1);
    S2 = S_nap2_V(dots2, V2);
    H = lenghtS(min_elm(Dot[Dot.size() - 2], 0), min_elm(Dot[Dot.size() - 1], 0));

    V += H * (S1 + sqrt(S1 * S2) + S2) / 3;

    return V;
}

template<typename T>
T Cave_V3elN1(vector < vector <vector<T>>> Dot)
{
    T S = 0;
    T C = 0;
    T V = 0;
    T coef = 2.0 / 3;
    vector < vector < vector <T> > >  dots2;
    vector < vector < vector <T> > >  Numbers;
    vector< vector <vector<T>>> Temp;
    vector<T> a, b, c, d;
    for (int i = 0; i < Dot.size() - 2; i++)
    {
        Temp = { Dot[i], Dot[i + 1] };
        S = Cave_S2_2(Temp);
        C = (sqrt((max_elm(Dot[i], 2)[2] - min_elm(Dot[i], 2)[2]) * (max_elm(Dot[i], 2)[2] - min_elm(Dot[i], 2)[2])) +
            sqrt((max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]) * (max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]))) *
            0.5;
        V += coef * S * C;

        if (Dot[i + 1].size() > 4)
        {
            a = max_elm(Dot[i + 1], 0);
            b = min_elm(Dot[i + 1], 0);
            c = max_elm(Dot[i + 1], 1);
            d = min_elm(Dot[i + 1], 1);
            dots2 = { {a}, {c}, {b}, {d} };
            dots2 = dots_original_S(dots2, Dot[i + 1]);
            a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
            Numbers = { {},{},{},{} };


            for (int j = 0; j < Dot[i + 1].size(); j++) 
            {
                if (Dot[i + 1][j] != a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i + 1][j], lenghtS)].push_back(Dot[i + 1][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) < lenghtS(Numbers[0][j], d))
                {
                    dots2[0].push_back(Numbers[0][j]);
                }
                else
                {
                    dots2[3].push_back(Numbers[0][j]);
                }
            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots2[1].push_back(Numbers[1][j]);
                }
                else
                {
                    dots2[0].push_back(Numbers[1][j]);
                }
            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
                {
                    dots2[2].push_back(Numbers[2][j]);
                }
                else
                    dots2[1].push_back(Numbers[2][j]);
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots2[3].push_back(Numbers[3][j]);
                }
                else
                    dots2[2].push_back(Numbers[3][j]);
                
               
            }

            Numbers.clear();

            dots2 = Sort_lenghts(dots2, lenghtS);
            S = S_nap(dots2);
        }
        else
        {
            dots2 = { {max_elm(Dot[i + 1], 1)}, {max_elm(Dot[i + 1], 0)} , {min_elm(Dot[i + 1], 1)}, {min_elm(Dot[i + 1], 0)} };
            dots2 = dots_original_S(dots2, Dot[i + 1]);
            S = area_of_triangle_Ger(dots2[1][0], dots2[2][0], dots2[3][0]) + area_of_triangle_Ger(dots2[1][0], dots2[0][0], dots2[3][0]);

        }
        C = (sqrt((max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]) * (max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]))) *
            0.5;
        V -= coef * S * C;
    }
    Temp = { Dot[Dot.size() - 2],Dot[Dot.size() - 1] };
    S = Cave_S2_2(Temp);
    C = (sqrt((max_elm(Dot[Dot.size() - 2], 2)[2] - min_elm(Dot[Dot.size() - 2], 2)[2]) * (max_elm(Dot[Dot.size() - 2], 2)[2] - min_elm(Dot[Dot.size() - 2], 2)[2])) +
        sqrt((max_elm(Dot[Dot.size() - 1], 2)[2] - min_elm(Dot[Dot.size() - 1], 2)[2]) * (max_elm(Dot[Dot.size() - 1], 2)[2] - min_elm(Dot[Dot.size() - 1], 2)[2]))) *
        0.5;
    V += coef * S * C;
    return V;
}


template<typename T>
T Cave_V3elN2(vector < vector <vector<T>>> Dot)
{
    T S = 0;
    T C = 0;
    T V = 0;
    T coef = 2.0 / 3;
    vector < vector < vector <T> > >  dots2;
    vector < vector < vector <T> > >  Numbers;
    vector< vector <vector<double>>> Temp;
    vector<T> a, b, c, d;
    for (int i = 0; i < Dot.size() - 2; i++)
    {
        Temp = { Dot[i], Dot[i + 1] };
        S = Cave_S2(Temp);
        C = (sqrt((max_elm(Dot[i], 2)[2] - min_elm(Dot[i], 2)[2]) * (max_elm(Dot[i], 2)[2] - min_elm(Dot[i], 2)[2])) +
            sqrt((max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]) * (max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]))) *
            0.5;
        V += coef * S * C;

        if (Dot[i + 1].size() > 4)
        {
            a = max_elm(Dot[i + 1], 0);
            b = min_elm(Dot[i + 1], 0);
            c = max_elm(Dot[i + 1], 1);
            d = min_elm(Dot[i + 1], 1);
            dots2 = { {a}, {c}, {b}, {d} };
            dots2 = dots_original_V(dots2, Dot[i + 1]);
            a = dots2[0][0]; c = dots2[1][0]; b = dots2[2][0]; d = dots2[3][0];
            Numbers = { {},{},{},{} };


            for (int j = 0; j < Dot[i + 1].size(); j++)
            {
                if (Dot[i + 1][j] != a && Dot[i + 1][j] != b && Dot[i + 1][j] != c && Dot[i + 1][j] != d)
                    Numbers[Sort_one(a, b, c, d, Dot[i + 1][j], lenghtS)].push_back(Dot[i + 1][j]);

            }

            for (int j = 0; j < Numbers[0].size(); j++)
            {
                if (lenghtS(Numbers[0][j], c) < lenghtS(Numbers[0][j], d))
                {
                    dots2[0].push_back(Numbers[0][j]);
                }
                else
                {
                    dots2[3].push_back(Numbers[0][j]);
                }
            }

            for (int j = 0; j < Numbers[1].size(); j++)
            {
                if (lenghtS(Numbers[1][j], b) < lenghtS(Numbers[1][j], a))
                {
                    dots2[1].push_back(Numbers[1][j]);
                }
                else
                {
                    dots2[0].push_back(Numbers[1][j]);
                }
            }

            for (int j = 0; j < Numbers[2].size(); j++)
            {
                if (lenghtS(Numbers[2][j], d) < lenghtS(Numbers[2][j], c))
                {
                    dots2[2].push_back(Numbers[2][j]);
                }
                else
                    dots2[1].push_back(Numbers[2][j]);
            }

            for (int j = 0; j < Numbers[3].size(); j++)
            {
                if (lenghtS(Numbers[3][j], a) < lenghtS(Numbers[3][j], b))
                {
                    dots2[3].push_back(Numbers[3][j]);
                }
                else
                    dots2[2].push_back(Numbers[3][j]);


            }

            Numbers.clear();

            dots2 = Sort_lenghts(dots2, lenghtS);
            for (int i = 0; i < dots2[1].size() - 1; i++)
            {
                S += area_of_triangle_Ger(c, dots2[1][i + 1], dots2[1][i]);
            }
            S += area_of_triangle_Ger(c, b, d);
            for (int i = 0; i < dots2[2].size() - 1; i++)
            {
                S += area_of_triangle_Ger(d, dots2[2][i + 1], dots2[2][i]);
            }
        }
        else
        {
            dots2 = { {max_elm(Dot[i + 1], 1)}, {max_elm(Dot[i + 1], 0)} , {min_elm(Dot[i + 1], 1)}, {min_elm(Dot[i + 1], 0)} };
            dots2 = dots_original_V(dots2, Dot[i + 1]);
            S = area_of_triangle_Ger(dots2[1][0], dots2[2][0], dots2[3][0]) + area_of_triangle_Ger(dots2[1][0], dots2[0][0], dots2[3][0]);

        }
        C = (sqrt((max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]) * (max_elm(Dot[i + 1], 2)[2] - min_elm(Dot[i + 1], 2)[2]))) *
            0.5;
        V -= coef * S * C;
    }
    Temp = { Dot[Dot.size() - 2],Dot[Dot.size() - 1] };
    S = Cave_S2(Temp);
    C = (sqrt((max_elm(Dot[Dot.size() - 2], 2)[2] - min_elm(Dot[Dot.size() - 2], 2)[2]) * (max_elm(Dot[Dot.size() - 2], 2)[2] - min_elm(Dot[Dot.size() - 2], 2)[2])) +
        sqrt((max_elm(Dot[Dot.size() - 1], 2)[2] - min_elm(Dot[Dot.size() - 1], 2)[2]) * (max_elm(Dot[Dot.size() - 1], 2)[2] - min_elm(Dot[Dot.size() - 1], 2)[2]))) *
        0.5;
    V += coef * S * C;
    return V;
}

template<typename T>
T Priv_V(vector < vector <vector<T>>> Dot)
{
    vector< vector <vector<T>>> Temp = { Dot[0], Dot[1] };
    T S = Cave_S(Temp);
    T h = (lenghtV(min_elm(Dot[0], 2), max_elm(Dot[0], 2))+ lenghtV(min_elm(Dot[1], 2), max_elm(Dot[1], 2)))*0.5;
    return S * h * 0.5;

}


int main()
{
    setlocale(LC_ALL, "Russian");
    vector< vector <vector<double>>> picket = {

{{23.39, 54.1, 86.2},
{22.33, 291.0, 56.2},
{10.42, 299.5, -0.4},
{11.24, 105.8, -4.1},
{26.71, 106.2, 54.3},
{18.68, 213.2, 41.7},
{36.44, 149.1, 50.8},
{39.86, 176.3, 25.0},
{26.59, 146.0, 29.0},
{17.87, 124.6, 49.7},
{22.72, 316.2, 36.6},
{17.86, 308.5, 7.7},
{20.92, 339.5,  11.8},
{27.01,  33.9, 36.5},
{22.51, 31.0 ,-4.7},
{11.43, 53.9 ,-17.5},
{18.22 ,111.7 , -27.5}},


{{23.99, 356.9, 11.9},
{1.24, 6.9, -79.7},
{24.16, 71.8, 82.2},
{18.63, 275.4, 0.1},
{9.54, 115.8, -1.2}},


{{19.99, 32.7, 14.3},
{0.55, 228.0, -83.5},
{2.11, 72.0, 85.9},
{2.51, 311.1, 0.0},
{2.39, 113.3, 15.2},
{8.73, 12.8 ,16.5},
{14.13, 29.7, 43.0},
{7.52, 59.2 ,36.3},
{11.82, 26.2, 13.1},
{2.21, 324.8 ,52.4},
{1.32, 115.6 ,-2.6},
{1.68, 120.2 ,-6.2},
{1.17 ,145.0, 85.5}},

{{14.06, 42.4, 6.9},
{0.61, 342.9, -84.8},
{0.22, 147.0, 82.9},
{0.16, 306.1, -2.4},
{1.18, 108.1, -2.4},
{1.16, 109.1, -43.9}},

{{0.30, 104.4, 35.2},
{7.71, 32.2, 29.7},
{1.10, 52.3, -84.8},
{5.05, 108.6, 87.0},
{3.34, 306.6, 2.1},
{3.02, 117.3, 6.1},
{5.89, 72.4, 44.1},
{4.70, 53.1, 20.8},
{4.99, 2.7, 29.9},
{5.47, 79.5, 72.4}},



{ {12.51,	18.7,	21.4	},
{0.16,	322.8,	-82.0},
{3.18,	94.2,	85.1},
{1.40,	314.6,	13.2},
{1.66	,119.8,	15.3},
{4.62,	63.1,	42.0},
{6.53,	353.1,	56.4},
{2.28,	327.8,	29.6},},


{ {14.83,	72.0,	13.7},
{0.35,	17.1,	-79.5},
{2.48,	15.9,	85.0},
{2.04,	338.1,	5.8},
{1.47,	155.9,	4.5},},


{{12.18, 	353.6,	 31.7},
{1.20,	356.5,	-78.1},
{3.03,	213.9,	84.8},
{8.07,	300.4,	7.0},
{2.33,	113.8,	20.2},
{6.31,	289.2,	50.0},
{2.39,	189.1,	11.5},
{6.52,	282.5,	61.5},},


{{24.51, 	25.4,	30.0},
{0.44,	112.9,	-84.0},
{3.80,	116.4,	85.1},
{3.80,	295.3,	14.2},
{1.70,	96.6,	0.4},
{3.80,	326.6,	73.8},},


{{16.38,	 50.7,	23.3},
{0.14,	25.9,	-80.9},
{2.54,	186.8,	78.4},
{2.90,	327.4,	8.3},
{2.36,	144.9,	5.2},
{9.98,	31.9,	4.2},},


{{10.02,	 4.4,	-6.9},
{0.13,	339.4,	-81.3},
{3.12,	115.3,	84.4},
{2.45,	290.1,	6.1},
{1.72,	107.0,	6.2},},


{{6.68,	32.5,	29.1},
{0.48,	306.5,	-82.4},
{2.33,	75.7,	81.3},
{2.04,	303.4,	14.2},
{2.03,	120.8,	2.0},},


{{5.07,	236.9,	24.3},
{1.28,	230.3,	-81.5},
{3.03,	251.7,	85.1},
{3.84,	211.4,	10.0},
{0.46,	19.9,	-2.7},},


{{3.50,	307.4,	24.3},
{0.14,	73.1,	-78.3},
{2.86,	176.5,	84.5},
{1.53,	165.5,	2.1},
{1.49,	348.0,	8.8},
{6.66,	332.9,	24.6},},


{{7.15,	333.0,	24.6},
{0.94,	53.5,	-84.9},
{2.39,	203.9,	89.3},
{4.24,	255.3,	-1.4},
{1.03,	75.0,	6.7},
{5.06,	308.1,	13.3},
{9.27,	274.5,	-4.0},
{6.72,	26.3,	45.5},},

{{7.43,	297.0,	1.5},
{1.55,	95.3,	-86.3},
{1.65,	56.1,	82.6},
{2.72,	229.4,	-0.5},
{0.66,	39.0,	12.2},},


{{6.59,	349.8,	-25.9},
{0.49,	126.0,	-85.1},
{5.05,	13.1,	82.1},
{2.82,	251.3,	9.2},
{1.44,	73.9,	11.7},
{11.71,	281.1,	68.3},
{12.42,	256.1,	54.0},
{9.54,	205.0,	78.1},
{4.79,	228.7,	23.7},
{5.45,	181.9,	82.0},
{2.64,	60.7,	53.3},
{3.71,	337.0,	11.5},
{3.54,	26.2,	14.3},},


{{7.00,	10.0,	-30.6},
{1.51,	250.9,	-79.9},
{1.08,	70.7,	82.9},
{0.14,	287.6,	8.8},
{0.76,	96.3,	-0.4},
{2.99,	9.4,	5.1},},


{{17.89,	 9.5,	5.1},
{0.75,	269.1,	-83.8},
{0.99,	83.5,	82.8},
{1.58,	290.4,	3.6},
{0.64,	97.0,	2.3},
{3.08,	351.5,	1.9},},











    };


    vector< vector <vector<double>>> picket2 = {




    { {9.66	,33.9,	40.5},
    {14.19	,55.1,	85.6},
    {11.13	,123.3,	9.3},
    {1.62	,81.9,	-85.4},
    {0.57	,277.4,	-0.2},
    {5.14	,354.6,	77.1},
    {14.19	,126.0,	84.4},
    {24.62	,109.2,	51.0} },



    { {7.23,	123.3,	-17.6},
    {0.97	,125.1,	-85.0},
    {0.94,	10.2,	86.1},
    {0.74,	262.3,	9.8},
    {6.76,	91.4,	5.7},
    {10.90,	92.6,	10.4},
    {1.04	,187.3,	-82.5},
    {0.99	,83.4,	86.4},
    {19.11	,98.2,	58.1},
    {0.94	,262.4,	-2.5},
    {10.63	,85.4,	10.7} },

    { {8.20	,27.6,	25.5},
    {1.04	,18.0,	-80.0},
    {16.94	,129.5,	83.5},
    {1.84	,258.9,	0.5},
    {9.31	,80.6,	8.9} },

    { {11.80,	354.8,	16.0},
    {0.49,	345.9,	-80.0},
    {4.58,	11.8,	87.1},
    {1.16,	275.4,	-0.6},
    {2.88,	90.7,	12.2} },


    { {13.20,	354.8,	-7.1},
    {1.58,	282.3,	-74.4},
    {9.02,	341.6,	84.7},
    {0.60,	257.5,	-4.2},
    {0.82,	79.2,	2.2} },

    { {59.2,	352.9,	16.6},
    {1.96,	245.4,	-73.7},
    {36.46,	46.0,	70.7},
    {1.42	,260.9,	3.9},
    {2.34	,76.0,	2.4} },


    { {2.77	,354.6,	-5.6},
    {0.76	,240.2,	-79.2},
    {3.96	,312.6,	86.1},
    {1.53,	241.9,	0.4},
    {5.02,	64.0,	15.2} },

    {{ 13.02,	321.2,	-28.9 },
    { 0.42,	215.6,	-85.1 },
    { 5.49,	335.5,	83.5 },
    { 2.83,	246.6,	-2.8 },
    { 4.74,	53.0,	13.3 },
    { 5.41,	180.0,	6.9 },
    { 14.67,	356.0,	7.3 },
    { 8.15	,11.6,	35.7 },
    { 6.32	,336.9,	24.9 },
    { 5.96	,320.1,	1.9 }},

    {{15.14,	327.2,	-9.2},
    {0.71,	229.9,	-83.5},
    {2.65,	340.4,	75.7},
    {4.30,	204.2,	-2.4},
    {10.14,	7.9	,4.4},
    {8.12,	64.0,	21.8},
    {7.45,	67.8,	2.0},
    {7.25,	262.0,	-6.6},
    {4.85,	272.7,	37.9} },

    {{20.18,	 303.6,	1.6},
    {0.14,	168.7,	-80.4},
    {1.94,	313.7,	86.2},
    {3.74,	242.2,	8.3},
    {5.05,	56.7,	-5.0},
    {4.48,	8.4,	24.5},
    {4.93,	341.6,	42.2},
    {5.16,	275.6,	38.8},
    {13.43,	279.9,	17.1} },

    /*{{14.30,	301.4,	15.9},
    {0.61,	330.9,	-81.9}	,
    {1.82,	154.8,	79.6}	,
    {3.68,	277.9,	7.6},
    {10.22,	88.1,	0.1},
    {9.71,	357.8,	9.2},
    {9.98,	333.8,	13.2},
    {28.93,	122.6,	-13.8} },

    { {10.39,	35.6,	3.3},
    {1.32,	280.1,	-76.5},
    {1.28,	307.7,	79.3},
    {0.96,	275.7,	-5.5},
    {0.67,	84.5,	4.2},
    {4.60,	323.1,	46.0},
    {2.39,	5.2,	12.7} },


    { {5.50,	9.3,	-38.9},
    {1.19,	217.9,	-80.2},
    {1.60,	7.2,	83.8},
    {3.56,	318.4,	17.4},
    {2.48,	113.0,	-7.0},
    {1.00,	39.0,	14.0},
    {5.11,	336.9,	29.3},
    {3.96,	240.9,	63.3} },

    { {7.95,	29.3,	-23.7},
    {0.33,	276.9,	-80.3},
    {2.36,	308.3,	88.1},
    {1.82,	282.4,	7.0},
    {2.51,	75.3,	9.5} },

    { {17.16, 	354.3,	6.3},
    {0.31,	219.2,	-68.3},
    {2.24,	323.1,	88.6},
    {3.51,	235.1,	11.1},
    {3.52,	50.5,	0.1},
    {4.57,	73.5,	-6.6},
    {11.93,	292.5,	34.7},
    {7.11, 314.2,	29.9} },

    { {4.33,	333.3,	21.7},
    {0.79,	259.0,	-84.3},
    {0.79,	253.6,	11.5},
    {2.84,	354.2,	81.8},
    {2.79,	62.1,	5.1},
    {3.94,	292.1,	23.5},
    {10.17,	308.5,	-1.0} },


    { {5.25,	354.0,	-8.4},
    {0.89,	194.7,	-87.4},
    {2.23,	5.8,	87.0},
    {6.59,	271.3,	0.9},
    {2.53,	78.0,	-11.3},
    {5.07,	29.1,	-18.1},
    {9.52,	325.4,	-14.7} },

    { {17.05, 	17.1,	-17.7},
    {0.15,	249.7,	-79.3},
    {9.03,	35.9,	84.1},
    {0.63,	275.8,	2.4},
    {3.17,	108.6,	9.9},
    {5.42,	325.8,	52.3} },


    { {4.55,	7.0,	-4.8},
    {1.27,	290.6,	-87.4},
    {0.16,	357.0,	75.4},
    {0.36,	309.7,	-12.2},
    {1.55,	120.5,	-0.3},
    {8.55,	110.6,	47.4},
    {6.29,	72.1,	19.4},
    {5.04,	33.2,	-12.5} },



    { {11.87, 	42.2,	-20.4},
    {1.98,	343.1,	-72.2},
    {6.04,	56.3,	77.2},
    {3.75,	340.9,	-7.7},
    {0.75,	142.2,	-3.4},
    {7.81,	60.2,	-8.6} },


    { {12.03, 	64.9,	-48.7},
    {1.24,	323.0,	-83.5},
    {6.53,	60.4,	78.7},
    {3.53,	290.5,	0.1},
    {0.16,	101.7,	-1.0} },


    { {10.59, 	336.7,	-3.7},
    {1.29,	184.4,	-79.1},
    {5.32,	39.2,	74.8},
    {0.15,	255.2,	-1.1},
    {4.47,	74.9,	2.2} },

    { {27.65, 	11.8,	-22.9},
    {0.55,	227.8,	-83.3},
    {10.56,	21.5,	73.9},
    {1.27,	265.4,	7.4},
    {3.49,	77.3,	-4.0} },


    { {26.95,  	1.4,	-11.9},
    {0.15,	265.1,	-81.3},
    {6.28,	123.6,	80.7},
    {3.22,	260.6,	5.5},
    {4.38,	53.3,	1.4},
    {4.85,	117.1,	10.0},
    {6.11,	90.2,	20.6},
    {6.00,	71.9,	4.7},
    {4.49,	7.8,	42.4},
    {4.52,	309.6,	11.5} },



    { {12.56, 	346.5,	3.9},
    {0.64,	306.6,	-68.1},
    {2.83,	12.9,	82.8},
    {2.46,	245.9,	9.2},
    {4.89,	69.4,	0.0},
    {10.46,	10.1,	-1.8},
    {7.61, 281.1,	9.3},
    {13.98,	281.8,	37.8},
    {16.13,	302.5,	19.1},
    {15.73,	322.2,	14.1},
    {13.62,	329.6,	27.6} },


    { {8.69,	330.9,	6.9},
    {0.94,	304.5,	-80.9},
    {0.88,	319.1,	81.4},
    {1.72,	254.8,	1.5},
    {0.17,	62.7,	0.7},
    {3.02,	354.7,	-66.4},
    {2.18,	123.8,	-55.2} },

    { {5.58,	349.8,	-29.9},
    {0.31,	291.1,	-71.8},
    {0.53,	128.8,	80.0},
    {0.56,	250.8,	16.7},
    {0.20,	82.8,	-7.7} },

    { {7.02,	346.2,	4.9},
    {1.32,	245.9,	-75.1},
    {6.53,	1.0,	80.4},
    {1.26,	275.8,	5.1},
    {2.58,	86.6,	2.1},
    {17.34,	53.7,	34.8},
    {4.82,	93.4,	39.7},
    {5.19,	298.6,	38.7},
    {37.73,	259.3,	-2.2} },

    { {8.68,	9.9,	-3.3},
    {1.58,	255.9,	-74.9},
    {6.90,	279.5,	81.2},
    {3.53,	211.8,	-2.8},
    {0.34,	354.2,	13.3} },


    { {11.50, 	259.2,	-2.1},
    {2.61,	59.5,	-88.9},
    {6.70,	281.8,	73.9},
    {0.56,	162.0,	5.4},
    {6.35,	340.3,	7.3} },

    { {3.85,	280.8,	-14.1},
    {2.85,	222.0,	-85.6},
    {6.05,	291.3,	77.6},
    {1.39,	182.1,	3.6},
    {3.55,	335.2,	-3.1},
    {6.76,	270.9,	7.2} },



    { {39.01, 	250.8,	-8.2},
    {0.27,	203.1,	-82.8},
    {1.87,	278.5,	75.4},
    {3.54,	155.2,	-0.7},
    {2.92,	316.0,	1.9},
    {1.83,	169.7,	-35.0},
    {9.09,	219.1,	-11.3},
    {3.42,	237.9,	14.1},
    {7.07,	258.7,	14.0},
    {6.19,	25.2,	20.6},
    {4.40,	113.1,	15.1},
    {9.52,	322.2,	49.4} },

    { {7.78,	269.7,	7.5},
    {1.20,	28.7,	-84.0},
    {0.49,	59.0,	86.5},
    {0.66,	230.7,	9.5},
    {0.52,	29.0,	10.0},
    {7.62,	245.1,	3.6} },

    { {12.74, 	325.7,	18.3},
    {0.76,	245.6,	-78.3},
    {2.92,	155.6,	87.0},
    {4.68,	172.9,	-4.2},
    {1.33,	335.6,	-4.0},
    {8.68,	356.7,	48.9},
    {8.07,	43.2,	36.6},
    {2.30,	71.0,	14.6},
    {5.54,	200.6,	7.0},
    {5.09,	258.6,	9.8} },

    { {7.16,	245.2,	3.7},
    {0.53,	205.0,	-86.2},
    {2.58,	44.6,	78.3},
    {5.29,	153.0,  -2.1},
    {0.75,	336.0,	9.0},
    {8.72,	164.3,	-24.7},
    {7.45,	214.8,	0.5},
    {5.16,	251.5,	38.2} },


    { {11.64, 	272.2,	7.1},
    {0.17,	157.4,	-71.1},
    {5.36,	282.4,	87.5},
    {1.98,	159.3,	3.8},
    {3.76,	345.0,	6.2},
    {7.39,	187.4,	6.4},
    {12.42,	187.8,	21.8},
    {10.96,	208.4,	19.9},
    {4.12,	227.3,	16.0},
    {7.91,	232.5,	31.9},
    {6.88,	266.1,	59.1},
    {11.15,	283.1,	70.9} },


    { {16.97,	282.0,	37.0},
    {0.20,	271.1,	-68.2},
    {11.83,	263.3,	88.5},
    {3.96,	216.7,	-2.1},
    {1.62,	44.0,	7.6} },

    { {11.43, 	279.7,	70.7},
    {1.83,	303.4,	-65.5},
    {0.26,	159.7,	77.1},
    {7.48,	143.0,	-3.0},
    {0.16,	327.1,	12.3},
    {7.05,	207.3,	-2.8},
    {6.47,	165.3,	-6.3},
    {7.56,	95.0,	-9.3}}*/
    };






    vector< vector <vector<double>>> Dots = Dots_interpenter(picket);
    vector< vector <vector<double>>> Temp;
    double S = 0;
    double V = 0;
    int a = 20;
    int b = 21;
    vector<double> Sa;
    vector<double> Vs;
    double S1_sum = 0;
    double S2_sum = 0;
    double V1_sum = 0;
    double V2_sum = 0;
    vector<double> h = { 33, 16, 7.5, 3, 4.3, 3, 4, 4, 5, 2, 4, 3, 4, 5, 5, 6, 12, 2 };
    vector<double> s = { 936, 416, 52, 30, 80, 38, 70, 132, 72, 60, 25, 24, 16, 56,  40 ,64 ,32, 28 };
    vector<double> S_p = { 0,0,0,0 };
    vector<double> V_p = { 0,0,0,0,0,0,0 };
    vector<double> s1;
    vector<double> s2;
    vector<double> s3;
    vector<double> s4;
    vector<double> v1;
    vector<double> v2;
    vector<double> v3;
    vector<double> v4;
    vector<double> v5;
    vector<double> v6;
    vector<double> v7;
    int n = 0;
    for (int i = 0; i < Dots.size() - 1; i++)
    {
        cout << " Пикеты  номер   " << a << "   и " << "   " << b << "   " << "\n";
        Temp = { Dots[i], Dots[i + 1] };
        S = Cave_S(Temp);
        V = Cave_V3elN1(Temp);
        V1_sum += V;
        cout << "Площадь участка пещеры : " << "\n";
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        s1.push_back(S);
        S_p[0] += S / s[i];
        S = Cave_S2(Temp);
        s2.push_back(S);
        S1_sum += s[i];
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        S_p[1] += S / s[i];
        S = Cave_S2_2(Temp);
        s3.push_back(S);
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        S_p[2] += S / s[i];
        S = Cave_S2_3(Temp);
        s4.push_back(S);
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        S_p[3] += S / s[i];

        cout << s[i] << " метров" << "   " << " в процентах  " << "\n" << "\n";
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n";
        V_p[0] += V / (s[i] * 0.78 * h[i]);
        v1.push_back(V);
        V = Cave_V3elN2(Temp);
        v2.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[1] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V2(Temp);
        v3.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[2] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V2_2(Temp);
        v4.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[3] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V(Temp);
        v5.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[4] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V_BN1(Temp);
        v6.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[5] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V_BN2(Temp);
        v7.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        cout << "В процентах от примерного значения:  " << (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[6] += V / (s[i] * 0.78 * h[i]);
        V2_sum += (s[i] * 0.78 * h[i]);
        a++; b++; n++;
    }
    cout << "\n";
    cout << "\n";
    cout << "Площадь пещеры: ";
    cout << S1_sum << " метров" << "    ";
    cout << "\n";
    S = Cave_S(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[0] / n << "   " << "\n";
    S = Cave_S2(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[1] / n << "   " << "\n";
    S = Cave_S2_2(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[2] / n << "   " << "\n";
    S = Cave_S2_3(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[3] / n << "   " << "\n";

    cout << "\n";
    cout << "\n";
    cout << "Объём пещеры: ";
    cout << V2_sum << " метров" << "    " << "\n";
    V = Cave_V3elN1(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[0] / n << "   " << "\n";
    V = Cave_V3elN2(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[1] / n << "   " << "\n";
    V = Cave_V2(Dots);
    cout << "Объём всей пещеры: ";
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[2] / n << "   " << "\n";
    V = Cave_V2_2(Dots);
    cout << "Объём всей пещеры: ";
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[3] / n << "   " << "\n";
    V = Cave_V(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[4] / n << "   " << "\n";
    V = Cave_V_BN1(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[5] / n << "   " << "\n";
    V = Cave_V_BN2(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[6] / n << "   " << "\n";

    for (int i = 0; i < s.size(); i++)
    {
        cout << s[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s1[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s2[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s3[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s4[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";


    for (int i = 0; i < s.size(); i++)
    {
        cout << s[i]*h[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v1[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v2[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v3[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v4[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v5[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v6[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v7[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    s1.clear();
    s2.clear();
    s3.clear();
    s4.clear();
    v1.clear();
    v2.clear();
    v3.clear();
    v4.clear();
    v5.clear();
    v6.clear();
    v7.clear();
    Dots = Dots_interpenter(picket2);
    h = { 14, 12, 12, 6, 16, 6, 6, 4, 4, 6, 8, 4, 3, 6, 8, };
    s = { 62, 90, 102, 42, 328, 42, 156, 180, 256 };
    a = 4;
    b = 5;
    S1_sum = 0;
    V2_sum = 0;
    n = 0;
    S_p = { 0,0,0,0 };
    V_p = { 0,0,0,0,0,0,0 };
    for (int i = 0; i < Dots.size() - 1; i++)
    {
        cout << " Пикеты  номер   " << a << "   и " << "   " << b << "   " << "\n";
        Temp = { Dots[i], Dots[i + 1] };
        S = Cave_S(Temp);
        V = Cave_V3elN1(Temp);
        V1_sum += V;
        cout << "Площадь участка пещеры : " << "\n";
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        s1.push_back(S);
        S_p[0] += S / s[i];
        S = Cave_S2(Temp);
        s2.push_back(S);
        S1_sum += s[i];
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        S_p[1] += S / s[i];
        S = Cave_S2_2(Temp);
        s3.push_back(S);
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        S_p[2] += S / s[i];
        S = Cave_S2_3(Temp);
        s4.push_back(S);
        cout << S << " метров" << "   " << " в процентах  " << S / s[i] << "   " << "\n";
        S_p[3] += S / s[i];

        cout << s[i] << " метров" << "   " << " в процентах  " << "\n" << "\n";
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n";
        V_p[0] += V / (s[i] * 0.78 * h[i]);
        v1.push_back(V);
        V = Cave_V3elN2(Temp);
        v2.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[1] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V2(Temp);
        v3.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[2] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V2_2(Temp);
        v4.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[3] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V(Temp);
        v5.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[4] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V_BN1(Temp);
        v6.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[5] += V / (s[i] * 0.78 * h[i]);
        V = Cave_V_BN2(Temp);
        v7.push_back(V);
        cout << "Объём участка пещеры: ";
        cout << V << " кубических метров." << "     ";
        cout << "В процентах от примерного значения:  " << V / (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        cout << "В процентах от примерного значения:  " << (s[i] * 0.78 * h[i]) << "   " << "\n" << "\n";
        V_p[6] += V / (s[i] * 0.78 * h[i]);
        V2_sum += (s[i] * 0.78 * h[i]);
        a++; b++; n++;
    }
    cout << "\n";
    cout << "\n";
    cout << "Площадь пещеры: ";
    cout << S1_sum << " метров" << "    ";
    cout << "\n";
    S = Cave_S(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[0] / n << "   " << "\n";
    S = Cave_S2(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[1] / n << "   " << "\n";
    S = Cave_S2_2(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[2] / n << "   " << "\n";
    S = Cave_S2_3(Dots);
    cout << "Площадь всей пещеры : " << "\n";
    cout << S << " метров" << "   " << " в процентах  " << S / (S1_sum) << "   " << "точность метода  " << S_p[3] / n << "   " << "\n";

    cout << "\n";
    cout << "\n";
    cout << "Объём пещеры: ";
    cout << V2_sum << " метров" << "    " << "\n";
    V = Cave_V3elN1(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[0] / n << "   " << "\n";
    V = Cave_V3elN2(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[1] / n << "   " << "\n";
    V = Cave_V2(Dots);
    cout << "Объём всей пещеры: ";
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[2] / n << "   " << "\n";
    V = Cave_V2_2(Dots);
    cout << "Объём всей пещеры: ";
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[3] / n << "   " << "\n";
    V = Cave_V(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[4] / n << "   " << "\n";
    V = Cave_V_BN1(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[5] / n << "   " << "\n";
    V = Cave_V_BN2(Dots);
    cout << V << " кубических метров" << "    " << " в процентах  " << V / V2_sum << "   " << "точность метода  " << V_p[6] / n << "   " << "\n";
    

    for (int i = 0; i < s.size(); i++)
    {
        cout << s[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s1[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s2[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s3[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << s4[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";


    for (int i = 0; i < s.size(); i++)
    {
        cout << s[i] * h[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v1[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v2[i] << "," << "    ";
    }
    cout << "\n"; cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v3[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v4[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v5[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v6[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";
    for (int i = 0; i < s.size(); i++)
    {
        cout << v7[i] << "," << "    ";
    }
    cout << "\n";
    cout << "\n";

}

