//  do_nm.hpp
//  Numerical_Method_Implicit
//  Created by Dongliang Yi on 4/15/16.
//  Copyright © 2016 Dongliang Yi. All rights reserved.

double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T);
double ** do_euler_implicitet(int N, int M, float r, double sigma, double K, double S_max, double T);
double ** do_crank_nicolson(int N, int M, float r, double sigma, double K, double S_max, double T);
