//
//  do_nm.cpp
//  Numerical_Method_Implicit
//
//  Created by Dongliang Yi on 4/15/16.
//  Copyright © 2016 Dongliang Yi. All rights reserved.
//


#include <iostream>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <gsl/gsl_linalg.h>


double ** do_euler_implicitet(int N, int M, float r, double sigma, double K, double S_max, double T)
{
    
    double ** result_table;
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    
    
     // Build Example Matrix
     gsl_matrix * t_m = gsl_matrix_alloc(N,N);
     gsl_matrix_set_zero(t_m);
     for(int i = 0; i < N; i++)
     {
     
     double B_i1 = (1+ sigma*sigma*(i+1)*(i+1)*T/M +r*T/M);
     gsl_matrix_set(t_m, i, i, B_i1);
     }
    
     for(int i = 0; i < N-1; i++)
     {
     double A_i1 = 0.5*r*(i+2)*T/M - 0.5*sigma*sigma*(i+2)*(i+2)*T/M;
     gsl_matrix_set(t_m, i+1, i, A_i1);
     }
    
     for(int i = 0; i < N-1; i++){
     
     double C_i1 = -0.5*r*(i+1)*T/M - 0.5*sigma*sigma*(i+1)*(i+1)*T/M;
     gsl_matrix_set(t_m, i, i+1,C_i1);
     }
    
    
     //Inverse
     std::cout << "::::::::::::::::::::::::::"<<std::endl;
     for(int i = 0; i < N; i++)
     {
     for(int j = 0; j < N; j++)
     std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
     std::cout<<std::endl;
     }
     
     std::cout << "::::::::::::::::::::::::::"<<std::endl;
     
     gsl_matrix *t_inverse = gsl_matrix_alloc(N,N);
     int s;
     gsl_permutation *p = gsl_permutation_alloc(N);
     gsl_linalg_LU_decomp(t_m, p, &s);
     gsl_linalg_LU_invert(t_m, p, t_inverse);
     
    
     for(int i = 0; i < N; i++)
     {
     for(int j = 0; j < N; j++)
     std::cout<<gsl_matrix_get(t_inverse, i, j)<<" ";
     std::cout<<std::endl;
     }
     
     // Solution Vector
     gsl_vector * s_v_old = gsl_vector_alloc(N);
     gsl_vector * s_v_new = gsl_vector_alloc(N);
     
     gsl_vector_set_zero(s_v_old);
     gsl_vector_set_zero(s_v_new);
     
     //example initial value
     
     for (int i = 0; i<N; i++) {
     
     double intrisic_value = ((i+1)*S_max/N - K) > 0? ((i+1)*S_max/N - K): 0;// this is the end value
     
     gsl_vector_set(s_v_old, i, intrisic_value);
     }
     
    
     std::cout<<"New\n";
     for(int i = 0; i < N; i++)
     std::cout<<gsl_vector_get(s_v_new, i)<<" ";
     std::cout<<std::endl;
     
     std::cout<<"Old\n";
     for(int i = 0; i < N; i++)
     std::cout<<gsl_vector_get(s_v_old, i)<<" ";
     std::cout<<std::endl;
     
     // printing out the initial condition
     std::cout<<"Starting\n";
     for(int i = 0; i < N; i++)
     std::cout<<gsl_vector_get(s_v_old, i)<<" ";
     std::cout<<std::endl;
     
     // this is the loop that steps through time
     for(int t_j = 0; t_j < M; t_j++)
     {
         double S_NP = S_max-K*exp(-r*(t_j+1)*T/M)+S_max/N;
         double C_i1 = -0.5*r*(N)*T/M - 0.5*sigma*sigma*(N)*(N)*T/M;
         double S_N = gsl_vector_get(s_v_old, N-1);
         double adjust_S_N = S_N - C_i1*S_NP;
         
         gsl_vector_set(s_v_old, N-1, adjust_S_N);

         
    // This is what actually does the matrix vector multiplication.
     gsl_blas_dgemv(CblasNoTrans, 1.0, t_inverse, s_v_old, 0.0, s_v_new);
     
     
     for(int i = 0; i < N; i++)
     std::cout<<gsl_vector_get(s_v_new, i)<<" ";
     std::cout<<std::endl;
     
     memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
     
     gsl_vector_swap(s_v_new, s_v_old);
     }
    
    
    //gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    
    
    return result_table;
};






double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T)
{
    
    double ** result_table;
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    
    // Build Example Matrix
    gsl_matrix * t_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m);
    for(int i = 0; i < N; i++)
    {
        
        double B_i1 = (1- sigma*sigma*(i+1)*(i+1)*T/M)/(1+r*T/M);
        gsl_matrix_set(t_m, i, i, B_i1);
    }
    
    
    
    for(int i = 0; i < N-1; i++)
    {
        double A_i1 = (-0.5*r*(i+2)*T/M + 0.5*sigma*sigma*(i+2)*(i+2)*T/M)/(1+r*T/M);
        gsl_matrix_set(t_m, i+1, i, A_i1);
    }
    
    
    for(int i = 0; i < N-1; i++){
        
        double C_i1 = (0.5*r*(i+1)*T/M + 0.5*sigma*sigma*(i+1)*(i+1)*T/M)/(1+r*T/M);
        gsl_matrix_set(t_m, i, i+1,C_i1);
    }
    
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
        std::cout<<std::endl;
    }
    
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    
    //example initial value
    
    for (int i = 0; i<N; i++) {
        
        double intrisic_value = ((i+1)*S_max/N - K) > 0? ((i+1)*S_max/N - K): 0;
        
        gsl_vector_set(s_v_old, i, intrisic_value);
    }
    
    
    //gsl_vector_set(s_v_old, N/2, 2.5);
    //gsl_vector_set(s_v_old, N/2 - 1, 1.5);
    //gsl_vector_set(s_v_old, N/2 + 1, 1.5);
    
    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // printing out the initial condition
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // this is the loop that steps through time
    for(int t_j = 0; t_j < M; t_j++)
    {
        // This is what actually does the matrix vector multiplication.
        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
        
        //adjust
        double S_NM1 = gsl_vector_get(s_v_new, N-2);
        double S_NM2 = gsl_vector_get(s_v_new, N-3);
        double S_N = (S_NM1*2-S_NM2) > 0? (S_NM1*2-S_NM2):0;
        gsl_vector_set(s_v_new, N-1, S_N);
        
        
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
    }
    
    gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    
    return result_table;

};


double ** do_crank_nicolson(int N, int M, float r, double sigma, double K, double S_max, double T)
{
    double ** result_table;
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    
    
    // Build Example Matrix
    gsl_matrix * t_m1 = gsl_matrix_alloc(N,N);
    gsl_matrix * t_m2 = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m1);
    gsl_matrix_set_zero(t_m2);
    
    for(int i = 0; i < N; i++)
    {
        
        double B_i1 = -0.5*T/M *(sigma*sigma*(i+1)*(i+1)+r);
        gsl_matrix_set(t_m1, i, i, 1-B_i1);
        gsl_matrix_set(t_m2, i, i, 1+B_i1);
    }
    
    
    
    for(int i = 0; i < N-1; i++)
    {
        double A_i1 = 0.25*T/M * (sigma*sigma*(i+2)*(i+2)-r*(i+2));
        gsl_matrix_set(t_m1, i+1, i, -A_i1);
        gsl_matrix_set(t_m2, i+1, i, A_i1);
    }
    
    
    for(int i = 0; i < N-1; i++){
        
        double C_i1 = 0.25*T/M * (sigma*sigma*(i+1)*(i+1)+r*(i+1));
        gsl_matrix_set(t_m1, i, i+1,-C_i1);
        gsl_matrix_set(t_m2, i, i+1,C_i1);
    }
    
    
    
    //Inverse
    std::cout << "M1 Matrix: "<<std::endl;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m1, i, j)<<" ";
        std::cout<<std::endl;
    }
    
    std::cout << "-----------------"<<std::endl<<std::endl;
    
    
    
    std::cout << "M2 Matrix: "<<std::endl;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m2, i, j)<<" ";
        std::cout<<std::endl;
    }
    
    std::cout << "-----------------"<<std::endl<<std::endl;
    
    
    
    gsl_matrix *t_inverse = gsl_matrix_alloc(N,N);
    int s;
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(t_m1, p, &s);
    gsl_linalg_LU_invert(t_m1, p, t_inverse);
    
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_int = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    gsl_vector_set_zero(s_v_int);
    
    //example initial value
    
    for (int i = 0; i<N; i++) {
        
        double intrisic_value = ((i+1)*S_max/N - K) > 0? ((i+1)*S_max/N - K): 0;// this is the end value
        
        gsl_vector_set(s_v_old, i, intrisic_value);
    }
    
    
    //gsl_vector_set(s_v_old, N/2, 2.5);
    //gsl_vector_set(s_v_old, N/2 - 1, 1.5);
    //gsl_vector_set(s_v_old, N/2 + 1, 1.5);
    
    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // printing out the initial condition
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // this is the loop that steps through time
    for(int t_j = 0; t_j < M; t_j++)
    {
        
        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m2, s_v_old, 0.0, s_v_int);
        
        
        
        
        double C_i1 = 0.25*T/M * (sigma*sigma*(N)*(N)+r*(N));
        
        double S_NP = S_max-K*exp(-r*(t_j+1)*T/M)+S_max/N + S_max-K*exp(-r*(t_j+1-1)*T/M)+S_max/N;;
        double S_N = gsl_vector_get(s_v_int, N-1);
        double adjust_S_N = S_N + C_i1*S_NP;
        
        
        
        gsl_vector_set(s_v_int, N-1, adjust_S_N);
        
        
        // This is what actually does the matrix vector multiplication.
        gsl_blas_dgemv(CblasNoTrans, 1.0, t_inverse, s_v_int, 0.0, s_v_new);
        
        //adjust
        //double S_NM1 = gsl_vector_get(s_v_new, N-2);
        //double S_NM2 = gsl_vector_get(s_v_new, N-3);
        
        
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
    }
    
    
    
    //gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    
    
    return result_table;

};
