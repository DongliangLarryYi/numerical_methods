//  main.cpp
//  Numerical_Method_Implicit
//  Created by Dongliang Yi on 4/15/16.
//  Copyright © 2016 Dongliang Yi. All rights reserved.

#include <iostream>
#include <fstream>
#include "do_nm.hpp"

int main(void)
{
    std::ofstream output_file,output_file2;
    //output_file.open("out.csv");
    output_file2.open("out2.csv");
    int N = 20; // stock price split
    int M = 20; //time split
    double r = 0.02;
    double sigma = .2;
    double K = 20;
    double S_max = 100;
    double T = 5; // for simplicity, I shortened the maturity. since there is no 20 year option in real world.
    
    double ** result_table2 = NULL;
    result_table2 = do_euler_implicitet(N, M, r, sigma, K, S_max, T);
    std::cout<<"The Result!\n";
    for(int i =0; i < M; i++)
    {
         for(int j = 0; j<N-1; j++)
        {
            std::cout<<result_table2[i][j]<<" ";
            output_file2<<result_table2[i][j]<<',';
        }
        std::cout<<result_table2[i][N-1]<<std::endl;
        output_file2<<result_table2[i][N-1]<<"\n";
    }
}
