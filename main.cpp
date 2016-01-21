//
//  main.cpp
//  MC_Layer
//
//  Created by Administrator on 04/12/15.
//  Copyright © 2015 Administrator. All rights reserved.
//


#include <iostream>
#include "Layer_Reflection.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    
    Layer_Reflection layer;
    
    //std::cout<<"Running..."<<"\n";
    //double start = omp_get_wtime( );
    layer.Calculate();
    //double end = omp_get_wtime( );
    //std::cout<<"Time cost £∫"<<end -start<<"\n";
    
    std::cout << "Hello, World!\n";
    return 0;
}
