// Layers_MC.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include "Layer_Reflection.h"

using namespace std;

int main()
{
		
	Layer_Reflection layer;	
		
	cout<<"Running..."<<"\n"; 
	//double start = omp_get_wtime( );
	layer.Calculate();
	//double end = omp_get_wtime( );
	//cout<<"Time cost ："<<end -start<<"\n";
	return 0;
}

