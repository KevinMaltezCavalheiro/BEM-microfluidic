#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
using namespace std;

int main(int argc, const char * argv[]) 
{
    #pragma omp parallel
    #pragma omp critical
    cout << "Greetings from thread " << omp_get_thread_num() << endl;
    return 0;
}

