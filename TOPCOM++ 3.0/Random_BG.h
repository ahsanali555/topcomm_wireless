#pragma once
#include <cstdlib>
// #include <iostream>

class Random{
 private:

  
 public:
  unsigned long int seed_u;
  unsigned long int seed;  //previously LONG INT

  Random(void) {
    this->seed=987654321u;
    this->seed_u=123456789lu;
  }
  ~Random(void){;}
  void bubbleSort(int a[], int size);
  double gauss(double sdev, double mean);
  double uniform(double a, double b);
  int uniform(int  a, int b); // [a, b)
  int nonUniform(int  a, int b);

}; 
