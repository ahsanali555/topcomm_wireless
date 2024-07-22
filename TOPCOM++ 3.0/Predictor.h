#pragma once
#include <stdio.h>
class Predictor
{
public:
	Predictor();
	~Predictor();

	void Display(FILE * file) const;

	//! Insert a new dat in table
	void Insert(const double x, const double y);
	void Reset() 
	{
		ltab = 0;
	} 

	//! Insert a new dat in table
	double Predict(const double x) const;

private:
	double* xt;
	double* yt;
	int ltab;
};

