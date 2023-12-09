#ifndef DRIVEN_FORCE_H
#define DRIVEN_FORCE_H


#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <linalg.h>




template <typename T, unsigned Dim>
using Coordinates = linalg::vec<T, Dim>;




//--------------------------------------------------DrivenForce-------------------------------------------------------------------

/**
 * @brief class DrivenForce - class for represent driven force with amplitude F and frequency W. 
 *              
 */
template <typename T>
class DrivenForce
{
	T F_, W_;
	T (*Func_)(T F, T W, Coordinates<T, 3> State);

public:
	DrivenForce(T F, T W, T (*Func)(T F, T W, Coordinates<T, 3> State)) : F_(F), W_(W), Func_(Func) {};
	T getF() const { return F_; }
	T getW() const { return W_; }
	T operator()(Coordinates<T, 3> State) const { return Func_(F_, W_, State); }
};



#endif // DRIVEN_FORCE_H
