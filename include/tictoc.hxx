// tic toc timing tools
// requires compiler flag -std=c++11
// A.Wingen						18.04.16

// Define
//--------
#ifndef TICTOC_INCLUDED
#define TICTOC_INCLUDED

// Include
//--------
#include <stack>
#include <chrono>
using namespace std::chrono;

// Global
//--------
stack<high_resolution_clock::time_point> tictoc_stack;

// Functions
//----------
inline void tic() {tictoc_stack.push(high_resolution_clock::now());}

inline double toc()
{
double T = duration<double, std::micro>(high_resolution_clock::now() - tictoc_stack.top()).count();
tictoc_stack.pop();
return T*1e-6;
}

#endif //  TICTOC_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
