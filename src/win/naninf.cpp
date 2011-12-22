#include <numeric>
#include <limits>
extern "C"
{
#include "../data.h"
}

extern "C"
int isnan(DOUBLE a)
{
    return a!=a;
}

extern "C"
int isinf(DOUBLE a)
{
    return std::numeric_limits<DOUBLE>::infinity() == abs(a);
}
