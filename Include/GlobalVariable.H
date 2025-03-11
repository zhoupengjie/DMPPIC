#ifndef __GLOBALVARIABLE_H__
#define __GLOBALVARIABLE_H__

#include <random>

#include "Random.H"

std::default_random_engine generator;
std::normal_distribution<scalar> distribution(0.0, 1.0);

Foam::Random rnd;

struct GlobalVariable{
    scalar duffusion_coefficient;
    const scalar kB = 1.38064852e-23;
    const scalar T = 298.15;
} GV;

#endif
