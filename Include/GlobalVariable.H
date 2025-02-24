#ifndef __GLOBALVARIABLE_H__
#define __GLOBALVARIABLE_H__

// #include <random>

#include "Random.H"

// std::default_random_engine generator;
// std::normal_distribution<scalar> distribution(0.0, 1.0);

Foam::Random rnd;  // Create a local random number generator

struct GlobalVariable{
    scalar duffusion_coefficient;
} GV;

#endif
