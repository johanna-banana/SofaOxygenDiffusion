//
// Created by Somers on 1/4/2021.
//
#include "SofaOxygenDiffusion.h"

int main()
{
    SofaOxygenDiffusion sofasim = SofaOxygenDiffusion();
    sofasim.init_simulation();
    sofasim.run_simulation();
}