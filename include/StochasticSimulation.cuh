/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#pragma once

#define MAX_SSA_REACTANTS 60    //  sized for CUDA shared memory limits
#define MAX_SSA_REACTIONS 96
#define SSA_BLOCK_SIZE 32

#include "ReactionNetwork.hpp"

/******************************************************************************
 *  SSASysInfo struct
 ******************************************************************************/

struct SSASysInfo {
    std::vector<double> reactionRates{};
    std::vector<int> initialConditions{};
    std::vector<int> reactantCoefficients{};
    std::vector<int> transitionCoefficients{};    
};

struct SSASimInfo {
    double tMax = 0;
    int warps = 0;
    
    double tGrid = 0;
    int savedPaths = 0;
    std::filesystem::path outputDir{};

    ReactionNetwork base{};
};

/******************************************************************************
 *  Stochastic simulation algorithm
 ******************************************************************************/

void SSA(const SSASimInfo& simInfo, const std::vector<SSASysInfo> &sysInfos);
