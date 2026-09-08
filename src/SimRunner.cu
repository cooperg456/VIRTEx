/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "SimRunner.hpp"

#define MAX_SSA_REACTANTS 60    //  sized for CUDA shared memory limits
#define MAX_SSA_REACTIONS 96

#define MAX_SSA_BOUNDARIES 84   //  sized for CUDA constant memory limits

#define SSA_BLOCK_SIZE 32       //  multiple of 32. "it works on my machine"

/******************************************************************************
 *  SimRunner class
 ******************************************************************************/