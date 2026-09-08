/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#pragma once

#include <SimBuilder.hpp>

#include <vector>

/******************************************************************************
 *  SimRunner class
 ******************************************************************************/

class SimRunner {
public:
    SimRunner(const std::vector<SimBuilder::Model>& models, const SimBuilder::Bounds& bounds,
              const SimBuilder::Analysis& analysis, int warps);
    ~SimRunner();

    volatile double* ssa(unsigned long long seed);
    volatile double* sde(unsigned long long seed);
    volatile double* tau(unsigned long long seed);

    void sync() const;

private:

};
