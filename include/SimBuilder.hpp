/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#pragma once

#include <string>
#include <vector>

namespace SimBuilder {

/******************************************************************************
 *  Sweep builder
 ******************************************************************************/

struct Model {
    std::vector<std::string> reactants{};
    std::vector<std::string> reactions{};

    std::vector<double> rates{};
    std::vector<int> ics{};
    std::vector<int> rcs{};
    std::vector<int> tcs{};
};

struct Modifiers {
    std::vector<std::tuple<std::string, std::vector<double>>> rates{};
    std::vector<std::tuple<std::string, std::vector<int>>> reactants{};
    std::vector<std::tuple<std::string, std::vector<int>>> products{};
    std::vector<int> ics{};
};

struct Bounds {
    double tMax = 0;
    std::vector<std::vector<std::tuple<std::string, std::string, int>>> bounds{};
};

struct Analysis {
    bool paths = false;
    double tGrid = 0;

    bool exits = false;
};

static size_t getIdx(const std::string& str, std::vector<std::string> vec);

std::vector<Model> BuildSweep(const Model& model, const Modifiers& mods);

}   //  SimBuilder
