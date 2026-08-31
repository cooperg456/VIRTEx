/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#pragma once

#include <filesystem>
#include <tuple>

/******************************************************************************
 *  CLIAdapter namespace
 ******************************************************************************/

namespace CLIAdapter {

    struct Overrides {
        std::vector<std::tuple<std::string, std::vector<double>>> rates{};
        std::vector<std::tuple<std::string, std::vector<int>>> reactants{};
        std::vector<std::tuple<std::string, std::vector<int>>> products{};
        std::vector<int> initialConditions{};
    };

    struct StoppingConditions {
        double tMax = 0;
        std::vector<std::vector<std::tuple<std::string, std::string, int>>> bounds{};
    };

    struct SavePaths {
        int saved = 0;
        double tGrid = 0;
    };

    struct Analysis {
        SavePaths savePaths{};
    };

    struct Arguments {
        int warps = 0;
        unsigned long long seed;

        Overrides overrides{};
        StoppingConditions stoppingConditions{};
        Analysis analysis{};

        std::filesystem::path inputFile{};
        std::filesystem::path outputDir{};
    };

    [[nodiscard]] Arguments parseArguments(int argc, char *argv[]);

};
