/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#pragma once

#include "SimBuilder.hpp"

#include <filesystem>

namespace CLITool {

/******************************************************************************
 *  ProgressBar class
 ******************************************************************************/

class ProgressBar {
public:
    ProgressBar(volatile double* progress, int nProgress, const std::string* title);
    ~ProgressBar();

    bool show();

private:
    volatile double* _progress;

    int _nProgress;
    int _barLength;

    std::string _title{};

    std::chrono::time_point<std::chrono::steady_clock> _start{};

    static inline const std::string _barChars[8] = {
        " ",      "\u258F",
        "\u258E", "\u258D",
        "\u258C", "\u258B",
        "\u258A", "\u2589"
    };
};

/******************************************************************************
 *  CLI argument parser
 ******************************************************************************/

struct Arguments {
    bool verbose = false;
    bool watch = false;

    std::filesystem::path inputFile{};
    std::filesystem::path outputDir{};

    int warps = 1;
    unsigned long long seed;

    SimBuilder::Bounds bounds{};
    SimBuilder::Modifiers mods{};
    SimBuilder::Analysis analysis{};
};

Arguments parseArguments(int argc, char *argv[]);

}   //  CLITool
