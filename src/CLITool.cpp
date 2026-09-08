/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "CLITool.hpp"

#include "CLI/CLI11.hpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <random>

namespace CLITool {

/******************************************************************************
 *  progressBar class
 ******************************************************************************/

ProgressBar::ProgressBar(volatile double* progress, const int nProgress, const std::string* title) {
    constexpr int offset = 20;
    constexpr int width  = 80;

    _progress  = progress;
    _nProgress = nProgress;

    if (title) {
        _title = *title + ": ";
    }

    _barLength = width - static_cast<int>(_title.length()) - offset;
}

ProgressBar::~ProgressBar() {
    std::cout << std::endl;
}

bool ProgressBar::show() {
    if (_start == std::chrono::steady_clock::time_point{}) {
        _start =  std::chrono::steady_clock::now();
    }

    double minProgress = *std::ranges::min_element(_progress, _progress + _nProgress);
    minProgress = std::min(1.0, std::max(static_cast<double>(minProgress), 0.0));
    const int minBarProgress = static_cast<int>(_barLength * minProgress);

    double avgProgress = 0;
    for (int i = 0; i < _nProgress; i++) {
        avgProgress += _progress[i] / _nProgress;
    }
    avgProgress = std::min(1.0, std::max(static_cast<double>(avgProgress), 0.0));
    int avgBarProgress = static_cast<int>(_barLength * avgProgress);
    const double avgBarFraction = _barLength * avgProgress - avgBarProgress;

    std::string bar;
    for (int i = 0; i < minBarProgress; i++) {
        bar += "\u2591";
    }
    for (int i = 0; i < avgBarProgress - minBarProgress; i++) {
        bar += "\u2588";
    }
    if (avgBarFraction > 0.0) {
        bar += _barChars[static_cast<size_t>(avgBarFraction * 8)];
        avgBarProgress++;
    }
    for (int i = 0; i < _barLength - avgBarProgress; i++) {
        bar += ' ';
    }

    const auto elapsed = duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - _start);
    const std::chrono::hh_mm_ss time{elapsed};

    std::cout << "\r" << _title
              << std::setw(5) << std::fixed << std::setprecision(1) << (avgProgress * 100.0)
              << "% |" << bar << "| [" << time << "]" << std::flush;

    return *std::ranges::min_element(_progress, _progress + _nProgress) < 1.0;
}

/******************************************************************************
 *  CLI argument parser
 ******************************************************************************/

Arguments parseArguments(int argc, char *argv[]) {
    CLI::App app{"VIRTEx — Viral Infection, Recovery, and Transmission Explorer"};

    Arguments args{};

    std::random_device rd;
    args.seed = rd();

    app.add_option("input", args.inputFile,
        "Specify the reaction network JSON file")
        ->required();

    auto outputDirOpt = app.add_option("-o,--output", args.outputDir,
        "Directory to place simulation output files");

    app.add_option("-W,--warps", args.warps,
        "Number of sample path warps (32 paths per warp)")
        ->check(CLI::NonNegativeNumber);

    app.add_option("-S,--seed", args.seed,
        "Base seed for initializing cuRAND")
        ->check(CLI::NonNegativeNumber);

    app.add_flag("-v,--verbose", args.verbose);

    app.add_flag("-w,--watch", args.watch);

    app.add_option("-T,--tMax", args.bounds.tMax,
        "Maximum simulation time")
        ->check(CLI::NonNegativeNumber);

    auto boundOpt = app.add_option("-B,--bound", args.bounds.bounds,
        "System exit bounds. Ex."
        "\n--bound H = 1000 V = 0"
        "\n--bound H < 1");

    app.add_option("--rate", args.mods.rates,
        "Override a reaction's rate. Ex."
        "\n--rate GammaV 15"
        "\n--rate GammaV $(seq 1.5 1.5 150)");

    app.add_option("--reactants", args.mods.reactants,
        "Override a reaction's reactants. Ex."
        "\n--reactants Gamma 1 0 1"
        "\n--reactants Gamma $(printf '%i 0 1' $(seq 1 10))");

    app.add_option("--products", args.mods.products,
        "Override a reaction's products. Ex."
        "\n--reactants aI 0 35 0 1"
        "\n--reactants aI $(printf '0 %i 0 1' $(seq 1 50))");

    app.add_option("--ic", args.mods.ics,
        "Override the initial conditions. Ex."
        "\n--ic 1000 3 0 0"
        "\n--ic $(printf '1000 %i 0 0' $(seq 1 30))")
        ->check(CLI::NonNegativeNumber);

    auto tGridOpt = app.add_option("-t,--tGrid", args.analysis.tGrid,
        "Time grid spacing for saved trajectories")
        ->check(CLI::NonNegativeNumber);

    app.add_flag("-s,--samplePaths", args.analysis.paths,
        "Number of paths to save trajectories for")
        ->needs(tGridOpt)
        ->needs(outputDirOpt);

    app.add_flag("-e,--extinctions", args.analysis.exits)
        ->needs(boundOpt)
        ->needs(outputDirOpt);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        app.exit(e);
        std::exit(0);
    }

    return args;
}

}   //  CLITool