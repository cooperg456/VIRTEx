/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "CLIAdapter.hpp"

#include "CLI/CLI11.hpp"

#include <random>

/******************************************************************************
 *  Parse arguments
 ******************************************************************************/

CLIAdapter::Arguments CLIAdapter::parseArguments(int argc, char *argv[]) {
    CLI::App app{"VIRTEx — Viral Infection, Recovery, and Transmission Explorer"};

    std::random_device rd;

    Arguments args{};
    args.seed = rd();

    app.add_option("input", args.inputFile,
        "Specify the reaction network JSON file")
        ->required();

    app.add_option("-W,--warps", args.warps,
        "Number of sample path warps (32 paths per warp)")
        ->required()
        ->check(CLI::NonNegativeNumber);

    app.add_option("-T,--tMax", args.stoppingConditions.tMax,
        "Maximum simulation time")
        ->required()
        ->check(CLI::NonNegativeNumber);

    app.add_option("-B,--bound", args.stoppingConditions.bounds,
        "System exit bounds"
        "\n--bound H = 1000 V = 0"
        "\n--bound H < 1");

    app.add_option("-S,--seed", args.seed,
        "Base seed for initializing cuRAND")
        ->check(CLI::NonNegativeNumber);

    app.add_option("--rate", args.overrides.rates,
        "Override a reaction's rate. Ex."
        "\n--rate GammaV 15"
        "\n--rate GammaV $(seq 1.5 1.5 150)");

    app.add_option("--reactants", args.overrides.reactants,
        "Override a reaction's reactants. Ex."
        "\n--reactants Gamma 1 0 1"
        "\n--reactants Gamma $(printf '%i 0 1' $(seq 1 10))");

    app.add_option("--products", args.overrides.products,
        "Override a reaction's products. Ex."
        "\n--reactants aI 0 35 0 1"
        "\n--reactants aI $(printf '0 %i 0 1' $(seq 1 50))");

    app.add_option("--ic", args.overrides.initialConditions,
        "Override the initial conditions. Ex."
        "\n--ic 1000 3 0 0"
        "\n--ic $(printf '1000 %i 0 0' $(seq 1 30))")
        ->check(CLI::NonNegativeNumber);

    auto outputDirOpt = app.add_option("-o,--output", args.outputDir,
        "Directory to place simulation output files");

    auto tGridOpt = app.add_option("-t,--tGrid", args.analysis.savePaths.tGrid,
        "Time grid spacing for saved trajectories")
        ->check(CLI::NonNegativeNumber);

    app.add_option("-s,--save", args.analysis.savePaths.saved,
        "Number of paths to save trajectories for")
        ->needs(tGridOpt)
        ->needs(outputDirOpt)
        ->check(CLI::NonNegativeNumber);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        app.exit(e);
        std::exit(0);
    }

    return args;
}
