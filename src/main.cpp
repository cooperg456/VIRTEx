/******************************************************************************
*  CuCTMC -- https://github.com/cooperg456/CuCTMC
 *
 *  A set of parallel tools to analyze Continuous-Time Markov Chains (CTMCs)
 *  of the Chemical Master Equation (CME) type
 ******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "StochasticSimulation.cuh"
#include "ReactionNetwork.hpp"

#include "CLI/CLI11.hpp"

#include <random>
#include <tuple>

/******************************************************************************
 *  Entry point
 ******************************************************************************/
 
int main(int argc, char *argv[]) {
    CLI::App app{"CuCTMC -- Parallel tools for analyzing Chemical Master Equation CTMCs"};

    /**************************************************************************
     *  CLI Arguments
     **************************************************************************/

    std::filesystem::path inputFile;
    app.add_option("input", inputFile,
        "Specify the reaction network JSON file")
        ->required();

    int warps;
    app.add_option("-W,--warps", warps,
        "Number of sample path warps (32 paths per warp)")
        ->required()
        ->check(CLI::NonNegativeNumber);

    double tMax;
    app.add_option("-T,--tMax", tMax,
        "Maximum simulation time")
        ->required()
        ->check(CLI::NonNegativeNumber);

    std::random_device rd;
    unsigned long long seed = rd();
    app.add_option("-S,--seed", seed,
        "Base seed for initializing cuRAND")
        ->check(CLI::NonNegativeNumber);

    std::vector<std::tuple<std::string, std::vector<double>>> rates;
    app.add_option("--rate", rates,
        "Override a reaction's rate. Ex. "
        "\n--rate GammaV 15"
        "\n--rate GammaV $(seq 1.5 1.5 150)");

    std::vector<std::tuple<std::string, std::vector<int>>> reacts;
    app.add_option("--reactants", reacts,
        "Override a reaction's reactants. Ex. "
        "\n--reactants Gamma 1 0 1"
        "\n--reactants Gamma $(printf '%i 0 1' $(seq 1 10))");

    std::vector<std::tuple<std::string, std::vector<int>>> prods;
    app.add_option("--products", prods,
        "Override a reaction's products. Ex. "
        "\n--reactants aI 0 35 0 1"
        "\n--reactants aI $(printf '0 %i 0 1' $(seq 1 50))");

    std::vector<int> initialConditions;
    app.add_option("--ic", initialConditions,
        "Override the initial conditions. Ex. "
        "\n--ic 1000 3 0 0"
        "\n--ic $(printf '1000 %i 0 0' $(seq 1 30))")
        ->check(CLI::NonNegativeNumber);

    std::filesystem::path outputDir{};
    auto outputDirOpt = app.add_option("-o,--output", outputDir,
        "Directory to place simulation output files");

    double tGrid = 0;
    auto tGridOpt = app.add_option("-t,--tGrid", tGrid,
        "Time grid spacing for saved trajectories")
        ->check(CLI::NonNegativeNumber);

    int saved = 0;
    app.add_option("-s,--save", saved,
        "Number of paths to save trajectories for")
        ->needs(tGridOpt)
        ->needs(outputDirOpt)
        ->check(CLI::NonNegativeNumber);

    /**************************************************************************
     *  Input checking
     **************************************************************************/

    CLI11_PARSE(app, argc, argv);

    ReactionNetwork ctmc(inputFile);

    std::vector<int> sweepSizes;
    size_t reactantsSize = ctmc.reactants.size();
    sweepSizes.reserve(rates.size());
    for (auto & rate : rates) {
           sweepSizes.push_back(static_cast<int>(std::get<1>(rate).size()));
        }

    for (auto & react : reacts) {
        int sweepSize = static_cast<int>(std::get<1>(react).size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error(
                "Length of required reactants must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    for (auto & prod : prods) {
        int sweepSize = static_cast<int>(std::get<1>(prod).size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error(
                "Length of reaction products must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    if (!initialConditions.empty()) {
        int sweepSize = static_cast<int>(initialConditions.size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error(
                "Length of initial conditions must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    if (sweepSizes.size() > 1 && std::adjacent_find(
        sweepSizes.begin(), sweepSizes.end(), std::not_equal_to<>()) != sweepSizes.end()) {
        throw std::runtime_error("Length of each sweep must be equal");
    }

    int sweepSize = 1;
    if (!sweepSizes.empty()) {
        sweepSize = sweepSizes.front();
    }

    /**************************************************************************
     *  Create Input Structs
     **************************************************************************/

    SSASimInfo simInfo;
    simInfo.seed = seed;
    simInfo.tMax = tMax;
    simInfo.warps = warps;
    simInfo.tGrid = tGrid;
    simInfo.savedPaths = saved;
    simInfo.outputDir = outputDir; 
    simInfo.base = ctmc;

    std::vector<SSASysInfo> sysInfos(sweepSize);
    for (size_t i = 0; i < sweepSize; i++) {
        sysInfos[i].reactionRates = ctmc.reactionRates;
        sysInfos[i].initialConditions = ctmc.initialConditions;
        sysInfos[i].reactantCoefficients = ctmc.reactantCoefficients;
        sysInfos[i].transitionCoefficients = ctmc.transitionCoefficients;

        for (auto & j : rates) {
            size_t reaction = ctmc.getReactionIdx(std::get<0>(j));
            double rate = std::get<1>(j)[i];

            sysInfos[i].reactionRates[reaction] = rate;
        }

        //  reacts loop must come before prods
        for (auto & j : reacts) {
            size_t reaction = ctmc.getReactionIdx(std::get<0>(j));
            std::vector<int>& react = std::get<1>(j);

            for (size_t k = 0; k < reactantsSize; k++) {
                size_t l = reaction * reactantsSize + k;
                size_t reactIdx = i * reactantsSize + k;

                int beta = ctmc.reactantCoefficients[l] + ctmc.transitionCoefficients[l];

                sysInfos[i].reactantCoefficients[l] = react[reactIdx];
                sysInfos[i].transitionCoefficients[l] = beta - react[reactIdx];
            }
        }

        for (auto & j : prods) {
            size_t reaction = ctmc.getReactionIdx(std::get<0>(j));
            std::vector<int>& prod = std::get<1>(j);

            for (size_t k = 0; k < reactantsSize; k++) {
                size_t l = reaction * reactantsSize + k;
                size_t prodIdx = i * reactantsSize + k;

                sysInfos[i].transitionCoefficients[l] = prod[prodIdx] - sysInfos[i].reactantCoefficients[l];
            }
        }

        if (!initialConditions.empty()) {
            for (size_t j = 0; j < reactantsSize; j++) {
                sysInfos[i].initialConditions[j] = initialConditions[i * reactantsSize + j];
            }
        }
    }

    /**************************************************************************
     *  Run and return
     **************************************************************************/

    SSA(simInfo, sysInfos);

    return 0;
}
    