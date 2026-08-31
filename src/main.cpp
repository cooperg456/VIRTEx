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

#include "CLIAdapter.hpp"
#include "ReactionNetwork.hpp"
#include "StochasticSimulation.cuh"

#include <random>
#include <tuple>

/******************************************************************************
 *  Entry point
 ******************************************************************************/
 
int main(int argc, char *argv[]) {
    CLIAdapter::Arguments args = CLIAdapter::parseArguments(argc, argv);

    /**************************************************************************
     *  Input checking
     **************************************************************************/

    ReactionNetwork ctmc(args.inputFile);

    std::vector<int> sweepSizes;
    size_t reactantsSize = ctmc.reactants.size();
    sweepSizes.reserve(args.overrides.rates.size());
    for (auto & rate : args.overrides.rates) {
           sweepSizes.push_back(static_cast<int>(std::get<1>(rate).size()));
        }

    for (auto & react : args.overrides.reactants) {
        int sweepSize = static_cast<int>(std::get<1>(react).size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error(
                "Length of required reactants must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    for (auto & prod : args.overrides.products) {
        int sweepSize = static_cast<int>(std::get<1>(prod).size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error(
                "Length of reaction products must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    if (!args.overrides.initialConditions.empty()) {
        int sweepSize = static_cast<int>(args.overrides.initialConditions.size());
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
    simInfo.seed = args.seed;
    simInfo.tMax = args.stoppingConditions.tMax;
    simInfo.warps = args.warps;
    simInfo.tGrid = args.analysis.savePaths.tGrid;
    simInfo.savedPaths = args.analysis.savePaths.saved;
    simInfo.outputDir = args.outputDir;
    simInfo.base = ctmc;

    for (auto& boundSet : args.stoppingConditions.bounds) {
        for (auto& bound : boundSet) {
            simInfo.boundVals.push_back(static_cast<int>(ctmc.getReactantIdx(std::get<0>(bound))));

            if (std::string& op = std::get<1>(bound); op == "=") {
                simInfo.boundVals.push_back(0);
            }
            else if (op == ">") {
                simInfo.boundVals.push_back(1);
            }
            else if (op == "<") {
                simInfo.boundVals.push_back(2);
            }
            else {
                throw std::runtime_error("Invalid comparison operator: " + op);
            }

            simInfo.boundVals.push_back(std::get<2>(bound));
        }
        simInfo.boundIdxs.push_back(static_cast<int>(simInfo.boundVals.size()));
    }

    std::vector<SSASysInfo> sysInfos(sweepSize);
    for (size_t i = 0; i < sweepSize; i++) {
        sysInfos[i].reactionRates = ctmc.reactionRates;
        sysInfos[i].initialConditions = ctmc.initialConditions;
        sysInfos[i].reactantCoefficients = ctmc.reactantCoefficients;
        sysInfos[i].transitionCoefficients = ctmc.transitionCoefficients;

        for (auto & j : args.overrides.rates) {
            size_t reaction = ctmc.getReactionIdx(std::get<0>(j));
            double rate = std::get<1>(j)[i];

            sysInfos[i].reactionRates[reaction] = rate;
        }

        //  reacts loop must come before prods
        for (auto & j : args.overrides.reactants) {
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

        for (auto & j : args.overrides.products) {
            size_t reaction = ctmc.getReactionIdx(std::get<0>(j));
            std::vector<int>& prod = std::get<1>(j);

            for (size_t k = 0; k < reactantsSize; k++) {
                size_t l = reaction * reactantsSize + k;
                size_t prodIdx = i * reactantsSize + k;

                sysInfos[i].transitionCoefficients[l] = prod[prodIdx] - sysInfos[i].reactantCoefficients[l];
            }
        }

        if (!args.overrides.initialConditions.empty()) {
            for (size_t j = 0; j < reactantsSize; j++) {
                sysInfos[i].initialConditions[j] = args.overrides.initialConditions[i * reactantsSize + j];
            }
        }
    }

    /**************************************************************************
     *  Run and return
     **************************************************************************/

    SSA(simInfo, sysInfos);

    return 0;
}
    