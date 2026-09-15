/******************************************************************************
 *  VIRTEx — Viral Infection, Recovery, and Transmission Explorer
 *  https://github.com/cooperg456/VIRTEx
 *
 *  A set of parallel tools for analyzing Continuous-Time Markov Chains (CTMCs)
 *  derived from the Chemical Master Equation (CME).
 *
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "CLITool.hpp"
#include "FileIO.hpp"
#include "SimBuilder.hpp"
#include "SimRunner.hpp"

#include <thread>

/******************************************************************************
 *  Entry point
 ******************************************************************************/
 
int main(int argc, char *argv[]) {
    CLITool::Arguments args = CLITool::parseArguments(argc, argv);

    auto model = FileIO::loadModel_JSON(args.inputFile);
    auto sim = SimRunner(SimBuilder::BuildSweep(model, args.mods), args.bounds, args.analysis, args.warps);
    sim.ssa(args.seed);

    if (args.watch) {
        const std::string title = "Progress";
        auto pb = CLITool::ProgressBar(progress, static_cast<int>(args.mods.rates.size()) * args.warps * 32, &title);
        while (pb.show()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
    }

    sim.sync();

    //  TODO:   process and output

    return 0;
}
    