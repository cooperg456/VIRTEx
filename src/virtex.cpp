////////////////////////////////////////////////////////////////////////////////
//  VIRTEx — Viral Infection, Recovery, and Transmission Explorer
//  https://github.com/cooperg456/VIRTEx
//
//  A set of parallel tools for analyzing Continuous-Time Markov Chains (CTMCs)
//  derived from the Chemical Master Equation (CME).
//
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#include <argparse/argparse.hpp>

#include <iostream>
#include <random>

#include <virtex/Model.hpp>
#include <virtex/SimObject.hpp>

#include "Helpers.hpp"

namespace {
    struct VIRTExArgs : argparse::Args {
        std::filesystem::path &model = arg("modelPath", "path to model xml file");
        std::filesystem::path &params = arg("sweepPath", "path to sweep xml file");
        std::filesystem::path &output = kwarg("o,output", "path to output directory").set_default(".");
        bool &verbose = flag("v,verbose", "toggle verbose output");
    };
}

int main(int argc, char *argv[]) {
    auto args = argparse::parse<VIRTExArgs>(argc, argv);

    if (args.verbose) {
        args.print();
    }

    Vx::Model::Model model = Vx::Helpers::parseModel(Vx::Helpers::loadFile(args.model));

    Vx::Model::Params params = Vx::Helpers::parseParams(Vx::Helpers::loadFile(args.params), model);

    if (!params.seed) {
        std::random_device rd;
        params.seed = rd();
    }

    if (args.verbose) {
        std::cout << "\n" << model << "\n" << params << "\n";
    }

    Vx::SimObject sim = Vx::SimObject(params, Vx::SimType::StochasticSimulation);

    sim.runSimulation();
    Vx::SimOutput* out = sim.deviceSynchronize();

    auto paths = out->getPaths();
    for (const auto &path: paths) {
        std::cout << path << "\n";
    }

    return 0;
}
