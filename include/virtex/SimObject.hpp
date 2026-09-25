////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#pragma once

#include <virtex/Model.hpp>

namespace Vx {
    enum class SimType {
        StochasticSimulation
    };

    class SimOutput {
    public:
        SimOutput(Model::Params const &params, int numTrials);

        ~SimOutput();

        SimOutput(const SimOutput &) = delete;

        SimOutput(SimOutput &&) noexcept;

        SimOutput &operator=(const SimOutput &) = delete;

        SimOutput &operator=(SimOutput &&) noexcept;

        std::vector<int> getPaths() const;

        std::vector<int> getExits() const;

        std::vector<double> getTimes() const;

    private:
        friend class SimObject;

        int _numTrials = 0;

        const Model::Params* _params;

        int* d_paths = nullptr;

        int* d_exits = nullptr;
        double* d_times = nullptr;
    };

    class SimObject {
    public:
        SimObject(Model::Params const &params, SimType type);

        ~SimObject();

        SimObject(const SimObject &) = delete;

        SimObject(SimObject &&) noexcept;

        SimObject &operator=(const SimObject &) = delete;

        SimObject &operator=(SimObject &&) noexcept;

        void runSimulation();

        SimOutput* deviceSynchronize();

    private:
        const Model::Params &_params;
        const SimType _type;

        int _numTrials = 1;

        double* d_rates = nullptr;
        int* d_initial = nullptr;
        int* d_alpha = nullptr;
        int* d_trans = nullptr;

        SimOutput* simOutput = nullptr;
    };
}
