////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#include <virtex/SimObject.hpp>

#include "cuda_runtime.h"
#include "curand_kernel.h"

#define MAX_SSA_REACTANTS 60    //  sized for CUDA shared memory limits
#define MAX_SSA_REACTIONS 96

#define MAX_SSA_BOUNDARIES 84   //  sized for CUDA constant memory limits

#define SSA_BLOCK_SIZE 32       //  multiple of 32. "it works on my machine"

__constant__ static int c_boundIdxs[MAX_SSA_BOUNDARIES];
__constant__ static int c_boundVals[MAX_SSA_REACTANTS * 3];

__device__
static double propensity(const int *x, const int nx, const double *rates, const int nr, const int *alpha, double *a) {
    double a0 = 0;

    for (int i = 0; i < nr; i++) {
        a[i] = rates[i];

        for (int m = 0; m < nx; m++) {
            if (const int alpha_jm = alpha[i * nx + m]; alpha_jm > 0) {
                if (x[m] < alpha_jm) {
                    a[i] = 0.0;
                    break;
                }

                for (int k = x[m] - alpha_jm + 1; k <= x[m]; k++) {
                    a[i] *= k;
                }
            }
        }
        a0 += a[i];
    }
    return a0;
}

__global__
static void stochasticSimulation() {

}

using namespace Vx;

SimObject::SimObject(Model::Params const &params, SimType type) : _params(params), _type(type) {
    int _numBlocksPer = static_cast<int>(params.size / 32) + 1;

    //  TODO:   LOOP

    const size_t reactions = _params.model->reactions.size();
    const size_t species = _params.model->species.size();

    cudaMalloc(&d_rates, _numBlocks * reactions * sizeof(double));
    cudaMalloc(&d_initial, _numBlocks * species * sizeof(int));
    cudaMalloc(&d_alpha, _numBlocks * reactions * species * sizeof(int));
    cudaMalloc(&d_trans, _numBlocks * reactions * species * sizeof(int));

    if (_params.analysis.savePaths.t != 0) {
        int times = static_cast<int>(_params.tMax / _params.analysis.savePaths.t) + 1;
        cudaMalloc(&d_paths, _params.size * species * times * sizeof(int));
    }

    if (!_params.analysis.stopping.bounds.empty()) {
        cudaMalloc(&d_exits, _params.size * species * sizeof(int));
        cudaMalloc(&d_times, _params.size * sizeof(double));

        cudaMemset(d_exits, -1, _params.size * species * sizeof(int));
    }

    //  TODO:   COPY
}

SimObject::~SimObject() {
    cudaFree(d_rates);
    cudaFree(d_initial);
    cudaFree(d_alpha);
    cudaFree(d_trans);
    cudaFree(d_initial);
    cudaFree(d_exits);
    cudaFree(d_times);
}

SimObject::SimObject(SimObject &&other) noexcept : _params(other._params), _type(other._type) {}

SimObject &SimObject::operator=(SimObject &&other) noexcept {
    return *this;
}

void SimObject::runSimulation() {
    switch (_type) {
        case SimType::StochasticSimulation:
            stochasticSimulation<<<_numBlocks, 32>>>();
            break;
        default:
            break;
    }
}

SimOutput SimObject::deviceSynchronize() {
    cudaDeviceSynchronize();

    auto result = SimOutput{};

    return result;
}
