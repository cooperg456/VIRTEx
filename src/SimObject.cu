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

__constant__
static double c_timeStep;

__constant__
static int c_boundIdxs[MAX_SSA_BOUNDARIES];

__constant__
static int c_boundVals[MAX_SSA_REACTANTS * 3];

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

__device__
static int boundary(const int *x, const int nBounds) {
    for (int s = 0; s < nBounds; s++) {
        const int start = (s == 0) ? 0 : c_boundIdxs[s - 1];
        const int end = c_boundIdxs[s];

        bool match = true;
        for (int k = start; k < end; k += 3) {
            const int m = c_boundVals[k];
            const int op = c_boundVals[k + 1];
            const int value = c_boundVals[k + 2];

            bool ok;
            switch (op) {
                case 0:  ok = (x[m] == value); break;
                case 1:  ok = (x[m] >  value); break;
                case 2:  ok = (x[m] <  value); break;
                case 3:  ok = (x[m] >= value); break;
                default: ok = (x[m] <= value); break;
            }

            if (!ok) {
                match = false;
                break;
            }
        }

        if (match) {
            return s;
        }
    }
    return -1;
}

__global__
static void stochasticSimulation(const int nx, const int nr,
                                 const int pathsPerTrial, const int blocksPerTrial,
                                 const int nBounds, const int nPoints,
                                 const double tMax, const unsigned long long seed,
                                 const double *d_rates, const int *d_initial,
                                 const int *d_alpha, const int *d_trans,
                                 int *d_paths, int *d_exits, double *d_times) {

    const unsigned int trial = blockIdx.x / blocksPerTrial;
    const unsigned int path = (blockIdx.x % blocksPerTrial) * blockDim.x + threadIdx.x;

    //  fill shared memory

    __shared__ double s_rates[MAX_SSA_REACTIONS];
    __shared__ int s_alpha[MAX_SSA_REACTANTS * MAX_SSA_REACTIONS];
    __shared__ int s_trans[MAX_SSA_REACTANTS * MAX_SSA_REACTIONS];

    for (int i = static_cast<int>(threadIdx.x); i < nr; i += static_cast<int>(blockDim.x)) {
        s_rates[i] = d_rates[trial * nr + i];
    }

    for (int i = static_cast<int>(threadIdx.x); i < nx * nr; i += static_cast<int>(blockDim.x)) {
        s_alpha[i] = d_alpha[trial * nx * nr + i];
        s_trans[i] = d_trans[trial * nx * nr + i];
    }

    __syncthreads();

    if (path >= static_cast<unsigned int>(pathsPerTrial)) {
        return;
    }

    const size_t sample = static_cast<size_t>(trial) * pathsPerTrial + path;

    //  initialize cuRAND

    curandStatePhilox4_32_10 state{};
    curand_init(seed, sample, 0, &state);

    //  Initialize t = t0 and x = x0, set n = 0.

    int x[MAX_SSA_REACTANTS];
    for (int i = 0; i < nx; i++) {
        x[i] = d_initial[trial * nx + i];
    }
    double t = 0.0;

    const bool savePath = (d_paths != nullptr) && (nPoints > 0);
    int *p = savePath ? d_paths + sample * nPoints * nx : nullptr;

    int saved = 0;
    if (savePath) {
        for (int i = 0; i < nx; i++) {
            p[i] = x[i];
        }
    }

    int exitCode = -1;

    while (true) {

        if (nBounds > 0) {
            exitCode = boundary(x, nBounds);
            if (exitCode >= 0) {
                break;
            }
        }

        if (t >= tMax) {
            break;
        }

        //  Compute a_j(x), j = 1, 2, ..., K and a_0(x).

        double a[MAX_SSA_REACTIONS];
        const double a0 = propensity(x, nx, s_rates, nr, s_alpha, a);

        if (a0 <= 0.0) {
            t = tMax;
            break;
        }

        //  Generate r_1, r_2 ~ U([0, 1]).

        const double2 r = curand_uniform2_double(&state);

        //  Use the Golden rule to transform r_1 into tau ~ Exp(a_0(x))

        const double tau = log(1.0 / r.x) / a0;
        const double tNext = t + tau;

        //  Let j be the smallest integer for which sum_{i <= j} a_i > r_2 * a_0

        int j = 0;
        for (double sum_a = a[0]; sum_a < r.y * a0 && j < nr - 1;) {
            j++;
            sum_a += a[j];
        }

        if (savePath) {
            while (saved + 1 < nPoints && (saved + 1) * c_timeStep < tNext) {
                saved++;
                for (int i = 0; i < nx; i++) {
                    p[saved * nx + i] = x[i];
                }
            }
        }

        if (tNext >= tMax) {
            t = tMax;
            break;
        }

        //  Increment t by tau and x by v_j

        for (int i = 0; i < nx; i++) {
            x[i] += s_trans[j * nx + i];
        }
        t = tNext;
    }

    //  record exit state and exit time

    if (d_exits != nullptr) {
        for (int i = 0; i < nx; i++) {
            d_exits[sample * nx + i] = x[i];
        }
    }

    if (d_times != nullptr) {
        d_times[sample] = t;
    }

    //  backfill remainder of saved path with final state

    if (savePath) {
        while (saved + 1 < nPoints) {
            saved++;
            for (int i = 0; i < nx; i++) {
                p[saved * nx + i] = x[i];
            }
        }
    }
}

using namespace Vx;

SimOutput::SimOutput(Model::Params const &params, int numTrials) {
    const size_t species = params.model->species.size();
    _numTrials = numTrials;
    _params = &params;

    if (params.analysis.savePaths.t != 0) {
        int times = static_cast<int>(params.tMax / params.analysis.savePaths.t) + 1;
        cudaMalloc(&d_paths, _params->size * _numTrials * species * times * sizeof(int));
    }

    if (!params.analysis.stopping.bounds.empty()) {
        cudaMalloc(&d_exits, _params->size * _numTrials * species * sizeof(int));
        cudaMalloc(&d_times, _params->size * _numTrials * sizeof(double));
    }
}

SimOutput::~SimOutput() {
    cudaFree(d_paths);
    cudaFree(d_exits);
    cudaFree(d_times);
}

SimOutput::SimOutput(SimOutput &&other) noexcept :
        d_paths(other.d_paths),
        d_exits(other.d_exits),
        d_times(other.d_times) {
    other.d_paths = nullptr;
    other.d_exits = nullptr;
    other.d_times = nullptr;
}

SimOutput &SimOutput::operator=(SimOutput &&other) noexcept {
    if (this != &other) {
        cudaFree(d_paths);
        cudaFree(d_exits);
        cudaFree(d_times);
        d_paths = other.d_paths;
        d_exits = other.d_exits;
        d_times = other.d_times;
        other.d_paths = nullptr;
        other.d_exits = nullptr;
        other.d_times = nullptr;
    }
    return *this;
}

std::vector<int> SimOutput::getPaths() const {
    if (d_paths) {
        const size_t species = _params->model->species.size();
        int times = static_cast<int>(_params->tMax / _params->analysis.savePaths.t) + 1;
        int n = _params->size * _numTrials * species * times;

        auto result = std::vector<int>(n, 0);
        cudaMemcpy(result.data(), d_paths, n * sizeof(int), cudaMemcpyDeviceToHost);
        return result;
    }
    return std::vector<int>(0);
}

std::vector<int> SimOutput::getExits() const {
    if (d_exits) {
        const size_t species = _params->model->species.size();
        int n = _params->size * _numTrials * species;

        auto result = std::vector<int>(n, 0);
        cudaMemcpy(result.data(), d_exits, n * sizeof(int), cudaMemcpyDeviceToHost);
        return result;
    }
    return std::vector<int>(0);
}

std::vector<double> SimOutput::getTimes() const {
    if (d_exits) {
        const size_t species = _params->model->species.size();
        int n = _params->size * _numTrials;

        auto result = std::vector<double>(n, 0);
        cudaMemcpy(result.data(), d_exits, n * sizeof(int), cudaMemcpyDeviceToHost);
        return result;
    }
    return std::vector<double>(0);
}

SimObject::SimObject(Model::Params const &params, SimType type) : _params(params), _type(type) {
    const Model::Model &model = *_params.model;
    const size_t reactions = model.reactions.size();
    const size_t species = model.species.size();
    const size_t nAxes = _params.axes.size();

    std::vector<size_t> axisLen(nAxes);
    for (size_t a = 0; a < nAxes; a++) {
        const Model::Axis &axis = _params.axes[a];
        if (!axis.initial.empty()) {
            axisLen[a] = axis.initial[0].values.size();
        }
        else if (!axis.rates.empty()) {
            axisLen[a] = axis.rates[0].values.size();
        }
        else if (!axis.products.empty()) {
            axisLen[a] = axis.products[0].values.size();
        }
        else {
            axisLen[a] = axis.reactants[0].values.size();
        }
        _numTrials *= axisLen[a];
    }

    std::vector<double> rates(_numTrials * reactions);
    std::vector<int> initial(_numTrials * species, 0);
    std::vector<int> alpha(_numTrials * reactions * species, 0);
    std::vector<int> trans(_numTrials * reactions * species, 0);

    std::vector<double> baseRates(reactions);
    std::vector<int> baseAlpha(reactions * species, 0);
    std::vector<int> baseProd(reactions * species, 0);
    for (size_t j = 0; j < reactions; j++) {
        baseRates[j] = model.reactions[j].rate;
        for (const auto &r: model.reactions[j].reactants) {
            baseAlpha[j * species + (r.species - model.species.data())] = static_cast<int>(r.n);
        }
        for (const auto &p: model.reactions[j].products) {
            baseProd[j * species + (p.species - model.species.data())] = static_cast<int>(p.n);
        }
    }

    std::vector<size_t> idx(nAxes);
    std::vector<int> prod(reactions * species);

    for (size_t trial = 0; trial < _numTrials; trial++) {

        size_t rem = trial;
        for (size_t a = nAxes; a-- > 0;) {
            idx[a] = rem % axisLen[a];
            rem /= axisLen[a];
        }

        double *r = &rates[trial * reactions];
        int *x0 = &initial[trial * species];
        int *al = &alpha[trial * reactions * species];

        std::copy(baseRates.begin(), baseRates.end(), r);
        std::copy(baseAlpha.begin(), baseAlpha.end(), al);
        prod = baseProd;

        for (size_t a = 0; a < nAxes; a++) {
            const Model::Axis &axis = _params.axes[a];

            for (const auto &s: axis.initial) {
                x0[s.species - model.species.data()] = static_cast<int>(s.values[idx[a]]);
            }
            for (const auto &s: axis.rates) {
                r[s.reaction - model.reactions.data()] = s.values[idx[a]];
            }
            for (const auto &s: axis.reactants) {
                al[(s.reaction - model.reactions.data()) * species + (s.species - model.species.data())] = static_cast<int>(s.values[idx[a]]);
            }
            for (const auto &s: axis.products) {
                prod[(s.reaction - model.reactions.data()) * species + (s.species - model.species.data())] = static_cast<int>(s.values[idx[a]]);
            }
        }

        for (size_t i = 0; i < reactions * species; i++) {
            trans[trial * reactions * species + i] = prod[i] - al[i];
        }
    }

    //  create device resources

    cudaMalloc(&d_rates, _numTrials * reactions * sizeof(double));
    cudaMalloc(&d_initial, _numTrials * species * sizeof(int));
    cudaMalloc(&d_alpha, _numTrials * reactions * species * sizeof(int));
    cudaMalloc(&d_trans, _numTrials * reactions * species * sizeof(int));

    //  copy to device resources

    cudaMemcpy(d_rates, rates.data(), rates.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_initial, initial.data(), initial.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_alpha, alpha.data(), alpha.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_trans, trans.data(), trans.size() * sizeof(int), cudaMemcpyHostToDevice);

    //  constant memory

    const double timeStep = _params.analysis.savePaths.t;
    cudaMemcpyToSymbol(c_timeStep, &timeStep, sizeof(double));

    std::vector<int> boundIdxs;
    std::vector<int> boundVals;

    for (const auto &bound: _params.analysis.stopping.bounds) {
        for (const auto &condition: bound.conditions) {
            boundVals.push_back(static_cast<int>(condition.species - model.species.data()));

            if (condition.op == Model::Comparison::Equal) {
                boundVals.push_back(0);
            }
            else if (condition.op == Model::Comparison::GreaterThan) {
                boundVals.push_back(1);
            }
            else if (condition.op == Model::Comparison::LessThan) {
                boundVals.push_back(2);
            }
            else if (condition.op == Model::Comparison::GreaterEqual) {
                boundVals.push_back(3);
            }
            else if (condition.op == Model::Comparison::LessEqual) {
                boundVals.push_back(4);
            }

            boundVals.push_back(static_cast<int>(condition.value));
        }
        boundIdxs.push_back(static_cast<int>(boundVals.size()));
    }

    if (!boundIdxs.empty()) {
        cudaMemcpyToSymbol(c_boundIdxs, boundIdxs.data(), boundIdxs.size() * sizeof(int));
        cudaMemcpyToSymbol(c_boundVals, boundVals.data(), boundVals.size() * sizeof(int));
    }
}

SimObject::~SimObject() {
    cudaFree(d_rates);
    cudaFree(d_initial);
    cudaFree(d_alpha);
    cudaFree(d_trans);

    delete simOutput;
}

SimObject::SimObject(SimObject &&other) noexcept :
        _params(other._params),
        _type(other._type),
        _numTrials(other._numTrials),
        d_rates(other.d_rates),
        d_initial(other.d_initial),
        d_alpha(other.d_alpha),
        d_trans(other.d_trans),
        simOutput(other.simOutput) {
    other.d_rates = nullptr;
    other.d_initial = nullptr;
    other.d_alpha = nullptr;
    other.d_trans = nullptr;
    other.simOutput = nullptr;
}

SimObject &SimObject::operator=(SimObject &&other) noexcept {
    return *this;
}

void SimObject::runSimulation() {
    const Model::Model &model = *_params.model;

    const int nx = static_cast<int>(model.species.size());
    const int nr = static_cast<int>(model.reactions.size());

    const int pathsPerTrial = static_cast<int>(_params.size);
    const int blocksPerTrial = (pathsPerTrial + SSA_BLOCK_SIZE - 1) / SSA_BLOCK_SIZE;

    const int nBounds = static_cast<int>(_params.analysis.stopping.bounds.size());
    const int nPoints = _params.analysis.savePaths.t != 0 ? static_cast<int>(_params.tMax / _params.analysis.savePaths.t) + 1 : 0;

    delete simOutput;
    simOutput = new SimOutput(_params, _numTrials);

    switch (_type) {
        case SimType::StochasticSimulation:
            stochasticSimulation<<<blocksPerTrial * _numTrials, SSA_BLOCK_SIZE>>>(
                nx, nr, pathsPerTrial, blocksPerTrial, nBounds, nPoints,
                _params.tMax, _params.seed, d_rates, d_initial, d_alpha, d_trans,
                simOutput->d_paths, simOutput->d_exits, simOutput->d_times);
            break;
        default:
            break;
    }
}

SimOutput* SimObject::deviceSynchronize() {
    cudaDeviceSynchronize();
    SimOutput* result = simOutput;
    simOutput = nullptr;
    return result;
}
