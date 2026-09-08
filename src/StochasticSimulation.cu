/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "StochasticSimulation.hpp"

#include "cuda_runtime.h"
#include "curand_kernel.h"

#include <fstream>
#include <iostream>

/******************************************************************************
 *  CUDA constant memory
 ******************************************************************************/

__constant__ static int c_boundIdxs[MAX_SSA_BOUNDARIES];
__constant__ static int c_boundVals[MAX_SSA_REACTANTS * 3];

/******************************************************************************
 *  Stochastic simulation algorithm
 ******************************************************************************/

__global__
static void cuSSA(const int reactants, const int reactions, const int warps, const int savedPaths, const int boundSets,
                  const double tGrid, const double tMax, const unsigned long long seed,
                  const int* d_conds, const double* d_rates, const int* d_alpha, const int* d_trans,
                  double* d_times, int* d_paths, double* d_exitTimes, int* d_exits) {

    const unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
    const unsigned int sim = blockIdx.x / warps;
    const unsigned int path = idx % (warps * SSA_BLOCK_SIZE);

    //  fill shared memory

    __shared__ double s_rates[MAX_SSA_REACTIONS];
    __shared__ int s_alpha[MAX_SSA_REACTANTS * MAX_SSA_REACTIONS];
    __shared__ int s_trans[MAX_SSA_REACTANTS * MAX_SSA_REACTIONS];

    for (unsigned int i = threadIdx.x; i < reactions; i += blockDim.x) {
        s_rates[i] = d_rates[blockIdx.x * reactions + i];
    }

    for (unsigned int i = threadIdx.x; i < reactants * reactions; i += blockDim.x) {
        s_alpha[i] = d_alpha[blockIdx.x * reactants * reactions + i];
        s_trans[i] = d_trans[blockIdx.x * reactants * reactions + i];
    }

    __syncthreads();

    //  initialize cuRAND

    curandStatePhilox4_32_10 state{};
    curand_init(seed, idx, 0, &state);

    //  Initialize t = t0 and x = x0, set n = 0.

    int x[MAX_SSA_REACTANTS];
    for (int i = 0; i < reactants; i++) {
        x[i] = d_conds[blockIdx.x * reactants + i];
    }
    double t = 0;

    //  record first n paths

    int tSaved = 0;
    int pointsPerPath = 0;
    if (path < savedPaths) {
        pointsPerPath = static_cast<int>(tMax / tGrid + 1);
        for (int i = 0; i < reactants; i++) {
            d_paths[(sim * savedPaths * pointsPerPath * reactants) + (path * pointsPerPath * reactants) + i] = x[i];
        }
        d_times[(sim * savedPaths * pointsPerPath) + (path * pointsPerPath)] = t;
    }

    int exitCode = -1;

    //  while (t < tMax);

    do {

            //  Compute a_j(x), j = 1, 2, ..., K and a_0(x).

            double a[MAX_SSA_REACTIONS];
            double a0 = 0;
            for (int j = 0; j < reactions; j++) {
                a[j] = s_rates[j];
                for (int m = 0; m < reactants; m++) {
                    if (const int alpha_jm = s_alpha[j * reactants + m]; alpha_jm > 0) {
                        if (x[m] < alpha_jm) {
                            a[j] = 0.0;
                            break;
                        }
                        for (int k = x[m] - alpha_jm + 1; k <= x[m]; k++) {
                            a[j] *= k;
                        }
                    }
                }
                a0 += a[j];
            }

            if (a0 != 0) {
                //  Generate r_1, r_2 ∼ U([0, 1]).

                double2 r = curand_uniform2_double(&state);

                //  Use the Golden rule to transform r_1 into τ ∼ Exp(a_0(x))

                double tau = log(1 / r.x) / a0;

                //  Let j be the smallest integer for which

                int j = 0;
                for (double sum_a = a[0]; sum_a < r.y * a0 && j < reactions - 1;) {
                    j++;
                    sum_a += a[j];
                }

                //  Increment t by τ and x by v_j

                for (int i = 0; i < reactants; i++) {
                    x[i] += s_trans[reactants * j + i];
                }
                t += tau;
            }
            else {
                t = tMax;
            }

            //  record first n paths

            if (path < savedPaths) {
                while (tSaved * tGrid < t && tSaved < pointsPerPath - 1) {
                    tSaved++;
                    for (int i = 0; i < reactants; i++) {
                        d_paths[(sim * savedPaths * pointsPerPath * reactants) + (path * pointsPerPath * reactants) + (tSaved * reactants) + i] = x[i];
                    }
                    d_times[(sim * savedPaths * pointsPerPath) + (path * pointsPerPath) + tSaved] = tSaved * tGrid;
                }
            }

            //  check boundary conditions

            for (int s = 0; s < boundSets; s++) {
                const int start = (s == 0) ? 0 : c_boundIdxs[s - 1];
                const int end = c_boundIdxs[s];
                bool match = true;
                for (int k = start; k < end; k += 3) {
                    const int reactantIdx = c_boundVals[k];
                    const int op = c_boundVals[k + 1];
                    const int value = c_boundVals[k + 2];
                    bool ok;
                    if (op == 0) ok = (x[reactantIdx] == value);
                    else if (op == 1) ok = (x[reactantIdx] > value);
                    else ok = (x[reactantIdx] < value);
                    if (!ok) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    exitCode = s;
                    break;
                }
            }
            if (exitCode >= 0) break;

    } while (t < tMax);

    //  record exit time and code

    d_exitTimes[idx] = t;
    d_exits[idx] = exitCode;

    //  backfill remainder of saved path with final state

    if (path < savedPaths) {
        while (tSaved < pointsPerPath - 1) {
            tSaved++;
            for (int i = 0; i < reactants; i++) {
                d_paths[(sim * savedPaths * pointsPerPath * reactants) + (path * pointsPerPath * reactants) + (tSaved * reactants) + i] = x[i];
            }
            d_times[(sim * savedPaths * pointsPerPath) + (path * pointsPerPath) + tSaved] = tSaved * tGrid;
        }
    }
}

void SSA(const SSASimInfo& simInfo, const std::vector<SSASysInfo> &sysInfos) {
    int numBlocks = simInfo.warps * static_cast<int>(sysInfos.size());

    //  ReactionNetwork info

    double *d_reactionRates;
    int *d_initialConditions;
    int *d_reactantCoefficients;
    int *d_transitionCoefficients;

    double *d_timePoints = nullptr;
    int *d_samplePaths = nullptr;

    double *d_exitTimes;
    int *d_exits;

    const int n_reactionRates = static_cast<int>(simInfo.base.reactions.size());
    const int n_initialConditions = static_cast<int>(simInfo.base.initialConditions.size());
    const int n_reactantCoefficients = static_cast<int>(simInfo.base.reactantCoefficients.size());
    const int n_transitionCoefficients = static_cast<int>(simInfo.base.transitionCoefficients.size());

    const int n_exits = static_cast<int>(sysInfos.size()) * simInfo.warps * SSA_BLOCK_SIZE;
    const int boundSets = static_cast<int>(simInfo.boundIdxs.size());

    int n_timePoints = 0;
    int n_samplePaths = 0;
    int pointsPerPath = 0;
    if (simInfo.savedPaths) {
        pointsPerPath = static_cast<int>(simInfo.tMax / simInfo.tGrid + 1);
        n_timePoints = simInfo.savedPaths * pointsPerPath;
        n_samplePaths = n_timePoints * n_initialConditions;
    }

    //  malloc buffers

    cudaMalloc(&d_reactionRates, numBlocks * n_reactionRates * sizeof(double));
    cudaMalloc(&d_reactantCoefficients, numBlocks * n_reactantCoefficients * sizeof(int));
    cudaMalloc(&d_transitionCoefficients, numBlocks * n_transitionCoefficients * sizeof(int));

    if (simInfo.savedPaths) {
        cudaMallocManaged(&d_timePoints, n_timePoints * sysInfos.size() * sizeof(double));
        cudaMallocManaged(&d_samplePaths, n_samplePaths * sysInfos.size() * sizeof(int));
    }
    cudaMalloc(&d_initialConditions, numBlocks * n_initialConditions * sizeof(int));

    cudaMallocManaged(&d_exitTimes, n_exits * sizeof(double));
    cudaMallocManaged(&d_exits, n_exits * sizeof(int));

    //  memcpy to gpu

    for (int sim = 0; sim < static_cast<int>(sysInfos.size()); sim++) {
        for (int w = 0; w < simInfo.warps; w++) {
            const int block = sim * simInfo.warps + w;
            cudaMemcpy(d_reactionRates + n_reactionRates * block, sysInfos[sim].reactionRates.data(), n_reactionRates * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(d_reactantCoefficients + n_reactantCoefficients * block, sysInfos[sim].reactantCoefficients.data(), n_reactantCoefficients * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_transitionCoefficients + n_transitionCoefficients * block, sysInfos[sim].transitionCoefficients.data(), n_transitionCoefficients * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_initialConditions + n_initialConditions * block, sysInfos[sim].initialConditions.data(), n_initialConditions * sizeof(int), cudaMemcpyHostToDevice);
        }
    }

    //  copy boundary info to constant memory

    cudaMemcpyToSymbol(c_boundIdxs, simInfo.boundIdxs.data(),
                        simInfo.boundIdxs.size() * sizeof(int));
    cudaMemcpyToSymbol(c_boundVals, simInfo.boundVals.data(),
                        simInfo.boundVals.size() * sizeof(int));

    //  prefetches

    int device;
    cudaGetDevice(&device);

    cudaMemLocation memlocDev{};
    memlocDev.id = device;
    memlocDev.type = cudaMemLocationTypeDevice;

    if (simInfo.savedPaths) {
        cudaMemPrefetchAsync(d_timePoints, n_timePoints * sysInfos.size() * sizeof(double), memlocDev, 0);
        cudaMemPrefetchAsync(d_samplePaths, n_samplePaths * sysInfos.size() * sizeof(int), memlocDev, 0);
    }

    cudaMemPrefetchAsync(d_exitTimes, n_exits * sizeof(double), memlocDev, 0);
    cudaMemPrefetchAsync(d_exits, n_exits * sizeof(int), memlocDev, 0);

    //  kernel launch

    cuSSA<<<numBlocks, SSA_BLOCK_SIZE>>>(n_initialConditions, n_reactionRates,
                                         simInfo.warps, simInfo.savedPaths, boundSets,
                                         simInfo.tGrid, simInfo.tMax, simInfo.seed,
                                         d_initialConditions, d_reactionRates,
                                         d_reactantCoefficients, d_transitionCoefficients,
                                         d_timePoints, d_samplePaths,
                                         d_exitTimes, d_exits);

    cudaDeviceSynchronize();

    //  free unneeded allocations

    cudaFree(d_reactionRates);
    cudaFree(d_reactantCoefficients);
    cudaFree(d_transitionCoefficients);
    cudaFree(d_initialConditions);

    //  prefetch results back to host

    cudaMemLocation memlocHost{};
    memlocHost.id = 0;
    memlocHost.type = cudaMemLocationTypeHost;

    cudaMemPrefetchAsync(d_exitTimes, n_exits * sizeof(double), memlocHost, 0);
    cudaMemPrefetchAsync(d_exits, n_exits * sizeof(int), memlocHost, 0);

    if (simInfo.savedPaths) {
        cudaMemPrefetchAsync(d_timePoints, n_timePoints * sysInfos.size() * sizeof(double), memlocHost, 0);
        cudaMemPrefetchAsync(d_samplePaths, n_samplePaths * sysInfos.size() * sizeof(int), memlocHost, 0);
    }

    cudaDeviceSynchronize();

    //  get data

    std::filesystem::create_directories(simInfo.outputDir);

    const int pathsPerSim = simInfo.warps * SSA_BLOCK_SIZE;

    for (size_t sim = 0; sim < sysInfos.size(); sim++) {
        std::filesystem::path exitFile = simInfo.outputDir;
        exitFile.append("ssa_exits_" + std::to_string(sim) + ".csv");

        std::ofstream efile(exitFile);
        if (!efile.is_open()) {
            std::cerr << "Failed to open csv file for writing: " << exitFile << "\n";
            continue;
        }

        efile << "path,exit_time,exit_code\n";

        for (int path = 0; path < pathsPerSim; path++) {
            const size_t exitIdx = (sim * pathsPerSim) + path;
            efile << path << "," << d_exitTimes[exitIdx] << "," << d_exits[exitIdx] << "\n";
        }

        efile.close();
    }

    if (simInfo.savedPaths) {
        for (size_t sim = 0; sim < sysInfos.size(); sim++) {
            std::filesystem::path outFile = simInfo.outputDir;
            outFile.append("ssa_trajectories_" + std::to_string(sim) + ".csv");

            std::ofstream file(outFile);
            if (!file.is_open()) {
                std::cerr << "Failed to open csv file for writing: " << outFile << "\n";
                continue;
            }

            file << "path,time";
            for (int r = 0; r < n_initialConditions; r++) {
                file << ",reactant_" << r;
            }
            file << "\n";

            for (int path = 0; path < simInfo.savedPaths; path++) {
                for (int step = 0; step < pointsPerPath; step++) {
                    size_t timeIdx = ((size_t)sim * simInfo.savedPaths * pointsPerPath)
                                    + (path * pointsPerPath) + step;
                    size_t pathBase = ((size_t)sim * simInfo.savedPaths * pointsPerPath * n_initialConditions)
                                    + (path * pointsPerPath * n_initialConditions)
                                    + (step * n_initialConditions);

                    file << path << "," << d_timePoints[timeIdx];
                    for (int r = 0; r < n_initialConditions; r++) {
                        file << "," << d_samplePaths[pathBase + r];
                    }
                    file << "\n";
                }
            }

            file.close();
        }

        cudaFree(d_timePoints);
        cudaFree(d_samplePaths);
    }

    cudaFree(d_exitTimes);
    cudaFree(d_exits);
}

/*
for (auto& boundSet : args.bounds.bounds) {
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
