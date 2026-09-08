/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "SimBuilder.hpp"

#include <algorithm>
#include <stdexcept>

namespace SimBuilder {

/******************************************************************************
 *  Sweep builder
 ******************************************************************************/

static size_t getIdx(const std::string& str, std::vector<std::string> vec) {
    const auto it = std::ranges::find(vec, str);
    if (it == vec.end()) {
        throw std::runtime_error("Unknown member: " + str);
    }
    return std::distance(vec.begin(), it);
}

std::vector<Model> BuildSweep(const Model& model, const Modifiers& mods) {
    std::vector<int> sweepSizes;
    size_t reactantsSize = model.reactants.size();
    sweepSizes.reserve(mods.rates.size());
    for (auto & rate : mods.rates) {
        sweepSizes.push_back(static_cast<int>(std::get<1>(rate).size()));
    }

    for (auto & react : mods.reactants) {
        int sweepSize = static_cast<int>(std::get<1>(react).size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error("Length of required reactants must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    for (auto & prod : mods.products) {
        int sweepSize = static_cast<int>(std::get<1>(prod).size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error("Length of reaction products must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    if (!mods.ics.empty()) {
        int sweepSize = static_cast<int>(mods.ics.size());
        if (sweepSize % reactantsSize != 0) {
            throw std::runtime_error("Length of initial conditions must be a multiple of the number of reactants");
        }
        sweepSizes.push_back(sweepSize / static_cast<int>(reactantsSize));
    }

    if (sweepSizes.size() > 1 && std::ranges::adjacent_find(sweepSizes, std::not_equal_to<>()) != sweepSizes.end()) {
        throw std::runtime_error("Length of each sweep must be equal");
    }

    int sweepSize = 1;
    if (!sweepSizes.empty()) {
        sweepSize = sweepSizes.front();
    }

    std::vector<Model> sweep(sweepSize);
    
    for (size_t i = 0; i < sweepSize; i++) {
        sweep[i] = model;

        for (auto & j : mods.rates) {
            size_t reaction = getIdx(std::get<0>(j), model.reactions);
            double rate = std::get<1>(j)[i];

            sweep[i].rates[reaction] = rate;
        }

        //  reacts loop must come before prods
        for (auto & j : mods.reactants) {
            size_t reaction = getIdx(std::get<0>(j), model.reactions);
            std::tuple_element_t<1, std::tuple<std::string, std::vector<int>>> react = std::get<1>(j);

            for (size_t k = 0; k < reactantsSize; k++) {
                size_t l = reaction * reactantsSize + k;
                size_t reactIdx = i * reactantsSize + k;

                int beta = model.rcs[l] + model.tcs[l];

                sweep[i].rcs[l] = react[reactIdx];
                sweep[i].tcs[l] = beta - react[reactIdx];
            }
        }

        for (auto & j : mods.products) {
            size_t reaction = getIdx(std::get<0>(j), model.reactions);
            std::tuple_element_t<1, std::tuple<std::string, std::vector<int>>> prod = std::get<1>(j);

            for (size_t k = 0; k < reactantsSize; k++) {
                size_t l = reaction * reactantsSize + k;
                size_t prodIdx = i * reactantsSize + k;

                sweep[i].tcs[l] = prod[prodIdx] - sweep[i].rcs[l];
            }
        }

        if (!mods.ics.empty()) {
            for (size_t j = 0; j < reactantsSize; j++) {
                sweep[i].ics[j] = mods.ics[i * reactantsSize + j];
            }
        }
    }
    
    return sweep;
}

}   //  SimBuilder