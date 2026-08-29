/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "ReactionNetwork.hpp"

#include "nlohmann/json.hpp"

#include <fstream>

/******************************************************************************
 *  ReactionNetwork functions
 ******************************************************************************/

ReactionNetwork::ReactionNetwork(const std::filesystem::path& jsonFile) {
    std::ifstream f(jsonFile);

    if (!f.is_open()) {
        std::string msg = "Failed to open file: " + jsonFile.string();
        if (f.bad()) msg += " (fatal stream error)";
        if (f.fail()) msg += std::string(" - ") + std::strerror(errno);

        f.close();
        throw std::runtime_error(msg);
    }

    nlohmann::json network = nlohmann::json::parse(f);
    f.close();

    size_t n_reactants = network["reactants"].size();
    size_t n_reactions = network["reactions"].size();

    reactants = network["reactants"];
    initialConditions = network["conditions"].get<std::vector<int>>();

    reactions.resize(n_reactions);
    reactionRates.resize(n_reactions);
    reactantCoefficients = std::vector<int>(n_reactants * n_reactions, 0);
    transitionCoefficients = std::vector<int>(n_reactants * n_reactions, 0);

    for (size_t i = 0; i < network["reactions"].size(); i++) {
        reactions[i] = network["reactions"][i]["name"];
        reactionRates[i] = network["reactions"][i]["rate"];

        for (auto& reactant : network["reactions"][i]["reactants"].items()) {
            size_t j = getReactantIdx(reactant.key());
            reactantCoefficients[i * n_reactants + j] += static_cast<int>(reactant.value());
            transitionCoefficients[i * n_reactants + j] -= static_cast<int>(reactant.value());
        }

        for (auto& product : network["reactions"][i]["products"].items()) {
            size_t j = getReactantIdx(product.key());
            transitionCoefficients[i * n_reactants + j] += static_cast<int>(product.value());
        }
    }
}

/******************************************************************************
 *  Helper functions
 ******************************************************************************/

size_t ReactionNetwork::getReactantIdx(const std::string& reactant) {
    auto it = std::find(reactants.begin(), reactants.end(), reactant);
    if (it == reactants.end()) {
        throw std::runtime_error("Unknown reactant: " + reactant);
    }
    return std::distance(reactants.begin(), it);
}

size_t ReactionNetwork::getReactionIdx(const std::string& reaction) {
    auto it = std::find(reactions.begin(), reactions.end(), reaction);
    if (it == reactions.end()) {
        throw std::runtime_error("Unknown reaction: " + reaction);
    }
    return std::distance(reactions.begin(), it);
}
