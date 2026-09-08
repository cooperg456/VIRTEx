/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#include "FileIO.hpp"

#include "nlohmann/json.hpp"

#include <fstream>

namespace FileIO {

/******************************************************************************
 *  JSON model parser
 ******************************************************************************/

SimBuilder::Model loadModel_JSON(const  std::filesystem::path& filePath) {
    std::ifstream f(filePath);

    if (!f.is_open()) {
        std::string msg = "Failed to open file: " + filePath.string();
        f.close();
        throw std::runtime_error(msg);
    }

    nlohmann::json network = nlohmann::json::parse(f);
    f.close();

    SimBuilder::Model m{};
    m.reactants = network["reactants"].get<std::vector<std::string>>();
    m.ics = network["conditions"].get<std::vector<int>>();

    size_t n_reactants = network["reactants"].size();
    size_t n_reactions = network["reactions"].size();

    m.reactions.resize(n_reactions);
    m.rates.resize(n_reactions);

    m.rcs = std::vector<int>(n_reactants * n_reactions, 0);
    m.tcs = std::vector<int>(n_reactants * n_reactions, 0);;

    for (size_t i = 0; i < network["reactions"].size(); i++) {
        m.reactions[i] = network["reactions"][i]["name"];
        m.rates[i] = network["reactions"][i]["rate"];

        for (auto& reactant : network["reactions"][i]["reactants"].items()) {
            size_t j = SimBuilder::getIdx(reactant.key(), m.reactants);
            m.rcs[i * n_reactants + j] += static_cast<int>(reactant.value());
            m.tcs[i * n_reactants + j] -= static_cast<int>(reactant.value());
        }

        for (auto& product : network["reactions"][i]["products"].items()) {
            size_t j = SimBuilder::getIdx(product.key(), m.reactants);
            m.tcs[i * n_reactants + j] += static_cast<int>(product.value());
        }
    }

    return m;
}

}   //  FileIO