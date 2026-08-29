/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#pragma once

#include <vector>
#include <string>
#include <filesystem>

/******************************************************************************
 *  ReactionNetwork class
 ******************************************************************************/

class ReactionNetwork {
public:
    ReactionNetwork() = default;
    explicit ReactionNetwork(const std::filesystem::path& jsonFile);

    /**************************************************************************
     *  Reaction network parameters
     **************************************************************************/
    
    std::vector<std::string> reactants{};
    std::vector<std::string> reactions{};

    std::vector<double> reactionRates{};
    std::vector<int> initialConditions{};
    std::vector<int> reactantCoefficients{};
    std::vector<int> transitionCoefficients{};

    /**************************************************************************
     *  Helper functions
     **************************************************************************/

    size_t getReactantIdx(const std::string& reactant);
    size_t getReactionIdx(const std::string& reaction);
};
