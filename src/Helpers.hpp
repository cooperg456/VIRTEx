////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#pragma once

#include <filesystem>

#include <virtex/Model.hpp>

namespace Vx::Helpers {
    std::vector<std::string> tokenize(std::string str, const std::string &sep);

    std::vector<double> linspace(double start, double stop, size_t n);

    std::vector<double> geomspace(double start, double stop, size_t n);

    std::string loadFile(const std::filesystem::path &filename);

    Vx::Model::Model parseModel(const std::string &xmlString);

    Vx::Model::Params parseParams(const std::string &xmlString, const Vx::Model::Model &model);
}
