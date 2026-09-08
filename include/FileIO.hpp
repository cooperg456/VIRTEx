/******************************************************************************
 *  Copyright (c) 2026 Cooper Gray
 *
 *  This application is free software, meaning you can redistribute it and/or
 *  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
 *  for details. You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 ******************************************************************************/

#pragma once

#include "SimBuilder.hpp"

#include <filesystem>

namespace FileIO {

/******************************************************************************
 *  JSON model parser
 ******************************************************************************/

SimBuilder::Model loadModel_JSON(const std::filesystem::path& filePath);

}   //  FileIO

