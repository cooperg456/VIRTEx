////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#pragma once

#include <string>
#include <vector>

namespace Vx::Model {
    struct Species {
        std::string id{};
        std::string symbol{};
        std::string desc{};

        friend std::ostream &operator<<(std::ostream &os, const Species &obj);
    };

    struct NSpecies {
        unsigned int n = 0;

        const Species *species = nullptr;

        friend std::ostream &operator<<(std::ostream &os, const NSpecies &obj);
    };

    struct Reaction {
        double rate = 0;

        std::string id{};
        std::string symbol{};
        std::string desc{};

        std::vector<NSpecies> reactants{};
        std::vector<NSpecies> products{};

        friend std::ostream &operator<<(std::ostream &os, const Reaction &obj);
    };

    struct Model {
        std::string desc{};

        std::vector<Species> species{};
        std::vector<Reaction> reactions{};

        friend std::ostream &operator<<(std::ostream &os, const Model &obj);
    };

    enum class Comparison {
        Equal,
        LessThan,
        GreaterThan,
        LessEqual,
        GreaterEqual
    };

    struct Condition {
        const Species *species = nullptr;
        Comparison op = Comparison::Equal;
        unsigned int value = 0;

        friend std::ostream &operator<<(std::ostream &os, const Condition &obj);
    };

    struct Bound {
        std::string id{};

        std::vector<Condition> conditions{};

        friend std::ostream &operator<<(std::ostream &os, const Bound &obj);
    };

    struct Stopping {
        std::vector<Bound> bounds{};

        friend std::ostream &operator<<(std::ostream &os, const Stopping &obj);
    };

    struct SavePaths {
        double t = 0;

        friend std::ostream &operator<<(std::ostream &os, const SavePaths &obj);
    };

    //  "it's... it's... Analysis!!!"
    struct Analysis {
        SavePaths savePaths{};
        Stopping stopping{};

        friend std::ostream &operator<<(std::ostream &os, const Analysis &obj);
    };

    struct USweep {
        std::vector<unsigned int> values{};

        const Species *species = nullptr;
        const Reaction *reaction = nullptr;

        friend std::ostream &operator<<(std::ostream &os, const USweep &obj);
    };

    struct DSweep {
        std::vector<double> values{};

        const Species *species = nullptr;
        const Reaction *reaction = nullptr;

        friend std::ostream &operator<<(std::ostream &os, const DSweep &obj);
    };

    struct Axis {
        std::vector<USweep> initial{};
        std::vector<USweep> reactants{};
        std::vector<USweep> products{};
        std::vector<DSweep> rates{};

        friend std::ostream &operator<<(std::ostream &os, const Axis &obj);
    };

    struct Params {
        double tMax = 0;
        unsigned int size = 0;
        unsigned long long seed = 0;

        std::string desc{};

        Analysis analysis{};

        std::vector<Axis> axes{};

        const Model *model = nullptr;

        friend std::ostream &operator<<(std::ostream &os, const Params &obj);
    };
}
