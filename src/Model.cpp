////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#include <virtex/Model.hpp>

#include <iostream>

using namespace Vx::Model;

namespace {
    const char *toString(Comparison op) {
        switch (op) {
            case Comparison::Equal: return "=";
            case Comparison::LessThan: return "<";
            case Comparison::GreaterThan: return ">";
            case Comparison::LessEqual: return "<=";
            case Comparison::GreaterEqual: return ">=";
        }
        return "";
    }
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Species &obj) {
    return os << obj.id << " (" << obj.symbol << ") - " << obj.desc;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const NSpecies &obj) {
    return os << (obj.n > 1 ? std::to_string(obj.n) + " " : "") << obj.species->symbol;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Reaction &obj) {
    os << obj.id << " (" << obj.symbol << ") : ";
    for (size_t i = 0; i + 1 < obj.reactants.size(); i++) {
        os << obj.reactants[i] << " + ";
    }
    if (!obj.reactants.empty()) {
        os << obj.reactants[obj.reactants.size() - 1];
    } else {
        os << "0";
    }
    os << " -> (rate " << obj.rate << ") -> ";
    for (size_t i = 0; i + 1 < obj.products.size(); i++) {
        os << obj.products[i] << " + ";
    }
    if (!obj.products.empty()) {
        os << obj.products[obj.products.size() - 1];
    } else {
        os << "0";
    }
    if (!obj.desc.empty()) {
        os << " - " << obj.desc;
    }
    return os;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Model &obj) {
    os << "Model: " << obj.desc << "\n";
    os << "Species:\n";
    for (auto &_species: obj.species) {
        os << "\t" << _species << "\n";
    }
    os << "Reactions:\n";
    for (auto &_reaction: obj.reactions) {
        os << "\t" << _reaction << "\n";
    }
    return os;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Condition &obj) {
    if (obj.species) {
        os << obj.species->symbol;
    } else {
        os << "?";
    }
    return os << " " << toString(obj.op) << " " << obj.value;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Bound &obj) {
    os << obj.id << ": ";
    for (size_t i = 0; i + 1 < obj.conditions.size(); i++) {
        os << obj.conditions[i] << ", ";
    }
    if (!obj.conditions.empty()) {
        os << obj.conditions[obj.conditions.size() - 1];
    }
    return os;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Stopping &obj) {
    os << "stopping:\n";
    for (auto &_bound: obj.bounds) {
        os << "\t\t" << _bound << "\n";
    }
    return os;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const SavePaths &obj) {
    return os << "savepaths: t = " << obj.t << "\n";
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Analysis &obj) {
    if (obj.savePaths.t != 0) {
        os << "\t" << obj.savePaths;
    }
    if (!obj.stopping.bounds.empty()) {
        os << "\t" << obj.stopping;
    }
    return os;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const USweep &obj) {
    if (obj.species) {
        os << "species " << obj.species->symbol;
    } else if (obj.reaction) {
        os << "reaction " << obj.reaction->symbol;
    }
    os << " : [";
    for (size_t i = 0; i + 1 < obj.values.size(); i++) {
        os << obj.values[i] << ", ";
    }
    if (!obj.values.empty()) {
        os << obj.values[obj.values.size() - 1];
    }
    return os << "]";
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const DSweep &obj) {
    if (obj.species) {
        os << "species " << obj.species->symbol;
    } else if (obj.reaction) {
        os << "reaction " << obj.reaction->symbol;
    }
    os << " : [";
    for (size_t i = 0; i + 1 < obj.values.size(); i++) {
        os << obj.values[i] << ", ";
    }
    if (!obj.values.empty()) {
        os << obj.values[obj.values.size() - 1];
    }
    return os << "]";
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Axis &obj) {
    if (!obj.initial.empty()) {
        os << "\tinitial:\n";
        for (auto &_u: obj.initial) {
            os << "\t\t" << _u << "\n";
        }
    }
    if (!obj.reactants.empty()) {
        os << "\treactants:\n";
        for (auto &_u: obj.reactants) {
            os << "\t\t" << _u << "\n";
        }
    }
    if (!obj.products.empty()) {
        os << "\tproducts:\n";
        for (auto &_u: obj.products) {
            os << "\t\t" << _u << "\n";
        }
    }
    if (!obj.rates.empty()) {
        os << "\trates:\n";
        for (auto &_d: obj.rates) {
            os << "\t\t" << _d << "\n";
        }
    }
    return os;
}

std::ostream &Vx::Model::operator<<(std::ostream &os, const Params &obj) {
    os << "Params: " << obj.desc << "\n";
    os << "tMax: " << obj.tMax << "\n";
    os << "size: " << obj.size << "\n";
    os << "seed: " << obj.seed << "\n";
    os << "Analysis:\n" << obj.analysis;
    os << "Axes:\n";
    for (size_t i = 0; i < obj.axes.size(); i++) {
        os << "\tAxis (" << (i + 1) << "):\n" << obj.axes[i] << "\n";
    }
    return os;
}
