////////////////////////////////////////////////////////////////////////////////
//  VIRTEx — Viral Infection, Recovery, and Transmission Explorer
//  https://github.com/cooperg456/VIRTEx
//
//  A set of parallel tools for analyzing Continuous-Time Markov Chains (CTMCs)
//  derived from the Chemical Master Equation (CME).
//
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#include <argparse/argparse.hpp>
#include <xml/xml.h>

#include <iostream>
#include <fstream>
#include <random>

#include <virtex/Model.hpp>
#include <virtex/SimObject.hpp>

namespace {
    std::vector<std::string> tokenize(std::string str, const std::string &sep) {
        std::vector<std::string> tokens;

        size_t pos = 0;
        while ((pos = str.find(sep)) != std::string::npos) {
            std::string token = str.substr(0, pos);
            tokens.push_back(token);
            str.erase(0, pos + sep.length());
        }
        tokens.push_back(str);

        return tokens;
    }

    std::vector<double> linspace(double start, double stop, size_t n) {
        double step = (stop - start) / static_cast<double>(n - 1);

        std::vector<double> result(n);
        for (size_t i = 0; i < n - 1; i++) {
            result[i] = (start + static_cast<double>(i) * step);
        }
        result[n - 1] = stop;

        return result;
    }

    std::vector<double> geomspace(double start, double stop, size_t n) {
        double step = (std::log(stop) - std::log(start)) / static_cast<double>(n - 1);

        std::vector<double> result(n);
        for (size_t i = 0; i < n - 1; i++) {
            result[i] = (start * std::exp(static_cast<double>(i) * step));
        }
        result[n - 1] = stop;

        return result;
    }

    std::string loadFile(const std::filesystem::path &filename) {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Failed to open " + filename.string());
        }
        std::stringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }

    Vx::Model::Model parseModel(const std::string &xmlString) {
        XMLNode *root = xml_parse_string(xmlString.c_str());
        XMLNode *model = xml_node_child_at(root, 0);

        Vx::Model::Model _model{};

        _model.desc = xml_node_attr(model, "desc");
        for (size_t i = 0; i < model->children->len; i++) {
            XMLNode *component = xml_node_child_at(model, i);

            if (!strcmp(component->tag, "species")) {
                Vx::Model::Species species;

                species.id = xml_node_attr(component, "id");
                species.symbol = xml_node_attr(component, "symbol");
                species.desc = xml_node_attr(component, "desc");

                _model.species.push_back(species);
            }

            if (!strcmp(component->tag, "reaction")) {
                Vx::Model::Reaction reaction;

                reaction.id = xml_node_attr(component, "id");
                reaction.symbol = xml_node_attr(component, "symbol");
                reaction.desc = xml_node_attr(component, "desc");
                reaction.rate = std::stod(xml_node_attr(component, "rate"));

                for (size_t j = 0; j < component->children->len; j++) {
                    XMLNode *species = xml_node_child_at(component, j);

                    Vx::Model::NSpecies nSpecies;

                    nSpecies.n = std::stoul(xml_node_attr(species, "n"));

                    std::string id = xml_node_attr(species, "species");
                    for (const auto &_species: _model.species) {
                        if (_species.id == id) {
                            nSpecies.species = &_species;
                            break;
                        }
                    }

                    if (!strcmp(species->tag, "reactant")) {
                        reaction.reactants.push_back(nSpecies);
                    }
                    if (!strcmp(species->tag, "product")) {
                        reaction.products.push_back(nSpecies);
                    }
                }
                _model.reactions.push_back(reaction);
            }
        }
        xml_node_free(root);
        return _model;
    }

    Vx::Model::Params parseParams(const std::string &xmlString, const Vx::Model::Model &model) {
        XMLNode *root = xml_parse_string(xmlString.c_str());
        XMLNode *params = xml_node_child_at(root, 0);

        Vx::Model::Params _params{};

        if (xml_node_attr(params, "desc")) {
            _params.desc = xml_node_attr(params, "desc");
        }
        if (xml_node_attr(params, "batchsize")) {
            _params.size = std::stoul(xml_node_attr(params, "batchsize"));
        }
        if (xml_node_attr(params, "T")) {
            _params.tMax = std::stod(xml_node_attr(params, "T"));
        }
        if (xml_node_attr(params, "seed")) {
            _params.seed = std::stoul(xml_node_attr(params, "seed"));
        }

        _params.model = &model;
        _params.desc = xml_node_attr(params, "desc");
        for (size_t i = 0; i < params->children->len; i++) {
            XMLNode *component = xml_node_child_at(params, i);

                        if (!strcmp(component->tag, "analysis")) {
                for (size_t j = 0; j < component->children->len; j++) {
                    XMLNode *analysis = xml_node_child_at(component, j);

                    if (!strcmp(analysis->tag, "savepaths")) {
                        _params.analysis.savePaths.t = std::stod(xml_node_attr(analysis, "t"));
                    }

                    if (!strcmp(analysis->tag, "stopping")) {
                        for (size_t k = 0; k < analysis->children->len; k++) {
                            XMLNode *bound = xml_node_child_at(analysis, k);

                            Vx::Model::Bound _bound{};

                            _bound.id = xml_node_attr(bound, "id");

                            auto conditions = tokenize(xml_node_attr(bound, "conditions"), " ");
                            for (const auto &condition: conditions) {
                                auto conditionStr = tokenize(condition, ":");

                                if (conditionStr.size() != 3) {
                                    throw std::runtime_error("Bad condition '" + condition + "', expected species:op:value");
                                }

                                Vx::Model::Condition _condition{};

                                for (const auto &_species: model.species) {
                                    if (_species.id == conditionStr[0]) {
                                        _condition.species = &_species;
                                        break;
                                    }
                                }
                                if (!_condition.species) {
                                    throw std::runtime_error("Unknown species '" + conditionStr[0] + "' in bound");
                                }

                                const char *comparison = conditionStr[1].c_str();
                                if (!strcmp(comparison, "=")) {
                                    _condition.op = Vx::Model::Comparison::Equal;
                                }
                                if (!strcmp(comparison, "<")) {
                                    _condition.op = Vx::Model::Comparison::LessThan;
                                }
                                if (!strcmp(comparison, ">")) {
                                    _condition.op = Vx::Model::Comparison::GreaterThan;
                                }
                                if (!strcmp(comparison, "<=")) {
                                    _condition.op = Vx::Model::Comparison::LessEqual;
                                }
                                if (!strcmp(comparison, ">=")) {
                                    _condition.op = Vx::Model::Comparison::GreaterEqual;
                                }

                                _condition.value = std::stoul(conditionStr[2]);

                                _bound.conditions.push_back(_condition);
                            }

                            _params.analysis.stopping.bounds.push_back(_bound);
                        }
                    }
                }
            }

            if (!strcmp(component->tag, "axis")) {
                Vx::Model::Axis axis{};

                for (size_t j = 0; j < component->children->len; j++) {
                    XMLNode *object = xml_node_child_at(component, j);

                    if (!strcmp(object->tag, "rate")) {
                        Vx::Model::DSweep sweep{};

                        if (xml_node_attr(object, "reaction")) {
                            std::string reactionId = xml_node_attr(object, "reaction");
                            for (const auto &_reactions: model.reactions) {
                                if (_reactions.id == reactionId) {
                                    sweep.reaction = &_reactions;
                                    break;
                                }
                            }
                        }

                        std::vector<std::string> tokens = tokenize(xml_node_attr(object, "values"), " ");
                        for (const auto &token: tokens) {
                            std::vector<std::string> _tokens = tokenize(token, ":");

                            if (_tokens.size() == 1) {
                                sweep.values.push_back(std::stod(_tokens[0]));
                            }
                            if (_tokens.size() == 3) {
                                std::vector<double> values;

                                if (!strcmp(xml_node_attr(object, "scale"), "log")) {
                                    values = geomspace(std::stod(_tokens[0]), std::stod(_tokens[1]),
                                                       std::stoul(_tokens[2]));
                                } else {
                                    values = linspace(std::stod(_tokens[0]), std::stod(_tokens[1]),
                                                      std::stoul(_tokens[2]));
                                }

                                sweep.values.insert(sweep.values.end(), values.begin(), values.end());
                            }
                        }
                        axis.rates.push_back(sweep);
                    } else {
                        Vx::Model::USweep sweep{};

                        if (xml_node_attr(object, "species")) {
                            std::string speciesId = xml_node_attr(object, "species");
                            for (const auto &_species: model.species) {
                                if (_species.id == speciesId) {
                                    sweep.species = &_species;
                                    break;
                                }
                            }
                        }

                        if (xml_node_attr(object, "reaction")) {
                            std::string reactionId = xml_node_attr(object, "reaction");
                            for (const auto &_reactions: model.reactions) {
                                if (_reactions.id == reactionId) {
                                    sweep.reaction = &_reactions;
                                    break;
                                }
                            }
                        }

                        std::vector<std::string> tokens = tokenize(xml_node_attr(object, "values"), " ");
                        for (const auto &token: tokens) {
                            std::vector<std::string> _tokens = tokenize(token, ":");

                            if (_tokens.size() == 1) {
                                sweep.values.push_back(std::stoul(_tokens[0]));
                            }
                            if (_tokens.size() == 3) {
                                std::vector<double> values;

                                if (!strcmp(xml_node_attr(object, "scale"), "log")) {
                                    values = geomspace(std::stod(_tokens[0]), std::stod(_tokens[1]),
                                                       std::stoul(_tokens[2]));
                                } else {
                                    values = linspace(std::stod(_tokens[0]), std::stod(_tokens[1]),
                                                      std::stoul(_tokens[2]));
                                }

                                sweep.values.insert(sweep.values.end(), values.begin(), values.end());
                            }
                        }

                        if (!strcmp(object->tag, "initial")) {
                            axis.initial.push_back(sweep);
                        }
                        if (!strcmp(object->tag, "product")) {
                            axis.products.push_back(sweep);
                        }
                        if (!strcmp(object->tag, "reactant")) {
                            axis.reactants.push_back(sweep);
                        }
                    }
                }
                _params.axes.push_back(axis);
            }
        }
        return _params;
    }

    struct VIRTExArgs : argparse::Args {
        std::filesystem::path &model = arg("modelPath", "path to model xml file");
        std::filesystem::path &params = arg("sweepPath", "path to sweep xml file");
        std::filesystem::path &output = kwarg("o,output", "path to output directory").set_default(".");
        bool &verbose = flag("v,verbose", "toggle verbose output");
    };
}

int main(int argc, char *argv[]) {
    auto args = argparse::parse<VIRTExArgs>(argc, argv);

    if (args.verbose) {
        args.print();
    }

    Vx::Model::Model model = parseModel(loadFile(args.model));

    Vx::Model::Params params = parseParams(loadFile(args.params), model);

    if (!params.seed) {
        std::random_device rd;
        params.seed = rd();
    }

    if (args.verbose) {
        std::cout << "\n" << model << "\n" << params << "\n";
    }

    Vx::SimObject sim = Vx::SimObject(params, Vx::SimType::StochasticSimulation);

    sim.runSimulation();
    Vx::SimOutput* out = sim.deviceSynchronize();

    auto paths = out->getPaths();
    for (const auto &path: paths) {
        std::cout << path << "\n";
    }

    return 0;
}
