////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2026 Cooper Gray
//
//  This application is free software, meaning you can redistribute it and/or
//  modify it under the terms of the Apache License Version 2.0. See `LICENSE`
//  for details. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0

#include <nanobind/nanobind.h>

#include <virtex/Model.hpp>
#include <virtex/SimObject.hpp>

#include "Helpers.hpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(virtex, m) {
    m.doc() = "Viral Infection, Recovery, and Transmission Explorer";

    nb::class_<Vx::Model::Model>(m, "Model")
        .def(nb::init<>());

    nb::class_<Vx::Model::Params>(m, "Params")
        .def(nb::init<>());

    m.def("load_file", &Vx::Helpers::loadFile, "filename"_a,
          "Read a file into a string.");

    m.def("parse_model", &Vx::Helpers::parseModel, "xml_string"_a,
          "Parse a model description from an XML string.");

    m.def("parse_params", &Vx::Helpers::parseParams, "xml_string"_a, "model"_a,
          nb::keep_alive<0, 2>(),
          "Parse simulation parameters from an XML string.");

    nb::enum_<Vx::SimType>(m, "SimType")
        .value("StochasticSimulation", Vx::SimType::StochasticSimulation);

    nb::class_<Vx::SimOutput>(m, "SimOutput")
        .def(nb::init<const Vx::Model::Params &, int>(),
             "params"_a, "num_trials"_a,
             nb::keep_alive<1, 2>())
        .def("get_paths", &Vx::SimOutput::getPaths)
        .def("get_exits", &Vx::SimOutput::getExits)
        .def("get_times", &Vx::SimOutput::getTimes);

    nb::class_<Vx::SimObject>(m, "SimObject")
        .def(nb::init<const Vx::Model::Params &, Vx::SimType>(),
             "params"_a, "type"_a,
             nb::keep_alive<1, 2>())
        .def("run_simulation", &Vx::SimObject::runSimulation,
             nb::call_guard<nb::gil_scoped_release>())
        .def("device_synchronize", &Vx::SimObject::deviceSynchronize,
             nb::rv_policy::reference_internal);
}
