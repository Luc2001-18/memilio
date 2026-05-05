/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "pybind_util.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#define OSEIRVECTOR_BINDINGS_SKIP_MAIN
#include "ode_seir_malaria.cpp"

namespace py = pybind11;

PYBIND11_MODULE(_simulation_oseirvector, m)
{
    m.def("simulate", &simulate, "Simulates the OSEIR Vector metapopulation model", py::arg("t0") = 0,
          py::arg("tmax") = 6940, py::arg("dt") = 0.1, py::arg("TimeExposed") = 15.0,
          py::arg("TimeInfectedAsymptomatic") = 100.0, py::arg("TimeInfectedSymptomatic") = 7.0,
          py::arg("TransmissionProbabilityOnContact") = 0.1, py::arg("AsymptomaticProbability") = 0.287,
          py::arg("TimeWaningImmunity") = 180.0, py::arg("BitingRateNorth") = 0.4, py::arg("BitingRateCenter") = 0.4,
          py::arg("BitingRateSouth") = 0.4, py::arg("TransmissionVectorToHuman") = 0.24, 
          py::arg("TransmissionHumanToVector") = 0.02,
          py::arg("ic_scale") = 1.0);
    m.attr("__version__") = "dev";
}
