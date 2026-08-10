/*
 * Copyright 2025 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "xc.hpp"
#include <filesystem>
#include <simde/integration_grids/molecular_grid.hpp>
#include <sstream>

namespace scf::xc {
namespace {
const auto desc = R"(

GridFromFile
------------

Reads a grid from a file. The file is expected to be a text file with four
columns: x, y, z, weight. Each row of the file corresponds to a grid point.
The columns should be space separated. An example of the file format is:

.. code-block:: text

    1.0 2.0 3.0 4.0
    5.0 6.0 7.0 8.0

This would create two grid points: (1.0, 2.0, 3.0) with weight 4.0 and
(5.0, 6.0, 7.0) with weight 8.0.
)";

}

using pt = simde::MolecularGrid;

MODULE_CTOR(GridFromFile) {
    satisfies_property_type<pt>();

    add_input<std::filesystem::path>("Path to Grid File");

    description(desc);
}

MODULE_RUN(GridFromFile) {
    // const auto& [molecule] = pt::unwrap_inputs(inputs);

    const auto& path =
      inputs.at("Path to Grid File").value<std::filesystem::path>();

    std::ifstream file(path);
    if(!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + path.string());
    }

    // Create the grid
    std::string line;
    std::vector<chemist::GridPoint> grid_points;
    while(std::getline(file, line)) {
        auto first_character = line.find_first_not_of(" \t\r\n");
        if(first_character == std::string::npos) continue;
        double x, y, z, weight;
        std::stringstream ss(line);
        if(!(ss >> x >> y >> z >> weight) || !(ss >> std::ws).eof()) {
            throw std::runtime_error("Malformed grid file line: '" + line +
                                     "'");
        }
        grid_points.emplace_back(weight, x, y, z);
    }
    chemist::Grid grid(grid_points.begin(), grid_points.end());

    file.close();

    // Return the grid
    auto rv = results();
    return pt::wrap_results(rv, std::move(grid));
}
} // namespace scf::xc
