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

#include "../../test_scf.hpp"
#include <iostream>
#include <pluginplay/pluginplay.hpp>

using namespace scf;

using pt = simde::MolecularGrid;

TEST_CASE("GridFromFile") {
    pluginplay::ModuleManager mm;
    scf::load_modules(mm);
    auto& mod = mm.at("Grid From File");

    int process_rank;
    MPI_Comm_rank(mm.get_runtime().mpi_comm(), &process_rank);
    auto rank_str = std::to_string(process_rank);

    auto h2 = test_scf::make_h2<chemist::Molecule>();

    SECTION("file doesn't exist") {
        std::filesystem::path p = "non_existent_file.txt";
        mod.change_input("Path to Grid File", p);
        REQUIRE_THROWS_AS(mod.run_as<pt>(h2), std::runtime_error);
    }

    SECTION("only three points") {
        std::filesystem::path p = "only_three_points" + rank_str + ".txt";
        std::ofstream file(p);
        file << "1.0 2.0 3.0\n";
        file << "5.0 6.0 7.0 8.0\n";
        file.close();
        mod.change_input("Path to Grid File", p);
        REQUIRE_THROWS_AS(mod.run_as<pt>(h2), std::runtime_error);
        std::filesystem::remove(p);
    }

    SECTION("five numbers points") {
        std::filesystem::path p = "five_numbers_points" + rank_str + ".txt";
        std::ofstream file(p);
        file << "1.0 2.0 3.0 4.0 5.0\n";
        file << "5.0 6.0 7.0 8.0\n";
        file.close();
        mod.change_input("Path to Grid File", p);
        REQUIRE_THROWS_AS(mod.run_as<pt>(h2), std::runtime_error);
        std::filesystem::remove(p);
    }

    SECTION("non-mumeric data") {
        std::filesystem::path p = "non-numeric_data" + rank_str + ".txt";
        std::ofstream file(p);
        file << "1.0 2.0 3.0 hello\n";
        file << "5.0 6.0 7.0 8.0\n";
        file.close();
        mod.change_input("Path to Grid File", p);
        REQUIRE_THROWS_AS(mod.run_as<pt>(h2), std::runtime_error);
        std::filesystem::remove(p);
    }

    SECTION("good file") {
        std::filesystem::path p = "good_file" + rank_str + ".txt";
        std::ofstream file(p);
        file << "1.0 2.0 3.0 4.0\n";
        file << "5.0 6.0 7.0 8.0\n";
        file.close();
        mod.change_input("Path to Grid File", p);
        auto grid = mod.run_as<pt>(h2);
        REQUIRE(grid.size() == 2);
        REQUIRE(grid.at(0).weight() == 4.0);
        REQUIRE(grid.at(0).point().x() == 1.0);
        REQUIRE(grid.at(0).point().y() == 2.0);
        REQUIRE(grid.at(0).point().z() == 3.0);
        REQUIRE(grid.at(1).weight() == 8.0);
        REQUIRE(grid.at(1).point().x() == 5.0);
        REQUIRE(grid.at(1).point().y() == 6.0);
        REQUIRE(grid.at(1).point().z() == 7.0);
        std::filesystem::remove(p);
    }

    SECTION("skip blank line") {
        std::filesystem::path p = "skip_blank_line" + rank_str + ".txt";
        std::ofstream file(p);
        file << "1.0 2.0 3.0 4.0\n";
        file << "\n";
        file << "5.0 6.0 7.0 8.0\n";
        file.close();
        mod.change_input("Path to Grid File", p);
        auto grid = mod.run_as<pt>(h2);
        REQUIRE(grid.size() == 2);
        REQUIRE(grid.at(0).weight() == 4.0);
        REQUIRE(grid.at(0).point().x() == 1.0);
        REQUIRE(grid.at(0).point().y() == 2.0);
        REQUIRE(grid.at(0).point().z() == 3.0);
        REQUIRE(grid.at(1).weight() == 8.0);
        REQUIRE(grid.at(1).point().x() == 5.0);
        REQUIRE(grid.at(1).point().y() == 6.0);
        REQUIRE(grid.at(1).point().z() == 7.0);
        std::filesystem::remove(p);
    }
}
