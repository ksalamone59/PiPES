#ifndef IO_H
#define IO_H

#include <filesystem>
#include <format>
#include <stdexcept>
#include <string_view>

namespace io
{
    /**
     * @brief Loads a table of values from a file into a 2D array
     * @tparam N_ROWS: Number of rows in the table
     * @tparam N_COLS: Number of columns in the table
     * @param file_path: Path to the file containing the table
     * @returns A 2D array containing the loaded table values
     */
    template<std::size_t N_ROWS, std::size_t N_COLS>
    inline std::array<std::array<double, N_COLS>, N_ROWS> load_table(const std::filesystem::path &file_path)
    {
        std::array<std::array<double, N_COLS>, N_ROWS> table{};
        std::ifstream file(file_path);
        if(!file.is_open())
        {
            throw std::runtime_error("Failed to open file: " + file_path.string());
        }
        std::string line;
        std::size_t row_idx{0};
        while(std::getline(file, line) && row_idx < N_ROWS)
        {
            if(line.empty() || line.starts_with("#")) continue;
            std::istringstream ss(line);
            for(std::size_t col_idx = 0; col_idx < N_COLS; col_idx++)
            {
                ss >> table[row_idx][col_idx];
            }
            row_idx++;
        }
        return table;
    }
}

#endif