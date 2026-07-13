#ifndef IO_H
#define IO_H

#include <filesystem>
#include <format>
#include <stdexcept>
#include <string_view>

namespace io
{
    /**
     * @brief Gets the file path for the phase shift data based on lab momentum.
     * @param lab_momentum The lab momentum.
     * @return The file path for the phase shift data.
     */
    inline std::filesystem::path get_phase_shift_file_path(const double lab_momentum)
    {
        #ifndef PHASE_SHIFT_DATADIR
        throw std::runtime_error("Phase shift data directory not defined! Please check CMake.");
        #endif 
        std::filesystem::path file_path = PHASE_SHIFT_DATADIR;
        return (file_path  / std::format("phase_shifts_{}.dat", lab_momentum));
    }
}

#endif