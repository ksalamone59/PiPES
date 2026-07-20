#ifndef PHASE_SHIFT_GEN_H
#define PHASE_SHIFT_GEN_H

#include "physics_helpers.h"
#include "io.h"
#include <filesystem>
#include "partial_wave.h"

class PhaseShiftLoader
{
    public:
        PhaseShiftLoader()
        {
            phase_shift_data = io::load_table<51, 7>(table_file);
        }
        ~PhaseShiftLoader() = default;
        /**
         * @brief Returns s and p phase shifts for a given momentum based on table 
         * @param momentum_lab: The lab momentum in MeV
         * @param verbose: Whether to print the phase shifts to the console
         * @returns An array of phase shifts in degrees: {S11, S31, P11, P13, P31, P33}
         */
        partial_wave set_s_p_phase_shifts(const double momentum_lab, bool verbose = false) const
        {
            partial_wave phase_shifts;
            phase_shifts.lmax = 1;
            auto results = physics_helpers::interpolate(momentum_lab, phase_shift_data);
            phase_shifts[static_cast<std::size_t>(wave::S11)] = physics_helpers::deg2rad(results[0]);
            phase_shifts[static_cast<std::size_t>(wave::S31)] = physics_helpers::deg2rad(results[1]);
            phase_shifts[static_cast<std::size_t>(wave::P11)] = physics_helpers::deg2rad(results[2]);
            phase_shifts[static_cast<std::size_t>(wave::P13)] = physics_helpers::deg2rad(results[3]);
            phase_shifts[static_cast<std::size_t>(wave::P31)] = physics_helpers::deg2rad(results[4]);
            phase_shifts[static_cast<std::size_t>(wave::P33)] = physics_helpers::deg2rad(results[5]);
            if(verbose)
            {
                std::cout << "Phase shifts for momentum " << momentum_lab << " MeV:" << std::endl;
                std::cout << "S11: " << physics_helpers::rad2deg(phase_shifts[static_cast<std::size_t>(wave::S11)]) << " deg" << std::endl;
                std::cout << "S31: " << physics_helpers::rad2deg(phase_shifts[static_cast<std::size_t>(wave::S31)]) << " deg" << std::endl;
                std::cout << "P11: " << physics_helpers::rad2deg(phase_shifts[static_cast<std::size_t>(wave::P11)]) << " deg" << std::endl;
                std::cout << "P13: " << physics_helpers::rad2deg(phase_shifts[static_cast<std::size_t>(wave::P13)]) << " deg" << std::endl;
                std::cout << "P31: " << physics_helpers::rad2deg(phase_shifts[static_cast<std::size_t>(wave::P31)]) << " deg" << std::endl;
                std::cout << "P33: " << physics_helpers::rad2deg(phase_shifts[static_cast<std::size_t>(wave::P33)]) << " deg" << std::endl;
            }
            return phase_shifts;
        }
    private:
        std::array<std::array<double, 7>, 51> phase_shift_data;
        const std::filesystem::path table_file = std::filesystem::path(PHASE_SHIFT_DATADIR) / "phase_shifts.dat";
};

#endif 