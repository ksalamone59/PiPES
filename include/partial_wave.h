/**
   * @brief Partial wave handling
   * @author Kyle Salamone 
   * @date July 2026 
 */

#ifndef PARTIAL_WAVE_H
#define PARTIAL_WAVE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string_view>
#include <algorithm>
#include <sstream>

#include "physics_helpers.h"

enum class wave {S11, S31, P11, P13, P31, P33, D13, D15, D33, D35, F15, F17, F35, F37};
/**
 * @brief Converts a wave index to a string
 * @param w: The wave index to convert
 * @returns A string representing the wave index
 */
constexpr std::string_view wave_idx_to_string(const wave w)
{
    switch(w)
    {
        case wave::S11: return "S11";
        case wave::S31: return "S31";
        case wave::P11: return "P11";
        case wave::P13: return "P13";
        case wave::P31: return "P31";
        case wave::P33: return "P33";
        case wave::D13: return "D13";
        case wave::D15: return "D15";
        case wave::D33: return "D33";
        case wave::D35: return "D35";
        case wave::F15: return "F15";
        case wave::F17: return "F17";
        case wave::F35: return "F35";
        case wave::F37: return "F37";
        default: throw std::invalid_argument("Unknown wave");
    }
    return "Unknown";
}

struct partial_wave
{
    // Phase shifts 
    double s11{0.}, s31{0.};
    double p11{0.}, p13{0.}, p31{0.}, p33{0.};
    double d13{0.}, d15{0.}, d33{0.}, d35{0.};
    double f15{0.}, f17{0.}, f35{0.}, f37{0.};
    double lmax = 0;
    partial_wave() = default;
    ~partial_wave() = default;
    /**
     * @brief Returns a reference to the phase shift at the specified index
     * @param idx: The index of the phase shift to return
     * @returns A reference to the phase shift at the specified index
     */
    double& operator[] (const std::size_t idx) 
    {
        switch(idx)
        {
            case 0: return s11; // I=1/2, J=1/2
            case 1: return s31; // I=3/2, J=1/2
            case 2: return p11; // I=1/2, J=1/2
            case 3: return p13; // I=1/2, J=3/2
            case 4: return p31; // I=3/2, J=1/2
            case 5: return p33; // I=3/2, J=3/2
            case 6: return d13; // I=1/2, J=3/2
            case 7: return d15; // I=1/2, J=5/2
            case 8: return d33; // I=3/2, J=3/2
            case 9: return d35; // I=3/2, J=5/2
            case 10: return f15; // I=1/2, J=5/2
            case 11: return f17; // I=1/2, J=7/2
            case 12: return f35; // I=3/2, J=5/2
            case 13: return f37; // I=3/2, J=7/2
            default: throw std::out_of_range("Invalid partial wave index.");
        }
    }
    /**
     * @brief Returns a constant reference to the phase shift at the specified index
     * @param idx: The index of the phase shift to return
     * @returns A constant reference to the phase shift at the specified index
     */
    const double& operator[](const std::size_t idx) const
    {
        switch(idx)
        {
            case 0: return s11; // I=1/2, J=1/2
            case 1: return s31; // I=3/2, J=1/2
            case 2: return p11; // I=1/2, J=1/2
            case 3: return p13; // I=1/2, J=3/2
            case 4: return p31; // I=3/2, J=1/2
            case 5: return p33; // I=3/2, J=3/2
            case 6: return d13; // I=1/2, J=3/2
            case 7: return d15; // I=1/2, J=5/2
            case 8: return d33; // I=3/2, J=3/2
            case 9: return d35; // I=3/2, J=5/2
            case 10: return f15; // I=1/2, J=5/2
            case 11: return f17; // I=1/2, J=7/2
            case 12: return f35; // I=3/2, J=5/2
            case 13: return f37; // I=3/2, J=7/2
            default: throw std::out_of_range("Invalid partial wave index.");
        }
    }
    /**
     * @brief Transforms all phase shifts from degrees to radians
     */
    void transform_to_radians() 
    {
        s11 = physics_helpers::deg2rad(s11);
        s31 = physics_helpers::deg2rad(s31);
        p11 = physics_helpers::deg2rad(p11);
        p13 = physics_helpers::deg2rad(p13);
        p31 = physics_helpers::deg2rad(p31);
        p33 = physics_helpers::deg2rad(p33);
        d13 = physics_helpers::deg2rad(d13);
        d15 = physics_helpers::deg2rad(d15);
        d33 = physics_helpers::deg2rad(d33);
        d35 = physics_helpers::deg2rad(d35);
        f15 = physics_helpers::deg2rad(f15);
        f17 = physics_helpers::deg2rad(f17);
        f35 = physics_helpers::deg2rad(f35);
        f37 = physics_helpers::deg2rad(f37);
    }
    /**
     * @brief Sets all phase shifts to zero
     */
    void set_all_zero()
    {
        s11 = 0.0;
        s31 = 0.0;
        p11 = 0.0;
        p13 = 0.0;
        p31 = 0.0;
        p33 = 0.0;
        d13 = 0.0;
        d15 = 0.0;
        d33 = 0.0;
        d35 = 0.0;
        f15 = 0.0;
        f17 = 0.0;
        f35 = 0.0;
        f37 = 0.0;
    }
    auto operator<=>(const partial_wave& other) const = default;
    /**
     * @brief Returns the number of phase shifts USED in the analysis, not the number loaded 
     * @returns The number of phase shifts USED in the analysis
     */
    inline std::size_t size() const noexcept 
    {
        return 2 + 4 * lmax; 
    }
    /**
     * @brief Returns number of phase shifts loaded, not necessarily the number used in the analysis 
     * @returns The number of phase shifts loaded, not necessarily the number used in the analysis
     */
    inline constexpr std::size_t loaded_size() const noexcept { return 14; }
};

/**
 * @brief Loads phase shifts from a data file
 * @param filename: The name of the file to load
 * @returns A partial_wave object containing the loaded phase shifts
 */
inline partial_wave load_phase_shifts(const std::string &filename, int lmax = 1)
{
    std::ifstream file(filename, std::ios_base::in);
    if(!file.is_open())
    {
        std::cerr << "Failed to open phase shift file: " << filename << "\n";
        throw std::runtime_error("Failed to open phase shift file.");
    }
    std::string line;
    while(true)
    {
        std::getline(file, line); 
        if(!line.starts_with('#') && !line.empty())
        {
            break;
        }
    }
    std::istringstream ss(line);
    partial_wave pw;
    pw.lmax = lmax;
    ss >> pw.s11 >> pw.s31 >> pw.p11 >> pw.p13 >> pw.p31 >> pw.p33 >> pw.d13 >> pw.d15 >> pw.d33 >> pw.d35 >> pw.f15 >> pw.f17 >> pw.f35 >> pw.f37;
    std::cout << "Loaded the following phase shifts (degrees) from " << filename << ":\n";
    for(std::size_t i=0;i<pw.size();i++) 
    {
        std::cout << "phase_shifts[" << wave_idx_to_string(static_cast<wave>(i)) << "]: " << pw[i] << "\n";
    }
    pw.transform_to_radians();
    return pw;
}

#endif 