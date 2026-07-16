/**
   * @brief Partial wave analysis header and functions 
   * @author Kyle Salamone 
   * @date July 2026  
 */

#ifndef PartialWaveGen_H
#define PartialWaveGen_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <array>
#include <random>
#include <cmath>
#include <numbers>
#include <memory>
#include <optional>

#include "pi_p_amplitude.h"
#include "io.h"


// TODO: only works for pi+p currently. Need to adjust for pi-p
class PartialWaveGen
{
    private:
        double momentum_lab{0.0}, momentum_cm{0.0};
        bool verbose{false};
        double gamma_boost{0.0};
        double alpha_kinematic{0.0};
        charge charge_polarity;
        double theta_min{0.}, theta_max{0.};
        int L_MAX{1};
        std::mt19937 mersenne_twister;
        std::uniform_real_distribution<double> uniform_dist, costheta_dist;
        std::mt19937::result_type seed {987654321};
        double bin_size{0.0001}; // Step in cos(theta)
        std::vector<double> dsigma_domega_costheta; // For comparison
        double umin{0.}, umax{0.};
        double max_cs{0.};
        bool lab_frame{false}; // !lab_frame = cm_frame
        PiPAmplitude pi_p_amplitude;
    public:
        PartialWaveGen() = default;
        ~PartialWaveGen() = default;
        /**
         * @brief Constructor for partial wave analysis
         * @param momentum_lab_ Lab momentum (in MeV)
         * @param ch Charge of the pion (+1 for pi+, -1 for pi-)
         * @param phase_shift_file std::optional path to the phase shift file. Defaults to stored ones in data/ based on lab momentum
         * @param theta_min_ Minimum scattering angle, in LAB frame
         * @param theta_max_ Maximum scattering angle, in LAB frame
         * @param L_MAX Maximum partial wave to calculate. Default: p-wave (L_MAX = 1)
         * @param verbose_ Whether to print verbose output
         * @param seed_ Seed for the random number generator
         * @param bin_size_ Size of the bins for caching 
         */
        PartialWaveGen(const double momentum_lab_, const charge ch, const double theta_min_, const double theta_max_, const std::optional<std::string> &phase_shift_file_path = std::nullopt, const int L_MAX_ = 1, bool verbose_ = false, std::mt19937::result_type seed_ = 987654321, double bin_size_ = 0.0001)
            : momentum_lab(momentum_lab_), verbose(verbose_), charge_polarity(ch), L_MAX(L_MAX_), seed(seed_), bin_size(bin_size_)
        {
            if(charge_polarity != charge::plus && charge_polarity != charge::minus)
            {
                throw std::invalid_argument("Invalid charge polarity");
            }
            set_seed(seed_);
            auto phase_shift_file = phase_shift_file_path.value_or(io::get_phase_shift_file_path(momentum_lab));
            pi_p_amplitude = PiPAmplitude(phase_shift_file, momentum_lab_, charge_polarity, verbose_);
            gamma_boost = pi_p_amplitude.get_gamma_cm_boost();
            alpha_kinematic = pi_p_amplitude.get_alpha_kinematic();
            momentum_cm = pi_p_amplitude.get_momentum_cm();
            theta_min = physics_helpers::theta_lab_to_cm(theta_min_, gamma_boost, alpha_kinematic);
            theta_max = physics_helpers::theta_lab_to_cm(theta_max_, gamma_boost, alpha_kinematic);
            umax = std::cos(physics_helpers::deg2rad(theta_min));
            umin = std::cos(physics_helpers::deg2rad(theta_max));
            uniform_dist = std::uniform_real_distribution<double>(0.0, 1.0);
            costheta_dist = std::uniform_real_distribution<double>(umin, umax);
            int n_bins = static_cast<int>(std::round((umax - umin) / bin_size)) + 1;
            dsigma_domega_costheta.resize(n_bins);
            for(int i=0; i<n_bins; i++)
            {
                double cos_theta = umin + i * bin_size;
                double theta = physics_helpers::rad2deg(std::acos(cos_theta));
                dsigma_domega_costheta[i] = pi_p_amplitude.dsigma_domega_cm(theta);
                if(dsigma_domega_costheta[i] > max_cs) max_cs = dsigma_domega_costheta[i];
            }
            max_cs *= 1.1;
        }
        /**
         * @brief Sets the seed for mersenne twister
         */
        void set_seed(const std::mt19937::result_type seed) noexcept
        {
            this->seed = seed;
            mersenne_twister.seed(this->seed);
        }
        /**
         * @brief Sets output to be lab frame values 
         */
        void set_lab_frame() 
        {
            this->lab_frame = true;
        }
        /**
         * @brief Sets output to be CM frame values 
         */
        void set_cm_frame()
        {
            this->lab_frame = false;
        }
        /**
         * @brief Returns whether generator is in lab frame or CM frame
         * @returns True if in lab frame, false if in CM frame
         */
        bool is_in_lab_frame() const noexcept
        {
            return lab_frame;
        }
        /**
         * @brief Gets minimum theta in correct frame 
         * @returns Minimum theta in the correct frame, in degrees
         */
        double get_min_theta() const noexcept
        {
            if(lab_frame)
            {
                return physics_helpers::theta_cm_to_lab(theta_min, gamma_boost, alpha_kinematic);
            }
            return theta_min;
        }
        /**
         * @brief Gets maximum theta in correct frame 
         * @returns Maximum theta in the correct frame, in degrees
         */
        double get_max_theta() const noexcept
        {
            if(lab_frame)
            {
                return physics_helpers::theta_cm_to_lab(theta_max, gamma_boost, alpha_kinematic);
            }
            return theta_max;
        }
        /**
         * @brief Returns CM frame minimum theta, in degrees, always
         * @returns CM frame minimum theta, in degrees
         */
        double get_min_theta_cm() const noexcept
        {
            return theta_min;
        }
        /**
         * @brief Returns CM frame maximum theta, in degrees, always
         * @returns CM frame maximum theta, in degrees
         */
        double get_max_theta_cm() const noexcept
        {
            return theta_max;
        }
        /**
         * @brief Returns reference to object to calculate PiPAmplitude object
         * @returns Reference to PiPAmplitude object
         */
        const PiPAmplitude& get_pi_p_amplitude() const noexcept
        {
            return pi_p_amplitude;
        }
        /**
         * @brief Samples random uniform phi in accordance with theory
         * @returns Random phi value, in radians 
         */
        double sample_phi() noexcept
        {
            return 2.0 * std::numbers::pi * uniform_dist(mersenne_twister);
        }
        /**
         * @brief Samples random theta according to SAID DCS
         * @brief \f[ p(\cos(\theta)) \propto \frac{d\sigma}{d\Omega} \f]
         * @returns Random theta value in the CM frame, in radians 
         */
        double sample_theta() noexcept 
        {
            double theta = 0.0;
            while(true)
            {
                double costheta_sample = costheta_dist(mersenne_twister);
                double offset = (costheta_sample - umin) / bin_size;
                std::size_t idx_lo = static_cast<std::size_t>(std::floor(offset + 1e-9));
                std::size_t idx_hi = idx_lo + 1;
                if(idx_hi >= dsigma_domega_costheta.size())
                {
                    idx_hi = dsigma_domega_costheta.size() - 1;
                }
                double costheta_lo = umin + idx_lo * bin_size;
                double frac = (costheta_sample - costheta_lo) / bin_size;
                double cs_lo = dsigma_domega_costheta[idx_lo];
                double cs_hi = dsigma_domega_costheta[idx_hi];
                double dsigma = cs_lo + frac * (cs_hi - cs_lo);
                double prob = std::fmin(1.0, dsigma / max_cs);
                double dice = uniform_dist(mersenne_twister);
                if(dice < prob)
                {
                    theta = physics_helpers::rad2deg(std::acos(costheta_sample));
                    break;
                }
            }
            return physics_helpers::deg2rad(theta);
        }
        /**
         * @brief Returns the outgoing pion momentum in the proper frame.
         * @param theta: CM scattering angle, in degrees
         * @param phi: Azimuthal angle, in degrees
         * @returns The outgoing pion four-vector in the given frame (in MeV)
         */
        fourVector get_pion_four_vector(const double theta, const double phi)
        {
            const double k = momentum_cm;
            const double theta_rad = physics_helpers::deg2rad(theta);
            const double phi_rad = physics_helpers::deg2rad(phi);
            const double px = k * std::sin(theta_rad) * std::cos(phi_rad);
            const double py = k * std::sin(theta_rad) * std::sin(phi_rad);
            const double pz = k * std::cos(theta_rad);
            const double E = std::sqrt(k*k + physics_helpers::m_pion_squared);
            fourVector pion_four_vector(E, px, py, pz);
            if(lab_frame)
            {
                pion_four_vector = physics_helpers::boost_cm_to_lab(pion_four_vector, gamma_boost);
            }
            return pion_four_vector;
        }   
        /**
         * @brief Returns the recoil proton 4 vector 
         * @param theta: CM scattering angle, in degrees
         * @param phi: Azimuthal angle, in degrees
         * @returns The recoil proton four-vector in the given frame (in MeV)
         */
        fourVector get_proton_four_vector(const double theta, const double phi)
        {
            const double k = momentum_cm;
            const double theta_rad = physics_helpers::deg2rad(theta);
            const double phi_rad = physics_helpers::deg2rad(phi);
            const double px = -k * std::sin(theta_rad) * std::cos(phi_rad);
            const double py = -k * std::sin(theta_rad) * std::sin(phi_rad);
            const double pz = -k * std::cos(theta_rad);
            const double E = std::sqrt(k*k + physics_helpers::m_proton_squared);
            fourVector proton_four_vector(E, px, py, pz);
            if(lab_frame)
            {
                proton_four_vector = physics_helpers::boost_cm_to_lab(proton_four_vector, gamma_boost);
            }
            return proton_four_vector;
        }
};  

#endif 