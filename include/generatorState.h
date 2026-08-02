#ifndef GENERATOR_STATE_H
#define GENERATOR_STATE_H

#include "pi_p_amplitude.h"
#include "physics_helpers.h"

#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <numbers>

/**
 * @brief Immutable momentum-dependent state for the partial wave generator.
 *
 * This class contains everything that depends on the beam momentum:
 *
 * - PiPAmplitude
 * - CM momentum
 * - Lorentz boost parameters
 * - angular limits
 * - differential cross section lookup table
 * - rejection sampling normalization
 *
 * A generatorState should be created once and then reused.
 * This makes it suitable for momentum caching:
 *
 * momentum -> generatorState
 */
class generatorState
{
    private:
        double momentum_lab{0.0};
        double momentum_cm{0.0};
        double pion_lab_energy{0.0};
        double beta_boost{0.0};
        double gamma_boost{0.0};
        double alpha_kinematic{0.0};
        double theta_lab_min{0.0};
        double theta_lab_max{0.0};
        double theta_cm_min{0.0};
        double theta_cm_max{0.0};
        double umin{0.0};
        double umax{0.0};
        double bin_size{0.0001};
        int L_MAX{0};
        PiPAmplitude pi_p_amplitude;
        std::vector<double> dsigma_domega_costheta;
        double max_cs{0.0};
        // In fm^2
        double total_cs{0.0};

    public:
        /**
         * @brief Construct a momentum state.
         *
         * This performs all expensive calculations:
         * - phase shift interpolation
         * - amplitude initialization
         * - differential cross section lookup generation
         */
        generatorState(const std::shared_ptr<const PhaseShiftLoader>& phase_shift_loader, int L_MAX_, const double momentum_lab_,const charge ch,const double theta_min_lab_,const double theta_max_lab_,const double bin_size_=0.0001)
            : momentum_lab(momentum_lab_), theta_lab_min(theta_min_lab_), theta_lab_max(theta_max_lab_), bin_size(bin_size_), L_MAX(L_MAX_),
            pi_p_amplitude(PiPAmplitude(phase_shift_loader, momentum_lab_, ch, true, L_MAX_))
        {
            initialize_kinematics();
            initialize_angles();
            build_lookup_table();
        }
    private:
        /**
         * @brief Calculate momentum-dependent kinematics.
         */
        inline void initialize_kinematics()
        {
            pion_lab_energy = physics_helpers::pion_lab_energy(momentum_lab);
            beta_boost = physics_helpers::beta_boost(momentum_lab, pion_lab_energy);
            momentum_cm = pi_p_amplitude.get_momentum_cm();
            gamma_boost = pi_p_amplitude.get_gamma_cm_boost();
            alpha_kinematic = pi_p_amplitude.get_alpha_kinematic();
        }
        /**
         * @brief Convert angular limits from lab to CM frame.
         */
        inline void initialize_angles()
        {
            theta_cm_min = physics_helpers::theta_lab_to_cm(theta_lab_min, gamma_boost, alpha_kinematic);
            theta_cm_max = physics_helpers::theta_lab_to_cm(theta_lab_max, gamma_boost, alpha_kinematic);
            umax = std::cos(physics_helpers::deg2rad(theta_cm_min));
            umin = std::cos(physics_helpers::deg2rad(theta_cm_max));
        }
        /**
         * @brief Build rejection sampling lookup table.
         */
        inline void build_lookup_table()
        {
            int n_bins = static_cast<int>(std::round((umax - umin) / bin_size)) + 1;
            dsigma_domega_costheta.resize(n_bins);
            max_cs = 0.0;
            for(int i=0; i<n_bins; i++)
            {
                double cos_theta = umin + i * bin_size;
                double theta = physics_helpers::rad2deg(std::acos(cos_theta));
                dsigma_domega_costheta[i] = pi_p_amplitude.dsigma_domega_cm(theta);
                if(dsigma_domega_costheta[i] > max_cs)
                {
                    max_cs = dsigma_domega_costheta[i];
                }
            }
            max_cs *= 1.1;
            auto dcs = [this](double cos_theta) -> double
            {
                double theta_deg = physics_helpers::rad2deg(std::acos(cos_theta));
                return pi_p_amplitude.dsigma_domega_cm(theta_deg);
            };
            total_cs = 2.0 * std::numbers::pi * physics_helpers::adaptive_simpson(umin, umax, dcs, 1e-4);
        }

    public:
        inline double get_momentum_lab() const noexcept
        {
            return momentum_lab;
        }
        inline double get_momentum_cm() const noexcept
        {
            return momentum_cm;
        }
        inline double get_beta_boost() const noexcept
        {
            return beta_boost;
        }
        inline double get_gamma_boost() const noexcept
        {
            return gamma_boost;
        }
        inline double get_alpha_kinematic() const noexcept
        {
            return alpha_kinematic;
        }
        inline double get_theta_lab_min() const noexcept
        {
            return theta_lab_min;
        }
        inline double get_theta_lab_max() const noexcept
        {
            return theta_lab_max;
        }
        inline double get_theta_cm_min() const noexcept
        {
            return theta_cm_min;
        }
        inline double get_theta_cm_max() const noexcept
        {
            return theta_cm_max;
        }
        inline double get_umin() const noexcept
        {
            return umin;
        }
        inline double get_umax() const noexcept
        {
            return umax;
        }
        inline double get_total_cs() const noexcept
        {
            return total_cs;
        }
        inline const std::vector<double> &get_dsigma_domega_costheta() const noexcept
        {
            return dsigma_domega_costheta;
        }
        inline double get_max_cs() const noexcept
        {
            return max_cs;
        }
        inline const PiPAmplitude &get_amplitude() const noexcept
        {
            return pi_p_amplitude;
        }
        inline double get_bin_size() const noexcept
        {
            return bin_size;
        }
};

#endif 