#ifndef GENERATOR_STATE_H
#define GENERATOR_STATE_H

#include "pi_p_amplitude.h"
#include "physics_helpers.h"

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
        // Momentum-dependent kinematics
        double momentum_lab{0.0};
        double momentum_cm{0.0};
        double pion_lab_energy{0.0};

        double beta_boost{0.0};
        double gamma_boost{0.0};
        double alpha_kinematic{0.0};
        // Angular limits
        double theta_lab_min{0.0};
        double theta_lab_max{0.0};
        double theta_cm_min{0.0};
        double theta_cm_max{0.0};
        // cos(theta) sampling limits
        double umin{0.0};
        double umax{0.0};
        // Lookup table resolution
        double bin_size{0.0001};
        // Physics amplitude
        PiPAmplitude pi_p_amplitude;
        // Sampling cache
        std::vector<double> dsigma_domega_costheta;
        double max_cs{0.0};

    public:
        /**
         * @brief Construct a momentum state.
         *
         * This performs all expensive calculations:
         * - phase shift interpolation
         * - amplitude initialization
         * - differential cross section lookup generation
         */
        generatorState(const std::shared_ptr<const PhaseShiftLoader>& phase_shift_loader,const double momentum_lab_,const charge ch,const double theta_min_lab_,const double theta_max_lab_,const double bin_size_=0.0001)
            : momentum_lab(momentum_lab_), theta_lab_min(theta_min_lab_), theta_lab_max(theta_max_lab_), bin_size(bin_size_), 
            pi_p_amplitude(PiPAmplitude(phase_shift_loader, momentum_lab_, ch))
        {
            initialize_kinematics();
            initialize_angles();
            build_lookup_table();
        }
    private:
        /**
         * @brief Calculate momentum-dependent kinematics.
         */
        void initialize_kinematics()
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
        void initialize_angles()
        {
            theta_cm_min = physics_helpers::theta_lab_to_cm(theta_lab_min, gamma_boost, alpha_kinematic);
            theta_cm_max = physics_helpers::theta_lab_to_cm(theta_lab_max, gamma_boost, alpha_kinematic);
            umax = std::cos(physics_helpers::deg2rad(theta_cm_min));
            umin = std::cos(physics_helpers::deg2rad(theta_cm_max));
        }
        /**
         * @brief Build rejection sampling lookup table.
         */
        void build_lookup_table()
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
        }

    public:
        double get_momentum_lab() const noexcept
        {
            return momentum_lab;
        }
        double get_momentum_cm() const noexcept
        {
            return momentum_cm;
        }
        double get_beta_boost() const noexcept
        {
            return beta_boost;
        }
        double get_gamma_boost() const noexcept
        {
            return gamma_boost;
        }
        double get_alpha_kinematic() const noexcept
        {
            return alpha_kinematic;
        }
        double get_theta_lab_min() const noexcept
        {
            return theta_lab_min;
        }
        double get_theta_lab_max() const noexcept
        {
            return theta_lab_max;
        }
        double get_theta_cm_min() const noexcept
        {
            return theta_cm_min;
        }
        double get_theta_cm_max() const noexcept
        {
            return theta_cm_max;
        }
        double get_umin() const noexcept
        {
            return umin;
        }
        double get_umax() const noexcept
        {
            return umax;
        }
        const std::vector<double> &get_dsigma_domega_costheta() const noexcept
        {
            return dsigma_domega_costheta;
        }
        double get_max_cs() const noexcept
        {
            return max_cs;
        }
        const PiPAmplitude &get_amplitude() const noexcept
        {
            return pi_p_amplitude;
        }
        double get_bin_size() const noexcept
        {
            return bin_size;
        }
};

#endif 