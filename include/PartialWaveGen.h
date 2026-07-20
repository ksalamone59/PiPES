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

#include "pi_p_amplitude.h"
#include "momentum_cache.h"
#include "io.h"
#include "generatorState.h"

/**
 * @brief Struct that holds all information for a given elastic scattering event
 */
struct elasticEvent
{
    double thetaCM{0.0};
    double phiCM{0.0};

    double thetaLab{0.0};
    double phiLab{0.0};

    fourVector pionLab;
    fourVector protonLab;

    fourVector pionCM;
    fourVector protonCM;
};

/**
 * @brief Partial wave event generator.
 *
 * This class only handles:
 * - random sampling
 * - beam direction
 * - frame transformations
 * - event construction
 *
 * All momentum-dependent physics lives in generatorState.
 */
class PartialWaveGen
{
    private:
        bool verbose{false};
        charge charge_polarity;
        int L_MAX{1};
        // Random number generation
        std::mt19937 mersenne_twister;
        std::uniform_real_distribution<double> uniform_dist;
        std::uniform_real_distribution<double> costheta_dist;
        std::mt19937::result_type seed{987654321};
        // Lookup resolution
        double bin_size{0.0001};
        // Output frame selection
        bool lab_frame{false};
        // Beam orientation
        threeVector beam_direction{0.0,0.0,1.0};
        /**
         * @brief Momentum-dependent physics state.
         *
         * This will later be replaced/looked up from a momentum cache.
         */
        std::shared_ptr<const generatorState> state;
    public:
        PartialWaveGen() = default;
        ~PartialWaveGen() = default;
        /**
         * @brief Construct generator.
         *
         * @param momentum_lab_ Beam momentum in MeV
         * @param ch Pion charge
         * @param theta_min_ Minimum lab scattering angle
         * @param theta_max_ Maximum lab scattering angle
         * @param L_MAX_ Maximum partial wave
         */
        PartialWaveGen(const charge ch, const int L_MAX_=1, bool verbose_= false, std::mt19937::result_type seed_= 987654321) : verbose(verbose_), charge_polarity(ch), L_MAX(L_MAX_), seed(seed_)
        {
            if(charge_polarity != charge::plus && charge_polarity != charge::minus)
            {
                throw std::invalid_argument(
                    "Invalid pion charge");
            }
            set_seed(seed_);
            uniform_dist = std::uniform_real_distribution<double>(0.0, 1.0);
        }
        /**
         * @brief Replace the current state.
         *
         * This is the function that will be used
         * by the momentum cache later.
         */
        void set_state(const std::shared_ptr<const generatorState>& new_state)
        {
            state = new_state;
            costheta_dist = std::uniform_real_distribution<double>(state->get_umin(), state->get_umax());
            set_beam_direction(beam_direction);
            this->bin_size = state->get_bin_size();
        }
        /**
         * @brief Get current generator state.
         */
        const generatorState &get_state() const noexcept
        {
            return *state;
        }
        /**
         * @brief Returns PiPAmplitude.
         */
        const PiPAmplitude &get_pi_p_amplitude() const noexcept
        {
            return state->get_amplitude();
        }

        /**
         * @brief Set RNG seed.
         */
        void set_seed(const std::mt19937::result_type seed_)noexcept
        {
            seed = seed_;
            mersenne_twister.seed(seed);
        }
        /**
         * @brief Set output in lab frame.
         */
        void set_lab_frame()
        {
            lab_frame = true;
        }
        /**
         * @brief Set output in CM frame.
         */
        void set_cm_frame()
        {
            lab_frame = false;
        }
        /**
         * @brief Check output frame.
         */
        bool is_in_lab_frame() const noexcept
        {
            return lab_frame;
        }
        /**
         * @brief Minimum angle in current frame.
         */
        double get_min_theta() const noexcept
        {
            return lab_frame ? state->get_theta_lab_min() : state->get_theta_cm_min();
        }
        /**
         * @brief Maximum angle in current frame.
         */
        double get_max_theta() const noexcept
        {
            return lab_frame ? state->get_theta_lab_max() : state->get_theta_cm_max();
        }
        /**
         * @brief Minimum CM angle.
         */
        double get_min_theta_cm() const noexcept
        {
            return state->get_theta_cm_min();
        }
        /**
         * @brief Maximum CM angle.
         */
        double get_max_theta_cm() const noexcept
        {
            return state->get_theta_cm_max();
        }
        /**
         * @brief Sample uniform azimuthal angle.
         */
        double sample_phi() noexcept
        {
            return 2.0 * std::numbers::pi * uniform_dist(mersenne_twister);
        }
        /**
         * @brief Samples theta according to differential cross section.
         *
         * Uses rejection sampling from the cached lookup table.
         */
        double sample_theta() noexcept
        {
            const auto& lookup = state->get_dsigma_domega_costheta();
            const double umin = state->get_umin();
            while(true)
            {
                double costheta_sample = costheta_dist(mersenne_twister);
                double offset = (costheta_sample - umin) / bin_size;
                std::size_t idx_lo = static_cast<std::size_t>(std::floor(offset + 1e-9));
                std::size_t idx_hi = idx_lo + 1;
                if(idx_hi >= lookup.size()) idx_hi = lookup.size()-1;
                double costheta_lo = umin + idx_lo * bin_size;
                double frac = (costheta_sample - costheta_lo) / bin_size;
                double cs = lookup[idx_lo] + frac * (lookup[idx_hi]-lookup[idx_lo]);
                double probability = std::min(1.0, cs / state->get_max_cs());
                if(uniform_dist(mersenne_twister) < probability)
                {
                    return std::acos(costheta_sample);
                }
            }
        }
        /**
         * @brief Construct pion CM four-vector.
         */
        fourVector get_pion_four_vector(const double theta, const double phi) const
        {
            const double k = state->get_momentum_cm();
            const double px = k * std::sin(theta) * std::cos(phi);
            const double py = k * std::sin(theta) * std::sin(phi);
            const double pz = k * std::cos(theta);
            const double E = std::sqrt(k*k + physics_helpers::m_pion_squared);
            return fourVector(E,px,py,pz);
        }
        /**
         * @brief Construct recoil proton CM four-vector.
         */
        fourVector get_proton_four_vector(const double theta, const double phi) const
        {
            const double k = state->get_momentum_cm();
            const double px = -k * std::sin(theta) * std::cos(phi);
            const double py = -k * std::sin(theta) * std::sin(phi);
            const double pz = -k * std::cos(theta);
            const double E = std::sqrt(k*k +physics_helpers::m_proton_squared);
            return fourVector(E,px,py,pz);
        }
        /**
         * @brief Generate complete elastic scattering event.
         */
        elasticEvent sample_event()
        {
            elasticEvent event;
            event.thetaCM = sample_theta();
            event.phiCM = sample_phi();
            event.pionCM = get_pion_four_vector(event.thetaCM, event.phiCM);
            event.protonCM = get_proton_four_vector(event.thetaCM, event.phiCM);
            threeVector beta_cm_vector = {0., 0., state->get_beta_boost()};
            event.pionLab = physics_helpers::boost_cm_to_lab(event.pionCM, beta_cm_vector);
            event.protonLab = physics_helpers::boost_cm_to_lab(event.protonCM, beta_cm_vector);
            physics_helpers::rotate_four_vector(event.pionLab, beam_direction);
            physics_helpers::rotate_four_vector(event.protonLab, beam_direction);
            std::tie(event.thetaLab,event.phiLab) = physics_helpers::spherical_angles_relative_to_axis(
                                    {event.pionLab.px,event.pionLab.py,event.pionLab.pz},beam_direction);
            return event;
        }
        /**
         * @brief Change incoming beam direction.
         */
        void set_beam_direction(const threeVector& direction)
        {
            beam_direction = direction.normalize();
        }
};

#endif