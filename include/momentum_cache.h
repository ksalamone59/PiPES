#ifndef MOMENTUM_CACHE_H
#define MOMENTUM_CACHE_H

#include <unordered_map>
#include <memory>
#include <mutex>

#include "generatorState.h"

class MomentumCacheState
{
    public:
        MomentumCacheState(const charge ch, const double theta_min_lab_, const double theta_max_lab_, const double bin_width = 1.0) 
            : bin_width(bin_width), pion_charge(ch), theta_min_lab(theta_min_lab_), theta_max_lab(theta_max_lab_) 
        {
            phase_shift_loader = std::make_shared<PhaseShiftLoader>();
        }
        ~MomentumCacheState() = default;
        std::shared_ptr<const generatorState> get_state(const double lab_momentum)
        {
            const int momentum_idx = static_cast<int>(std::round(lab_momentum / bin_width));
            std::lock_guard<std::mutex> lock(mutex);
            auto it = cross_section_cache.find(momentum_idx);
            if(it != cross_section_cache.end())
            {
                return it->second;
            }
            // Update the momentum based on the cache state
            auto new_state = std::make_shared<const generatorState>(phase_shift_loader, momentum_idx, pion_charge, theta_min_lab, theta_max_lab);
            new_state->get_dsigma_domega_costheta();
            cross_section_cache[momentum_idx] = new_state;
            return new_state;   
        }
        void clear()
        {
            std::lock_guard<std::mutex> lock(mutex);
            cross_section_cache.clear();
        }
        std::size_t size() const
        {
            std::lock_guard<std::mutex> lock(mutex);
            return cross_section_cache.size();
        }
    private:
        std::unordered_map<int, std::shared_ptr<const generatorState>> cross_section_cache; // Example cache, wrt cos theta
        double bin_width{1.};
        mutable std::mutex mutex;
        std::shared_ptr<const PhaseShiftLoader> phase_shift_loader;
        charge pion_charge;
        double theta_min_lab{0.0}, theta_max_lab{0.0};
};

#endif 