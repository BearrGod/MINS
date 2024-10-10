#ifndef MINS_OPTIONS_ZUPT_H
#define MINS_OPTIONS_ZUPT_H
#include <map>
#include <memory>
#include <vector>

namespace ov_core {
class YamlParser;
struct FeatureInitializerOptions;
} // namespace ov_core

namespace mins {
    struct OptionsZeroVelocityUpdate {
        void load(const std::shared_ptr<ov_core::YamlParser>& parser = nullptr) ; 
         
        void print() ;  

        double zupt_chi2_multipler = 0 ; 
        double zupt_max_velocity = 0.1 ; 
        double zupt_noise_multiplier = 10 ; 
        double zupt_max_disparity = 0.5 ; 
        bool zupt_only_at_beginning = false ; 
        bool try_zupt = false ; 
        double sigma_pix = 1;
        double sigma_pix_sq = 1;
    } ; 
} // namespace mins
#endif // MINS_OPTIONS_ZUPT_H 