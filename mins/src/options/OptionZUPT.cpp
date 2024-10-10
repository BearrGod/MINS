#include "OptionZUPT.h"
#include "utils/Print_Logger.h"
#include "utils/opencv_yaml_parse.h"


void mins::OptionsZeroVelocityUpdate::load(const std::shared_ptr<ov_core::YamlParser>& parser){
    if(parser!=nullptr){
        std::string f = "config_estimator";
        parser->parse_external(f,"est","zupt_chi2_multipler", zupt_chi2_multipler);
        parser->parse_external(f,"est","zupt_max_velocity", zupt_max_velocity);
        parser->parse_external(f,"est","zupt_noise_multiplier", zupt_noise_multiplier);
        parser->parse_external(f,"est","zupt_max_disparity", zupt_max_disparity);
        parser->parse_external(f,"est","zupt_only_at_beginning", zupt_only_at_beginning);
        parser->parse_external(f,"est","try_zupt", try_zupt);
        parser->parse_external(f,"est","sigma_pix", sigma_pix,false);
    }
}

void mins::OptionsZeroVelocityUpdate::print() {
  PRINT2(BOLDBLUE "Options - ZeroVelocityUpdate\n" RESET);
  PRINT2("\t- try_zupt: %s\n", try_zupt ? "true" : "false");
  PRINT2("\t- zupt_chi2_multipler: %.9f\n", zupt_chi2_multipler);
  PRINT2("\t- zupt_max_velocity: %.9f\n", zupt_max_velocity);
  PRINT2("\t- zupt_noise_multiplier: %.9f\n", zupt_noise_multiplier);
  PRINT2("\t- zupt_only_at_beginning: %s\n", zupt_only_at_beginning ? "true" : "false");
  PRINT2("\t- zupt_max_disparity: %.9f\n", zupt_max_disparity);
  PRINT2("\t- sigma_pix: %.1f\n", sigma_pix);

}