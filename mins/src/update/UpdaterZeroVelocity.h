/*
 * OpenVINS: An Open Platform for Visual-Inertial Research
 * Copyright (C) 2018-2023 Patrick Geneva
 * Copyright (C) 2018-2023 Guoquan Huang
 * Copyright (C) 2018-2023 OpenVINS Contributors
 * Copyright (C) 2018-2019 Kevin Eckenhoff
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef MINS_UPDATER_ZEROVELOCITY_H
#define MINS_UPDATER_ZEROVELOCITY_H

#include <memory> 
#include "options/OptionsEstimator.h"
#include "state/Propagator.h"
#include "state/State.h"
#include "utils/sensor_data.h"
#include "options/OptionsIMU.h"
#include "options/OptionZUPT.h"




namespace mins {

class UpdaterStatistics;
struct OptionsZeroVelocityUpdate;
struct OptionsIMU ; // Imu noises here

std::vector<ov_core::ImuData> select_imu_readings(const std::vector<ov_core::ImuData> &imu_data, double time0, double time1,
                                                           bool warn = true)  ; 

class UpdaterZeroVelocity {
    public:
        UpdaterZeroVelocity(){} ; // Default constructur made explicit here
        UpdaterZeroVelocity(OptionsZeroVelocityUpdate& options , OptionsIMU& noises 
                            ,std::shared_ptr<Propagator> prop ,double gravity_mag) ;

        /**
        * @brief Feed function for inertial data
        * @param message Contains our timestamp and inertial information
        * @param oldest_time Time that we can discard measurements before
        */
        void feed_imu(const ov_core::ImuData &message, double oldest_time = -1) ;

        /**
         * @brief This will remove any IMU measurements that are older then the given measurement time
         * @param oldest_time Time that we can discard measurements before (in IMU clock)
         */
        void clean_old_imu_measurements(double oldest_time) ;

         /**
         * @brief Will first detect if the system is zero velocity, then will update.
         * @param state State of the filter
         * @param timestamp Next camera timestamp we want to see if we should propagate to.
         * @return True if the system is currently at zero velocity
         */
        bool try_update(std::shared_ptr<State> state, double timestamp);

    protected : 
        /// Our history of IMU messages (time, angular, linear)
        std::vector<ov_core::ImuData> imu_data ;
        
        /// Options used during update (chi2 multiplier)
        OptionsZeroVelocityUpdate _options ;
        
        /// Container for imu noises value 
        OptionsIMU _noises ; 
        
        /// Feature tracker database with all features in it
        // std::shared_ptr<ov_core::FeatureDatabase> _db;

        /// Our propagator!
        std::shared_ptr<Propagator> _prop;

        /// Gravity vector
        Eigen::Vector3d _gravity;

        /// Chi squared 95th percentile table (lookup would be size of residual)
        std::map<int, double> chi_squared_table;

        /// Estimate for time offset at last propagation time
        double last_prop_time_offset = 0.0;
        bool have_last_prop_time_offset = false;

        /// Last timestamp we did zero velocity update with
        double last_zupt_state_timestamp = 0.0;

        /// Number of times we have called update
        int last_zupt_count = 0;

} ; 

} // namespace  mins

#endif // OV_MSCKF_UPDATER_ZEROVELOCITY_H
