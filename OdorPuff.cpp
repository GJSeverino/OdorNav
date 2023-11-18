//------------------------------------------------------------------------
//                       1D Odor Puff Simulation 
//
//                    Gabriel J. Severino     |     2023
// 
//
//
//--------------------------------------------------------------------------

#include "OdorPuff.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>

    OdorPuff::OdorPuff(double rate, 
                   double dw_speed, 
                   double init_intensity, 
                   double puff_D,
                   double max_x,
                   int delay_steps,
                   double dt,
                   double decay_rate)
: rate(rate), dw_speed(dw_speed), init_intensity(init_intensity), puff_D(puff_D),
  max_x(max_x), dt(dt), r0(4), delay_steps(delay_steps), decay_rate(decay_rate) {}


    // Set the wind speed.

    void OdorPuff::setWindSpeed(double speed) {
        dw_speed = speed;
    }

    // Set the position of the odor source and initialize the Puff.
    
    void OdorPuff::setSource(double source_x) {
        source_position = source_x;

        // Clear any existing Puff.
        puff_xs.clear();
        puff_durations.clear();
        puff_sizes.clear();

        puff_xs.push_back(source_position);
        puff_durations.push_back(0);
        puff_sizes.push_back(r0);

        // // Generate initial Puff.
        // for(int i = 0; i < delay_steps; i++) {
        //     evolve_Puff();
        // }
    }

    // Update the position and properties of Puff.
    void OdorPuff::evolve_Puff() {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double rand_unif_1 = distribution(generator);
        double prob = 1 - exp(-rate * dt);

        // Based on a probability, create a new puff at the source position.
        if(rand_unif_1 < prob) {
            puff_xs.push_back(source_position); 
            puff_durations.push_back(0);
            puff_sizes.push_back(r0);
        }

        // Update each puff's position and properties.
        for(size_t i = 0; i < puff_xs.size(); i++) {
            puff_xs[i] += dt * dw_speed;
            puff_durations[i] += dt;
            puff_sizes[i] = sqrt(r0*r0 + 4 * puff_D * puff_durations[i]);
        }

        // Remove Puff that have traveled beyond the specified limit.
        for(size_t i = 0; i < puff_xs.size(); ) {
            if(puff_xs[i] > max_x) {
                puff_xs.erase(puff_xs.begin() + i);
                puff_durations.erase(puff_durations.begin() + i);
                puff_sizes.erase(puff_sizes.begin() + i);
            } else {
                i++;
            }
        }
        // Update the initial intensity based on the decay rate and the duration of the Puff.
        for(size_t i = 0; i < puff_xs.size(); i++) {
            init_intensity *= exp(-decay_rate * puff_durations[i]);
        }
    }

  // Calculate the odor concentration at a given position.
    double OdorPuff::compute_concentration(double x) {
        double signal = 0;

        // Sum the contributions from each puff to the overall concentration.
        for(size_t i = 0; i < puff_xs.size(); i++) {
            double distance = fabs(x - puff_xs[i]);
            double scaled_distance = distance / puff_sizes[i];
            double gaussian_part = exp(-scaled_distance * scaled_distance);
            double puff_prefactor = init_intensity / (M_PI * puff_sizes[i] * puff_sizes[i]);
            signal += gaussian_part * puff_prefactor;
        }

        return signal;
    }



    // Generate a data file with odor concentrations at specified locations over time.
    void OdorPuff::generate_data_file(std::vector<double>& locations, double run_time, const std::string& filename) {
        int num_steps = static_cast<int>(round(run_time / dt));
        std::ofstream out(filename);

        for(int i = 0; i < num_steps; i++) {
            evolve_Puff();
            for(double loc : locations) {
                double concentration = compute_concentration(loc);
                out << concentration << " ";
            }
            out << "\n";
        }

        out.close();
    }





// int main() {


//     OdorPuff simulation;

//     // Define the source position and set it in the simulation.
//     double source_position = 200;
//     simulation.setSource(source_position);
//     simulation.setWindSpeed(-150);
    
//     // Specify the locations for which we want to compute concentrations.
//     std::vector<double> locations;
//     for(int i = 0; i <= 600; i++) {
//         locations.push_back(i);
//     }

//     simulation.generate_data_file(locations, 60, "output.dat");



//     return 0;
// }
