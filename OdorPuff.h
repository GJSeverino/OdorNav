// OdorPuff.h
#ifndef ODORPUFF_H
#define ODORPUFF_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>

class OdorPuff {
private:
      // Simulation parameters.
    double rate;                           // Poisson release rate 
    double dw_speed;                       // downwind speed
    double init_intensity;                 // Initial intensity of Puff
    double puff_D;                         // puff growth diffusivity 
    double max_x;                          // Maximum distance to track Puff
    double dt;                             // Time step
    double r0;                             // Initial puff radius
    double source_position;                // Position of the odor source
    double decay_rate;                     // Decay rate of odor
    int delay_steps;                       // Number of steps to delay the release of Puff

    // Containers to track the properties of Odor puffs.
    std::vector<double> puff_xs;           // center of puff along x axis 
    std::vector<double> puff_durations;    // duration of puff
    std::vector<double> puff_sizes;        // radius of puff (starts with initial size of r0)

    // Random number generator.
    std::default_random_engine generator;

public:
    OdorPuff(double rate=10, 
             double dw_speed=150, 
             double init_intensity=1074.7, 
             double puff_D=10,
             double max_x=600,
             int delay_steps=240,
             double dt=0.01,
             double decay_rate=0.0);

    void setWindSpeed(double speed);
    void setSource(double source_x);
    void evolve_Puff();
    double compute_concentration(double x);
    void generate_data_file(std::vector<double>& locations, double run_time, const std::string& filename);
};

#endif // ODORPUFF_H
