
#include "Sniffer.h"
#include <cmath>  // for sine and exponential functions
#include "random.h"

// Constants
const double SpaceHeight = 100.0;
const double SpaceWidth = 100.0;

// Constants for motion
const double MaxAngle = M_PI / 12.0; // Pi/12
const double MaxThrust = 1.5;
const double Friction = 0.9;

// Contructor
Sniffer::Sniffer(int networksize) {
    Set(networksize);
}

// destructor
Sniffer::~Sniffer() {
    Set(0);
}


// Initialize the agent
void Sniffer::Set(int networksize) {
    size = networksize;
    // gain = 2.0;
    sensorweights.SetBounds(1, size);
    sensorweights.FillContents(0.0);
    posX = 0.0;
    posY = 0.0;
    pastposX = 0.0;
    pastposY = 0.0;
    velocity = 0.0;
    pastTheta = 0.0;
    theta = 0.0;

    sensor = 0.0;
    
}

// Reset the state of the agent
void Sniffer::Reset(double initposX, double initposY, double initTheta) {
    posX = initposX;
    posY = initposY;
    pastposX = initposX;
    pastposY = initposY;
    sensor = 0.0;
    velocity = 0.0;
    theta = initTheta;
    NervousSystem.RandomizeCircuitState(0.0, 0.0);
}

// Map the output of neuron 3 to the breathing rate
double Sniffer::MapBreathingRate(double neuronOutput){
    double minRate = 0.1;
    double maxRate = 2.0;

    return minRate + (maxRate - minRate) * neuronOutput;

}

// // Respiration sense function
// void Sniffer::Sense(double chemical_concentration, double current_time) {
    
//     // Use the output of neuron 3 to control the breathing rate
//     double breathingRate = MapBreathingRate(NervousSystem.NeuronOutput(3));

//     // Phase of the sin wave determines inhale/exhale cycle 
//     double phase = sin(current_time * 2 * M_PI * breathingRate);
    
//     // Agent perceives the chemical concentration only inhalation (positive phase of the sine wave)
//     if (phase > 0) {
//         sensor = chemical_concentration * phase; // * phase for a more continuous signal. May need to multiply by derivative of phase. 
//     } else {
//         sensor = 0.0;
//     }
// }


// Normal sense function
void Sniffer::Sense(double chemical_concentration, double current_time){
    if (current_time > 0.0) {

        
        sensor = chemical_concentration;
    } else {
        sensor = 0.0;
    }
}

// Step in time
void Sniffer::Step(double StepSize) {
    
    pastposX = posX;
    pastposY = posY;
    pastTheta = theta;

    for (int i = 1; i <= size; i++) {
        NervousSystem.SetNeuronExternalInput(i, sensor * sensorweights[i]);
    }

	// Update the nervous system
    NervousSystem.EulerStep(StepSize);

        // Constants for motion
    const double MaxAngle = M_PI / 12.0; // Pi/12
    const double MaxThrust = 1.5;
    const double Friction = 0.9;

    double outputMotorRight = NervousSystem.NeuronOutput(1); 
    double outputMotorLeft = NervousSystem.NeuronOutput(2); 

    // Calculate the torque and thrust based on the neural outputs
    double torque = (outputMotorRight - outputMotorLeft) * MaxAngle;
    double thrust = (outputMotorRight + outputMotorLeft) * MaxThrust;

    // Update velocity and angle
    velocity = velocity * Friction + StepSize * thrust;
    theta += StepSize * torque;

    // Calculate the new position based on velocity and angle
    posX += StepSize * velocity * cos(theta);
    posY += StepSize * velocity * sin(theta);

    // // Wrap around environment 
    // if (posX >= SpaceWidth) posX -= SpaceWidth;
    // if (posX < 0.0) posX += SpaceWidth;
    // if (posY >= SpaceHeight) posY -= SpaceHeight;
    // if (posY < 0.0) posY += SpaceHeight;

    // Zero boundary conditions
    if (posX >= SpaceWidth) posX = SpaceWidth;
    if (posX < 0.0) posX = 0.0;
    if (posY >= SpaceHeight) posY = SpaceHeight;
    if (posY < 0.0) posY = 0.0;


}




