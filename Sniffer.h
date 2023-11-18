#pragma once

#include "CTRNN.h"

// The Sniffer Agent class declaration
class Sniffer {
public:
    // The constructor
    Sniffer(int networksize);

    // The destructor
   	~Sniffer();

    // Accessors
    double GetPosX() const { return posX; }
    double GetPosY() const { return posY; }
    double GetPastPosX() const { return pastposX; }
    double GetPastPosY() const { return pastposY; }
    double GetVelocity() const { return velocity; }
    double GetTheta() const { return theta; }
    void SetPosX(double newPosX) { posX = newPosX; }
    void SetPosY(double newPosY) { posY = newPosY; }
    void SetSensorWeight(int index, double value) { sensorweights[index] = value; }

    // Control methods
    void Set(int networksize);
    void Reset(double initposX, double initposY);
    double MapBreathingRate(double neuronOutput);
    void Sense(double chemical_concentration, double current_time);
    void Step(double StepSize);

    // Properties
    int size;
    double posX, posY, pastposX, pastposY, velocity, theta, pastTheta, gain, sensor;
    TVector<double> sensorweights;
    CTRNN NervousSystem;
};
