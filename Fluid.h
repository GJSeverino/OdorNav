#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

const int spaceN = 100;
const int iter = 16;

inline int IX(int x, int y) {
    return x + y * spaceN;
}




class Fluid {
public:
    int size;
    float dt;
    float diff;
    float visc;
    std::vector<float> s;
    std::vector<float> odor;
    std::vector<float> Vx;
    std::vector<float> Vy;
    std::vector<float> Vx0;
    std::vector<float> Vy0;

    Fluid(float dt, float diffusion, float viscosity);

    void step();
    void addOdor(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);
    float getOdorConcentration(float x, float y);
    void saveodor(std::ofstream& file);



private:
    void set_bnd(int b, std::vector<float>& x);
    void lin_solve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c);
    void diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff, float dt);
    void project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div);
    void advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY, float dt);
};

#endif // FLUID_H
