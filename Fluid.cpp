#include <vector>
#include <cmath>
#include <algorithm> // For std::fill
#include <iostream>
#include <fstream>
#include "Fluid.h"



Fluid::Fluid(float dt, float diffusion, float viscosity)
    : size(spaceN), dt(dt), diff(diffusion), visc(viscosity),
      s(spaceN * spaceN, 0), odor(spaceN * spaceN, 0),
      Vx(spaceN * spaceN, 0), Vy(spaceN * spaceN, 0),
      Vx0(spaceN * spaceN, 0), Vy0(spaceN * spaceN, 0) {}



void Fluid::saveodor(std::ofstream& file) {
    for (int j = 0; j < spaceN; j++) {
        for (int i = 0; i < spaceN; i++) {
            file << odor[IX(i, j)] << " ";
        }
        file << "\n";
    }
    file << "\n\n"; // Delimiter between frames
}

float Fluid::getOdorConcentration(float x, float y) {
    // Check bounds and return -1.0f if out of range
    if (x < 0.0f || x > static_cast<float>(spaceN) || y < 0.0f || y > static_cast<float>(spaceN)) {
        std::cerr << "Coordinates out of bounds!" << std::endl;
        return -1.0f;
    }

    // Get the integer part of the coordinates
    int x_int = static_cast<int>(x);
    int y_int = static_cast<int>(y);

    // Get the fractional part of the coordinates
    float x_frac = x - static_cast<float>(x_int);
    float y_frac = y - static_cast<float>(y_int);

    // Check bounds for top-right corner of the interpolation square
    if (x_int >= spaceN - 1 || y_int >= spaceN - 1) {
        return odor[IX(x_int, y_int)]; // Return the value at the bottom-left corner
    }

    // Perform bilinear interpolation
    float bl = odor[IX(x_int, y_int)];     // Bottom Left
    float br = odor[IX(x_int + 1, y_int)]; // Bottom Right
    float tl = odor[IX(x_int, y_int + 1)]; // Top Left
    float tr = odor[IX(x_int + 1, y_int + 1)]; // Top Right

    float b = bl + x_frac * (br - bl); // Interpolate along the bottom edge
    float t = tl + x_frac * (tr - tl); // Interpolate along the top edge
    return b + y_frac * (t - b); // Interpolate between top and bottom
}

void Fluid::set_bnd(int b, std::vector<float>& x) {
    for (int i = 1; i < spaceN - 1; i++) {
        x[IX(i, 0)] = (b == 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, spaceN - 1)] = (b == 2) ? -x[IX(i, spaceN - 2)] : x[IX(i, spaceN - 2)];
    }

    for (int j = 1; j < spaceN - 1; j++) {
        x[IX(0, j)] = (b == 1) ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(spaceN - 1, j)] = (b == 1) ? -x[IX(spaceN - 2, j)] : x[IX(spaceN - 2, j)];
    }

    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, spaceN - 1)] = 0.5 * (x[IX(1, spaceN - 1)] + x[IX(0, spaceN - 2)]);
    x[IX(spaceN - 1, 0)] = 0.5 * (x[IX(spaceN - 2, 0)] + x[IX(spaceN - 1, 1)]);
    x[IX(spaceN - 1, spaceN - 1)] = 0.5 * (x[IX(spaceN - 2, spaceN - 1)] + x[IX(spaceN - 1, spaceN - 2)]);
}

void Fluid::lin_solve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < spaceN - 1; j++) {
            for (int i = 1; i < spaceN - 1; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)] +
                     a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] +
                          x[IX(i, j + 1)] + x[IX(i, j - 1)])) *
                    cRecip;
            }
        }
        set_bnd(b, x);
    }
}

void Fluid::diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff, float dt) {
    float a = dt * diff * (spaceN - 2) * (spaceN - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}


void Fluid::project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div) {
    for (int j = 1; j < spaceN - 1; j++) {
        for (int i = 1; i < spaceN - 1; i++) {
            div[IX(i, j)] = (-0.5 * (velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] +
                                     velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)])) / spaceN;
            p[IX(i, j)] = 0;
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for (int j = 1; j < spaceN - 1; j++) {
        for (int i = 1; i < spaceN - 1; i++) {
            velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * spaceN;
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * spaceN;
        }
    }

    set_bnd(1, velocX);
    set_bnd(2, velocY);
}


void Fluid::advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY, float dt) {
    float i0, i1, j0, j1;

    float dtx = dt * (spaceN - 2);
    float dty = dt * (spaceN - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    for (int j = 1; j < spaceN - 1; j++) {
        for (int i = 1; i < spaceN - 1; i++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x = i - tmp1;
            y = j - tmp2;

            if (x < 0.5) x = 0.5;
            if (x > spaceN + 0.5) x = spaceN + 0.5;
            i0 = floor(x);
            i1 = i0 + 1.0;

            if (y < 0.5) y = 0.5;
            if (y > spaceN + 0.5) y = spaceN + 0.5;
            j0 = floor(y);
            j1 = j0 + 1.0;

            s1 = x - i0;
            s0 = 1.0 - s1;
            t1 = y - j0;
            t0 = 1.0 - t1;

            d[IX(i, j)] =
                s0 * (t0 * d0[IX((int)i0, (int)j0)] + t1 * d0[IX((int)i0, (int)j1)]) +
                s1 * (t0 * d0[IX((int)i1, (int)j0)] + t1 * d0[IX((int)i1, (int)j1)]);
        }
    }

    set_bnd(b, d);
}

void Fluid::addOdor(int x, int y, float amount) {
    int index = IX(x, y);
    this->odor[index] += amount;
}

void Fluid::addVelocity(int x, int y, float amountX, float amountY) {
    int index = IX(x, y);
    this->Vx[index] += amountX;
    this->Vy[index] += amountY;
}




void Fluid::step() {

    // Diffuse velocity
    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);

    // Project velocity
    project(Vx0, Vy0, Vx, Vy);

    // Advect velocity
    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);

    // Project velocity
    project(Vx, Vy, Vx0, Vy0);

    // Diffuse and advect odor (odor concentration)
    diffuse(0, s, odor, diff, dt);
    advect(0, odor, s, Vx, Vy, dt);
}

// int main() {
//     std::ofstream file("odor_data.dat");

//     Fluid fluid(0.01, 0.0, 0.0000001);

//     int sourceX = 50; // Example x-coordinate for odor source
//     int sourceY = 50;  // Example y-coordinate for odor source
//     float odorAmount = 10; // Amount of odor to add
//     float velocityX = 2; // Velocity in the x-direction
//     float velocityY = 2; // Velocity in the y-direction

//     for (int i = 0; i < 1000; i++) {
//         fluid.addOdor(sourceX, sourceY, odorAmount);
//         fluid.addVelocity(sourceX, sourceY, velocityX, velocityY);
//         fluid.step();
//         fluid.saveodor(file); // Save odor data to file

//         // float concentration = fluid.getOdorConcentration(50, 50);
//         // std::cout << "Odor concentration at (50,50): " << concentration << std::endl;
//     }

//     file.close();
//     return 0;
// }
