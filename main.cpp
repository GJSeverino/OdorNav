// *******************************************************************************************
//       Gabriel Juliano Severino 
//               2023
//
//       Active Respiratory Agent 
// 
//
// TO DO:
// 
// 
// *******************************************************************************************

#include "OdorPuff.h"
#include "TSearch.h"
#include "Sniffer.h"
#include "CTRNN.h"
#include "random.h"
#include <random>

#define PRINTOFILE

// Task params
const double StepSize = 0.01;
const double RunDuration = 5000; // 800.0;
const double TransDuration = 3000.0; // Transient duration 
const double EvalDuration = RunDuration - TransDuration; // Evaluation duration

const double SpaceHeight = 100.0; // Size of the space
const double SpaceWidth = 100.0; 

const double peakPositionX = 50.0; // Position of the chemical source
const double peakPositionY = 50.0; 

// EA params
const int POPSIZE = 300; //500 
const int GENS = 100; //100
const double MUTVAR = 0.15;
const double CROSSPROB = 0.05;
const double EXPECTED = 1.1;
const double ELITISM = 0.02;

// Nervous system params
const int N = 3;
const double WR = 8.0;
const double SR = 8.0;
const double BR = 8.0;
const double TMIN = 1.0;
const double TMAX = 10.0;

int	VectSize = N*N + 2*N + N;

// ------------------------------------
// Genotype-Phenotype Mapping Functions
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	int k = 1;
	// Time-constants
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), TMIN, TMAX);
		k++;
	}
	// Bias
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), -BR, BR);
		k++;
	}
	// Weights
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			phen(k) = MapSearchParameter(gen(k), -WR, WR);
			k++;
		}
	}
	// Sensor Weights
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), -SR, SR);
		k++;
	}
}

double DistanceGradient(double posX, double posY, double peakPosX, double peakPosY, double steepness = 1.5) {
    // Calculate direct Euclidean distance along x and y axes
    double dx = std::abs(posX - peakPosX);
    double dy = std::abs(posY - peakPosY);

    // Euclidean distance in a 2D space
    double effective_distance = sqrt(dx * dx + dy * dy);

    // Normalize the distance to the maximum possible distance
    const double max_distance = sqrt((SpaceWidth * SpaceWidth) + (SpaceHeight * SpaceHeight));
    double normalized_distance = effective_distance / max_distance;

    // Calculate the gradient concentration using the normalized distance and steepness
    return 1.0 - std::abs(normalized_distance) * steepness;
}

// ------------------------------------
// Fitness function 
// ------------------------------------
// double FitnessFunction(TVector<double> &genotype, RandomState &rs)
// {
// 	// Map genotype to phenotype
// 	TVector<double> phenotype;
// 	phenotype.SetBounds(1, VectSize);
// 	GenPhenMapping(genotype, phenotype);

// 	// Create the agent
// 	Sniffer Agent(N);

// 	// Instantiate the nervous systems
// 	Agent.NervousSystem.SetCircuitSize(N);
	
// 	int k = 1;

// 	// Time-constants
// 	for (int i = 1; i <= N; i++) {
// 		Agent.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
// 		k++;
// 	}
// 	// Biases
// 	for (int i = 1; i <= N; i++) {
// 		Agent.NervousSystem.SetNeuronBias(i,phenotype(k));
// 		k++;
// 	}
// 	// Weights
// 	for (int i = 1; i <= N; i++) {
// 		for (int j = 1; j <= N; j++) {
// 			Agent.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
// 			k++;
// 		}
// 	}
// 	// Sensor Weights
// 	for (int i = 1; i <= N; i++) {
// 		Agent.SetSensorWeight(i,phenotype(k));
// 		k++;
// 	}

//     double totalfit = 0.0, totaldist = 0.0, dist = 0.0;
//     double totaltime = 0, totaltrials = 0;

//     // Vary the steepness of the gradient
//     const double minSteepness = 0.5;  
//     const double maxSteepness = 2.0;  
//     const double steepnessStep = 0.5;

// 		for (double steepness = minSteepness; steepness <= maxSteepness; steepness += steepnessStep) {
// 			// for (double x = 1.0; x < SpaceWidth; x += 50) {  // Trials across starting positions
// 			// 	for(double y = 1.0; y < SpaceHeight; y += 50) {
// 				double x = rs.UniformRandom(0.0, SpaceWidth);
// 				double y = rs.UniformRandom(0.0, SpaceHeight);

// 				// Peak position of chemical gradient
// 				const double peakPositionX = rs.UniformRandom(0.0, SpaceWidth);
// 				const double peakPositionY = rs.UniformRandom(0.0, SpaceHeight);

// 				// Set agent's position
// 				Agent.Reset(x,y);

// 				totaldist = 0.0;
// 				totaltime = 0;

// 				for (double time = 0; time < RunDuration; time += StepSize) {

// 					// Calculate chemical gradient at Sniffer position
// 					double gradientValue = DistanceGradient(Agent.posX, Agent.posY, peakPositionX, peakPositionY, steepness);

// 					// Sense the gradient
// 					Agent.Sense(gradientValue, time);
					
// 					// Move based on sensed gradient
// 					Agent.Step(StepSize);

// 					if (time > TransDuration) {
   						 
// 						 // Calculate the distance to the peak considering the 2D environment and wrap-around boundaries
//    					 	double dx = std::abs(Agent.posX - peakPositionX);
//     					double dy = std::abs(Agent.posY - peakPositionY);
//     					double directDist = sqrt(dx * dx + dy * dy);

//     					// Ensure that the minimum distance is at least 2.0 to avoid division by zero or negative fitness
//     					if (directDist < 2.0)
//         					directDist = 2.0;

//     					totaldist += directDist;
//     					totaltime += 1;
// 					}
// 				}

// 					double maxDist = sqrt((SpaceWidth * SpaceWidth) + (SpaceHeight * SpaceHeight)); 
// 					double averageDist = totaldist / totaltime;
// 					double fitnessForThisTrial = 0.0;
// 					double thresholdDistance = 4.0; // Define a threshold distance

// 					if (averageDist < thresholdDistance) {
// 						// Provide a high fitness score if within threshold distance
// 						fitnessForThisTrial = 1.0 - exp(-averageDist);
// 					} else {
// 						// Penalize heavily if not within threshold distance
// 						fitnessForThisTrial = exp(-averageDist / maxDist);
// 					}
// 					totalfit += fitnessForThisTrial;
// 					totaltrials += 1;
// 				}

//     // Average fitness across trials
//     return totalfit / totaltrials;
// }


double FitnessFunctionChemoIndex(TVector<double> &genotype, RandomState &rs)
{
	// Map genotype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);

	// Create the agent
	Sniffer Agent(N);

	// Instantiate the nervous systems
	Agent.NervousSystem.SetCircuitSize(N);
	
	int k = 1;

	// Time-constants
	for (int i = 1; i <= N; i++) {
		Agent.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
		k++;
	}
	// Biases
	for (int i = 1; i <= N; i++) {
		Agent.NervousSystem.SetNeuronBias(i,phenotype(k));
		k++;
	}
	// Weights
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			Agent.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
			k++;
		}
	}
	// Sensor Weights
	for (int i = 1; i <= N; i++) {
		Agent.SetSensorWeight(i,phenotype(k));
		k++;
	}

    double totalFit = 0.0;
    int trials = 0;

    // Vary the steepness of the gradient
    const double minSteepness = 0.5;
    const double maxSteepness = 2.0;
    const double steepnessStep = 0.5;

    for (double steepness = minSteepness; steepness <= maxSteepness; steepness += steepnessStep) {
        for (int reps = 1; reps <= 2; reps++) {  

            double x = rs.UniformRandom(0.0, SpaceWidth);
            double y = rs.UniformRandom(0.0, SpaceHeight);

            // Peak position of chemical gradient
            const double peakPositionX = rs.UniformRandom(0.0, SpaceWidth);
            const double peakPositionY = rs.UniformRandom(0.0, SpaceHeight);

            // Calculate initial distance
            double initialDist = sqrt(pow(x - peakPositionX, 2) + pow(y - peakPositionY, 2));
            if (initialDist < 1.0) initialDist = 1.0; // Avoid division by zero

            // Set agent's position
            Agent.Reset(x, y);

            double totalDist = 0.0;
			int evaltimesteps = 0;

            for (double time = 0; time < RunDuration; time += StepSize) {
                
				// Calculate chemical gradient at Sniffer position
				double gradientValue = DistanceGradient(Agent.posX, Agent.posY, peakPositionX, peakPositionY, steepness);
				// Sense the gradient
				Agent.Sense(gradientValue, time);
					
				// Move based on sensed gradient
				Agent.Step(StepSize);

				if (time > TransDuration) {
					double dx = std::abs(Agent.posX - peakPositionX);
					double dy = std::abs(Agent.posY - peakPositionY);
					totalDist += sqrt(dx * dx + dy * dy);

					evaltimesteps++;

					}
							}

							// Only calculate the average distance if there were evaluation time steps
							double averageDist = evaltimesteps > 0 ? totalDist / evaltimesteps : initialDist;
							double normalizedDist = averageDist / initialDist;
							double fitnessForThisTrial = 1.0 - normalizedDist;
							fitnessForThisTrial = fitnessForThisTrial < 0.0 ? 0.0 : fitnessForThisTrial; // Ensure non-negative fitness

							totalFit += fitnessForThisTrial;
							trials++;
						}
					}
					return totalFit / trials;
				}



//                 if (time > TransDuration) {
//                     double dx = std::abs(Agent.posX - peakPositionX);
//                     double dy = std::abs(Agent.posY - peakPositionY);
//                     totalDist += sqrt(dx * dx + dy * dy);
//                 }
//             }

//             double averageDist = totalDist / (RunDuration / StepSize); // run or eval duration? 
//             double fitnessForThisTrial = 1.0 - (averageDist / initialDist);
//             fitnessForThisTrial = fitnessForThisTrial < 0.0 ? 0.0 : fitnessForThisTrial; // Ensure non-negative fitness

//             totalFit += fitnessForThisTrial;
//             trials++;
//         }
//     }
//     return totalFit / trials;
// }

// ================================================
// FUNCTIONS FOR ANALYZING A SUCCESFUL CIRCUIT
// ================================================

void BehavioralTraces_Specific(TVector<double> &genotype,double x1, double y1, double chemicalsourceX, double chemicalsourceY, double steepness)
{
    // Start output file for positions
    ofstream positionFile("AgentPositions.dat");

    // Map genotype to phenotype
    TVector<double> phenotype;
    phenotype.SetBounds(1, VectSize);
    GenPhenMapping(genotype, phenotype);

    // Create the agent
    Sniffer Agent(N);

    // Instantiate the nervous systems
    Agent.NervousSystem.SetCircuitSize(N);
    int k = 1;

    // Time-constants
    for (int i = 1; i <= N; i++) {
        Agent.NervousSystem.SetNeuronTimeConstant(i, phenotype(k));
        k++;
    }
    // Biases
    for (int i = 1; i <= N; i++) {
        Agent.NervousSystem.SetNeuronBias(i, phenotype(k));
        k++;
    }
    // Weights
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            Agent.NervousSystem.SetConnectionWeight(i, j, phenotype(k));
            k++;
        }
    }
    // Sensor Weights
    for (int i = 1; i <= N; i++) {
        Agent.SetSensorWeight(i, phenotype(k));
        k++;
    }

    // Run the simulation
    // for (double y1 = 0.0; y1 < SpaceHeight; y1 += 50.0) {
    //     for (double x1 = 0.0; x1 < SpaceWidth; x1 += 50.0) {
            
            const double peakPositionX = chemicalsourceX;
            const double peakPositionY = chemicalsourceY;

            // Set agent's position
            Agent.Reset(x1, y1);

            for (double time = 0; time < RunDuration; time += StepSize) {
                
                // Sense the gradient
                double gradientValue = DistanceGradient(Agent.posX, Agent.posY, peakPositionX, peakPositionY, steepness);

                Agent.Sense(gradientValue, time);

                // Move based on sensed gradient
                Agent.Step(StepSize);

                // Save the position at each timestep
                positionFile << Agent.posX << " " << Agent.posY << endl;
            }
    //     }
    // }

    positionFile.close();
}

void BehavioralTraces_Across_Conditions(TVector<double> &genotype, double steepness)
{
    // Start output file for positions
    ofstream positionFile("AgentPositions.dat");

    // Map genotype to phenotype
    TVector<double> phenotype;
    phenotype.SetBounds(1, VectSize);
    GenPhenMapping(genotype, phenotype);

    // Create the agent
    Sniffer Agent(N);

    // Instantiate the nervous systems
    Agent.NervousSystem.SetCircuitSize(N);
    int k = 1;

    // Time-constants
    for (int i = 1; i <= N; i++) {
        Agent.NervousSystem.SetNeuronTimeConstant(i, phenotype(k));
        k++;
    }
    // Biases
    for (int i = 1; i <= N; i++) {
        Agent.NervousSystem.SetNeuronBias(i, phenotype(k));
        k++;
    }
    // Weights
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            Agent.NervousSystem.SetConnectionWeight(i, j, phenotype(k));
            k++;
        }
    }
    // Sensor Weights
    for (int i = 1; i <= N; i++) {
        Agent.SetSensorWeight(i, phenotype(k));
        k++;
    }

	const int NumPeakPositions = 6; // Change this as needed

    // Initialize the agent's position
    double initPosX = 50.0;
    double initPosY = 50.0;
    
    // Array to hold the peak positions
    double peakPositionsX[NumPeakPositions] = { /* x coordinates */ };
    double peakPositionsY[NumPeakPositions] = { /* y coordinates */ };

    
    double angleStep = (2 * M_PI) / NumPeakPositions;
    double radius = 25; // Change the radius as needed


    for (int i = 0; i < NumPeakPositions; ++i) {
        peakPositionsX[i] = initPosX + radius * std::cos(i * angleStep);
        peakPositionsY[i] = initPosY + radius * std::sin(i * angleStep);
    }

    // Iterate over all peak positions
    for (int i = 0; i < NumPeakPositions; ++i) {
        // Reset agent to starting position
        Agent.Reset(initPosX, initPosY);

        // Peak position of chemical gradient
        double peakPositionX = peakPositionsX[i];
        double peakPositionY = peakPositionsY[i];

        for (double time = 0; time < RunDuration; time += StepSize) {
            // Sense the gradient
            double gradientValue = DistanceGradient(Agent.posX, Agent.posY, peakPositionX, peakPositionY, steepness);

            Agent.Sense(gradientValue, time);

            // Move based on sensed gradient
            Agent.Step(StepSize);

            // Save the position at each timestep
            positionFile << Agent.posX << " " << Agent.posY << std::endl;
        }
    }

    positionFile.close();
}

// ================================================
// C. ADDITIONAL EVOLUTIONARY FUNCTIONS
// ================================================
int TerminationFunction(int Generation, double BestPerf, double AvgPerf, double PerfVar) {
	if (BestPerf > 0.99) return 1;
	else return 0;
}

// ------------------------------------
// Display functions
// ------------------------------------
void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
	cout << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}

void ResultsDisplay(TSearch &s)
{
	TVector<double> bestVector;
	ofstream BestIndividualFile;
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);

	// Save the genotype of the best individual
	bestVector = s.BestIndividual();
	BestIndividualFile.open("best.gen.dat");
	BestIndividualFile << bestVector << endl;
	BestIndividualFile.close();

	// Also show the best individual in the Circuit Model form
	BestIndividualFile.open("best.ns.dat");
	GenPhenMapping(bestVector, phenotype);
	Sniffer Agent(N);

	// Instantiate the nervous system
	Agent.NervousSystem.SetCircuitSize(N);
	int k = 1;
	// Time-constants
	for (int i = 1; i <= N; i++) {
		Agent.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
		k++;
	}
	// Bias
	for (int i = 1; i <= N; i++) {
		Agent.NervousSystem.SetNeuronBias(i,phenotype(k));
		k++;
	}
	// Weights
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			Agent.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
			k++;
		}
	}
		// Sensor Weights
	for (int i = 1; i <= N; i++) {
		Agent.SetSensorWeight(i,phenotype(k));
		k++;
	}
	BestIndividualFile << Agent.NervousSystem << endl;
	BestIndividualFile << Agent.sensorweights << "\n" << endl;
	BestIndividualFile.close();
}

// ------------------------------------
// THE MAIN PROGRAM 
// ------------------------------------
int main (int argc, const char* argv[]) 
{

//-------------------------------------
// EVOLUTION 
//-------------------------------------
	
	long randomseed = static_cast<long>(time(NULL));
	if (argc == 2)
		randomseed += atoi(argv[1]);

	TSearch s(VectSize);
	
	#ifdef PRINTOFILE

	ofstream file;
	file.open("evol.dat");
	cout.rdbuf(file.rdbuf());
	
	// save the seed to a file
	ofstream seedfile;
	seedfile.open ("seed.dat");
	seedfile << randomseed << endl;
	seedfile.close();
	
	RandomState rs(randomseed); 

	#endif
	
	// Configure the search
	s.SetRandomSeed(randomseed);
	s.SetSearchResultsDisplayFunction(ResultsDisplay);
	s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
	s.SetSelectionMode(RANK_BASED);
	s.SetReproductionMode(GENETIC_ALGORITHM);
	s.SetPopulationSize(POPSIZE);
	s.SetMaxGenerations(GENS);
	s.SetCrossoverProbability(CROSSPROB);
	s.SetCrossoverMode(UNIFORM);
	s.SetMutationVariance(MUTVAR);
	s.SetMaxExpectedOffspring(EXPECTED);
	s.SetElitistFraction(ELITISM);
	s.SetSearchConstraint(1);
	
	/* Stage 1 */
	s.SetSearchTerminationFunction(TerminationFunction);
	s.SetEvaluationFunction(FitnessFunctionChemoIndex); 

	s.ExecuteSearch();


// ================================================
// B. MAIN FOR ANALYZING A SUCCESFUL CIRCUIT
// ================================================
	// ifstream genefile;
	// genefile.open("best.gen.dat");
	// TVector<double> genotype(1, VectSize);
	// genefile >> genotype;

	// // //     agent x position, agent y position, Peak x position, Peak y position, steepness
	// BehavioralTraces_Specific(genotype, 10.0, 10.0, 50.0, 50.0, 1.5);
	//  BehavioralTraces_Across_Conditions(genotype, 1.5);
	// // BehavioralTraces(genotype, 455.0, 2.5); 


	return 0;
}
