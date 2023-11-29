// *******************************************************************************************
//                             Gabriel Juliano Severino 
//                                      2023
//
//            Chemoreceptive Navigation in Static and Fluctuating Environments
//
// TO DO:
// 
// 
// *******************************************************************************************

#include "Fluid.h"
#include "TSearch.h"
#include "Sniffer.h"
#include "CTRNN.h"
#include "random.h"
#include <random>

#define PRINTOFILE

// Task params
const double StepSize = 0.01;
const double RunDuration = 100; // 6000 and 5000 for transient 
const double TransDuration = 80.0; // Transient duration 
const double EvalDuration = RunDuration - TransDuration; // Evaluation duration

const double SpaceHeight = 100.0; // Size of the space
const double SpaceWidth = 100.0; 

const double peakPositionX = 50.0; // Position of the chemical source
const double peakPositionY = 50.0; 

// Plume params

const double nu = 0.0000001;             // Kinematic viscosity
const double diffusionRate = 0.0;
const double dt_fluid = 0.1; // Time step for NS simulation


const int nsStepsPerAgentStep = static_cast<int>(1); // 10 steps per agent step StepSize / deltaTimeNS
const int totalAgentSteps = static_cast<int>(RunDuration / StepSize); 





// EA params
const int POPSIZE = 1; //500 
const int GENS = 1; //100
const double MUTVAR = 0.2;
const double CROSSPROB = 0.05;
const double EXPECTED = 1.1;
const double ELITISM = 0.02;

// Nervous system params
const int N = 5;
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


double FitnessFunctionChemoIndexFLUX(TVector<double> &genotype, RandomState &rs)
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

///////////////////////


    Fluid fluid(dt_fluid, diffusionRate, nu);

    int sourceX = 50;      // Example x-coordinate for odor source
    int sourceY = 50;      // Example y-coordinate for odor source
    float odorAmount = 10; // Amount of odor to add
    float velocityX = 2;   // Velocity in the x-direction
    float velocityY = 2;   // Velocity in the y-direction

    double totalFit = 0.0;
    int trials = 0;

        for (double theta = 0.0; theta < 2*M_PI; theta += M_PI/2) {  

            double x = rs.UniformRandom(0.0, SpaceWidth);
            double y = rs.UniformRandom(0.0, SpaceHeight);

            // Calculate initial distance
            double initialDist = sqrt(pow(x - sourceX, 2) + pow(y - sourceY, 2));
            if (initialDist < 1.0) initialDist = 1.0; // Avoid division by zero
			// printf("Initial distance: %f\n", initialDist);

            // Set agent's position
            Agent.Reset(x, y, theta);

            double dist = 0.0;
            int agentStepCounter = 0;
        for (double time = 0; time < RunDuration; time += StepSize) {

				// Calculate chemical gradient at Sniffer position
				double concentration = fluid.getOdorConcentration(Agent.posX, Agent.posY);

                        // Check if it's time to update the fluid simulation (Step every 10 agent steps)
                if (agentStepCounter % 10 == 0) {
                    fluid.addOdor(sourceX, sourceY, odorAmount);
                    fluid.addVelocity(sourceX, sourceY, velocityX, velocityY);
                    fluid.step(); 
                }

				// Sense the gradient
				Agent.Sense(concentration, time);
					
				// Move based on sensed gradient
				Agent.Step(StepSize);
                                                                
            // Evaluate fitness after transient duration
            if (time > TransDuration) {
                double dx = std::abs(Agent.posX - sourceX);
                double dy = std::abs(Agent.posY - sourceY);
                dist += sqrt(dx * dx + dy * dy);
                agentStepCounter++;
            }
        }

        double totaldist = (dist / (EvalDuration / StepSize));
        double fitnessForThisTrial = (initialDist - totaldist) / initialDist;
        fitnessForThisTrial = fitnessForThisTrial < 0.0 ? 0.0 : fitnessForThisTrial; // Ensure non-negative fitness

        totalFit += fitnessForThisTrial;
        trials++;
    }

    return totalFit / trials;
}


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
        for (double theta = 0.0; theta < 2*M_PI; theta += M_PI/2) {  

            double x = rs.UniformRandom(0.0, SpaceWidth);
            double y = rs.UniformRandom(0.0, SpaceHeight);

            // Peak position of chemical gradient
            const double peakPositionX = rs.UniformRandom(0.0, SpaceWidth);
            const double peakPositionY = rs.UniformRandom(0.0, SpaceHeight);

            // Calculate initial distance
            double initialDist = sqrt(pow(x - peakPositionX, 2) + pow(y - peakPositionY, 2));
            if (initialDist < 1.0) initialDist = 1.0; // Avoid division by zero
			// printf("Initial distance: %f\n", initialDist);

            // Set agent's position
            Agent.Reset(x, y, theta);

            double dist = 0.0;
            double totalDist = 0.0;

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
                    dist += sqrt(dx * dx + dy * dy);
                }
            }

            double totaldist = (dist / (EvalDuration / StepSize));
            double fitnessForThisTrial = (initialDist - totaldist)/initialDist;
            fitnessForThisTrial = fitnessForThisTrial < 0.0 ? 0.0 : fitnessForThisTrial; // Ensure non-negative fitness

            totalFit += fitnessForThisTrial;
            trials++;
        }
    }
    return totalFit / trials;
}

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
            Agent.Reset(x1, y1, 0.0);

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
        Agent.Reset(initPosX, initPosY, 0.0);

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
	// s.SetEvaluationFunction(FitnessFunctionChemoIndex); 
    s.SetEvaluationFunction(FitnessFunctionChemoIndexFLUX);

	s.ExecuteSearch();


// ================================================
// B. MAIN FOR ANALYZING A SUCCESFUL CIRCUIT
// ================================================
	// ifstream genefile;
	// genefile.open("best.gen.dat");
	// TVector<double> genotype(1, VectSize);
	// genefile >> genotype;

	 // Get the index from the argument
    // int index = std::stoi(argv[1]);
	
	    // Construct the filename
    // std::string filename = "Evolved_Agents/StaticSensorStaticEnv/5N/Genotypes/best.gen_" + std::to_string(index) + ".dat";
    
    
    
    
    // std::string filename = "best.gen.dat";


    // ifstream genefile;
	// genefile.open(filename);
	// TVector<double> genotype(1, VectSize);
	// genefile >> genotype;

	// // //     agent x position, agent y position, Peak x position, Peak y position, steepness
	// BehavioralTraces_Specific(genotype, 10.0, 10.0, 50.0, 50.0, 1.5);
	// BehavioralTraces_Across_Conditions(genotype, 1.5);
	// // // // BehavioralTraces(genotype, 455.0, 2.5); 








	return 0;
}
