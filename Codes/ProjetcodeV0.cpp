/**
 * @file Projet_Version_V6.cpp
 * @author Olivier Oribes
 * @brief Newtonian N-body Solar System simulation with energy conservation
 *        analysis and force computation benchmarking.
 *
 * This program simulates the evolution of a simplified Solar System model
 * using Newton's law of universal gravitation.
 *
 * Gravitational accelerations can be computed using two independent methods:
 *   - A classical double-loop implementation
 *   - A matrix-based implementation using the Eigen linear algebra library
 *
 * Time integration is performed using a fourth-order Yoshida symplectic
 * integrator, while a Leapfrog integrator is also available for comparison.
 *
 * During the simulation, the program exports the Cartesian trajectories of
 * every celestial body and continuously evaluates the kinetic, potential,
 * and total mechanical energy of the system in order to verify the numerical
 * conservation properties of the chosen integration scheme.
 *
 * Initial conditions correspond to the NASA JPL Horizons ephemerides at the
 * reference epoch of 13/03/2025.
 *
 * Generated output files:
 *   - data_Yoshida*.txt : individual body trajectories
 *   - positions.txt     : Cartesian coordinates for visualization
 *   - data_Energy.txt   : kinetic, potential and total energy history
 *   - settingdata.txt   : body names and radii for visualization
 *
 * Main features:
 *   - Newtonian N-body gravitational simulation
 *   - Matrix-based force computation using Eigen
 *   - Classical double-loop force computation
 *   - Fourth-order Yoshida symplectic integrator
 *   - Alternative Leapfrog integrator
 *   - Mechanical energy conservation analysis
 *   - Trajectory export for visualization
 *   - Execution time measurement
 *
 * Target system:
 *   - Standard C++17 compiler
 *   - Eigen 3 library
 *
 * @version 6.0
 * @date 2026-07-06
 *
 * @copyright Copyright (c) 2026 Olivier Oribes
 *
 * @license BSD-3-Clause
 */


#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <thread>
#include <chrono>
#include <vector>
#include "Eigen/Dense"

using namespace Eigen;

//-----------------Constants-------------------------------------
//---------------------------------------------------------------
constexpr double DT    = 1;  // Time step (days)
constexpr double SOFT  = 0.01;  // Softening parameter
constexpr int    MAXSTEP = 365*165; //  1 year of steps
constexpr double UA  = 1.495978707e11;      // m
constexpr double DAY = 86400;               // s
constexpr double G_SI = 6.67430e-11;        // m³/kg/s²
constexpr double G = G_SI * DAY * DAY / (UA * UA * UA); // Gravitational constant  UA³/kg/jour²

//-----------------Body structure--------------------------------
//---------------------------------------------------------------
struct BODY
{
    std::string name;
    double mass;         // kg
    double radius;       // km
    Vector3d position;   // (x, y, z) in AU
    Vector3d velocity;   // (vx, vy, vz) in AU/day
    Vector3d half_velocity;
    Vector3d acceleration;

    BODY(const std::string &name_,
         double mass_,
         double radius_,
         const Vector3d &pos_,
         const Vector3d &vel_,
         const Vector3d & /* half_vel_ */,
         const Vector3d & /* acc_ */)
      : name(name_),
        mass(mass_),
        radius(radius_),
        position(pos_),
        velocity(vel_),
        half_velocity(Vector3d::Zero()),
        acceleration(Vector3d::Zero())
    {}
    // second constructor with 5 parameters
    BODY(const std::string &name_,
        double mass_,
        double radius_,
        const Vector3d &position_,
        const Vector3d &velocity_)
    : name(name_),
    mass(mass_),
    radius(radius_),
    position(position_),
    velocity(velocity_),
    half_velocity(Vector3d::Zero()),
    acceleration(Vector3d::Zero())
    {}
};

//-----------------Prototypes------------------------------------
//---------------------------------------------------------------
void leapfrog_method(std::vector<BODY> &bodies, int N);
void yoshida_step(std::vector<BODY>& bodies, double dt, int N);
void compute_forces_matrix(std::vector<BODY> &bodies, int N);
void compute_forces_loop(std::vector<BODY>& bodies);
void display_position(const std::vector<BODY>& bodies, int i);
void energy(const std::vector<BODY> &bodies, double &E_kin, double &E_pot, double &E_tot, int N, int i);
void print(const std::vector<BODY> &bodies, int i);
void save_settings(const std::vector<BODY>& bodies);
std::vector<BODY> initialization();


//-----------------main------------------------------------------
//---------------------------------------------------------------
int main()
{
    auto starttime = std::chrono::high_resolution_clock::now();
    std::vector<BODY> bodies = initialization(); // Initial positions, velocities, radius

    int N = static_cast<int>(bodies.size());    // N bodies

    double E_kin, E_pot, E_tot;                 // Energy

    compute_forces_matrix(bodies, N); // Acceleration initialization Matrix method
    //compute_forces_loop(bodies);    // Acceleration initialization Loop method

    save_settings(bodies);

    // Integration loop
    for(int i = 0; i < MAXSTEP; i++)
    {   

        print(bodies,i);

        //leapfrog_method(bodies, N);
        yoshida_step(bodies, DT, N);

        display_position(bodies,i);

        energy(bodies, E_kin, E_pot, E_tot, N, i);

    }

    auto endtime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endtime - starttime;
    std::cout << "Total time duration : " << elapsed.count() << " s\n";
    
    return 0;
}

//-----------------Function definitions---------------------

//-----------------------------------------------------------------------
// Compute accelerations by building a matrix M and X_x, X_y, X_z
// M(i, j) = G*m_j / (|x_j - x_i|^2 + SOFT^2)^(3/2)  (i != j), else 0
//-----------------------------------------------------------------------
void compute_forces_matrix(std::vector<BODY> &bodies, int N)
{
    //auto starttimematrix = std::chrono::high_resolution_clock::now();

    MatrixXd M(N,N);
    MatrixXd X_x(N,N), X_y(N,N), X_z(N,N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i==j)
            {
                M(i,j) = X_x(i,j) = X_y(i,j) = X_z(i,j) = 0.0;
            }
            else
            {
                Vector3d diff = bodies[j].position - bodies[i].position;
                double dist2 = diff.squaredNorm() + SOFT*SOFT;
                double dist3 = dist2*std::sqrt(dist2);
                M(i,j)   = G*bodies[j].mass/dist3;
                X_x(i,j) = diff.x();
                X_y(i,j) = diff.y();
                X_z(i,j) = diff.z();
            }
        }
    }

    VectorXd acc_x = (M.cwiseProduct(X_x)).rowwise().sum();
    VectorXd acc_y = (M.cwiseProduct(X_y)).rowwise().sum();
    VectorXd acc_z = (M.cwiseProduct(X_z)).rowwise().sum();

    for (int i = 0; i < N; ++i)
    {
        bodies[i].acceleration = Vector3d(acc_x(i), acc_y(i), acc_z(i));


        


    }

    //auto endtimematrix = std::chrono::high_resolution_clock::now();
    /*
    std::chrono::duration<double> elapsedmatrix = endtimematrix - starttimematrix;
    std::cout << "Matrix computation duration : " << elapsedmatrix.count() << " s\n";
    std::this_thread::sleep_for(std::chrono::seconds(1));
    */

}

//-----------------------------------------------------------------------
// Compute accelerations by building a matrix M and X_x, X_y, X_z
// with a double loop
//-----------------------------------------------------------------------
void compute_forces_loop(std::vector<BODY>& bodies) 
{
    //auto starttimeloop = std::chrono::high_resolution_clock::now();

    int N = bodies.size();
    double force;
    double r2, r3;

    // Réinitialisation des accélérations
    for (int i = 0; i < N; i++) 
    {
        for (int d = 0; d < 3; d++) 
        {
            bodies[i].acceleration[d] = 0.0;
        }
    }

    // Calcul des forces gravitationnelles entre chaque paire
    for (int i = 0; i < N; i++) 
    {

        for (int j = i+1; j < N; j++) 
        { // On ne recalcule pas deux fois la même interaction
            std::vector<double> r(3);

            for (int d = 0; d < 3; d++) 
            {
                r[d] = bodies[j].position[d] - bodies[i].position[d];
            }

            // Calcul du dénominateur
            r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]; // distance square
            r3 = pow(r2, 1.5); // we removed softening parameter

            for (int d = 0; d < 3; d++) 
            {
                force = G * bodies[i].mass * bodies[j].mass * r[d] / r3;
                bodies[i].acceleration[d] += force / bodies[i].mass;
                bodies[j].acceleration[d] -= force / bodies[j].mass;
            }
        }
    }

    /*
    auto endtimeloop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedloop = endtimeloop - starttimeloop;
    std::cout << "Loop computation duration : " << elapsedloop.count() << " s\n";
    std::this_thread::sleep_for(std::chrono::seconds(1));
    */
}

//-----------------------------------------------------------------------
// Leapfrog integrator
//-----------------------------------------------------------------------
void leapfrog_method(std::vector<BODY> &bodies,int N)
{   

    //auto starttimeleapfrog = std::chrono::high_resolution_clock::now();

    // half-step velocity
    for (BODY &b : bodies)
    {
        b.half_velocity = b.velocity + 0.5 * b.acceleration * DT;

        // position update
        b.position += b.half_velocity * DT;
    }

    // Compute accelerations
    compute_forces_matrix(bodies,N);

    // final velocity
    for (BODY &b : bodies)
    {
        b.velocity = b.half_velocity + 0.5 * b.acceleration * DT;
    }

    /*
    auto endtimeleapfrog = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedleapfrog = endtimeleapfrog - starttimeleapfrog;
    std::cout << "Leapfrog method duration : " << elapsedleapfrog.count() << " s\n";
    std::this_thread::sleep_for(std::chrono::seconds(1));
    */
}


//-----------------------------------------------------------------------
// Yoshida 4 integrator
//-----------------------------------------------------------------------
void yoshida_step(std::vector<BODY>& bodies, double dt, int N) 
{
    //auto starttimeyoshida = std::chrono::high_resolution_clock::now();

    // Yoshida coeff

    const double w1 = 1.0 / (2.0 * (2.0 - pow(2.0, 1.0 / 3.0)));
    const double w2 = (1.0 - pow(2.0, 1.0 / 3.0)) / (2.0 * (2.0 - pow(2.0, 1.0 / 3.0)));
    const double w3 = w2;
    const double w4 = w1;

    const double c1 = 1.0 / (2.0 - pow(2.0, 1.0 / 3.0));
    const double c2 = -pow(2.0, 1.0 / 3.0) / (2.0 - pow(2.0, 1.0 / 3.0));
    const double c3 = c1;

    // first integration step
    for (auto& b : bodies) {
        b.position += w1 * dt * b.velocity;
    }

    compute_forces_matrix(bodies, N);

    for (auto& b : bodies) {
        b.velocity += c1 * dt * b.acceleration;
    }

    // second integration step
    for (auto& b : bodies) {
        b.position += w2 * dt * b.velocity;
    }

    compute_forces_matrix(bodies, N);

    for (auto& b : bodies) {
        b.velocity += c2 * dt * b.acceleration;
    }

    // Third integration step
    for (auto& b : bodies) {
        b.position += w3 * dt * b.velocity;
    }

    compute_forces_matrix(bodies, N);

    for (auto& b : bodies) {
        b.velocity += c3 * dt * b.acceleration;
    }

    // Update position
    for (auto& b : bodies) {
        b.position += w4 * dt * b.velocity;
    }
    
    /*
    auto endtimeyoshida = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedyoshida = endtimeyoshida - starttimeyoshida;
    std::cout << "Yoshida 4 method duration : " << elapsedyoshida.count() << " s\n";
    std::this_thread::sleep_for(std::chrono::seconds(1));
    */

}
//-----------------------------------------------------------------------
// Print data to files
//-----------------------------------------------------------------------
void print(const std::vector<BODY> &bodies, int i)
{
    // File stay open during simulation
    static std::vector<std::ofstream> files;

    if(i == 0)
    {
        files.resize(bodies.size());

        for(size_t k = 0; k < bodies.size(); ++k)
        {
            std::string filename = "data_Yoshida" + bodies[k].name + ".txt";
            files[k].open(filename);

            if (!files[k].is_open())
            {
                std::cerr << "Error opening file " << filename << " !" << std::endl;
                exit(1);
            }

            // Header
            files[k] << std::left
                     << std::setw(30) << "X[AU]"
                     << std::setw(30) << "Y[AU]"
                     << std::setw(30) << "Z[AU]"<< "\n";
                    /*<< std::setw(30) << "Vx[AU/day]"
                     << std::setw(30) << "Vy[AU/day]"
                     << std::setw(30) << "Vz[AU/day]"
                     << std::setw(30) << "Ax[AU/day]"
                     << std::setw(30) << "Ay[AU/day]"
                     << std::setw(30) << "Az[AU/day]" << "\n";*/
        }
    }

    for (size_t k = 0; k < bodies.size(); ++k)
    {   
        files[k] << std::left << std::fixed << std::setprecision(16)
                 << std::setw(30) << bodies[k].position.x()
                 << std::setw(30) << bodies[k].position.y()
                 << std::setw(30) << bodies[k].position.z()
                 /*<< std::setw(30) << bodies[k].velocity.x()
                 << std::setw(30) << bodies[k].velocity.y()
                 << std::setw(30) << bodies[k].velocity.z()
                 << std::setw(30) << std::scientific << bodies[k].acceleration.x()
                 << std::setw(30) << std::scientific << bodies[k].acceleration.y()
                 << std::setw(30) << std::scientific << bodies[k].acceleration.z()*/
                 << "\n";
    }

    // Close file after last step
    if(i == MAXSTEP - 1)
    {
        for(auto &file : files)
            file.close();
        
        std::cout << "Data printed successfully." << std::endl;
    }
}
void save_settings(const std::vector<BODY>& bodies)
{
    std::ofstream file("settingdata.txt");

    if (!file.is_open())
    {
        std::cerr << "Error opening file settingdata.txt!\n";
        exit(1);
    }

    // Write a header
    file << "# Planet settings\n";
    file << std::left << std::setw(20) << "#Name" 
         << std::setw(20) << "Radius (km)" << "\n";

    // Write each planet's name and radius
    for (const auto& body : bodies)
    {
        file << std::left << std::setw(20) << body.name
             << std::setw(20) << body.radius << "\n";
    }

    file.close();
    std::cout << "Settings saved to settingdata.txt \n" << std::endl;
}
//-----------------------------------------------------------------------
// Print position to create an animation
//-----------------------------------------------------------------------
void display_position(const std::vector<BODY>& bodies, int i)
{
    static std::ofstream file2("positions.txt");

    if (!file2.is_open()) 
    {
        std::cerr << "Error opening file positions.txt !\n";
        exit(1);
    }

    // Write the step number
    file2 << "# Step " << i << "\n";

    for (const auto& b : bodies) 
    {
        // Write name and position (x, y, z) with fixed precision
        file2 << std::left << std::setw(10) << b.name
              << std::fixed << std::setprecision(6)
              << std::setw(15) << b.position.x()
              << std::setw(15) << b.position.y()
              << std::setw(15) << b.position.z()
              << "\n";
    }

    file2 << "\n";

    // Close file after the last step
    if (i == MAXSTEP - 1) 
    {
        file2.close();
    }
}

//-----------------------------------------------------------------------
// Calculate Energy and print them in a file
//-----------------------------------------------------------------------
void energy(const std::vector<BODY>& bodies, double& E_kin, double& E_pot, double& E_tot, int N, int step)
{
    static std::ofstream file1("data_Energy.txt");

    if (!file1.is_open())
    {
        std::cerr << "Error opening file data_Energy.txt!\n";
        exit(1);
    }

    // HEADER
    if (step == 0)
    {
        file1 << std::left << std::setw(10) << "Step"
              << std::setw(20) << "E_kin"
              << std::setw(20) << "E_pot"
              << std::setw(20) << "E_tot" << "\n";
    }

    E_kin = 0.0;
    E_pot = 0.0;

    // Kinetic energy
    for (const auto& b : bodies)
    {
        double V_squared = b.velocity.squaredNorm();
        E_kin += 0.5 * b.mass * V_squared;
    }

    // Potential energy
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            Vector3d diff = bodies[i].position - bodies[j].position;
            double dist = std::sqrt(diff.squaredNorm() + SOFT * SOFT);
            E_pot -= G * bodies[i].mass * bodies[j].mass / dist;
        }
    }

    E_tot = E_kin + E_pot;

    file1 << std::left << std::scientific << std::setprecision(10)
          << std::setw(10) << step
          << std::setw(25) << E_kin
          << std::setw(25) << E_pot
          << std::setw(25) << E_tot << "\n";

    // Close in the end
    if (step == MAXSTEP - 1)
        file1.close();
}


//-----------------------------------------------------------------------
// Initialization
//-----------------------------------------------------------------------
std::vector<BODY> initialization() //Body position to date 13/03/2025
{
    std::vector<BODY> bodies =
    {
        {
            "Sun",
            1988410e24,
            695700.0,
            Vector3d(-5.216296657882820e-03, -5.139436247090928e-03,  1.708593667133784e-04),
            Vector3d( 7.280213226431683e-06, -2.789836694350716e-06, -1.238683299932505e-07)
        },
        {
            "Mercury",
            3.302e23,
            2439.4,
            Vector3d(-2.084520935855677e-01,  2.497359688077958e-01,  3.964071972873255e-02),
            Vector3d(-2.767745812565210e-02, -1.643508208949270e-02,  1.196278856355838e-03)
        },
        {
            "Venus",
            48.685e23,
            6051.84,
            Vector3d(-7.019861769318310e-01,  1.683091426476831e-01,  4.275646977214004e-02),
            Vector3d(-4.979898899535033e-03, -1.972132190568706e-02,  1.679721332602159e-05)
        },
        {
            "Earth",
            5.97219e24,
            6378.137,
            Vector3d(-9.901283494896050e-01,  1.278355489635163e-01,  1.679920267418484e-04),
            Vector3d(-2.571576631925831e-03, -1.711046177466572e-02,  1.472251630093178e-06)
        },
        {
            "Mars",
            6.4171e23,
            3389.92,
            Vector3d(-1.292422638475013e+00,  1.041698122892262e+00,  5.367433101452403e-02),
            Vector3d(-8.295770650586537e-03, -9.665330256906203e-03,  9.994132590083996e-07)
        },
        {
            "Jupiter",
            1.8982e27,
            71492.0,
            Vector3d( 5.150522486733893e-01,  5.074583148712015e+00, -3.257009678167901e-02),
            Vector3d(-7.595280620836376e-03,  1.121229025888626e-03,  1.652490291241097e-04)
        },
        {
            "Saturn",
            5.6834e26,
            60268.0,
            Vector3d( 9.498260015567574e+00, -1.379633365740638e+00, -3.541863996348656e-01),
            Vector3d( 4.922723455925417e-04,  5.509936081894204e-03, -1.157511325862617e-04)
        },
        {
            "Uranus",
            86.813e24,
            3389.92,
            Vector3d( 1.086486894342725e+01,  1.623349864770252e+01, -8.046535777081534e-02),
            Vector3d(-3.297459488770539e-03,  2.004283097715195e-03,  5.019408215753177e-05)
        },
        {
            "Neptune",
            102.409e24,
            24624,
            Vector3d( 2.987665737327115e+01, -4.149336646400792e-01, -6.799934441949956e-01),
            Vector3d( 2.280989371154266e-05,  3.157029843809297e-03, -6.586457086506603e-05)
        },
    };

    return bodies;
}
