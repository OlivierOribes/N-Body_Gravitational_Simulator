/**
 * @file Projet_Version_V8.cpp
 * @author Olivier Oribes
 * @brief N-body Solar System simulation using Newtonian gravity.
 *
 * This program simulates the evolution of a 23-body Solar System model using
 * Newton's law of universal gravitation.
 *
 * The equations of motion are integrated using a fourth-order symplectic
 * Yoshida integrator (with an alternative Leapfrog implementation available
 * for comparison). Gravitational accelerations are computed using matrix
 * operations provided by the Eigen linear algebra library in order to
 * efficiently evaluate all pairwise interactions.
 *
 * Initial positions and velocities are extracted from the NASA JPL Horizons
 * ephemerides (reference date: 13/03/2025). During the simulation, the program
 * periodically exports the Cartesian coordinates of every body for
 * post-processing and orbital visualization.
 *
 * The simulation also generates a configuration file containing the physical
 * characteristics (name and radius) of each celestial body, allowing external
 * visualization tools to correctly scale rendered objects.
 *
 * Main features:
 *   - Newtonian N-body gravitational simulation
 *   - Matrix-based force computation using Eigen
 *   - Fourth-order Yoshida symplectic integrator
 *   - Alternative Leapfrog integrator
 *   - 3D Cartesian trajectory export
 *   - NASA Horizons initial conditions
 *   - Execution time measurement
 *
 * Target system:
 *   - Standard C++17 compiler
 *   - Eigen 3 library
 *
 * @version 8.0
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

//=========================Constants=====================================
//=======================================================================
constexpr double DT    = 1.0/24;  // Time step (Days)
constexpr int    MAXSTEP = (365*24)*165; //  number of steps
constexpr double UA  = 1.495978707e11;      // m
constexpr double DAY = 86400;               // s
constexpr double G_SI = 6.67430e-11;        // m³/kg/s²
constexpr double G = G_SI * DAY * DAY / (UA * UA * UA); // Gravitational constant  UA³/kg/jour²

//==========================Body structure===============================
//=======================================================================
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

//=========================Prototypes====================================
//=======================================================================
void leapfrog_method(std::vector<BODY> &bodies, int N);
void yoshida_step(std::vector<BODY>& bodies, double dt, int N);
void compute_forces_matrix(std::vector<BODY> &bodies, int N);
void save_settings(const std::vector<BODY>& bodies);
std::vector<BODY> initialization();
void display_position(const std::vector<BODY>& bodies, int i);

//============================  MAIN  ===================================
//=======================================================================
int main()
{
    auto starttime = std::chrono::high_resolution_clock::now();

    std::vector<BODY> bodies = initialization(); // Initial positions, velocities, radius (From Nasa program)
    int N = static_cast<int>(bodies.size());    // N bodies

    compute_forces_matrix(bodies, N); // Acceleration initialization Matrix method
    //compute_forces_loop(bodies);    // Acceleration initialization Loop method

    save_settings(bodies);

    // Integration loop
    for(int i = 0; i < MAXSTEP; i++)
    {   
        
        display_position(bodies,i);
        //leapfrog_method(bodies, N);
        yoshida_step(bodies, DT, N);
        
    }

    auto endtime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endtime - starttime;
    std::cout << "Total time duration : " << elapsed.count() << " s\n";
    
    return 0;
}

//=======================Function definitions============================

//=======================================================================
// Compute accelerations by building a matrix M and X_x, X_y, X_z
// M(i, j) = G*m_j / (|x_j - x_i|^2 + SOFT^2)^(3/2)
//=======================================================================
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
                double dist2 = diff.squaredNorm();
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
//=======================================================================
// Leapfrog integrator
//=======================================================================
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


//=======================================================================
// Yoshida 4 integrator
//=======================================================================
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
    for (BODY& b : bodies) 
    {              
        b.position += w1 * dt * b.velocity;
    }

    compute_forces_matrix(bodies, N);

    for (BODY& b : bodies) 
    {
        b.velocity += c1 * dt * b.acceleration;
    }

    // second integration step
    for (BODY& b : bodies) 
    {
        b.position += w2 * dt * b.velocity;
    }

    compute_forces_matrix(bodies, N);

    for (BODY& b : bodies) 
    {
        b.velocity += c2 * dt * b.acceleration;
    }

    // Third integration step
    for (BODY& b : bodies) 
    {
        b.position += w3 * dt * b.velocity;
    }

    compute_forces_matrix(bodies, N);

    for (BODY& b : bodies) 
    {
        b.velocity += c3 * dt * b.acceleration;
    }

    // Update position
    for (BODY& b : bodies) 
    {
        b.position += w4 * dt * b.velocity;
    }
    
    /*
    auto endtimeyoshida = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedyoshida = endtimeyoshida - starttimeyoshida;
    std::cout << "Yoshida 4 method duration : " << elapsedyoshida.count() << " s\n";
    std::this_thread::sleep_for(std::chrono::seconds(1));
    */

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
//=======================================================================
// Print position to create a gif
//=======================================================================
// Function to write the positions of simulation bodies and their corresponding trajectory data
void display_position(const std::vector<BODY>& bodies,int i)
{
    // File stay open during simulation
    static std::ofstream file2("positions_lowstep.txt");

    if (!file2.is_open()) 
    {
        std::cerr << "Error opening file positions_lowstep.txt!\n";
        exit(1);
    }

    if((i%24==0))
        {
        // Write the current simulation step number
        file2 << "# Step " << i << "\n";


        // Loop over each body and write its current cartesian position
        for (const auto& b : bodies) 
        {
            file2 << std::left << std::setw(10) << b.name
            << std::fixed << std::setprecision(6)
            << std::setw(15) << b.position.x()
            << std::setw(15) << b.position.y()
            << std::setw(15) << b.position.z()
            << "\n";
        }
    

        file2 << "\n";
    }

    if (i == MAXSTEP - 1) 
    {
        file2.close();
        std::cout << "File saved in positions_lowstep.txt!\n";
    }
}


//=======================================================================
//=======================================================================
// Initialization
//=======================================================================
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
            "Moon",
            7.342e22,
            1737.4,
            Vector3d(-9.925938711509329e-01, 1.288367006296058e-01,  2.486981337651200e-04),
            Vector3d( -2.804579923753571e-03, -1.763131248304115e-02, -4.775310362211178e-05)
        },
        {
            "Mars",
            6.4171e23,
            3389.92,
            Vector3d(-1.292422638475013e+00,  1.041698122892262e+00,  5.367433101452403e-02),
            Vector3d(-8.295770650586537e-03, -9.665330256906203e-03,  9.994132590083996e-07)
        },
        {
            "Cérès",
            9.46e20,
            469.7,
            Vector3d(2.565019648633334e+00,  -1.432703630657016e+00,  -5.185122518609194e-01),
            Vector3d(4.577003423133978e-03, 8.366352227268607e-03,  -5.771505205964025e-04)
        },
        {
            "Vesta",
            2.7e20,
            262.7,
            Vector3d(-1.980477452873451e+00,  -9.651812081728403e-01,  2.691986187706480e-01),
            Vector3d(5.859289990811343e-03, -1.038296780705886e-02,  -4.045169231304599e-04)
        },
        {
            "Pallas",
            2.11e20,
            256.5,
            Vector3d(1.174270033605817e+00,  -2.652286271056669e+00,  1.732438835228203e+00),
            Vector3d(7.932086012594429e-03, 1.578119089756420e-03,  -1.779623439180384e-03)
        },
        {
            "Hygiea",
            1.0e20,
            203.56,
            Vector3d(2.261784015369848e+00,  2.635659253724826e+00,  1.882420562169424e-01),
            Vector3d(-6.494801882476291e-03, 5.763872507002603e-03,  -3.364925699815282e-04)
        },
        {
            "Jupiter",
            1.8982e27,
            71492.0,
            Vector3d( 5.150522486733893e-01,  5.074583148712015e+00, -3.257009678167901e-02),
            Vector3d(-7.595280620836376e-03,  1.121229025888626e-03,  1.652490291241097e-04)
        },
        {
            "Ganymède",
            1.481e23,
            2631.0,
            Vector3d( 5.194941870098311e-01, 5.068782456316153e+00, -3.273463365578770e-02),
            Vector3d(-2.522526215489928e-03,  4.836863610819146e-03,  3.796650417889915e-04)
        },
        {
            "Europa",
            4.5e22,
            1560.8,
            Vector3d( 5.197119708520060e-01, 5.073893900927814e+00, -3.249969893705672e-02),
            Vector3d(-6.364878610546718e-03,  8.958980865491096e-03,  4.445929613139840e-04)
        },
        {
            "Io",
            8.93e22,
            1821.49,
            Vector3d( 5.147533234396415e-01, 5.071769835587237e+00, -3.268096251354626e-02),
            Vector3d(2.203784785972939e-03,  -7.243787499494322E-04,  2.389488981407589E-04)
        },
        {
            "Callisto",
            1.075938e23,
            2410.3,
            Vector3d( 5.064601717059863e-01, 5.065474525741683e+00, -3.297738591680897e-02),
            Vector3d(-4.200113375600487e-03,  -2.144164169362075e-03,  1.087603124696974e-04)
        },
        {
            "Chiron",
            2e18,
            83.0,
            Vector3d( 1.715639111460142e+01, 7.049099432341470e+00, 2.730008991046606e-01),
            Vector3d(-1.603643521422809e-03,  2.744877268171411e-03,  -3.866914287914001e-04)
        },
        {
            "Chariklo",
            7e18,
            151.0,
            Vector3d( 1.282122374241851e+01, -1.166506550153236e+01, 2.228145193931309e+00),
            Vector3d(2.635684058937510e-03,  2.421134872526948e-03,  1.513813962920798e-03)
        },
        {
            "Saturn",
            5.6834e26,
            60268.0,
            Vector3d( 9.498260015567574e+00, -1.379633365740638e+00, -3.541863996348656e-01),
            Vector3d( 4.922723455925417e-04,  5.509936081894204e-03, -1.157511325862617e-04)
        },
        {
            "Titan",
            1.345e23,
            2575.5,
            Vector3d( 9.489870412550033e+00, -1.379719163129494e+00, -3.534606054832699e-01),
            Vector3d(5.633531534934078e-04, 2.726225455999667e-03,  1.312817922747998e-03)
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
        {
            "Pluto",
            1.307e22,
            1188.3,
            Vector3d( 1.841869962011130e+01, -2.994616254175570e+01, -2.123356584602299e+00),
            Vector3d(2.762654795734760e-03,  9.601212984816840e-04,  -8.906515241853556e-04)
        },
        {
            "Eris",
            1.646e22,
            1163.0,
            Vector3d( 8.534906490011645e+01, 3.920501585419132e+01, -1.777547216715094e+01),
            Vector3d(-4.626371396454441e-04,   8.615744578508323e-04,  9.355372265221904e-04)
        },

    };


    return bodies;
}
