//=======================================================================
// Compute the mechanical energy of the N-body system
//
// Computes the kinetic, potential and total mechanical energy of the
// simulated system.
//
// The kinetic energy is evaluated as:
//
//      E_kin = Σ (1/2 m v²)
//
// while the gravitational potential energy is obtained by summing the
// pairwise interactions between all bodies:
//
//      E_pot = - Σ G m_i m_j / r_ij
//
// The total mechanical energy is finally computed as:
//
//      E_tot = E_kin + E_pot
//
// These quantities are used to assess the numerical conservation of
// energy throughout the simulation.
//=======================================================================
void compute_energies(const vector<Body>& bodies,
                      double& kinetic,
                      double& potential,
                      double& total){
    int N = bodies.size();
    double r, r2, v2, diff;
    kinetic = 0.0;
    potential = 0.0;

    // Kinetic energy
    for (const auto& b : bodies) {
        v2 = b.velocity[0] * b.velocity[0] + b.velocity[1] * b.velocity[1] + b.velocity[2] * b.velocity[2];
        kinetic += 0.5 * b.mass * v2;
    }

    // Potential energy
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            r2 = 0.0;
            for (int d = 0; d < 3; d++) {
                diff = bodies[j].position[d] - bodies[i].position[d];
                r2 += diff * diff;
            }
            r = sqrt(r2); // + theta * theta [add softening parameter]
            potential -= G * bodies[i].mass * bodies[j].mass / r;
        }
    }

    total = kinetic + potential;
}

//=======================================================================
// Save normalized mechanical energies
//
// Writes the normalized kinetic, potential and total mechanical
// energies to the output file.
//
// Each energy is divided by its corresponding initial value in order
// to monitor the relative conservation of energy during the numerical
// integration.
//
// Output format:
//
//      Step  Ekin/Ekin₀  Epot/Epot₀  Etot/Etot₀
//=======================================================================
void save_energies(ofstream& energyFile,
                   double step,
                   double kinetic,
                   double potential,
                   double total,
                   double kinetic0,
                   double potential0,
                   double total0){
    energyFile << step << " " << kinetic / kinetic0 << " " << potential / potential0 << " " << total / total0 << endl;
}
