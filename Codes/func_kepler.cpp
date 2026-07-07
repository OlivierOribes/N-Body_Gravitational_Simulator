//=======================================================================
// Verify Kepler's Third Law
//
// Computes the theoretical and numerical values of the ratio T²/a³ for
// the Earth-Sun system.
//
// The theoretical value is obtained from:
//
//      T² / a³ = 4π² / G(M☉ + M⊕)
//
// while the numerical value is estimated from the simulated orbit using
// the current orbital radius and velocity:
//
//      T ≈ 2πa / v
//
// Both ratios are accumulated over one complete Earth revolution
// (365 simulation steps). Their average values are then written to the
// output file for validation of Kepler's Third Law.
//=======================================================================

void check_kepler_third_law(const vector<Body>& bodies,
                            ofstream& keplerFile,
                            double step){

    const Body& sun = bodies[0]; // The Sun is the central body
    const Body& planet = bodies[3]; // Select Earth
    static double sum_kepler_ratio_theoretical = 0.0;
    static double sum_kepler_ratio_numerical = 0.0;
    static int count = 0;
    static int n = 0;

    // Compute theoretical value
    double kepler_ratio_theoretical = (4 * M_PI * M_PI) / (G * (sun.mass + planet.mass));

    // Compute numerical semi-major axis
    double a = planet.position.norm();

    // Compute numerical orbital period
    double v = planet.velocity.norm();
    double T = 2 * M_PI * a / v; // T ≈ 2πa / v (since T = 2π sqrt(r³ / GM))

    // Compute numerical ratio
    double kepler_ratio_numerical = (T * T) / (a * a * a);

    sum_kepler_ratio_theoretical += kepler_ratio_theoretical;
    sum_kepler_ratio_numerical += kepler_ratio_numerical;
    count++;

    // Average over one complete Earth revolution
    if (count == 365) {
        double avg_kepler_ratio_theoretical = sum_kepler_ratio_theoretical / 365.0;
        double avg_kepler_ratio_numerical = sum_kepler_ratio_numerical / 365.0;

        // Save the averages to the .res file
        keplerFile << n << " " << avg_kepler_ratio_theoretical << " " << avg_kepler_ratio_numerical << endl;

        // Reset accumulators and counter
        sum_kepler_ratio_theoretical = 0.0;
        sum_kepler_ratio_numerical = 0.0;
        count = 0;
        n++;
    }
}

