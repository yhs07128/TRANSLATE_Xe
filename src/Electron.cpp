#include <array>
#include <chrono>
#include <fstream>
#include <string>

#include "Electron.h"

/*
 * Converts electronvolts to joules
 * 
 * @param eV Energy (in eV)
 * @return Energy (in J)
 */
inline double eV_to_J(double eV)
{
    return eV * 1.60218e-19;
}

/*
 * Converts joules to electronvolts
 * 
 * @param J Energy (in J)
 * @return Energy (in eV)
 */
inline double J_to_eV(double J)
{
    return J * 6.242e18;
}

/*
 * Given an electron's kinetic energy, returns its velocity (classically)
 * 
 * @param eV The electron's kinetic energy (in eV)
 * @return The velocity of the electron (in m/s)
 */
inline double eV_to_v(double eV)
{
    double J = eV_to_J(eV);
    return sqrt(2 * J / m_e);
}

/*
 * Produces a random unit vector
 * 
 * @param gen The random number generator to be used
 * @return Random unit vector
 */
Vec random_unit_vector(std::mt19937& gen)
{
    std::normal_distribution<double> dist(0.0, 1.0);

    Vec v(dist(gen), dist(gen), dist(gen));

    if (norm(v) != 0)
        v /= norm(v);

    return v;
}

/*
 * Return a velocity whose magnitude is chosen from the Maxwell-Boltzmann distribution for
 * argon / the electron. The velocity's direction is randomly chosen
 * 
 * @param gen The random number generator to be used
 * @param electron Whether the velocity should be chosen from the electron distribution. If false, it draws from the argon distribution
 * @return The velocity vector (in m/s)
 */
Vec random_velocity(std::mt19937& gen, bool electron=false)
{
    if (!interactions)
        return Vec(0, 0, 0);

    double a;

    if (!electron)
        a = sqrt(k_b * T / M_a);
    else
        a = sqrt(k_b * T / m_e);

    std::gamma_distribution<double> maxwell(1.5, 1 / (2 * a * a));
    Vec v = random_unit_vector(gen) * sqrt(maxwell(gen));

    return v;
}

/*
 * Chooses the starting position of simulated electron randomly from a disc at the
 * top of the simulated volume. If running the simulation in a uniform field, this returns
 * the origin.
 * 
 * @param gen The random number generator to be used
 * @return The position vector (in m)
 */
Vec starting_pos(std::mt19937& gen)
{
    if (uniform_field)
        return Vec(0, 0, 0);

    std::normal_distribution<double> dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    Vec random_disc(0, dist(gen), dist(gen));

    if (norm(random_disc) != 0)
        random_disc /= norm(random_disc);

    random_disc *= uniform(gen);
    random_disc *= (r_max - r_min) * spawn_width_scale;

    return Vec((z_max - z_min) * spawn_height_scale + z_min, random_disc.y, random_disc.z);
}

/*
 * Given a parallel coordinate to the tip array plane (y or z), this shifts it to 
 * within the bounds of a single tip volume. Effectively, this tiles the simulated volume
 * like a chessboard.
 * 
 * @param x The coordinate to be shifted (in m)
 * @return The shifted coordinate (in m)
 */
double shift_coord(double x)
{
    bool negative = signbit(x);
    double width = r_max - r_min;
    
    if (negative) {
        while ((x + 2 * width) < r_max)
            x += 2 * width;
    } else {
        while ((x - 2 * width) > -r_max)
            x -= 2 * width;
    }

    return x;
}

/*
 * Shifts the y and z coordinates of a position vector to lie within the bounds of
 * a single tip, for the purposes of calculating the force. This function is used
 * only when simulating a tip array.
 * 
 * @param pos The position vector to be altered (in m)
 * @return The shifted position vector
 */
Vec local_coords(Vec pos)
{
    return Vec(pos.x, shift_coord(pos.y), shift_coord(pos.z));
}

/*
 * Gives the radial coordinate of the simulated electron. If the electron lies outside
 * the simulated volume, its radial coordinate is chosen to be the maximum of the volume.
 * 
 * @param local_pos The position vector of the electron relative to the nearest tip (in the single tip case, this is the origin) (in m)
 * @return The electron's radial coordinate
 */
double r(Vec local_pos)
{
    double r_val = sqrt(local_pos.y * local_pos.y + local_pos.z * local_pos.z);

    if (r_val > r_max)
        r_val = r_max;

    return r_val;
}

/*
 * Gives the height cordinate of the simulated electron. This is just x, so this function
 * effectively acts like a clamp.
 * 
 * @param local_pos The position vector of the electron relative to the nearest tip (in m) (in the single tip case, this is the origin)
 * @return The electron's height coordinate
 */
double z(Vec local_pos)
{
    double z_val = local_pos.x;

    if (z_val > z_max)
        z_val = z_max;

    if (z_val < 0)
        z_val = 0;

    return z_val;
}

/*
 * Given a position vector, returns whether it lies within or on the surface of the
 * simulated tip (or tip array).
 * 
 * @param pos The position vector the check (in m)
 * @return Whether the vector touches the tip
 */
bool hit_check(Vec pos)
{
    Vec local = local_coords(pos);

    double r_val;
    double z_val;
    
    if (!single_tip) {
        r_val = r(local);
        z_val = z(local);
    } else {
        r_val = r(pos);
        z_val = z(pos);
    }

    if (z_val <= z_min) 
        return true;

    // For rounded tips
    if (rounded_tip) {
        // For cylindrical sides
        if (top_radius == bot_radius) {
            if ( ((z_val - top_z) * (z_val - top_z) + r_val * r_val <= top_radius * top_radius) || (z_val <= top_z && r_val <= top_radius) )
                return true;
        // For tapered sides
        } else {
            // For cylindrical sides after a tapered section
            if (z_val < bot_z) {
                if (r_val < bot_radius)
                    return true;
                else
                    return false;
            }
            if ( ((z_val - top_z) * (z_val - top_z) + r_val * r_val <= top_radius * top_radius) || (z_val <= top_z && (z_val <= (bot_z - top_z) / (bot_radius - top_radius) * (r_val - bot_radius) + bot_z)) )
                return true;
        }
    // For flat tips
    } else {
        // For cylindrical sides
        if (top_radius == bot_radius) {
            if (z_val <= top_z && r_val <= top_radius)
                return true;
        // For tapered sides
        } else {
            if (z_val <= top_z && (z_val <= (bot_z - top_z) / (bot_radius - top_radius) * (r_val - bot_radius) + bot_z))
                return true;
        }
    }

    return false;
}

/*
 * Gives the radial index to access the electric field lookup table.
 * 
 * @param r_val The radial coordinate (in m)
 * @return The associated index
 */
int r_index(double r_val)
{
    return floor(r_conversion_factor * (r_val - r_min));
}

/*
 * Gives the height index to access the electric field lookup table.
 * 
 * @param z_val The height coordinate (in m)
 * @return The associated index
 */
int z_index(double z_val)
{
    return floor(z_conversion_factor * (z_val - z_min));
}

/*
 * Gives the acceleration felt by the electron in a non-uniform field.
 * 
 * @param pos The position vector of the electron (in m)
 * @param volts The bulk field (in V/cm)
 * @return The acceleration of the electron (in m/s^2)
 */
Vec accel_from_E(Vec pos, double volts)
{
    Vec local = local_coords(pos);

    double r_val;
    double z_val;
    
    if (!single_tip) {
        r_val = r(local);
        z_val = z(local);
    } else {
        r_val = r(pos);
        z_val = z(pos);
    }

    int i = r_index(r_val);
    int j = z_index(z_val);
    
    double E_r_val = E_r[num_z_steps * i + j];
    double E_z_val = E_z[num_z_steps * i + j];

    Vec inward_vec;
    
    if (!single_tip)
        inward_vec = Vec(0, -local.y, -local.z);
    else
        inward_vec = Vec(0, -pos.y, -pos.z);

    if (norm(inward_vec) != 0)
        inward_vec /= norm(inward_vec);

    inward_vec *= E_r_val;

    Vec e_field(-E_z_val, inward_vec.y, inward_vec.z);

    return (volts / bulk_field) * (e / m_e) * e_field * 1e2;
}

/*
 * The constructor for the electron class. Sets up the basic initial conditions.
 * 
 * @param ions A parameter passed to child electrons to keep track of the total number of ionizations
 * @param volts The bulk field (in V/cm)
 * @param position The initial position vector (in m)
 * @param velocity The initial velocity vector (in m)
 * @param gen The random number generator to be used
 */
Electron::Electron(int& ions, double volts, Vec position, Vec velocity, std::mt19937& gen):
    ionized(ions), child_ions(0), x(position), v(velocity), accel((e / m_e) * volts * 1e2, 0, 0), volts_per_cm(volts), total_time(0), generator(gen)
{
    if (!uniform_field)
        accel = accel_from_E(x, volts_per_cm);

    energy = J_to_eV(0.5 * m_e * dot(v, v));

}

/*
 * Gives the index to be used when accessing cross-section lookup tables.
 * 
 * @param v The velocity of the electron (in m/s) (this is sometimes the velocity of the electron relative to a chosen argon atom)
 * @return The associated index
 */
int Electron::index(double v) const
{
    double E = J_to_eV(0.5 * m_e * v * v);

    if (E < eV_min)
        E = eV_min;

    int ind = round(((eV_steps - 1) / log10(eV_max / eV_min)) * (log10(E) - log10(eV_min)));
    assert(ind < 1000);

    return ind;
}

/*
 * Removes an amount of energy from the electron, then updates its velocity to match.
 * 
 * @param eV The energy to be removed (in eV)
 */
void Electron::remove_energy(double eV)
{
    energy -= eV;

    if (energy < 0)
        energy = 0;
    
    v = sqrt((2 * eV_to_J(energy)) / m_e) * (v / norm(v));

    return;
}

/*
 * Removes the ionization energy of argon from the electron, then, if using the full
 * algorithm, spawns a new electron at the current one's position. The ionization counter is
 * increased by one.
 */
void Electron::ionization()
{
    if (use_recursive_ionization) {
        std::uniform_real_distribution<double> rand_eV(1.0, 5.0);
        Vec near_therm_vel = random_unit_vector(generator) * eV_to_v(rand_eV(generator));
        ionization_electrons.push_back(new Electron(ionized, volts_per_cm, x, near_therm_vel, generator));
    }

    remove_energy(15.76);

    ionized++;
    child_ions++;

    return;
}

/*
 * Updates the electron's velocity, assuming an elastic collision with an argon
 * atom.
 * 
 * @param u The magnitude of the electron's velocity relative to a chosen argon atom (in m/s)
 * @param vm The chosen argon atom's velocity vector (in m/s)
 */
void Electron::elastic_collision(double u, Vec vm)
{
    Vec unit_vec = random_unit_vector(generator);

    Vec a = (M_a * u) * unit_vec;
    Vec b = (m_e) * v;
    Vec c = (M_a) * vm;

    Vec v1 = (a + b + c) / (m_e + M_a);

    double r = x_sec_p(u) / std::max(x_sec_e(u), x_sec_p(u));

    std::uniform_real_distribution<double> rand_momentum_change(0.0, 1.0);

    if (rand_momentum_change(generator) < r)
        v = v1;
    else
        v = norm(v1) / norm(v) * v;

    return;
}

/*
 * Classically updates the electron's position and velocity.
 */
void Electron::update_pos_vel()
{
    Vec x_old = x;
    
    if (!uniform_field) {
        x += v * time_to_collision + 0.5 * accel_from_E(x, volts_per_cm) * time_to_collision * time_to_collision;
        v += 0.5 * (accel_from_E(x_old, volts_per_cm) + accel_from_E(x, volts_per_cm)) * time_to_collision;
    } else {
        x += v * time_to_collision + 0.5 * accel * time_to_collision * time_to_collision;
        v += accel * time_to_collision;
    }

    return;
}

/*
 * Draws the next timestep to use.
 * 
 * @param u The magnitude of the electron's velocity relative to a chosen argon atom (in m/s)
 * @return The drawn timestep (in s)
 */
double Electron::next_collision(double u)
{
    std::exponential_distribution<double> exponential(lambda);
    double time_step = exponential(generator);

    if (time_step > 3 * beta)
        time_step = 3 * beta;

    return time_step;
}

/*
 * Gives the probability for an interaction.
 * 
 * @param u The magnitude of the electron's velocity relative to a chosen argon atom (in m/s)
 * @param x_sec The value of the interaction's cross-section at the kinetic energy associated with the above relative velocity (in cm^2)
 * @return The associated interaction probability
 */
double Electron::probability(double u, double x_sec)
{
    double prob = (u * 1e2 * x_sec) / K_max;
    assert (prob < 1);
    return prob;
}

/*
 * Advances the simulation by one step. In this order, this function...
 *      1. Updates all child electrons
 *      2. Draws a timestep
 *      3. Updates the parent electron's position and velocity
 *      4. Determines the probability of all interactions
 *      5. Simulates the interaction that passes the probability check
 *      6. Updates the parent electron's energy and total simulated time 
 */
void Electron::update()
{
    // Updates all child electrons
    for (auto elec : ionization_electrons) {
        elec->update();
        // If a child electron hits a tip, end its simulation
        if (!uniform_field)
            ionization_electrons.erase(std::remove_if(ionization_electrons.begin(), ionization_electrons.end(), [](auto const& i){ return hit_check(i->x); }), ionization_electrons.end());
    }
    
    // Draw a timestep and update the electron's position and velocity
    time_to_collision = next_collision(norm(v));
    update_pos_vel();

    // If interactions are turned off, implement a terminal velocity
    if (!interactions) {
        if (norm(v) > 1000) {
            v /= norm(v);
            v *= 1000;
        }
    }
    
    // Draw an argon vector and determine interaction probabilities
    energy = J_to_eV(0.5 * m_e * dot(v, v));

    Vec vm = random_velocity(generator);
    double u = norm(v - vm);
    
    double col_prob = probability(u, x_sec_e(u));
    double ion_prob = probability(u, x_sec_i(u));

    double ex_11_prob = probability(u, x_sec_ex(u, 11));
    double ex_13_prob = probability(u, x_sec_ex(u, 13));
    double ex_14_prob = probability(u, x_sec_ex(u, 14));
    double ex_15_prob = probability(u, x_sec_ex(u, 15));

    assert(col_prob + ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob <= 1);

    // If interactions are turned off, set all probabilities to 0
    if (!interactions) {
        col_prob = 0;
        ion_prob = 0;
        ex_11_prob = 0;
        ex_13_prob = 0;
        ex_14_prob = 0;
        ex_15_prob = 0;
    }
    
    std::uniform_real_distribution<double> rand_event(0.0, 1.0);
    double prob = rand_event(generator);
    
    // Simulate the interaction that passes the check
    if (prob < ion_prob) {
        ionization();
    } else if (prob < (ion_prob + ex_11_prob)) {
        remove_energy(11.68);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob)) {
        remove_energy(13.21);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob)) {
        remove_energy(14.10);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob)) {
        remove_energy(15.23);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob + col_prob)) {
        elastic_collision(u, vm);
    }
    
    // Update the electron's energy and total time simulated
    energy = J_to_eV(0.5 * m_e * dot(v, v));
    total_time += time_to_collision;
}

/*
 * Generates a set of simulations, saving the files to py/studies/simulation-runs/
 * The files generated are plain text, with each row being a simulation step. The columns
 * are...
 * 
 * time (in ns), x-position (in um), y-position (in um), z-position (in um), energy (in eV), drift velocity (in m/s), number of ionizations
 * 
 * @param volts The bulk field (in V/cm)
 * @param cutoff The cutoff time (in s) (this is superseded by hit checks when using a non-uniform field)
 * @param cores The number of available processor cores (used in naming files)
 * @param write_every Save one simulated step for every write_every steps
 * @param k The core that the simulation runs on (used in naming files)
 * @param batches The number of batches to be generated
 * @param bar The progress bar to be used
 */
void generate_plot(int volts, double cutoff, int cores, int write_every, int k, int batches, ProgressBar& bar)
{
    for (int i = 0; i < batches; i++) {
        bar.new_file(cores * i + (k + 1), k);
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());

        int total_ions = 0;
        Electron elec(total_ions, volts, starting_pos(generator), random_velocity(generator, true), generator);

        double starting_x = 0;
        if (!uniform_field)
            starting_x = elec.position().x;

        std::ofstream file("../py/simulation-runs/" + std::to_string(volts) + "V - " + std::to_string(cores * i + (k + 1)) + ".txt");
        assert(file.is_open());
        
        int simulation_step = 1;

        while (elec.elapsed_time() < cutoff) {
            elec.update();

            if (!uniform_field && hit_check(elec.position()))
                break;

            if (simulation_step++ % write_every != 0)
                continue;

            if (uniform_field)
                bar.update(elec.elapsed_time() / cutoff, k);
            else
                bar.update(1 - (elec.position().x - z_min) / ((z_max - z_min) * spawn_height_scale), k);

            if (bar.min_prog(k))
                bar.display();
            
            file << elec.elapsed_time() * 1e9 << "," << elec.position().x * 1e6 << "," << elec.position().y * 1e6 << "," << elec.position().z * 1e6 << "," << elec.ke() << "," << ((elec.position().x - starting_x) / elec.elapsed_time()) << "," << total_ions << "\n";
        }

        elec.update();
        bar.update(1, k);

        file << elec.elapsed_time() * 1e9 << "," << elec.position().x * 1e6 << "," << elec.position().y * 1e6 << "," << elec.position().z * 1e6 << "," << elec.ke() << "," << ((elec.position().x - starting_x) / elec.elapsed_time()) << "," << total_ions << "\n";
        file.close();
        assert(!file.is_open());
    }

    return;
}