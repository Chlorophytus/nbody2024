#include "../include/particle.hpp"
using namespace nbody;

static std::unique_ptr<std::vector<particle>> particles{nullptr};
static U64 width;
static U64 height;

// This is a gravitational coefficient
// Positive G-force attracts particles, but negative G-force repels them
constexpr static F32 G_force = 8.0f;

// Clamp minimum distance
constexpr static F32 minimum = 8.0f;

decltype(particles) &nbody::system_create(U32 seed, U64 num_particles, U64 x,
                                          U64 y) {
  std::cerr << "With seed " << seed << std::endl;
  std::cerr << "With " << num_particles << " particles" << std::endl;
  width = x;
  height = y;
  std::mt19937 randomizer(seed);

  particles = std::make_unique<std::vector<particle>>();
  particles->resize(num_particles);
  // Generates the position ordinates
  std::uniform_real_distribution gen_ordinate(-0.5f, 0.5f);
  // Acceleration and velocity default to 0.0
  // Mass is an exponential distribution
  std::exponential_distribution gen_mass(0.15f);

  U32 id = 0;
  for (particle &p : *particles) {
    p.x_position = width * gen_ordinate(randomizer);
    p.y_position = height * gen_ordinate(randomizer);
    p.x_velocity = 0.0f;
    p.y_velocity = 0.0f;
    p.x_acceleration = 0.0f;
    p.y_acceleration = 0.0f;
    p.mass = std::clamp(gen_mass(randomizer), 4.0f, 64.0f);
    p.id = id++;
  }

  return particles;
}

void nbody::particle::gravitate(particle &other) {
  // Don't gravitate to ourselves
  if (id == other.id) {
    return;
  }

  // Calculate the difference of the X positions
  F32 x_delta = other.x_position - x_position;

  // Calculate the difference of the Y positions
  F32 y_delta = other.y_position - y_position;

  // When calculating, clamp the absolute values to a certain amount so a
  // division by 0 does not occur.
  if (std::abs(x_delta) < minimum) {
    x_delta = std::signbit(x_delta) ? -minimum : minimum;
  }
  if (std::abs(y_delta) < minimum) {
    y_delta = std::signbit(y_delta) ? -minimum : minimum;
  }

  // One component of Newton's gravitation formula is to multiply both masses
  const F32 m1_m2 = other.mass * mass;

  // Calculate the final Newton's gravitation formula values

  // X gravitation force value
  F32 x_force = G_force * (m1_m2 / (x_delta * x_delta));
  // Y gravitation force value
  F32 y_force = G_force * (m1_m2 / (y_delta * y_delta));

  // Squaring the distances removed the signs, let's bring the signedness back
  if (std::signbit(x_delta)) {
    x_force *= -1;
  }
  if (std::signbit(y_delta)) {
    y_force *= -1;
  }

  // Newton's second law, F = ma. In this case, a = (F/m) based on solving the
  // literal equation.

  // Self X acceleration is attraction in one direction
  x_acceleration += x_force / mass;
  // Self Y acceleration is attraction in one direction
  y_acceleration += y_force / mass;

  // Neighbor X acceleration is attraction in other direction
  other.x_acceleration -= x_force / other.mass;
  // Neighbor Y acceleration is attraction in other direction
  other.y_acceleration -= y_force / other.mass;
}
void nbody::particle::iterate(const U64 particles_size) {
  // X acceleration is a change in X velocity
  x_velocity += x_acceleration / particles_size;
  // Y acceleration is a change in Y velocity
  y_velocity += y_acceleration / particles_size;

  // Reset X acceleration accumulator
  x_acceleration = 0.0f;
  // Reset Y acceleration accumulator
  y_acceleration = 0.0f;

  // X velocity can bounces off the window edges
  if (std::abs(x_position) > (width / 2.0f)) {
    x_velocity *= -1.0f;
  }

  // X velocity is a change in X position
  x_position += x_velocity;

  // Y velocity can bounce off the window edges
  if (std::abs(y_position) > (height / 2.0f)) {
    y_velocity *= -1.0f;
  }

  // Y velocity is a change in Y position
  y_position += y_velocity;
}

void nbody::system_tick(const U64 num_threads) {
  const U64 particles_size = particles->size();

  // A semaphore that makes sure all threads complete before the next iteration
  std::atomic<U64> thread_semaphore(num_threads);

  // Loop that dispatches threads to each core
  for (U64 thread_i = 0; thread_i < num_threads; thread_i++) {
    // One thread handles one region of particles and detaches
    std::thread{[&thread_semaphore, &num_threads, thread_i, particles_size]() {
      // Particles cache
      std::vector<particle> particles_stride{};

      // Emplace particles into particle cache
      for (U64 particle_i = thread_i; particle_i < particles_size;
           particle_i += num_threads) {
        particles_stride.emplace_back((*particles)[particle_i]);
      }

      // Calculate the accelerations for every particle
      for (particle &callee : particles_stride) {
        std::for_each(particles->begin(), particles->end(),
                      [&callee](particle &p) { callee.gravitate(p); });
      }

      // Calculate velocities and positions then
      std::for_each(
          particles_stride.begin(), particles_stride.end(),
          [particles_size](particle &p) { p.iterate(particles_size); });

      // Strided synchronize threaded particle cache into particles
      U64 stride_i = 0;
      for (U64 particle_i = thread_i; particle_i < particles_size;
           particle_i += num_threads) {
        (*particles)[particle_i] = particles_stride[stride_i];
        stride_i++;
      }
      // Decrement semaphore, we're done on this end
      thread_semaphore--;
    }}.detach();
  }

  // Spin until semaphore is 0, which signals all threads did their jobs
  do {
    std::this_thread::sleep_for(std::chrono::microseconds(1));
  } while (thread_semaphore > 0);
}
