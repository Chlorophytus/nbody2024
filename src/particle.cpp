#include "../include/particle.hpp"
using namespace nbody;

static std::unique_ptr<std::vector<particle>> particles{nullptr};
// static std::unique_ptr<std::vector<particle_cluster8>> clusters{nullptr};
constexpr static F32 G_force = -8.0f;
constexpr static F32 minimum = 2.0f;
static U64 width;
static U64 height;

decltype(particles) &nbody::system_create(U32 seed, U64 num_particles, U64 x,
                                          U64 y) {
  // if (particles || clusters) {
  //   std::cerr << "Already initialized the system" << std::endl;
  //   return particles;
  // }
  std::cerr << "Non-vectorized with seed " << seed << std::endl;
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
  std::exponential_distribution gen_mass(0.1f);

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
  if (id == other.id) {
    return;
  }
  F32 x_delta = other.x_position - x_position;
  if (std::abs(x_delta) < minimum) {
    x_delta = std::signbit(x_delta) ? -minimum : minimum;
  }

  F32 y_delta = other.y_position - y_position;
  if (std::abs(y_delta) < minimum) {
    y_delta = std::signbit(y_delta) ? -minimum : minimum;
  }
  const F32 m1_m2 = other.mass * mass;

  // calculate G forces
  F32 x_force = G_force * (m1_m2 / (x_delta * x_delta));
  F32 y_force = G_force * (m1_m2 / (y_delta * y_delta));
  // squaring removes the negativity in the magnitude
  if (std::signbit(x_delta)) {
    x_force *= -1;
  }
  if (std::signbit(y_delta)) {
    y_force *= -1;
  }

  // Newton's second law, F = ma. In this case, a = (F/m) based on solving the
  // literal equation.
  x_acceleration -= x_force / mass;
  y_acceleration -= y_force / mass;
  other.x_acceleration += x_force / other.mass;
  other.y_acceleration += y_force / other.mass;
}
void nbody::particle::iterate(const U64 particles_size) {
  // Acceleration is delta velocity
  x_velocity += x_acceleration / particles_size;
  y_velocity += y_acceleration / particles_size;
  // Reset acceleration accumulators
  x_acceleration = 0.0f;
  y_acceleration = 0.0f;

  if (std::abs(x_position) > (width / 2.0f)) {
    x_velocity *= -1.0f;
  }
  x_position += x_velocity;

  if (std::abs(y_position) > (height / 2.0f)) {
    y_velocity *= -1.0f;
  }
  y_position += y_velocity;
}

void nbody::system_tick(const U64 num_threads) {
  const U64 particles_size = particles->size();
  std::atomic<U64> thread_semaphore(num_threads);
  for (U64 thread_i = 0; thread_i < num_threads; thread_i++) {
    std::thread{[&thread_semaphore, &num_threads, thread_i, particles_size]() {
      thread_local std::vector<particle> particles_stride{};
      // Emplace particles into threaded particle cache
      for (U64 particle_i = thread_i; particle_i < particles_size;
           particle_i += num_threads) {
        particles_stride.emplace_back((*particles)[particle_i]);
      }
      // Calculate the pull for every particle
      for (particle &callee : particles_stride) {
        std::for_each(particles->begin(), particles->end(),
                      [&callee](particle &p) { callee.gravitate(p); });
      }
      // Calculate velocity and position then
      std::for_each(
          particles_stride.begin(), particles_stride.end(),
          [particles_size](particle &p) { p.iterate(particles_size); });
      // Strided sync threaded particle cache into particles
      U64 stride_i = 0;
      for (U64 particle_i = thread_i; particle_i < particles_size;
           particle_i += num_threads) {
        (*particles)[particle_i] = particles_stride[stride_i];
        stride_i++;
      }
      // Decrement semaphore
      thread_semaphore--;
    }}.detach();
  }
  // Spin until semaphore is 0
  do {
    std::this_thread::sleep_for(std::chrono::microseconds(1));
  } while (thread_semaphore > 0);
}