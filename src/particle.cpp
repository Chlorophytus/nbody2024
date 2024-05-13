#include "../include/particle.hpp"
using namespace nbody;

static std::unique_ptr<std::vector<particle>> particles{nullptr};
constexpr static F32 G_force = -0.0001f;
static U64 width;
static U64 height;

std::unique_ptr<std::vector<particle>> &nbody::system_get() {
  return particles;
}

void nbody::system_create(U32 seed, U64 num_particles, U64 x, U64 y) {
  std::cerr << "with seed " << seed << std::endl;
  width = x;
  height = y;
  std::mt19937 randomizer(seed);

  particles = std::make_unique<std::vector<particle>>();
  particles->resize(num_particles);
  // Generates the position ordinates
  std::uniform_real_distribution gen_ordinate(-0.5f, 0.5f);
  // Acceleration and velocity default to 0.0
  // Mass is an exponential distribution
  std::exponential_distribution gen_mass(0.125f);

  U32 id = 0;
  for (particle &p : *particles) {
    p.x_position = width * gen_ordinate(randomizer);
    p.y_position = height * gen_ordinate(randomizer);
    p.x_velocity = 0.0f;
    p.y_velocity = 0.0f;
    p.x_acceleration = 0.0f;
    p.y_acceleration = 0.0f;
    p.mass = std::clamp(gen_mass(randomizer), 2.0f, 64.0f);
    p.id = id++;
  }
}

void nbody::particle::gravitate(particle &other) {
  if (id == other.id) {
    return;
  }
  // calculate G forces
  const F32 x_delta = std::sqrt(x_position - other.x_position);
  const F32 y_delta = std::sqrt(y_position - other.y_position);
  const F32 x_delta_squared = std::pow(x_delta, 2.0f);
  const F32 y_delta_squared = std::pow(y_delta, 2.0f);
  const F32 m1_m2 = other.mass * mass;
  const F32 x_force = G_force * (m1_m2 / std::max(512.0f, x_delta_squared));
  const F32 y_force = G_force * (m1_m2 / std::max(512.0f, y_delta_squared));

  // Newton's second law, F = ma. In this case, a = (F/m) based on solving the
  // literal equation.
  x_acceleration -= x_force / mass;
  y_acceleration -= y_force / mass;
  other.x_acceleration += x_force / other.mass;
  other.y_acceleration += y_force / other.mass;
}
void nbody::particle::iterate() {
  // Acceleration is delta velocity
  x_velocity += x_acceleration;
  y_velocity += y_acceleration;

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
  std::atomic<U64> thread_semaphore(num_threads);
  std::mutex particles_access;

  for (U64 thread_i = 0; thread_i < num_threads; thread_i++) {
    std::thread{[&thread_semaphore, &particles_access, &num_threads,
                 thread_i]() {
      thread_local std::vector<particle> particles_stride{};
      {
        std::lock_guard<std::mutex> lock(particles_access);
        const U64 particles_size = particles->size();
        for (U64 particle_i = thread_i; particle_i < particles_size;
             particle_i += num_threads) {
          particles_stride.emplace_back((*particles)[particle_i]);
        }
      }
      for (particle &callee : particles_stride) {
        std::lock_guard<std::mutex> lock(particles_access);
        std::for_each(particles->begin(), particles->end(),
                      [&callee](particle &p) { callee.gravitate(p); });
      }
      std::for_each(particles_stride.begin(), particles_stride.end(),
                    [](particle &p) { p.iterate(); });
      {
        // std::lock_guard<std::mutex> lock(particles_access);
        const U64 particles_size = particles->size();
        U64 stride_i = 0;
        for (U64 particle_i = thread_i; particle_i < particles_size;
             particle_i += num_threads) {
          (*particles)[particle_i] = particles_stride[stride_i];
          stride_i++;
        }
      }
      thread_semaphore--;
    }}.detach();
  }
  do {
    std::this_thread::sleep_for(std::chrono::microseconds(1));
  } while (thread_semaphore > 0);
}