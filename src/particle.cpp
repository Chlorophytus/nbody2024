#include "../include/particle.hpp"
#include "raymath.h"
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
    p.position = Vector2{.x = width * gen_ordinate(randomizer), .y = height * gen_ordinate(randomizer)};
    p.velocity = Vector2Zero();
    p.acceleration = Vector2Zero();
    p.mass = std::clamp(gen_mass(randomizer), 4.0f, 64.0f);
    p.color = Color{.r = 255,
                    .g = static_cast<U8>(256 - (p.mass * 4.0f)),
                    .b = 0,
                    .a = 255};
    p.id = id++;
  }

  return particles;
}

void nbody::particle::gravitate(particle &other) {
  // Don't gravitate to ourselves
  if (id == other.id) {
    return;
  }

  // Calculate the difference of the positions
  Vector2 delta = other.position - position;

  // When calculating, clamp the absolute values to a certain amount so a
  // division by 0 does not occur.
  if (std::abs(delta.x) < minimum) {
    delta.x = std::signbit(delta.x) ? -minimum : minimum;
  }
  if (std::abs(delta.y) < minimum) {
    delta.y = std::signbit(delta.y) ? -minimum : minimum;
  }

  // One component of Newton's gravitation formula is to multiply both masses
  const Vector2 m1_m2 = Vector2{.x = other.mass * mass, .y = other.mass * mass};

  // Calculate the final Newton's gravitation formula values

  // Gravitation force value
  Vector2 force = (m1_m2 / (delta * delta)) * G_force;

  // Squaring the distances removed the signs, let's bring the signedness back
  if (std::signbit(delta.x)) {
    force.x *= -1;
  }
  if (std::signbit(delta.y)) {
    force.y *= -1;
  }
  // Newton's second law, F = ma. In this case, a = (F/m) based on solving the
  // literal equation.

  // Self acceleration is attraction in one direction
  acceleration += force / mass;

  // Neighbor acceleration is attraction in other direction
  other.acceleration -= force / other.mass;
}
void nbody::particle::iterate(const U64 particles_size) {
  // acceleration is a change in velocity
  velocity += acceleration / particles_size;

  // Reset acceleration accumulator
  acceleration = Vector2Zero();

  // velocity can bounce off the window edges
  if (std::abs(position.x) > (width / 2.0f)) {
    velocity.x *= -1.0f;
  }
  // Y velocity can bounce off the window edges
  if (std::abs(position.y) > (height / 2.0f)) {
    velocity.y *= -1.0f;
  }

  // Y velocity is a change in Y position
  position += velocity;
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
