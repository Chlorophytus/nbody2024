#pragma once
#include "main.hpp"

namespace nbody {
/// Store 1 particle
struct particle {
  F32 x_position;
  F32 y_position;
  F32 x_velocity;
  F32 y_velocity;
  F32 x_acceleration;
  F32 y_acceleration;
  F32 mass;
  U32 id;

  /// @brief Calculate gravity of one particle to another
  /// @param other Another particle to gravitate toward
  void gravitate(particle &);

  /// @brief Calculate positions and velocities
  /// @param num_particles Number of particles in this system, to divide the
  /// acceleration accumulator
  void iterate(const U64);
};

/// @brief Initialize an N-body simulation
/// @param seed Random Number Generator seed value
/// @param num_particles Particle count
/// @param x Simulation width
/// @param y Simulation height
/// @return N-body simulation singleton
std::unique_ptr<std::vector<particle>> &system_create(U32, U64, U64, U64);

/// @brief Initialize a *vectorized* N-body simulation
/// @param seed Random Number Generator seed value
/// @param num_particles 8-particle cluster count
/// @param x Simulation width
/// @param y Simulation height
/// @return N-body simulation singleton
// std::unique_ptr<std::vector<particle_cluster8>> &
//     system_create_vectorized(U32, U64, U64, U64);

/// @brief Use regular non-vectorized calculations for the N-body simulation
/// @param num_threads Worker thread count
void system_tick(const U64);

/// @brief Use AVX2 and FMA3 intrinsics to calculate the N-body simulation
/// @param num_threads Worker thread count
// void system_tick_vectorized(const U64);
} // namespace nbody
