#pragma once
#include "main.hpp"

namespace nbody {
struct particle {
  F32 x_position;
  F32 y_position;
  F32 x_velocity;
  F32 y_velocity;
  F32 x_acceleration;
  F32 y_acceleration;
  F32 mass;
  U32 id;

  void gravitate(particle &);
  void iterate();
};

void system_create(U32, U64, U64, U64);
void system_tick(const U64);
void system_tick_vectorized(const U64);
std::unique_ptr<std::vector<particle>> &system_get();
} // namespace nbody
