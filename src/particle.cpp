#include "../include/particle.hpp"
using namespace nbody;

static std::unique_ptr<std::vector<particle>> particles{nullptr};
// static std::unique_ptr<std::vector<particle_cluster8>> clusters{nullptr};
constexpr static F32 G_force = -8.0f;
constexpr static F32 minimum = 2.0f;
static U64 width;
static U64 height;
// const static __m256 minimum_vec_hi = _mm256_set1_ps(minimum);
// const static __m256 minimum_vec_lo = _mm256_set1_ps(-minimum);
// const static __m256 G_force_vec = _mm256_set1_ps(G_force);

// ============================================================================
// PARTICLE NON VECTORIZED
// ============================================================================
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
#if 0
// ============================================================================
// PARTICLE VECTORIZED
// ============================================================================
decltype(clusters) &nbody::system_create_vectorized(U32 seed, U64 num_particles,
                                                    U64 x, U64 y) {
  if (particles || clusters) {
    std::cerr << "Already initialized the system" << std::endl;
    return clusters;
  }
  std::cerr << "Vectorized with seed " << seed << std::endl;
  std::cerr << "With " << num_particles << " clusters of 8 particles"
            << std::endl;
  width = x;
  height = y;
  std::mt19937 randomizer(seed);

  clusters = std::make_unique<std::vector<particle_cluster8>>();
  clusters->resize(num_particles);

  // Generates the position ordinates
  std::uniform_real_distribution gen_ordinate(-0.5f, 0.5f);
  // Acceleration and velocity default to 0.0
  // Mass is an exponential distribution
  std::exponential_distribution gen_mass(0.1f);

  U32 id = 0;
  for (particle_cluster8 &c : *clusters) {
    c.x_position = _mm256_set_ps(
        width * gen_ordinate(randomizer), width * gen_ordinate(randomizer),
        width * gen_ordinate(randomizer), width * gen_ordinate(randomizer),
        width * gen_ordinate(randomizer), width * gen_ordinate(randomizer),
        width * gen_ordinate(randomizer), width * gen_ordinate(randomizer));
    c.y_position = _mm256_set_ps(
        height * gen_ordinate(randomizer), height * gen_ordinate(randomizer),
        height * gen_ordinate(randomizer), height * gen_ordinate(randomizer),
        height * gen_ordinate(randomizer), height * gen_ordinate(randomizer),
        height * gen_ordinate(randomizer), height * gen_ordinate(randomizer));
    c.x_velocity = _mm256_setzero_ps();
    c.y_velocity = _mm256_setzero_ps();
    c.x_acceleration = _mm256_setzero_ps();
    c.y_acceleration = _mm256_setzero_ps();
    c.mass = _mm256_set_ps(std::clamp(gen_mass(randomizer), 4.0f, 64.0f),
                           std::clamp(gen_mass(randomizer), 4.0f, 64.0f),
                           std::clamp(gen_mass(randomizer), 4.0f, 64.0f),
                           std::clamp(gen_mass(randomizer), 4.0f, 64.0f),
                           std::clamp(gen_mass(randomizer), 4.0f, 64.0f),
                           std::clamp(gen_mass(randomizer), 4.0f, 64.0f),
                           std::clamp(gen_mass(randomizer), 4.0f, 64.0f),
                           std::clamp(gen_mass(randomizer), 4.0f, 64.0f));
    c.id = id++;
  }

  return clusters;
}

void nbody::system_tick_vectorized(const U64 num_threads) {
  const U64 particles_size = clusters->size();
  std::atomic<U64> thread_semaphore(num_threads);

  for (U64 thread_i = 0; thread_i < num_threads; thread_i++) {
    std::thread{[&thread_semaphore, &num_threads, thread_i, particles_size]() {
      thread_local std::vector<particle_cluster8> particles_stride{};
      // Emplace particles into threaded particle cache
      for (U64 particle_i = thread_i; particle_i < particles_size;
           particle_i += num_threads) {
        particles_stride.emplace_back((*clusters)[particle_i]);
      }
      // Calculate the pull for every particle
      for (particle_cluster8 &callee : particles_stride) {
        std::for_each(clusters->begin(), clusters->end(),
                      [&callee](particle_cluster8 &p) {
                        callee.gravitate_vectorized(p);
                      });
      }
      // Calculate velocity and position then
      std::for_each(particles_stride.begin(), particles_stride.end(),
                    [particles_size](particle_cluster8 &p) {
                      p.iterate_vectorized(particles_size);
                    });
      // Strided sync threaded particle cache into particles
      U64 stride_i = 0;
      for (U64 particle_i = thread_i; particle_i < particles_size;
           particle_i += num_threads) {
        (*clusters)[particle_i] = particles_stride[stride_i];
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
void nbody::particle_cluster8::gravitate_vectorized(particle_cluster8 &other) {
  if (id == other.id) {
    return;
  }
  // Unroller dummy numbers, can be any float that isn't zero.
  const __m256 test_unroll[8]{
      _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f),
      _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f),
      _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f),
      _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f),
      _mm256_set_ps(0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f),
      _mm256_set_ps(0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f),
      _mm256_set_ps(0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f),
      _mm256_set_ps(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f)};

  // X deltas
  __m256 x_delta = _mm256_sub_ps(other.x_position, x_position);
  const __m256 x_compare_hi =
      _mm256_cmp_ps(x_delta, minimum_vec_hi, _CMP_LT_OQ);
  const __m256 x_compare_lo =
      _mm256_cmp_ps(x_delta, minimum_vec_lo, _CMP_GT_OQ);

  // Y deltas
  __m256 y_delta = _mm256_sub_ps(other.y_position, y_position);
  const __m256 y_compare_hi =
      _mm256_cmp_ps(x_delta, minimum_vec_hi, _CMP_LT_OQ);
  const __m256 y_compare_lo =
      _mm256_cmp_ps(y_delta, minimum_vec_lo, _CMP_GT_OQ);

  // NOTE: Comparisons will be inverted
  F32 x_delta_clamped[8]{0.0f};
  F32 y_delta_clamped[8]{0.0f};

  _mm256_storeu_ps(x_delta_clamped, x_delta);
  _mm256_storeu_ps(y_delta_clamped, y_delta);

  // Two birds with one stone
  // HACK: vectorization bottleneck!
  for (auto i = 0; i < 8; i++) {
    // Test is inverted, compare high into bit 1
    U32 x_compare = _mm256_testz_ps(x_compare_hi, test_unroll[i]);
    x_compare <<= 1;
    // Compare low into bit 0
    x_compare = _mm256_testz_ps(x_compare_lo, test_unroll[i]);

    U32 y_compare = _mm256_testz_ps(y_compare_hi, test_unroll[i]);
    y_compare <<= 1;
    y_compare = _mm256_testz_ps(y_compare_lo, test_unroll[i]);

    switch (x_compare) {
    case 0: {
      // Do nothing
      break;
    }
    case 1 << 0: {
      // Hit
      x_delta_clamped[i] = -minimum;
      break;
    }
    case 1 << 1: {
      // Hit
      x_delta_clamped[i] = minimum;
      break;
    }
    default: {
      throw std::runtime_error{"X compare domain"};
    }
    }
    switch (y_compare) {
    case 0: {
      break;
    }
    case 1 << 0: {
      y_delta_clamped[i] = -minimum;
      break;
    }
    case 1 << 1: {
      y_delta_clamped[i] = minimum;
      break;
    }
    default: {
      throw std::runtime_error{"Y compare domain"};
    }
    }
    // Negative?
    const auto x_negative = std::signbit(x_delta_clamped[i]);
    const auto y_negative = std::signbit(y_delta_clamped[i]);

    // Reusing numbers here...
    x_delta_clamped[i] *= x_delta_clamped[i] * (x_negative ? -1 : 1);
    y_delta_clamped[i] *= y_delta_clamped[i] * (y_negative ? -1 : 1);
  }
  x_delta = _mm256_loadu_ps(x_delta_clamped);
  y_delta = _mm256_loadu_ps(y_delta_clamped);

  const __m256 m1_m2 = _mm256_mul_ps(other.mass, mass);

  // calculate G forces
  // NOTE: We squared the deltas beforehand
  __m256 x_force = _mm256_mul_ps(G_force_vec, _mm256_div_ps(m1_m2, x_delta));
  __m256 y_force = _mm256_mul_ps(G_force_vec, _mm256_div_ps(m1_m2, y_delta));

  // Newton's second law, F = ma. In this case, a = (F/m) based on solving the
  // literal equation.
  x_acceleration = _mm256_sub_ps(x_acceleration, _mm256_div_ps(x_force, mass));
  y_acceleration = _mm256_sub_ps(y_acceleration, _mm256_div_ps(y_force, mass));
  other.x_acceleration =
      _mm256_add_ps(other.x_acceleration, _mm256_div_ps(x_force, other.mass));
  other.y_acceleration =
      _mm256_add_ps(other.y_acceleration, _mm256_div_ps(y_force, other.mass));
}

void nbody::particle_cluster8::iterate_vectorized(const U64 particles_size) {
  // Acceleration is delta velocity
  x_velocity = _mm256_add_ps(
      x_velocity,
      _mm256_div_ps(x_acceleration, _mm256_set1_ps(particles_size * 8)));
  y_velocity = _mm256_add_ps(
      y_velocity,
      _mm256_div_ps(y_acceleration, _mm256_set1_ps(particles_size * 8)));
  // Reset acceleration accumulators
  x_acceleration = _mm256_setzero_ps();
  y_acceleration = _mm256_setzero_ps();

  F32 x_positions[8]{0.0f};
  F32 y_positions[8]{0.0f};
  F32 x_multiply_velocities[8]{1.0f};
  F32 y_multiply_velocities[8]{1.0f};

  _mm256_storeu_ps(x_positions, x_position);
  _mm256_storeu_ps(y_positions, y_position);

  for(auto i = 0; i < 8; i++) {
    if(std::abs(x_positions[i]) > (width / 2.0f)) {
      x_multiply_velocities[i] = -1.0f;
    }
    if(std::abs(y_positions[i]) > (height / 2.0f)) {
      y_multiply_velocities[i] = -1.0f;
    }
  }

  x_velocity = _mm256_mul_ps(x_velocity, _mm256_loadu_ps(x_multiply_velocities));
  y_velocity = _mm256_mul_ps(y_velocity, _mm256_loadu_ps(y_multiply_velocities));

  x_position = _mm256_add_ps(x_position, x_velocity);
  y_position = _mm256_add_ps(y_position, y_velocity);
}
#endif