#include "../include/main.hpp"
#include "../include/particle.hpp"

// constexpr auto use_avx2 = true;
constexpr U64 particles_count = 1536;
constexpr U64 width = 1280;
constexpr U64 height = 720;

int main(int argc, char **argv) {
  // Fail safe
  int status = EXIT_FAILURE;
  try {

    InitWindow(width, height, "nbody2024 " nbody_VSTRING_FULL);
    SetTargetFPS(60);
    std::random_device rd;
    S64 simulation_time = 0;
#if 0
    if (use_avx2) {
      auto &&avx2_system =
          nbody::system_create_vectorized(rd(), particles_count / 8, width, height);
      F32 x[8]{0.0f};
      F32 y[8]{0.0f};
      F32 mass[8]{0.0f};
      while (!WindowShouldClose()) {
        auto t0 = std::chrono::steady_clock::now();
        nbody::system_tick_vectorized(std::thread::hardware_concurrency());
        auto t1 = std::chrono::steady_clock::now();

        BeginDrawing();
        ClearBackground(BLACK);
        for (const nbody::particle_cluster8 &c : *avx2_system) {
          _mm256_storeu_ps(x, c.x_position);
          _mm256_storeu_ps(y, c.y_position);
          _mm256_storeu_ps(mass, c.mass);

          for (auto i = 0; i < 8; i++) {
            DrawCircleLinesV(
                Vector2{.x = x[i] + (width / 2), .y = y[i] + (height / 2)},
                mass[i] / 2,
                Color{.r = 255,
                      .g = static_cast<U8>(256 - (mass[i] * 4.0f)),
                      .b = 0,
                      .a = 255});
          }
        }
        DrawFPS(20, 20);
        std::stringstream timing_stream{};
        timing_stream << "SIM ";
        timing_stream << std::chrono::duration_cast<std::chrono::milliseconds>(
            t1 - t0);
        DrawText(timing_stream.str().c_str(), 20, 50, 20, GRAY);

        EndDrawing();
      }
    } else {
#endif
    auto &&system = nbody::system_create(rd(), particles_count, width, height);

    while (!WindowShouldClose()) {
      auto t0 = std::chrono::steady_clock::now();
      nbody::system_tick(std::thread::hardware_concurrency());
      auto t1 = std::chrono::steady_clock::now();

      BeginDrawing();
      ClearBackground(BLACK);
      for (const nbody::particle &p : *system) {
        DrawCircleLinesV(Vector2{.x = p.x_position + (width / 2),
                                 .y = p.y_position + (height / 2)},
                         p.mass / 2,
                         Color{.r = 255,
                               .g = static_cast<U8>(256 - (p.mass * 4.0f)),
                               .b = 0,
                               .a = 255});
      }
      DrawFPS(20, 20);
      std::stringstream timing_stream{};
      timing_stream << "SIM ";
      timing_stream << std::chrono::duration_cast<std::chrono::milliseconds>(
          t1 - t0);
      DrawText(timing_stream.str().c_str(), 20, 50, 20, GRAY);

      EndDrawing();
    }
    // }
    // Success is assumed when everything finishes with no exceptions
    status = EXIT_SUCCESS;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
  }
  return status;
}