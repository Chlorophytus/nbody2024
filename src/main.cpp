#include "../include/main.hpp"
#include "../include/particle.hpp"

constexpr U64 particles_count = 768;
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
      std::string particle_text{};
      timing_stream << "Simulation ";
      timing_stream << std::chrono::duration_cast<std::chrono::milliseconds>(
          t1 - t0);
      particle_text = std::to_string(particles_count) + " particles";
      DrawText(timing_stream.str().c_str(), 20, 50, 20, GRAY);
      DrawText(particle_text.c_str(), 20, 80, 20, GRAY);

      EndDrawing();
    }
    // Success is assumed when everything finishes with no exceptions
    status = EXIT_SUCCESS;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
  }
  return status;
}
