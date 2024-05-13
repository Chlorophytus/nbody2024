#include "../include/main.hpp"
#include "../include/particle.hpp"

int main(int argc, char **argv) {
  // Fail safe
  int status = EXIT_FAILURE;
  try {
    constexpr U64 width = 1280;
    constexpr U64 height = 720;

    InitWindow(width, height, "nbody2024 " nbody_VSTRING_FULL);
    SetTargetFPS(60);
    std::random_device rd;

    nbody::system_create(rd(), 1536, width, height);
    auto &&system = nbody::system_get();

    while (!WindowShouldClose()) {
      BeginDrawing();
      ClearBackground(BLACK);
      for (const nbody::particle &p : *system) {
        DrawCircleLinesV(Vector2{.x = p.x_position + (width / 2),
                            .y = p.y_position + (height / 2)},
                    p.mass / 2,
                    Color{.r = 255, .g = static_cast<U8>(256 - (p.mass * 4.0f)), .b = 0, .a = 255});
      }
      DrawFPS(20, 20);
      EndDrawing();
      auto t0 = std::chrono::steady_clock::now();
      nbody::system_tick(8);
      auto t1 = std::chrono::steady_clock::now();
      std::cerr << "Simulation time is " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0) << std::endl;
    }
    // Success is assumed when everything finishes with no exceptions
    status = EXIT_SUCCESS;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
  }
  return status;
}