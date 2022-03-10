#include "renderer.hpp"

int main() {
    Renderer renderer(480, 360);
    renderer.SetDrawColor(Color4{1, 0, 0, 1});
    renderer.SetClearColor(Color4{0.2, 0.2, 0.2, 1});
    renderer.Clear();
    for (int i = 100; i < 400; i++) {
        renderer.DrawPixel(i, 180);
    }
    // renderer.DrawLine(10, 10, 450, 300);
    renderer.Save("test_renderer.bmp");
    return 0;
}
