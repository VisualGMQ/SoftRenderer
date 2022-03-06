#include "renderer.hpp"

int main() {
    Renderer renderer(480, 360);
    renderer.SetDrawColor(Color4{1, 0, 0, 1});
    renderer.SetClearColor(Color4{0.1, 0.1, 0.1, 1});
    renderer.Clear();
    renderer.DrawLine(10, 10, 450, 300);
    renderer.Save("test_renderer.bmp");
    return 0;
}
