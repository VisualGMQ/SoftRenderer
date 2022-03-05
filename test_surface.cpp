#include "renderer.hpp"

int main() {
    Surface surface(480, 320);
    surface.Clear(Color4{0.5, 0.5, 0.5, 1});

    for (int i = 0 ; i < 400; i++) {
        surface.PutPixel(i+ 40, 160, Color4{1, 0, 0, 1});
    }

    surface.Save("test_surface.bmp");
    return 0;
}
