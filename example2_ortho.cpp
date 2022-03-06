#include "renderer.hpp"

enum UniformVar {
    Color = 0,
};

int main() {
    Renderer renderer(480, 320);
    renderer.SetClearColor(Color4{0.1, 0.1, 0.1, 1});
    renderer.Clear();
    renderer.SetViewport(0, 0, 480, 320);

    struct { Vec4 pos; Vec4 color; } vs_input[3] = {
        {Vec4{80, 20, 0, 1}, Vec4{1, 0, 0, 1}},
        {Vec4{400, 20, 0, 1}, Vec4{0, 1, 0, 1}},
        {Vec4{240, 300, 0, 1}, Vec4{0, 0, 1, 1}},
    };

    auto orthoMat = CreateOrtho(0, 480, 320, 0, 1, -1);

    renderer.SetVertexShader([&](int index, ShaderContext& output) {
        output.varyingVec4[Color] = vs_input[index].color;
        return orthoMat * vs_input[index].pos;        
    });

    renderer.SetFragmentShader([&](ShaderContext& input) {
        return input.varyingVec4[Color];
    });

    renderer.DrawPrimitive();
    renderer.Save("ortho_triangle.bmp");
    return 0;
}
