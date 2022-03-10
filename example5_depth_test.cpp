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
        {Vec4{0.5, 0.5, -1, 1}, Vec4{1, 0, 0, 1}},
        {Vec4{0.5, -0.5, -1, 1}, Vec4{0, 1, 0, 1}},
        {Vec4{-0.5, -0.5, -1, 1}, Vec4{0, 0, 1, 1}},
    };

    // auto orthoMat = CreateOrtho(-1, 1, 1, -1, 2, 0);
    auto orthoMat = CreatePersp(M_PI * 0.5, 480.0f/320.0f, -0.1, -100);

    renderer.SetVertexShader([&](int index, ShaderContext& output) {
        output.varyingVec4[Color] = vs_input[index].color;
        return orthoMat * vs_input[index].pos;        
    });

    renderer.SetFragmentShader([&](ShaderContext& input) {
        return input.varyingVec4[Color];
    });

    renderer.DrawPrimitive();

    renderer.SetVertexShader([&](int index, ShaderContext& output) {
        output.varyingVec4[Color] = vs_input[index].color;
        auto pos = vs_input[index].pos;
        pos.x -= 0.3;
        pos.z = -2;
        return orthoMat * pos;
    });

    renderer.DrawPrimitive();

    renderer.Save("depth_test.bmp");
    return 0;
}
