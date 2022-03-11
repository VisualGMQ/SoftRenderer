#include "renderer.hpp"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

enum UniformVar {
    Color = 0,
};

int main() {
    glm::ortho(0, 480, 320, 0);
    Renderer renderer(480, 320);
    renderer.SetClearColor(Color4{0.1, 0.1, 0.1, 1});
    renderer.Clear();
    renderer.SetViewport(0, 0, 480, 320);
    renderer.SetFaceCull(CW);

    struct { Vec4 pos; Vec4 color; } vs_input[3] = {
        {Vec4{0.5, 0.5, 1, 1}, Vec4{1, 0, 0, 1}},
        {Vec4{0.5, -0.5, 1, 1}, Vec4{0, 1, 0, 1}},
        {Vec4{-0.5, -0.5, 1, 1}, Vec4{0, 0, 1, 1}},
    };

    auto orthoMat = CreateOrtho(-1, 1, 1, -1, 1, -1);

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
