#include "renderer.hpp"
#include "interactive.hpp"
#include <string>

enum UniformVar {
    Color = 0,
};

struct { Vec4 pos; Vec4 color; } const vs_input[3] = {
    {Vec4{80, 20, 0, 1}, Vec4{1, 0, 0, 1}},
    {Vec4{400, 20, 0, 1}, Vec4{0, 1, 0, 1}},
    {Vec4{240, 300, 0, 1}, Vec4{0, 0, 1, 1}},
};

constexpr int WindowWidth = 480;
constexpr int WindowHeight = 320;

const Mat44 OrthoMat = CreateOrtho(0, 480, 320, 0, 1, -1);
Mat44 ModelMat = Mat44::Eye();

class MoveTriangleApp: public App {
public:
    MoveTriangleApp(): App("a - move left; b - move right", WindowWidth, WindowHeight) {}

    void OnInit() override {
        pos_.x = 0;
        pos_.y = 0;

        renderer_.reset(new Renderer(WindowWidth, WindowHeight));

        renderer_->SetClearColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();
        renderer_->SetViewport(0, 0, WindowWidth, WindowHeight);

        renderer_->SetVertexShader([&](int index, ShaderContext& output) {
            output.varyingVec4[Color] = vs_input[index].color;
            return OrthoMat * ModelMat * vs_input[index].pos;        
        });

        renderer_->SetFragmentShader([&](ShaderContext& input) {
            return input.varyingVec4[Color];
        });
    }

    void OnKeyDown(const SDL_KeyboardEvent& e) override {
        if (e.keysym.sym == SDLK_a) {
            pos_.x -= 5;
        }
        if (e.keysym.sym == SDLK_d) {
            pos_.x += 5;
        }
        if (e.keysym.sym == SDLK_w) {
            pos_.y -= 5;
        }
        if (e.keysym.sym == SDLK_s) {
            pos_.y += 5;
        }

        ModelMat = CreateTranslate(pos_.x, pos_.y, 0);
    }

    void OnRender() override {
        renderer_->SetDrawColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();
        renderer_->DrawPrimitive();
        SwapBuffer(renderer_->GetFramebuffer()->GetRaw());
    }

private:
    std::unique_ptr<Renderer> renderer_;
    Vec2 pos_;
};

int main() {
    MoveTriangleApp app;
    app.Run();
    return 0;
}
