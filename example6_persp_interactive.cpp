#include "renderer.hpp"
#include "interactive.hpp"
#include <string>

enum UniformVar {
    Color = 0,
};

struct { Vec4 pos; Vec4 color; } const vs_input[3] = {
    {Vec4{0.5, 0.5, 0, 1}, Vec4{1, 0, 0, 1}},
    {Vec4{0.5, -0.5, 0, 1}, Vec4{0, 1, 0, 1}},
    {Vec4{-0.5, -0.5, 0, 1}, Vec4{0, 0, 1, 1}},
};

constexpr int WindowWidth = 480;
constexpr int WindowHeight = 320;

const Mat44 PerspMat = CreatePersp(Radians(90), real(WindowWidth) / WindowHeight, -0.1, -100);
Mat44 ModelMat = Mat44::Eye();
Mat44 ViewMat = CreateTranslate(0, 0, -1);

class MoveTriangleApp: public App {
public:
    MoveTriangleApp(): App("a - move left; b - move right", WindowWidth, WindowHeight) {}

    void OnInit() override {
        pos_.x = 0;
        pos_.y = 0;

        renderer_ = new Renderer(WindowWidth, WindowHeight);

        renderer_->SetClearColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();
        renderer_->SetViewport(0, 0, WindowWidth, WindowHeight);

        std::cout << "ModelMat = " << ModelMat << std::endl;
        std::cout << "ViewMat * ModelMat = " << ViewMat * ModelMat << std::endl;

        renderer_->SetVertexShader([&](int index, ShaderContext& output) {
            output.varyingVec4[Color] = vs_input[index].color;
            return PerspMat * ViewMat * ModelMat * vs_input[index].pos;        
        });

        renderer_->SetFragmentShader([&](ShaderContext& input) {
            return input.varyingVec4[Color];
        });
    }

    void OnKeyDown(const SDL_KeyboardEvent& e) override {
        if (e.keysym.sym == SDLK_a) {
            pos_.x -= 0.1;
        }
        if (e.keysym.sym == SDLK_d) {
            pos_.x += 0.1;
        }
        if (e.keysym.sym == SDLK_w) {
            pos_.y += 0.1;
        }
        if (e.keysym.sym == SDLK_s) {
            pos_.y -= 0.1;
        }
        if (e.keysym.sym == SDLK_j) {
            pos_.z += 0.1;
        }
        if (e.keysym.sym == SDLK_u) {
            pos_.z -= 0.1;
        }
        if (e.keysym.sym == SDLK_t) {
            rotation_.x -= 5;
        }
        if (e.keysym.sym == SDLK_g) {
            rotation_.x += 5;
        }
        if (e.keysym.sym == SDLK_SPACE) {
            renderer_->SaveDepthBuf("depthbuf.bmp");
        }

        ModelMat = CreateTranslate(pos_.x, pos_.y, pos_.z) *
                   CreateRotate(Radians(rotation_.x), Radians(rotation_.y), Radians(rotation_.z));
    }

    void OnRender() override {
        renderer_->SetDrawColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();
        renderer_->DrawPrimitive();
        SwapBuffer(renderer_->GetFramebuffer()->GetRaw());
    }

    void OnQuit() override {
        delete renderer_;
    }

private:
    Renderer* renderer_;
    Vec3 pos_;
    Vec3 rotation_;
};

int main(int, char**) {
    MoveTriangleApp app;
    app.Run();
    return 0;
}
