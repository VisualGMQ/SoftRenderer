#include "renderer.hpp"
#include "interactive.hpp"
#include <string>

enum UniformVar {
    Texcoord = 0,
};

struct { Vec4 pos; Vec2 texcoord; } const vs_input[4] = {
    {Vec4{0.5, 0.5, 0, 1},   Vec2{1, 1}},
    {Vec4{0.5, -0.5, 0, 1},  Vec2{1, 0}},
    {Vec4{-0.5, -0.5, 0, 1}, Vec2{0, 0}},
    {Vec4{-0.5, 0.5, 0, 1},  Vec2{0, 1}},
};

// simulate Element Indices Buffer
inline unsigned int indices[] = {
    0, 1, 2,
    0, 2, 3,
};

inline unsigned int begin_index = 0;

constexpr int WindowWidth = 480;
constexpr int WindowHeight = 320;

const Mat44 PerspMat = CreatePersp(Radians(90), real(WindowWidth) / WindowHeight, -0.1, -100);
Mat44 ModelMat = Mat44::Eye();
Mat44 ViewMat = CreateTranslate(0, 0, -1);

Surface* WallTexture = nullptr;

class MoveTriangleApp: public App {
public:
    MoveTriangleApp(): App("a - move left; b - move right", WindowWidth, WindowHeight) {}

    void OnInit() override {
        WallTexture = new Surface("./assets/wall.bmp");
        pos_.x = 0;
        pos_.y = 0;

        renderer_.reset(new Renderer(WindowWidth, WindowHeight));

        renderer_->SetClearColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();
        renderer_->SetViewport(0, 0, WindowWidth, WindowHeight);

        std::cout << "ModelMat = " << ModelMat << std::endl;
        std::cout << "ViewMat * ModelMat = " << ViewMat * ModelMat << std::endl;

        renderer_->SetVertexShader([&](int index, ShaderContext& output) {
            output.varyingVec2[Texcoord] = vs_input[indices[begin_index + index]].texcoord;
            return PerspMat * ViewMat * ModelMat * vs_input[indices[begin_index + index]].pos;        
        });

        renderer_->SetFragmentShader([&](ShaderContext& input) {
            return TextureSample(WallTexture, input.varyingVec2[Texcoord]);
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

        ModelMat = CreateTranslate(pos_.x, pos_.y, pos_.z) *
                   CreateRotate(Radians(rotation_.x), Radians(rotation_.y), Radians(rotation_.z));
    }

    void OnRender() override {
        renderer_->SetDrawColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();

        begin_index = 0;
        renderer_->DrawPrimitive();
        begin_index = 3;
        renderer_->DrawPrimitive();

        SwapBuffer(renderer_->GetFramebuffer()->GetRaw());
    }

    void OnQuit() override {
        delete WallTexture;
    }

private:
    std::unique_ptr<Renderer> renderer_;
    Vec3 pos_;
    Vec3 rotation_;
};

int main(int, char**) {
    MoveTriangleApp app;
    app.Run();
    return 0;
}
