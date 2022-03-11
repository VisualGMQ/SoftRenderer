#include "renderer.hpp"
#include "interactive.hpp"
#include <string>

enum UniformVar {
    Texcoord = 0,
};

const struct {Vec4 pos; Vec2 texcoord;} vs_input[36] = {
    {Vec4{-0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 0.0f}}, // Bottom-left
    {Vec4{ 0.5f,  0.5f, -0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-right
    {Vec4{ 0.5f, -0.5f, -0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // bottom-right         
    {Vec4{ 0.5f,  0.5f, -0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-right
    {Vec4{-0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 0.0f}}, // bottom-left
    {Vec4{-0.5f,  0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // top-left

    {Vec4{-0.5f, -0.5f,  0.5f, 1.0f},  Vec2{0.0f, 0.0f}}, // bottom-left
    {Vec4{ 0.5f, -0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // bottom-right
    {Vec4{ 0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-right
    {Vec4{ 0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-right
    {Vec4{-0.5f,  0.5f,  0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // top-left
    {Vec4{-0.5f, -0.5f,  0.5f, 1.0f},  Vec2{0.0f, 0.0f}}, // bottom-left

    {Vec4{-0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // top-right
    {Vec4{-0.5f,  0.5f, -0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-left
    {Vec4{-0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // bottom-left
    {Vec4{-0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // bottom-left
    {Vec4{-0.5f, -0.5f,  0.5f, 1.0f},  Vec2{0.0f, 0.0f}}, // bottom-right
    {Vec4{-0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // top-right

    {Vec4{ 0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // top-left
    {Vec4{ 0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // bottom-right
    {Vec4{ 0.5f,  0.5f, -0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-right         
    {Vec4{ 0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // bottom-right
    {Vec4{ 0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // top-left
    {Vec4{ 0.5f, -0.5f,  0.5f, 1.0f},  Vec2{0.0f, 0.0f}}, // bottom-left     

    {Vec4{-0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // top-right
    {Vec4{ 0.5f, -0.5f, -0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-left
    {Vec4{ 0.5f, -0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // bottom-left
    {Vec4{ 0.5f, -0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // bottom-left
    {Vec4{-0.5f, -0.5f,  0.5f, 1.0f},  Vec2{0.0f, 0.0f}}, // bottom-right
    {Vec4{-0.5f, -0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // top-right

    {Vec4{-0.5f,  0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // top-left
    {Vec4{ 0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // bottom-right
    {Vec4{ 0.5f,  0.5f, -0.5f, 1.0f},  Vec2{1.0f, 1.0f}}, // top-right     
    {Vec4{ 0.5f,  0.5f,  0.5f, 1.0f},  Vec2{1.0f, 0.0f}}, // bottom-right
    {Vec4{-0.5f,  0.5f, -0.5f, 1.0f},  Vec2{0.0f, 1.0f}}, // top-left
    {Vec4{-0.5f,  0.5f,  0.5f, 1.0f},  Vec2{0.0f, 0.0f}}  // bottom-left        
};

constexpr int WindowWidth = 1024;
constexpr int WindowHeight = 720;

const Mat44 PerspMat = CreatePersp(Radians(90), real(WindowWidth) / WindowHeight, -0.1, -100);
Mat44 ModelMat = Mat44::Eye();
Mat44 ViewMat = CreateTranslate(0, 0, -1.5);

Surface* BoxTexture = nullptr;
inline unsigned int BeginIndex = 0;

class MoveTriangleApp: public App {
public:
    MoveTriangleApp(): App("a - move left; b - move right", WindowWidth, WindowHeight) {}

    void OnInit() override {
        BoxTexture = new Surface("./assets/container.jpg");
        pos_.x = 0;
        pos_.y = 0;

        renderer_.reset(new Renderer(WindowWidth, WindowHeight));
        renderer_->SetFaceCull(CW);

        renderer_->SetClearColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();
        renderer_->SetViewport(0, 0, WindowWidth, WindowHeight);

        std::cout << "ModelMat = " << ModelMat << std::endl;
        std::cout << "ViewMat * ModelMat = " << ViewMat * ModelMat << std::endl;

        renderer_->SetVertexShader([&](int index, ShaderContext& output) {
            output.varyingVec2[Texcoord] = vs_input[BeginIndex + index].texcoord;
            return PerspMat * ViewMat * ModelMat * vs_input[BeginIndex + index].pos;        
        });

        renderer_->SetFragmentShader([&](ShaderContext& input) {
            return TextureSample(BoxTexture, input.varyingVec2[Texcoord]);
        });
    }

    void OnKeyDown(const SDL_KeyboardEvent& e) override {
        if (e.keysym.sym == SDLK_a) {
            rotation_.y -= 5;
        }
        if (e.keysym.sym == SDLK_d) {
            rotation_.y += 5;
        }
        if (e.keysym.sym == SDLK_w) {
            rotation_.x += 5;
        }
        if (e.keysym.sym == SDLK_s) {
            rotation_.x -= 5;
        }
        if (e.keysym.sym == SDLK_u) {
            rotation_.z += 5;
        }
        if (e.keysym.sym == SDLK_j) {
            rotation_.z -= 5;
        }

        ModelMat = CreateRotate(Radians(rotation_.x), Radians(rotation_.y), Radians(rotation_.z));
    }

    void OnRender() override {
        renderer_->SetDrawColor(Color4{0.1, 0.1, 0.1, 1});
        renderer_->Clear();

        BeginIndex = 0;
        for (int i = 0; i < 12; i++) {
            renderer_->DrawPrimitive();
            BeginIndex += 3;
        }

        SwapBuffer(renderer_->GetFramebuffer()->GetRaw());
    }

    void OnQuit() override {
        delete BoxTexture;
    }

private:
    std::unique_ptr<Renderer> renderer_;
    Vec3 pos_;
    Vec3 rotation_;
};

int main(int, char**) {
    Renderer::Init();
    MoveTriangleApp app;
    app.Run();
    Renderer::Quit();
    return 0;
}
