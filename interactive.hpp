#pragma once

#include "SDL.h"
#include <string>
#include <chrono>

class App {
public:
    App(const char* title, int w, int h): title_(title) {
        SDL_Init(SDL_INIT_EVERYTHING);
        window_ = SDL_CreateWindow(title,
                                   SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                   w, h,
                                   SDL_WINDOW_SHOWN);
        if (!window_) {
            SDL_Log("can't create window");
        }
        renderer_ = SDL_CreateRenderer(window_, -1, 0);
        if (!renderer_) {
            SDL_Log("can't create SDL renderer");
        }
        SDL_Log("SDL init OK");
    }

    virtual ~App() {
        SDL_DestroyRenderer(renderer_);
        SDL_DestroyWindow(window_);
        SDL_Quit();
    }

    void Run() {
        OnInit();
        auto t = std::chrono::high_resolution_clock::now();
        SDL_Log("start app");
        while (!ShouldExit()) {
            SDL_Event event;
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    Exit();
                }
                if (event.type == SDL_KEYDOWN) {
                    OnKeyDown(event.key);
                }
                if (event.type == SDL_KEYUP) {
                    OnKeyUp(event.key);
                }
                if (event.type == SDL_MOUSEBUTTONDOWN) {
                    OnMouseDown(event.button);
                }
                if (event.type == SDL_MOUSEBUTTONUP) {
                    OnMouseUp(event.button);
                }
                if (event.type == SDL_MOUSEMOTION) {
                    OnMotion(event.motion);
                }
                if (event.type == SDL_WINDOWEVENT) {
                    if (event.window.event == SDL_WINDOWEVENT_RESIZED) {
                        OnWindowResize(event.window.data1, event.window.data2);
                    }
                }
            }
            auto elapse = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t);
            t = std::chrono::high_resolution_clock::now();
            SDL_SetWindowTitle(window_, (title_ + "fps: " + std::to_string(int(1000.0 / elapse.count()))).c_str());
            OnRender();
        }
        OnQuit();
    }

    void Exit() { isQuit_ = true; }
    bool ShouldExit() const { return isQuit_; }

    SDL_Window* GetWindow() const { return window_; }

    void SwapBuffer(SDL_Surface* surface) {
        SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer_, surface);
        if (!texture) {
            SDL_Log("swap buffer failed");
        } else {
            SDL_RenderCopy(renderer_, texture, nullptr, nullptr);
            SDL_DestroyTexture(texture);
        }
        SDL_RenderPresent(renderer_);
    }

    virtual void OnInit() {}
    virtual void OnQuit() {}
    virtual void OnRender() {}

    virtual void OnKeyDown(const SDL_KeyboardEvent&) {}
    virtual void OnKeyUp(const SDL_KeyboardEvent&) {}
    virtual void OnMouseDown(const SDL_MouseButtonEvent&) {}
    virtual void OnMouseUp(const SDL_MouseButtonEvent&) {}
    virtual void OnMotion(const SDL_MouseMotionEvent&) {}
    virtual void OnWindowResize(int, int) {}

private:
    bool isQuit_ = false;
    SDL_Window* window_;
    SDL_Renderer* renderer_;
    std::string title_;
};
