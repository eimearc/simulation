#pragma once

struct WinParams
{
  int spinXFace = 0;
  int spinYFace = 0;
  bool rotate = false;
  bool translate = false;
  int origX = 0;
  int origY = 0;
  int origXPos = 0;
  int origYPos = 0;
  int width = 1024;
  int height = 720;
};

constexpr float INCREMENT = 0.01f;
constexpr float ZOOM = 0.1f;
