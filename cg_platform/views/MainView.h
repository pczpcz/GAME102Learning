#pragma once

#include "MainUIManager.h"

//主窗口：通过osg的帧循环将事件消息传递给imgui
class MainView : public OsgView
{
public:
    MainView(const std::string &name)  : OsgView(name)
    {
    }

    ~MainView() 
    {
    }

	void render_imgui() override 
    {
        UIManager::instance().render_imgui();
    }

	void render_osg() override {}

};

