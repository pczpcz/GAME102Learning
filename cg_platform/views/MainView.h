#pragma once

#include "MainUIManager.h"

//�����ڣ�ͨ��osg��֡ѭ�����¼���Ϣ���ݸ�imgui
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

