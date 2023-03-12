#pragma once

#include <osgViewer/View>
#include <osgViewer/CompositeViewer>
#include "OsgView.h"
#include "imgui_impl_opengl3.h"

class ImGuiInitOperation : public osg::Operation
{
public:
    ImGuiInitOperation()
        : osg::Operation("ImGuiInitOperation", false)
    {
    }

    void operator()(osg::Object* object) override
    {
        osg::GraphicsContext* context = dynamic_cast<osg::GraphicsContext*>(object);
        if (!context)
            return;

        if (!ImGui_ImplOpenGL3_Init())
        {
            //
        }
    }
};

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
        ImGui::ShowDemoWindow();
    }

	void render_osg() override {}

private:

};

