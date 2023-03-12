#pragma once

#include <osgViewer/View>
#include <osgViewer/CompositeViewer>
#include "OsgView.h"
#include "imgui_impl_opengl3.h"

#include "OvUI/Modules/Canvas.h"
#include "OvUI/Panels/PanelMenuBar.h"

#include "PanelsManager.h"
#include "MenuBar.h"

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

//主窗口：通过osg的帧循环将事件消息传递给imgui
class MainView : public OsgView
{
public:
    MainView(const std::string &name)  : OsgView(name), m_panelsManager(m_canvas)
    {
        setup_ui();
    }

    ~MainView() 
    {
    }

    void setup_ui() 
    {
        //创建菜单栏
        m_panelsManager.CreatePanel<OvEditor::Panels::MenuBar>("Menu Bar");

        //创建属性窗口

        m_canvas.MakeDockspace(true);
    }

	void render_imgui() override 
    {
        //ImGui::ShowDemoWindow();
        m_canvas.Draw();
    }

	void render_osg() override {}

private:
    OvUI::Modules::Canvas m_canvas;
    OvEditor::Core::PanelsManager m_panelsManager;
};

