#pragma once
//
#include "OsgView.h"

#include "PanelsManager.h"
#include "MenuBar.h"
#include "Console.h"

//主窗口：通过osg的帧循环将事件消息传递给imgui
class MainView : public OsgView
{
public:
    MainView(const std::string& name) : OsgView(name)
    {
    }

    void create_imgui()  override
    {
        /*OvUI::Settings::PanelWindowSettings settings;*/
        //settings.closable = true;
        //settings.collapsable = true;
        //settings.dockable = true;
        //创建菜单栏
        UIManager::instance().m_panelsManager.CreatePanel<OvEditor::Panels::MenuBar>("Menu Bar");
        ////创建控制台窗口
        //m_panelsManager.CreatePanel<OvEditor::Panels::Console>("Console", true, settings);
     }

    ~MainView()
    {
    }

    void size(int& x, int& y, int& width, int& height) override
    {
        osg::GraphicsContext::WindowingSystemInterface* wsi = osg::GraphicsContext::getWindowingSystemInterface();
        if (!wsi) return;

        unsigned int w, h;
        osg::GraphicsContext::ScreenIdentifier main_screen_id;
        main_screen_id.readDISPLAY();
        main_screen_id.setUndefinedScreenDetailsToDefaultScreen();
        wsi->getScreenResolution(main_screen_id, w, h);

        x = 0;
        y = 0;
        width = w;
        height = h;
    }

    void drawUi() override
    {
        UIManager::instance().draw();
    }
};

