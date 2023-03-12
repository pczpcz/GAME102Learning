#pragma once

#include <osg/GraphicsContext>
#include <osgViewer/View>
#include <osgViewer/CompositeViewer>

#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "OsgImGuiHandler.hpp"

#include "OsgView.h"
#include "SceneView.h"

#include "OvUI/Modules/Canvas.h"
#include "OvUI/Panels/PanelMenuBar.h"

#include "PanelsManager.h"
#include "MenuBar.h"
#include "Console.h"

class Context
{
public:
    static Context& instance()
    {
        static Context context;
        return context;
    }

    ~Context()
    {
        int a = 0;
    }

    osg::ref_ptr<osg::GraphicsContext> get_main_gc()
    {
        if (!m_maingc)
        {
            osg::GraphicsContext::WindowingSystemInterface* wsi = osg::GraphicsContext::getWindowingSystemInterface();
            if (!wsi)
            {
                return nullptr;
            }

            unsigned int width, height;
            osg::GraphicsContext::ScreenIdentifier main_screen_id;

            main_screen_id.readDISPLAY();
            main_screen_id.setUndefinedScreenDetailsToDefaultScreen();
            wsi->getScreenResolution(main_screen_id, width, height);

            osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
            traits->x = 100;
            traits->y = 100;
            traits->width = width - 100;
            traits->height = height - 100;
            traits->windowDecoration = true;
            traits->doubleBuffer = true;
            traits->sharedContext = 0;
            traits->readDISPLAY();
            traits->setUndefinedScreenDetailsToDefaultScreen();

            m_maingc = osg::GraphicsContext::createGraphicsContext(traits.get());
        }
        return m_maingc;
    }

private:
    Context() {}
    Context& operator=(const Context&) = delete;

private:
    osg::ref_ptr<osg::GraphicsContext> m_maingc;
};

class UIManager
{
public:
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

            UIManager::instance().setup_ui();
        }
    };

    static UIManager& instance()
    {
        static UIManager m_uimanager;
        return m_uimanager;
    }

    ~UIManager()
    {
    }

    int create_osgview(osg::GraphicsContext* gc, osg::ref_ptr<OsgView>imguiview)
    {
        if (nullptr == gc || nullptr == imguiview)
            return -1;

        int tid;
        osgViewer::View* view = get_view(imguiview->get_name(), tid);
        if (view) return -1;

        // id++;
        view = new osgViewer::View;   //会不会泄漏？？
        view->setName(imguiview->get_name());
        m_mainViewer.addView(view);
        m_mapViews.insert(std::pair<std::string, int>(imguiview->get_name(), 0));

        view->getCamera()->setGraphicsContext(gc);
        const osg::GraphicsContext::Traits* traits = gc->getTraits();
        if (traits)
        {
            view->getCamera()->setProjectionMatrixAsPerspective(30.0, double(traits->width) / double(traits->height), 1.0, 1000.0);
            view->getCamera()->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
            view->getCamera()->setGraphicsContext(gc);
            view->setCameraManipulator(new osgGA::TrackballManipulator);

            //嵌合imgui
            view->addEventHandler(imguiview);
        }
    }

    osgViewer::View* get_view(const std::string& name, int& id)
    {
        auto ittr = m_mapViews.find(name);
        if (ittr != m_mapViews.end())
        {
            id = ittr->second;
            return m_mainViewer.getView(id);
        }
        return nullptr;
    }

    int run() { return m_mainViewer.run(); }

    void setup_ui()
    {
        OvUI::Settings::PanelWindowSettings settings;
        settings.closable = true;
        settings.collapsable = true;
        settings.dockable = true;

        //创建菜单栏
        m_panelsManager.CreatePanel<OvEditor::Panels::MenuBar>("Menu Bar");

        //创建属性窗口
        m_panelsManager.CreatePanel<OvEditor::Panels::Console>("Console", true, settings);

        m_canvas.MakeDockspace(true);
    }

    void render_imgui()
    {
        ImGui::ShowDemoWindow();
        m_canvas.Draw();
    }

private:
    UIManager() : m_panelsManager(m_canvas)
    {
        init_imgui();
        m_mainViewer.setRealizeOperation(new ImGuiInitOperation);
    }

    UIManager& operator=(const UIManager&) = delete;

    void init_imgui()
    {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();

        // Keyboard mapping. ImGui will use those indices to peek into the io.KeyDown[] array.
        io.KeyMap[ImGuiKey_Tab] = ImGuiKey_Tab;
        io.KeyMap[ImGuiKey_LeftArrow] = ImGuiKey_LeftArrow;
        io.KeyMap[ImGuiKey_RightArrow] = ImGuiKey_RightArrow;
        io.KeyMap[ImGuiKey_UpArrow] = ImGuiKey_UpArrow;
        io.KeyMap[ImGuiKey_DownArrow] = ImGuiKey_DownArrow;
        io.KeyMap[ImGuiKey_PageUp] = ImGuiKey_PageUp;
        io.KeyMap[ImGuiKey_PageDown] = ImGuiKey_PageDown;
        io.KeyMap[ImGuiKey_Home] = ImGuiKey_Home;
        io.KeyMap[ImGuiKey_End] = ImGuiKey_End;
        io.KeyMap[ImGuiKey_Delete] = ImGuiKey_Delete;
        io.KeyMap[ImGuiKey_Backspace] = ImGuiKey_Backspace;
        io.KeyMap[ImGuiKey_Enter] = ImGuiKey_Enter;
        io.KeyMap[ImGuiKey_Escape] = ImGuiKey_Escape;
        io.KeyMap[ImGuiKey_A] = osgGA::GUIEventAdapter::KeySymbol::KEY_A;
        io.KeyMap[ImGuiKey_C] = osgGA::GUIEventAdapter::KeySymbol::KEY_C;
        io.KeyMap[ImGuiKey_V] = osgGA::GUIEventAdapter::KeySymbol::KEY_V;
        io.KeyMap[ImGuiKey_X] = osgGA::GUIEventAdapter::KeySymbol::KEY_X;
        io.KeyMap[ImGuiKey_Y] = osgGA::GUIEventAdapter::KeySymbol::KEY_Y;
        io.KeyMap[ImGuiKey_Z] = osgGA::GUIEventAdapter::KeySymbol::KEY_Z;
    }

private:
    osgViewer::CompositeViewer m_mainViewer;
    std::unordered_map<std::string, int> m_mapViews;
    OvUI::Modules::Canvas m_canvas;
    OvEditor::Core::PanelsManager m_panelsManager;

    //static int id;
};

