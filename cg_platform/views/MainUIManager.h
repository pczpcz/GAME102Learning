#pragma once

#include <osg/GraphicsContext>
#include <osgViewer/View>
#include <osgViewer/CompositeViewer>
#include <osgGA/TrackballManipulator>

#include "imgui.h"
#include "imgui_impl_opengl3.h"

#include "OvUI/Modules/Canvas.h"
#include "OvUI/Panels/PanelMenuBar.h"
#include "PanelsManager.h"

#include "ViewBase.h"

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

    int run() { return m_mainViewer.run(); }

    void setup_ui()
    {
        //创建imgui控件
        for (int i = 0; i < m_vecViews.size(); ++i) 
            m_vecViews.at(i)->create_imgui();

        //创建scene view -- imgui部分
        m_canvas.MakeDockspace(true);
    }

    void draw()
    {
        ImGui::ShowDemoWindow();
        m_canvas.Draw();
    }

    osgViewer::View* get_view(const std::string& name)
    {
        auto ittr = m_mapViews.find(name);
        if (ittr != m_mapViews.end())
        {
            int id = ittr->second;
            return m_mainViewer.getView(id);
        }
        return nullptr;
    }

    void register_view(ViewBase* imguiview)
    {
        if (!imguiview) return;
        int x, y, width, height;
        imguiview->size(x, y, width, height);

        osgViewer::View* view = get_view(imguiview->get_name());
        if (view) return;

        osg::GraphicsContext::WindowingSystemInterface* wsi = osg::GraphicsContext::getWindowingSystemInterface();
        if (!wsi) return;

        //unsigned int width, height;
        //osg::GraphicsContext::ScreenIdentifier main_screen_id;
        //main_screen_id.readDISPLAY();
        //main_screen_id.setUndefinedScreenDetailsToDefaultScreen();
        //wsi->getScreenResolution(main_screen_id, width, height);

        osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
        traits->x = x;
        traits->y = y;
        traits->width = width;
        traits->height = height;
        traits->windowDecoration;
        traits->doubleBuffer = true;
        traits->sharedContext = 0;
        traits->readDISPLAY();
        traits->setUndefinedScreenDetailsToDefaultScreen();
        if (traits)
        {
            osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());
            osgViewer::View* view = new osgViewer::View;
            view->setName(imguiview->get_name());
            m_mainViewer.addView(view);

            view->getCamera()->setGraphicsContext(gc);
            view->getCamera()->setProjectionMatrixAsPerspective(30.0, double(traits->width) / double(traits->height), 1.0, 1000.0);
            view->getCamera()->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
            view->getCamera()->setGraphicsContext(gc);
            view->setCameraManipulator(new osgGA::TrackballManipulator);

            //嵌合imgui
            imguiview->addEventHandler(view);
            imguiview->setCamera(view->getCamera());
            m_mapViews.insert(std::pair<std::string, int>(imguiview->get_name(), 0));
            m_vecViews.push_back(imguiview);
        }
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

public:
    OvEditor::Core::PanelsManager m_panelsManager;
private:
    //imgui
    OvUI::Modules::Canvas m_canvas;

    //osg
    std::unordered_map<std::string, int> m_mapViews;
    osgViewer::CompositeViewer m_mainViewer;
    std::vector<ViewBase*> m_vecViews;
};

