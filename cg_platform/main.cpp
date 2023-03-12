#include <iostream>

#include <osgViewer/Viewer>
#include <osgViewer/config/SingleWindow>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/GraphicsContext>
#include <osgDB/Registry>
#include <osgDB/ReadFile>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_opengl3.h"
#include "imgui/OsgImGuiHandler.hpp"

#include <osg/Geometry>
#include <vector>
#include <unordered_map>
#include <deque>
#include <iterator>
#include <algorithm>

#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#undef max
#undef min
#include "algorithm/Triangulation/Delauney/Delaunay_Mutithread.h"
//#include "utils/USRefl/include/USRefl_99.h"

#include "views/MainView.h"
#include "views/SceneView.h"

#include "tests/ClockGeometry/ClockGeometry.h"

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
            traits->width = width-100;
            traits->height = height-100;
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

    osgViewer::View* get_view(const std::string& name, int &id)
    {
        auto ittr = m_mapViews.find(name);
        if (ittr != m_mapViews.end())
        {
           id = ittr->second;
           return m_mainViewer.getView(id);
        }
        return nullptr;
    }

    int run()  { return m_mainViewer.run();  }

private:
    UIManager()
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
    
    //static int id;
};

int main(int argc, char** argv)
{
    //----------------------------------------todo list--------------------------------------------------------//
    //1. asset导入


    //----------------------------------------关于设计--------------------------------------------------------//
    //1. 主窗口使用compositeViewer，这样内部的osg子窗口View可以共享同一帧循环，且内部的子窗口可以有不同的上下文，也就可能/可以各自独立的响应鼠标事件
    //2. 部分模块可能要设计成单例，这样界面之间通信会方便一点


    //-----------------------------------------主程序---------------------------------------------------------//
    //窗口初始化逻辑
    //考虑：内部osg窗口，要不要和主osg窗口共享同一个帧循环
    //1. 主窗口：（osg包含imgui的界面）--MainView
    osg::ref_ptr<MainView> mainview = new MainView("MainView");
    UIManager::instance().create_osgview(Context::instance().get_main_gc(), mainview);

    //2. 加载资源--仅仅测试用，后面计划由asset模块进行导入，并考虑加载效率，分页LOD。。。？？？
    osg::ref_ptr<osg::Group> root = new osg::Group;
    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    osg::ref_ptr<ClockGeometry> geo = new ClockGeometry;
    geo->createGeometry();
    geode->addChild(geo.get());
    root->addChild(geode.get());
    int id;
    osgViewer::View* view = UIManager::instance().get_view(mainview->get_name(), id);
    if (!view) return 0;
    view->setSceneData(root.get());

    //3. 主窗口内部有两种子窗口
    //（1）单纯imgui窗口（也包括imgui窗口包含imgui子窗口的情况）
    //（2）包含有osg窗口的imgui窗口 -- 如， SceneView
    //SceneView sceneview;                //todo: scene使用自己独立的gc，而如果窗口内还要显示imgui控件，则这部分的imgui控件使用这个独立的gc？？？
    //sceneview.register_view(&mainview);

    //4. 开启渲染循环
    return UIManager::instance().run();
}
