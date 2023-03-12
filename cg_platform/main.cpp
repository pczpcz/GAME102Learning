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
        view = new osgViewer::View;   //�᲻��й©����
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

            //Ƕ��imgui
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
    //1. asset����


    //----------------------------------------�������--------------------------------------------------------//
    //1. ������ʹ��compositeViewer�������ڲ���osg�Ӵ���View���Թ���ͬһ֡ѭ�������ڲ����Ӵ��ڿ����в�ͬ�������ģ�Ҳ�Ϳ���/���Ը��Զ�������Ӧ����¼�
    //2. ����ģ�����Ҫ��Ƴɵ�������������֮��ͨ�Ż᷽��һ��


    //-----------------------------------------������---------------------------------------------------------//
    //���ڳ�ʼ���߼�
    //���ǣ��ڲ�osg���ڣ�Ҫ��Ҫ����osg���ڹ���ͬһ��֡ѭ��
    //1. �����ڣ���osg����imgui�Ľ��棩--MainView
    osg::ref_ptr<MainView> mainview = new MainView("MainView");
    UIManager::instance().create_osgview(Context::instance().get_main_gc(), mainview);

    //2. ������Դ--���������ã�����ƻ���assetģ����е��룬�����Ǽ���Ч�ʣ���ҳLOD������������
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

    //3. �������ڲ��������Ӵ���
    //��1������imgui���ڣ�Ҳ����imgui���ڰ���imgui�Ӵ��ڵ������
    //��2��������osg���ڵ�imgui���� -- �磬 SceneView
    //SceneView sceneview;                //todo: sceneʹ���Լ�������gc������������ڻ�Ҫ��ʾimgui�ؼ������ⲿ�ֵ�imgui�ؼ�ʹ�����������gc������
    //sceneview.register_view(&mainview);

    //4. ������Ⱦѭ��
    return UIManager::instance().run();
}
