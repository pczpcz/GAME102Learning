#include <iostream>

#include <osgViewer/Viewer>
#include <osgViewer/config/SingleWindow>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/GraphicsContext>
#include <osgDB/Registry>
#include <osgDB/ReadFile>

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

#include "tests/ClockGeometry/ClockGeometry.h"


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
