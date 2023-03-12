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
