#include <iostream>

#include <osgViewer/Viewer>
#include <osg/Geometry>
#include <vector>
#include <unordered_map>

//#undef max
//#undef min
//#include "utils/USRefl/include/USRefl_99.h"

#include "tests/ClockGeometry/ClockGeometry.h"

#include "views/MainView.h"
#include "views/SceneView.h"
#include "views/MainUIManager.h"

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
    UIManager::instance();
    //ע�����
    UIManager::instance().register_view(new MainView("Main View"));
    //register_view(new SceneView("Scene View"));
    //
    //


    //2. ������Դ--���������ã�����ƻ���assetģ����е��룬�����Ǽ���Ч�ʣ���ҳLOD������������
    osg::ref_ptr<osg::Group> root = new osg::Group;
    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    osg::ref_ptr<ClockGeometry> geo = new ClockGeometry;
    geo->createGeometry();
    geode->addChild(geo.get());
    root->addChild(geode.get());
    osgViewer::View* view = UIManager::instance().get_view("Main View");
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
