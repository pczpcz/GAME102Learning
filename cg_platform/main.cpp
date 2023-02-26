#include <iostream>

#include <osgViewer/Viewer>
#include <osgViewer/config/SingleWindow>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>

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
#include "1_Delauney/Delaunay_Mutithread.h"

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
            std::cout << "ImGui_ImplOpenGL3_Init() failed\n";
        }
    }
};

class ImGuiDemo : public OsgImGuiHandler
{
protected:
    void drawUi() override
    {
        // ImGui code goes here...
        ImGui::ShowDemoWindow();
    }
};

int main(int argc, char** argv)
{
    osgViewer::Viewer viewer;
    viewer.apply(new osgViewer::SingleWindow(100, 100, 640, 480));
    viewer.setRealizeOperation(new ImGuiInitOperation);
    viewer.addEventHandler(new ImGuiDemo);

    Delaunay_Mutithread dm;
    typedef Delaunay_Mutithread::Point_2 Point_2;
    std::vector<Point_2> vecPoints;
    //dm.readInputFromFile("D:\\inputs\\CapitalA.txt");
    //dm.readInputFromFile("D:\\inputs\\cdt.txt", vecPoints);
    //dmc.readInputFromFile("E:\\inputs\\cornercases.txt");
    //dmc.readInputFromFile("E:\\inputs\\crossing-edges.txt");
    //dm.readInputFromFile("D:\\inputs\\ditch.txt");
    //dm.readInputFromFile("D:\\inputs\\gh_issue.txt");
    //dmc.readInputFromFile("E:\\inputs\\guitar_no_box.txt");
    //dmc.readInputFromFile("E:\\inputs\\Hanging.txt");
    //dmc.readInputFromFile("E:\\inputs\\Hanging2.txt");
    dm.readInputFromFile("D:\\inputs\\island.txt", vecPoints);
    //dm.readInputFromFile("D:\\inputs\\island.txt");
    //dmc.readInputFromFile("E:\\inputs\\issue-42-full-boundary-overlap.txt");
    //dmc.readInputFromFile("E:\\inputs\\issue-42-hole-overlaps-bondary.txt");
    //dmc.readInputFromFile("E:\\inputs\\issue-42-multiple-boundary-overlaps.txt");
    //dmc.readInputFromFile("E:\\inputs\\issue-42-multiple-boundary-overlaps-conform-to-edge.txt");
    //dmc.readInputFromFile("E:\\inputs\\issue-65-wrong-edges.txt");
    //dm.readInputFromFile("D:\\inputs\\kidney.txt");
    //dmc.readInputFromFile("E:\\inputs\\Letter-u.txt");
    //dmc.readInputFromFile("E:\\inputs\\OnEdge.txt");
    //dmc.readInputFromFile("E:\\inputs\\overlapping-constraints.txt");
    //dmc.readInputFromFile("E:\\inputs\\overlapping-constraints2.txt");
    //dmc.readInputFromFile("E:\\inputs\\points_on_constraint_edge.txt");
    //dmc.readInputFromFile("E:\\inputs\\ProblematicCase1.txt");
    //dmc.readInputFromFile("E:\\inputs\\regression_issue_38_wrong_hull_small.txt");
    //dmc.readInputFromFile("E:\\inputs\\square-with-crack.txt");
    //dmc.readInputFromFile("E:\\inputs\\unit-square.txt");

    osg::ref_ptr<osg::Group> root = new osg::Group;
    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    dm.delaunay(vecPoints);
    osg::ref_ptr<osg::Geometry> geometry = dm.createGeometry();
    geode->addDrawable(geometry.get());
    root->addChild(geode.get());
    viewer.setSceneData(root.get());

    viewer.setCameraManipulator(new osgGA::TrackballManipulator);
    viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
    viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));	//ÏÔÊ¾Íø¸ñ
    viewer.addEventHandler(new osgViewer::StatsHandler);

    return viewer.run();
}
