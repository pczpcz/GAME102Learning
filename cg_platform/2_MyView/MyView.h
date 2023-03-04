#pragma once

#include "../imgui/imgui.h"

#include <windows.h>

#include <osg/Geometry>
#include <osgViewer/CompositeViewer>
#include <osgGA/TrackballManipulator>
#include <osgViewer/ViewerEventHandlers>
#include <osgGA/TerrainManipulator>
#include <osgGA/StateSetManipulator>

#define _USE_MATH_DEFINES
#undef min
#undef max
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

class MyView
{
public:

	~MyView() 
	{
		m_viewer.setDone(true);
		Sleep(1000);
		m_viewer.stopThreading();
	}

	MyView(int x = 0, int y = 0, int width = 500, int height = 500)  : m_bShowNormal(true)
	{
		//gc
		osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
		traits->x = x;
		traits->y = y;
		traits->width = width;
		traits->height = height;
		traits->windowDecoration = true;
		traits->doubleBuffer = true;
		traits->sharedContext = 0;
		traits->readDISPLAY();
		traits->setUndefinedScreenDetailsToDefaultScreen();

		osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());
		if (!gc.valid()) return;
		gc->setClearColor(osg::Vec4f(0.2f, 0.2f, 0.6f, 1.0f));
		gc->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// main camera
		osg::ref_ptr<osg::Camera> main_camera = m_viewer.getCamera();
		main_camera->setGraphicsContext(gc);
		main_camera->setViewport(new osg::Viewport(traits->x, traits->y, traits->width, traits->height));
		main_camera->setDrawBuffer(GL_BACK);
		main_camera->setReadBuffer(GL_BACK);
		main_camera->setClearMask(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		main_camera->setClearColor(osg::Vec4f(0.2f, 0.2f, 0.4f, 1.0f));
		main_camera->setProjectionMatrixAsPerspective(60.0f, static_cast<double>(traits->width) / static_cast<double>(traits->height), 1.0, 1000.0);
		m_viewer.setCamera(main_camera.get());

		osg::ref_ptr<osgGA::StateSetManipulator> statesetManipulator = new osgGA::StateSetManipulator;
		statesetManipulator->setStateSet(m_viewer.getCamera()->getOrCreateStateSet());
		m_viewer.setCameraManipulator(new osgGA::TerrainManipulator);
		m_viewer.addEventHandler(statesetManipulator.get());
		m_viewer.addEventHandler(new osgViewer::StatsHandler);
	
		//normal camera
		m_normalCamera = new osg::Camera;
		m_normalCamera->setGraphicsContext(gc.get());
		m_normalCamera->setViewport(new osg::Viewport(x, y, width / 4, height / 4));
		GLenum buffer = traits->doubleBuffer ? GL_BACK : GL_FRONT;
		m_normalCamera->setDrawBuffer(buffer);
		m_normalCamera->setReadBuffer(buffer);
		m_normalCamera->setClearMask(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		m_normalCamera->setClearColor(osg::Vec4f(0.0f, 0.0f, 0.0f, 1.0f));
		m_viewer.addSlave(m_normalCamera.get(), osg::Matrixd(), osg::Matrixd::scale(1.0, 1.0, 1.0));
	}

	osg::ref_ptr<osg::Geometry> create_mesh_geometry() 
	{
		osg::ref_ptr<osg::Geometry> geo = new osg::Geometry;

		//生成顶点
		if (!m_vecPointsRef)
			m_vecPointsRef = new osg::Vec3Array;
		m_vecPointsRef->clear();
		m_vecPointsRef->reserve(m_mesh.n_vertices());
		MyMesh::VertexIter v_it, v_end(m_mesh.vertices_end());
		for (v_it = m_mesh.vertices_begin(); v_it != v_end; ++v_it) 
		{
			MyMesh::Point &p = m_mesh.point(*v_it);
			m_vecPointsRef->push_back(osg::Vec3(p[0], p[1], p[2]));

			if (abs(p[0] - 9142.640129) < 10e-3 && abs(p[1] - 1049.542249) < 10e-3) 
			{
				int a = 0;
			}
		}
		geo->setVertexArray(m_vecPointsRef.get());

		//生成三角面片及法线
		osg::DrawElementsUInt* primiset = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
		MyMesh::FaceIter f_it, f_end(m_mesh.faces_end());
		for (f_it = m_mesh.faces_begin(); f_it != f_end; ++f_it) 
		{
			for (auto vh : f_it->vertices()) 
			{
				primiset->push_back(vh.idx());
			}
		}
		geo->addPrimitiveSet(primiset);
		return geo.release();
	}

	void mesh_from_file(const std::string& path) 
	{
		std::ifstream ifs(path);
		if (!ifs.is_open()) return;

		int iVertexNum, iFaceNum;
		ifs >> iVertexNum >> iFaceNum;
		std::vector<MyMesh::VertexHandle> vec_vhs;
		vec_vhs.resize(iVertexNum);
		for (std::size_t i = 0; i < iVertexNum; ++i)
		{
			double x, y, z;
			ifs >> x >> y >> z;
			vec_vhs[i] = m_mesh.add_vertex(MyMesh::Point(x, y, z));
			if (z > 0.0)
				int a = 0;
		}

		std::vector<MyMesh::VertexHandle>  vec_fhs;
		for (std::size_t i = 0; i < iFaceNum; ++i)
		{
			int i1, i2, i3;
			ifs >> i1 >> i2 >> i3;
			vec_fhs.clear();
			vec_fhs.push_back(vec_vhs[i1]);
			vec_fhs.push_back(vec_vhs[i2]);
			vec_fhs.push_back(vec_vhs[i3]);
			m_mesh.add_face(vec_fhs);
		}
	}

	osg::ref_ptr<osg::Geometry> create_test_geomery()
	{
		osg::ref_ptr<osg::Geometry> geo = new osg::Geometry;
		if (!m_vecPointsRef)
			m_vecPointsRef = new osg::Vec3Array;
		geo->setVertexArray(m_vecPointsRef.get());
		m_vecPointsRef->push_back(osg::Vec3(-1.f, 0.f, -1.f));
		m_vecPointsRef->push_back(osg::Vec3(1.f, 0.f, -1.f));
		m_vecPointsRef->push_back(osg::Vec3(1.f, 0.f, 1.f));
		m_vecPointsRef->push_back(osg::Vec3(-1.f, 0.f, 1.f));
		geo->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS, 0, 4));
		return geo.release();
	}

	int run() 
	{
		m_vecNodes.clear();
		if (!m_root)
			m_root = new osg::Group;

		mesh_from_file("E:\\demo.txt");

		osg::ref_ptr<osg::Geode> geode = new osg::Geode;
		osg::ref_ptr<osg::Geometry> geo = create_mesh_geometry();
		m_vecGeometrys.push_back(geo);
		geode->addDrawable(geo.get());
		m_vecNodes.push_back(geode);

		for (int i = 0; i < m_vecNodes.size(); ++i)
			m_root->addChild(m_vecNodes[i]);
		
		m_viewer.setSceneData(m_root.get());
		return m_viewer.run();
	}

private:
	osgViewer::Viewer m_viewer;

	osg::ref_ptr<osg::Camera> m_normalCamera;

	osg::ref_ptr<osg::Group> m_root;
	std::vector<osg::ref_ptr<osg::Node>> m_vecNodes;
	std::vector<osg::ref_ptr<osg::Geometry>> m_vecGeometrys;
	osg::ref_ptr<osgGA::TrackballManipulator> trackball;
	osg::ref_ptr<osg::Vec3Array> m_vecPointsRef;

	MyMesh m_mesh;

	bool m_bShowNormal;
};

