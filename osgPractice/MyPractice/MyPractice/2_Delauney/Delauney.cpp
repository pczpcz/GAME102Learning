#include "../pch.h"
#include "Delauney.h"

#include <random>
#include <algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <fstream>

void DelauneyTriangulation::addPoints(osg::Vec3 &point)
{
	if (!m_vecPointsRef) 
		m_vecPointsRef = new osg::Vec3Array;
	m_vecPointsRef->push_back(point);
}

void DelauneyTriangulation::addPoints(osg::Vec3Array &vecArrayVertex)
{
	if (!m_vecPointsRef)
		m_vecPointsRef = new osg::Vec3Array;
	m_vecPointsRef->insert(m_vecPointsRef->end(), vecArrayVertex.begin(), vecArrayVertex.end());
}

void DelauneyTriangulation::delauneyTriangulation_OSG()
{
	if (!m_vecPointsRef)
		return;

	m_delaunayTriRef = new osgUtil::DelaunayTriangulator(m_vecPointsRef.get());
	m_delaunayTriRef->triangulate();
	m_PrimitiveRef = m_delaunayTriRef->getTriangles();
}

void DelauneyTriangulation::delauneyTriangulation_CGAL()
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Projection_traits_xy_3<K>  Gt;
	typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
	typedef K::Point_3   Point;
	typedef Delaunay::Finite_face_handles Face_handles;

	/*
	std::ifstream in("2_Delauney/data/triangulation_prog1.cin");
	std::istream_iterator<Point> begin(in);
	std::istream_iterator<Point> end;
	Delaunay dt(begin, end);*/

	std::vector<Point> vecPoints;
	for (auto ittr = m_vecPointsRef->begin(); ittr != m_vecPointsRef->end(); ++ittr) 
	{
		vecPoints.push_back(Point(ittr->x(), ittr->y(), ittr->z()));
	}
	Delaunay dt(vecPoints.begin(), vecPoints.end());

	if (!m_vecPointsRef)
		m_vecPointsRef = new osg::Vec3Array;

	/*
	for (auto ittr = dt.finite_faces_begin(); ittr != dt.finite_faces_end(); ++ittr) 
	{
		for (int i = 0; i < 3; ++i) 
		{
			float x = dt.triangle(ittr).vertex(i).hx();
			float y = dt.triangle(ittr).vertex(i).hy();
			float z = dt.triangle(ittr).vertex(i).hz();
			m_vecPointsRef->push_back(osg::Vec3(x, y, z));
		}
	}*/
	m_vecPointsRef->clear();
	for (auto handle : dt.finite_face_handles())
	{
		for (int i = 0; i < 3; ++i)
		{
			float x = handle->vertex(i)->point().x();
			float y = handle->vertex(i)->point().y();
			float z = handle->vertex(i)->point().z();
			m_vecPointsRef->push_back(osg::Vec3(x, y, z));
		}
	}

	m_PrimitiveRef = new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, m_vecPointsRef->size());

	std::cout << dt.number_of_vertices() << std::endl;
}

void DelauneyTriangulation::delauneyTriangulation_User()
{
	return ;
}

osg::ref_ptr<osg::Vec3Array> DelauneyTriangulation::getRandomExampleVertexSets_1()
{
	osg::ref_ptr<osg::Vec3Array> vecArrayVertex = new osg::Vec3Array;

	std::normal_distribution<float> distribution(0.0, 1.0);
	std::default_random_engine engine(time(0));
	std::vector<float> vecHeight;

	for (int i = 0; i < 10; ++i)
	{
		float z = distribution(engine) * 100.0;
		vecHeight.push_back(z);
	}
	std::sort(vecHeight.begin(), vecHeight.end(), std::greater<float>());
	vecArrayVertex->push_back(osg::Vec3(0.0, 0.0, vecHeight[0]));
	for (int iHeightIndex = 1; iHeightIndex < vecHeight.size(); ++iHeightIndex)
	{
		float r = std::pow(2, iHeightIndex);
		for (float theta = 0; theta < 360; theta += 10)
		{
			float fRand = rand() % 100 / 100.0;
			float x = r * cos(osg::DegreesToRadians(theta));
			float y = r * sin(osg::DegreesToRadians(theta));
			float z = vecHeight[iHeightIndex];
			vecArrayVertex->push_back(osg::Vec3(x, y, z));
		}
	}

	return vecArrayVertex.release();
}

osg::ref_ptr<osg::Vec3Array> DelauneyTriangulation::getRandomExampleVertexSets_2()
{
	osg::ref_ptr<osg::Vec3Array> vecArrayVertex = new osg::Vec3Array;

	float iMax = 500.0;
	srand(time(0));

	for (int i = 0; i < 30; ++i) 
	{
		float x = rand() % 100 / 100.0 * iMax;
		float y = rand() % 100 / 100.0 * iMax;
		float z = 0.0;
		vecArrayVertex->push_back(osg::Vec3(x, y, z) * osg::Matrix::translate(osg::Vec3(500.0, 0.0, 0.0)));
	}

	return vecArrayVertex.release();
}

osg::ref_ptr<osg::Vec3Array> DelauneyTriangulation::getRandomExampleVertexSets_3()
{
	osg::ref_ptr<osg::Vec3Array> vecArrayVertex = new osg::Vec3Array;

	float iMax = 500.0;
	srand(time(0));

	for (int i = 0; i < 30; ++i)
	{
		float x = rand() % 100 / 100.0 * iMax;
		float y = rand() % 100 / 100.0 * iMax;
		float z = 0.0;
		vecArrayVertex->push_back(osg::Vec3(x, y, z) * osg::Matrix::translate(osg::Vec3(-500.0, 0.0, 0.0)));
	}

	/*
	vecArrayVertex->push_back(osg::Vec3(-50.0, -50.0, 0.0));
	vecArrayVertex->push_back(osg::Vec3(-50.0, 50.0, 0.0));
	vecArrayVertex->push_back(osg::Vec3(50.0, 50.0, 0.0));
	vecArrayVertex->push_back(osg::Vec3(50.0, -50.0, 0.0));*/

	return vecArrayVertex.release();
}

osg::ref_ptr<osg::Geometry> DelauneyTriangulation::createGeometry()
{
	osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;

	if (m_vecPointsRef)
		geometry->setVertexArray(m_vecPointsRef.get());

	if (m_PrimitiveRef)
		geometry->addPrimitiveSet(m_PrimitiveRef.get());

	return geometry.release();
}
