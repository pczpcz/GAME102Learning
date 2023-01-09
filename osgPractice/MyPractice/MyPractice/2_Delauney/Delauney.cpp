#include "CDT.h"
#include "remove_at.hpp"

#include "../pch.h"

#include "Delauney.h"

#include <random>
#include <algorithm>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
DelauneyTriangulation::DelauneyTriangulation()
{
	m_vecPointsRef = new osg::Vec3Array;
	m_ContraintsRef = new osgUtil::DelaunayConstraint;
	m_delaunayTriRef = new osgUtil::DelaunayTriangulator;

	m_logFileStream.open("D:\\log.txt");
}

void DelauneyTriangulation::delauneyTriangulation_OSG()
{
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

	/*遍历方法一：
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

	//遍历方法二：
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

	//m_PrimitiveRef = new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, m_vecPointsRef->size());

	std::cout << dt.number_of_vertices() << std::endl;
}

//还有问题：100w数据点，有一条很长的划线
void DelauneyTriangulation::delauneyTriangulation_OGRECDT()
{
	CDT::Triangulation<double> cdt;
	try 
	{
		cdt.insertVertices(
			m_vecPointsRef.get()->begin(),
			m_vecPointsRef.get()->end(),
			[](const osg::Vec3& p) { return p.x(); },
			[](const osg::Vec3& p) { return p.y(); }
		);
		cdt.eraseSuperTriangle();
	}
	catch (std::runtime_error e)
	{
		std::cout << e.what() << std::endl;
	}

	osg::ref_ptr<osg::DrawElementsUInt> primitiveRef = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
	for (auto ittr = cdt.triangles.begin(); ittr != cdt.triangles.end(); ++ittr)
	{
		unsigned int p1 = ittr->vertices[0] /* - unsigned int(3)*/;
		unsigned int p2 = ittr->vertices[1] /* - unsigned int(3)*/;
		unsigned int p3 = ittr->vertices[2] /* - unsigned int(3)*/;

		if (p1 < m_vecPointsRef->size() && p2 < m_vecPointsRef->size() && p3 < m_vecPointsRef->size()) 
		{
			float md1 = abs(m_vecPointsRef->at(p1).x() - m_vecPointsRef->at(p2).x()) + abs(m_vecPointsRef->at(p1).y() - m_vecPointsRef->at(p2).y());
			float md2 = abs(m_vecPointsRef->at(p1).x() - m_vecPointsRef->at(p3).x()) + abs(m_vecPointsRef->at(p1).y() - m_vecPointsRef->at(p3).y());
			float md3 = abs(m_vecPointsRef->at(p2).x() - m_vecPointsRef->at(p3).x()) + abs(m_vecPointsRef->at(p2).y() - m_vecPointsRef->at(p3).y());
			if (md1 > 100.0 || md2 > 100.0 || md3 > 100.0)
			{
				//continue;
			}

			primitiveRef->push_back(p1);
			primitiveRef->push_back(p2);
			primitiveRef->push_back(p3);
		}
		else 
		{
			m_logFileStream << p1 << ", " << p2 << "," << p3 << std::endl;
		}
	}

	m_PrimitiveRef = primitiveRef;
}

void DelauneyTriangulation::delauneyConstrainTriangulation_OGRECDT()
{
	CDT::Triangulation<double> cdt;
	cdt.insertVertices(
		m_vecPointsRef.get()->begin(),
		m_vecPointsRef.get()->end(),
		[](const osg::Vec3& p) { return p.x(); },
		[](const osg::Vec3& p) { return p.y(); }
	);
	cdt.insertEdges(
		m_vecEdges.begin(),
		m_vecEdges.end(),
		[](std::pair<int, int>& p) { return p.first; },
		[](std::pair<int, int>& p) { return p.second; }
	);
	cdt.eraseOuterTrianglesAndHoles();

	osg::ref_ptr<osg::DrawElementsUInt> primitiveRef = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, cdt.triangles.size() * 3);
	for (int i = 0; i < cdt.triangles.size(); ++i)
	{
		primitiveRef->addElement(cdt.triangles[i].vertices[0]);
		primitiveRef->addElement(cdt.triangles[i].vertices[1]);
		primitiveRef->addElement(cdt.triangles[i].vertices[2]);
	}

	m_PrimitiveRef = primitiveRef;
}

void DelauneyTriangulation::delauneyConstrainTriangulation_CGAL()
{
}

void DelauneyTriangulation::delaunayConstrainTriangulation_OSG()
{
	m_ContraintsRef = new osgUtil::DelaunayConstraint;
	osg::ref_ptr<osg::DrawElementsUInt> primitiveRef= new osg::DrawElementsUInt(osg::PrimitiveSet::LINES);
	for (int i = 0; i < m_vecEdges.size(); ++i) 
	{
		primitiveRef->addElement(m_vecEdges[i].first);
		primitiveRef->addElement(m_vecEdges[i].second);
	}
	m_ContraintsRef->addPrimitiveSet(primitiveRef.get());

	m_delaunayTriRef->addInputConstraint(m_ContraintsRef.get());
	m_delaunayTriRef->setInputPointArray(m_vecPointsRef.get());

	int iTimeStart = ::GetTickCount();
	int iTimeEnd = iTimeStart;

	m_delaunayTriRef->triangulate();

	iTimeEnd = ::GetTickCount();
	std::cout << "delaunay triangulation time used: " << (iTimeEnd - iTimeStart) / 1000.0 / 60.0 << "(min)" << std::endl;

	m_delaunayTriRef->removeInternalTriangles(m_ContraintsRef.get());
	iTimeEnd = ::GetTickCount();
	std::cout << "delaunay removeInternalTriangles time used: " << (iTimeEnd - iTimeStart) / 1000.0 / 60.0 << "(min)" << std::endl;

	m_PrimitiveRef = m_delaunayTriRef->getTriangles();

}

/* 3D测试数据
osg::ref_ptr<osg::Vec3Array> DelauneyTriangulation::getRandomExampleVertexSets_1(osg::Matrix &transform)
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

			//当把y, z对调后，由于对于（x, z）, 可能有2个y值与之对应，这种情况下的三维空间点的三角化出现错乱，
			//由此推测：在osg中的该接口都是先投影到x，y平面，在进行的三角化，之后再加上高程信息，恢复到三维空间中
			//所谓的2D和3D三角化的区别：一个是对表面的划分（通过三角形），一个是对三维空间的划分（通过四面体）
			//三维空间的表面应该没有3D dylaunay triangulation的提法，除非三维实体的点云进行3D dylaunay triangulation，找到它的凸包表面？？
			//其他的扩展：（猜的）
			//(1) CGAL的Tutorials的Gis那个示例：DSM到DTM过程中，要移除植被和建筑物，如果有这些点的存在，是不是会导致网格划分失败，这个示例后面还要好好看看，包括补充一些DSM/DTM gis相关的知识
			//(2) 可能应用一：利用这个限制，可以把限制线围成的那个面剔除，比如把建筑物识别出来，再剔除，对应位置就形成了一个大洞
			//(2) 可能应用二：利用这个限制，可以限制网格的形成方式，比如先把多边形的中轴线找到，把中轴线作为限制，进行网格划分

			//float z = r * sin(osg::DegreesToRadians(theta));
			//float y = vecHeight[iHeightIndex];

			vecArrayVertex->push_back(osg::Vec3(x, y, z)* transform);
		}
	}

	return vecArrayVertex.release();
}*/

int DelauneyTriangulation::getExmapleVertexSetsFromFile(const std::string & filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open()) return -1;

	std::string str;
	std::vector<std::string> vecStr;
	bool bHasPritiveSetToAdd = false;
	bool bNewLine = false;
	int iConstrainPointCount = 0;

	int iEdgeStartIndex = 0;
	int iEdgeEndIndex = 0;

	int iTimeStart = ::GetTickCount();
	int iTimeEnd = iTimeStart;
	while(std::getline(ifs, str))
	{
		boost::split(vecStr, str, boost::is_any_of(","));
		if (!vecStr.empty())
		{
			if (vecStr.size() == 1)
			{
				if (vecStr[0] == "") continue;
				bNewLine = true;	//linex or lint count
			}

			if (bNewLine)
			{
				if (bHasPritiveSetToAdd && iConstrainPointCount==2)
				{
					m_vecEdges.push_back(std::make_pair(iEdgeStartIndex, iEdgeEndIndex));
					iEdgeStartIndex = iEdgeEndIndex;
				}
				iConstrainPointCount = 0;
				iEdgeStartIndex = iEdgeEndIndex;
				bHasPritiveSetToAdd = false;
			}

			if (vecStr.size() == 3)
			{
				++iConstrainPointCount;
				bNewLine = false;
				bHasPritiveSetToAdd = true;
				m_vecPointsRef->push_back(std::move(osg::Vec3(std::stof(vecStr[0]), std::stof(vecStr[1]), std::stof(vecStr[2]))));
			}
		}
	}

	iTimeEnd = ::GetTickCount();
	std::cout << "data reading and pretreament time used: " << (iTimeEnd - iTimeStart) / 1000.0 / 60.0 << "(min)" << std::endl;

	//处理最后一行的数据
	if (bHasPritiveSetToAdd && iConstrainPointCount == 2)
	{
		m_vecEdges.push_back(std::make_pair(iEdgeStartIndex, iEdgeEndIndex));
	}

	//去重
	CDT::DuplicatesInfo di;
	di = CDT::FindDuplicates<float>(
		m_vecPointsRef.get()->begin(),
		m_vecPointsRef.get()->end(),
		[](const osg::Vec3& p) { return p.x(); },
		[](const osg::Vec3& p) { return p.y(); }
	);
	osg::Vec3Array &vecVertex = *m_vecPointsRef.get();
	vecVertex.erase(
		remove_at(vecVertex.begin(), vecVertex.end(), di.duplicates.begin(), di.duplicates.end()),
		vecVertex.end()
	);
}

int DelauneyTriangulation::getRegularTeatExmapleData(int iVertexNums, float fRatio)
{
	//目的：生成100%无重复数据，约束线无相交的数据集
	int iRowTotal = 10;
	int iColTotal = 10;
	float fConstrainRatio = fRatio;
	if (fConstrainRatio < 0 || fConstrainRatio >= 0.9) 
	{
		fConstrainRatio = 0.5;
	}

	if (iVertexNums > 0)
	{
		float num = sqrt(iVertexNums);
		iRowTotal = num + 0.5;
		iColTotal = num + 0.5;
	}

	if (iRowTotal < 10) iRowTotal = 10;
	if (iColTotal < 10) iColTotal = 10;

	int index = 0;
	float fDilute = 0.01;
	for (int iRow = 0; iRow < iRowTotal; iRow+=1/*iRowTotal*fDilute*/)
	{
		for (int iCol = 0; iCol < iColTotal; iCol+=1/*iColTotal*fDilute*/)
		{
			index = iCol + iRow * iColTotal;
			m_vecPointsRef->push_back(osg::Vec3((float)iRow, (float)iCol, 0.0));
		}
	}
	
	int iConstrainEdgesNum = iRowTotal * iColTotal * fConstrainRatio;
	int iCount = 0;
	int iEdgeStartIndex = 0;
	int iEdgeEndIndex = 0;
	for (int iRow = 0; iRow < iRowTotal-1; ++iRow)
	{
		for (int iCol = 0; iCol < iColTotal - 1; iCol+=2)
		{
			iEdgeStartIndex = iRow + iCol * iRowTotal;
			iEdgeEndIndex = iRow + (iCol + 1) * iRowTotal;
			
			++iCount;
			m_vecEdges.push_back(std::pair<int, int>(iEdgeStartIndex, iEdgeEndIndex));

			if (iCount >= iConstrainEdgesNum) 
				return 0;
		}
	}

	if (iCount < iConstrainEdgesNum) 
	{
		iEdgeStartIndex = 0;
		iEdgeEndIndex = 0;
		for (int iCol = 0; iCol < iColTotal - 1; iCol += 2)
		{
			for (int iRow = 0; iRow < iRowTotal -1; ++iRow)
			{
				iEdgeStartIndex = iCol + iRow * iColTotal;
				iEdgeEndIndex = iCol + (iRow + 1) * iColTotal;
				
				++iCount;
				m_vecEdges.push_back(std::pair<int, int>(iEdgeStartIndex, iEdgeEndIndex));

				if (iCount >= iConstrainEdgesNum)
					return 0;
			}
		}
	}

	return 0;
}

osg::ref_ptr<osg::Geometry> DelauneyTriangulation::createGeometry()
{
	osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;
	geometry->setVertexArray(m_vecPointsRef.get());
	geometry->addPrimitiveSet(m_PrimitiveRef.get());
	return geometry.release();
}
