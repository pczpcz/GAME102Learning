//#include "..\pch.h"

//#include "Delaunay_Mutithread.h"

//#include <iostream>
//#include <fstream>

///*
//Delaunay_Mutithread::Delaunay_Mutithread()
//{
//	CLog::getInstance();
//}
//
//int Delaunay_Mutithread::getRegularTeatExmapleData(int iVertexNums, float fRatio)
//{
//	//目的：生成100%无重复数据，约束线无相交的数据集
//	int iRowTotal = 5;
//	int iColTotal = 5;
//	float fConstrainRatio = fRatio;
//	if (fConstrainRatio < 0 || fConstrainRatio >= 0.9)
//	{
//		fConstrainRatio = 0.5;
//	}
//
//	if (iVertexNums > 0)
//	{
//		float num = sqrt(iVertexNums);
//		iRowTotal = num + 0.5;
//		iColTotal = num + 0.5;
//	}
//
//	if (iRowTotal < 3) iRowTotal = 3;
//	if (iColTotal < 3) iColTotal = 3;
//
//	int index = 0;
//	float fDilute = 0.01;
//	for (int iRow = 0; iRow < iRowTotal; iRow += 1/*iRowTotal*fDilute*/)
//	{
//		for (int iCol = 0; iCol < iColTotal; iCol += 1/*iColTotal*fDilute*/)
//		{
//			index = iCol + iRow * iColTotal;
//			insert((float)iRow, (float)iCol);
//		}
//	}
//
//	return 0;
//}
//
//void Delaunay_Mutithread::readInputFromFile(const std::string & fileName)
//{
//	std::ifstream f(fileName);
//	if (!f.is_open())
//	{
//		return;
//	}
//	int nVerts;
//	int nEdges;
//	f >> nVerts >> nEdges;
//
//	std::vector<double> vv;
//	vv.reserve(nVerts);
//	for (std::size_t i = 0; i < nVerts; ++i)
//	{
//		double x, y, z;
//		f >> x >> y;
//		insert(x, y);
//	}
//}
//
//void Delaunay_Mutithread::delaunay()
//{
//	if (m_vecMyPoint.empty())
//		return;
//
//	IDelaBella2<double>* idb = IDelaBella2<double>::Create();
//	size_t size = m_vecMyPoint.size();
//
//	int verts = idb->Triangulate(size, &m_vecMyPoint[0].x, &m_vecMyPoint[0].y, sizeof(MyPoint));
//
//	if (verts > 0)
//	{
//		int tris = idb->GetNumPolygons();
//		const IDelaBella2<double>::Simplex* dela = idb->GetFirstDelaunaySimplex();
//		for (int i = 0; i < tris; i++)
//		{
//			if (dela) {
//				m_vecIndexs.push_back(dela->v[0]->i);
//				m_vecIndexs.push_back(dela->v[1]->i);
//				m_vecIndexs.push_back(dela->v[2]->i);
//			}
//			dela = dela->next;
//		}
//	}
//
//	idb->Destroy();
//}
//
//void Delaunay_Mutithread::insert(double x, double y)
//{
//	m_vecMyPoint.push_back(MyPoint(x, y));
//}
//
//osg::ref_ptr<osg::Geometry> Delaunay_Mutithread::createGeometry() 
//{
//	osg::ref_ptr<osg::Geometry> geo = new osg::Geometry;
//	m_vecPointsRef = new osg::Vec3Array;
//	m_PrimiRef = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
//
//	for (auto &data : m_vecMyPoint) 
//	{
//		m_vecPointsRef->push_back(osg::Vec3(data.x, data.y, 0.0));
//	}
//
//	for (auto &index : m_vecIndexs)
//	{
//		m_PrimiRef->push_back(index);
//	}
//
//	geo->setVertexArray(m_vecPointsRef.get());
//	geo->addPrimitiveSet(m_PrimiRef.get());
//
//	return geo.release();
//}
