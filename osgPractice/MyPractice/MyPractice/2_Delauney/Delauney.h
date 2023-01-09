#pragma once

#include <osg/Geometry>
#include <osgUtil/DelaunayTriangulator>
#include <fstream>

//TODO��
//�������ݼ���û�н��棬���������ȡǰ1000�����ݣ���ÿ���ȸ�ֵһ����ɫ�Աȿ�����
//��ʹ�ù�����Զ����������ǻ�����Ч����
//�����ݽ������뵽�߳���

class DelauneyTriangulation : public osg::Referenced
{
public:
	DelauneyTriangulation();

	osg::ref_ptr<osg::Geometry> createGeometry();

	int getRegularTeatExmapleData(int iVertexNums, float fRatio = 0.5);
	int getExmapleVertexSetsFromFile(const std::string &filename);
	
	void delauneyTriangulation_OSG();
	void delauneyTriangulation_OGRECDT();
	void delauneyTriangulation_CGAL();

	void delaunayConstrainTriangulation_OSG();
	void delauneyConstrainTriangulation_OGRECDT();
	void delauneyConstrainTriangulation_CGAL();

private:
	osg::ref_ptr<osg::Vec3Array> m_vecPointsRef;
	std::vector<std::pair<int, int>> m_vecEdges;			//Լ����
	osg::ref_ptr<osg::Vec4Array> m_vecPointsColorRef;

	osg::ref_ptr<osg::DrawElementsUInt> m_PrimitiveRef;

	osg::ref_ptr<osgUtil::DelaunayTriangulator> m_delaunayTriRef;
	osg::ref_ptr<osgUtil::DelaunayConstraint> m_ContraintsRef;

	std::ofstream m_logFileStream;
};