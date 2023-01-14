#pragma once

#include <osg/Geometry>
#include <osgUtil/DelaunayTriangulator>
#include <fstream>

//TODO：
//看看数据集有没有交叉，比如随机截取前1000组数据，给每条先赋值一个颜色对比看看；
//先使用规则的自定义数据三角化看看效果。
//将数据解析放入到线程中

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
	void delauneyTriangulation_Triangle();

	void delaunayConstrainTriangulation_OSG();
	void delauneyConstrainTriangulation_OGRECDT();
	void delauneyConstrainTriangulation_CGAL();

private:
	osg::ref_ptr<osg::Vec3Array> m_vecPointsRef;
	std::vector<std::pair<int, int>> m_vecEdges;			//约束边
	osg::ref_ptr<osg::Vec4Array> m_vecPointsColorRef;

	osg::ref_ptr<osg::DrawElementsUInt> m_PrimitiveRef;

	osg::ref_ptr<osgUtil::DelaunayTriangulator> m_delaunayTriRef;
	osg::ref_ptr<osgUtil::DelaunayConstraint> m_ContraintsRef;

	std::ofstream m_logFileStream;
};