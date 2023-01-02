#pragma once

//计划实现/调用以下接口：
//1. osg接口
//2. CGAL接口
//3. 自定义算法接口

//界面功能需求
//1. 点击离散点图标，从文件导入离散点
//2. 执行投影到2D平面
//3. 对2D平面上的点进行Delauney Triangulation
//4. 将2D网格投影回3D曲面，并渲染显示
//5. 参数显示窗口，修改2D点的坐标，点击重新划分按钮可重新划分
//6. 后期可以实现导入CAD的等高线，生成三维的三角网格曲面
#include <osg/Geometry>
#include <osgUtil/DelaunayTriangulator>

class DelauneyTriangulation : public osg::Referenced
{
public:
	void addPoints(osg::Vec3 &point);
	void addPoints(osg::Vec3Array &vecArrayVertex);

	void delauneyTriangulation_OSG();
	void delauneyTriangulation_CGAL();
	void delauneyTriangulation_User();

	static osg::ref_ptr<osg::Vec3Array> getRandomExampleVertexSets_1();
	static osg::ref_ptr<osg::Vec3Array> getRandomExampleVertexSets_2();
	static osg::ref_ptr<osg::Vec3Array> getRandomExampleVertexSets_3();

	osg::ref_ptr<osg::Geometry> createGeometry();

private:
	osg::ref_ptr<osg::Vec3Array> m_vecPointsRef;
	osg::ref_ptr<osg::PrimitiveSet> m_PrimitiveRef;
	osg::ref_ptr<osgUtil::DelaunayTriangulator> m_delaunayTriRef;
};

