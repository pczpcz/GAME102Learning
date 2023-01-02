#pragma once

//�ƻ�ʵ��/�������½ӿڣ�
//1. osg�ӿ�
//2. CGAL�ӿ�
//3. �Զ����㷨�ӿ�

//���湦������
//1. �����ɢ��ͼ�꣬���ļ�������ɢ��
//2. ִ��ͶӰ��2Dƽ��
//3. ��2Dƽ���ϵĵ����Delauney Triangulation
//4. ��2D����ͶӰ��3D���棬����Ⱦ��ʾ
//5. ������ʾ���ڣ��޸�2D������꣬������»��ְ�ť�����»���
//6. ���ڿ���ʵ�ֵ���CAD�ĵȸ��ߣ�������ά��������������
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

