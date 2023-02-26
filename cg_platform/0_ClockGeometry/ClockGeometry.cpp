#include <atlstr.h>
#include <atltime.h>
#include "ClockGeometry.h"

ClockGeometry::ClockGeometry()
	: m_fRadius(50)
	, m_fCircleSegmentLevel(32)
{
	SYSTEMTIME systm;
	GetLocalTime(&systm);
	m_curTime = CTime(systm);

	m_fLenRatio[0] = 0.4;
	m_fLenRatio[1] = 0.5;
	m_fLenRatio[2] = 0.6;

	m_fWidthRatio[0] = 0.02;
	m_fWidthRatio[1] = 0.01;
	m_fWidthRatio[2] = 0.005;

	m_vecPointerIndex.resize(6, 0);

	osg::ref_ptr<TimeUpdataCallback> updateCallback = new TimeUpdataCallback;
	setUpdateCallback(updateCallback.get());
	setUseDisplayList(false);
	setUseVertexBufferObjects(true);
}

void ClockGeometry::createGeometry()
{
	CString strTime = m_curTime.Format(_T("%Y-%m-%d %H:%M:%S %A"));

	int iCircleStartIndex = 0;
	int iCircleEndIndex = 0;
	int iscaleStartIndex = 0;
	int iscaleEndIndex = 0;

	osg::ref_ptr<osg::PrimitiveSet> circlePrimSetRef;
	osg::ref_ptr<osg::PrimitiveSet> scalePrimSetRef;
	std::vector<osg::ref_ptr<osg::PrimitiveSet>> vecTimePointerPrimSetsRef;

	//1. 表盘(顺逆时针，有什么影响，还没有添加颜色，表盘颜色和线条颜色如何控制，和其他部件绘制的先后顺序有什么影响？？)
	float delta = 360.0 / m_fCircleSegmentLevel;
	osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
	for (int i = 0; i < m_fCircleSegmentLevel; ++i)
	{
		float theta = osg::DegreesToRadians(i * delta);
		vertices->push_back(osg::Vec3(m_fRadius * cos(theta), m_fRadius * sin(theta), 0.0));
	}
	iCircleEndIndex = vertices->size();
	circlePrimSetRef = new osg::DrawArrays(osg::DrawArrays::LINE_LOOP, iCircleStartIndex, iCircleEndIndex - iCircleStartIndex);

	//2. 刻度
	iscaleStartIndex = vertices->size();
	delta = 360.0 / 60;
	for (int i = 0; i < 60; ++i)
	{
		float theta = osg::DegreesToRadians(i * delta);
		if (0 == i % 5)
		{
			vertices->push_back(osg::Vec3(0.95 * m_fRadius * cos(theta), 0.95 * m_fRadius * sin(theta), 0.0));
			vertices->push_back(osg::Vec3(0.65 * m_fRadius * cos(theta), 0.65 * m_fRadius * sin(theta), 0.0));
		}
		else 
		{
			vertices->push_back(osg::Vec3(0.95 * m_fRadius * cos(theta), 0.95 * m_fRadius * sin(theta), 0.0));
			vertices->push_back(osg::Vec3(0.75 * m_fRadius * cos(theta), 0.75 * m_fRadius * sin(theta), 0.0));
		}
	}
	iscaleEndIndex = vertices->size();
	scalePrimSetRef = new osg::DrawArrays(osg::DrawArrays::LINES, iscaleStartIndex, iscaleEndIndex - iscaleStartIndex);

	//3. 时针、分针、秒针（运动的，如何优化效率？？？, 要不要独立出去，单独一个drawable）
	int iStart = 0;
	int iEnd = iStart;
	osg::ref_ptr<osg::Vec3Array> verticesTmp;
	for (int iPointer = 0; iPointer < 3; ++iPointer)
	{
		verticesTmp = getRawPointer(iPointer);
		if (verticesTmp->empty()) continue;

		int iStart = vertices->size();
		osg::Matrix matrix = getPointerMatrix(iPointer);
		for (int iVertex = 0; iVertex < verticesTmp->size(); ++iVertex)
		{
			osg::Vec3 point = verticesTmp->at(iVertex) * matrix;
			vertices->push_back(point);
		}
		iEnd = vertices->size();
		m_vecPointerIndex[iPointer * 2] = iStart;
		m_vecPointerIndex[iPointer * 2 + 1] = iEnd;
		vecTimePointerPrimSetsRef.push_back(new osg::DrawArrays(osg::DrawArrays::QUADS, iStart, iEnd - iStart));
	}

	setVertexArray(vertices.get());

	addPrimitiveSet(circlePrimSetRef.get());
	addPrimitiveSet(scalePrimSetRef.get());
	for (auto &primset : vecTimePointerPrimSetsRef)
		addPrimitiveSet(primset.get());
}

void ClockGeometry::updateCurTime()
{
	SYSTEMTIME systm;
	GetLocalTime(&systm);
	m_curTime = CTime(systm);

	osg::ref_ptr<osg::Vec3Array> vertices = dynamic_cast<osg::Vec3Array *>(getVertexArray());
	if (!vertices || ePointerIndex_SecondEnd >= m_vecPointerIndex.size()) return;

	std::vector<int> vecIndex;
	vecIndex.push_back(m_vecPointerIndex[ePointerIndex_HourStart]);
	vecIndex.push_back(m_vecPointerIndex[ePointerIndex_MinuteStart]);
	vecIndex.push_back(m_vecPointerIndex[ePointerIndex_SecondStart]);

	//更新顶点坐标
	osg::ref_ptr<osg::Vec3Array> verticesTmp;
	for (int iPointer = 0; iPointer < 3; ++iPointer)
	{
		verticesTmp = getRawPointer(iPointer);
		if (verticesTmp->empty()) continue;

		osg::Matrix matrix = getPointerMatrix(iPointer);
		for (int iVertex = 0; iVertex < verticesTmp->size(); ++iVertex)
		{
			osg::Vec3 point = verticesTmp->at(iVertex) * matrix;
			int index = vecIndex[iPointer] + iVertex;
			if (index < vertices->size()) 
			{
				vertices->at(index) = point;
			}
		}
	}

	vertices->dirty();
}

osg::ref_ptr<osg::Vec3Array> ClockGeometry::getRawPointer(int index)
{
	if (index >= 3) return new osg::Vec3Array();

	osg::ref_ptr<osg::Vec3Array> verticesTmp = new osg::Vec3Array;
	verticesTmp->push_back(osg::Vec3(0.0, 0.0, 0.0));
	verticesTmp->push_back(osg::Vec3(m_fRadius*m_fWidthRatio[index], 0.0, 0.0));
	verticesTmp->push_back(osg::Vec3(m_fRadius*m_fWidthRatio[index], m_fLenRatio[index] * m_fRadius, 0.0));
	verticesTmp->push_back(osg::Vec3(0.0, m_fLenRatio[index] * m_fRadius, 0.0));

	return verticesTmp.release();
}

osg::Matrix ClockGeometry::getPointerMatrix(int iPointer)
{
	if (iPointer >= 3) return osg::Matrix();
	
	float fRotateDelta[3];
	fRotateDelta[0] = 360.0 / 12;
	fRotateDelta[1] = 360.0 / 60;
	fRotateDelta[2] = 360.0 / 60;

	float fRotate[3];
	fRotate[0] = m_curTime.GetHour() % 12;
	fRotate[1] = m_curTime.GetMinute() % 60;
	fRotate[2] = m_curTime.GetSecond() % 60;

	osg::Matrixf matrix;
	float theta = osg::DegreesToRadians(fRotateDelta[iPointer] * fRotate[iPointer]);
	float xTrans = m_fRadius * m_fWidthRatio[iPointer] * 0.5;
	matrix = osg::Matrixf::translate(osg::Vec3(-1 * xTrans, 0.0, 0.0)) * osg::Matrixf::rotate(osg::Quat(-1 * theta, osg::Vec3(0.0, 0.0, 1.0)));
	return matrix;
}

void TimeUpdataCallback::update(osg::NodeVisitor * nv, osg::Drawable * dw)
{
	osg::ref_ptr<ClockGeometry> clockGeometry = dynamic_cast<ClockGeometry*>(dw);
	if (!clockGeometry) return;

	clockGeometry->updateCurTime();
}
