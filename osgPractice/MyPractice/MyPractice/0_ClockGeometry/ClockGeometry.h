#pragma once

#include <osg/Geometry>

class TimeUpdataCallback : public osg::Drawable::UpdateCallback 
{
	virtual void update(osg::NodeVisitor *nv, osg::Drawable* dw);
};

class ClockGeometry : public osg::Geometry
{
public:
	enum EPointerIndex
	{
		ePointerIndex_HourStart,
		ePointerIndex_HourEnd,
		ePointerIndex_MinuteStart,
		ePointerIndex_MinuteEnd,
		ePointerIndex_SecondStart,
		ePointerIndex_SecondEnd,
	};

	ClockGeometry();

	//øΩ±¥ππ‘Ï£ø£ø

	void createGeometry();
	void updateCurTime();

protected:
	osg::ref_ptr<osg::Vec3Array> getRawPointer(int index);
	osg::Matrix getPointerMatrix(int iPointer);

private:
	float m_fRadius;
	float m_fCircleSegmentLevel;

	std::vector<int> m_vecPointerIndex;

	float m_fLenRatio[3];
	float m_fWidthRatio[3];

	CTime m_curTime;
};

