#pragma once

#include "ViewBase.h"
#include "OsgImGuiHandler.hpp"
#include <osg/camera>
#include <osgGA/TrackballManipulator>

//osg�а���imgui�� osg���¼���Ϣ���ݸ�imgui
class OsgView : public ViewBase, public OsgImGuiHandler
{
public:
	OsgView(const std::string& name) : ViewBase(name) {}

	~OsgView(){}

	void render_imgui() override //���������imgui�ؼ�
	{
	}

	void render_osg() override {}

	std::string& get_name() { return m_strName; }

protected:
	void drawUi() override { render_imgui(); }

};

