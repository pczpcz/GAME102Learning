#pragma once

#include "ViewBase.h"
#include "OsgImGuiHandler.hpp"
#include <osg/camera>
#include <osgGA/TrackballManipulator>

//osg中包含imgui， osg将事件消息传递给imgui
class OsgView : public ViewBase, public OsgImGuiHandler
{
public:
	OsgView(const std::string& name) : ViewBase(name) {}

	~OsgView(){}

	void render_imgui() override //在这里绘制imgui控件
	{
	}

	void render_osg() override {}

	std::string& get_name() { return m_strName; }

protected:
	void drawUi() override { render_imgui(); }

};

