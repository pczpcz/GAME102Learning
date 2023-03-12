#pragma once
#include "ViewBase.h"

//imgui包含osg的窗口, imgui选项被操作-》查询并更新osg节点状态-》osg渲染
class SceneView : public ViewBase
{
public:
	SceneView(const std::string &name) : ViewBase(name) {}
	virtual ~SceneView() {}

	void render_imgui() override {}
	void render_osg() override {}
};

