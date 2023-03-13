#pragma once
#include "ViewBase.h"

//imgui包含osg的窗口, imgui选项被操作-》查询并更新osg节点状态-》osg渲染
class SceneView : public ViewBase
{
public:
	SceneView(const std::string &name) : ViewBase(name) {}
	virtual ~SceneView() {}

	//重点：对于imgui内部包含osg窗口的窗口，要把imgui的窗口事件传递给osg窗口，如imgui窗口被关闭了，那通过imgui窗口的回调，osg窗口也应该被关闭（停止view的渲染）

private:

};

