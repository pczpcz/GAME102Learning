#pragma once
#include "ViewBase.h"

//imgui����osg�Ĵ���, imguiѡ�����-����ѯ������osg�ڵ�״̬-��osg��Ⱦ
class SceneView : public ViewBase
{
public:
	SceneView(const std::string &name) : ViewBase(name) {}
	virtual ~SceneView() {}

	void render_imgui() override {}
	void render_osg() override {}
};

