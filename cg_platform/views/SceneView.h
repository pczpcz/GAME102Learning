#pragma once
#include "ViewBase.h"

//imgui����osg�Ĵ���, imguiѡ�����-����ѯ������osg�ڵ�״̬-��osg��Ⱦ
class SceneView : public ViewBase
{
public:
	SceneView(const std::string &name) : ViewBase(name) {}
	virtual ~SceneView() {}

	//�ص㣺����imgui�ڲ�����osg���ڵĴ��ڣ�Ҫ��imgui�Ĵ����¼����ݸ�osg���ڣ���imgui���ڱ��ر��ˣ���ͨ��imgui���ڵĻص���osg����ҲӦ�ñ��رգ�ֹͣview����Ⱦ��

private:

};

