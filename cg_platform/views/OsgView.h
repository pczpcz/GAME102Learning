#pragma once

#include "ViewBase.h"
#include "OsgImGuiHandler.hpp"
#include <osg/camera>

#include "MainUIManager.h"

//osg中包含imgui， osg将事件消息传递给imgui
class OsgView : public ViewBase, public OsgImGuiHandler
{
public:
	OsgView(const std::string& name) : ViewBase(name) 
	{
	}

	~OsgView(){}

	void create_imgui()  override 
	{
		/*OvUI::Settings::PanelWindowSettings settings;*/
		//settings.closable = true;
		//settings.collapsable = true;
		//settings.dockable = true;
		////创建菜单栏
		//m_panelsManager.CreatePanel<OvEditor::Panels::MenuBar>("Menu Bar");
		////创建控制台窗口
		//m_panelsManager.CreatePanel<OvEditor::Panels::Console>("Console", true, settings);
	 }

	void setCamera(osg::ref_ptr<osg::Camera> camera) 
	{
		m_camera = camera;
	}

	void setGraphicContext(osg::ref_ptr<osg::GraphicsContext> gc) 
	{
		if (m_camera)
			m_camera->setGraphicsContext(gc.get());
	}

	void updateOSGWnd(int x, int y, int width, int height)
	{
		osgViewer::GraphicsWindow* gw = dynamic_cast<osgViewer::GraphicsWindow*>(m_camera->getGraphicsContext());
		gw->setWindowRectangle(x, y, width, height);
	}

	void addEventHandler(osgViewer::View* view) 
	{
		if (view)
			view->addEventHandler(this);
	}

	void setSceneData(osg::Group* root) override  
	{
		osgViewer::View *view = UIManager::instance().get_view(m_strName);
		if (view)  view->setSceneData(root);
	}

protected:
	osg::ref_ptr<osg::Camera> m_camera;
};

