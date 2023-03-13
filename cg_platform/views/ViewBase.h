#pragma once

#include <osgViewer/View>
#include <osgGA/GUIEventAdapter>
#include "imgui.h"

class ViewBase
{
public:
    ViewBase(const std::string& name) : m_strName(name) {}

    ~ViewBase()
    {
    }

    virtual void create_imgui() {}

    virtual void size(int& x, int &y, int &width, int &height) 
    {
        x = 0;
        y = 0;
        width = 1000;
        height = 1000;
    }

    virtual void setCamera(osg::ref_ptr<osg::Camera>) {}

    std::string& get_name() { return m_strName; }

    virtual void addEventHandler(osgViewer::View* view) {}

    virtual void setSceneData(osg::Group* root) {}

protected:
	std::string m_strName;
};

/**
 * Imporant Note: Dear ImGui expects the control Keys indices not to be
 * greater thant 511. It actually uses an array of 512 elements. However,
 * OSG has indices greater than that. So here I do a conversion for special
 * keys between ImGui and OSG.
 */
static int ConvertFromOSGKey(int key)
{
    using KEY = osgGA::GUIEventAdapter::KeySymbol;

    switch (key)
    {
    case KEY::KEY_Tab:
        return ImGuiKey_Tab;
    case KEY::KEY_Left:
        return ImGuiKey_LeftArrow;
    case KEY::KEY_Right:
        return ImGuiKey_RightArrow;
    case KEY::KEY_Up:
        return ImGuiKey_UpArrow;
    case KEY::KEY_Down:
        return ImGuiKey_DownArrow;
    case KEY::KEY_Page_Up:
        return ImGuiKey_PageUp;
    case KEY::KEY_Page_Down:
        return ImGuiKey_PageDown;
    case KEY::KEY_Home:
        return ImGuiKey_Home;
    case KEY::KEY_End:
        return ImGuiKey_End;
    case KEY::KEY_Delete:
        return ImGuiKey_Delete;
    case KEY::KEY_BackSpace:
        return ImGuiKey_Backspace;
    case KEY::KEY_Return:
        return ImGuiKey_Enter;
    case KEY::KEY_Escape:
        return ImGuiKey_Escape;
    default: // Not found
        return -1;
    }
}

