#pragma once

/*目的：
 * 1. 统一资源ID的分配
 * 2. 方便类属性和UI的绑定
 */

#include <unordered_set>
#include <Windows.h>

template <typename Wnd, typename T, typename ParentWnd>
class MyWnd
{
	enum CtrlType
	{
		Enum_Type,
		bool_Type,
		value_type,
		string_type
	};

	//MyView的资源ID可以唯一，但是里面包裹CWnd可以相同
public:
	MyWnd(enum CtrlType type, const std::string& name, std::vector<T>& vecNames, std::vector<T>&& vecValues, ParentWnd* pWnd)
	{
		m_pWnd = pWnd;
	}

protected:
	virtual void Create()
	{
		for (int i = 0; i < m_wnds.size(); ++i)
		{
			m_wnds[i]->Update();
		}
	}

private:
	int m_id;
	std::string m_name;

	std::vector<MyCtrl<ParentWnd, T>*> m_wnds;
	std::vector<std::string> m_vecNames;
	std::vector<T> m_vecValues;
	ParentWnd* m_pWnd;
};

