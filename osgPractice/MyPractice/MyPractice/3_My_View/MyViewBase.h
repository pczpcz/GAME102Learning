#pragma once

/*
#define f("IDC_", )
{
#if

}
*/

#include <Windows.h>

template <typename ParentWnd>
class MyCheckBox1 : public CButton, public MyCtrl<ParentWnd, ParentWnd<T>::type>
{
	DECLARE_DYNAMIC(kkk);
public:
	MyCheckBox(int id):
	{
		m_id = id;
	}

	void Create() override
	{
		if (m_pWnd)
		{
			CButton::Create(m_Name.c_str(), BS_AUTOCHECKBOX | WS_TABSTOP, &m_rect, m_pWnd, m_id);
		}
	}

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV Ö§³Ö

	DECLARE_MESSAGE_MAP()

public:
	bool bChecked;
	afx_msg void OnBnClickedCheck();
};

