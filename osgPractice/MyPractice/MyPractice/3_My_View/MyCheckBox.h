#pragma once
#include "MyView.h"
#include "afxwin.h"

class MyCtrl
{
public:
	virtual void Create() {}

	void SetCallback(MyCallback callback)
	{
		m_callback = callback;
	}

public:
	MyCallback m_callback;

public:
	static int m_id;
	CRect m_rect;
	std::string m_Name;
	CWnd* m_pWnd;
};

class MyCallback
{
	virtual void operator()()
	{
	}
};

class MyCheckBox : public CButton, public MyCtrl
{
	DECLARE_DYNAMIC(MyCheckBox);
public:

	MyCheckBox(int id)
	{
		m_id = id;
	}

	void Create() override
	{
		if (m_pWnd)
		{
			CA2T wt(m_Name.c_str());
			CButton::Create(wt, BS_AUTOCHECKBOX|WS_TABSTOP, m_rect, m_pWnd, m_id);
		}
	}

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV Ö§³Ö

	DECLARE_MESSAGE_MAP()

public:
	afx_msg void OnBnClickedCheck();

public:
	bool m_bChecked;
};
