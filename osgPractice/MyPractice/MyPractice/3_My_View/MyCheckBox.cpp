#include "..\pch.h"
#include "..\resource.h"
#include "MyCheckBox.h"

IMPLEMENT_DYNAMIC(MyCheckBox, CButton)

BEGIN_MESSAGE_MAP(MyCheckBox, CButton)
	ON_BN_CLICKED(IDC_CHECK1, &MyCheckBox::OnBnClickedCheck)
END_MESSAGE_MAP()

void MyCheckBox::DoDataExchange(CDataExchange* pDX)
{
	CButton::DoDataExchange(pDX);
	//DDX_Control(pDX, IDC_CHECK1, chek);
}

void MyCheckBox::OnBnClickedCheck()
{




}
