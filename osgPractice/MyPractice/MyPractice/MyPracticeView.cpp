
// MyPracticeView.cpp: CMyPracticeView 类的实现
//

#include "pch.h"
#include "framework.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "MyPractice.h"
#endif

#include "MyPracticeDoc.h"
#include "MyPracticeView.h"

#include <osgViewer/ViewerEventHandlers>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osgDB/ReadFile>
#include <osgViewer/Viewer>
#include <osgViewer/api/Win32/GraphicsWindowWin32>
#include <osg/Geode>
#include <osg/Group>

#include "0_ClockGeometry/ClockGeometry.h"
#include "2_Delauney/Delauney.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

osg::ref_ptr<osgViewer::Viewer> g_viewer;
bool g_finished;

void render(void *)
{
	while (!g_viewer->done())
	{
		g_viewer->frame();
	}
	g_finished = true;
}

// CMyPracticeView

IMPLEMENT_DYNCREATE(CMyPracticeView, CView)

BEGIN_MESSAGE_MAP(CMyPracticeView, CView)
	// 标准打印命令
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CMyPracticeView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_CREATE()
END_MESSAGE_MAP()

// CMyPracticeView 构造/析构

CMyPracticeView::CMyPracticeView() noexcept
{
	// TODO: 在此处添加构造代码

}

CMyPracticeView::~CMyPracticeView()
{
	g_viewer->setDone(true);
	Sleep(1000);
	g_viewer->stopThreading();
}

BOOL CMyPracticeView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: 在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式

	return CView::PreCreateWindow(cs);
}

// CMyPracticeView 绘图

void CMyPracticeView::OnDraw(CDC* /*pDC*/)
{
	CMyPracticeDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: 在此处为本机数据添加绘制代码
}


// CMyPracticeView 打印


void CMyPracticeView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CMyPracticeView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 默认准备
	return DoPreparePrinting(pInfo);
}

void CMyPracticeView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 添加额外的打印前进行的初始化过程
}

void CMyPracticeView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 添加打印后进行的清理过程
}

void CMyPracticeView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CMyPracticeView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CMyPracticeView 诊断

#ifdef _DEBUG
void CMyPracticeView::AssertValid() const
{
	CView::AssertValid();
}

void CMyPracticeView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CMyPracticeDoc* CMyPracticeView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMyPracticeDoc)));
	return (CMyPracticeDoc*)m_pDocument;
}
#endif //_DEBUG


// CMyPracticeView 消息处理程序


int CMyPracticeView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	return 0;
}


void CMyPracticeView::OnInitialUpdate()
{
	CView::OnInitialUpdate();

	osg::ref_ptr<osg::Referenced> windata = new osgViewer::GraphicsWindowWin32::WindowData(m_hWnd);

	CRect rect;
	GetWindowRect(&rect);

	osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
	traits->x = 0;
	traits->y = 0;
	traits->width = rect.Width();
	traits->height = rect.Height();
	traits->windowDecoration = false;
	traits->doubleBuffer = true;
	traits->inheritedWindowData = windata;

	osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());

	g_viewer = new osgViewer::Viewer;
	osg::ref_ptr<osg::Camera> camera = g_viewer->getCamera();
	camera->setGraphicsContext(gc);
	camera->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
	camera->setProjectionMatrixAsPerspective(60.0f, static_cast<double>(traits->width)/ static_cast<double>(traits->height), 1.0f, 1000.0f);

	g_viewer->setCamera(camera.get());

	osg::ref_ptr<osg::Group> root = new osg::Group;
	osg::ref_ptr<osg::Geode> geode = new osg::Geode;

	//test_0: MFC环境搭建
	//g_viewer->setSceneData(osgDB::readNodeFile("E:\\osgData\\cow.osg"));

	//test_1: 时钟绘制
	//osg::ref_ptr<ClockGeometry> clockGeometryRef = new ClockGeometry;
	//clockGeometryRef->createGeometry();
	//geode->addDrawable(clockGeometryRef.get());
	//root->addChild(geode.get());

	//test_2: Delaunay Triangulation
	osg::ref_ptr<osg::Geode> geode6 = new osg::Geode;
	osg::ref_ptr<DelauneyTriangulation> delaunayTriangulation6 = new DelauneyTriangulation;
	//delaunayTriangulation6->getRegularTeatExmapleData(1000000, 0.5);
	//delaunayTriangulation6->getExmapleVertexSetsFromFile("D:\\demo.dat");
	//delaunayTriangulation6->delauneyTriangulation_OSG();
	//delaunayTriangulation6->delauneyTriangulation_OGRECDT();
	//osg::ref_ptr<osg::Geometry> geometry6 = delaunayTriangulation6->createGeometry();
	//geode6->addDrawable(geometry6.get());

	root->addChild(geode6.get());

	g_viewer->setSceneData(root.get());

	g_viewer->setCameraManipulator(new osgGA::TrackballManipulator);
	g_viewer->setThreadingModel(osgViewer::Viewer::SingleThreaded);

	g_viewer->addEventHandler(new osgGA::StateSetManipulator(g_viewer->getCamera()->getOrCreateStateSet()));	//显示网格
	g_viewer->addEventHandler(new osgViewer::StatsHandler);														//状态信息

	g_finished = false;
	_beginthread(render, 0, NULL);
}
