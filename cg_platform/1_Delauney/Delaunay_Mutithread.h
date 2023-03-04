#pragma once

#include <osg/Geometry>

#include <vector>
#include <unordered_map>
#include <deque>
#include <iterator>
#include <algorithm>
#include <random>

#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#undef max
#undef min
#include "delabella_modified.h"

#define MAX_DOUBLE		((std::numeric_limits<double>::max)());
#define MIN_DOUBLE		((std::numeric_limits<double>::min)());
#define MAX_FLOAT		((std::numeric_limits<double>::max)());
#define MIN_FLOAT		((std::numeric_limits<double>::min)());
#define MAX_INT			((std::numeric_limits<int>::max)());
#define MIN_INT			((std::numeric_limits<int>::min)());
#define MAX_UINT		((std::numeric_limits<unsigned int>::max)());
#define MIN_UINT		((std::numeric_limits<unsigned int>::min)());

//-----------------------------version_20230214-----------------------------------------------------------------------------
//关于使用(zhangchuanpan, 2023.02.14)
//1. 多次重复调用delaunay(),可能会出问题，暂时只支持一次性调用,请一次性插入insert()所有的数据，并调用一次delaunay()
//2. 目前有两种多线程实现方案（实际实现的是方案二），它们的区别如下：
//方案一：对每个分区单独进行三角化，合并操作也是多线程进行，分区间的合并操作互不影响，
//		  三角形可以在每次合并过程中单独输出，而不需要考虑其他正在进行的合并操作。但是该方案，
//        在外围会有一些三角形相交的问题，需要特殊处理；由于合并操作互不相关，单个合并结束后，
//		  内存可以提前释放，该方案对内存要求较低，适合超大数据量，比如上亿？？
//方案二：每个分区也是单独进行三角化，不同的是真正的合并操作是在主线程进行；
//        主线程搜集每个分区中识别出来的风险三角形 / 点，统一做一次三角化；然后和每个分区的网格进行比对，
//        合并恢复出完整的网格；这个方法不需要对外围的三角形相交做特殊处理，但是最终的结果，
//        需要等到最后一步通过主线程来输出。这导致分区部分中间数据内存无法提前释放，对内存的需求可能也会比较大
//目前考虑准确性要求，选择【方案二】（暂不考虑数据量非常大，以致内存不足的情况）
//3. 方案二性能测试情况：
//3.1 测试条件：781W数据点， 主频：3.2GHz， CPU：Intel(R) Core(TM) i5 - 6500 CPU  4核心/4线程，   内存：8GB
//3.2 性能对比：（使用CTime输出日志，粗略测算）
//(1) 3线程：25s左右，其中， 三角化占：13s左右，数据导入导出占：12s左右；
//(2) 单线程：90s左右，其中， 三角化占: 87s左右，数据导出时间占：3s左右；
//4. 关于合并方向
//目前只支持水平左右两个方向的合并，暂不支持上下左右四个方向的合并方式
//5. 关于约束
//暂不支持
//6. 关于面的数量
//不同的分区数量划分，可能导致最终结果的三角形数有略微的差异，几十个？可能原因是由于不同的分区
//划分影响到了分区的上下切线的识别，主要影响外围边界的三角形，对内部无影响（也有可能是其他原因，还在确认排查）

//-----------------------------version_20230301-----------------------------------------------------------------------------
//1. 简化代码（3000行-->1800行），去掉了对dela算法的不必要封装（之前做封装主要是为了可以方便接入其他的三角化算法）
//2. 删除了上下切线的查找逻辑（用三角网格的所有边界点代替）
//3. 修改多线程逻辑（删除了线程池），优化了内存使用和释放（尽量使用已有的预分配连续内存vector或者dela内部已有的内存，去掉了map/set成员变量，在必要时通过局部变量使用），减少数据拷贝操作（在主流程直接调用dela算法接口并将数据索引传递给接口，而不是把数据拷贝传给dela封装类需要时再调用接口）
//4. 增加上下左右四个合并方向（之前只有左右两个方向）的功能
//5. 增加限制多边形区域进行三角划分的功能(只支持一个外轮廓，支持多个内轮廓--内轮廓接口功能已实现，但是界面没有内轮廓识别功能，所以还未接入该功能)
//6. 删除了预采样负载均衡功能，直接根据总区域范围按面积平均划分分区（分区划分：1x3、3x3，5x4,.......，nxm）
//7. 针对JN_3.las点云418W分区（其他两个分区468W、1200W）dela中的split算法效率低的问题，参考CDT对数据进行随机化处理
//8. 增加三角形顶点及索引数据的分区提取功能

//可能存在的问题：
// 1. 还是存在一些重复面
// 2. 未实现带约束的功能
// 3. 如果使用限制区域功能，在限制区域的边界是有一些变形（三角化算法凸包特性决定的），该问题的可能解决办法：
//	3.1 方法一：添加限制区域与登高线的插值点一起进行三角化，必要时还可以在最后删除掉插值点所涉及的面
//	3.2 方法二：先三角化后，再使用限制区域与三角形的位置关系，剔除掉限制区域以外的三角形
//	如果使用方法二，效率不会比现有三角化更快，只会更慢，慢多少要看采用的方法
//	单从效率考虑，方法一更好（看后面实际需求再实现）

//点云性能测试结果：
//（1）JN_3.las(2100W-单线程-限制区域未开启)：用时240s
//（2）JN_3.las(2100W-3线程-限制区域未开启)： 用时60s
//（3）grounds.las(2500W-3线程-限制区域未开启)： 用时90s
//（4）grounds.las(2500W-单线程-限制区域未开启)： 用时195s
//结论：效率提升2-4倍，由于对分区数据进行随机化处理，同时线程的负载可能不均衡，所以不同的数据集提升幅度有差异

//关于限制区域功能：
//限制区域屏蔽掉了部分点，使得后续分区参与三角化的点变少，三角化速度可能提升；
//但也和限制区域和线程的负载均衡有关系，有些线程点数大量减少，而有些线程点数可能基本不变，merge操作是在线程join之后调用，所以也有可能不提升，限制区域屏蔽逻辑也有一定开销，不排除效率降低的可能
//另一个问题：限制区域有可能导致某个线程的点数<3, 此时会将这部分点加入到merge中

//关于分区逻辑：
//1. 根据当前可用线程数，和点云的长宽比，自动计算分区排布

//--------------------------------------------------------------------------------------------

struct Face
{
	unsigned int index1;
	unsigned int index2;
	unsigned int index3;

	Face()
	{
		index1 = MAX_INT;
		index2 = MAX_INT;
		index3 = MAX_INT;
	}

	Face(const Face& face)
	{
		if (face.index1 > face.index2) {
			index1 = face.index1;
			index2 = face.index2;
		}
		else {
			index1 = face.index2;
			index2 = face.index1;
		}
		index3 = face.index3;

		unsigned int tmp;
		if (index3 > index1) {
			tmp = index1;
			index1 = index3;
			index3 = index2;
			index2 = tmp;
		}
		else {
			if (index3 > index2) {
				tmp = index3;
				index3 = index2;
				index2 = tmp;
			}
		}
	}

	Face(unsigned int index1, unsigned int index2, unsigned int index3, bool bRisk = true)
	{
		if (index1 > index2) {
			this->index1 = index1;
			this->index2 = index2;
		}
		else {
			this->index1 = index2;
			this->index2 = index1;
		}
		this->index3 = index3;

		unsigned int tmp;
		if (this->index3 > this->index1) {
			tmp = this->index1;
			this->index1 = this->index3;
			this->index3 = this->index2;
			this->index2 = tmp;
		}
		else {
			if (this->index3 > this->index2) {
				tmp = this->index3;
				this->index3 = this->index2;
				this->index2 = tmp;
			}
		}
	}

	Face& operator=(const Face& face)
	{
		index1 = face.index1;
		index2 = face.index2;
		index3 = face.index3;

		return *this;
	}

	bool operator<(const Face& face) const
	{
		if (index1 != face.index1)
			return index1 < face.index1;

		if (index2 != face.index2)
			return index2 < face.index2;

		if (index3 != face.index3)
			return index3 < face.index3;

		return false;
	}
};

template <typename T = double, typename I = int>
struct CDelaMerge : public CDelaBella2<T, I>
{
	typedef CDelaMerge<T, I> Self;

	struct Point_2;
	typedef struct MyVert : CDelaBella2<T, I>::Vert
	{
		Point_2 point()
		{
			return Point_2(this->x, this->y);
		}
	} MyVert;

	typedef struct MyFace : CDelaBella2<T, I>::Face
	{
		MyVert *vertex(I vi)
		{
			if (vi >= 0 && vi < 3)
				return (MyVert *)this->v[vi];
			return nullptr;
		}
	} MyFace;

	//todo: 封装成迭代器
	typedef MyVert*					Vertex_handle;
	typedef MyFace*					Face_handle;

	struct Point_2
	{
		Point_2()
		{
			x = (T)0;
			y = (T)0;
			m_info = MAX_UINT;
		}

		Point_2(T tx, T ty)
		{
			x = tx;
			y = ty;
			m_info = MAX_UINT;
		}

		Point_2(T tx, T ty, unsigned int info)
		{
			x = tx;
			y = ty;
			m_info = info;
		}

		Point_2(const Point_2& p)
		{
			x = p.x;
			y = p.y;
			m_info = p.m_info;
		}

		Point_2(Vertex_handle vh, unsigned int info)
		{
			x = vh->x;
			y = vh->y;
			m_info = info;
		}
		T x;
		T y;
		unsigned int m_info;
	};

	struct Segment_2
	{
		Segment_2() {}

		Segment_2(Vertex_handle vh1, Vertex_handle vh2)
		{
			if (vh1)
			{
				x1 = vh1->x;
				y1 = vh1->y;
			}

			if (vh2)
			{
				x2 = vh2->x;
				y2 = vh2->y;
			}
		}

		Segment_2(const Point_2 &p1, const Point_2& p2)
		{
			x1 = p1.x;
			y1 = p1.y;

			x2 = p2.x;
			y2 = p2.y;
		}

		Segment_2(const Segment_2& seg)
		{
			x1 = seg.x1;
			y1 = seg.y1;
			x2 = seg.x2;
			y2 = seg.y2;
		}

		Segment_2(T xx1, T yy1, T xx2, T yy2) 
		{
			x1 = xx1;
			y1 = yy1;
			x2 = xx2;
			y2 = yy2;
		}

		T x1;
		T y1;
		T x2;
		T y2;
	};

	struct Circle_2
	{
		Circle_2(Vertex_handle vh1, Vertex_handle vh2, Vertex_handle vh3)
		{
			if (vh1)
			{
				x1 = vh1->x;
				y1 = vh1->y;
			}

			if (vh2)
			{
				x2 = vh2->x;
				y2 = vh2->y;
			}

			if (vh3)
			{
				x3 = vh3->x;
				y3 = vh3->y;
			}

			//检查共线??

			center();
		}

		T x1, y1, x2, y2, x3, y3;

		Point_2 cp;
		T r2;

	private:
		void center()
		{
			//todo: 检查共线??

			T d12 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
			T d23 = (x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3);
			T d31 = (x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1);

			//找到最大角，以降低误差
			int iMaxTheta = 3;
			if (d23 > d12 && d23 > d31)
				iMaxTheta = 1;
			else if (d31 > d12 && d31 > d23)
				iMaxTheta = 2;

			if (iMaxTheta == 1)
			{
				T xx0 = (x1 + x2) / 2;
				T yy0 = (y1 + y2) / 2;
				T m0 = (y2 - y1) * (-1);
				T n0 = (x2 - x1);

				T xx1 = (x1 + x3) / 2;
				T yy1 = (y1 + y3) / 2;
				T m1 = (y3 - y1) * (-1);
				T n1 = (x3 - x1);

				//分母为0？？
				T t = (m1 * yy0 - m1 * yy1 - (xx0 * n1 - xx1 * n1)) / (m0 * n1 - m1 * n0);
				cp.x = xx0 + m0 * t;
				cp.y = yy0 + n0 * t;
			}
			else if (iMaxTheta == 2)
			{
				T xx0 = (x2 + x1) / 2;
				T yy0 = (y2 + y1) / 2;
				T m0 = (y1 - y2) * (-1);
				T n0 = (x1 - x2);

				T xx1 = (x2 + x3) / 2;
				T yy1 = (y2 + y3) / 2;
				T m1 = (y3 - y2) * (-1);
				T n1 = (x3 - x2);

				//分母为0？？
				T t = (m1 * yy0 - m1 * yy1 - (xx0 * n1 - xx1 * n1)) / (m0 * n1 - m1 * n0);
				cp.x = xx0 + m0 * t;
				cp.y = yy0 + n0 * t;
			}
			else
			{
				T xx0 = (x3 + x2) / 2;
				T yy0 = (y3 + y2) / 2;
				T m0 = (y2 - y3) * (-1);
				T n0 = (x2 - x3);

				T xx1 = (x3 + x1) / 2;
				T yy1 = (y3 + y1) / 2;
				T m1 = (y1 - y3) * (-1);
				T n1 = (x1 - x3);

				//分母为0？？
				T t = (m1 * yy0 - m1 * yy1 - (xx0 * n1 - xx1 * n1)) / (m0 * n1 - m1 * n0);
				cp.x = xx0 + m0 * t;
				cp.y = yy0 + n0 * t;
			}

			r2 = (cp.x - x1) * (cp.x - x1) + (cp.y - y1) * (cp.y - y1);
		}
	};

	struct Edge_2
	{
		Edge_2() {}

		Edge_2(const Edge_2 &edge)
		{
			b = edge.b;
			e = edge.e;
		}

		Edge_2(I bb, I ee)
			: b(bb)
			, e(ee)
		{
		}

		I b;
		I e;
	};

	struct Polygon_2 
	{
		Polygon_2() {}

		Polygon_2(std::vector<Point_2>& vecPoints)			//&&  std::move   ????
		{
			m_vecPoints.insert(m_vecPoints.end(), vecPoints.begin(), vecPoints.end());
		}
		std::vector<Point_2> m_vecPoints;
	};

	CDelaMerge()
	{
		m_face_risk_head = nullptr;
		m_vertex_risk_head = nullptr;
	}

	~CDelaMerge() {}

	static CDelaMerge<T, I>* Create();

	I Triangulate(I points, const T* x, const T* y, size_t advance_bytes, I stop = -1) override
	{
		return CDelaBella2<T, I>::Triangulate(points, x, y, advance_bytes, stop);
	}
	
	Face_handle m_face_risk_head;
	Vertex_handle m_vertex_risk_head;
};

struct DMContainer
{
	typedef typename CDelaMerge<double, int>							Dela;
	typedef typename Dela::Face_handle									Face_handle;
	typedef typename Dela::Vertex_handle								Vertex_handle;
	typedef typename Dela::Point_2										Point_2;
	typedef typename Dela::Segment_2									Segment_2;
	typedef typename Dela::Circle_2										Circle_2;
	typedef typename Dela::Polygon_2									Polygon_2;

	typedef CGAL::Exact_predicates_inexact_constructions_kernel			K;

	struct Block
	{
		Block() 
		{
			m_iRow = -1;
			m_iCol = -1;

			unsigned int m_ib_raw = MAX_UINT;
			unsigned int m_ie_raw = MAX_UINT;
			unsigned int m_ib = MAX_UINT;
			unsigned int m_ie = MAX_UINT;

			m_dt = nullptr;

			nei[0] = true;
			nei[1] = true;
			nei[2] = true;
			nei[3] = true;
		}

		void reserve(int num)
		{
		}
		void release()
		{
		}
		size_t size() 
		{
			return m_ie - m_ib;
		}

		size_t raw_size()
		{
			return m_ie_raw - m_ib_raw;
		}

		int m_iRow;
		int m_iCol;

		unsigned int m_ib_raw;		//执行限制区域前的初始分区
		unsigned int m_ie_raw;		//执行限制区域前的初始分区

		unsigned int m_ib;			//执行限制区域后的初始分区
		unsigned int m_ie;			//执行限制区域后的初始分区

		bool nei[4];
		double limits[4];

		Dela *m_dt;
	};

	DMContainer() 
		: m_iRows(-1)
		, m_iCols(-1)
	{
	}

	int insert(std::vector<Point_2> &vecPoints) 
	{
		m_dMaxX =  (-1) * MAX_DOUBLE;
		m_dMinX = MAX_DOUBLE;
		m_dMaxY = (-1) *MAX_DOUBLE;
		m_dMinY = MAX_DOUBLE;

		m_vecPoints.swap(vecPoints);
		for (auto &p : m_vecPoints)
		{
			if (p.x > m_dMaxX)
				m_dMaxX = p.x;
			if (p.x < m_dMinX)
				m_dMinX = p.x;
			if (p.y > m_dMaxY)
				m_dMaxY = p.y;
			if (p.y < m_dMinY)
				m_dMinY = p.y;
		}
		return m_vecPoints.size();
	}

	std::vector<std::pair<unsigned int, unsigned int>> partition(unsigned int ib, unsigned int ie, double* dx, int nx, double* dy, int ny, int offset)
	{
		class CompareX
		{
		public:
			CompareX(double x) :dx(x) {}
			bool operator()(Point_2& p) const
			{
				return p.x <= dx;
			}
		private:
			double dx;
		};

		class CompareY
		{
		public:
			CompareY(double y) :dy(y) {}
			bool operator()(Point_2& p) const
			{
				return p.y <= dy;
			}
		private:
			double dy;
		};

		std::vector<std::pair<unsigned int, unsigned int>> vecPartsX;
		unsigned int ibegin = ib;
		unsigned int iend = ie;
		for (int i = 0; i < nx; ++i) {
			if (!dx) break;

			auto ittrX = std::partition(m_vecPoints.begin() + ibegin, m_vecPoints.begin() + iend, CompareX(*(double*)((char*)dx + i * offset)));
			unsigned int sub_ib = ibegin;
			unsigned int sub_im = std::distance(m_vecPoints.begin(), ittrX);
			unsigned int sub_ie = iend;
			vecPartsX.push_back(std::pair<unsigned int, unsigned int>(sub_ib, sub_im));
			if (i == nx - 1) {
				vecPartsX.push_back(std::pair<unsigned int, unsigned int>(sub_im, sub_ie));
			}
			ibegin = sub_im;
			iend = sub_ie;
		}

		std::vector<std::pair<unsigned int, unsigned int>> vecPartsY;
		if (vecPartsX.empty()) 
		{
			for (int j = 0; j < ny; ++j) {
				if (!dy) break;

				auto ittrY = std::partition(m_vecPoints.begin() + ibegin, m_vecPoints.begin() + iend, CompareY(*(double*)((char*)dy + j * offset)));
				unsigned int sub_ib = ibegin;
				unsigned int sub_im = std::distance(m_vecPoints.begin(), ittrY);
				unsigned int sub_ie = iend;
				vecPartsY.push_back(std::pair<unsigned int, unsigned int>(sub_ib, sub_im));
				if (j == ny - 1) {
					vecPartsY.push_back(std::pair<unsigned int, unsigned int>(sub_im, sub_ie));
				}
				ibegin = sub_im;
				iend = sub_ie;
			}
		}
		else 
		{
			for (int i = 0; i < vecPartsX.size(); ++i) {
				if (!dy) break;

				ibegin = vecPartsX[i].first;
				iend = vecPartsX[i].second;
				for (int j = 0; j < ny; ++j) {
					auto ittrY = std::partition(m_vecPoints.begin() + ibegin, m_vecPoints.begin() + iend, CompareY(*(double*)((char*)dy + j * offset)));
					unsigned int sub_ib = ibegin;
					unsigned int sub_im = std::distance(m_vecPoints.begin(), ittrY);
					unsigned int sub_ie = iend;
					vecPartsY.push_back(std::pair<unsigned int, unsigned int>(sub_ib, sub_im));
					if (j == ny - 1) {
						vecPartsY.push_back(std::pair<unsigned int, unsigned int>(sub_im, sub_ie));
					}
					ibegin = sub_im;
					iend = sub_ie;
				}
			}
		}

		if (vecPartsY.empty())
			return vecPartsX;
		return vecPartsY;
	};

	bool partition(int &iRows, int &iCols, unsigned int ib, unsigned int ie)
	{
		if (iRows <= 0 || iCols <= 0)
			return false;

		m_iRows = iRows;
		m_iCols = iCols;

		double XRange = m_dMaxX - m_dMinX;
		double YRange = m_dMaxY - m_dMinY;
		if (XRange < 10e-6 || YRange < 10e-6)
			return -1;

		double XStep = XRange / m_iCols;
		double YStep = YRange / m_iRows;

		std::vector<double> vecSepXs;
		for (int iCol = 1; iCol < m_iCols; iCol++) {
			double dSepX = m_dMinX + iCol * XStep;
			vecSepXs.push_back(dSepX);
		}

		std::vector<double> vecSepYs;
		for (int iRow = 1; iRow < m_iRows; iRow++) {
			double dSepY = m_dMinY + iRow * YStep;
			vecSepYs.push_back(dSepY);
		}

		if (vecSepXs.empty() && vecSepYs.empty())
			return false;

		//1. 线程分区block 调用partition
		typedef std::pair<unsigned int, unsigned int> RangePair;
		std::vector<RangePair> vecBlocks;
		double *pSepX = nullptr;
		double *pSepY = nullptr;
		if (!vecSepXs.empty())
			pSepX = &vecSepXs[0];
		if (!vecSepYs.empty())
			pSepY = &vecSepYs[0];
		if (m_vecPoints.empty())
			return false;
		vecBlocks = partition(ib, ie, pSepX, vecSepXs.size(), pSepY, vecSepYs.size(), sizeof(double));
		/*std::sort(vecBlocks.begin(), vecBlocks.end(), [](RangePair &p1, RangePair &p2)->bool {
			if (p1.first < p2.first) {
				return true;
			} else if (p1.first == p2.first) {
				if (p1.second < p2.second) {
					return true;
				}
			}
			return false;
		});*/

		if (vecBlocks.size() != m_iCols * m_iRows) 
		{
			//划分异常
			m_vecBlocks.clear();
			return false;
		}

		for (int iCol = 0; iCol < m_iCols; ++iCol) {
			for (int iRow = 0; iRow < m_iRows; ++iRow) {
				int iBlock = iRow + iCol * m_iRows;
				Block newBlock;
				newBlock.m_ib_raw = vecBlocks[iBlock].first;
				newBlock.m_ie_raw = vecBlocks[iBlock].second;
				newBlock.m_ib = vecBlocks[iBlock].first;
				newBlock.m_ie = vecBlocks[iBlock].second;
				get_limit(iRow, iCol, newBlock.limits);

				if (iCol == 0)
					newBlock.nei[2] = false;
				if (iCol == m_iCols - 1)
					newBlock.nei[3] = false;
				if (iRow == 0)
					newBlock.nei[1] = false;
				if (iRow == m_iRows - 1)
					newBlock.nei[0] = false;
				m_vecBlocks.push_back(newBlock);
			}
		}

		//处理空分区和合并点数太少的分区
		/*
		for (int iCol = 0; iCol < m_iCols; ++iCol) {
			for (int iRow = 0; iRow < m_iRows; ++iRow) {
				int iBlock = iRow + iCol * m_iRows;
				if (iBlock < vecBlocks.size()) {
					if (vecBlocks[iBlock].first >= vecBlocks[iBlock].second)
						continue;
					int curSize = vecBlocks[iBlock].second - vecBlocks[iBlock].first;
					bool bNewBlock = false;
					if (iRow == 0) 
						bNewBlock = true;
					if (!m_vecBlocks.empty()) {
						Block &lastBlock = m_vecBlocks.back();
						int lastSize = lastBlock.m_ie - lastBlock.m_ib;
						if (lastSize >= 100000 && curSize >= 100000)
							bNewBlock = true;
					}

					if (bNewBlock) {
						Block newBlock;
						newBlock.m_ib = vecBlocks[iBlock].first;
						newBlock.m_ie = vecBlocks[iBlock].second;
						get_limit(iRow, iCol, newBlock.limits);

						if (iCol == 0)
							newBlock.nei[2] = false;
						if (iCol == m_iCols - 1)
							newBlock.nei[3] = false;
						if (iRow == 0)
							newBlock.nei[1] = false;
						if (iRow == m_iRows - 1)
							newBlock.nei[0] = false;
						m_vecBlocks.push_back(newBlock);
						continue;
					}
					if (!m_vecBlocks.empty()) {
						Block &lastBlock = m_vecBlocks.back();
						lastBlock.m_ie = vecBlocks[iBlock].second;
						double expandLimits[4];
						get_limit(iRow, iCol, expandLimits);
						lastBlock.limits[1] = expandLimits[1];
						if (iRow == m_iRows - 1)
							lastBlock.nei[1] = false;
					}
				}
			}
		}*/

		m_vecvecIndexs.resize(m_vecBlocks.size());

		return true;
	}

	bool hasNeighbor(int iRow, int iCol, int i)
	{
		if (i == 0) 
		{
			return iRow != m_iRows - 1;
		}
		else if (i == 1)
		{
			return iRow != 0;
		}
		else if (i == 2)
		{
			return iCol != 0;
		}
		else if (i == 3)
		{
			return iCol != m_iCols - 1;
		}
		return false;
	}

	void saveIndex(int iBlockIndex, unsigned int index)
	{
		m_vecvecIndexs[iBlockIndex].push_back(index);
	}

	void saveMergeIndex(unsigned int index) 
	{
		m_vecMergeTriIndexs.push_back(index);
	}

	void reserve(int num)
	{
		m_vecPoints.reserve(num);
	}
	void release() 
	{
		std::vector<Point_2>().swap(m_vecPoints);
		for (int i = 0; i < m_vecBlocks.size(); ++i) 
		{
			if (m_vecBlocks[i].m_dt) {
				m_vecBlocks[i].m_dt->Destroy();
				m_vecBlocks[i].m_dt = nullptr;
			}
		}
	}

	bool do_circle_intersect(Face_handle fh, double *limits, bool *nei)
	{
		if (!fh)
			return false;

		Circle_2 cir(fh->vertex(0), fh->vertex(1), fh->vertex(2));
		double dis2 = MAX_UINT;
		for (int i = 0; i < 4; ++i) {
			if ((i == 2 || i == 3) && nei[i]) {
				dis2 = (cir.cp.x - limits[i]) * (cir.cp.x - limits[i]);
			} else if (i == 0 || i == 1 && nei[i]) {
				dis2 = (cir.cp.y - limits[i]) * (cir.cp.y - limits[i]);
			}

			if (cir.r2 >= dis2) {
				return true;
			}
		}
		return false;
	}

	bool do_face_intersect(Face_handle handle, Segment_2 seg)
	{
		Segment_2 seg1(handle->vertex(0), handle->vertex(1));
		Segment_2 seg2(handle->vertex(0), handle->vertex(2));
		Segment_2 seg3(handle->vertex(1), handle->vertex(2));
		if (do_intersect(seg, seg1)
			|| do_intersect(seg, seg2)
			|| do_intersect(seg, seg3))
		{
			return true;
		}
		return false;
	}

#define eps 1.0e-6
	bool do_intersect(const Segment_2 &l1, const Segment_2 &l2)
	{
		//快速排斥实验
		if ((l1.x1 > l1.x2 ? l1.x1 : l1.x2) < (l2.x1 < l2.x2 ? l2.x1 : l2.x2) ||
			(l1.y1 > l1.y2 ? l1.y1 : l1.y2) < (l2.y1 < l2.y2 ? l2.y1 : l2.y2) ||
			(l2.x1 > l2.x2 ? l2.x1 : l2.x2) < (l1.x1 < l1.x2 ? l1.x1 : l1.x2) ||
			(l2.y1 > l2.y2 ? l2.y1 : l2.y2) < (l1.y1 < l1.y2 ? l1.y1 : l1.y2)   )
		{
			return false;
		}
		//跨立实验
		if (   (((l1.x1 - l2.x1)*(l2.y2 - l2.y1) - (l1.y1 - l2.y1)*(l2.x2 - l2.x1))* ((l1.x2 - l2.x1)*(l2.y2 - l2.y1) - (l1.y2 - l2.y1)*(l2.x2 - l2.x1))) > eps
			|| (((l2.x1 - l1.x1)*(l1.y2 - l1.y1) - (l2.y1 - l1.y1)*(l1.x2 - l1.x1))* ((l2.x2 - l1.x1)*(l1.y2 - l1.y1) - (l2.y2 - l1.y1)*(l1.x2 - l1.x1))) > eps )
		{
			return false;
		}
		return true;
	}

	bool do_intersect_cgal(Segment_2 seg1, Segment_2 seg2) 
	{
		typedef CGAL::Segment_2<K>	CGSegment_2;
		typedef K::Point_2			CGPoint_2;

		CGSegment_2 segment1(CGPoint_2(seg1.x1, seg1.y1), CGPoint_2(seg1.x2, seg1.y2));
		CGSegment_2 segment2(CGPoint_2(seg2.x1, seg2.y1), CGPoint_2(seg2.x2, seg2.y2));

		return CGAL::do_intersect(segment1, segment2);
	}

	int getOutlineSegs(Polygon_2 &polygon, std::vector<Segment_2> &vecSegs)
	{
		int size = polygon.m_vecPoints.size();
		if (size <= 2)
			return vecSegs.size();
		for (int j = 0; j < size - 1; ++j) {
			int k = j + 1;
			Segment_2 seg(polygon.m_vecPoints[j], polygon.m_vecPoints[k]);
			vecSegs.push_back(seg);
		}
		Segment_2 seg(polygon.m_vecPoints[size - 1], polygon.m_vecPoints[0]);
		vecSegs.push_back(seg);
		return vecSegs.size();
	}

	//(==0， >0,  <0)  ---> (相交,  内部,  外部)
	int do_quad_intersect(double x1, double y1, double x2, double y2, std::vector<Segment_2> &vecSegs, Polygon_2 &polygon)
	{
		if (vecSegs.empty())
			return 1;

		Segment_2 seg1(x1, y1, x2, y1);
		Segment_2 seg2(x2, y1, x2, y2);
		Segment_2 seg3(x2, y2, x1, y2);
		Segment_2 seg4(x1, y2, x1, y1);

		bool bIntersect = false;
		for (auto &seg : vecSegs)
		{
			if (!bIntersect)
				bIntersect = do_intersect(seg, seg1);
			if (!bIntersect)
				bIntersect = do_intersect(seg, seg2);
			if (!bIntersect)
				bIntersect = do_intersect(seg, seg3);
			if (!bIntersect)
				bIntersect = do_intersect(seg, seg4);

			if (bIntersect)
				return 0;
		}

		//不相交, 仍有几种情况: 分区在轮廓外、分区在轮廓内、 分区和轮廓重叠、分区包含了轮廓，
		//但是由于特殊的划分使包含的情况不成立，所以不考虑这种情况
		if (false == do_quad_point_inside(x1, y1, polygon))
			return -1;
		if (false == do_quad_point_inside(x1, y2, polygon))
			return -1;
		if (false == do_quad_point_inside(x2, y1, polygon))
			return -1;
		if (false == do_quad_point_inside(x2, y2, polygon))
			return -1;

		return 1;		//分区四个点都在轮廓上且不相交--》在轮廓内部
	}

	bool do_quad_point_inside(double x, double y, Polygon_2& poly)
	{
		//射线法：  (点数太多，效率？？)
		// ref1： https://wrfranklin.org/Research/Short_Notes/pnpoly.html
		// ref2： https://erich.realtimerendering.com/ptinpoly/
		auto pnpoly = [&](float testx, float testy)->bool {
			int i, j, c = 0;
			std::vector<Point_2> &vecPoints = poly.m_vecPoints;
			int nvert = vecPoints.size();
			for (i = 0, j = nvert - 1; i < nvert; j = i++) {
				if (((vecPoints[i].y > testy) != (vecPoints[j].y > testy)) &&
					(testx < (vecPoints[j].x - vecPoints[i].x) * (testy - vecPoints[i].y) / (vecPoints[j].y - vecPoints[i].y) + vecPoints[i].x)) {
					c = !c;
				}
			}
			return c;
		};
		return pnpoly(x, y);
	}

	bool get_limit(int iRow, int iCol, double *pLimit) 
	{
		if (!pLimit) 
			return false;

		double XRange = m_dMaxX - m_dMinX;
		double YRange = m_dMaxY - m_dMinY;

		if (XRange < 10e-6 || YRange < 10e-6)
			return false;

		double XStep = XRange / m_iCols;
		double YStep = YRange / m_iRows;

		double dMinX = m_dMinX + iCol * XStep;
		double dMaxX = m_dMinX + (iCol + 1) * XStep;
		double dMinY = m_dMinY + iRow * YStep;
		double dMaxY = m_dMinY + (iRow + 1) * YStep;

		if (iRow == m_iRows - 1)
			dMaxY = m_dMaxY;

		if (iCol == m_iCols - 1)
			dMaxX = m_dMaxX;

		pLimit[0] = dMaxY;
		pLimit[1] = dMinY;
		pLimit[2] = dMinX;
		pLimit[3] = dMaxX;

		return true;
	}

	bool delaunay(int iBlockIndex)
	{
		if (iBlockIndex < 0 || iBlockIndex >= m_vecBlocks.size())
			return false;
		Block &block = m_vecBlocks[iBlockIndex];
		if (!block.m_dt) {
			block.m_dt = Dela::Create();
		}

		if (m_vecPoints.empty())
			return false;

		//JN_3.las点云418W分区（其他两个分区468W、1200W）dela中的split算法效率低的问题，参考CDT对数据进行随机化处理
		static std::mt19937 randGenerator(9001);
		auto random_shuffle = [&](unsigned int ib, unsigned int ie)
		{
			if (ib >= ie) return;
			unsigned int n = ie - ib;
			for (unsigned int i = n; i >= 1; --i)
			{
				std::swap(m_vecPoints[ib + i - 1], m_vecPoints[ib + randGenerator() % (i)]);
			}
		};
		random_shuffle(block.m_ib, block.m_ie);

		//处理限制区域，调用partition
		struct SQuadrans
		{
			unsigned int ib;
			unsigned int ie;

			double x1;
			double x2;
			double y1;
			double y2;
		};
		typedef std::vector<std::pair<unsigned int, unsigned int>> VecBlockRange;
		std::deque<SQuadrans> deQuads;

		std::vector<unsigned int> vecMaskIndexs;
		//从轮廓的外围矩形的中心开始划分，避免出现分区包含整个轮廓的情况
		auto quad_div = [&](SQuadrans& quad, std::vector<Segment_2> &vecSegs, Polygon_2 &polygon, int OutlineType)
		{
			vecMaskIndexs.clear();
			if (quad.ib >= quad.ie)
				return false;

			if (vecSegs.empty() || vecSegs.size() <= 2)
				return false;

			double rectMinx = MAX_DOUBLE;
			double rectMaxx = (-1) * MAX_DOUBLE;
			double rectMiny = MAX_DOUBLE;
			double rectMaxy = (-1) * MAX_DOUBLE;

			for (auto& seg : vecSegs) 
			{
				if (seg.x1 < rectMinx)
					rectMinx = seg.x1;
				if (seg.x2 < rectMinx)
					rectMinx = seg.x2;
				if (seg.x2 > rectMaxx)
					rectMaxx = seg.x2;
				if (seg.x1 > rectMaxx)
					rectMaxx = seg.x1;
				if (seg.y1 < rectMiny)
					rectMiny = seg.y1;
				if (seg.y2 < rectMiny)
					rectMiny = seg.y2;
				if (seg.y2 > rectMaxy)
					rectMaxy = seg.y2;
				if (seg.y1 > rectMaxy)
					rectMaxy = seg.y1;
			}

			if (rectMinx < quad.x1)
				quad.x1 = rectMinx;
			if (rectMaxx > quad.x2)
				quad.x2 = rectMaxx;
			if (rectMiny < quad.y1)
				quad.y1 = rectMiny;
			if (rectMaxy > quad.y2)
				quad.y2 = rectMaxy;

			VecBlockRange vecBR;
			double sepx = (rectMinx + rectMaxx) / 2;
			double sepy = (rectMiny + rectMaxy) / 2;

			vecBR = partition(quad.ib, quad.ie, &sepx, 1, &sepy, 1, sizeof(double));
			if (vecBR.size() != 4)
				return false;

			if (vecBR[0].first < vecBR[0].second) {
				SQuadrans quad1;
				quad1.x1 = quad.x1;
				quad1.y1 = quad.y1;
				quad1.x2 = sepx;
				quad1.y2 = sepy;
				quad1.ib = vecBR[0].first;
				quad1.ie = vecBR[0].second;
				deQuads.push_back(quad1);
			}

			if (vecBR[1].first < vecBR[1].second) {
				SQuadrans quad2;
				quad2.x1 = quad.x1;
				quad2.y1 = sepy;
				quad2.x2 = sepx;
				quad2.y2 = quad.y2;
				quad2.ib = vecBR[1].first;
				quad2.ie = vecBR[1].second;
				deQuads.push_back(quad2);
			}

			if (vecBR[2].first < vecBR[2].second) {
				SQuadrans quad3;
				quad3.x1 = sepx;
				quad3.y1 = quad.y1;
				quad3.x2 = quad.x2;
				quad3.y2 = sepy;
				quad3.ib = vecBR[2].first;
				quad3.ie = vecBR[2].second;
				deQuads.push_back(quad3);
			}

			if (vecBR[3].first < vecBR[3].second) {
				SQuadrans quad4;
				quad4.x1 = sepx;
				quad4.y1 = sepy;
				quad4.x2 = quad.x2;
				quad4.y2 = quad.y2;
				quad4.ib = vecBR[3].first;
				quad4.ie = vecBR[3].second;
				deQuads.push_back(quad4);
			}

			while (!deQuads.empty())
			{
				quad = deQuads.front();
				deQuads.pop_front();
				int retPos = do_quad_intersect(quad.x1, quad.y1, quad.x2, quad.y2, vecSegs, polygon);
				if (retPos == 0) {	//相交
					if ((quad.ie - quad.ib) <= (unsigned int)10) {
						for (unsigned int ui = quad.ib; ui < quad.ie; ++ui) {
							if (false == do_quad_point_inside(m_vecPoints[ui].x, m_vecPoints[ui].y, polygon)) {
								if (OutlineType == 0) {
									vecMaskIndexs.push_back(ui);	//目前只支持一个外轮廓，一个以上的话就不行了
								}
							} else {
								if (OutlineType == 1) {
									vecMaskIndexs.push_back(ui);
								}
							}
						}
					} else {
						VecBlockRange vecBR;
						double sepx = (quad.x1 + quad.x2) / 2;
						double sepy = (quad.y1 + quad.y2) / 2;
						if (quad.ib >= quad.ie)
							return false;

						vecBR = partition(quad.ib, quad.ie, &sepx, 1, &sepy, 1, sizeof(double));
						if (vecBR.size() != 4)
							continue;

						if (vecBR[0].first < vecBR[0].second) {
							SQuadrans quad1;
							quad1.x1 = quad.x1;
							quad1.y1 = quad.y1;
							quad1.x2 = sepx;
							quad1.y2 = sepy;
							quad1.ib = vecBR[0].first;
							quad1.ie = vecBR[0].second;
							deQuads.push_back(quad1);
						}

						if (vecBR[1].first < vecBR[1].second) {
							SQuadrans quad2;
							quad2.x1 = quad.x1;
							quad2.y1 = sepy;
							quad2.x2 = sepx;
							quad2.y2 = quad.y2;
							quad2.ib = vecBR[1].first;
							quad2.ie = vecBR[1].second;
							deQuads.push_back(quad2);
						}

						if (vecBR[2].first < vecBR[2].second) {
							SQuadrans quad3;
							quad3.x1 = sepx;
							quad3.y1 = quad.y1;
							quad3.x2 = quad.x2;
							quad3.y2 = sepy;
							quad3.ib = vecBR[2].first;
							quad3.ie = vecBR[2].second;
							deQuads.push_back(quad3);
						}

						if (vecBR[3].first < vecBR[3].second) {
							SQuadrans quad4;
							quad4.x1 = sepx;
							quad4.y1 = sepy;
							quad4.x2 = quad.x2;
							quad4.y2 = quad.y2;
							quad4.ib = vecBR[3].first;
							quad4.ie = vecBR[3].second;
							deQuads.push_back(quad4);
						}
					}
				} else if (retPos < 0) {//全部在区域外面
					if (OutlineType == 0) {						
						for (unsigned int ui = quad.ib; ui < quad.ie; ++ui) {
							vecMaskIndexs.push_back(ui);
						}
					}
				} else {
					if (OutlineType == 1) {
						for (unsigned int ui = quad.ib; ui < quad.ie; ++ui) {
							vecMaskIndexs.push_back(ui);
						}
					}
				}
			}//while (!deQuads.empty()) end

			//将被屏蔽的点移动到分区最后
			std::sort(vecMaskIndexs.begin(), vecMaskIndexs.end());
			unsigned int iBack = block.m_ie - 1;
			if (!vecMaskIndexs.empty())
			{
				for (int i = vecMaskIndexs.size(); i >= 1; i--)
				{
					unsigned int iMask = vecMaskIndexs[i - 1];
					if (iBack > iMask) {
						std::swap(m_vecPoints[iMask], m_vecPoints[iBack]);
					}
					block.m_ie = iBack;
					iBack--;
				}
			}
			return true;
		};// quad_div definition end
		
		//只支持一个外轮廓
		if (!m_polyOutLines.empty())
		{
			Polygon_2 &polygon = m_polyOutLines[0];
			SQuadrans quadOutline;
			quadOutline.ib = block.m_ib;
			quadOutline.ie = block.m_ie;
			quadOutline.x1 = block.limits[2];
			quadOutline.x2 = block.limits[3];
			quadOutline.y1 = block.limits[1];
			quadOutline.y2 = block.limits[0];
			std::vector<Segment_2> vecSegs;
			getOutlineSegs(polygon, vecSegs);
			quad_div(quadOutline, vecSegs, polygon, 0);
		}

		//多个内轮廓
		for (auto &polygon : m_polyHoles)
		{
			SQuadrans quadOutline;
			quadOutline.ib = block.m_ib;
			quadOutline.ie = block.m_ie;
			quadOutline.x1 = block.limits[2];
			quadOutline.x2 = block.limits[3];
			quadOutline.y1 = block.limits[1];
			quadOutline.y2 = block.limits[0];
			std::vector<Segment_2> vecSegs;
			getOutlineSegs(polygon, vecSegs);
			quad_div(quadOutline, vecSegs, polygon, 1);
		}

		if (block.size() <= 2)
			return true; //这部分点要加入到merge中，见merge

		int iRet = block.m_dt->Triangulate(block.size(), &m_vecPoints[block.m_ib].x, &m_vecPoints[block.m_ib].y, sizeof(Point_2));
		if (iRet <= 0)
			return false;

		//单线程
		if (m_vecBlocks.size() == 1) 
		{
			size_t size = block.m_dt->GetNumPolygons();
			if (size <= 0)
				return false;
			Face_handle fh = static_cast<Face_handle>(const_cast<IDelaBella2<double, int>::Simplex*>(block.m_dt->GetFirstDelaunaySimplex()));
			for (int i = 0; i < size && fh; ++i) 
			{
				unsigned int p1 = block.m_ib + fh->vertex(0)->i;
				unsigned int p2 = block.m_ib + fh->vertex(1)->i;
				unsigned int p3 = block.m_ib + fh->vertex(2)->i;

				saveIndex(iBlockIndex, p1);
				saveIndex(iBlockIndex, p2);
				saveIndex(iBlockIndex, p3);

				fh = (Face_handle)fh->next;
			}
			return true;
		}

		//识别风险面
		size_t size = block.m_dt->GetNumPolygons();
		if (size <= 0)
			return false;
		Face_handle fh = static_cast<Face_handle>(const_cast<IDelaBella2<double, int>::Simplex*>(block.m_dt->GetFirstDelaunaySimplex()));
		for (int i = 0; i < size && fh; ++i)
		{
			Face_handle fh_tmp = (Face_handle)fh->next;
			if (do_circle_intersect(fh, block.limits, block.nei)) {
				if (block.m_dt->m_face_risk_head && block.m_dt->m_face_risk_head != fh) {
					Face_handle tmp = block.m_dt->m_face_risk_head;
					block.m_dt->m_face_risk_head = fh;
					block.m_dt->m_face_risk_head->next = tmp;
				}
				else if (nullptr == block.m_dt->m_face_risk_head) {
					block.m_dt->m_face_risk_head = fh;
					block.m_dt->m_face_risk_head->next = nullptr;
				}
			} else {
				unsigned int p1 = block.m_ib + fh->vertex(0)->i;
				unsigned int p2 = block.m_ib + fh->vertex(1)->i;
				unsigned int p3 = block.m_ib + fh->vertex(2)->i;

				saveIndex(iBlockIndex, p1);
				saveIndex(iBlockIndex, p2);
				saveIndex(iBlockIndex, p3);

				if ((abs(m_vecPoints[p1].x - 1236.929658) < 10e-4 && abs(m_vecPoints[p1].y - 4316.058791) < 10e-4)
					|| (abs(m_vecPoints[p1].x - 1158.188142) < 10e-4 && abs(m_vecPoints[p1].y - 4268.846663) < 10e-4)
					|| (abs(m_vecPoints[p1].x - 1167.637176) < 10e-4 && abs(m_vecPoints[p1].y - 4391.598187) < 10e-4))
				{
					int a = 0;
				}

				if ((abs(m_vecPoints[p2].x - 1236.929658) < 10e-4 && abs(m_vecPoints[p2].y - 4316.058791) < 10e-4)
					|| (abs(m_vecPoints[p2].x - 1158.188142) < 10e-4 && abs(m_vecPoints[p2].y - 4268.846663) < 10e-4)
					|| (abs(m_vecPoints[p2].x - 1167.637176) < 10e-4 && abs(m_vecPoints[p2].y - 4391.598187) < 10e-4))
				{
					int a = 0;
				}

				if ((abs(m_vecPoints[p3].x - 1236.929658) < 10e-4 && abs(m_vecPoints[p3].y - 4316.058791) < 10e-4)
					|| (abs(m_vecPoints[p3].x - 1158.188142) < 10e-4 && abs(m_vecPoints[p3].y - 4268.846663) < 10e-4)
					|| (abs(m_vecPoints[p3].x - 1167.637176) < 10e-4 && abs(m_vecPoints[p3].y - 4391.598187) < 10e-4))
				{
					int a = 0;
				}
			}
			fh = fh_tmp;
		}

		//识别边界点
		int num = block.m_dt->GetNumBoundaryVerts();
		Vertex_handle vh = static_cast<Vertex_handle>(const_cast<IDelaBella2<double, int>::Vertex *>(block.m_dt->GetFirstBoundaryVertex()));
		for (int i = 0; i < num && vh; ++i) {
			Vertex_handle vh_tmp = (Vertex_handle)vh->next;
			if (block.m_dt->m_vertex_risk_head && vh != block.m_dt->m_vertex_risk_head) {
				Vertex_handle tmp = block.m_dt->m_vertex_risk_head;
				block.m_dt->m_vertex_risk_head = vh;
				block.m_dt->m_vertex_risk_head->next = tmp;
			}
			else if (nullptr == block.m_dt->m_vertex_risk_head) {
				block.m_dt->m_vertex_risk_head = vh;
				block.m_dt->m_vertex_risk_head->next = nullptr;
			}
			vh = vh_tmp;
		}
		return true;
	}

	bool merge()
	{
		//merge区域三角化
		Dela *delaTmp = Dela::Create();
		std::vector<Point_2> vecPointsTmp;
		std::map<Face, bool> facesmap;

		for (int iBlock = 0; iBlock < m_vecBlocks.size(); ++iBlock)
		{
			Block &block = m_vecBlocks[iBlock];

			if (block.size() <= 2) 
			{
				unsigned int i = block.m_ib;
				while (i < block.m_ie) {
					Point_2 p(m_vecPoints[i].x, m_vecPoints[i].y, i);
					vecPointsTmp.push_back(p);
					m_vecMergeVertexIndexs.push_back(std::pair<unsigned int, int>(i, iBlock));
					i++;
				}
			}

			if (!block.m_dt)
				continue;

			Face_handle fh = block.m_dt->m_face_risk_head;
			while (fh)
			{
				//构建风险面的索引
				unsigned int pi1 = block.m_ib + fh->vertex(0)->i;
				unsigned int pi2 = block.m_ib + fh->vertex(1)->i;
				unsigned int pi3 = block.m_ib + fh->vertex(2)->i;

				Point_2 p1(fh->vertex(0), pi1);
				Point_2 p2(fh->vertex(1), pi2);
				Point_2 p3(fh->vertex(2), pi3);

				Face face(pi1, pi2, pi3);
				facesmap.insert(std::pair<Face, bool>(face, true));

				vecPointsTmp.push_back(p1);
				vecPointsTmp.push_back(p2);
				vecPointsTmp.push_back(p3);
				m_vecMergeVertexIndexs.push_back(std::pair<unsigned int, int>(pi1, iBlock));
				m_vecMergeVertexIndexs.push_back(std::pair<unsigned int, int>(pi2, iBlock));
				m_vecMergeVertexIndexs.push_back(std::pair<unsigned int, int>(pi3, iBlock));
				fh = (Face_handle)fh->next;
			}

			Vertex_handle vh = block.m_dt->m_vertex_risk_head;
			while (vh)
			{
				unsigned int pi = block.m_ib + vh->i;
				Point_2 p(vh, pi);
				vecPointsTmp.push_back(p);
				m_vecMergeVertexIndexs.push_back(std::pair<unsigned int, int>(pi, iBlock));
				vh = (Vertex_handle)vh->next;
			}

			block.m_dt->Destroy();
			block.m_dt = nullptr;
		}

		if (vecPointsTmp.empty()) 
		{
			delaTmp->Destroy();
			delaTmp = nullptr;
			return true;
		}

		int iRet = delaTmp->Triangulate(vecPointsTmp.size(), &vecPointsTmp[0].x, &vecPointsTmp[0].y, sizeof(Point_2));
		if (iRet <= 0) 
		{
			delaTmp->Destroy();
			delaTmp = nullptr;
			return false;
		}
		
		auto do_face_duplicate = [&](unsigned int info1, unsigned int info2, unsigned int info3)->bool
		{
			Face face(info1, info2, info3);
			auto ittr = facesmap.find(face);
			if (ittr != facesmap.end()) {
				return true;
			}
			return false;
		};

		//输出剩余面
		std::vector<Segment_2> vecSegments;
		getAllSegs(vecSegments);
		size_t size = delaTmp->GetNumPolygons();
		Face_handle fh = static_cast<Face_handle>(const_cast<IDelaBella2<double, int>::Simplex*>(delaTmp->GetFirstDelaunaySimplex()));
		for (int i = 0; i < size; ++i)
		{
			if (fh->vertex(0)->i >= vecPointsTmp.size()
				|| fh->vertex(1)->i >= vecPointsTmp.size()
				|| fh->vertex(2)->i >= vecPointsTmp.size()) 
			{
				fh = (Face_handle)fh->next; 
				continue;
			}

			unsigned int info1 = vecPointsTmp[fh->vertex(0)->i].m_info;
			unsigned int info2 = vecPointsTmp[fh->vertex(1)->i].m_info;
			unsigned int info3 = vecPointsTmp[fh->vertex(2)->i].m_info;

			//saveIndex(0, info1);	//for debug
			//saveIndex(0, info2);	//for debug
			//saveIndex(0, info3);	//for debug

			if ((abs(m_vecPoints[info1].x - 1236.929658) < 10e-4 && abs(m_vecPoints[info1].y - 4316.058791) < 10e-4)
				|| (abs(m_vecPoints[info1].x - 1158.188142) < 10e-4 && abs(m_vecPoints[info1].y - 4268.846663) < 10e-4)
				|| (abs(m_vecPoints[info1].x - 1167.637176) < 10e-4 && abs(m_vecPoints[info1].y - 4391.598187) < 10e-4))
			{
				int a = 0;
			}

			if ((abs(m_vecPoints[info2].x - 1236.929658) < 10e-4 && abs(m_vecPoints[info2].y - 4316.058791) < 10e-4)
				|| (abs(m_vecPoints[info2].x - 1158.188142) < 10e-4 && abs(m_vecPoints[info2].y - 4268.846663) < 10e-4)
				|| (abs(m_vecPoints[info2].x - 1167.637176) < 10e-4 && abs(m_vecPoints[info2].y - 4391.598187) < 10e-4))
			{
				int a = 0;
			}

			if ((abs(m_vecPoints[info3].x - 1236.929658) < 10e-4 && abs(m_vecPoints[info3].y - 4316.058791) < 10e-4)
				|| (abs(m_vecPoints[info3].x - 1158.188142) < 10e-4 && abs(m_vecPoints[info3].y - 4268.846663) < 10e-4)
				|| (abs(m_vecPoints[info3].x - 1167.637176) < 10e-4 && abs(m_vecPoints[info3].y - 4391.598187) < 10e-4))
			{
				int a = 0;
			}

			if (do_face_duplicate(info1, info2, info3)) {
				saveMergeIndex(info1);
				saveMergeIndex(info2);
				saveMergeIndex(info3);
				fh = (Face_handle)fh->next;
				continue;
			}

			for (auto &seg : vecSegments) {
				if (do_face_intersect(fh, seg)) {
					saveMergeIndex(info1);
					saveMergeIndex(info2);
					saveMergeIndex(info3);
					break;
				}
			}
			fh = (Face_handle)fh->next;
		}

		delaTmp->Destroy();
		delaTmp = nullptr;
		return true;
	}

	int getAllSegs(std::vector<Segment_2> &vecSegs) 
	{
		for (int iBlock = 0; iBlock < m_vecBlocks.size(); ++iBlock) 
		{
			Block &block = m_vecBlocks[iBlock];
			if (block.nei[0]) 
			{
				Segment_2 segTop;
				segTop.x1 = block.limits[2];
				segTop.x2 = block.limits[3];
				segTop.y1 = block.limits[0];
				segTop.y2 = block.limits[0];
				vecSegs.push_back(segTop);
			}

			if (block.nei[1])
			{
				Segment_2 segBottom;
				segBottom.x1 = block.limits[2];
				segBottom.x2 = block.limits[3];
				segBottom.y1 = block.limits[1];
				segBottom.y2 = block.limits[1];
				vecSegs.push_back(segBottom);
			}

			if (block.nei[2])
			{
				Segment_2 segLeft;
				segLeft.x1 = block.limits[2];
				segLeft.x2 = block.limits[2];
				segLeft.y1 = block.limits[0];
				segLeft.y2 = block.limits[1];
				vecSegs.push_back(segLeft);
			}

			if (block.nei[3])
			{
				Segment_2 segRight;
				segRight.x1 = block.limits[3];
				segRight.x2 = block.limits[3];
				segRight.y1 = block.limits[0];
				segRight.y2 = block.limits[1];
				vecSegs.push_back(segRight);
			}
		}
		return vecSegs.size();
	}

	int insert_outlines(const std::vector<Polygon_2> &vecPoints)
	{
		m_polyOutLines.insert(m_polyOutLines.end(), vecPoints.begin(), vecPoints.end());
		return vecPoints.size();
	}

	int insert_holes(const std::vector<Polygon_2> &vecPoints)
	{
		m_polyHoles.insert(m_polyHoles.end(), vecPoints.begin(), vecPoints.end());
		return vecPoints.size();
	}

	//data:
	std::vector<Point_2> m_vecPoints;
	std::vector<Block> m_vecBlocks;
	std::vector<Polygon_2> m_polyOutLines;
	std::vector<Polygon_2> m_polyHoles;
	std::vector<std::vector<unsigned int>> m_vecvecIndexs;

	std::vector<std::pair<unsigned int,  int>> m_vecMergeVertexIndexs;		//std::vector<std::pair<vertex_index,  block_index>>
	std::vector<unsigned int> m_vecMergeTriIndexs;
	
	int m_iRows;
	int m_iCols;

	double m_dMaxX;
	double m_dMinX;
	double m_dMaxY;
	double m_dMinY;
};

class Delaunay_Mutithread
{
public:
	typedef typename DMContainer::Point_2		Point_2;
	typedef typename DMContainer::Polygon_2		Polygon_2;
	typedef			 DMContainer::Block			Block;

	Delaunay_Mutithread() {}

	~Delaunay_Mutithread() 
	{
		m_container.release();
	}

	//测试用
	void readInputFromFile(const std::string& fileName, std::vector<Point_2> &vecPoints)
	{
		std::ifstream f(fileName);
		if (!f.is_open())
		{
			return;
		}
		int nVerts;
		int nEdges;
		f >> nVerts >> nEdges;

		std::vector<double> vv;
		vv.reserve(nVerts);
		for (std::size_t i = 0; i < nVerts; ++i)
		{
			double x, y, z;
			f >> x >> y;
			vecPoints.push_back(Point_2(x, y, i));
		}
	}

	osg::ref_ptr<osg::Geometry> createGeometry()
	{
		osg::ref_ptr<osg::Geometry> geo = new osg::Geometry;
		m_vecPointsRef = new osg::Vec3Array;
		m_PrimiRef = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);

		for (auto &data : m_container.m_vecPoints)
		{
			m_vecPointsRef->push_back(osg::Vec3(data.x, data.y, 0.0));
		}

		for (auto &vecIndexs : m_container.m_vecvecIndexs)
		{
			for (auto &index : vecIndexs)
				m_PrimiRef->push_back(index);




		}

		if (!m_vecPointsRef->empty())
			geo->setVertexArray(m_vecPointsRef.get());

		if (!m_PrimiRef->empty())
			geo->addPrimitiveSet(m_PrimiRef.get());

		return geo.release();
	}

	int insert(std::vector<Point_2> &vecPoints) 
	{
		return m_container.insert(vecPoints);
	}

	int insert_outlines(const std::vector<Polygon_2> &vecPoints) 
	{
		return m_container.insert_outlines(vecPoints);
	}

	int insert_holes(const std::vector<Polygon_2> &vecPoints) 
	{
		return m_container.insert_holes(vecPoints);
	}

	int getThreadNum(size_t size)
	{
		int iMaxThread = std::thread::hardware_concurrency();
		if (iMaxThread <= 0)
		{
			iMaxThread = 2;
		}

		//最合适数量: 如果少于1000000，则只用单线程就好了，少数据量用多线程加速效果不一定明显
		int iVertexPerThread = 1000000;

		//测试
		//iVertexPerThread = 3;

		int iThreadNum = 2;
		iThreadNum = size / iVertexPerThread;
		if (iThreadNum <= 1)
			iThreadNum = 2;
		if (iThreadNum >= iMaxThread)
			iThreadNum = iMaxThread;

		return iThreadNum - 1;
	}

	bool getBlockParam(int &iThread, int &iRows, int &iCols, unsigned int ib, unsigned int ie)
	{
		if (iThread == 1) 
		{
			iRows = 1;
			iCols = 1;
			m_container.m_vecBlocks.resize(1);
			m_container.m_vecBlocks[0].m_iCol = 0;
			m_container.m_vecBlocks[0].m_iRow = 0;
			m_container.m_vecBlocks[0].m_ib_raw = ib;
			m_container.m_vecBlocks[0].m_ie_raw = ie;
			m_container.m_vecBlocks[0].m_ib = ib;
			m_container.m_vecBlocks[0].m_ie = ie;
			m_container.m_vecBlocks[0].nei[0] = false;
			m_container.m_vecBlocks[0].nei[1] = false;
			m_container.m_vecBlocks[0].nei[2] = false;
			m_container.m_vecBlocks[0].nei[3] = false;
			m_container.m_vecBlocks[0].limits[0] = m_container.m_dMaxY;
			m_container.m_vecBlocks[0].limits[1] = m_container.m_dMinY;
			m_container.m_vecBlocks[0].limits[2] = m_container.m_dMinX;
			m_container.m_vecBlocks[0].limits[3] = m_container.m_dMaxX;
			m_container.m_vecvecIndexs.resize(1);
			return true;
		}

		unsigned int size = ie - ib;
		double xrange = m_container.m_dMaxX - m_container.m_dMinX;
		double yrange = m_container.m_dMaxY - m_container.m_dMinY;
		double xyratio = xrange / yrange;

		if (xyratio < 0.2) 
		{
			iCols = 1;
			iRows = iThread;
		}
		else if (xyratio > 5) 
		{
			iCols = iThread;
			iRows = 1;
		}
		else 
		{
			int iRowNum = 1;
			int iColNum = 1;
			while (iRowNum *  iColNum < iThread)
			{
				if (xyratio > 1.0)
				{
					iColNum++;
					iRowNum = iColNum / xyratio;
				}
				else
				{
					iRowNum++;
					iColNum = iRowNum * xyratio;
				}
			}
			iCols = iColNum;
			iRows = iRowNum;
		}

		return m_container.partition(iRows, iCols, ib, ie);
	}

	bool delaunay(std::vector<Point_2> &vecPoints)
	{
		//memory enough ??

		auto first_tri = [](DMContainer *pDMC, int iBlockIndex) {
			if (pDMC) {
				pDMC->delaunay(iBlockIndex);
			}
		};

		if (insert(vecPoints) <= 0)
			return false;

		int iThread = getThreadNum(m_container.m_vecPoints.size());

		//测试
		iThread = 10;
		//generate_random_outline(0.7, 0.7, 0.3, 0.3);

		//生成分区
		int iRows = 1;
		int iCols = 1;
		if (!getBlockParam(iThread, iRows, iCols, 0, m_container.m_vecPoints.size()))
		{
			return false;
		}

		//多线程三角化
		iThread = m_container.m_vecBlocks.size();
		if (iThread == 1) 
		{
			first_tri(&m_container, 0);
			return true;
		}
		std::vector<std::thread> vecThreads(iThread);
		for (int iBlock = 0; iBlock < m_container.m_vecBlocks.size(); ++iBlock)
		{
			vecThreads[iBlock] = std::move(std::thread(first_tri, &m_container, iBlock));
		}
		std::for_each(vecThreads.begin(), vecThreads.end(), std::mem_fn(&std::thread::join));
		return m_container.merge();
	}

	//测试用
	void generate_random_outline(double outline_ratio_x, double outline_ratio_y, double hole_ratio_x, double hole_ratio_y)
	{
		double center_x = (m_container.m_dMaxX - m_container.m_dMinX) / 2 - 1000.0;
		double center_y = (m_container.m_dMaxY - m_container.m_dMinY) / 2;
		double ouline_x_range = (m_container.m_dMaxX - m_container.m_dMinX) * outline_ratio_x;
		double ouline_y_range = (m_container.m_dMaxY - m_container.m_dMinY) * outline_ratio_y;
		double hole_x_range = (m_container.m_dMaxX - m_container.m_dMinX) * hole_ratio_x;
		double hole_y_range = (m_container.m_dMaxY - m_container.m_dMinY) * hole_ratio_y;

		double po_x1 = center_x - ouline_x_range / 2;
		double po_y1 = center_y - ouline_y_range / 2;
		double po_x2 = center_x + ouline_x_range / 2;
		double po_y2 = center_y + ouline_y_range / 2;

		double ph_x1 = center_x - hole_x_range / 2;
		double ph_y1 = center_y - hole_y_range / 2;
		double ph_x2 = center_x + hole_x_range / 2;
		double ph_y2 = center_y + hole_y_range / 2;

		Polygon_2 poly_outline;
		std::vector<Polygon_2> vecPolyOutlines;
		poly_outline.m_vecPoints.push_back(Point_2(po_x2, po_y1));
		poly_outline.m_vecPoints.push_back(Point_2(po_x2, po_y2));
		poly_outline.m_vecPoints.push_back(Point_2(po_x1, po_y2));
		poly_outline.m_vecPoints.push_back(Point_2(po_x1, po_y1));

		vecPolyOutlines.push_back(poly_outline);
		insert_outlines(vecPolyOutlines);

		Polygon_2 poly_hole;
		std::vector<Polygon_2> vecPolyHoles;
		poly_hole.m_vecPoints.push_back(Point_2(ph_x2, ph_y1));
		poly_hole.m_vecPoints.push_back(Point_2(ph_x2, ph_y2));
		poly_hole.m_vecPoints.push_back(Point_2(ph_x1, ph_y2));
		poly_hole.m_vecPoints.push_back(Point_2(ph_x1, ph_y1));

		vecPolyHoles.push_back(poly_hole);
		insert_holes(vecPolyHoles);
	}

public:
	DMContainer m_container;

	osg::ref_ptr<osg::Vec3Array> m_vecPointsRef;
	osg::ref_ptr<osg::DrawElementsUInt> m_PrimiRef;
};

template<typename T, typename I>
inline CDelaMerge<T, I>* CDelaMerge<T, I>::Create()
{
	CDelaMerge<T, I> *ret = 0;
	try
	{
		ret = new CDelaMerge<T, I>;
	}
	catch (...)
	{
		ret = 0;
	}
	return ret;
}

