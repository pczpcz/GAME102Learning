#pragma once

#include <osg/Geometry>

#include <vector>
#include <unordered_map>
#include <deque>
#include <iterator>
#include <algorithm>

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
//����ʹ��(zhangchuanpan, 2023.02.14)
//1. ����ظ�����delaunay(),���ܻ�����⣬��ʱֻ֧��һ���Ե���,��һ���Բ���insert()���е����ݣ�������һ��delaunay()
//2. Ŀǰ�����ֶ��߳�ʵ�ַ�����ʵ��ʵ�ֵ��Ƿ������������ǵ��������£�
//����һ����ÿ�����������������ǻ����ϲ�����Ҳ�Ƕ��߳̽��У�������ĺϲ���������Ӱ�죬
//		  �����ο�����ÿ�κϲ������е��������������Ҫ�����������ڽ��еĺϲ����������Ǹ÷�����
//        ����Χ����һЩ�������ཻ�����⣬��Ҫ���⴦�����ںϲ�����������أ������ϲ�������
//		  �ڴ������ǰ�ͷţ��÷������ڴ�Ҫ��ϵͣ��ʺϳ������������������ڣ���
//��������ÿ������Ҳ�ǵ����������ǻ�����ͬ���������ĺϲ������������߳̽��У�
//        ���߳��Ѽ�ÿ��������ʶ������ķ��������� / �㣬ͳһ��һ�����ǻ���Ȼ���ÿ��������������бȶԣ�
//        �ϲ��ָ������������������������Ҫ����Χ���������ཻ�����⴦���������յĽ����
//        ��Ҫ�ȵ����һ��ͨ�����߳���������⵼�·��������м������ڴ��޷���ǰ�ͷţ����ڴ���������Ҳ��Ƚϴ�
//Ŀǰ����׼ȷ��Ҫ��ѡ�񡾷����������ݲ������������ǳ��������ڴ治��������
//3. ���������ܲ��������
//3.1 ����������781W���ݵ㣬 ��Ƶ��3.2GHz�� CPU��Intel(R) Core(TM) i5 - 6500 CPU  4����/4�̣߳�   �ڴ棺8GB
//3.2 ���ܶԱȣ���ʹ��CTime�����־�����Բ��㣩
//(1) 3�̣߳�25s���ң����У� ���ǻ�ռ��13s���ң����ݵ��뵼��ռ��12s���ң�
//(2) ���̣߳�90s���ң����У� ���ǻ�ռ: 87s���ң����ݵ���ʱ��ռ��3s���ң�
//4. ���ںϲ�����
//Ŀǰֻ֧��ˮƽ������������ĺϲ����ݲ�֧�����������ĸ�����ĺϲ���ʽ
//5. ����Լ��
//�ݲ�֧��
//6. �����������
//��ͬ�ķ����������֣����ܵ������ս����������������΢�Ĳ��죬��ʮ��������ԭ�������ڲ�ͬ�ķ���
//����Ӱ�쵽�˷������������ߵ�ʶ����ҪӰ����Χ�߽�������Σ����ڲ���Ӱ�죨Ҳ�п���������ԭ�򣬻���ȷ���Ų飩

//-----------------------------version_20230223-----------------------------------------------------------------------------
//1. �򻯴��루3000��-->1600�У���ȥ���˶�dela�㷨�Ĳ���Ҫ��װ��֮ǰ����װ��Ҫ��Ϊ�˿��Է���������������ǻ��㷨��
//2. ɾ�����������ߵĲ����߼�����������������б߽����棩
//3. �޸Ķ��߳��߼���ɾ�����̳߳أ����Ż����ڴ�ʹ�ú��ͷţ�����ʹ�����е�Ԥ���������ڴ�vector����dela�ڲ����е��ڴ棬ȥ����map/set��Ա�������ڱ�Ҫʱͨ���ֲ�����ʹ�ã����������ݿ�����������������ֱ�ӵ���dela�㷨�ӿڲ��������������ݸ��ӿڣ������ǰ����ݿ�������dela��װ����Ҫʱ�ٵ��ýӿڣ�
//4. �������������ĸ��ϲ�����֮ǰֻ�������������򣩵Ĺ���
//5. �������ƶ��������������ǻ��ֵĹ���(ֻ֧��һ����������֧�ֶ��������)
//6. ɾ����Ԥ�������ؾ��⹦�ܣ�ֱ�Ӹ���������Χ�����ƽ�����ַ������������֣�1x3��3x3��5x4,.......��nxm��

//���ܴ��ڵ����⣺
// 1. ���Ǵ���һЩ�ظ���
// 2. δʵ�ִ�Լ���Ĺ���
// 3. ���ʹ�����������ܣ�����������ı߽�����һЩ���Σ����ǻ��㷨͹�����Ծ����ģ���������Ŀ��ܽ���취��
//	3.1 ����һ���������������Ǹ��ߵĲ�ֵ��һ��������ǻ�����Ҫʱ�����������ɾ������ֵ�����漰����
//	3.2 �������������ǻ�����ʹ�����������������ε�λ�ù�ϵ���޳����������������������
//	���ʹ�÷�������Ч�ʲ�����������ǻ����죬ֻ�������������Ҫ�����õķ���
//	����Ч�ʿ��ǣ�����һ���ã�������ʵ��������ʵ�֣�

//���ܲ��Խ����
//(1) 3�̣߳�781W���ݵ㣬���������ܣ�����������							��ʱ15s����
//(2) 3�̣߳�781W���ݵ㣬���������ܣ� ������0.7�� 0.7�� 0.3�� 0.3������	��ʱ10s����
//�����������ε��˲��ֵ㣬ʹ�ú��������������ǻ��ĵ���٣����ǻ��ٶȿ���������
//��Ҳ������������̵߳ĸ��ؾ����й�ϵ����Щ�̵߳����������٣�����Щ�̵߳������ܻ������䣬merge���������߳�join֮����ã�����Ҳ�п��ܲ��������������������߼�Ҳ��һ�����������ų�Ч�ʽ��͵Ŀ���
//��һ�����⣺���������п��ܵ���ĳ���̵߳ĵ���<3, ��ʱ�Ὣ�ⲿ�ֵ���뵽merge��

//���ڷ����߼���
//1. ���ĳ�е����������������е�ĳ����������С��100000�����ϲ������е������������棨���ײ���������ϲ������ݲ�֧��ˮƽ�����ĺϲ�
//2. ���ݵ�ǰ�����߳������͵��Ƶĳ���ȣ��Զ���������Ų�

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

	//todo: ��װ�ɵ�����
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

			//��鹲��??

			center();
		}

		T x1, y1, x2, y2, x3, y3;

		Point_2 cp;
		T r2;

	private:
		void center()
		{
			//todo: ��鹲��??

			T d12 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
			T d23 = (x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3);
			T d31 = (x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1);

			//�ҵ����ǣ��Խ������
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

				//��ĸΪ0����
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

				//��ĸΪ0����
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

				//��ĸΪ0����
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

		int m_iRow;
		int m_iCol;
		unsigned int m_ib;
		unsigned int m_ie;

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
		m_dMaxX = MIN_DOUBLE;
		m_dMinX = MAX_DOUBLE;
		m_dMaxY = MIN_DOUBLE;
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

		//1. �̷߳���block ����partition
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
			//�����쳣
			m_vecBlocks.clear();
			return false;
		}

		//����շ����ͺϲ�����̫�ٵķ���
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
		}

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

	bool do_intersect(Segment_2 seg1, Segment_2 seg2) 
	{
		typedef CGAL::Segment_2<K>	CGSegment_2;
		typedef K::Point_2			CGPoint_2;

		CGSegment_2 segment1(CGPoint_2(seg1.x1, seg1.y1), CGPoint_2(seg1.x2, seg1.y2));
		CGSegment_2 segment2(CGPoint_2(seg2.x1, seg2.y1), CGPoint_2(seg2.x2, seg2.y2));

		if (CGAL::do_intersect(segment1, segment2))
		{
			return true;
		}
		return false;
	}

	int getOutlineSegs(std::vector<Segment_2> &vecSegs) 
	{
		for (int i = 0; i < m_polyOutLines.size(); ++i)
		{
			int size = m_polyOutLines[i].m_vecPoints.size();
			if (size <= 2)
				continue;
			for (int j = size - 1; j >= 1 ; --j) {
				int k = j - 1;
				Segment_2 seg_ccw(m_polyOutLines[i].m_vecPoints[j], m_polyOutLines[i].m_vecPoints[k]);
				vecSegs.push_back(seg_ccw);
			}
			Segment_2 seg_ccw(m_polyOutLines[i].m_vecPoints[0], m_polyOutLines[i].m_vecPoints[size-1]);
			vecSegs.push_back(seg_ccw);
		}
		return vecSegs.size();
	}

	int getHoleSegs(std::vector<Segment_2> &vecSegs) 
	{
		for (int i = 0; i < m_polyHoles.size(); ++i)
		{
			int size = m_polyHoles[i].m_vecPoints.size();
			if (size <= 2)
				continue;
			for (int j = size - 1; j >= 1; --j) {
				int k = j - 1;
				Segment_2 seg_ccw(m_polyHoles[i].m_vecPoints[j], m_polyHoles[i].m_vecPoints[k]);
				vecSegs.push_back(seg_ccw);
			}
			Segment_2 seg_ccw(m_polyHoles[i].m_vecPoints[0], m_polyHoles[i].m_vecPoints[size - 1]);
			vecSegs.push_back(seg_ccw);
		}
		/*
		for (int i = 0; i < m_polyHoles.size(); ++i)
		{
			int size = m_polyHoles[i].m_vecPoints.size();
			if (size <= 2)
				continue;
			for (int j = 0; j < size - 1; ++j) {
				int k = j + 1;
				Segment_2 seg_cw(m_polyHoles[i].m_vecPoints[j], m_polyHoles[i].m_vecPoints[k]);
				vecSegs.push_back(seg_cw);
			}
			Segment_2 seg_cw(m_polyHoles[i].m_vecPoints[size - 1], m_polyHoles[i].m_vecPoints[0]);
			vecSegs.push_back(seg_cw);
		}
		*/
		return vecSegs.size();
	}

	//(==0�� >0,  <0)  ---> (�ཻ,  �ڲ�,  �ⲿ)
	int do_quad_intersect(double x1, double y1, double x2, double y2, std::vector<Segment_2> &vecSegs)
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

		//���ཻ�����м�������������������⡢�����������ڡ� �����������ص�������������������
		//������������Ļ���ʹ��������������������Բ������������
		for (auto& seg : vecSegs)
		{
			if (false == do_quad_point_inside(x1, y1, vecSegs))
			{
				return -1;
			}
		}

		return 1;
	}

	bool do_quad_point_inside(double x, double y, std::vector<Segment_2> &vecSegs)
	{
		for (auto &seg : vecSegs)
		{
			//int ori = predicates::adaptive::orient2d(seg.x1, seg.y1, seg.x2, seg.y2, x, y);
			//if (ori <= 0)
				//return false;

			typedef K::Point_2			CGPoint_2;
			CGPoint_2 p1(seg.x1, seg.y1);
			CGPoint_2 p2(seg.x2, seg.y2);
			CGPoint_2 p(x, y);

			if (CGAL::RIGHT_TURN == CGAL::orientation(p1, p2, p)) 
			{
				return false;
			}
		}

		return true;
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

		//�����������򣬵���partition
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

		std::vector<Segment_2> vecOutlineSegs;
		std::vector<Segment_2> vecHoleSegs;
		std::vector<unsigned int> vecMaskIndexs;
		getOutlineSegs(vecOutlineSegs);
		getHoleSegs(vecHoleSegs);   //??

		//����������Χ���ε����Ŀ�ʼ���֣�������ַ��������������������
		auto quad_div = [&](SQuadrans& quad, std::vector<Segment_2>& vecSegs, int OutlineType)
		{
			if (quad.ib >= quad.ie)
				return false;

			double rectMinx = MAX_DOUBLE;
			double rectMaxx = MIN_DOUBLE;
			double rectMiny = MAX_DOUBLE;
			double rectMaxy = MIN_DOUBLE;

			for (auto& seg : vecSegs) 
			{
				if (seg.x1 < rectMinx)
					rectMinx = seg.x1;
				if (seg.x2 < rectMinx)
					rectMinx = seg.x1;
				if (seg.x2 > rectMaxx)
					rectMaxx = seg.x1;
				if (seg.x1 > rectMaxx)
					rectMaxx = seg.x1;
				if (seg.y1 < rectMinx)
					rectMiny = seg.x1;
				if (seg.y2 < rectMinx)
					rectMiny = seg.x1;
				if (seg.y2 < rectMinx)
					rectMaxy = seg.x1;
				if (seg.y1 < rectMinx)
					rectMaxy = seg.x1;
			}

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
				int retPos = do_quad_intersect(quad.x1, quad.y1, quad.x2, quad.y2, vecSegs);
				if (retPos == 0) {	//�ཻ
					if ((quad.ie - quad.ib) <= (unsigned int)10) {
						for (unsigned int ui = quad.ib; ui < quad.ie; ++ui) {
							if (false == do_quad_point_inside(m_vecPoints[ui].x, m_vecPoints[ui].y, vecOutlineSegs)) {
								if (OutlineType == 0)
									vecMaskIndexs.push_back(ui);	//Ŀǰֻ֧��һ����������һ�����ϵĻ��Ͳ�����
							} else {
								if (OutlineType == 1)
									vecMaskIndexs.push_back(ui);
							}
						}
					}
					else {
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
				} else if (retPos < 0) {//ȫ������������
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

			unsigned int iBack = block.m_ie - 1;
			if (!vecMaskIndexs.empty())
			{
				for (int i = vecMaskIndexs.size() - 1; i >= 0; i--)
				{
					unsigned int iMask = vecMaskIndexs[i];
					if (iBack > iMask) {
						std::swap(m_vecPoints[iMask], m_vecPoints[iBack]);
					}
					block.m_ie = iBack;
					iBack--;
				}
			}
			vecMaskIndexs.clear();
			return true;
		};// quad_div definition end
		
		//ֻ֧��һ��������
		SQuadrans quadOutline;
		quadOutline.ib = block.m_ie;
		quadOutline.ie = block.m_ib;
		quadOutline.x1 = block.limits[2];
		quadOutline.x2 = block.limits[3];
		quadOutline.y1 = block.limits[1];
		quadOutline.y2 = block.limits[0];
		quad_div(quadOutline, vecOutlineSegs, 0);

		//todo: ���������
		SQuadrans quadHole;
		quadHole.ib = block.m_ib;
		quadHole.ie = block.m_ie;
		quadHole.x1 = block.limits[2];
		quadHole.x2 = block.limits[3];
		quadHole.y1 = block.limits[1];
		quadHole.y2 = block.limits[0];
		quad_div(quadHole, vecHoleSegs, 1);

		if (block.size() <= 2)
			return true; //�ⲿ�ֵ�Ҫ���뵽merge�У���merge

		int iRet = block.m_dt->Triangulate(block.size(), &m_vecPoints[block.m_ib].x, &m_vecPoints[block.m_ib].y, sizeof(Point_2));
		if (iRet <= 0)
			return false;

		//ʶ�������
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
			}
			fh = fh_tmp;
		}

		//ʶ��߽��
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
		//merge�������ǻ�
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
					vecPointsTmp.push_back(m_vecPoints[i]);
					i++;
				}
			}

			if (!block.m_dt)
				continue;

			Face_handle fh = block.m_dt->m_face_risk_head;
			while (fh)
			{
				//���������������
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
				fh = (Face_handle)fh->next;
			}

			Vertex_handle vh = block.m_dt->m_vertex_risk_head;
			while (vh)
			{
				unsigned int pi = block.m_ib + vh->i;
				Point_2 p(vh, pi);
				vecPointsTmp.push_back(p);
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

		//���ʣ����
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

			if (do_face_duplicate(info1, info2, info3)) {
				saveIndex(0, info1);
				saveIndex(0, info2);
				saveIndex(0, info3);
				fh = (Face_handle)fh->next;
				continue;
			}

			for (auto &seg : vecSegments) {
				if (do_face_intersect(fh, seg)) {
					saveIndex(0, info1);
					saveIndex(0, info2);
					saveIndex(0, info3);
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

	Delaunay_Mutithread() {}

	~Delaunay_Mutithread() 
	{
		m_container.release();
	}

	//������
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

		//���������: �������300000����ֻ�õ��߳̾ͺ��ˣ����������ö��̼߳���Ч����һ������
		int iVertexPerThread = 300000;

		//����
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
		iRows = 1;		//for debug
		iCols = 3;		//for debug
		iThread = 3;	//for debug

		return m_container.partition(iRows, iCols, ib, ie);
	}

	bool delaunay(std::vector<Point_2> &vecPoints)
	{
		//memory enough ??

		if (insert(vecPoints) <= 0)
			return false;

		int iThread = getThreadNum(m_container.m_vecPoints.size());

		//����
		generate_random_outline(0.7, 0.7, 0.3, 0.3);

		//���̷߳������ǻ�
		auto first_tri = [](DMContainer *pDMC, int iBlockIndex) {
			if (pDMC) {
				pDMC->delaunay(iBlockIndex);
			}
		};

		int iRows = 1;
		int iCols = 1;
		if (!getBlockParam(iThread, iRows, iCols, 0, m_container.m_vecPoints.size()))
		{
			return false;
		}

		iThread = m_container.m_vecBlocks.size();
		std::vector<std::thread> vecThreads(iThread);
		for (int iBlock = 0; iBlock < m_container.m_vecBlocks.size(); ++iBlock) {
			vecThreads[iBlock] = std::move(std::thread(first_tri, &m_container, iBlock));
		}
		std::for_each(vecThreads.begin(), vecThreads.end(), std::mem_fn(&std::thread::join));

		return m_container.merge();
	}

	//������
	void generate_random_outline(double outline_ratio_x, double outline_ratio_y, double hole_ratio_x, double hole_ratio_y)
	{
		double center_x = (m_container.m_dMaxX - m_container.m_dMinX) / 2;
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
		poly_outline.m_vecPoints.push_back(Point_2(po_x1, po_y1));
		poly_outline.m_vecPoints.push_back(Point_2(po_x1, po_y2));
		poly_outline.m_vecPoints.push_back(Point_2(po_x2, po_y2));
		poly_outline.m_vecPoints.push_back(Point_2(po_x2, po_y1));
		vecPolyOutlines.push_back(poly_outline);
		insert_outlines(vecPolyOutlines);

		Polygon_2 poly_hole;
		std::vector<Polygon_2> vecPolyHoles;
		poly_hole.m_vecPoints.push_back(Point_2(ph_x1, ph_y1));
		poly_hole.m_vecPoints.push_back(Point_2(ph_x1, ph_y2));
		poly_hole.m_vecPoints.push_back(Point_2(ph_x2, ph_y2));
		poly_hole.m_vecPoints.push_back(Point_2(ph_x2, ph_y1));
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
