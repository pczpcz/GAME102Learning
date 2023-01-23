#pragma once

#include <iostream>
#include <thread>
#include <queue>
#include <limits>
#include <mutex>

#include <unordered_map>
#include <map>

#include <osg/Geometry>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Circle_2.h>

#include "ThreadPool.h"

//#include "..\ntm_tin_producing_context.hpp"

#define MAX_DOUBLE		((std::numeric_limits<double>::max)());
#define MIN_DOUBLE		((std::numeric_limits<double>::min)());
#define MAX_FLOAT		((std::numeric_limits<double>::max)());
#define MIN_FLOAT		((std::numeric_limits<double>::min)());
#define MAX_INT			((std::numeric_limits<double>::max)());
#define MIN_INT			((std::numeric_limits<double>::min)());

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
		} else {
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
		} else {
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
		} else {
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
		} else {
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

class Delaunay_Mutithread_CGAL
{
public:
	typedef CGAL::Exact_predicates_inexact_constructions_kernel						K;
	typedef CGAL::Delaunay_mesh_vertex_base_2<K>									Vertex_base;
	typedef CGAL::Delaunay_mesh_face_base_2<K>										Face_base;
	typedef Face_base																Fb;
	typedef CGAL::Triangulation_vertex_base_2<K>									Vb;
	typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K, Vb>		Info;
	typedef CGAL::Triangulation_data_structure_2<Info, Fb>							TDS;
	typedef CGAL::Exact_predicates_tag												Itag;
	//typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>				CDT2D;
	typedef CGAL::Delaunay_triangulation_2<K, TDS>									CDT2D;
	typedef CDT2D::Finite_face_handles												Face_handles;
	typedef CDT2D::Finite_faces_iterator											Face_ittr;
	typedef CDT2D::Finite_vertices_iterator											Vertex_ittr;
	typedef CDT2D::Face_circulator													Circulator;
	typedef CDT2D::Vertex_circulator												VCirculator;
	typedef CDT2D::Face_handle														Face_handle;
	typedef CDT2D::Vertex_handle													Vertex_handle;
	typedef CDT2D::Point															Point2D;
	typedef K::Point_2																Point_2;
	typedef CGAL::Segment_2<K>														Segment_2;
	typedef CGAL::Circle_2<K>														Circle_2;
	/*
	typedef CGAL::Projection_traits_xy_3<K>  Gt;
	typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
	typedef K::Point_3   Point;
	typedef Delaunay::Finite_face_handles Face_handles;
	typedef Delaunay::Face_handle face_handle;
	typedef Delaunay::Vertex_handle Vertex_handle;

	typedef Delaunay::Finite_face_handles Face_handles;
	typedef Delaunay::Face_handle face_handle;
	typedef Delaunay::Vertex_handle Vertex_handle;*/

	Delaunay_Mutithread_CGAL();

	~Delaunay_Mutithread_CGAL();

	void read_data();

	void Delaunay_Mutithread_CGAL::readInputFromFile(const std::string& fileName);

	int getRegularTeatExmapleData(int iVertexNums, float fRatio);

	osg::ref_ptr<osg::Geometry> createGeometry();

	void delaunay();

	int getThreadNum(size_t size);

	void delaunay(osg::Vec3Array::iterator begin, osg::Vec3Array::iterator end);

	void saveFaceIndexs(Face_handle fh);

	void saveFaceIndexs(unsigned int index1, unsigned int index2, unsigned int index3);

	struct SSectionData
	{
		enum Orientation
		{
			Near_Left,
			Near_Right,
			Not_Near
		};

		enum Status
		{
			Status_Idle,
			Status_Tri,
			Status_Merging,
			Status_Finish
		};

		//分区编号，用于判断是否全部合并结束，或者用于确定分区范围（由采样决定）
		int m_iSampleBeginX;
		int m_iSampleEndX;

		Status m_status;

		bool m_bLeftFinishMerge;
		bool m_bRightFinishMerge;

		//boundingBox
		double m_dLeft;
		double m_dRight;
		double m_dTop;
		double m_dBottom;

		//左右极限的顶点
		Vertex_handle vertex_lr[2];

		typedef osg::Vec3Array::iterator DataIttr;
		std::vector<DataIttr> m_vecDataIttrs;

		std::set<Face_handle> m_setRiskHandles_left;
		std::set<Face_handle> m_setRiskHandles_right;

		std::deque<Face_handle> m_deqBarrierFaces;

		std::unordered_map<unsigned int, Vertex_handle> m_vertexHashtable;
		std::map<Face, bool> m_facesmap;

		Delaunay_Mutithread_CGAL::CDT2D m_dt;

		SSectionData()
		{
			m_bLeftFinishMerge = false;
			m_bRightFinishMerge = false;

			m_iSampleBeginX = -1;
			m_iSampleEndX = -1;

			m_status = Status_Idle;

			m_dLeft = MAX_DOUBLE;
			m_dRight = (-1) * MAX_DOUBLE;
			m_dTop = (-1) * MAX_DOUBLE;
			m_dBottom = MAX_DOUBLE;
		}

		SSectionData(const SSectionData &data)
		{
			m_bLeftFinishMerge = false;
			m_bRightFinishMerge = false;

			m_iSampleBeginX = data.m_iSampleBeginX;
			m_iSampleEndX = data.m_iSampleEndX;

			m_status = Status_Idle;

			m_dLeft = MAX_DOUBLE;
			m_dRight = (-1) * MAX_DOUBLE;
			m_dTop = (-1) * MAX_DOUBLE;
			m_dBottom = MAX_DOUBLE;
		}

		SSectionData(const SSectionData& data1, const SSectionData& data2)
		{
			m_bLeftFinishMerge = false;
			m_bRightFinishMerge = false;

			m_iSampleBeginX = data1.m_iSampleBeginX;
			m_iSampleEndX = data2.m_iSampleEndX;

			m_status = Status_Idle;

			m_dLeft = MAX_DOUBLE;
			m_dRight = (-1) * MAX_DOUBLE;
			m_dTop = (-1) * MAX_DOUBLE;
			m_dBottom = MAX_DOUBLE;

			m_dLeft = (data1.m_dLeft < data2.m_dLeft) ? data1.m_dLeft : data2.m_dLeft;
			m_dRight = (data1.m_dRight > data2.m_dRight) ? data1.m_dRight : data2.m_dRight;
			m_dTop = (data1.m_dTop > data2.m_dTop) ? data1.m_dTop : data2.m_dTop;
			m_dBottom = (data1.m_dBottom < data2.m_dBottom) ? data1.m_dBottom : data2.m_dBottom;
		}

		~SSectionData()
		{
		}

		void update_right_most_vertex(Face_handle fh)
		{
			Vertex_handle vh;
			if (fh->vertex(0)->point().x() > fh->vertex(1)->point().x()) {
				vh = fh->vertex(0);
			} else {
				vh = fh->vertex(1);
			}

			if (vh->point().x() < fh->vertex(2)->point().x()) {
				vh = fh->vertex(2);
			}

			if (vertex_lr[1] == Vertex_handle() 
				|| vertex_lr[1]->point().x() < vh->point().x()) 
			{
				vertex_lr[1] = vh;
			}
		}

		void update_left_most_vertex(Face_handle fh)
		{
			Vertex_handle vh;
			if (fh->vertex(0)->point().x() < fh->vertex(1)->point().x()) {
				vh = fh->vertex(0);
			} else {
				vh = fh->vertex(1);
			}

			if (vh->point().x() > fh->vertex(2)->point().x()) {
				vh = fh->vertex(2);
			}

			if (vertex_lr[0] == Vertex_handle()
				|| vertex_lr[0]->point().x() > vh->point().x())
			{
				vertex_lr[0] = vh;
			}
		}

		bool is_border(SSectionData& section, Vertex_handle vh)
		{
			Circulator fcc = section.m_dt.incident_faces(vh);
			Circulator fcc_end(fcc);
			do {
				if (section.m_dt.is_infinite(fcc)){
					return true;
				}
			} while (++fcc != fcc_end);
			return false;
		}

		Vertex_handle get_ccw_first_vertex(Vertex_handle vh_base, std::vector<Vertex_handle> vecVCandidates)
		{
			if (vecVCandidates.empty())
				return Vertex_handle();

			Vertex_handle vh_target = vecVCandidates[0];
			for (int i = 1; i < vecVCandidates.size(); ++i) 
			{
				CGAL::Orientation ori = CGAL::orientation(vh_base->point(), vh_target->point(), vecVCandidates[i]->point());
				if (CGAL::RIGHT_TURN == ori) {
					vh_target = vecVCandidates[i];
				} else if (CGAL::COLLINEAR == ori) 
				{
					//如果共线选择最近的那个
					CGAL::Comparison_result ret = CGAL::compare_distance_to_point(vh_base->point(), vh_target->point(), vecVCandidates[i]->point());
					if (CGAL::LARGER == ret) {
						vh_target = vecVCandidates[i];
					}
				}
			}
			return vh_target;
		}

		Vertex_handle get_cw_first_vertex(Vertex_handle vh_base, std::vector<Vertex_handle> &vecVCandidates)
		{
			if (vecVCandidates.empty())
				return Vertex_handle();

			Vertex_handle vh_target = vecVCandidates[0];
			for (int i = 1; i < vecVCandidates.size(); ++i)
			{
				CGAL::Orientation ori = CGAL::orientation(vh_base->point(), vh_target->point(), vecVCandidates[i]->point());
				if (CGAL::LEFT_TURN == ori) {
					vh_target = vecVCandidates[i];
				}
				else if (CGAL::COLLINEAR == ori)
				{
					//如果共线选择最近的那个
					CGAL::Comparison_result ret = CGAL::compare_distance_to_point(vh_base->point(), vh_target->point(), vecVCandidates[i]->point());
					if (CGAL::LARGER == ret) {
						vh_target = vecVCandidates[i];
					}
				}
			}
			return vh_target;
		}

		//遍历到上切线，去重后，将顶点加入到三角化中(zig_zag)
		void get_up_tangent_vetex(SSectionData& leftGraph, SSectionData& rightGraph, std::vector<Vertex_handle>& vecVHandles)
		{
			Vertex_handle hvl = leftGraph.vertex_lr[1];
			Vertex_handle hvr = rightGraph.vertex_lr[0];
			if (hvl == Vertex_handle() || hvr == Vertex_handle())
				return;

			std::set<Vertex_handle> VertexIndexsSet;
			VertexIndexsSet.insert(hvl);
			VertexIndexsSet.insert(hvr);
			vecVHandles.push_back(hvl);
			vecVHandles.push_back(hvr);
			std::vector<Vertex_handle> vecVCandidates;
			bool bLeftUpdate = true;
			bool bRightUpdate = true;
			while (bLeftUpdate || bRightUpdate)
			{
				if (bLeftUpdate || bRightUpdate)
				{
					if (leftGraph.m_dt.is_infinite(hvl))
						continue;

					VCirculator ccl = leftGraph.m_dt.incident_vertices(hvl);
					VCirculator ccl_end(ccl);
					vecVCandidates.clear();
					do {
						if (!leftGraph.m_dt.is_infinite(ccl->handle()) && is_border(leftGraph, ccl->handle())) {
							if (VertexIndexsSet.find(ccl->handle()) == VertexIndexsSet.end()) {
								if (CGAL::RIGHT_TURN != CGAL::orientation(hvl->point(), hvr->point(), ccl->handle()->point())) {
									vecVCandidates.push_back(ccl->handle());
								}
							}
						}
					} while (++ccl != ccl_end);

					Vertex_handle vh_next_border = get_ccw_first_vertex(hvl, vecVCandidates);
					if (Vertex_handle() != vh_next_border) {
						vecVHandles.push_back(vh_next_border);
						VertexIndexsSet.insert(vh_next_border);
						hvl = vh_next_border;
						bLeftUpdate = true;
					} else {
						bLeftUpdate = false;
					}
				}

				if (bLeftUpdate || bRightUpdate)
				{
					if (rightGraph.m_dt.is_infinite(hvr))
						continue;

					vecVCandidates.clear();
					VCirculator ccr = rightGraph.m_dt.incident_vertices(hvr);
					VCirculator ccr_end(ccr);
					do {
						if (!rightGraph.m_dt.is_infinite(ccr->handle()) && is_border(rightGraph, ccr->handle())) {
							if (VertexIndexsSet.find(ccr->handle()) == VertexIndexsSet.end()) {
								if (CGAL::RIGHT_TURN != CGAL::orientation(hvl->point(), hvr->point(), ccr->handle()->point())) {
									vecVCandidates.push_back(ccr->handle());
								}
							}
						}
					} while (++ccr != ccr_end);

					Vertex_handle vh_next_border = get_cw_first_vertex(hvr, vecVCandidates);
					if (Vertex_handle() != vh_next_border) {
						vecVHandles.push_back(vh_next_border);
						VertexIndexsSet.insert(vh_next_border);
						hvr = vh_next_border;
						bRightUpdate = true;
					} else {
						bRightUpdate = false;
					}
				}
			}
		}

		//遍历到下切线，去重后，将顶点加入到三角化中(zig_zag)
		void get_bottom_tangent_vetex(SSectionData& leftGraph, SSectionData& rightGraph, std::vector<Vertex_handle>& vecVHandles)
		{
			Vertex_handle hvl = leftGraph.vertex_lr[1];
			Vertex_handle hvr = rightGraph.vertex_lr[0];
			if (hvl == Vertex_handle() || hvr == Vertex_handle())
				return;

			std::set<Vertex_handle> VertexIndexsSet;
			VertexIndexsSet.insert(hvl);
			VertexIndexsSet.insert(hvr);
			vecVHandles.push_back(hvl);
			vecVHandles.push_back(hvr);
			std::vector<Vertex_handle> vecVCandidates;
			bool bLeftUpdate = true;
			bool bRightUpdate = true;
			while (bLeftUpdate || bRightUpdate)
			{
				if (bLeftUpdate || bRightUpdate)
				{
					if (leftGraph.m_dt.is_infinite(hvl))
						continue;

					VCirculator ccl = leftGraph.m_dt.incident_vertices(hvl);
					VCirculator ccl_end(ccl);
					vecVCandidates.clear();
					do {
						if (!leftGraph.m_dt.is_infinite(ccl->handle()) && is_border(leftGraph, ccl->handle())) {
							if (VertexIndexsSet.find(ccl->handle()) == VertexIndexsSet.end()) {
								if (CGAL::LEFT_TURN != CGAL::orientation(hvl->point(), hvr->point(), ccl->handle()->point())) {
									vecVCandidates.push_back(ccl->handle());
								}
							}
						}
					} while (++ccl != ccl_end);

					Vertex_handle vh_next_border = get_cw_first_vertex(hvl, vecVCandidates);
					if (Vertex_handle() != vh_next_border) {
						vecVHandles.push_back(vh_next_border);
						VertexIndexsSet.insert(vh_next_border);
						hvl = vh_next_border;
						bLeftUpdate = true;
					} else  {
						bLeftUpdate = false;
					}
				}

				if (bLeftUpdate || bRightUpdate)
				{
					if (rightGraph.m_dt.is_infinite(hvr))
						continue;

					vecVCandidates.clear();
					VCirculator ccr = rightGraph.m_dt.incident_vertices(hvr);
					VCirculator ccr_end(ccr);
					do {
						if (!rightGraph.m_dt.is_infinite(ccr->handle()) && is_border(rightGraph, ccr->handle())) {
							if (VertexIndexsSet.find(ccr->handle()) == VertexIndexsSet.end()) {
								if (CGAL::LEFT_TURN != CGAL::orientation(hvl->point(), hvr->point(), ccr->handle()->point())) {
									vecVCandidates.push_back(ccr->handle());
								}
							}
						}
					} while (++ccr != ccr_end);

					Vertex_handle vh_next_border = get_ccw_first_vertex(hvr, vecVCandidates);
					if (Vertex_handle() != vh_next_border) {
						vecVHandles.push_back(vh_next_border);
						VertexIndexsSet.insert(vh_next_border);
						hvr = vh_next_border;
						bRightUpdate = true;
					} else {
						bRightUpdate = false;
					}
				}
			}
		}

		void boundingBox(double x, double y) 
		{
			if (x < m_dLeft)
				m_dLeft = x;
			if (x > m_dRight)
				m_dRight = x;

			if (y < m_dBottom)
				m_dBottom = y;
			if (y > m_dTop)
				m_dTop = y;
		}

		bool check_circle_intersect(Face_handle handle, Segment_2& seg, Orientation ori)
		{
			Point2D &p1 = handle->vertex(0)->point();
			Point2D &p2 = handle->vertex(1)->point();
			Point2D &p3 = handle->vertex(2)->point();
			Circle_2 cir(p1, p2, p3);

			if (CGAL::do_intersect(cir, seg))
			{
				if (ori == Near_Left) {
					m_setRiskHandles_left.insert(handle);
				} else {
					m_setRiskHandles_right.insert(handle);
				}
				
				m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(handle->vertex(0)->info(), handle->vertex(0)));
				m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(handle->vertex(1)->info(), handle->vertex(1)));
				m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(handle->vertex(2)->info(), handle->vertex(2)));
				Face face(handle->vertex(0)->info(), handle->vertex(1)->info(), handle->vertex(2)->info());
				m_facesmap.insert(std::pair<Face, bool>(face, false));
				return true;
			}

			update_right_most_vertex(handle);
			update_left_most_vertex(handle);

			m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(handle->vertex(0)->info(), handle->vertex(0)));
			m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(handle->vertex(1)->info(), handle->vertex(1)));
			m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(handle->vertex(2)->info(), handle->vertex(2)));

			Face face(handle->vertex(0)->info(), handle->vertex(1)->info(), handle->vertex(2)->info());
			m_facesmap.insert(std::pair<Face, bool>(face, true));
			return false;
		}
		
		Orientation orientation(SSectionData& other)
		{
			if (m_iSampleBeginX < 0 || m_iSampleEndX < 0)
				return Not_Near;

			if (m_iSampleBeginX >= m_iSampleEndX || other.m_iSampleBeginX >= other.m_iSampleEndX)
				return Not_Near;
			
			if (m_iSampleBeginX == other.m_iSampleEndX)
				return Near_Right;
			if (m_iSampleEndX == other.m_iSampleBeginX)
				return Near_Left;
		}

		void insert_delaunay(Face_handle& handle)
		{
			Vertex_handle h1 = m_dt.insert(handle->vertex(0)->point());
			Vertex_handle h2 = m_dt.insert(handle->vertex(1)->point());
			Vertex_handle h3 = m_dt.insert(handle->vertex(2)->point());
			h1->info() = handle->vertex(0)->info();
			h2->info() = handle->vertex(1)->info();
			h3->info() = handle->vertex(2)->info();
		}

		void insert_risk_delaunay(SSectionData& other, Orientation ori)
		{
			if (ori == Near_Left) {
				for (auto handle : other.m_setRiskHandles_left) {
					insert_delaunay(handle);
				}
			} else {
				for (auto handle : other.m_setRiskHandles_right) {
					insert_delaunay(handle);
				}
			}
		}

		void insert_tangent_delaunay(SSectionData& graph1, SSectionData& graph2) 
		{
			std::vector<Vertex_handle> vecVHandles;
			get_up_tangent_vetex(graph1, graph2, vecVHandles);
			get_bottom_tangent_vetex(graph1, graph2, vecVHandles);
		
			for (auto& vh : vecVHandles)
			{
				auto ittr = m_vertexHashtable.find(vh->info());
				if (ittr == m_vertexHashtable.end())
				{
					Vertex_handle h1 = m_dt.insert(vh->point());
					h1->info() = vh->info();
					m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(h1->info(), h1));
				}
			}
		}

		bool check_face_intersect(Face_handle& handle, Segment_2 seg)
		{
			Point2D p1 = handle->vertex(0)->point();
			Point2D p2 = handle->vertex(1)->point();
			Point2D p3 = handle->vertex(2)->point();

			Segment_2 seg1(p1, p2);
			Segment_2 seg2(p1, p3);
			Segment_2 seg3(p2, p3);

			if (CGAL::do_intersect(seg, seg1)
				|| CGAL::do_intersect(seg, seg2)
				|| CGAL::do_intersect(seg, seg3))
			{
				Face face(handle->vertex(0)->info(), handle->vertex(1)->info(), handle->vertex(2)->info());
				auto ittr = m_facesmap.find(face);
				if (ittr == m_facesmap.end()) {
					m_facesmap.insert(std::pair<Face, bool>(face, true));
				} else {
					ittr->second = true;
				}
				return true;
			}
			return false;
		}

		bool check_face_duplicate(Face_handle& handle)
		{
			auto ittr = m_facesmap.find(Face(handle->vertex(0)->info(), handle->vertex(1)->info(), handle->vertex(2)->info()));
			if (ittr != m_facesmap.end())
			{
				ittr->second = true;
				return true;
			}
			return false;
		}

		int merge(std::vector<Face_handle> &vecMergeResult, SSectionData& leftSection, SSectionData& rightSection)
		{
			std::vector<Face_handle> vecFaceCandidates;
			std::vector<unsigned int> vecDuliInfos;
			int iCount = 0;
			while (!m_deqBarrierFaces.empty())
			{
				Face_handle cur_face = m_deqBarrierFaces.front();
				m_deqBarrierFaces.pop_front();
			
				//bfs由已确定的面为基础，确定四周3个邻居面
				for (int i = 0; i < 3; ++i)
				{
					Face_handle nfh = cur_face->neighbor(i);

					if (m_dt.is_infinite(nfh->vertex(0)) || m_dt.is_infinite(nfh->vertex(1)) || m_dt.is_infinite(nfh->vertex(2))) 
						continue;
					if (leftSection.m_dt.is_infinite(nfh->vertex(0)) || leftSection.m_dt.is_infinite(nfh->vertex(1)) || leftSection.m_dt.is_infinite(nfh->vertex(2)))
						continue;
					if (rightSection.m_dt.is_infinite(nfh->vertex(0)) || rightSection.m_dt.is_infinite(nfh->vertex(1)) || rightSection.m_dt.is_infinite(nfh->vertex(2)))
						continue;

					//已经确定的就跳过
					auto ittr = leftSection.m_facesmap.find(Face(nfh->vertex(0)->info(), nfh->vertex(1)->info(), nfh->vertex(2)->info()));
					if (ittr != leftSection.m_facesmap.end() && ittr->second)
						continue;
					ittr = rightSection.m_facesmap.find(Face(nfh->vertex(0)->info(), nfh->vertex(1)->info(), nfh->vertex(2)->info()));
					if (ittr != rightSection.m_facesmap.end() && ittr->second)
						continue;

					Vertex_handle vh0 = cur_face->vertex(i);
					Vertex_handle vh1 = cur_face->vertex(CDT2D::ccw(i));
					Vertex_handle vh2 = cur_face->vertex(CDT2D::ccw(CDT2D::ccw(i)));

					//走到对边边上的其中一个点上面(vh1),根据顶点所处的分区，捞出该顶点周围所有的面
					unsigned int vi = vh1->info();
					bool bLeft = true;
					vecFaceCandidates.clear();
					auto ittr_vertex = leftSection.m_vertexHashtable.find(vi);
					if (ittr_vertex != leftSection.m_vertexHashtable.end()){
						bLeft = true;
						get_faces_around_vertex(leftSection, ittr_vertex->second, vecFaceCandidates);
					} else {
						auto ittr_vertex = rightSection.m_vertexHashtable.find(vi);
						if (ittr_vertex != rightSection.m_vertexHashtable.end()) {
							bLeft = false;
							get_faces_around_vertex(rightSection, ittr_vertex->second, vecFaceCandidates);
						} else {
							continue;
						}
					}

					//根据顶点的周围的面，判断哪个才是真正的邻居
					for (auto &facehandle : vecFaceCandidates)
					{
						iCount = 0;
						int other_index = -1;
						vecDuliInfos.clear();
						vecDuliInfos.push_back(facehandle->vertex(0)->info());
						vecDuliInfos.push_back(facehandle->vertex(1)->info());
						vecDuliInfos.push_back(facehandle->vertex(2)->info());
						for (auto& index : vecDuliInfos) 
						{
							if (index == vh1->info() || index == vh2->info()) {
								++iCount;
							} else {
								other_index = index;
							}
						}

						if (iCount != 2)
							continue;
						if (other_index == vh0->info())
							continue;

						//已经有两个点在 重复面和跨边面组成的连续区域 上，如果剩下的最后一个点还在这个区域上，那这个面就不是我们要找点的面）
						//这个能够工作，还是不是太能解释，猜测：这个连续区域是个凸包，如果用bfs往外一层，还是个凸包
						//要不然没有办法解释如果这个连续区域的外围是凹的情况（好像没有办法正常工作）
						auto ittr = m_vertexHashtable.find(other_index);
						if (ittr != m_vertexHashtable.end())
							continue;

						bool bNew = false;
						if (bLeft){
							auto ittr = leftSection.m_facesmap.find(Face(facehandle->vertex(0)->info(), facehandle->vertex(1)->info(), facehandle->vertex(2)->info()));
							if (ittr != leftSection.m_facesmap.end() && !ittr->second) {
								ittr->second = true;
								bNew = true;
							}
						} else {
							auto ittr = rightSection.m_facesmap.find(Face(facehandle->vertex(0)->info(), facehandle->vertex(1)->info(), facehandle->vertex(2)->info()));
							if (ittr != rightSection.m_facesmap.end() && !ittr->second) {
								ittr->second = true;
								bNew = true;
							}
						}

						if (bNew) {
							vecMergeResult.push_back(facehandle);
							m_deqBarrierFaces.push_back(facehandle);
							m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(facehandle->vertex(0)->info(), facehandle->vertex(0)));
							m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(facehandle->vertex(1)->info(), facehandle->vertex(1)));
							m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(facehandle->vertex(2)->info(), facehandle->vertex(2)));
						}
					}
				}
			}
			return vecMergeResult.size();
		}

		int get_infinite_vertex_index(SSectionData& section, Circulator cc, std::vector<int> &vecIndexs)
		{
			if (section.m_dt.is_infinite(cc->vertex(0)))
				vecIndexs.push_back(0);

			if (section.m_dt.is_infinite(cc->vertex(1)))
				vecIndexs.push_back(1);

			if (section.m_dt.is_infinite(cc->vertex(2)))
				vecIndexs.push_back(2);

			return vecIndexs.size();
		}

		int get_faces_around_vertex(SSectionData& section, Vertex_handle& vhandle, std::vector<Face_handle> &vecFaceCandidates)
		{
			TDS& tds = section.m_dt.tds();
			Circulator cc = tds.incident_faces(vhandle);
			Circulator cc_done(cc);
			Face_handle fh;

			std::vector<int> infiniteIndexs;
			do{
				if (get_infinite_vertex_index(section, cc, infiniteIndexs) > 0)
					continue;

				Face face(cc->vertex(0)->info(), cc->vertex(1)->info(), cc->vertex(2)->info());
				auto ittr = section.m_facesmap.find(face);
				if (ittr != section.m_facesmap.end() && !ittr->second)
					vecFaceCandidates.push_back(cc);
			} while (++cc != cc_done);

			return vecFaceCandidates.size();
		}
	};

	struct merging_task : public task
	{
		merging_task(Delaunay_Mutithread_CGAL *obj, SSectionData *data1, SSectionData *data2) 
		{
			m_obj = obj;
			m_data1 = data1;
			m_data2 = data2;
		}

		void operator()()
		{
			if (m_obj)
				m_obj->merge(m_data1, m_data2);
		}
		Delaunay_Mutithread_CGAL *m_obj;
		SSectionData* m_data1;
		SSectionData* m_data2;
	};

	struct delaunay_task : public task
	{
		delaunay_task(Delaunay_Mutithread_CGAL *obj, SSectionData *data)
		{
			m_obj = obj;
			m_data = data;
		}

		void operator()()
		{
			if (m_obj)
				m_obj->delaunay_CGAL(m_data);
		}
		Delaunay_Mutithread_CGAL *m_obj;
		SSectionData* m_data;
	};

private:
	std::recursive_mutex m_mutex;
	CThreadPool *m_pThreadPool;

	osg::ref_ptr<osg::Vec3Array> m_vecPointsRef;
	osg::ref_ptr<osg::DrawElementsUInt> m_PrimitveSetRef;

	std::vector<task*> m_vecTasks;
	std::vector<SSectionData> m_vecSectionData;  //

protected:
	SSectionData::Orientation orientation(SSectionData &section1, SSectionData &section2);

	void checkMergingAndFinish();
	int hasMerging(std::set<std::pair<int, int>>& setMergingIndexs);
	bool checkAllFinish();

	void merge(SSectionData* section1, SSectionData* section2);
	void delaunay_CGAL(SSectionData *sectionData);
	void sample(float fSampleRatio, int iThreadNum, osg::Vec3Array::iterator begin, osg::Vec3Array::iterator end);
};
