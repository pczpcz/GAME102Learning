#pragma once
//先自己写一个数据结构，后面看看把BGL或者CGAL用起来，或者用其他的， 如OpenMesh（模板使用还是不太熟悉，用起来比较费时间）

#include <vector>
#include <limits>
#include <algorithm>

struct Face;
struct HalfEdge;
struct Vertex;

struct Face
{
	Face():pBelongToEdge(nullptr){}

	HalfEdge *pBelongToEdge;
};

struct HalfEdge
{
	HalfEdge():index(-1),pPrevEdge(nullptr),pNextEdge(nullptr),pOwnFace(nullptr){}
	HalfEdge(int iBegin, int iEnd, HalfEdge* pEdge1 = nullptr, HalfEdge* pEdge2 = nullptr,Face* pFace = nullptr)
		:index(iBegin),pPrevEdge(pEdge1),pNextEdge(pEdge2),pOwnFace(pFace) {}

	int index;

	HalfEdge *pPrevEdge;
	HalfEdge *pNextEdge;
	HalfEdge* pOpsiteEdge;

	Face *pOwnFace;
};

struct Vertex
{
	Vertex():px(0),pz(0),pOutEdge(nullptr){}
	Vertex(int index, HalfEdge *pEdge = nullptr):px(-1),py(-1),pz(-1),pOutEdge(pEdge) {}

	float px;
	float py;
	float pz;
	HalfEdge *pOutEdge;
};

class Graph
{
};

//仅使用于2D的情况
class HalfEdgeDataStruct : public Graph
{
public:
	Vertex* vertex(int index) 
	{
		if (index >= 0 && index < m_vecVertexs.size()) 
		{
			return &m_vecVertexs[index];
		}
		return nullptr;
	}
	HalfEdge* halfEdge(int iBegin, int iEnd) 
	{
		Vertex* pVertex = vertex(iBegin);
		if (!pVertex || pVertex->pOutEdge) return nullptr;

		if (pVertex->pOutEdge->pNextEdge || pVertex->pOutEdge->pNextEdge->index == iEnd)
		{
			return pVertex->pOutEdge;
		}
		return nullptr;
	}
	Face* face(int iBegin, int iEnd) 
	{
		HalfEdge* pEdge = halfEdge(iBegin, iEnd);
		if (!pEdge) return nullptr;

		return pEdge->pOwnFace;
	}
	int AdjanceFace(int iBegin, int iEnd, std::vector<Face *> &vecFaces)
	{
		HalfEdge* pEdge = halfEdge(iBegin, iEnd);
		if (!pEdge) return -1;

		getAdjanceFaces(pEdge, vecFaces);
		return vecFaces.size();
	}

	bool addEdge(int iBegin, int iEnd) 
	{
		//写不下去了。。。。。。
		//写不下去了。。。。。。

		//尴尬了，还是直接用库吧

	}

	bool removeEdge(int iBegin, int iEnd)
	{
		//写不下去了。。。。。。
		//写不下去了。。。。。。

		//尴尬了，还是直接用库吧
	}

	bool exchangeEdge(HalfEdge* pEdge1, HalfEdge* pEdge2)
	{
		//写不下去了。。。。。。
		//写不下去了。。。。。。

		//尴尬了，还是直接用库吧
	}

	int getFaceVertexs(Face* face, std::vector<int> &vecVertex) 
	{
		HalfEdge* pEdge = face->pBelongToEdge;
		if (!pEdge) 
		{
			return -1;
		}

		HalfEdge* pNextEdge = pEdge;
		do 
		{
			vecVertex.push_back(pNextEdge->index);
			pNextEdge = pNextEdge->pNextEdge;
		} while (pEdge != pNextEdge && !pNextEdge);

		return vecVertex.size();
	}
	int getBoundEdges(Face* face, std::vector<HalfEdge *> &vecEdges) 
	{
		if (!face) return -1;
		return getFaceEdges(face->pBelongToEdge, vecEdges);
	}
	int getAdjanceFaces(HalfEdge* pEdge, std::vector<Face *> &vecFaces)
	{
		if (!pEdge || !pEdge->pOwnFace) return -1;

		vecFaces.push_back(pEdge->pOwnFace);
		if (pEdge->pOpsiteEdge && pEdge->pOpsiteEdge->pOwnFace)
			vecFaces.push_back(pEdge->pOpsiteEdge->pOwnFace);
		return vecFaces.size();
	}
	int getFaceEdges(HalfEdge* edge, std::vector<HalfEdge *> &vecEdges)
	{
		if (!edge) return -1;
		vecEdges.push_back(edge);
		HalfEdge* pEdge = edge->pNextEdge;
		while(!pEdge && pEdge->index>=0 && pEdge->index != edge->index)
		{
			vecEdges.push_back(pEdge);
			pEdge = pEdge->pNextEdge;
		}

		return vecEdges.size();
	}
	int getOutEdges(Vertex* vertex, std::vector<HalfEdge *> &vecEdges)
	{
		if (!vertex && !vertex->pOutEdge) return -1;
	
		HalfEdge* pOutEdge = vertex->pOutEdge;
		HalfEdge* pOpsiteEdge = pOutEdge->pOpsiteEdge;
		vecEdges.push_back(pOutEdge);
		while (!pOpsiteEdge && !pOpsiteEdge->pNextEdge && pOpsiteEdge->pNextEdge != vertex->pOutEdge) 
		{
			pOutEdge = pOpsiteEdge->pNextEdge;
			pOpsiteEdge = pOutEdge->pOpsiteEdge;
			vecEdges.push_back(pOutEdge);
		}
		return vecEdges.size();
	}
	int getAdjanceFaces(Vertex* vertex, std::vector<Face*> &vecFaces)
	{
		if (!vertex) return -1;

		std::vector<HalfEdge*> vecEdges;
		getOutEdges(vertex, vecEdges);
		for (int i = 0; i < vecEdges.size(); ++i)
		{
			if (vecEdges[i]->pOwnFace)
				vecFaces.push_back(vecEdges[i]->pOwnFace);
		}
		return vecFaces.size();
	}

private:
	std::vector<Vertex> m_vecVertexs;
	std::vector<HalfEdge> m_vecEdges;
	std::vector<Face> m_vecFaces;
};


