#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/HalfedgeDS_default.h>

#include <osg/Geometry>

#include <limits>

//Ŀ�ģ�
//1. ��ϰʹ��CGAL
//2. ʵ�ַ����㷨
//3. ʵ�ֶ��̣߳����ǳ������ݼ�
//4. ����Ĳ������ݼ��ĵ��롢��ʾ���ӻ��������������osg���жԽ�
//5. ����Ĳ�ͬ�㷨��ȷ�ԡ�Ч�ʱȽϷ�������ζԱȲ�ͬ�㷨���Ա�ÿ���������Ƿ�һ�£����߼򵥿������������Ƿ��غϣ�
//6. ����Ľ���CGAL�Ľӿڲ��ԣ�������Ϥ������ݽṹ

//��ʼǰ���Ա�һ��CDT��CGAL�����⣺
//1. CGAL�����ݽṹ����Compact_container<base>���ο�STL���ڴ�ķ��������Ż�����������ݿ�������ÿ��node�ڵ㶼Ҫ�洢
// ��һ���ڵ�ָ�룬ռ����8���ֽ�(x64)�� CGAL�п���ʹ��������Union�ķ��������ָ������ݣ����ڴ��������ͷŹ����У��л��Ṳ��ͬһ���ڴ��ַ���µģ�Դ��û�п�̫����
// CDT�Ļ�ֱ��ʹ��std::vector, std::map���൱�ڵ����洢��ָ�������
//2. ׼ȷ�����棬���ڶԱ���֤��
//3. �����ⶼʹ�����������㷨��������point location������㷨�����ϲ�ͬ
//(1)CGAL�����˱��������ݽṹ���ص㣬��һ���������ͨ��zig_zag�ҵ�Ŀ���棬
// ����ÿ�ο��Դ��ĸ���λ��ȷ��һ�����о�Ӧ����log(n)��Ч��Ӧ�ò���̫�ͣ�����point location��֮��Ĳ�����������ǻ���
// ͬһ��������ݽṹ��û��CDT�е�kdtreeά���������������������kd���������������Ŀ�����
//(2)CDTʹ����kd������point location���̣���Ȼ��ÿ�ε����������У�ֱ�ӽ��������ݵ����ô��룬������һ���ֵĿ���������kd������Ҳ�洢���ݣ�
//�ò��ִ洢�����ڶ������ʱ�İ�����ݽṹ֮�⣻�о�kd�����ԽϿ�Ķ�λ������㸽�����򣬵������Ҫ��ά��kd��Ϊ���ۣ�
//���ҵ��������֮�󣬻���û�а취100%ȷ������Ŀ�����򣬻�Ҫ�����򸽽�����һ�������жϣ����������Ч�ʻ�Ҫ����ȷ�ϣ�
//�о�kd����ĳЩ��������ݵ㼯������£�����־ֲ����е����⣨����㼯����һ��С�����ڣ�
//4. CGAL�����˱Ƚ϶�ķ�װ��ʹ������Ƚ��ѣ���ͬ�㷨�������˲�ͬ�����ݽṹ��CDTʹ�ñȽϼ�
//5. ��Ŀǰ��ʹ������������ߺ���Ҫ����һ�����ݵĿ�����ȥ

//ΪʲôҪ��devide and conquer��
//1. �Լ��ĸо�����Ϊ�ṹ����Ψһ�ԣ����һ��С�����ڣ�����100���˴˿�����͹���ڣ��Ѿ�����delaunay���ǻ��֣�
//��ô�����͹����Χ����㻰�о����Խ�flip edge�Ĵ���������ͣ�����ʹ���������㷨���о��ƻ��������ľֲ����ԣ�
//���ܵ��ºܶ಻��Ҫ��flip edge(�����ڲ���ǰ�������Ѿ�����Ԥ�У�������ĳ�����ƿռ�������)��
//(1)���ڴ����ݼ���һ�����з�������ÿ�������Ӽ��ֲ����������루�����ö��̣߳���Ȼ���֮����Ӽ��ٽ���һ��devide and conquer
//�ں��Ӽ���������Ӧ�����������������������ף�
//2. devide and conquer�о�������Ȼ�оֲ������ƣ��������ڴ��ʹ��������Ͳ����ı���Ѷ�������񶼱Ƚ�ռ����

//���ݽṹ��ѡ��
//1. �ų���Triangulation_data_structure_2
//CGAL��Triangulation_data_structure_2��������ǻ����������ݽṹ��ÿ�β���һ���㣬������е����������Σ�������Ҫ�ҵ���������ߵ�λ��
//�������㷨�ǰ�һЩ�����Ƭ����������������������У�����һЩ�����ĵ㡢�ߡ���
//����������ݽṹ�о��������䣬��������ʵ�ַ������Խ����������  ����ֲ���_tds����������+���ν�ϣ�������ǲ���
//�����ˣ�����
//2. ���ǿ���ʹ�ó���Ŀ⣬ͨ�õ����ݽṹ��BGL�м������������ݽṹ��
//3. ��򵥵��ڽӱ�Ӧ�ò��У���Ϊ�򵥵�CCW��CW���������鷳
//���ۣ���halfedge���԰ɣ�CGAL��һ�����ݽṹ��HalfedgeDS
//�ο�1��https://doc.cgal.org/latest/HalfedgeDS/index.html#Chapter_Halfedge_Data_Structures
//�ο�2��https://doc.cgal.org/latest/BGL/group__PkgBGLTraits.html

//surface_mesh



//TODO: devide_and_conqure�ı߽���������(�˻����)
//1. ���㹲��
//2. �ĵ㹲Բ
//3. ���Ƕ��̰߳汾�� ��������Ӽ����������У���ʱû����������

//TODO: devide_and_conqure�㷨���еĿ�Դ��Ŀ

//���϶���Ŀǰ���������⣬��������������������ҪĿ�Ļ������Ȱ�CGAL������
//TODO: ��زο����ϻ���
// [0]https://www.doc88.com/p-694154354521.html?r=1
// [1]
// [2]
// 
//author: zhanghcuanpan
//date: 2023.01.14

template <typename T>
class MyTrait
{
    typedef T Point_2;
};

template <>
class MyTrait<osg::Vec3Array>
{
    typedef osg::Vec3Array::iterator Point_2;
};

typedef MyTrait<osg::Vec3Array> OSGTrait;

#define INFINITE  std::numeric_limits<unsigned int>::max()

//ģ���̸㲻����ֱ��д����
template <class traits = OSGTrait>
class Delaunay_devide_and_conqure
{
public:
    typedef typename traits::Point_2                                         Point;
    typedef typename CGAL::Exact_predicates_inexact_constructions_kernel     K;
    typedef typename CGAL::HalfedgeDS_default<traits>                        HDS;
    typedef typename CGAL::HalfedgeDS_decorator<HDS>                         Decorator;
    typedef typename HDS::Halfedge                                           Halfedge;
    typedef typename HDS::Halfedge_handle                                    Halfedge_handle;
    typedef typename HDS::Vertex                                             Vertex;
    typedef typename HDS::Vertex_handle                                      Vertex_handle;
    typedef typename HDS::Face                                               Face;
    typedef typename HDS::Face_handle                                        Face_handle;

    typedef typename CGAL::Projection_traits_xy_3<K>                         Gt;    //����������м���

    typedef typename traits::Point_2                             OutContainerIttr;
    typedef typename std::pair<Vertex_handle, Vertex_handle>     LimitHandlesPair;
    typedef typename std::pair<Vertex_handle, Vertex_handle>     TangentLine;
    typedef typename std::pair<TangentLine, TangentLine>         TopAndBottomTangent;

    Delaunay_devide_and_conqure() {}

    Delaunay_devide_and_conqure(OutContainerIttr begin, OutContainerIttr end)
    {
        _pHds = new HDS;
        _pDecorator = new Decorator(*_pHds);

        unsigned int size = std::distance(begin, end);
        if (size > 0)
            insert(begin, end);

        //v, h, f
        _pHds->reserve(size, size, size);
    }

    ~Delaunay_devide_and_conqure()
    {
        delete _pDecorator;
        delete _pHds;
    }

    class LexiSortOpr
    {
        bool operator()(const OutContainerIttr& ittr1, const OutContainerIttr& ittr2)
        {
            if (ittr1->x() < ittr2->x()) {
                return true;
            }
            else {
                if (ittr1->y() < ittr2->y())
                    return true;
                return false;
            }
        }
    };

private:
    std::set<OutContainerIttr, LexiSortOpr> m_setPointsIttrs;  //�����ݵ��������lexicographically ascending��
    HDS* _pHds;
    Decorator* _pDecorator;
    Gt _gt;

public:
    void insert(OutContainerIttr begin, OutContainerIttr end)
    {
        //�Բ�������ݵ㼯�������򣬷�������devide_and_conqure
        for (OutContainerIttr ittr = begin; ittr != end; ++ittr)
        {
            //����ܹ�˳��ȥ�ؾͺ���
            m_setPointsIttrs.insert(ittr);
        }
    }

    bool is_valid_handle(Vertex_handle vh) { return true; }
    bool is_valid_handle(Face_handle fh) { return true; }
    bool is_valid_handle(Halfedge_handle hh) { return true; }
    bool is_valid_TangentLine(TangentLine line)
    {
        if (is_valid_handle(line.first) && is_valid_handle(line.second))
            return true;
        return false;
    }

    bool is_same_tangent(TangentLine line1, TangentLine line2)
    {
        if (!is_valid_TangentLine(line1) || !is_valid_TangentLine(line2))
            return false;

        if (line1.first == line2.first && line1.second == line2.second)
            return true;

        return false;
    }

    //tool_2: Ѱ������͹������������
    TopAndBottomTangent findTangent(LimitHandlesPair leftHull, LimitHandlesPair rightHull)
    {
    }

    //���������Σ���ֻ��һ������ʱ��������������ľ����Ч��������Ȼ�������β��뵽���ݽṹ�У�һ����ʱ��ͬ��
    LimitHandlesPair insert_hds(OutContainerIttr pr1, OutContainerIttr pr2, OutContainerIttr pr3)
    {
        std::vector<OutContainerIttr> vecValidIttrs;
        if (pr1 != OutContainerIttr())
            vecValidIttrs.push_back(pr1);
        if (pr2 != OutContainerIttr())
            vecValidIttrs.push_back(pr2);
        if (pr3 != OutContainerIttr())
            vecValidIttrs.push_back(pr3);

        Vertex_handle h1;
        Vertex_handle h2;
        Vertex_handle h3;
        Vertex_handle LeftMostVertexHandle;
        Vertex_handle RightMostVertexHandle;

        //������ʱ��
        int size = vecValidIttrs.size();
        if (size == 3)
        {
            Gt::Point_2 p1(vecValidIttrs[0]->x(), vecValidIttrs[0]->y());
            Gt::Point_2 p2(vecValidIttrs[1]->x(), vecValidIttrs[1]->y());
            Gt::Point_2 p3(vecValidIttrs[2]->x(), vecValidIttrs[2]->y());
            CGAL::Orientation orientation = _gt.orientation_2_object(p1, p2, p3);
            if (orientation == CGAL::CLOCKWISE)
            {
                //p1->p3->p2
                Vertex v1(vecValidIttrs[0]);
                Vertex v2(vecValidIttrs[2]);
                Vertex v3(vecValidIttrs[1]);
                //vh1 = _pDecorator->vertices_push_back(v1);
                //vh2 = _pDecorator->vertices_push_back(v2);
                //vh3 = _pDecorator->vertices_push_back(v3);

                Halfedge_handle hf1 = _pHds->edges_push_back(Halfedge());
                Halfedge_handle hf2 = _pHds->edges_push_back(Halfedge());
                Halfedge_handle hf3 = _pHds->edges_push_back(Halfedge());
                //hf1->HBase::set_prev(hf2);
                //_pDecorator->set_face();
                //_pDecorator->set_vertex();
                //_pDecorator->set_prev();


                Face face(edge1);
                Face_handle fh = _pDecorator->faces_push_back(face);

                LeftMostVertexHandle = (vecValidIttrs[0].x() < vecValidIttrs[1].x()) ? vh1 : vh2;
                LeftMostVertexHandle = (LeftMostVertexHandle->point()->x() < vecValidIttrs[2].x()) ? LeftMostVertexHandle : vh3;

                RightMostVertexHandle = (vecValidIttrs[0].y() < vecValidIttrs[1].y()) ? vh2 : vh1;
                RightMostVertexHandle = (RightMostVertexHandle->point()->y() < vecValidIttrs[2].y()) ? vh3 : RightMostVertexHandle;
            }
            else
            {
                //p1->p2->p3
                Vertex v1(vecValidIttrs[0]);
                Vertex v2(vecValidIttrs[1]);
                Vertex v3(vecValidIttrs[2]);
                vh1 = _pDecorator->vertices_push_back(v1);
                vh2 = _pDecorator->vertices_push_back(v2);
                vh3 = _pDecorator->vertices_push_back(v3);

                Halfedge edge1;
                Halfedge edge2;
                Halfedge edge3;

                Face face;

                LeftMostVertexHandle = (vecValidIttrs[0].x() < vecValidIttrs[1].x()) ? vh1 : vh2;
                LeftMostVertexHandle = (LeftMostVertexHandle->point()->x() < vecValidIttrs[2].x()) ? LeftMostVertexHandle : vh3;

                RightMostVertexHandle = (vecValidIttrs[0].y() < vecValidIttrs[1].y()) ? vh2 : vh1;
                RightMostVertexHandle = (RightMostVertexHandle->point()->y() < vecValidIttrs[2].y()) ? vh3 : RightMostVertexHandle;
            }
        }
        else if (size == 2) 
        {
        
        }
        else if (size == 1) 
        {
        
        }
        else 
        {
            return LimitHandlesPair();
        }
        return LimitHandlesPair(LeftMostVertexHandle, RightMostVertexHandle);
    }

    //devide_and_conqure
    //�����������߼���λ�õľ��
    LimitHandlesPair devide_and_conqure(OutContainerIttr begin, OutContainerIttr end)
    {
        unsigned int size = std::distance(begin, end);
        if (size == 1 || size <= 0)
        {
            OutContainerIttr pr1 = OutContainerIttr();
            OutContainerIttr pr2 = OutContainerIttr();
            OutContainerIttr pr3 = OutContainerIttr();
            return insert_hds(p1, p2, p3);
        }
        else if (size == 2)
        {
            //step1 : ����_tds������������
            OutContainerIttr pr1 = begin;
            OutContainerIttr pr2 = begin + 1;
            OutContainerIttr pr3 = OutContainerIttr();
            return insert_hds(p1, p2, p3);
        }
        else if (size == 3)
        {
            //����_tds������������
            OutContainerIttr pr1 = begin;
            OutContainerIttr pr2 = begin + 1;
            OutContainerIttr pr3 = begin + 2;
            return insert_hds(p1, p2, p3);
        }

        //dv_callback

        //conditon: size > 3
        int mid = size/2;
        OutContainerIttr begin_1 = begin;
        OutContainerIttr end_1 = begin + mid;
        OutContainerIttr begin_2 = end1 + 1;
        OutContainerIttr end_2 = end;

        LimitHandlesPair limitHandles_L = devide_and_conqure(begin_1, end_1);
        LimitHandlesPair limitHandles_R = devide_and_conqure(begin_2, end_2);
        return merge(limitHandles_L, limitHandles_R);
    }

    //merge function for devide_and_conqure

    LimitHandlesPair merge(LimitHandlesPair limitLeft, LimitHandlesPair limitRight)
    {
        //step1: Ѱ������͹�����������ߣ�͹��������1�����������Ρ�Ҳ�����ǵ�����һ���㣬Ҳ�����ǵ�����һ���ߣ�Ҳ������delaunay����   
        TopAndBottomTangent tangents = findTangent(LimitHandlesPair limitLeft, LimitHandlesPair limitRight)
        bool bLeft = is_valid_TangentLine(tangents.first);
        bool bRight = is_valid_TangentLine(tangents.first);
        if (!bLeft || bRight)
        {
            if (bLeft) return limitLeft;  
            if (bRight) return limitRight;
            return LimitHandlesPair();
        }

        //step2: ��ʼ�ݽ������ߣ�������delaunay���ǻ���
        //TODO: merge_callback
        while(!is_same_tangent(tangents.first, tangents.second))
        {
            //step2_1���������ߵġ��󡿶��㣨λ����͹����Ϊ���㣬�ڡ��ҡ�͹���У�Ѱ�������ߵġ���ѡ�㡿
            //Ѱ�ҹ����У���Ҫʱ�Ͽ����ҡ�͹����ĳЩ�ߣ���ߺ��治�����������ˣ�


            //step2_2���������ߵġ��ҡ����㣨λ����͹����Ϊ���㣬�ڡ���͹���У�Ѱ�������ߵġ���ѡ�㡿
            //Ѱ�ҹ����У���Ҫʱ�Ͽ�����͹����ĳЩ�ߣ���ߺ��治�����������ˣ�


            //step2_3: 
            //�����ϻ���ͺ�ѡ���У�ȷ��������, ��_tds�в���ñ�
            //ѭ���ݽ��������µ�������

            
        }//step2_4: continue the while circle

        Vertex_handle LeftMostVertexHandle = limitLeft.first;
        Vertex_handle RightMostVertexHandle = limitRight.second;
        return LimitHandlesPair(LeftMostVertexHandle, RightMostVertexHandle);
    }

private:


};

