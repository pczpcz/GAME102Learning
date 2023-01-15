#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/HalfedgeDS_default.h>

#include <osg/Geometry>

#include <limits>

//目的：
//1. 练习使用CGAL
//2. 实现分治算法
//3. 实现多线程，考虑超大数据集
//4. 方便的测试数据集的导入、显示可视化方案，方便地与osg进行对接
//5. 方便的不同算法正确性、效率比较方案（如何对比不同算法，对比每个三角形是否一致，或者简单看看重心坐标是否重合）
//6. 方便的进行CGAL的接口测试，用于熟悉半边数据结构

//开始前，对比一下CDT和CGAL两个库：
//1. CGAL的数据结构容器Compact_container<base>，参考STL对内存的分配做了优化，如果把数据看成链表，每个node节点都要存储
// 下一个节点指针，占用了8个字节(x64)， CGAL有可能使用了类似Union的方法，这个指针和数据，在内存的申请和释放过程中，有机会共用同一个内存地址（猜的，源码没有看太懂）
// CDT的话直接使用std::vector, std::map，相当于单独存储了指针和数据
//2. 准确性上面，还在对比验证中
//3. 两个库都使用了逐点插入算法，但是在point location上面的算法处理上不同
//(1)CGAL利用了本身半边数据结构的特点，由一个面出发，通过zig_zag找到目标面，
// 由于每次可以从四个方位中确定一个，感觉应该是log(n)，效率应该不会太低，他的point location和之后的插入操作，都是基于
// 同一个半边数据结构，没有CDT中的kdtree维护开销（插入操作开销，kd树本身搜索最近点的开销）
//(2)CDT使用了kd树加速point location过程，虽然在每次的搜索过程中，直接将顶点数据的引用传入，避免了一部分的开销；但是kd树本身也存储数据，
//该部分存储独立于顶点插入时的半边数据结构之外；感觉kd树可以较快的定位到最近点附近区域，但是这个要以维护kd树为代价；
//在找到了最近点之后，还是没有办法100%确定就是目标区域，还要对区域附近进行一次搜索判断，这个搜索的效率还要加以确认；
//感觉kd树在某些特殊的数据点集的情况下，会出现局部集中的问题（多个点集中在一个小区域内）
//4. CGAL库做了比较多的封装，使用上面比较难，不同算法可能用了不同的数据结构，CDT使用比较简单
//5. 从目前的使用情况看，二者好像都要传入一份数据的拷贝进去

//为什么要用devide and conquer？
//1. 自己的感觉：因为结构具有唯一性，如果一个小区域内（比如100个彼此靠近的凸包内）已经满足delaunay三角划分；
//那么在这个凸包外围插入点话感觉可以将flip edge的次数降到最低；但是使用逐点插入算法，感觉破坏了这样的局部特性，
//可能导致很多不必要的flip edge(除非在插入前对数据已经有了预判，进行了某种类似空间排序处理)：
//(1)对于大数据集的一种折中方法：对每个数据子集分布进行逐点插入（可以用多线程）；然后对之后的子集再进行一次devide and conquer
//融合子集结果（这个应该有人做过，后面找找文献）
//2. devide and conquer感觉好像天然有局部性优势，而且在内存的使用率上面和并发的编程难度上面好像都比较占优势

//数据结构的选型
//1. 排除了Triangulation_data_structure_2
//CGAL的Triangulation_data_structure_2是针对三角化的特殊数据结构，每次插入一个点，会把所有点连成三角形，所以需要找到插入面或线的位置
//而分治算法是把一些破碎的片段逐渐连接起来，允许过程中，存在一些孤立的点、线、面
//所以这个数据结构感觉很难适配，不过它的实现方法可以借鉴，本来想  定义局部的_tds，用逐点插入+分治结合，结果还是不行
//尴尬了！！！
//2. 还是考虑使用成熟的库，通用的数据结构，BGL有几个，如半边数据结构，
//3. 最简单的邻接表应该不行，因为简单的CCW，CW操作都很麻烦
//结论：用halfedge试试吧，CGAL有一个数据结构：HalfedgeDS
//参考1：https://doc.cgal.org/latest/HalfedgeDS/index.html#Chapter_Halfedge_Data_Structures
//参考2：https://doc.cgal.org/latest/BGL/group__PkgBGLTraits.html

//surface_mesh



//TODO: devide_and_conqure的边界条件整理(退化情况)
//1. 三点共线
//2. 四点共圆
//3. 考虑多线程版本： 如果两个子集是上下排列，此时没有上下切线

//TODO: devide_and_conqure算法现有的开源项目

//以上都是目前不成熟的理解，后面慢慢深入修正，主要目的还是想先把CGAL用起来
//TODO: 相关参考资料汇总
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

//模板编程搞不定，直接写类了
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

    typedef typename CGAL::Projection_traits_xy_3<K>                         Gt;    //先用这个进行计算

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
    std::set<OutContainerIttr, LexiSortOpr> m_setPointsIttrs;  //对数据点进行排序（lexicographically ascending）
    HDS* _pHds;
    Decorator* _pDecorator;
    Gt _gt;

public:
    void insert(OutContainerIttr begin, OutContainerIttr end)
    {
        //对插入的数据点集进行排序，方便后面的devide_and_conqure
        for (OutContainerIttr ittr = begin; ittr != end; ++ittr)
        {
            //如果能够顺便去重就好了
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

    //tool_2: 寻找两个凸包的上下切线
    TopAndBottomTangent findTangent(LimitHandlesPair leftHull, LimitHandlesPair rightHull)
    {
    }

    //插入三角形（当只有一个顶点时，其他两个顶点的句柄无效，但是仍然当三角形插入到数据结构中；一条边时，同理）
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

        //三点逆时针
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
    //返回左右两边极限位置的句柄
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
            //step1 : 操作_tds，将两点连线
            OutContainerIttr pr1 = begin;
            OutContainerIttr pr2 = begin + 1;
            OutContainerIttr pr3 = OutContainerIttr();
            return insert_hds(p1, p2, p3);
        }
        else if (size == 3)
        {
            //操作_tds，将两点连线
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
        //step1: 寻找两个凸包的上下切线（凸包可能是1个单独三角形、也可能是单独的一个点，也可能是单独的一条边，也可能是delaunay网格）   
        TopAndBottomTangent tangents = findTangent(LimitHandlesPair limitLeft, LimitHandlesPair limitRight)
        bool bLeft = is_valid_TangentLine(tangents.first);
        bool bRight = is_valid_TangentLine(tangents.first);
        if (!bLeft || bRight)
        {
            if (bLeft) return limitLeft;  
            if (bRight) return limitRight;
            return LimitHandlesPair();
        }

        //step2: 开始递进下切线，以满足delaunay三角划分
        //TODO: merge_callback
        while(!is_same_tangent(tangents.first, tangents.second))
        {
            //step2_1：以下切线的【左】顶点（位于左凸包）为基点，在【右】凸包中，寻找下切线的【候选点】
            //寻找过程中，必要时断开【右】凸包的某些边（这边后面不会再连起来了）


            //step2_2：以下切线的【右】顶点（位于左凸包）为基点，在【左】凸包中，寻找下切线的【候选点】
            //寻找过程中，必要时断开【左】凸包的某些边（这边后面不会再连起来了）


            //step2_3: 
            //从以上基点和候选点中，确定下切线, 在_tds中插入该边
            //循环递进到这条新的下切线

            
        }//step2_4: continue the while circle

        Vertex_handle LeftMostVertexHandle = limitLeft.first;
        Vertex_handle RightMostVertexHandle = limitRight.second;
        return LimitHandlesPair(LeftMostVertexHandle, RightMostVertexHandle);
    }

private:


};

