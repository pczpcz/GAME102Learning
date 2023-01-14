#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/TDS_2/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>

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


//TODO: devide_and_conqure的边界条件整理




//TODO: devide_and_conqure算法现有的开源项目




//以上都是目前不成熟的理解，后面慢慢深入修正，主要目的还是想先把CGAL用起来
//TODO: 相关参考资料汇总
// [0]
// [1]
// [2]
// 
//author: zhanghcuanpan
//date: 2023.01.14

using namespace CGAL;

template < class Gt, class Tds = Triangulation_data_structure_2 <Triangulation_vertex_base_2<Gt>,Triangulation_face_base_2<Gt> > >
    class Delaunay_devide_and_conqure : public Triangulation_2
{
public:
    typedef Gt Geom_traits;
    typedef typename Geom_traits::Point_2                       Point;
    typedef typename Geom_traits::Segment_2                     Segment;
    typedef typename Geom_traits::Triangle_2                    Triangle;
    typedef typename Geom_traits::Orientation_2                 Orientation_2;
    typedef typename Geom_traits::Compare_x_2                   Compare_x;
    typedef typename Geom_traits::Compare_y_2                   Compare_y;
    typedef typename Geom_traits::Side_of_oriented_circle_2     Side_of_oriented_circle;

    typedef Triangulation_2<Gt, Tds>                            Triangulation;
    typedef typename Triangulation::size_type                   size_type;
    typedef typename Triangulation::Locate_type                 Locate_type;
    typedef typename Triangulation::Face_handle                 Face_handle;
    typedef typename Triangulation::Vertex_handle               Vertex_handle;
    typedef typename Triangulation::Edge                        Edge;
    typedef typename Triangulation::Edge_circulator             Edge_circulator;
    typedef typename Triangulation::Face_circulator             Face_circulator;
    typedef typename Triangulation::Vertex_circulator           Vertex_circulator;
    typedef typename Triangulation::Finite_edges_iterator       Finite_edges_iterator;
    typedef typename Triangulation::Finite_faces_iterator       Finite_faces_iterator;
    typedef typename Triangulation::Finite_vertices_iterator    Finite_vertices_iterator;
    typedef typename Triangulation::All_faces_iterator          All_faces_iterator;

    Delaunay_devide_and_conqure(Delaunay_devide_and_conqure&&) = default;
    Delaunay_devide_and_conqure& operator=(const Delaunay_devide_and_conqure&) = default;
    Delaunay_devide_and_conqure& operator=(Delaunay_devide_and_conqure&&) = default;
    ~Delaunay_devide_and_conqure() = default;

    template <class InputIterator>
    Delaunay_triangulation_2(InputIterator first, InputIterator last, const Gt& gt = Gt())
        : Triangulation_2<Gt, Tds>(gt)
    {
        //insert(first, last);
    }

    /*
    void insert(iterator begin, iterator end)
    {
        //对插入的数据点集进行排序，方便后面的devide_and_conqure
        //数据点集可以直接插入到osg::Vec3Array, 方便和osg衔接
        _tds.insert(begin, end);
        sort(begin, end);
    }
    */

    /* 一些工具：
    //tool_1: 顶点排序
    void sort(iterator begin, iterator end);

    //tool_2: 寻找两个凸包的上下切线
    std::pair<EdgeHandle, EdgeHandle> findTangent(iterator begin1, iterator end1, iterator begin2, iterator end2);

    //tool_3: ...


    */

    //devide_and_conqure
    /*伪代码
    void devide_and_conqure(iterator begin, iterator end) 
    {
        //默认数据点集已经做过排序了
        int size = std::distance(begin, end);
        if (size == 1 || size <= 0)
        {
            return;
        }
        else if (size == 2)
        {
            step1: 操作_tds，将两点连线
            return;
        }
        else if (size == 2)
        {
            //操作_tds，将两点连线
            _tds.colinear(begin, end)
            return;
        }
        else if (size == 3)
        {
            //操作_tds，将三点连成三角形
            _tds.coltriangal(begin, end)
            return;
        }

        //dv_callback

        //conditon: size > 3
        int mid = size/2;
        iterator begin1 = begin;
        iterator end1 = begin + mid;
        iterator begin2 = end1 + 1;
        iterator end2 = end;

        devide_and_conqure(begin1, end1);
        devide_and_conqure(begin2, end2);
        merge(begin1, end1, begin2, end2);
    }*/

    //merge function for devide_and_conqure
    /*伪代码
    void merge(iterator begin1, iterator end1, iterator begin2, iterator end2)
    {
        //就不对iterator的有效性做判断了，默认是ok的
        //step1: 寻找两个凸包的上下切线（凸包可能是1个单独三角形、也可能是单独的一个点，也可能是单独的一条边，也可能是delaunay网格）

        //step2: 开始递进下切线，以满足delaunay三角划分
        //TODO: merge_callback
        while(下切线 != 上切线)
        {
            //step2_1：以下切线的【左】顶点（位于左凸包）为基点，在【右】凸包中，寻找下切线的【候选点】
            //寻找过程中，必要时断开【右】凸包的某些边（这边后面不会再连起来了）


            //step2_2：以下切线的【右】顶点（位于左凸包）为基点，在【左】凸包中，寻找下切线的【候选点】
            //寻找过程中，必要时断开【左】凸包的某些边（这边后面不会再连起来了）


            //step2_3: 
            //从以上基点和候选点中，确定下切线, 在_tds中插入该边
            //循环递进到这条新的下切线

            
        }//step2_4: continue the while circle
    }
    */

private:
    //继承的基类成员：
    //Gt _gt;
    //Tds _tds;
    //Vertex_handle _infinite_vertex;

};

