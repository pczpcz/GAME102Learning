#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/TDS_2/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>

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


//TODO: devide_and_conqure�ı߽���������




//TODO: devide_and_conqure�㷨���еĿ�Դ��Ŀ




//���϶���Ŀǰ���������⣬��������������������ҪĿ�Ļ������Ȱ�CGAL������
//TODO: ��زο����ϻ���
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
        //�Բ�������ݵ㼯�������򣬷�������devide_and_conqure
        //���ݵ㼯����ֱ�Ӳ��뵽osg::Vec3Array, �����osg�ν�
        _tds.insert(begin, end);
        sort(begin, end);
    }
    */

    /* һЩ���ߣ�
    //tool_1: ��������
    void sort(iterator begin, iterator end);

    //tool_2: Ѱ������͹������������
    std::pair<EdgeHandle, EdgeHandle> findTangent(iterator begin1, iterator end1, iterator begin2, iterator end2);

    //tool_3: ...


    */

    //devide_and_conqure
    /*α����
    void devide_and_conqure(iterator begin, iterator end) 
    {
        //Ĭ�����ݵ㼯�Ѿ�����������
        int size = std::distance(begin, end);
        if (size == 1 || size <= 0)
        {
            return;
        }
        else if (size == 2)
        {
            step1: ����_tds������������
            return;
        }
        else if (size == 2)
        {
            //����_tds������������
            _tds.colinear(begin, end)
            return;
        }
        else if (size == 3)
        {
            //����_tds������������������
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
    /*α����
    void merge(iterator begin1, iterator end1, iterator begin2, iterator end2)
    {
        //�Ͳ���iterator����Ч�����ж��ˣ�Ĭ����ok��
        //step1: Ѱ������͹�����������ߣ�͹��������1�����������Ρ�Ҳ�����ǵ�����һ���㣬Ҳ�����ǵ�����һ���ߣ�Ҳ������delaunay����

        //step2: ��ʼ�ݽ������ߣ�������delaunay���ǻ���
        //TODO: merge_callback
        while(������ != ������)
        {
            //step2_1���������ߵġ��󡿶��㣨λ����͹����Ϊ���㣬�ڡ��ҡ�͹���У�Ѱ�������ߵġ���ѡ�㡿
            //Ѱ�ҹ����У���Ҫʱ�Ͽ����ҡ�͹����ĳЩ�ߣ���ߺ��治�����������ˣ�


            //step2_2���������ߵġ��ҡ����㣨λ����͹����Ϊ���㣬�ڡ���͹���У�Ѱ�������ߵġ���ѡ�㡿
            //Ѱ�ҹ����У���Ҫʱ�Ͽ�����͹����ĳЩ�ߣ���ߺ��治�����������ˣ�


            //step2_3: 
            //�����ϻ���ͺ�ѡ���У�ȷ��������, ��_tds�в���ñ�
            //ѭ���ݽ��������µ�������

            
        }//step2_4: continue the while circle
    }
    */

private:
    //�̳еĻ����Ա��
    //Gt _gt;
    //Tds _tds;
    //Vertex_handle _infinite_vertex;

};

