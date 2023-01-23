#include "Delaunay_Mutithread_CGAL.h"
#include <cstdlib>

Delaunay_Mutithread_CGAL::Delaunay_Mutithread_CGAL()
{
	m_vecPointsRef = new osg::Vec3Array;
	//m_PrimitveSetRef = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
}

Delaunay_Mutithread_CGAL::~Delaunay_Mutithread_CGAL()
{
	for (int i = 0; i < m_vecTasks.size(); ++i)
		delete m_vecTasks[i];

	if (m_pThreadPool) 
	{
		m_pThreadPool->stop();
		delete m_pThreadPool;
	}
}

void Delaunay_Mutithread_CGAL::read_data()
{
}

int Delaunay_Mutithread_CGAL::getThreadNum(size_t size)
{
	int iMaxThread = std::thread::hardware_concurrency() * 2;
	if (iMaxThread <= 0)
	{
		iMaxThread = 2;
	}

	//���������
	int iVertexPerThread = 3;

	int iThreadNum = 2;
	iThreadNum = size / iVertexPerThread;
	if (iThreadNum <= 0)
		iThreadNum = 2;
	if (iThreadNum >= iMaxThread)
		iThreadNum = iMaxThread;

	return iThreadNum - 1;
}

/*�������Ҫ��ȷһ������������������ʹ�ó�����ununiform data_set�ȵ�
//����sample��
//1. ���������Χ���֣�ÿ����Χ�ڵĵ��������ܲ�ͬ�����䵽ÿ���̵߳Ķ������Ͳ�ͬ��ÿ���̵߳ĸ��ؾͲ�ͬ
//2. ����������򻮷֣�����ϲ���ʱ�򣬺���������ᶼ�����룬��ô�ϲ�Ҳ�����鷳
//3. �ⲻ����ʱ����Ҳ�и���ķ������߽�߽����û�а취��֤һֱ������ֱ�ģ�
//������������Կ��ǽ�����Ե����������̲��䣨ֻҪ�߽��Ƕ���ģ�������Ǳ��룩�����濴��������ʵ��
//����������ͨ��Ԥ�Ȳ����ķ�ʽ���Ժ�������ʷ֣���ֻ�Ժ��ᣬ����ֳɾŹ������ʽ�����ؾ���Ͷ�����Ѽ�ˣ���

//1. �������е�chunk���룬���ǣ�����Ǹ��߱������ضϻ������ڷ��������棬��ʲôӰ��
//2. ���Ƿ�����Լ���ߵ���ʲôӰ��
*/
void Delaunay_Mutithread_CGAL::delaunay(osg::Vec3Array::iterator begin, osg::Vec3Array::iterator end)
{
	size_t size = std::distance(begin, end);
	if (size <= 0)
		return;

	int iThreadNum = getThreadNum(size);

	//����
	iThreadNum = 2;

	unsigned int vertexPerthread = size / iThreadNum;
	std::thread *pThreads = new std::thread[iThreadNum];

	if (!m_pThreadPool)
		m_pThreadPool = new CThreadPool;
	m_pThreadPool->setThreadNum(iThreadNum);

	//����Ԥ��������ʵ�ָ��ؾ���
	//sample(0.5, iThreadNum, begin, end);
	sample(1.0, iThreadNum, begin, end);

	if (m_vecSectionData.empty())
		return;

	//��� ��ʼ���� �� �̳߳ص��������
	for (int i = 0; i < iThreadNum; ++i)
	{
		if (m_vecSectionData[i].m_vecDataIttrs.empty())
			continue;
		delaunay_task  *pTask = new delaunay_task(this, &m_vecSectionData[i]);
		m_vecTasks.push_back(pTask);
		m_pThreadPool->addTask(pTask);
	}
	m_pThreadPool->start();
	m_pThreadPool->join();
}

//ԭ�������
// 
// ��ע�����������Ƶ��������Ķ����γɵķ����������Ƚ϶�Ĳ²�ɷ�
// 
// �ο�����1��https://dsa.cs.tsinghua.edu.cn/~deng/cg/project/2021s/2021s-i.pdf
// �ο�����2��D. Funke and P. Sanders. Parallel d-d delaunay trian-gulations in sharedand distributed memory.In 2017 Proceedings of the Ninteenth Workshop on Algorithm Engineering and Experiments(ALENEX), pages 207�C217. SIAM, 2017.
// 
// ������̣�
//1. �����õ��Ӽ������������+��ԭ������ͬ�������Σ�
//2. �����õ��Ӽ������������+��ԭ������ͬ�������Σ�
//3. ��ԭ���֣��޳��ǰ�ȫ�棩�У�ͨ���Ӽ��Ķ������ѭ�����������뻹ȱ�ٵ�������
//�е㿴�������������ԣ�

// ����һ����û��һ�ֿ��ܣ���������Σ�p1,p2,p4������֮ǰԭ���е������Σ�p1,p2,p3��;�����������β����ھӹ�ϵ�����ǽ����ϵ
// ���������,p2,p3,p4���γ���һ����������Σ���ԭ��������p1,p2���ھ�ʱ���Ƿ���ܻ��(p1,p2,p3��������ӽ���)��
// ���ۣ�������ԭ����������Χ�У���Ӧ�ð�(p1, p2, p3)���޳��������������ʵ���Ƿǰ�ȫ��
/*
	p1	#
		|  \
		|    #p4
		|  /
	p2	#  #p3
*/

//���������Ƿ���ڱ���������ΰ�Χ����������
//����p00,p01,p1,p2,p3,p4λ��һ����ߣ���p5,p6��p7λ������һ����ߣ���ʱ��p1,p2,p5����p3��p4, p6���ǿ��������
//p2��p3��p5��p6֮��һ�����л��п�������Σ�����ԭ��͹������Χ�߽磬һ����һϵ�л���һȦ��������Σ��������ڣ�
//1. (p1,p2,p5)��������β�û��Ӱ�쵽ԭ����p00, p1, p2�����ԭ���������Σ�û������
//2. (p01,p7,p2)�ƻ���ԭ����������(p01,p03,p02),�������ڲ�
//���ۣ�
//1. ��������κͿ��������֮��һ���������ſ�������Σ�����������͹����Ե�����Ų������������û�����뵽ԭ�����ڲ�
//2. �ο�����1��û�����������ֻ������/û�����룬����У���������һ����������Σ���������������һ�����������ķ��ߣ����������߷ָ�����
//3. û������ʱ���ϲ����������ǻ��������Σ�p00, p1, p2����ԭ��һ�£�����ԭ���б�����������Σ���ȻӦ����ԭ����������Χ�б��޳�

/*
	   p1 # --------# p5
       /   \       /
   p00*		\     /
	    \    \   /
		    p2 #
             /  \
		p01# ---- #p7
		    \   /
		   p3 # 
      		 / \	
		 p4 # - #p6
*/

//��������ͨ���⽻Բ��Աߵ�bounding box����ʶ��ǰ�ȫ��ķ������Ƿ�׼ȷ������ϲ�ʱ���Ѱ�ȫ��Ҳ������ȥ����ʲôӰ��û�У�
//1. ʶ������ķǰ�ȫ�����Χ�������»���֮ǰ��Ӧ�ò���͹���������������ǻ�֮��һ����һ��͹�������������۱�֤�ģ����������������ǻ�֮���͹����Χ����Щ�����Σ�������ǲ��Ե�
//2. Ҫ��һ���޳������淽����Ҫ��Ȼmerging��û�취�����ˣ�������û�а�������İ�ȫ���ȥ�������Ǿ�ȷʶ���˷ǰ�ȫ�棬���»��ֺ����Χ�����ζ������ǲ�׼ȷ��
//3. ���Զ���merging��˵��������������Ӧ�ò�������
//4. �ɴ˲����������⣺merging��ʱ������޳�����Щ����ȷ����

/*
		   ��ʶ������ķ�ȫ��߽硿                                    ���ϲ����͹���߽硿
			--------------------			
			\		           /
			 \		          /
			  \              /            re_triangulation()            �϶���һ��͹�����ԣ�
     safe	  |    risk     |  safe              --->
			  |             |
			  /             \
			 /               \
		   /                  \
		   ---------------------
*/

//�����ģ�
//1. ���÷������ĵĽ��ۣ��ڿ�Խ�ֽ����ϵĿ���������������Ų��ģ����ԾͿ�����������߽�Ϊ��㣬���������������ָ�����ȷ�Ľ��
//2. ��λָ������Իָ���߽�Ϊ����
/* α�������£���ʵ��ʵ�ֲ��������ģ�ֻ�����ڷ�����
	for (�����еĿ����������ȡ��һ�����������)
	{
		for(�ҵ����������ھ�)
		{
			//�������������Ƿ�����ģ�����һ����ԭ�����о��Ѿ���������������Σ�ʵ��д�����ʱ����������Ѵ�����������Ϊ�Ƿ�������ж�����
			if (��ԭ�����о��Ѿ�����--�������)
			{
				//����ھ��ǰ�ȫ��
			}
			else
			{
				if (����ھ�����һ�����������)
				{
					//�ǰ�ȫ��
				}
				else
				{
					//situation_1������ھ�����Ϊ��������ʶ��׼ȷ���������뵽�������У����������������ǻ�͹�����ƣ�
					//�����ԭ����������ȷ�Ļ��֣��������ǻ��б��𿪽��д���Ļ��֣������ǰ�ȫ�棬Ҳ�����Ƿǰ�ȫ�棩
					//situation_2: ����ھ��������»��ֺ��͹���ı߽��ϣ�������Ƿ���ȷ�棨�����ǰ�ȫ�棬Ҳ�����Ƿǰ�ȫ�棩

					//��ԭ�����еļ��Ͻ����˷�����
					//1. ԭ�����漯�� = �����漯�� + �Ƿ����漯��
					//2. �����漯�� = �������漯�� + δ�������漯�� + ����漯��

					//������Ҫ���漯�� = ԭ���ּ��� - �����漯�� + δ�������漯�� + ����漯��
					//���У�ԭ���ּ��ϣ������漯�ϣ�����漯�϶�����֪�ģ�����δ�������漯�ϣ�ò�ƿ���ͨ���жϷ����漯���е�����ԭ�������Ƿ���֣���ȷ���Ƿ���δ�������棬
					//�Ⲣδ�������е�δ�������棬��Ϊsituation_1��situation_2�Ĵ��ڡ�

					//�����ƵĽ�� = ԭ���ּ��� - �����漯�� + �����漯�����»��ֺ��ԭ�����漯���ظ����� + ����漯��

					//����ѶԷ����漯�ϵ��������ǻ����̣�Ҳ����ɴ�͹����Χ��������̣�����������Ȳ���̫���Ϊ�������漯�����»��ֺ��ԭ�����漯���ظ����� + ����漯�ϣ�����Ѿ�
					//ȷ���˴󲿷���ȷ�����ǻ��֣��������������������Χ������µĵ㣬�Ӷ����µ����룬�ɴ˿��Բ²�����ۣ�û�������Ƶ���ֻ�ܽв��ۣ���
					//�����ۡ�������situation_1����situation_2�б����󻮷ֵ������ζ�Ӧ�ú������Ѿ�ʶ������������Σ������漯�����»��ֺ��ԭ�����漯���ظ����� + ����漯�ϣ��й����ߣ�
					//������Ӧ���Ǵ�ģ�Ӧ���ǣ���ʶ�������ȷ���֣������ģ�������Χ���󻮷֣������ģ�֮�䣬�������ı߽磬�����Ϳ����������߽�����Χ�ָ�����ȷ����

					//�������ȷ���ֵ���Ͳ�����ԭ�����У���ôͨ����������߰����ҳ�����һ��������Ӧ�����棬����һ����ֻ�����ֿ��ܣ���1����ʶ����ظ��棬��2����ʶ��Ŀ���棻
					//����һ�����������Ҫ�ҵ��棬������׼ȷ�Ľ�����λ���������ǻ�ǰ�ķ����漯���У�������������ְ�����һЩ����ȫ�棬������������ѭ��������

					//���Ǻ���û�У�
					//1. ����������ڷ����漯�����������汻���»��֣�flip��,�����Ǳ�flip�������ߣ��������߾Ͳ����������ʶ��Ŀ���������ʶ����ظ����У�
					//2. �����������û�б�flip�����Ǳ����룬��ô������Ҳ�����������ʶ��Ŀ���������ʶ����ظ����У�
					//3. �������������Ϊ����������߱�flip��������ߣ���ô������һ���� ��ʶ����ظ����������ʶ��Ŀ�����ϵı�
					//4. �����ظ����ϵıߣ�����ʶ����ظ����ϵı߳�������ԭ���ּ������������������ߵ��棬Ҫô�ҵ����Լ���Ҫô���ҵ����������ھ�
					//5. ���ڿ�����ϵıߣ�����ʶ��Ŀ�����ϵı߳�������ԭ���ּ������������������ߵ��棬�����������������1��֮��Ҫ��������棨2���������ھӣ� Ҫ��������أ�
					//5.1 ��������Χ�еĽ�Ҫ���������������ǰɾ��������ɾ��������̫��ɾ��
					//5.2 ������һ���Ƕȣ�Ҫ������������ζ��ǿ�������εĶ��㣬�����еĶ�����ɼ��ϣ��κ���������λ���������֮�ڵ������Σ�����Ӧ�����ҵ��������Ƿ�Χ֮��
				}
			}
		}
	}

	//�����壺���ڷǰ�ȫ���ѡ����������һ�������ˣ�
	//��͹������Χ���ټ���һ���㣬���Ӱ��
	//1. ��p1,p2,p3,p4�����͹���߽磬���붥��5�󣬻�Ӱ�쵽�����棺����(p3,p4), Ϊ���γ�͹����(p2,p3)���ܻ���Ӱ��
	//���⣺��͹������Χ�ټ���һ���㣬���뵽����û�к�ȣ�����˵ֻӰ�쵽��һȦ???

	p1 *
		 \
		p2 *
		   |
		p3 *
		  /  * p5
	  p4 *
*/

void Delaunay_Mutithread_CGAL::saveFaceIndexs(Face_handle fh) 
{
	std::lock_guard<std::recursive_mutex> lock(m_mutex);
	if (!m_PrimitveSetRef)
		m_PrimitveSetRef = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
	m_PrimitveSetRef->push_back(fh->vertex(0)->info());
	m_PrimitveSetRef->push_back(fh->vertex(1)->info());
	m_PrimitveSetRef->push_back(fh->vertex(2)->info());
}

void Delaunay_Mutithread_CGAL::saveFaceIndexs(unsigned int index1, unsigned int index2, unsigned int index3) 
{
	std::lock_guard<std::recursive_mutex> lock(m_mutex);
	if (!m_PrimitveSetRef)
		m_PrimitveSetRef = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
	m_PrimitveSetRef->push_back(index1);
	m_PrimitveSetRef->push_back(index2);
	m_PrimitveSetRef->push_back(index3);
}

Delaunay_Mutithread_CGAL::SSectionData::Orientation Delaunay_Mutithread_CGAL::orientation(SSectionData& section1, SSectionData& section2)
{
	return section1.orientation(section2);
}

void Delaunay_Mutithread_CGAL::merge(SSectionData* section1, SSectionData* section2)
{
	if (!section1 || !section2)
		return;

	SSectionData deplicateSection(*section1, *section2);
	Point2D p1((section1->m_dRight + section2->m_dLeft) / 2, deplicateSection.m_dTop);
	Point2D p2((section1->m_dRight + section2->m_dLeft) / 2, deplicateSection.m_dBottom);
	Segment_2 seg(p1, p2);

	//1. ������ڷ������ཻ���յ������ʶ��
	for (auto face_handle : section1->m_dt.finite_face_handles())
	{
		//�������������Բ �� �����������������Ƿ��н��棬�����, ���������Ϊ���з��յ���
		if (!section1->check_circle_intersect(face_handle, seg, SSectionData::Near_Right)) {
			//saveFaceIndexs(face_handle);	//for debug
		}
	}
	for (auto face_handle : section2->m_dt.finite_face_handles())
	{
		//�������������Բ �� �����������������Ƿ��н��棬�����, ���������Ϊ���з��յ���
		if (!section2->check_circle_intersect(face_handle, seg, SSectionData::Near_Left)) {
			//saveFaceIndexs(face_handle);  //for debug
		}
	}

	//�Էǰ�ȫ�����½������ǻ�
	deplicateSection.insert_risk_delaunay(*section1, SSectionData::Near_Right);
	deplicateSection.insert_risk_delaunay(*section2, SSectionData::Near_Left);
	deplicateSection.insert_tangent_delaunay(*section1, *section2);

	//��ѡ����ȫ�棬�����������ظ���Ϊ�ϲ���׼��
	for (Face_handle face_handle : deplicateSection.m_dt.finite_face_handles())
	{
		//saveFaceIndexs(face_handle);       //for debug

		//����Ƿ��ߵ���, ��ϲ�ǰ�����ظ�����
		if (section1->check_face_intersect(face_handle, seg)
			|| section2->check_face_intersect(face_handle, seg)
			|| section1->check_face_duplicate(face_handle)
			|| section2->check_face_duplicate(face_handle))
		{
			deplicateSection.m_deqBarrierFaces.push_back(face_handle);
			deplicateSection.m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(face_handle->vertex(0)->info(), face_handle->vertex(0)));
			deplicateSection.m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(face_handle->vertex(1)->info(), face_handle->vertex(1)));
			deplicateSection.m_vertexHashtable.insert(std::pair<unsigned int, Vertex_handle>(face_handle->vertex(2)->info(), face_handle->vertex(2)));
			saveFaceIndexs(face_handle);
		}
	}

	//�ϲ������������ǻ�������ϲ�ǰ�Ľ��бȶԣ�
	std::vector<Face_handle> vecMergeResult;
	deplicateSection.merge(vecMergeResult, *section1, *section2);
	for (auto& handle : vecMergeResult)
	{
		saveFaceIndexs(handle);
	}

	section1->m_bRightFinishMerge = true;
	section2->m_bLeftFinishMerge = true;
	section1->m_status = SSectionData::Status_Finish;
	section2->m_status = SSectionData::Status_Finish;

	if (section1->m_bLeftFinishMerge && section1->m_bRightFinishMerge) 
	{
		for (auto face_handle : section1->m_dt.finite_face_handles()) {
			auto ittr_left = section1->m_setRiskHandles_left.find(face_handle);
			auto ittr_right = section1->m_setRiskHandles_right.find(face_handle);
			if (ittr_left == section1->m_setRiskHandles_left.end() && ittr_right == section1->m_setRiskHandles_right.end()) {
				saveFaceIndexs(face_handle);
			}
		}
	}

	if (section2->m_bLeftFinishMerge && section2->m_bRightFinishMerge)
	{
		for (auto face_handle : section2->m_dt.finite_face_handles()) {
			auto ittr_left = section2->m_setRiskHandles_left.find(face_handle);
			auto ittr_right = section2->m_setRiskHandles_right.find(face_handle);
			if (ittr_left == section2->m_setRiskHandles_left.end() && ittr_right == section2->m_setRiskHandles_right.end()) {
				saveFaceIndexs(face_handle);
			}
		}
	}

	//����Ƿ����µĺϲ�������Ӻϲ������̳߳أ�������Ƿ����ȫ��������ֹͣ�̳߳أ�
	checkMergingAndFinish();
}

void Delaunay_Mutithread_CGAL::delaunay_CGAL(SSectionData *sectionData)
{
	if (!sectionData || sectionData->m_vecDataIttrs.empty())
		return;

	//�Ƚ���һ�����ǻ�
	for (auto &ittr : sectionData->m_vecDataIttrs)
	{
		Vertex_handle vh = sectionData->m_dt.insert(Point2D(ittr->x(), ittr->y()));
		vh->info() = std::distance(m_vecPointsRef->begin(), ittr);  //
	}

	if (m_vecSectionData.size() == 1)
	{
		for (auto& face_handle : m_vecSectionData[0].m_dt.finite_face_handles()) 
		{
			saveFaceIndexs(face_handle);
		}
	}

	sectionData->m_status = SSectionData::Status_Tri;

	//����Ƿ����µĺϲ�������Ӻϲ������̳߳أ�������Ƿ����ȫ��������ֹͣ�̳߳أ�
	checkMergingAndFinish();

	return;
}

void Delaunay_Mutithread_CGAL::checkMergingAndFinish()
{
	//����Ƿ��п��Ժϲ��ķ�����>=2��,����У�����ӵ��������
	std::set<std::pair<int, int>> setMergeIndexs;
	if (hasMerging(setMergeIndexs))
	{
		for (auto& pair : setMergeIndexs)
		{
			merging_task* pTask = new merging_task(this, &m_vecSectionData[pair.first], &m_vecSectionData[pair.second]);
			m_vecTasks.push_back(pTask);
			m_pThreadPool->addTask(pTask);
		}
	}

	//����Ƿ������ȫ������
	if (checkAllFinish() && m_pThreadPool->empty())
		m_pThreadPool->stop();
}

int Delaunay_Mutithread_CGAL::hasMerging(std::set<std::pair<int, int>>& setMergingIndexs)
{
	std::lock_guard<std::recursive_mutex> lock(m_mutex);

	if (m_vecSectionData.size() < 2)
		return 0;

	std::set<int> setIndexs;

	for (int i = 0; i < m_vecSectionData.size() - 1; ++i)
	{
		for (int j = i + 1; j < m_vecSectionData.size(); ++j)
		{
			if (m_vecSectionData[i].m_status == SSectionData::Status_Merging
				|| m_vecSectionData[j].m_status == SSectionData::Status_Merging
				|| m_vecSectionData[i].m_status < SSectionData::Status_Tri
				|| m_vecSectionData[j].m_status < SSectionData::Status_Tri)
			{
				continue;
			}

			SSectionData::Orientation ori = orientation(m_vecSectionData[i], m_vecSectionData[j]);
			if (ori == SSectionData::Near_Left)
			{
				if (!m_vecSectionData[i].m_bRightFinishMerge && !m_vecSectionData[j].m_bLeftFinishMerge)
				{
					if (setIndexs.insert(i).second && setIndexs.insert(j).second)
					{
						setMergingIndexs.insert(std::make_pair(i, j));
						m_vecSectionData[i].m_status = SSectionData::Status_Merging;
						m_vecSectionData[j].m_status = SSectionData::Status_Merging;
					}
				}
				break;
			}
			else if (ori == SSectionData::Near_Right)
			{
				if (!m_vecSectionData[j].m_bRightFinishMerge && !m_vecSectionData[i].m_bLeftFinishMerge)
				{
					if (setIndexs.insert(i).second && setIndexs.insert(j).second)
					{
						setMergingIndexs.insert(std::make_pair(j, i));
						m_vecSectionData[i].m_status = SSectionData::Status_Merging;
						m_vecSectionData[j].m_status = SSectionData::Status_Merging;
					}
				}
				break;
			}
		}
	}
	return setMergingIndexs.size();
}


bool Delaunay_Mutithread_CGAL::checkAllFinish() 
{
	std::lock_guard<std::recursive_mutex> lock(m_mutex);

	//m_vecSectionData�������������߳�����ǰ���Ѿ�������
	int iFinishCount = 0;
	for (int i = 0; i < m_vecSectionData.size(); ++i) 
	{
		if (!m_vecSectionData[i].m_bLeftFinishMerge)
			++iFinishCount;
		if (!m_vecSectionData[i].m_bRightFinishMerge) 
			++iFinishCount;
	}
	return iFinishCount <= 0;
}

void Delaunay_Mutithread_CGAL::sample(float fSampleRatio, int iThreadNum, osg::Vec3Array::iterator begin, osg::Vec3Array::iterator end)
{
	//���ڵ�ȱ�ݣ�
	//1. fSampleRatio, iThreadNum����ѡ��������ڵ���ǰ���Լ�����ȷ�ϣ����������������Ż������ֲ����㣬ֱ�ӷ���
	//2. ����ܹ���ǰ֪�������Сֵ��Χ�����ֻ���Ӿ��⣬Ŀǰ�����СֵҲ���в����Ļ�õ�
	//3. ��������������һ���߳�ֻ��һ���㣬����û������

	unsigned int size = std::distance(begin, end);
	if (size <= 0)
		return;

	if (fSampleRatio <= 0 || iThreadNum <= 0)
		return;

	if (fSampleRatio > 1.0)
		fSampleRatio = 1.0;

	m_vecSectionData.resize(iThreadNum);
	if (iThreadNum == 1)
	{
		for (auto ittr = begin; ittr != end; ++ittr)
		{
			m_vecSectionData[0].m_iSampleBeginX = 0;
			//m_vecSectionData[0].m_status = SSectionData::Status_Idle;
			m_vecSectionData[0].m_iSampleEndX = 1;	 //��ʱֻ���Ǻ���Ļ���
			m_vecSectionData[0].boundingBox(ittr->x(), ittr->y());
			m_vecSectionData[0].m_vecDataIttrs.push_back(ittr);
			m_vecSectionData[0].m_bLeftFinishMerge = true;
			m_vecSectionData[0].m_bRightFinishMerge = true;
		}
		return;
	}

	unsigned int vertexPerthread = size / iThreadNum;
	if (vertexPerthread <= 0)
		return;

	for (int i = 0; i < m_vecSectionData.size(); ++i)
	{
		m_vecSectionData[i].m_vecDataIttrs.reserve(vertexPerthread + 1);
		m_vecSectionData[i].m_iSampleBeginX = i;
		m_vecSectionData[i].m_iSampleEndX = (i + 1);
		//m_vecSectionData[i].m_status = SSectionData::Status_Idle;
	}

	//1. �������ȷ�����ݷֲ�(��η�ֹ�����)
	std::set<osg::Vec3Array::iterator> setSamples;
	unsigned int totalSampleNum = fSampleRatio * size;
	srand(time(0));
	unsigned int randOffset = 0;
	double dMax_x = MIN_DOUBLE;
	double dMin_x = MAX_DOUBLE;
	for (; setSamples.size() < totalSampleNum; )
	{
		randOffset = rand() % size;
		//randOffset = i;
		auto ittr = begin + randOffset;
		setSamples.insert(ittr);

		if (ittr->x() > dMax_x)
			dMax_x = ittr->x();
		if (ittr->x() < dMin_x)
			dMin_x = ittr->x();
	}

	totalSampleNum = setSamples.size();

	unsigned int slot = 10000;
	double range = dMax_x - dMin_x;
	double slot_distance = range / slot;
	if (range < 0 || slot_distance < 0 || abs(slot_distance) < 10e-6)
		return;

	std::vector<unsigned int> vecSlots(slot, unsigned int(0));
	for (auto &ittr : setSamples)
	{
		unsigned int slotIndex = ((ittr->x() - dMin_x) / slot_distance - 0.5);
		if (slotIndex < 0 || slotIndex >= slot)
			return;

		++vecSlots[slotIndex];
	}

	//2. ���ݷֲ�����������򻮷�
	typedef std::pair<double, double> DivRange;
	std::vector<DivRange> vecDivRanges;
	unsigned int average = (totalSampleNum * 1.0 / (unsigned int)iThreadNum) + 0.5;
	if (average <= unsigned int(0))
		return;

	unsigned int total = unsigned int(0);
	unsigned int rangeIndex = 0;
	for (int slotIndex = 0; slotIndex < vecSlots.size(); ++slotIndex)
	{
		if (vecSlots[slotIndex] <= unsigned int(0))
			continue;

		if (rangeIndex >= iThreadNum)
			break;

		total += vecSlots[slotIndex];
		for (; total >= rangeIndex * average; )
		{
			double low = dMin_x + slotIndex * slot_distance; //��������
			double high = dMin_x + (slotIndex + 1) * slot_distance;
			if (!vecDivRanges.empty()) {
				vecDivRanges[rangeIndex-1].second = low;
			}
			vecDivRanges.push_back(DivRange(low, high));
			++rangeIndex;
			break;
		}
	}

	//3. ��ʽ����
	if (vecDivRanges.size() != iThreadNum)
		return;

	auto getDivNo = [&](double value)->int
	{
		for (int i = 0; i < vecDivRanges.size()-1; ++i)
		{
			if (value < vecDivRanges[i].second)
				return i;
		}
		
		return vecDivRanges.size()-1;
	};

	for (auto ittr = begin; ittr != end; ++ittr) 
	{
		int divNo = getDivNo(ittr->x());
		if (divNo < 0 || divNo >= m_vecSectionData.size()) 
		{
			m_vecSectionData.clear();
			return;
		}

		m_vecSectionData[divNo].boundingBox(ittr->x(), ittr->y());
		m_vecSectionData[divNo].m_vecDataIttrs.push_back(ittr);
	}

	if (!m_vecSectionData.empty()) 
	{
		if (m_vecSectionData.size() == 1) 
		{
			m_vecSectionData[0].m_bLeftFinishMerge = true;
			m_vecSectionData[0].m_bRightFinishMerge = false;
		}
		m_vecSectionData[0].m_bLeftFinishMerge = true;
		m_vecSectionData[m_vecSectionData.size()-1].m_bRightFinishMerge = true;
	}

	return ;
}

void Delaunay_Mutithread_CGAL::readInputFromFile(const std::string& fileName)
{
	std::ifstream f(fileName);
	if (!f.is_open())
	{
		return;
	}
	int nVerts;
	int nEdges;
	f >> nVerts >> nEdges;

	// Read vertices
	std::vector<double> vv;
	vv.reserve(nVerts);
	for (std::size_t i = 0; i < nVerts; ++i)
	{
		double x, y;
		f >> x >> y;
		m_vecPointsRef->push_back(osg::Vec3(x, y, 0.0));
	}
	/*
	// Read edges
	std::vector<int> ee;
	for (std::size_t i = 0; i < nEdges; ++i)
	{
		int v1, v2;
		f >> v1 >> v2;
		ee.emplace_back(v1, v2);
	}
	*/
	return ;
}

osg::ref_ptr<osg::Geometry> Delaunay_Mutithread_CGAL::createGeometry()
{
	osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;
	if (m_PrimitveSetRef) 
	{
		geometry->setVertexArray(m_vecPointsRef.get());
		geometry->addPrimitiveSet(m_PrimitveSetRef.get());
	}
	return geometry.release();
}

void Delaunay_Mutithread_CGAL::delaunay()
{
	delaunay(m_vecPointsRef->begin(), m_vecPointsRef->end());
}

int Delaunay_Mutithread_CGAL::getRegularTeatExmapleData(int iVertexNums, float fRatio)
{
	//Ŀ�ģ�����100%���ظ����ݣ�Լ�������ཻ�����ݼ�
	int iRowTotal = 5;
	int iColTotal = 5;
	float fConstrainRatio = fRatio;
	if (fConstrainRatio < 0 || fConstrainRatio >= 0.9)
	{
		fConstrainRatio = 0.5;
	}

	if (iVertexNums > 0)
	{
		float num = sqrt(iVertexNums);
		iRowTotal = num + 0.5;
		iColTotal = num + 0.5;
	}

	if (iRowTotal < 5) iRowTotal = 5;
	if (iColTotal < 5) iColTotal = 5;

	int index = 0;
	float fDilute = 0.01;
	for (int iRow = 0; iRow < iRowTotal; iRow += 1/*iRowTotal*fDilute*/)
	{
		for (int iCol = 0; iCol < iColTotal; iCol += 1/*iColTotal*fDilute*/)
		{
			index = iCol + iRow * iColTotal;
			m_vecPointsRef->push_back(osg::Vec3((float)iRow, (float)iCol, 0.0));
		}
	}

	int iConstrainEdgesNum = iRowTotal * iColTotal * fConstrainRatio;
	int iCount = 0;
	int iEdgeStartIndex = 0;
	int iEdgeEndIndex = 0;
	for (int iRow = 0; iRow < iRowTotal - 1; ++iRow)
	{
		for (int iCol = 0; iCol < iColTotal - 1; iCol += 2)
		{
			iEdgeStartIndex = iRow + iCol * iRowTotal;
			iEdgeEndIndex = iRow + (iCol + 1) * iRowTotal;

			++iCount;
			//m_vecEdges.push_back(std::pair<int, int>(iEdgeStartIndex, iEdgeEndIndex));

			if (iCount >= iConstrainEdgesNum)
				return 0;
		}
	}

	if (iCount < iConstrainEdgesNum)
	{
		iEdgeStartIndex = 0;
		iEdgeEndIndex = 0;
		for (int iCol = 0; iCol < iColTotal - 1; iCol += 2)
		{
			for (int iRow = 0; iRow < iRowTotal - 1; ++iRow)
			{
				iEdgeStartIndex = iCol + iRow * iColTotal;
				iEdgeEndIndex = iCol + (iRow + 1) * iColTotal;

				++iCount;
				//m_vecEdges.push_back(std::pair<int, int>(iEdgeStartIndex, iEdgeEndIndex));

				if (iCount >= iConstrainEdgesNum)
					return 0;
			}
		}
	}

	return 0;
}

/*
void Delaunay_Mutithread_CGAL::DoEx(mono::tool::SElevationDomainData * dataset, const std::string & path, bool bValidBoundaryHeight, double dMaxTriangleEdgeLength, bool bNeedElementConvex)
{
}
*/
