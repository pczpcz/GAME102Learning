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

	//最合适数量
	int iVertexPerThread = 3;

	int iThreadNum = 2;
	iThreadNum = size / iVertexPerThread;
	if (iThreadNum <= 0)
		iThreadNum = 2;
	if (iThreadNum >= iMaxThread)
		iThreadNum = iMaxThread;

	return iThreadNum - 1;
}

/*后面可能要明确一下需求：最大的数据量，使用场景，ununiform data_set等等
//关于sample：
//1. 如果按规则范围划分，每个范围内的点数量可能不同，分配到每个线程的顶点数就不同，每个线程的负载就不同
//2. 如果按不规则划分，如果合并的时候，横轴或者纵轴都不对齐，那么合并也会有麻烦
//3. 外不调用时可能也有更大的分区，边界边界可能没有办法保证一直都是竖直的，
//此种情况，可以考虑将坐标对调，处理流程不变（只要边界是对齐的，看情况非必须），后面看看在哪里实现
//初步方案：通过预先采样的方式，对横轴进行剖分（且只对横轴，如果分成九宫格的形式，负载均衡和对齐很难兼顾）；

//1. 根据已有的chunk代码，考虑：如果登高线被分区截断或者落在分区线上面，有什么影响
//2. 考虑分区对约束边的有什么影响
*/
void Delaunay_Mutithread_CGAL::delaunay(osg::Vec3Array::iterator begin, osg::Vec3Array::iterator end)
{
	size_t size = std::distance(begin, end);
	if (size <= 0)
		return;

	int iThreadNum = getThreadNum(size);

	//测试
	iThreadNum = 3;

	unsigned int vertexPerthread = size / iThreadNum;
	std::thread *pThreads = new std::thread[iThreadNum];

	if (!m_pThreadPool)
		m_pThreadPool = new CThreadPool;
	m_pThreadPool->setThreadNum(iThreadNum);

	//进行预采样，以实现负载均衡
	//sample(0.5, iThreadNum, begin, end);
	sample(1.0, iThreadNum, begin, end);

	if (m_vecSectionData.empty())
		return;

	//添加 初始任务 到 线程池的任务队列
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

//原理分析：
// 
// 备注：不是理论推导，资料阅读后形成的分析，包含比较多的猜测成分
// 
// 参考资料1：https://dsa.cs.tsinghua.edu.cn/~deng/cg/project/2021s/2021s-i.pdf
// 参考资料2：D. Funke and P. Sanders. Parallel d-d delaunay trian-gulations in sharedand distributed memory.In 2017 Proceedings of the Ninteenth Workshop on Algorithm Engineering and Experiments(ALENEX), pages 207–217. SIAM, 2017.
// 
// 大概流程：
//1. 首先拿到子集（跨边三角形+与原划分相同的三角形）
//2. 首先拿到子集（跨边三角形+与原划分相同的三角形）
//3. 在原划分（剔除非安全面）中，通过子集的顶点进行循环搜索，补齐还缺少的三角形
//有点看不懂，举例试试：

// 分析一：有没有一种可能，跨边三角形（p1,p2,p4）代替之前原集中的三角形（p1,p2,p3）;而两个三角形不是邻居关系，而是交叉关系
// 这种情况下,p2,p3,p4会形成另一个跨边三角形，在原集中搜索p1,p2的邻居时，是否可能会把(p1,p2,p3给错误添加进来)？
// 结论：所以在原集的搜索范围中，就应该把(p1, p2, p3)给剔除掉，而这个面其实就是非安全面
/*
	p1	#
		|  \
		|    #p4
		|  /
	p2	#  #p3
*/

//分析二：是否存在被跨边三角形包围的特殊区域？
//假设p00,p01,p1,p2,p3,p4位于一个半边，而p5,p6，p7位于另外一个半边，此时（p1,p2,p5）（p3、p4, p6）是跨边三角形
//p2、p3、p5、p6之间一定还有还有跨边三角形，沿着原集凸包的外围边界，一定是一系列或者一圈跨边三角形，区别在于：
//1. (p1,p2,p5)跨边三角形并没有影响到原集（p00, p1, p2）这个原来的三角形，没有侵入
//2. (p01,p7,p2)破坏了原集的三角形(p01,p03,p02),侵入了内部
//结论：
//1. 跨边三角形和跨边三角形之间一定还存在着跨边三角形，而且是沿着凸包边缘连续排布，区别就是有没有侵入到原集的内部
//2. 参考结论1，没有特殊的区域，只有侵入/没有侵入，如果有，那它还是一个跨边三角形，跨边三角形组成了一道【连续】的防线，将左右两边分隔开来
//3. 没有侵入时，合并后重新三角化的三角形（p00, p1, p2）与原集一致；而在原集中被侵入的三角形，显然应该在原集的搜索范围中被剔除

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

//分析三：通过外交圆与对边的bounding box求交来识别非安全面的方法，是否准确，如果合并时，把安全面也包含进去了有什么影响没有？
//1. 识别出来的非安全面的外围，在重新划分之前，应该不是凸包，而在重新三角化之后，一定是一个凸包（这是有理论保证的），所以在重新三角化之后的凸包外围的那些三角形，大概率是不对的
//2. 要有一个剔除错误面方法，要不然merging就没办法工作了，无论有没有包含多余的安全面进去，就算是精确识别了非安全面，重新划分后的外围三角形都可能是不准确的
//3. 所以对于merging来说，分析三的问题应该不是问题
//4. 由此产生的新问题：merging的时候如何剔除掉这些不正确的面

/*
		   【识别出来的非全面边界】                                    【合并后的凸包边界】
			--------------------			
			\		           /
			 \		          /
			  \              /            re_triangulation()            肯定是一个凸包（略）
     safe	  |    risk     |  safe              --->
			  |             |
			  /             \
			 /               \
		   /                  \
		   ---------------------
*/

//分析四：
//1. 利用分析二的的结论：在跨越分界线上的跨边三角形是连续排布的，所以就可以这个连续边界为起点，向左右两边搜索恢复出正确的结果
//2. 如何恢复：（以恢复左边界为例）
/* 伪代码如下：（实际实现不是这样的，只是用于分析）
	for (从所有的跨边三角形中取出一个跨边三角形)
	{
		for(找到它的三个邻居)
		{
			//如果这个三角形是非侵入的，那它一定在原划分中就已经存在了这个三角形，实际写代码的时候是用这个已存在条件，作为是否侵入的判断条件
			if (在原划分中就已经存在--非侵入的)
			{
				//这个邻居是安全面
			}
			else
			{
				if (这个邻居是另一个跨边三角形)
				{
					//是安全面
				}
				else
				{
					//situation_1：这个邻居是因为：风险面识别不准确，而被加入到风险面中，而又由于重新三角化凸包限制，
					//这个在原划分中是正确的划分，在新三角化中被拆开进行错误的划分（可能是安全面，也可能是非安全面）
					//situation_2: 这个邻居属于重新划分后的凸包的边界上，大概率是非正确面（可能是安全面，也可能是非安全面）

					//对原划分中的集合进行了分析：
					//1. 原划分面集合 = 风险面集合 + 非风险面集合
					//2. 风险面集合 = 被入侵面集合 + 未被入侵面集合 + 跨边面集合

					//我们想要的面集合 = 原划分集合 - 风险面集合 + 未被入侵面集合 + 跨边面集合
					//其中，原划分集合，风险面集合，跨边面集合都是已知的，至于未被入侵面集合，貌似可以通过判断风险面集合中的面在原划分中是否出现，来确定是否是未被入侵面，
					//这并未包含所有的未被入侵面，因为situation_1，situation_2的存在。

					//不完善的结果 = 原划分集合 - 风险面集合 + 风险面集合重新划分后和原划分面集合重复的面 + 跨边面集合

					//如果把对风险面集合的重新三角化过程，也想象成从凸包外围的侵入过程，这个侵入的深度不会太深，因为（风险面集合重新划分后和原划分面集合重复的面 + 跨边面集合）这个已经
					//确定了大部分正确的三角划分，侵入过程是由于又在外围插入的新的点，从而导致的侵入，由此可以猜测出结论（没有能力推导，只能叫猜论）：
					//【猜论】：无论situation_1还是situation_2中被错误划分的三角形都应该和我们已经识别出来的三角形（风险面集合重新划分后和原划分面集合重复的面 + 跨边面集合）有公共边，
					//这个理解应该是错的，应该是：已识别出来正确划分（连续的），和外围错误划分（连续的）之间，有连续的边界，这样就可以以这条边界向外围恢复出正确划分

					//而这个正确划分的面就藏在在原划分中，怎么通过这个公共边把它找出来？一个边最多对应两个面，其中一个面只有两种可能：（1）已识别的重复面，（2）已识别的跨边面；
					//而另一个面就是我们要找的面，这个面更准确的讲，是位于重新三角化前的风险面集合中，但是这个集合又包含了一些不安全面，好像陷入了死循环。。。

					//但是好像并没有：
					//1. 如果这条边在风险面集合中所处的面被重新划分（flip）,而且是被flip的那条边，那这条边就不会出现在已识别的跨边面或者已识别的重复面中；
					//2. 而如果这条边没有被flip，而是被侵入，那么这条边也不会出现在已识别的跨边面或者已识别的重复面中；
					//3. 如果这条边是作为被侵入面或者被flip面的其他边，那么这条边一定是 已识别的重复面或者是已识别的跨边面上的边
					//4. 对于重复面上的边，从已识别的重复面上的边出发，在原划分集合中搜索具有这条边的面，要么找到它自己，要么就找到它真正的邻居
					//5. 对于跨边面上的边，从已识别的跨边面上的边出发，在原划分集合中搜索具有这条边的面，可能有两种情况，（1）之后要被拆掉的面（2）真正的邻居， 要如何区分呢？
					//5.1 将搜索范围中的将要被拆掉的三角形提前删掉，单独删除，好像不太好删除
					//5.2 从另外一个角度，要被拆掉的三角形都是跨边三角形的顶点，把所有的顶点组成集合，任何三个顶点位于这个集合之内的三角形，都不应该在我的搜索考虑范围之内
				}
			}
		}
	}

	//分析五：关于非安全面的选择（真的是最后一个分析了）
	//在凸包的外围，再加入一个点，如何影响
	//1. （p1,p2,p3,p4）组成凸包边界，加入顶点5后，会影响到几个面：除了(p3,p4), 为了形成凸包，(p2,p3)可能会受影响
	//问题：在凸包的外围再加入一个点，侵入到底有没有厚度，还是说只影响到第一圈???

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

	//1. 对与隔壁分区有相交风险的面进行识别
	for (auto face_handle : section1->m_dt.finite_face_handles())
	{
		//计算这个面的外接圆 和 它的左右两个分区是否有交叉，如果有, 将这个面认为是有风险的面
		if (!section1->check_circle_intersect(face_handle, seg, SSectionData::Near_Right)) {
			//saveFaceIndexs(face_handle);
		}
	}
	for (auto face_handle : section2->m_dt.finite_face_handles())
	{
		//计算这个面的外接圆 和 它的左右两个分区是否有交叉，如果有, 将这个面认为是有风险的面
		if (!section2->check_circle_intersect(face_handle, seg, SSectionData::Near_Left)) {
			//saveFaceIndexs(face_handle);
		}
	}

	//对非安全面重新进行三角化
	deplicateSection.insert_risk_delaunay(*section1, SSectionData::Near_Right);
	deplicateSection.insert_risk_delaunay(*section2, SSectionData::Near_Left);
	deplicateSection.insert_tangent_delaunay(*section1, *section2);

	//挑选出安全面，并检查跨边面和重复，为合并做准备
	for (Face_handle face_handle : deplicateSection.m_dt.finite_face_handles())
	{
		//saveFaceIndexs(face_handle);

		//检查是否跨边的面, 与合并前的面重复的面
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

	//合并（将重新三角化的面与合并前的进行比对）
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

	//检查是否有新的合并任务（添加合并任务到线程池），检查是否完成全部工作（停止线程池）
	checkMergingAndFinish();
}

void Delaunay_Mutithread_CGAL::delaunay_CGAL(SSectionData *sectionData)
{
	if (!sectionData || sectionData->m_vecDataIttrs.empty())
		return;

	//1. 先进行一次三角化
	for (auto &ittr : sectionData->m_vecDataIttrs)
	{
		//也可以封装到sectionData里面
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

	//检查是否有新的合并任务（添加合并任务到线程池），检查是否完成全部工作（停止线程池）
	checkMergingAndFinish();

	return;
}

void Delaunay_Mutithread_CGAL::checkMergingAndFinish()
{
	//检查是否有可以合并的分区（>=2）,如果有，就添加到任务队列
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

	//检查是否完成了全部工作
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
				|| m_vecSectionData[j].m_status == SSectionData::Status_Merging)
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

	//m_vecSectionData分区的数据是线程启动前就已经加入了
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
	//存在的缺陷：
	//1. fSampleRatio, iThreadNum参数选择合理性在调用前，自己处理确认，不在这里做调整优化，发现不满足，直接返回
	//2. 如果能够提前知道最大、最小值范围，划分会更加均衡，目前最大最小值也是有采样的获得的
	//3. 对于少数据量，一个线程只有一个点，这样没有意义

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
			m_vecSectionData[0].m_status = SSectionData::Status_Tri;
			m_vecSectionData[0].m_iSampleEndX = 1;	 //暂时只考虑横向的划分
			m_vecSectionData[0].boundingBox(ittr->x(), ittr->y());
			m_vecSectionData[0].m_vecDataIttrs.push_back(ittr);
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
		m_vecSectionData[i].m_status = SSectionData::Status_Tri;
	}

	//1. 随机采样确认数据分布(如何防止溢出啊)
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

	//2. 根据分布计算合适区域划分
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
			double low = dMin_x + slotIndex * slot_distance; //就这样吧
			double high = dMin_x + (slotIndex + 1) * slot_distance;
			if (!vecDivRanges.empty()) {
				vecDivRanges[rangeIndex-1].second = low;
			}
			vecDivRanges.push_back(DivRange(low, high));
			++rangeIndex;
			break;
		}
	}

	//3. 正式划分
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
	//目的：生成100%无重复数据，约束线无相交的数据集
	int iRowTotal = 10;
	int iColTotal = 10;
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

	if (iRowTotal < 10) iRowTotal = 10;
	if (iColTotal < 10) iColTotal = 10;

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
