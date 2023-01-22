#pragma once

#include <thread>
#include <deque>
#include <mutex>

struct task
{
	enum EPriority 
	{
		ePriority_High,
		ePriority_Low
	};

	virtual void operator()()
	{
		int a = 1;
	}
};

class CThreadPool
{
public:
	CThreadPool(int iThreadNum = 1)
		: m_bDone(false)
	{
		if (m_iThreadNum <= 0)
			m_iThreadNum = 1;
		int iMaxThread = std::thread::hardware_concurrency()  * 2;
		if (m_iThreadNum > iMaxThread - 1)
			m_iThreadNum = iMaxThread;
	}

	~CThreadPool() 
	{
		m_bDone = true;
	}

	bool empty() 
	{
		std::lock_guard<std::recursive_mutex> lock(m_mutex);
		return m_qTasks.empty();
	}

	void start() 
	{
		try
		{
			for (int i = 0; i < m_iThreadNum; ++i) {
				m_vecThreads.push_back(std::thread(&CThreadPool::worker_thread, this));
			}
		}
		catch (...)
		{
			m_bDone = true;
			//throw;
		}
	}

	void stop() //在task中调用
	{
		m_bDone = true;
	}

	void join()
	{
		for (auto &td : m_vecThreads){
			if (td.joinable())
				td.join();
		}
	}

	void setThreadNum(int iThreamNum) 
	{
		m_iThreadNum = iThreamNum;
	}

	void addTask(task* pTask, task::EPriority pri = task::ePriority_Low)
	{
		if (!pTask)
			return;
		
		std::lock_guard<std::recursive_mutex> lock(m_mutex);

		if (pri == task::ePriority_High){
			m_qTasks.push_front(pTask);
		} else {
			m_qTasks.push_back(pTask);
		}
	}

protected:

	bool popTask(task **pTask)
	{
		std::lock_guard<std::recursive_mutex> lock(m_mutex);
		if (m_qTasks.empty())
			return false;

		*pTask = m_qTasks.front();
		m_qTasks.pop_front();
		return true;
	}

	void worker_thread()
	{
		while (!m_bDone){
			task* pTask[1];
			if (popTask(pTask)){
				if (pTask && *pTask) 
					(*(*pTask))();
			}
			else{
				std::this_thread::yield();
			}
		}
	}

private:
	std::vector<std::thread> m_vecThreads;
	std::deque<task*> m_qTasks;			//后面再考虑智能指针
	int m_iThreadNum;
	std::recursive_mutex m_mutex;
	bool m_bDone;
};