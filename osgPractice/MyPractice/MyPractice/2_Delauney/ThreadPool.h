#pragma once

#include <thread>
#include <deque>
#include <mutex>

#include "../pch.h"

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

	void detach() 
	{
		for (auto& td : m_vecThreads) {
			td.detach();
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
		if (m_qTasks.empty())
			return false;

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

class CLog
{
public:
	static CLog& getInstance()
	{
		static CLog log;
		return log;
	}

	~CLog()
	{
		for (int i = 0; i < m_vecTasks.size(); ++i)
		{
			delete m_vecTasks[i];
			m_vecTasks[i] = nullptr;
		}
		ofs.close();
	}

	void LOG(const CString& strInfo)
	{
		task* pTask = new LogTask(strInfo);
		m_vecTasks.push_back(pTask);
		m_ThreadPool.addTask(pTask);
	}

	bool Logout(const CString& strInfo)
	{
		if (!ofs.is_open())
			return false;

		ofs << (LPCTSTR)strInfo << (LPCTSTR)(_T("\n"));
	}

	struct LogTask : public task
	{
		LogTask(const CString& strInfo = _T(""))
		{
			m_strInfo = strInfo;
		}

		void operator()() override
		{
			CLog::getInstance().Logout(m_strInfo);
		}

		CString m_strInfo;
	};

private:
	CLog() : m_filepath(_T("D:\\log.txt"))	//先写死了
	{
		ofs.open(m_filepath);
		m_ThreadPool.setThreadNum(1);
		m_ThreadPool.start();
		m_ThreadPool.detach();
	}

	CString m_filepath;
	std::wofstream ofs;
	CThreadPool m_ThreadPool;	//还是用的轮询，效率比较低，后面再看看，先用来排查bug
	std::vector<task*> m_vecTasks;
};