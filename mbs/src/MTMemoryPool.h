#ifdef WIN32
#include <windows.h>
typedef CRITICAL_SECTION MUTEXTYPE;
#define INITMUTEX(m_cs) InitializeCriticalSection(&m_cs)
#define DELMUTEX(m_cs) DeleteCriticalSection(&m_cs)
#define LOCK(m_cs) EnterCriticalSection(&m_cs)
#define UNLOCK(m_cs) LeaveCriticalSection(&m_cs)
#else
#include <pthread.h>
typedef pthread_mutex_t MUTEXTYPE;
#define INITMUTEX(m_cs) pthread_mutex_init(&m_cs, NULL)
#define DELMUTEX(m_cs) pthread_mutex_destroy(&m_cs)
#define LOCK(m_cs) pthread_mutex_lock(&m_cs)
#define UNLOCK(m_cs) pthread_mutex_unlock(&m_cs)
#endif

class CCriticalSection
{
  public:
    CCriticalSection()
    {
        INITMUTEX(m_cs);
    }

    ~CCriticalSection()
    {
        DELMUTEX(m_cs);
    }

    void Lock()
    {
        LOCK(m_cs);
    }

    void Unlock()
    {
        UNLOCK(m_cs);
    }

  protected:
    MUTEXTYPE m_cs;
};

template <typename POOLTYPE, typename LOCKTYPE> class CMTMemoryPool
{
  public:
    void *Alloc(unsigned int size)
    {
        void *p = NULL;
        m_lock.Lock();
        p = m_pool.Alloc(size);
        m_lock.Unlock();

        return p;
    }

    void Free(void *p)
    {
        m_lock.Lock();
        m_pool.Free(p);
        m_lock.Unlock();
    }

  private:
    POOLTYPE m_pool;
    LOCKTYPE m_lock;
};
