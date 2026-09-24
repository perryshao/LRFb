//---------------------------------------------------------------------------

#ifndef Unit_recoH
#define Unit_recoH

typedef int TypeReco;
class TForm1; // test c++ builder
#include "mini_math_kernel.h"
#include "MemoryPool.h"
#include "MTMemoryPool.h"

#define MAXPOINT 1000

class Reco
{
    friend TForm1;

  protected:
    Q Ep(int i_u, int i_l);
    point U[MAXPOINT];
    point L[MAXPOINT];
    int nU, nL;
    int iU, iL;
    int lastX;

    void InsertPoint(Q &curEp, point &NP);
    void Rebuild(Lpoint LL, int &n, int miny, int maxy);
    bool back;

  public:
    Reco();
    void Init(); // new sequence of recognition
    void Normal(int &a, int &b);

    // we prefer an horizontal segment of 3 pixels and thickness 1 instead of
    // a vertical segment of thickness 3

    bool valid();
    virtual point Sym(point &P) = 0;
    void Insert(Q &curEp, point &NP);
    void GiveP(Lpoint LL, int &n, int miny, int maxy);

    // define the new allocation and memory pool
    void *operator new(size_t size)
    {
        void *p = s_pool->Alloc(size);
        return p;
    }

    void operator delete(void *p, size_t size)
    {
        s_pool->Free(p);
    }

    static void NewPool()
    {
        //s_pool = new CMemoryPool<CTest>;
        s_pool = new CMTMemoryPool<CMemoryPool<Reco>, CCriticalSection>;
    }

    static void DeletePool()
    {
        delete s_pool;
        s_pool = NULL;
    }

    //static CMemoryPool<CTest>* s_pool;
    static CMTMemoryPool<CMemoryPool<Reco>, CCriticalSection> *s_pool;
};

class RecoLR : public Reco
{
  public:
    point Sym(point &P)
    {
        return point(P.x, P.y);
    }
};
class RecoRL : public Reco
{
  public:
    point Sym(point &P)
    {
        return point(-P.x, P.y);
    }
};
class RecoBT : public Reco
{
  public:
    point Sym(point &P)
    {
        return point(P.y, P.x);
    }
};
class RecoTB : public Reco
{
  public:
    point Sym(point &P)
    {
        return point(-P.y, -P.x);
    }
};

#endif
